'''
Water Model Force Matching by Variable Projection

Author: Jan Kessler
Originally created in 07/2013
Rebuilt in 06/2020
'''
import sys
import argparse
import numpy as np
import numpy.linalg as npla
import scipy.linalg as spla
import scipy.optimize as spo
import lmfit as lmf
import random as rd
import subprocess

                    
'''                ----                '''
'''    ----------------------------------------------------    '''
'''        Variable Projection for Force Field Matching        '''
'''    ----------------------------------------------------    '''

def read_xyzdat(filename, nout):
    xyzfile = open(filename, 'r')
    xyzdat = np.empty(nout, dtype=np.float64)
    it = 0
    for line in xyzfile.readlines():
        linedat=line.split()
        if linedat[0] == "O" or linedat[0] == "H":
            xyzdat[it]   = float(linedat[1])
            xyzdat[it+1] = float(linedat[2])
            xyzdat[it+2] = float(linedat[3])
            it += 3

    if nout != len(xyzdat):
        print("Error: XYZ file did not contain the expected amount of positions/forces.")
        quit()

    return xyzdat
    

class varpro_ffm_res(object):
    '''        Variable Projection Force Field Matching Residual Class (callable)
    Fit function (Buckingham): qo^2*F_Coul(alpha) + ooeps*F_BH(oo_sig, oo_gam) + apot*F_stretch(alp, reoh) + bpot*F_bend(thetad)
    Fit function (LJ):         qo^2*F_Coul(alpha) + oo_eps*F_LJ(oo_sig) + apot*F_stretch(alp, reoh) + bpot*F_bend(thetad) '''
    def __init__(self, natom, ndata, traj_fname, aifrc_fname, linpar_bounds):
        ndim = 3
        self.narr = ndim*natom*ndata  # number of values in force & trajectory files

        self.funstore = np.zeros(self.narr)  # stores the total forces
        self.phistore = np.zeros((self.narr, 4))  # stores the force components (divided by respective linear factors)

        self.traj_fname = traj_fname  # the file containing positions for which the reference forces are given
        self.aifrc = read_xyzdat(aifrc_fname, self.narr)  # the reference forces

        self.linpars = np.array([1.0, 0.0005, 0.15, 0.07])  # stores linear parameters, use some initial values for the first force evaluation (values have no influence on optimization!!)

        # handle linear parameter boundaries
        flags_nobounds = [bpair[0] == 0. and bpair[1] == 0. for bpair in linpar_bounds]
        self.linpar_lbounds = np.empty(4)
        self.linpar_ubounds = np.empty(4)
        for it, flag in enumerate(flags_nobounds):
            if flag: # linear variable may have any non-zero value
                self.linpar_lbounds[it] = 0.
                self.linpar_ubounds[it] = np.inf
            else: # linear variable stays within given boundaries
                self.linpar_lbounds[it] = linpar_bounds[it][0]
                self.linpar_ubounds[it] = linpar_bounds[it][1]
        self.flag_linpar_anybounds = np.any(flags_nobounds)

        self.counteval = 0  # counts number of force program executions

        # print header of output file
        stdoutfile = open("stdout", 'w')
        stdoutfile.write("")
        stdoutfile.write("-------        FFM Script Output        -------")
        stdoutfile.write("")
        stdoutfile.close()

    def __call__(self, params):
        '''params is expected to be of the Parameters type (lmfit)'''

        # execute force computation and linear optimization
        self._recompute_force_components(params)
        self._opt_linpars()
        self._recompute_force_superposition()
        self._report(params)

        # return residuals
        return self.aifrc - self.funstore
   
    def _recompute_force_components(self, params):
        '''Takes lmfit params (+ self.linpars) and calls force code to compute corresponding force components'''
    
        outfile = open('params_fit', 'w')
        outfile.write('&param')
        outfile.write('\nwmass = %.16f' %      32831.2525000000000003e0)
        outfile.write('\nomass = %.16f' %      29156.9471000000000001e0)
        outfile.write('\nhmass = %.16f' %      1837.1527000000000001e0)
        outfile.write('\nqo = %.16f' %       - np.sqrt(self.linpars[0]))
        outfile.write('\nalpha = %.16f' %      params['alpha'].value)
        outfile.write('\noo_sig = %.16f' %     params['oo_sig'].value)
        outfile.write('\noo_eps = %.16f' %     self.linpars[1])
        if "oo_gam" in params:  # then use Buckingham potential instead of standard Lennard-Jones
            outfile.write('\noo_gam = %.16f' % params['oo_gam'].value)
        else:
            outfile.write('\noo_gam = %.16f' % 0.e0)
        outfile.write('\nthetad = %.16f' %     params['thetad'].value)
        outfile.write('\nreoh = %.16f' %       params['reoh'].value)
        outfile.write('\napot = %.16f' %       self.linpars[2])
        outfile.write('\nbpot = %.16f' %       self.linpars[3])
        outfile.write('\nalp = %.16f' %        params['alp'].value)
        # UNUSED PARAMETERS
        outfile.write('\nalpb = %.16f' %       0.e0)
        outfile.write('\nwm = %.16f' %         0.e0)
        outfile.write('\nwh = %.16f' %         0.e0)
        outfile.write('\n&end\n')
        outfile.close()
    
        args = ("./wff_forces.x", "input", "params_fit", self.traj_fname)
        popen = subprocess.Popen(args, stdout=subprocess.PIPE)
        popen.wait()
        output = popen.stdout.read()
        
        self.counteval += 1
        
        stdoutfile = open("stdout", 'a')
        stdoutfile.write("---    Output of Evaluation No. " + str(self.counteval) + "    ---")
        stdoutfile.write(str(output))
        stdoutfile.write("")
        stdoutfile.close()
        
        self.phistore[:, 0] = read_xyzdat("vmd.frc1", self.narr)/self.linpars[0]
        self.phistore[:, 1] = read_xyzdat("vmd.frc2", self.narr)/self.linpars[1]
        self.phistore[:, 2] = read_xyzdat("vmd.frc3", self.narr)/self.linpars[2]
        self.phistore[:, 3] = read_xyzdat("vmd.frc4", self.narr)/self.linpars[3]

    def _opt_linpars_nonzero(self):
        '''Optimize for a new set of non-zero linear parameters'''
        linfit = spo.nnls(self.phistore, self.aifrc)
        self.linpars = linfit[0]  # the actual fit result

    def _opt_linpars_bound(self):
        '''Optimize for a new set of bound linear parameters'''
        self.linpars = spo.lsq_linear(self.phistore, self.aifrc, (self.linpar_lbounds, self.linpar_ubounds))
        self.linpars = self.linpars.x # the actual fit result

    def _opt_linpars(self):
        '''Optimize for a new set of linear parameters'''
        if self.flag_linpar_anybounds:
            self._opt_linpars_bound()
        else:
            self._opt_linpars_nonzero()

        if np.any(self.linpars == 0):
            print(self.linpars)
            print("Error: One or more linear parameters reached 0 in the linear optimization step. Exiting.")
            quit()

    def _recompute_force_superposition(self):
        '''Used to recompute the forces when only the linear factors are changed.'''

        self.funstore[:] = 0.0
        for it in range(4):
            self.funstore[:] += self.linpars[it] * self.phistore[:, it]

    def report_linear(self):
        '''Print current set of linear parameters'''
        print("Linear Parameters:")
        print("qo^2:", self.linpars[0])
        print("oo_eps:", self.linpars[1])
        print("apot:", self.linpars[2])
        print("bpot:", self.linpars[3])

    def _report(self, params):
        '''Report current state of parameters / residual'''

        print()
        print("Optimization Step Report:")
        print()
        self.report_linear()
        print()
        print("Non-linear Parameters:")
        print("alpha:", params['alpha'].value)
        print("oo_sig:", params['oo_sig'].value)
        if "oo_gam" in params:
            print("oo_gam:", params['oo_gam'].value)
        print("alp:", params['alp'].value)
        print("reoh", params['reoh'].value)
        print("thetad", params['thetad'].value)
        print()
        print("Residual Sum:", sum(np.power(self.aifrc - self.funstore, 2)))
        print()
        
'''                      --------                      '''

# --- Setup CMDLine Parser ---
parser = argparse.ArgumentParser(prog='FFM', description='Fit 4-Site Water Model using Variable Projection, given a dataset of positions and corresponding ab-initio forces.')

parser.add_argument('--natom', required=True, type=int, help='Number of atoms in each of the provided snapshots (required)')
parser.add_argument('--ndata', required=True, type=int, help='Number of provided snapshots in the trajectory/force files (required)')

parser.add_argument('--traj', required=True, help='File containing the trajectory snapshots (required)')
parser.add_argument('--frcs', required=True, help='File containing the force snapshots (required)')

# flags
parser.add_argument('--useBuck', action='store_true', help='Use Buckingham potential instead of Lennard-Jones')

# non-linear model parameter init, min and max
parser.add_argument('--alpha', type=float, nargs=3, metavar=('init', 'min', 'max'), default = [0.7, 0.5, 1.0], help='Init, Min and Max for alpha parameter')
parser.add_argument('--oo_sig', type=float, nargs=3, metavar=('init', 'min', 'max'), default = [6.0, 1.0, 10.0], help='Init, Min and Max for oo_sig parameter')
parser.add_argument('--oo_gam', type=float, nargs=3, metavar=('init', 'min', 'max'), default = [12.0, 0.1, 50.0], help='Init, Min and Max for oo_gam parameter')
parser.add_argument('--alp', type=float, nargs=3, metavar=('init', 'min', 'max'), default = [1.35, 0.01, 100.], help='Init, Min and Max for alp parameter')
parser.add_argument('--reoh', type=float, nargs=3, metavar=('init', 'min', 'max'), default = [1.8, 1.6, 2.], help='Init, Min and Max for reoh parameter')
parser.add_argument('--thetad', type=float, nargs=3, metavar=('init', 'min', 'max'), default = [104.5, 100., 115.], help='Init, Min and Max for thetad parameter')

# linear model parameter min and max
parser.add_argument('--qo2', type=float, nargs=2, metavar=('min', 'max'), default = [0., 0.], help='Min and Max for qo2 (qo*qo) parameter (default: unbound, non-negative)')
parser.add_argument('--oo_eps', type=float, nargs=2, metavar=('min', 'max'), default = [0., 0.], help='Min and Max for oo_eps parameter (default: unbound, non-negative)')
parser.add_argument('--apot', type=float, nargs=2, metavar=('min', 'max'), default = [0., 0.], help='Min and Max for apot parameter (default: unbound, non-negative)')
parser.add_argument('--bpot', type=float, nargs=2, metavar=('min', 'max'), default = [0., 0.], help='Min and Max for bpot parameter (default: unbound, non-negative)')

# fitting configuration
parser.add_argument('--fit_method', default='least_squares', choices=['least_squares', 'trf', 'dogbox'], help='Method used for the bounded non-linear fitting (defaults to Trust-Region-Reflective least squares optimizer)')
parser.add_argument('--fit_max_nfev', type=int, default=None, help='Maximum number of function evaluations for non-linear fit (stopping criterium, only working with least_squares method)')
parser.add_argument('--fit_ftol', type=float, default=1.e-8, help='Fit function tolerance for non-linear fit (stopping criterium, only working with least_squares method)')

parser.print_help()
args = parser.parse_args()

# Setup lmf parameter object
fitparams = lmf.Parameters()
fitparams.add('alpha',      value=args.alpha[0], min=args.alpha[1], max=args.alpha[2])
fitparams.add('oo_sig',     value=args.oo_sig[0], min=args.oo_sig[1], max=args.oo_sig[2])
if args.useBuck:
    fitparams.add('oo_gam', value=args.oo_gam[0], min=args.oo_gam[1], max=args.oo_gam[2])
fitparams.add('alp',        value=args.alp[0], min=args.alp[1], max=args.alp[2])
fitparams.add('reoh',       value=args.reoh[0], min=args.reoh[1], max=args.reoh[2])
fitparams.add('thetad',     value=args.thetad[0], min=args.thetad[1], max=args.thetad[2])

# Setup variable projection (VP)
varpro_ffm_obj = varpro_ffm_res(args.natom, args.ndata, args.traj, args.frcs, [args.qo2, args.oo_eps, args.apot, args.bpot])
residual = varpro_ffm_obj(fitparams)

# Execute VP
if args.fit_method=='least_squares':
    fitout = lmf.minimize(varpro_ffm_obj, fitparams, method=args.fit_method, nan_policy='raise', max_nfev=args.fit_max_nfev, ftol=args.fit_ftol)
else:
    fitout = lmf.minimize(varpro_ffm_obj, fitparams, method=args.fit_method, nan_policy='raise')

# Report
print('Final optimization result:')
print()
varpro_ffm_obj.report_linear()
print()
lmf.report_fit(fitout)
