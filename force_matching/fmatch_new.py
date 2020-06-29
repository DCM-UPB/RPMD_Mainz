'''
Water Model Force Matching by Variable Projection

Author: Jan Kessler
Originally created in 07/2013
Rebuilt in 06/2020
'''
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

    return xyzdat
    

class varpro_ffm_res(object):
    '''        Variable Projection Force Field Matching Residual Class (callable)
    Fit function (Buckingham): qo^2*F_Coul(alpha) + ooeps*F_BH(oo_sig, oo_gam) + apot*F_stretch(alp, reoh) + bpot*F_bend(thetad)
    Fit function (LJ):         qo^2*F_Coul(alpha) + oo_eps*F_LJ(oo_sig) + apot*F_stretch(alp, reoh) + bpot*F_bend(thetad) '''
    def __init__(self, natom, ndata, traj_fname, aifrc_fname, parout_fname):
        ndim = 3
        self.narr = ndim*natom*ndata  # number of values in force & trajectory files

        self.linpars = np.array([1.0, 0.0005, 0.15, 0.07])  # stores linear parameters, use some initial values for the first force evaluation (values have no influence on optimization!!)
        self.funstore = np.zeros(self.narr)  # stores the total forces
        self.phistore = np.zeros((self.narr, 4))  # stores the force components (divided by respective linear factors)

        self.traj_fname = traj_fname  # the file containing positions for which the reference forces are given
        self.aifrc = read_xyzdat(aifrc_fname, self.narr)  # the reference forces
        self.parout_fname = parout_fname  # the file used to pass the FF parameters to the force code

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
    
        outfile = open(self.parout_fname, 'w')
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
    
        args = ("./wff_forces.x", "input", self.parout_fname, self.traj_fname)
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

    def _opt_linpars(self):
        '''Optimize for a new set of linear parameters'''
        # left alternative bounded version commented out:
        # self.linpars = spo.lsq_linear(self.phistore, self.aifrc, (self.linpar_bounds[:, 0], self.linpar_bounds[:, 1]))
        # self.linpars = self.linpars.x # the actual fit result
        linfit = spo.nnls(self.phistore, self.aifrc)
        self.linpars = linfit[0]  # the actual fit result

        if np.any(self.linpars == 0):
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


useBuckingham = False

fitparams = lmf.Parameters()
fitparams.add('alpha',      value=0.7,   min=0.500, max=1.000)
fitparams.add('oo_sig',     value=6.0,   min=1.0,   max=10.0)
if useBuckingham:
    fitparams.add('oo_gam', value=13.0,  min=0.010, max=50.000)
fitparams.add('alp',        value=1.35,  min=0.010, max=100.0)
fitparams.add('reoh',       value=1.8,   min=1.600, max=2.000)
fitparams.add('thetad',     value=104.5, min=100.00, max=115.0)

natom = 375
ndata = 1500
traj_fname = "all_tray.xyz"
aifrc_fname = "FORCES-PBE-ALL.frc"
parout_fname = "param_PBE-new"
fitftol = 1.e-12

varpro_ffm_obj = varpro_ffm_res(natom, ndata, traj_fname, aifrc_fname, parout_fname)
residual = varpro_ffm_obj(fitparams)

fitout = lmf.minimize(varpro_ffm_obj, fitparams, ftol=fitftol)

print('Final optimization result:')
print()
varpro_ffm_obj.report_linear()
print()
lmf.report_fit(fitout)
