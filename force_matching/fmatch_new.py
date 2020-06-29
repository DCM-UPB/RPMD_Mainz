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
    xyzdat=np.empty(nout, dtype=np.float64)
    it=0
    for line in xyzfile.readlines():
        linedat=line.split()
        if linedat[0]=="O" or linedat[0]=="H":         
            xyzdat[it]   = float(linedat[1])
            xyzdat[it+1] = float(linedat[2])
            xyzdat[it+2] = float(linedat[3])
            it+=3

    return xyzdat
    
'''        VP FFM Residual Class (callable)
Fit function (Buckingham): lp[0]*F_Coul(nlp[0]) + lp[1]*F_LJ(nlp[1, 5]) + lp[2]*F_stretch(nlp[2:4]) + lp[3]*F_bend(nlp[4])
Fit function (LJ):         lp[0]*F_Coul(nlp[0]) + lp[1]*F_LJ(nlp[1]) + lp[2]*F_stretch(nlp[2:4]) + lp[3]*F_bend(nlp[4]) '''
class vpres_ffm( object ):  
    def __init__( self, natom, ndata, useBuckingham, linp0, traj_fname, aifrc_fname, parout_fname):
        ndim = 3
        self.ndim = ndim
        self.natom = natom
        self.ndata = ndata
        self.narr = ndim*natom*ndata
        self.counteval = 0

        self.useBuckingham = useBuckingham
        if useBuckingham:                   # nonlinear parameters
            self.nolpars = np.zeros(6)
        else:
            self.nolpars = np.zeros(5)
        self.funstore = np.zeros(self.narr)  # stores the total forces
        self.phistore = np.zeros((self.narr, 4))  # stores the force components

        self.linpars = np.array(linp0)              # linear parameters (expected length 4!)
        self.parout_fname = parout_fname
        self.traj_fname = traj_fname

        self.aifrc = read_xyzdat(aifrc_fname, self.narr)

        stdoutfile = open("stdout", 'w')
        stdoutfile.write("")
        stdoutfile.write("-------        FFM Script Output        -------")
        stdoutfile.write("")
        stdoutfile.close()

    def __call__( self, params):
        '''params is expected to be of the Parameters type (lmfit)'''
        hparams = np.empty(len(self.nolpars))
        hparams[0] = params["alpha"].value
        hparams[1] = params["oo_sig"].value
        hparams[2] = params["alp"].value
        hparams[3] = params["reoh"].value
        hparams[4] = params["thetad"].value
        if self.useBuckingham:
            hparams[5] = params["oo_gam"].value

        if not np.all(hparams == self.nolpars):  # recalculate force terms
            self.nolpars = hparams
            self.fitfun_ffm_phi()

        #self.linpars = spo.lsq_linear(self.phistore, self.aifrc, (self.linpar_bounds[:, 0], self.linpar_bounds[:, 1]))
        #self.linpars = self.linpars.x # the actual fit result
        self.linpars = spo.nnls(self.phistore, self.aifrc)
        self.linpars = self.linpars[0]  # the actual fit result

        print("non-linear params:", self.nolpars)
        print("linear params:", self.linpars)
        if np.any(self.linpars == 0):
            print("Error: One or more linear parameters reached 0 in the linear optimization step. Exiting.")
            quit()

        self.fitfun_ffm()

        print("res:", sum(np.power(self.aifrc - self.funstore, 2)))
        return self.aifrc - self.funstore
   
    def fitfun_ffm_phi(self):
        '''Expects lp[0:4],nlp[0:6],R[ndata,0:3],ndata,na'''
    
        outfile = open(self.parout_fname, 'w')
        outfile.write('&param')
        outfile.write('\nwmass = %.16f'%    32831.2525000000000003e0)
        outfile.write('\nomass = %.16f'%    29156.9471000000000001e0)
        outfile.write('\nhmass = %.16f'%    1837.1527000000000001e0)
        outfile.write('\nqo = %.16f'%     - np.sqrt(self.linpars[0]))
        outfile.write('\nalpha = %.16f'%    self.nolpars[0])
        outfile.write('\noo_sig = %.16f'%   self.nolpars[1])
        outfile.write('\noo_eps = %.16f'%   self.linpars[1])
        if self.useBuckingham:
            outfile.write('\noo_gam = %.16f'%   self.nolpars[5])
        else:
            outfile.write('\noo_gam = %.16f' % 0.e0)
        outfile.write('\nthetad = %.16f'%   self.nolpars[4])
        outfile.write('\nreoh = %.16f'%     self.nolpars[3])
        outfile.write('\napot = %.16f'%     self.linpars[2])
        outfile.write('\nbpot = %.16f'%     self.linpars[3])
        outfile.write('\nalp = %.16f'%      self.nolpars[2])
        # UNUSED PARAMETERS
        outfile.write('\nalpb = %.16f'%     0.e0)
        outfile.write('\nwm = %.16f'%       0.e0)
        outfile.write('\nwh = %.16f'%       0.e0)
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
        
    def fitfun_ffm(self):
        ''' Call fitfun_ffm_phi before calling this (unless phistore is already up to date)!'''
        
        self.funstore[:] = 0.0
        for it in range(4):
            self.funstore[:] += self.linpars[it]*self.phistore[:, it]
        
'''                      --------                      '''


'''        
ndata=10       
exlinpars=np.array((0.3,0.7))
exnolpars=np.array((-0.3,-0.6))

Rtest = np.empty((ndata,6))
for it1 in range(ndata):
    for it2 in range(6):
        Rtest[it1,it2]=rd.random()
        
dataf=fitfun_test(exlinpars,exnolpars,Rtest)
meanfabs=np.mean(abs(dataf))
print(dataf)
print(np.size(dataf))
noisefac=0.05*meanfabs
for it1 in range(ndata*6):
    dataf[it1]+=noisefac*rd.random()
print(dataf)
    
testvp=vpres_test(2,ndata)    
residual=testvp(exnolpars,dataf,Rtest)
print(residual)

fitpars=spo.leastsq(testvp,[-0.1,1.0],args=(dataf, Rtest))
print(fitpars)
'''

linp0 = np.array([1.5, 0.0001, 0.15, 0.06])
useBuckingham = False

fitparams = lmf.Parameters()
fitparams.add('alpha',      value=0.7,   min=0.500, max=1.000)
fitparams.add('oo_sig',     value=5.5,   min=5.0,   max=7.0)
if useBuckingham:
    fitparams.add('oo_gam', value=13.0,  min=0.010, max=50.000)
fitparams.add('alp',        value=1.35,  min=0.010, max=100.0)
fitparams.add('reoh',       value=1.8,   min=1.600, max=2.000)
fitparams.add('thetad',     value=104.5, min=100.00, max=115.0)

natom = 375
ndata = 1500
traj_fname = "all_tray.xyz"
aifrc_fname = "FORCES-TPSS-D3-ALL.frc"
parout_fname = "param_TPSS-D3-new"
fitftol = 1.e-12

vpres_ffm_obj = vpres_ffm(natom, ndata, useBuckingham, linp0, traj_fname, aifrc_fname, parout_fname)
residual = vpres_ffm_obj(fitparams)

fitout = lmf.minimize(vpres_ffm_obj, fitparams, ftol=fitftol)
print(vpres_ffm_obj.linpars)
print(vpres_ffm_obj.nolpars)

lmf.report_fit(fitout)
