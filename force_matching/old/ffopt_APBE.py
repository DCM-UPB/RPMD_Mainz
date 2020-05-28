'''
Created on 24.07.2013

@author: Jan Kessler

Water Model Force Matching by Variable Projection 
'''
import numpy as np
import numpy.linalg as npla
import scipy.linalg as spla
import scipy.optimize as spo
import lmfit as lmf
import random as rd
import subprocess




from scipy.linalg import qr
''' Copied out of source from https://bitbucket.org/tillsten/pymls/ '''

# -*- coding: utf-8 -*-
"""
Created on Thu May 05 20:02:00 2011

@author: Tillsten
"""

eps = np.finfo(float).eps

def mls(B, v, umin, umax, Wv=None, Wu=None, ud=None, u=None, W=None, imax=100):
    """
mls - Control allocation using minimal least squares.

[u,W,iter] = mls_alloc(B,v,umin,umax,[Wv,Wu,ud,u0,W0,imax])

 Solves the bounded sequential least-squares problem

   min ||Wu(u-ud)||   subj. to   u in M

 where M is the set of control signals solving

   min ||Wv(Bu-v)||   subj. to   umin <= u <= umax

 using a two stage active set method. Wu must be diagonal since the
 problem is reformulated as a minimal least squares problem. The
 implementation does not handle the case of coplanar controls.

  Inputs:
  -------
 B     control effectiveness matrix (k x m)
 v     commanded virtual control (k x 1)
 umin  lower position limits (m x 1)
 umax  upper position limits (m x 1)
 Wv    virtual control weighting matrix (k x k) [I]
 Wu    control weighting matrix (m x m), diagonal [I]
 ud    desired control (m x 1) [0]
 u0    initial point (m x 1)
 W0    initial working set (m x 1) [empty]
 imax  max no. of iterations [100]

  Outputs:
  -------
 u     optimal control
 W     optimal active set
 iter  no. of iterations (= no. of changes in the working set + 1)

                           0 if u_i not saturated
 Active set syntax: W_i = -1 if u_i = umin_i
                          +1 if u_i = umax_i

Directly Based on the code from:
    Ola Haerkegard, www.control.isy.liu.se/~ola
    see licsence.

"""

    #k = number of virtual controls
    #m = number of variables (actuators)
    k, m = B.shape
    if u == None:
        u = np.mean(umin + umax, 0)[:, None]
    if W == None:
        W = np.zeros((m, 1))
    if ud == None:
        ud = np.zeros((m, 1))
    if Wu == None:
        Wu = np.eye(m)
    if Wv == None:
        Wv = np.eye(k)
    phase = 1
    #Reformulate as a minimal least squares problem. See 2002-03-08 (1).
    A = Wv.dot(B).dot(np.linalg.pinv(Wu))
    b = Wv.dot(v - B.dot(ud))

    xmin = (umin - ud).flatten()
    xmax = (umax - ud).flatten()
    # Compute initial point and residual.
    x = Wu.dot(u - ud)
    r = np.atleast_2d(A.dot(x) - b)

    #Determine indeces of free variables
    i_free = (W == 0).flatten()
    m_free = np.sum(i_free)

    for i in range(imax):
    #print 'Iter: ', i
        if phase == 1:
            A_free = A[:, i_free]
            if m_free <= k:
                if m_free > 0:
                    p_free = np.linalg.lstsq(-A_free, r)[0]
            else:
                q1, r1 = qr(A_free.T)
                p_free = -q1.dot(np.solve(r1.T, r))
            p = np.zeros((m, 1))
            if A.shape[1] > 1:
                p[i_free] = p_free
            else:
                p[i_free] = p_free.flatten()

        else:
            i_fixed = np.logical_not(i_free)
            m_fixed = m - m_free

            if m_fixed > 0:
                HT = U[i_fixed.squeeze(), :].T
                V, Rtot = qr(np.atleast_2d(HT))
                V1 = V[:, :m_fixed]
                V2 = V[:, m_fixed + 1:]
                R = Rtot[:, m_fixed]
            else:
                V, Rtot = np.array([[]]), np.array([[]])
                V1 = V2 = R = V.T

            s = -V2.T.dot(z)
            pz = V2.dot(s)
            p = U.dot(pz)

        x_opt = x + p
        infeasible = np.logical_or(x_opt < xmin, x_opt > xmax)
        if not np.any(infeasible[i_free]):
            x = x_opt

            if phase == 1:
                r = r + A.dot(p)
            else:
                z = z + pz
            if phase == 1 and m_free >= k:
                phase = 2
                Utot, Stot = qr(A.T)
                U = Utot[:, k:]
                z = U.T.dot(x)
            else:
                lam = np.zeros((m, 1))
                if m_free < m:
                    if phase == 1:
                        g = A.T.dot(r)
                        lam = -W * g
                    else:
                        lam[i_fixed] = -W[i_fixed] * np.linalg.solve(R, V1.T.dot(z))
                    if np.all(lam >= -eps):
                        u = np.linalg.solve(Wu, x) + ud
                        return u

                lambda_neg, i_neg = np.min(lam), np.argmin(lam)
                W[i_neg] = 0
                i_free[i_neg] = True
                m_free += 1
        else:
            dist = np.ones(m)
            i_min = np.logical_and(i_free, p.flat < 0).flatten()
            i_max = np.logical_and(i_free, p.flat > 0).flatten()
            dist[i_min] = (xmin[i_min] - x[i_min]) / p[i_min]
            dist[i_max] = (xmax[i_max] - x[i_max]) / p[i_max]
            alpha, i_alpha = np.min(dist), np.argmin(dist)

            x = x + alpha * p
            if phase == 1:
                r = r + A.dot(alpha * p) #!!
            else:
                z = z + alpha * pz
            W[i_alpha] = np.sign(p[i_alpha])

            if i_free[i_alpha]:
                i_free[i_alpha] = False
                m_free -= 1

    u = np.linalg.solve(Wu, x) + ud
    return u

def bounded_lsq(A, b, lower_lim, upper_lim):
    """
    Minimizes:

    |Ax-b|_2

    for lower_lim<x<upper_lim.
    """
    return mls(A, b, lower_lim, upper_lim)




'''    Variable Projection Testing    '''


def fitfun_test_phi(alph,R):
    '''Expects alph[0:2],R[ndata,0:3]'''
    ndata=np.shape(R)[0]
    nout=ndata*6
    dvec  = np.empty((ndata,3))
    Rnorm = np.empty(ndata)
    
    for it in range(ndata):

        dvec[it]  = R[it,3:6]-R[it,0:3]
        Rnorm[it] = npla.norm(dvec[it,0:3])
    
    phi=np.zeros((nout,2))
    itr=0
    itrm=0
    for it in range(ndata):
        itr=it*6
        itrm=itr + 3
        
        phi[itr:itrm,0]   =  np.exp(alph[0]*Rnorm[it]) * dvec[it]
        phi[itr:itrm,1]   =  np.exp(alph[1]*Rnorm[it]) * dvec[it]
        phi[itrm:itrm+3,0] =  -phi[itr:itrm,0]      
        phi[itrm:itrm+3,1] =  -phi[itr:itrm,1]
        
    return phi

def fitfun_test(cs,alph,R):
    phi=fitfun_test_phi(alph,R)
    return cs[0]*phi[:,0] + cs[1]*phi[:,1]
    
'''        VP Test Residual Class (callable)
Fit function: [ c1 * exp(a1*|(r2-r1)|)
+ c2 * exp(a2*|(r2-r1)|) ] * (r2-r1)/|(r2-r1)| '''
class vpres_test( object ):  
    def __init__( self, natom, ndata):
        ndim = 3
        self.ndim = ndim
        self.natom = natom
        self.ndata = ndata
        
        self.linpars = np.empty(2)                  #linear parameters
        self.nolpars = np.empty(2)                  #nonlinear parameters
        self.funstore = np.zeros(ndim*natom*ndata)  #stores the total function
        self.phistore = np.zeros((ndim*natom*ndata,2))#stores the basis functions

    def __call__( self, params, datap, R):
        if all(params[0:2]==self.nolpars[0:2]):
            return datap-self.funstore
        else:
            self.nolpars = params
            self.phistore = fitfun_test_phi(self.nolpars,R)
            self.linpars=np.linalg.lstsq(self.phistore,datap)
            self.linpars=self.linpars[0]
            self.funstore=fitfun_test(self.linpars,self.nolpars,R)
            
            return datap-self.funstore
  
  
                    
'''                ----                '''
'''    ----------------------------------------------------    '''
'''        Variable Projection for Force Field Matching        '''
'''    ----------------------------------------------------    '''


def read_xyzdat(filename,nout):
    xyzfile = file(filename, 'r')
    xyzdat=np.empty(nout)
    it=0
    for line in xyzfile.readlines():
        linedat=line.split()
        if linedat[0]=="O" or linedat[0]=="H":         
            xyzdat[it]   = float(linedat[1])
            xyzdat[it+1] = float(linedat[2])
            xyzdat[it+2] = float(linedat[3])
            it+=3
    xyzout=np.array(xyzdat)    
    return xyzout   
    
'''        VP FFM Residual Class (callable)
Fit function: cs[0]*F_Coul(alph[0]) + cs[1]*F_LJ(alph[1:3]) + cs[2]*F_stretch(alph[3:5]) + cs[3]*F_bend(alph[5]) '''
class vpres_ffm( object ):  
    def __init__( self, natom, ndata, linp0, linp_bounds, traj_fname, aifrc_fname, parout_fname):
        ndim = 3
        self.ndim = ndim
        self.natom = natom
        self.ndata = ndata
        self.narr = ndim*natom*ndata
        self.counteval = 0
        
        self.linpars = np.empty(4)                  #linear parameters
        self.nolpars = np.empty(6)                  #nonlinear parameters
        self.funstore = np.zeros(self.narr)  #stores the total function
        self.phistore = np.zeros((self.narr,4))#stores the basis functions
        
        self.linpars=linp0
        self.linpar_bounds=linp_bounds
        self.parout_fname=parout_fname
        self.traj_fname=traj_fname

        self.aifrc=read_xyzdat(aifrc_fname,self.narr)
        
        stdoutfile = file("stdout", 'w')
        stdoutfile.write("")
        stdoutfile.write("-------        FFM Script Output        -------")
        stdoutfile.write("")
        stdoutfile.close()

    def __call__( self, params, respower):
        '''params is expected to be of the Parameters type (lmfit)'''
        hparams=np.empty(6)
        hparams[0]=params["alpha"].value
        hparams[1]=params["oo_sig"].value
        hparams[2]=params["oo_gam"].value
        hparams[3]=params["alp"].value
        hparams[4]=params["reoh"].value
        hparams[5]=params["thetad"].value
        if not all(hparams[0:6]==self.nolpars[0:6]):
            self.nolpars=hparams
            self.fitfun_ffm_phi()
            
            #self.linpars=bounded_lsq(self.phistore,self.aifrc,self.linpar_bounds[:,0],self.linpar_bounds[:,1])
            
            self.linpars=spla.lstsq(self.phistore,self.aifrc)
            self.linpars=self.linpars[0]
            
            self.fitfun_ffm()
            
        return np.power(np.absolute(self.aifrc-self.funstore),respower)
   
    def fitfun_ffm_phi(self):
        '''Expects cs[0:4],alph[0:6],R[ndata,0:3],ndata,na'''
    
        outfile = file(self.parout_fname, 'w')
        outfile.write('&param')
        outfile.write('\nwmass = %.16f'%    32831.2525000000000003e0)
        outfile.write('\nomass = %.16f'%    29156.9471000000000001e0)
        outfile.write('\nhmass = %.16f'%    1837.1527000000000001e0)
        outfile.write('\nqo = %.16f'%       - np.sqrt(self.linpars[0]))
        outfile.write('\nalpha = %.16f'%    self.nolpars[0])
        outfile.write('\noo_sig = %.16f'%   self.nolpars[1])
        outfile.write('\noo_eps = %.16f'%   self.linpars[1]) 
        outfile.write('\noo_gam = %.16f'%   self.nolpars[2])
        outfile.write('\nthetad = %.16f'%   self.nolpars[5])
        outfile.write('\nreoh = %.16f'%     self.nolpars[4])
        outfile.write('\napot = %.16f'%     self.linpars[2])
        outfile.write('\nbpot = %.16f'%     self.linpars[3])
        outfile.write('\nalp = %.16f'%      self.nolpars[3])
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
        
        stdoutfile = file("stdout", 'a')
        stdoutfile.write("---    Output of Evaluation No. " + str(self.counteval) + "    ---")
        stdoutfile.write(output)
        stdoutfile.write("")
        stdoutfile.close()
        
        self.phistore[:,0]=read_xyzdat("vmd.frc1",self.narr)/self.linpars[0]
        self.phistore[:,1]=read_xyzdat("vmd.frc2",self.narr)/self.linpars[1]
        self.phistore[:,2]=read_xyzdat("vmd.frc3",self.narr)/self.linpars[2]
        self.phistore[:,3]=read_xyzdat("vmd.frc4",self.narr)/self.linpars[3]
        
    def fitfun_ffm(self):
        ''' Call fitfun_ffm_phi before calling this!'''
        #self.fitfun_ffm_phi()
        
        self.funstore[:]=0.0
        for it in range(4):
            self.funstore[:] += self.linpars[it]*self.phistore[:,it]
        
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
print dataf
print np.size(dataf)
noisefac=0.05*meanfabs
for it1 in range(ndata*6):
    dataf[it1]+=noisefac*rd.random()
print dataf
    
testvp=vpres_test(2,ndata)    
residual=testvp(exnolpars,dataf,Rtest)
print residual    

fitpars=spo.leastsq(testvp,[-0.1,1.0],args=(dataf, Rtest))
print fitpars
'''

ndata=1500
linpars0=np.array([1.21,0.00025,0.135,0.062])
linp_bounds=np.array([[0.5,1.5],[0.0001,0.001],[0.05,0.2],[0.01,0.1]])
nolpars0=np.array([0.7,5.8,13.0,1.3,1.84,107.0])
fitparams = lmf.Parameters()
fitparams.add('alpha',     value = nolpars0[0],  min = 0.500,  max = 1.000)
fitparams.add('oo_sig',    value = nolpars0[1],  min = 0.010,  max = 10.000)
fitparams.add('oo_gam',    value = nolpars0[2],  min = 0.010,  max = 50.000)
fitparams.add('alp',       value = nolpars0[3],  min = 0.010,  max = 100.0)
fitparams.add('reoh',      value = nolpars0[4],  min = 1.500,  max = 2.000)
fitparams.add('thetad',    value = nolpars0[5],  min = 90.00,   max = 120.0)

natom=375
ndata=1500
linp0=linpars0
traj_fname="all_tray.xyz"
aifrc_fname="FORCES-PBE-ALL.frc"
parout_fname="param_PBE-new"
respower=1.0
fitftol=1.e-12
lskwords={'ftol':fitftol}

vpres_ffm_obj=vpres_ffm(natom, ndata, linpars0, linp_bounds, traj_fname, aifrc_fname, parout_fname)
residual=vpres_ffm_obj(fitparams,respower)

fitout=lmf.minimize(vpres_ffm_obj,fitparams,args=(respower,),ftol=fitftol)
print vpres_ffm_obj.linpars
print vpres_ffm_obj.nolpars

lmf.report_fit(fitparams)