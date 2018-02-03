import time
import os
import numpy as np
import matplotlib.pyplot as plt
import h5py
plt.rcParams["font.family"] = "Times New Roman"
class exp:
    ###################################################################
    # Global functions
    def __init__(self,
                            exp_dir= './',
                            R11_dir= './',
                            direction='u',
                            manufactured=False,
                            xi=0,
                            yi=0,
                            method=1,
                            option=2,
                            t=0,
                            X=0,
                            Y=0,
                            Z=0,
                            U=0,
                            homogeneous=(0,1),
                            average_z=True,
                            skip_last=True,
                            moving_mean=False,
                            nwindow=5,
                            labelfont={'fontsize':14,'fontname':'Times New Roman'},
                            D_ref=0.00745,
                            U_ref=27.5):
        ''' calc_options needs to be dictionary containing values defined by set_values function'''
        # calc_options
        self.exp_dir = exp_dir
        self.R11_dir = R11_dir
        self.direction = direction
        self.manufactured = manufactured
        self.xi = xi
        self.yi = yi
        self.method = method
        self.option = option
        self.t,self.X,self.Y,self.Z,self.U = t,X,Y,Z,U
        self.homogeneous=homogeneous
        self.average_z=average_z
        self.skip_last=skip_last
        if self.average_z:
            self.homogeneous=(0,1) # time and z direction
        else:
            self.homogeneous=(0,)
        self.moving_mean=moving_mean
        self.nwindow=nwindow
        self.labelfont=labelfont
        # non-dimensionalize stuff
        self.D_ref=D_ref
        self.U_ref=U_ref

        # load data
        self.f1 = h5py.File(self.exp_dir+'T093013_01_Velocity_1-500.mat')
        self.f2 = h5py.File(self.exp_dir+'T093013_01_Velocity_501-1000.mat')
        self.f3 = h5py.File(self.exp_dir+'T093013_01_Velocity_1001-1500.mat')
        self.f4 = h5py.File(self.exp_dir+'T093013_01_Velocity_1501-2000.mat')
        self.f5 = h5py.File(self.exp_dir+'T093013_01_Velocity_2001-2500.mat')
        self.f6 = h5py.File(self.exp_dir+'T093013_01_Velocity_2501-3000.mat')
        self.f7 = h5py.File(self.exp_dir+'T093013_01_Velocity_3001-3500.mat')
        self.f8 = h5py.File(self.exp_dir+'T093013_01_Velocity_3501-4000.mat')
        self.f9 = h5py.File(self.exp_dir+'T093013_01_Velocity_4001-4500.mat')
        self.f10 = h5py.File(self.exp_dir+'T093013_01_Velocity_4501-5000.mat')
        self.f_all=[self.f1,self.f2,self.f3,self.f4,self.f5,self.f6,self.f7,self.f8,self.f9,self.f10]

    def read_data_time(self,t=0,h5py_file=None):
        if h5py_file==None:
            h5py_file=self.f1
        if self.direction=='u':
            self.Us=np.array(h5py_file['Us'][t])
        elif self.direction=='v':
            self.Us=np.array(h5py_file['Vs'][t])
        elif self.direction=='w':
            self.Us=np.array(h5py_file['Ws'][t])
        self.X=np.array(h5py_file['X'][0,])
        self.Y=np.array(h5py_file['Y'][0,])
        self.Z=np.array(h5py_file['Z'][0,])
        # recenter and change to m
        self.Y=self.Y/(1000.)+15.*self.D_ref # mm to m and 15 diam. downstream
        self.X=self.X/1000.-(1.15*self.D_ref) # centered between nozzles
        self.Z=self.Z/1000.
        # non-dimensionalize everything
        self.X=self.X/self.D_ref
        self.Y=self.Y/self.D_ref
        self.Z=self.Z/self.D_ref
        self.Us=self.Us/self.U_ref
        self.Zmesh,self.Xmesh,self.Ymesh=np.meshgrid(self.Z,self.X,self.Y,indexing='ij')
    def read_data_U(self):
        self.Us=np.array(self.f1['Us'])
        self.Us=self.Us/self.U_ref
    def read_data_V(self):
        self.Us=np.array(self.f1['Vs'])
        self.Us=self.Us/self.U_ref
    def read_data_W(self):
        self.Us=np.array(self.f1['Ws'])
        self.Us=self.Us/self.U_ref
    def plot_time_show(self,time=0,z=3):
        self.read_data_time(t=time)
        plt.figure()
        plt.colorbar(plt.contourf(self.X[z],self.Y[z],V[z],30,cmap='jet'))
        plt.title('time=%.4d'%time+' z=%.4d'%z)
        plt.xlabel('x',fontsize=24)
        plt.ylabel('y',fontsize=24)
        plt.show()

    ###################################################################
    ###################################################################
    # for length scales
    def read_data_r(self):
        if self.average_z:
            if self.direction=='u':
                self.Us=np.concatenate([np.array(f['Us'][:,:,self.xi,:]) for f in self.f_all],axis=0)
            elif self.direction=='v':
                self.Us=np.concatenate([np.array(f['Vs'][:,:,self.xi,:]) for f in self.f_all],axis=0)
            elif self.direction=='w':
                self.Us=np.concatenate([np.array(f['Ws'][:,:,self.xi,:]) for f in self.f_all],axis=0)
            else:
                print 'direction is wrong at ',self.direction
        else:
            if self.direction=='u':
                self.Us=np.concatenate([np.array(f['Us'][:,10,self.xi,:]) for f in self.f_all],axis=0)
            elif self.direction=='v':
                self.Us=np.concatenate([np.array(f['Vs'][:,10,self.xi,:]) for f in self.f_all],axis=0)
            elif self.direction=='w':
                self.Us=np.concatenate([np.array(f['Ws'][:,10,self.xi,:]) for f in self.f_all],axis=0)
            else:
                print 'direction is wrong at ',self.direction
        print '  experiment reading data in ',self.exp_dir
        if self.skip_last:
            self.Us=self.Us[:-1]
        self.t=np.arange(0,0.5,0.0001)
        self.X=np.array(f['X'][0,self.xi])
        self.Y=np.array(f['Y'][0,:])
        self.Z=np.array(f['Z'][0,:])
        # switch Y dimension ordering
        self.Us=self.Us[:,:,::-1]
        self.Y=self.Y[::-1]
        # recenter and change to m
        #D = 0.00745 # diameter of jet in m
        self.Y=self.Y/(1000.)+15.*self.D_ref # mm to m and 15 diam. downstream
        self.X=self.X/1000.-(1.15*self.D_ref) # centered between nozzles
        self.Z=self.Z/1000.

        # non-dimensionalize everything
        self.X=self.X/self.D_ref
        self.Y=self.Y/self.D_ref
        self.Z=self.Z/self.D_ref
        self.Us=self.Us/self.U_ref
        self.t=self.t/(self.D_ref/self.U_ref)
        if self.skip_last:
            self.t=self.t[:-1]
    def read_data_r_XYZ(self):
        self.X=np.array(self.f1['X'][0,self.xi])
        self.Y=np.array(self.f1['Y'][0,:])
        self.Z=np.array(self.f1['Z'][0,:])
        #

        self.Y=self.Y[::-1]
        # recenter and change to m
        self.Y=self.Y/(1000.)+15.*self.D_ref # mm to m and 15 diam. downstream
        self.X=self.X/1000.-(1.15*self.D_ref) # centered between nozzles
        self.Z=self.Z/1000.
        # non-dimensionalize everything
        self.X=self.X/self.D_ref
        self.Y=self.Y/self.D_ref
        self.Z=self.Z/self.D_ref
    def read_data_r_R11(self):
        if self.direction=='u':
            R11title=self.R11_dir+'R11_output/XYR11/ijXYR11_%i_%i_%.4f_%.4f'%(self.xi,self.yi,self.X,self.Y[self.yi])+'.npy'
        elif self.direction=='v':
            R11title=self.R11_dir+'R22_output/XYR22/ijXYR22_%i_%i_%.4f_%.4f'%(self.xi,self.yi,self.X,self.Y[self.yi])+'.npy'
        else:
            print 'something is wrong with your directions here'
        self.R11_norm=np.load(R11title)
    def set_uprime_r(self):
        # calc u' v' and w'  velocity prime values
        if self.moving_mean:
            print 'error, moving mean not implemented yet'
            k=np.arange(self.nwindow,len(t)-self.nwindow-1,1)
            umean=np.empty(k.shape)
            for ki in k:
                umean[ki-self.nwindow]=self.Us[ki-self.nwindow:ki+self.nwindow+1].mean(axis=self.homogeneous)
            if self.average_z:
                self.u=self.Us[self.nwindow:-1*self.nwindow-1]-umean[:,np.newaxis]
            else:
                self.u=self.Us[self.nwindow:-1*self.nwindow-1]-umean
            
        else:
            u_avg=self.Us.mean(axis=self.homogeneous)[np.newaxis,np.newaxis,:]
            self.u=self.Us-u_avg
    def calc_R11_r_exp(self):
        if self.method==0:
            nr=len(self.Y) # space axial domain
            self.R11=np.zeros(nr)
            r=np.arange(0,nr)-self.yi
            for ri in r:
                self.R11[ri+self.yi]=(self.u[:,:,self.yi]*self.u[:,:,self.yi+ri]).mean(axis=self.homogeneous)
            self.r=self.Y[:nr]-self.Y[self.yi]
        elif self.method==1:
            self.yi=0
            nr=len(self.Y) # space axial domain
            self.R11=np.zeros(nr)
            #print 'error method 1 not implemented yet'
            #self.R11[0]= (self.u[:,:,0]*self.u[:,:,0]).mean(axis=self.homogeneous)
            r=np.arange(0,nr)
            for ri in r: 
                self.R11[ri]=((self.u[:,:,0]*self.u[:,:,ri]).mean(axis=self.homogeneous))
            self.r=self.Y[:nr]-self.Y[0]
        elif self.method==2:
            print 'error method 2 not implemented yet'
            r=np.arange(0,nr)
            for rispace in r:
                for rime in r:
                    self.R11[rispace]=self.R11[rispace]+(self.u[rime]*self.u[(rime+rispace)]).mean()
            self.R11=self.R11/float(nr)
        self.R11_norm= self.R11/self.R11[self.yi]
    def calc_integral_length_L(self):
        # calc tau integral scale
        self.L=np.empty(self.R11_norm.shape)
        #for i in range(len(self.L)):
            #self.L[i]=np.trapz(self.R11_norm[:,i],x=self.r)
        self.L=np.trapz(self.R11_norm,x=self.r)
        print 'L= ',self.L
    def calc_taylor_length_lambda(self,ax=None):
        r=self.r
        R11=self.R11_norm
        # set edges to mirror
        #  n_edge=10
        # edge_array=np.arange(1,n_edge)
        # dr=np.mean(np.diff(self.r))
        # rprevious = (self.r[0]-dr*edge_array)[::-1]
        # rafter = (self.r[-1]+dr*edge_array)
        # r = np.concatenate([rprevious,self.r,rafter])
        # print 'r',r
        # if self.yi>=len(self.r)/2:
        #     negative_yi=self.yi-len(self.r)
        #     R11previous = self.R11_norm[:n_edge] # bogus data to fill first few spots
        #     R11after = self.R11_norm[(-1*edge_array+negative_yi*2+1)]
        # else:
        #     R11previous = self.R11_norm[((edge_array+self.yi*2)[::-1])]
        #     R11after = self.R11_norm[-1:-1*n_edge:-1] # bogus data to fill last few spots

        # #R11 = np.concatenate([self.R11_norm[:n_edge][::-1],self.R11_norm,self.R11_norm[-1*n_edge:]])
        # R11 = np.concatenate([R11previous,self.R11_norm,R11after])
        self.set_plot_label_r()
        #rnew=np.linspace(self.r.min(),self.r.max(),300)
        #interp=interp1d(r,R11,kind='quadratic')
        #ax.plot(rnew,interp(rnew),'sg')
        dR=np.gradient(R11,edge_order=self.option)
        dr=np.gradient(r,edge_order=self.option)
        #dRdr=dR/dr
        d2R=np.gradient(dR,edge_order=self.option)
        #d2r=np.gradient(dr,edge_order=self.option)
        #d2Rdr2=d2R/dr**2
        d2Rdr2=d2R/dr**2
        f_prime_prime=d2Rdr2[self.yi]
        #lambdaf=(-0.5*d2Rdr2[self.yi])**-0.5+self.r[self.yi]
        #ax.plot(r,1-(r-r[self.yi])**2/lambdaf**2,'s',label='taylor '+self.label)
        #P=1-(r-r[self.yi])**2/lambdaf**2
        self.lambdaf=(-f_prime_prime/2.)**-0.5+r[self.yi] # r[self.yi] should equal 0
        #print 'lambdaf= ',self.lambdaf
        #print self.lambdaf,r[self.yi],f_prime_prime
        if ax!=None:
            P=1.+0.5*f_prime_prime*(r-r[self.yi])**2
            ax.plot(r,P,'s',label=self.label+' taylor poly')
            ax.plot(self.lambdaf,0.,'o',label=self.label+r' $\lambda_f$')
        print 'TaylorL=',self.lambdaf


    def plot_R11_r_exp(self,ax,plot_options,label=''):
        self.set_plot_label_r()
        ax.plot(self.r,self.R11_norm,label=self.label+','+label,**plot_options)
    def set_axis_labels_r(self):
        ax.set_xlabel(r'$r/D$',fontsize=24)
        if self.direction=='u':
            ax.set_ylabel(r'$R_{11}(r)/R_{11}(0)$',**self.labelfont)
        elif self.direction=='v':
            ax.set_ylabel(r'$R_{22}(r)/R_{22}(0)$',**self.labelfont)
        elif self.direction=='w':
            ax.set_ylabel(r'$R_{33}(r)/R_{33}(0)$',**self.labelfont)
        plt.tight_layout()
        ax.legend(loc='best',numpoints=1)
    def set_plot_label_r(self):
        self.label='exp,%s,x=%.4f,y=%.4f,m=%i,o=%i'%(
                self.direction,
                self.X ,
                self.Y[self.yi] ,
                self.method ,
                self.option ,)
    def save_R11_r(self,f_handle,):
        np.savetxt(f_handle,[[self.X,self.Y[self.yi],self.L,self.lambdaf]],delimiter=' ')
        if self.direction=='u':
            np.save('R11_output/XYR11/ijXYR11_%i_%i_%.4f_%.4f'%(self.xi,self.yi,self.X,self.Y[self.yi]),self.R11_norm)
        elif self.direction=='v':
            np.save('R22_output/XYR22/ijXYR22_%i_%i_%.4f_%.4f'%(self.xi,self.yi,self.X,self.Y[self.yi]),self.R11_norm)
        elif self.direction=='w':
            np.save('R33_output/XYR33/ijXYR33_%i_%i_%.4f_%.4f'%(self.xi,self.yi,self.X,self.Y[self.yi]),self.R11_norm)
        print 'saved ijxy=',self.xi,self.yi,self.X,self.Y[self.yi]
    def Manufactured_R11_r(self):
        # manufactured data (to check accuracy of method)
        L=self.r[-1]-self.r[0]
        # option 1
        if self.option==1:
            self.R11=(
                1./(L*np.exp(self.r)) 
                * ( 
                    (np.cos(self.r)/5. 
                    * (
                    np.exp(-1.*L) * (-1.*np.cos(L)**2 + 2.*np.sin(L)*np.cos(L)-2.)
                    + 3.))
                - (np.sin(self.r)/2.  * ( np.sin(L)**2))
                )
                ) 
        # option 2
        elif self.option==2:
            self.R11=(1./L)*(
                    np.cos(self.r)
                    * (( np.sin(L)      *np.cos(L)+L )  /   2.)
                    - ((  np.sin(self.r)   *np.sin(L)**2)  /   2.))
        #print 'manufactured R11        =',R11
        R11_norm=self.R11/self.R11[0]
        self.R11expMNF=R11_norm
        #print 'manufactured R11/R11[0] =',R11_norm
        # return normalized R11
        return R11_norm
    def set_manufactured_r(self):
        # manufactured data (to check accuracy of method)
        self.u=np.empty((10,6,len(self.r)))
        if self.manufactured==True:
            if self.option==1:
                self.u=(np.cos(self.r)*np.exp(-1.*self.r))[np.newaxis,np.newaxis,:] + 0.08*(2.*np.random.random(self.u.shape)-1.0)
            elif self.option==2:
                self.u=np.cos(self.r)[np.newaxis,np.newaxis,:] + 0.08*(2.*np.random.random(self.u.shape)-1.0)

    def R11_r_calc_plot(self,ax,plot_options,label=''):
        # Experimental data set
        if self.manufactured:
            #self.r=np.linspace(4,36,1260)
            self.r=np.linspace(0,0.005/(self.D_ref/self.U_ref),5000.)
            self.Y=self.r
            self.set_manufactured_r()
            self.Manufactured_R11_r()
            ax.plot(self.r,self.R11expMNF,marker='x',linewidth=0,label='exp Manufactured')
        else:
            self.read_data_r()
            self.set_uprime_r()
        self.calc_R11_r_exp()
        self.plot_R11_r_exp(ax,plot_options,label=label)
        self.calc_integral_length_L()
        self.calc_taylor_length_lambda(ax)
    def R11_r_calc_save(self,f_handle):
        # Experimental data set
        if self.manufactured:
            #self.r=np.linspace(4,36,1260)
            self.r=np.linspace(0,0.005/(self.D_ref/self.U_ref),5000.)
            self.Y=self.r
            self.set_manufactured_r()
            self.Manufactured_R11_r()
            ax.plot(self.r,self.R11expMNF,marker='x',linewidth=0,label='exp Manufactured')
        else:
            self.read_data_r()
            self.set_uprime_r()
        self.calc_R11_r_exp()
        #self.plot_R11_r_exp(ax,plot_options,label=label)
        self.calc_integral_length_L()
        self.calc_taylor_length_lambda()
        self.save_R11_r(f_handle)
    def R11_r_loadR11_calc_plot(self,ax,plot_options,label=''):
        # Experimental data set and precalculated R11 values
        self.read_data_r_XYZ()
        self.read_data_r_R11()
        self.plot_R11_r_exp(ax,plot_options,label=label)
        self.calc_integral_length_L()
        self.calc_taylor_length_lambda(ax)
    def R11_r_loadR11_calc_save(self,f_handle):
        self.read_data_r_XYZ()
        self.read_data_r_R11()
        self.calc_integral_length_L()
        self.calc_taylor_length_lambda()
        np.savetxt(f_handle,[[self.X,self.Y[self.yi],self.L,self.lambdaf]],delimiter=' ')
        #self.save_R11_r(f_handle)
    ###################################################################
    ###################################################################
    # for time scales
    def read_data_pt(self):
        if self.average_z:
            if self.direction=='u':
                self.Us=np.concatenate([np.array(f['Us'][:,:,self.xi,-1-self.yi]) for f in self.f_all],axis=0)
            elif self.direction=='v':
                self.Us=np.concatenate([np.array(f['Vs'][:,:,self.xi,-1-self.yi]) for f in self.f_all],axis=0)
            elif self.direction=='w':
                self.Us=np.concatenate([np.array(f['Ws'][:,:,self.xi,-1-self.yi]) for f in self.f_all],axis=0)
            else:
                print 'direction is wrong at ',self.direction
        else:
            if self.direction=='u':
                self.Us=np.concatenate([np.array(f['Us'][:,10,self.xi,-1-self.yi]) for f in self.f_all],axis=0)
            elif self.direction=='v':
                self.Us=np.concatenate([np.array(f['Vs'][:,10,self.xi,-1-self.yi]) for f in self.f_all],axis=0)
            elif self.direction=='w':
                self.Us=np.concatenate([np.array(f['Ws'][:,10,self.xi,-1-self.yi]) for f in self.f_all],axis=0)
            else:
                print 'direction is wrong at ',self.direction
        print '  experiment reading data in ',self.exp_dir
        if self.skip_last:
            self.Us=self.Us[:-1]
        self.t=np.arange(0,0.5,0.0001)
        self.X=np.array(f['X'][0,self.xi])
        self.Y=np.array(f['Y'][0,-1-self.yi])
        self.Z=np.array(f['Z'][0,:])
        # recenter and change to m
        #D = 0.00745 # diameter of jet in m
        self.Y=self.Y/(1000.)+15.*self.D_ref # mm to m and 15 diam. downstream
        self.X=self.X/1000.-(1.15*self.D_ref) # centered between nozzles
        self.Z=self.Z/1000.

        # non-dimensionalize everything
        self.X=self.X/self.D_ref
        self.Y=self.Y/self.D_ref
        self.Z=self.Z/self.D_ref
        self.Us=self.Us/self.U_ref
        self.t=self.t/(self.D_ref/self.U_ref)
        if self.skip_last:
            self.t=self.t[:-1]
    def set_uprime(self):
        # calc u' v' and w'  velocity prime values
        if self.moving_mean:
            k=np.arange(self.nwindow,len(t)-self.nwindow-1,1)
            umean=np.empty(k.shape)
            for ki in k:
                umean[ki-self.nwindow]=self.Us[ki-self.nwindow:ki+self.nwindow+1].mean(axis=self.homogeneous)
            if self.average_z:
                self.u=self.Us[self.nwindow:-1*self.nwindow-1]-umean[:,np.newaxis]
            else:
                self.u=self.Us[self.nwindow:-1*self.nwindow-1]-umean
            
        else:
            u_avg=self.Us.mean(axis=self.homogeneous)
            self.u=self.Us-u_avg
    def calc_R11_tau_exp(self):
        nt=len(self.t)/2 # half of time domain
        self.R11=np.zeros(nt)
        if self.method==0:
            tau=np.arange(0,nt)
            for taui in tau:
                self.R11[taui]=(self.u[:nt]*self.u[taui:nt+taui]).mean()
        elif self.method==1:
            self.R11[0]= (self.u*self.u).mean()
            tau=np.arange(1,nt)
            for taui in tau: 
                self.R11[taui]=((self.u[:(-1*taui)]*self.u[taui:]).mean())
        elif self.method==2:
            tau=np.arange(0,nt)
            for tauispace in tau:
                for tauime in tau:
                    self.R11[tauispace]=self.R11[tauispace]+(self.u[tauime]*self.u[(tauime+tauispace)]).mean()
            self.R11=self.R11/float(nt)
        self.tau=self.t[:len(self.t)/2]
        self.R11_norm= self.R11/self.R11[0]
    def calc_integral_time_T(self):
        # calc tau integral scale
        self.T=np.trapz(self.R11_norm,x=self.tau)
        print 'T= ',self.T
    def plot_R11_tau_exp(self,ax,plot_options,label=''):
        self.set_plot_label()
        ax.plot(self.tau,self.R11_norm,label=self.label+','+label,**plot_options)
    def set_axis_labels_tau(self):
        ax.set_xlabel(r'$\tau/\tau_\mathrm{ref}$',**self.labelfont)
        if self.direction=='u':
            ax.set_ylabel(r'$R_{11}(\tau)$',**self.labelfont)
        elif self.direction=='v':
            ax.set_ylabel(r'$R_{22}(\tau)$',**self.labelfont)
        elif self.direction=='w':
            ax.set_ylabel(r'$R_{33}(\tau)$',**self.labelfont)
        #ax.set_xlabel(r'$\tau/\tau_\mathrm{ref}$',**self.labelfont)
        #ax.set_ylabel(r'$R_{11}(\tau)$',**self.labelfont)
        plt.tight_layout()
        ax.legend(loc='best',numpoints=1)
    def set_plot_label(self):
            self.label='exp,%s,x=%.4f,y=%.4f,m=%i,o=%i'%(
            self.direction,
            self.X ,
            self.Y ,
            self.method ,
            self.option ,)
    def save_R11_tau(self,f_handle,):
        self.calc_R11_tau_exp(**calc_options)
        self.tau=self.t[:len(self.t)/2]
        # normalize by variance
        self.R11_norm= self.R11/self.R11[0]
        # calc tau integral scale
        self.T=np.trapz(R11_norm,x=tau)
        np.savetxt(f_handle,[[self.X,self.Y,self.T]],delimiter=' ')
        if self.direction=='u':
            np.save('R11_output/XYR11/ijXYR11_'+str(self.xi)+'_'+str(self.yi)+'_%.4f_%.4f'%(self.X,self.Y),R11_norm)
        elif self.direction=='v':
            np.save('R22_output/XYR22/ijXYR22_'+str(self.xi)+'_'+str(self.yi)+'_%.4f_%.4f'%(self.X,self.Y),R11_norm)
        elif self.direction=='w':
            np.save('R33_output/XYR33/ijXYR33_'+str(self.xi)+'_'+str(self.yi)+'_%.4f_%.4f'%(self.X,self.Y),R11_norm)
        print 'saved ijxy=',self.xi,self.yi,self.X,self.Y
    def Manufactured_R11(self):
        # manufactured data (to check accuracy of method)
        self.tau=self.t[:len(self.t)/2]
        L=self.t[-1]-self.t[0]
        # option 1
        if self.option==1:
            self.R11=(
                1./(L*np.exp(self.tau)) 
                * ( 
                    (np.cos(self.tau)/5. 
                    * (
                    np.exp(-1.*L) * (-1.*np.cos(L)**2 + 2.*np.sin(L)*np.cos(L)-2.)
                    + 3.))
                - (np.sin(self.tau)/2.  * ( np.sin(L)**2))
                )
                ) 
        # option 2
        elif self.option==2:
            self.R11=(1./L)*(
                    np.cos(self.tau)
                    * (( np.sin(L)      *np.cos(L)+L )  /   2.)
                    - ((  np.sin(self.tau)   *np.sin(L)**2)  /   2.))
        #print 'manufactured R11        =',R11
        R11_norm=self.R11/self.R11[0]
        self.R11expMNF=R11_norm
        #print 'manufactured R11/R11[0] =',R11_norm
        # return normalized R11
        return R11_norm
    def set_manufactured_tau(self):
        # manufactured data (to check accuracy of method)
        self.u=np.empty((len(self.t),)+(20,))
        if self.manufactured==True:
            if self.option==1:
                self.u=(np.cos(self.t)*np.exp(-1.*self.t)).reshape((len(self.t),1)) + 0.08*(2.*np.random.random((len(self.t),)+self.u.shape[1:])-1.0)
            elif self.option==2:
                self.u=np.cos(self.t).reshape((len(self.t),1)) + 0.08*(2.*np.random.random((len(self.t),)+self.u.shape[1:])-1.0)

    def R11_tau_calc_plot(self,ax,plot_options,label=''):
        # Experimental data set
        if self.manufactured:
            self.t=np.linspace(0,0.05/(self.D_ref/self.U_ref),5000.)
            self.set_manufactured_tau()
            self.Manufactured_R11()
            ax.plot(self.tau,self.R11expMNF,marker='x',linewidth=0,label='exp Manufactured')
        else:
            self.read_data_pt()
            self.set_uprime()
        self.calc_R11_tau_exp()
        self.plot_R11_tau_exp(ax,plot_options,label=label)
        self.calc_integral_time_T()
        #plt.show()
        #return R11exp1,R11exp2,R11expMNF
        #return R11exp1,_tau,R11expMNF,T
    ###################################################################
if __name__=="__main__":
    import sys
    direction=sys.argv[1]
    print 'direction=',direction
    # user inputs
    calc_options1={
            'exp_dir':'../DA_ExperimentalData/From_JFrank_5000_cold_SandiaC/',
            'R11_dir':'/home/shaun/Desktop/DA/DA_ExperimentalData/From_JFrank_5000_cold_SandiaC/R11_output/LengthScales_Taylor4/R11_lengthScales_taylor4/',
            'direction':direction,
            'manufactured':False,
            'method':0,
            'option':2,             # order of taylor derivative difference scheme (option 2 means second order accurate)
            'skip_last':True,
            'moving_mean':False,
            'average_z':True,
            'xi':-3,
            'yi':48,
            }
    plot_options1={
            'linewidth':0,
            'marker':'+',
            }
    label=''
    exp1=exp(**calc_options1)
    ###############################
    # for length scales
    ###### calc and plot
    fig=plt.figure()
    ax=plt.subplot(111)
    exp1.xi=47
    exp1.yi=0
    #exp1.R11_r_calc_plot(ax,plot_options1,label=label)
    #exp1.R11_r_loadR11_calc_plot(ax,plot_options1,label='loaded')
    exp1.yi=1
    exp1.R11_r_calc_plot(ax,plot_options1,label=label)
    plot_options1['marker']='x'
    exp1.R11_r_loadR11_calc_plot(ax,plot_options1,label='loaded')
    exp1.yi=2
    #exp1.R11_r_calc_plot(ax,plot_options1,label=label)
    exp1.yi=3
    #exp1.R11_r_calc_plot(ax,plot_options1,label=label)
    exp1.set_axis_labels_r()
    plt.show()

    ###### calc and save
    #  if direction=='u':
    #      with open('R11_output/R11.txt','a') as f_handle:
    #          for xi in range(169):
    #              for yi in range(126):
    #                  exp1.xi=xi
    #                  exp1.yi=yi
    #                  exp1.R11_r_calc_save(f_handle)
    #  elif direction=='v':
    #      with open('R22_output/R22.txt','a') as f_handle:
    #          for xi in range(169):
    #              for yi in range(126):
    #                  exp1.xi=xi
    #                  exp1.yi=yi
    #                  exp1.R11_r_calc_save(f_handle)
    #  elif direction=='w':
    #      with open('R33_output/R33.txt','a') as f_handle:
    #          for xi in range(169):
    #              for yi in range(126):
    #                  exp1.xi=xi
    #                  exp1.yi=yi
    #                  exp1.R11_r_calc_save(f_handle)

    ###############################
    # for time scales
    #fig=plt.figure()
    #ax=plt.subplot(111)
    #exp1.R11_tau_calc_plot(ax,plot_options1,label=label)
    #exp1.xi=-63
    #exp1.R11_tau_calc_plot(ax,plot_options1,label=label)
    #exp1.xi=-122
    #exp1.R11_tau_calc_plot(ax,plot_options1,label=label)
    #exp1.set_axis_labels_tau()
    #plt.show()
