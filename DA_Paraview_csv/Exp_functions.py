import time
import os
import numpy as np
import matplotlib.pyplot as plt
import h5py
plt.rcParams["font.family"] = "Times New Roman"
# where experimental data is stored (subdirectory or total location)
exp_dir='../DA_ExperimentalData/From_JFrank_5000_cold_SandiaC/'

f1 = h5py.File(exp_dir+'T093013_01_Velocity_1-500.mat')
f2 = h5py.File(exp_dir+'T093013_01_Velocity_501-1000.mat')
f3 = h5py.File(exp_dir+'T093013_01_Velocity_1001-1500.mat')
f4 = h5py.File(exp_dir+'T093013_01_Velocity_1501-2000.mat')
f5 = h5py.File(exp_dir+'T093013_01_Velocity_2001-2500.mat')
f6 = h5py.File(exp_dir+'T093013_01_Velocity_2501-3000.mat')
f7 = h5py.File(exp_dir+'T093013_01_Velocity_3001-3500.mat')
f8 = h5py.File(exp_dir+'T093013_01_Velocity_3501-4000.mat')
f9 = h5py.File(exp_dir+'T093013_01_Velocity_4001-4500.mat')
f10 = h5py.File(exp_dir+'T093013_01_Velocity_4501-5000.mat')
f_all=[f1,f2,f3,f4,f5,f6,f7,f8,f9,f10]

def read_data_pt(i=0,j=0,direction='u',skip_last=False,average_z=True):
    if average_z:
        if direction=='u':
            Us=np.concatenate([np.array(f['Us'][:,:,i,-1-j]) for f in f_all],axis=0)
        elif direction=='v':
            Us=np.concatenate([np.array(f['Vs'][:,:,i,-1-j]) for f in f_all],axis=0)
        elif direction=='w':
            Us=np.concatenate([np.array(f['Ws'][:,:,i,-1-j]) for f in f_all],axis=0)
        else:
            print 'direction is wrong at ',direction
    else:
        if direction=='u':
            Us=np.concatenate([np.array(f['Us'][:,10,i,-1-j]) for f in f_all],axis=0)
        elif direction=='v':
            Us=np.concatenate([np.array(f['Vs'][:,10,i,-1-j]) for f in f_all],axis=0)
        elif direction=='w':
            Us=np.concatenate([np.array(f['Ws'][:,10,i,-1-j]) for f in f_all],axis=0)
        else:
            print 'direction is wrong at ',direction
    print '  experiment reading data in ',exp_dir
    if skip_last:
        Us=Us[:-1]
    t=np.arange(0,0.5,0.0001)
    X=np.array(f['X'][0,i])
    Y=np.array(f['Y'][0,-1-j])
    Z=np.array(f['Z'][0,:])
    # recenter and change to m
    D = 0.00745 # diameter of jet in m
    Y=Y/(1000.)+15.*D # mm to m and 15 diam. downstream
    X=X/1000.-(1.15*D) # centered between nozzles
    Z=Z/1000.

    # non-dimensionalize everything
    D_ref=0.00745
    U_ref=27.5
    X=X/D_ref
    Y=Y/D_ref
    Z=Z/D_ref
    Us=Us/U_ref
    t=t/(D_ref/U_ref)
    if skip_last:
        t=t[:-1]
    return [t,X,Y,Z,Us]
def read_data_time(t=0,h5py_file=f1,direction='u'):
    if direction=='u':
        Us=np.array(h5py_file['Us'][t])
    elif direction=='v':
        Us=np.array(h5py_file['Vs'][t])
    elif direction=='w':
        Us=np.array(h5py_file['Ws'][t])
    X=np.array(h5py_file['X'][0,])
    Y=np.array(h5py_file['Y'][0,])
    Z=np.array(h5py_file['Z'][0,])
    # recenter and change to m
    D = 0.00745 # diameter of jet in m
    Y=Y/(1000.)+15.*D # mm to m and 15 diam. downstream
    X=X/1000.-(1.15*D) # centered between nozzles
    Z=Z/1000.
    # non-dimensionalize everything
    D_ref=0.00745
    U_ref=27.5
    X=X/D_ref
    Y=Y/D_ref
    Z=Z/D_ref
    Us=Us/U_ref
    Zmesh,Xmesh,Ymesh=np.meshgrid(Z,X,Y,indexing='ij')
    return [Xmesh,Ymesh,Zmesh,Us]
def read_data_U():
    Us=np.array(f['Us'])
    U_ref=27.5
    Us=Us/U_ref
    return Us
def read_data_V():
    Vs=np.array(f['Vs'])
    U_ref=27.5
    Vs=Vs/U_ref
    return Vs
def read_data_W():
    Ws=np.array(f['Ws'])
    U_ref=27.5
    Ws=Ws/U_ref
    return Ws
def plot_time_show(time=0,z=3):
    X,Y,Z,U,V,W=read_data_time(t=time)
    plt.figure()
    plt.colorbar(plt.contourf(X[z],Y[z],V[z],30,cmap='jet'))
    plt.title('time=%.4d'%time+' z=%.4d'%z)
    plt.xlabel('x',fontsize=24)
    plt.ylabel('y',fontsize=24)
    plt.show()
def calc_R11_tau_exp(i=0,j=0,direction='u',method=0,manufactured=False,option=1,skip_last=False,moving_mean=False,average_z=True,nwindow=5):
    t,X,Y,Z,Us=read_data_pt(i=i,j=j,direction=direction,skip_last=skip_last,average_z=average_z)
    # calc u' v' and w'  velocity prime values
    if average_z:
        homogeneous=(0,1) # time and z direction
    else:
        homogeneous=(0,)
    if moving_mean:
        #nwindow=5
        k=np.arange(nwindow,len(t)-nwindow-1,1)
        #uwindow=np.empty(k.shape)
        umean=np.empty(k.shape)
        for ki in k:
            umean[ki-nwindow]=Us[ki-nwindow:ki+nwindow+1].mean(axis=homogeneous)
        if average_z:
            u=Us[nwindow:-1*nwindow-1]-umean[:,np.newaxis]
            #u=uwindow-Us.mean(axis=homogeneous)
        else:
            u=Us[nwindow:-1*nwindow-1]-umean
            #u=uwindow-Us.mean(axis=homogeneous)
        
    else:
        u_avg=Us.mean(axis=homogeneous)
        u=Us-u_avg
    
    # manufactured data (to check accuracy of method)
    if manufactured==True:
        # option 1
        if option==1:
            u=(np.cos(t)*np.exp(-1.*t)).reshape((len(t),1)) + 0.08*(2.*np.random.random((len(t),)+u.shape[1:])-1.0)
        # option 2
        elif option==2:
            u=np.cos(t).reshape((len(t),1)) + 0.08*(2.*np.random.random((len(t),)+u.shape[1:])-1.0)
    nt=len(t)/2 # half of time domain
    R11=np.zeros(nt)
    if method==0:
        tau=np.arange(0,nt)
        for taui in tau:
            R11[taui]=(u[:nt]*u[taui:nt+taui]).mean(axis=homogeneous)
    elif method==1:
        R11[0]= (u*u).mean()
        tau=np.arange(1,nt)
        for taui in tau: 
            R11[taui]=((u[:(-1*taui)]*u[taui:]).mean())
    elif method==2:
        tau=np.arange(0,nt)
        for tauispace in tau:
            for tauime in tau:
                R11[tauispace]=R11[tauispace]+(u[tauime]*u[(tauime+tauispace)]).mean()
        R11=R11/float(nt)
    return t,X,Y,Z,R11
def plot_R11_tau_exp(ax,i,j,calc_options,plot_options):
    t,X,Y,Z,R11=calc_R11_tau_exp(i=i,j=j,**calc_options)
    tau=t[:len(t)/2]
    # normalize by variance
    R11_norm= R11/R11[0]
    ax.plot(tau,R11_norm,**plot_options)
    ax.set_xlabel(r'$\tau/\tau_\mathrm{ref}$',fontsize=24)
    direction=calc_options['direction']
    if direction=='u':
        ax.set_ylabel(r'$R_{11}(\tau)$',fontsize=24)
    elif direction=='v':
        ax.set_ylabel(r'$R_{22}(\tau)$',fontsize=24)
    elif direction=='w':
        ax.set_ylabel(r'$R_{33}(\tau)$',fontsize=24)
    # calc tau integral scale
    T=np.trapz(R11_norm,x=tau)
    print 'T= ',T
    return R11_norm,tau,T
def save_R11_tau(f_handle,i=0,j=0,calc_options={'direction':'u'}):
    t,X,Y,Z,R11=calc_R11_tau_exp(i=i,j=j,**calc_options)
    tau=t[:len(t)/2]
    # normalize by variance
    R11_norm= R11/R11[0]
    # calc tau integral scale
    T=np.trapz(R11_norm,x=tau)
    np.savetxt(f_handle,[[X,Y,T]],delimiter=' ')
    direction=calc_options['direction']
    if direction=='u':
        np.save('R11_output/XYR11/ijXYR11_'+str(i)+'_'+str(j)+'_%.4f_%.4f'%(X,Y),R11_norm)
    elif direction=='v':
        np.save('R22_output/XYR22/ijXYR22_'+str(i)+'_'+str(j)+'_%.4f_%.4f'%(X,Y),R11_norm)
    elif direction=='w':
        np.save('R33_output/XYR33/ijXYR33_'+str(i)+'_'+str(j)+'_%.4f_%.4f'%(X,Y),R11_norm)
    print 'saved ijxy=',i,j,X,Y
def Manufactured_R11(t,tau,option=1):
    L=t[-1]-t[0]
    # option 1
    if option==1:
        R11=(
            1./(L*np.exp(tau)) 
            * ( 
                (np.cos(tau)/5. 
                * (
                np.exp(-1.*L) * (-1.*np.cos(L)**2 + 2.*np.sin(L)*np.cos(L)-2.)
                + 3.))
            - (np.sin(tau)/2.  * ( np.sin(L)**2))
            )
            ) 
    # option 2
    elif option==2:
        R11=(1./L)*(
                np.cos(tau)
                * (( np.sin(L)      *np.cos(L)+L )  /   2.)
                - ((  np.sin(tau)   *np.sin(L)**2)  /   2.))
    #print 'manufactured R11        =',R11
    R11_norm=R11/R11[0]
    #print 'manufactured R11/R11[0] =',R11_norm
    # return normalized R11
    return R11_norm
def R11_tau_calc_plot(ax,calc_options,plot_options):
    # Experimental data set
    #calc_options1={
            #'direction':direction,
            #'manufactured':manufactured,
            #'method':method,
            #'option':option,
            #'skip_last':skip_last,
            #'moving_mean':moving_mean,
            #'average_z':average_z,
            #}
    #plot_options2=calc_options1.copy()
    #plot_options2['method']=1
    #plot_options2['label']='exp method=1,option=%i'%option
    manufactured=calc_options['manufactured']
    if manufactured:
        #t=np.linspace(0,200*np.pi,5000)
        t=np.arange(0,0.5,0.0001)/(0.00745/27.5)
        #t=np.linspace(0,0.5/(0.00745/27.5),5000)
        #R11exp1=plot_R11_tau_exp(ax,xi,yi,direction=direction,method=0,manufactured=manufactured,option=option,skip_last=skip_last,moving_mean=moving_mean,marker='+',linewidth=0,label='exp method=0,option=%i'%option)
        #R11exp2=plot_R11_tau_exp(ax,xi,yi,direction=direction,method=1,manufactured=manufactured,option=option,skip_last=skip_last,moving_mean=moving_mean,marker='+',linewidth=0,label='exp method=1,option=%i'%option)
        R11exp1,_tau,T=plot_R11_tau_exp(ax,xi,yi,calc_options,plot_options)
        #R11exp2=plot_R11_tau_exp(ax,xi,yi,**plot_options2)
        #plot_R11_tau_exp(ax,-1,0,direction=direction,method=2,manufactured=manufactured,option=option,skip_last=skip_last,moving_mean=moving_mean,marker='+',label='method=1,option=%i'%option)
        R11expMNF=Manufactured_R11(t,t[:len(t)/2],option=option)
        ax.plot(t[:len(t)/2],R11expMNF,marker='x',linewidth=0,label='exp Manufactured')
    else:
        #R11exp1=plot_R11_tau_exp(ax,xi,yi,direction=direction,method=0,manufactured=manufactured,option=option,skip_last=skip_last,moving_mean=moving_mean,marker='+',linewidth=0,label='exp method=0,option=%i'%option)
        #R11exp2=plot_R11_tau_exp(ax,xi,yi,direction=direction,method=1,manufactured=manufactured,option=option,skip_last=skip_last,moving_mean=moving_mean,marker='+',linewidth=0,label='exp method=1,option=%i'%option)
        R11exp1,_tau,T=plot_R11_tau_exp(ax,xi,yi,calc_options,plot_options)
        #R11exp2=plot_R11_tau_exp(ax,xi,yi,**plot_options2)
        #plot_R11_tau_exp(ax,-1,0,direction=direction,method=2,manufactured=manufactured,option=option,skip_last=skip_last,marker='x',label='method=1,option=%i'%option)
        R11expMNF=0.
    labelfont={
            'fontsize':14,
            'fontname':'Times New Roman'}
    ax.set_xlabel(r'$\tau/\tau_\mathrm{ref}$',**labelfont)
    ax.set_ylabel(r'$R_{11}(\tau)$',**labelfont)
    plt.tight_layout()
    ax.legend(loc='best',numpoints=1)
    #plt.show()
    #return R11exp1,R11exp2,R11expMNF
    return R11exp1,_tau,R11expMNF,T
if __name__=="__main__":
    # user inputs
    direction='u'   # ['u','v','w']
    option=2        # [1,2]
    method=0
    xi,yi=-3,48     # [0-169,0-126]
    manufactured=False
    skip_last=False
    moving_mean=False
    average_z=True

    fig=plt.figure()
    ax=plt.subplot(111)
    R11_tau_calc_plot(ax)
    plt.show()
