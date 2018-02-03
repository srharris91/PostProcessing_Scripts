import time
import os
import numpy as np
import matplotlib.pyplot as plt
import h5py
plt.rcParams["font.family"] = "Times New Roman"

f1 = h5py.File('../DA_ExperimentalData/From_JFrank_5000_cold_SandiaC/T093013_01_Velocity_1-500.mat')
f2 = h5py.File('../DA_ExperimentalData/From_JFrank_5000_cold_SandiaC/T093013_01_Velocity_501-1000.mat')
f3 = h5py.File('../DA_ExperimentalData/From_JFrank_5000_cold_SandiaC/T093013_01_Velocity_1001-1500.mat')
f4 = h5py.File('../DA_ExperimentalData/From_JFrank_5000_cold_SandiaC/T093013_01_Velocity_1501-2000.mat')
f5 = h5py.File('../DA_ExperimentalData/From_JFrank_5000_cold_SandiaC/T093013_01_Velocity_2001-2500.mat')
f6 = h5py.File('../DA_ExperimentalData/From_JFrank_5000_cold_SandiaC/T093013_01_Velocity_2501-3000.mat')
f7 = h5py.File('../DA_ExperimentalData/From_JFrank_5000_cold_SandiaC/T093013_01_Velocity_3001-3500.mat')
f8 = h5py.File('../DA_ExperimentalData/From_JFrank_5000_cold_SandiaC/T093013_01_Velocity_3501-4000.mat')
f9 = h5py.File('../DA_ExperimentalData/From_JFrank_5000_cold_SandiaC/T093013_01_Velocity_4001-4500.mat')
f10 = h5py.File('../DA_ExperimentalData/From_JFrank_5000_cold_SandiaC/T093013_01_Velocity_4501-5000.mat')
f_all=[f1,f2,f3,f4,f5,f6,f7,f8,f9,f10]
f1.keys()

def read_data_pt(i=0,j=0,direction='u'):
    if direction=='u':
        Us=np.concatenate([np.array(f['Us'][:,:,i,-1-j]) for f in f_all],axis=0)
    elif direction=='v':
        Us=np.concatenate([np.array(f['Vs'][:,:,i,-1-j]) for f in f_all],axis=0)
    elif direction=='w':
        Us=np.concatenate([np.array(f['Ws'][:,:,i,-1-j]) for f in f_all],axis=0)
    else:
        print 'direction is wrong at ',direction
    X=np.array(f['X'][0,i])
    Y=np.array(f['Y'][0,-1-j])
    Z=np.array(f['Z'][0,:])
    # recenter and change to m
    D = 0.00745 # diameter of jet in m
    Y=Y/(1000.)+15.*D # mm to m and 15 diam. downstream
    X=X/1000.-(1.15*D) # centered between nozzles
    Z=Z/1000.
    return [X,Y,Z,Us]
def read_data_time(t=0):
    Us=np.array(f['Us'][t])
    Vs=np.array(f['Vs'][t])
    Ws=np.array(f['Ws'][t])
    X=np.array(f['X'][0,])
    Y=np.array(f['Y'][0,])
    Z=np.array(f['Z'][0,])
    # recenter and change to m
    D = 0.00745 # diameter of jet in m
    Y=Y/(1000.)+15.*D # mm to m and 15 diam. downstream
    X=X/1000.-(1.15*D) # centered between nozzles
    Z=Z/1000.
    Zmesh,Xmesh,Ymesh=np.meshgrid(Z,X,Y,indexing='ij')
    return [Xmesh,Ymesh,Zmesh,Us,Vs,Ws]
def read_data_U():
    Us=np.array(f['Us'])
    return Us
def read_data_V():
    Vs=np.array(f['Vs'])
    return Vs
def read_data_W():
    Ws=np.array(f['Ws'])
    return Ws
def plot_time_show(time=0,z=3):
    X,Y,Z,U,V,W=read_data_time(t=time)
    plt.figure()
    plt.colorbar(plt.contourf(X[z],Y[z],V[z],30,cmap='jet'))
    plt.title('time=%.4d'%time+' z=%.4d'%z)
    plt.xlabel('x',fontsize=24)
    plt.ylabel('y',fontsize=24)
    #plt.savefig('Continuity/continuity_time=%.4i'%(time)+'_z=%.4i'%(z)+'.png', dpi=200,bbox_inches='tight',)
    #plt.close()
    plt.show()
#plot_time_show()
def calc_R11_tau_exp(t=np.arange(0,0.49991,0.0001),i=0,j=0,direction='u',method=0,manufactured=False,option=1):
    X,Y,Z,Us=read_data_pt(i=i,j=j,direction=direction)
    # calc u' v' and w'  velocity prime values
    u_avg=Us.mean(axis=(0,1))
    u=Us-u_avg
    
    # manufactured data (to check accuracy of method)
    #if manufactured==True:
        #length=len(u)
        #theta=np.linspace(0.,2.*np.pi,length)
        #u=np.cos(theta)
    if manufactured==True:
        # option 1
        if option==1:
            u=(np.cos(t)*np.exp(-1.*t)).reshape((len(t),1))*np.ones(u.shape)
        # option 2
        elif option==2:
            u=np.cos(t).reshape((len(t),1))*np.ones(u.shape)
        # plot manufactured solution
        fig=plt.figure()
        ax=plt.subplot(111)
        ax.plot(t,u[:,0],'.')
        ax.set_xlabel('t')
        ax.set_ylabel("u'")
        #plt.show()
    
    nt=2500 # half of time domain
    R11=np.zeros(nt)
    if method==0:
        tau=np.arange(0,nt)
        for taui in tau:
            R11[taui]=(u[:nt]*u[taui:nt+taui]).mean(axis=(0,1))
    elif method==1:
        R11[0]= (u*u).mean(axis=(0,1))
        tau=np.arange(1,nt)
        for taui in tau: 
            R11[taui]=((u[:(-1*taui)]*u[taui:]).mean(axis=(0,1)))
    elif method==2:
        tau=np.arange(0,nt)
        for tauispace in tau:
            for tauime in tau:
                R11[tauispace]=R11[tauispace]+(u[tauime]*u[(tauime+tauispace)]).mean()
        R11=R11/float(nt)
    return X,Y,Z,R11
def plot_R11_tau_exp(ax,i=0,j=0,t=np.arange(0,0.49991,0.0001),direction='u',method=0,manufactured=False,option=1,**kwargs):
    #tau=np.arange(0,0.24991,0.0001)
    tau=t[:len(t)/2]
    X,Y,Z,R11=calc_R11_tau_exp(i=i,j=j,t=t,direction=direction,method=method,manufactured=manufactured,option=option)
    #fig=plt.figure()
    #ax=plt.subplot(111)
    # normalize by variance
    R11_norm= R11/R11[0]
    #ax.plot(tau,R11_norm,'k.')
    ax.plot(tau,R11_norm,**kwargs)
    ax.set_xlabel(r'$\tau[s]$',fontsize=24)
    if direction=='u':
        ax.set_ylabel(r'$R_{11}(\tau)$',fontsize=24)
    elif direction=='v':
        ax.set_ylabel(r'$R_{22}(\tau)$',fontsize=24)
    elif direction=='w':
        ax.set_ylabel(r'$R_{33}(\tau)$',fontsize=24)
    if manufactured==True:
        ax.set_title('cos wave in time,dir='+direction+',ij=(%d,%d) and xyz=(%.4f,%.4f)'%(i,j,X,Y))
    else:
        ax.set_title('dir='+direction+',ij=(%d,%d) and xy=(%.4f,%.4f)'%(i,j,X,Y))

    # calc tau integral scale
    T=np.trapz(R11_norm,x=tau)
    print 'T= ',T
    print 'Calc R11=',R11
    print 'Call R11/R11[0]=',R11_norm
    #plt.show()
def save_R11_tau(f_handle,i=0,j=0,direction='u',manufactured=False):
    X,Y,Z,R11=calc_R11_tau_exp(i=i,j=j,direction=direction,manufactured=manufactured)
    tau=np.arange(0,0.24991,0.0001)
    # normalize by variance
    R11_norm= R11/R11[0]
    # calc tau integral scale
    T=np.trapz(R11_norm,x=tau)
    np.savetxt(f_handle,[[X,Y,T]],delimiter=' ')
    np.save('R11_output/XYR11/ijXYR11_'+str(i)+'_'+str(j)+'_%.4f_%.4f'%(X,Y),R11_norm)
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
    print 'manufactured R11        =',R11
    R11_norm=R11/R11[0]
    print 'manufactured R11/R11[0] =',R11_norm
    # return normalized R11
    return R11_norm
###################################################
###################################################
###################################################

def read_box_probe(direction=0):
    if direction==0:
        ux='U-X'
    elif direction==1:
        ux='U-Y'
    elif direction==2:
        ux='U-Z'
    else:
        print 'Problem with direction in read_box_probe, direction=',direction
    fn=[
    'box_probe.'+ux,
    'box_probe.X',
    'box_probe.Y',
    'box_probe.Z',
    ]
    box_probe={}
    for fi in fn:
        print "reading ",fi.split(".")[-1]
        box_probe[fi.split(".")[-1]]=np.genfromtxt(fi,delimiter=' ')
    u=box_probe[ux][:,3:]
    x=box_probe['X'][3:]
    y=box_probe['Y'][3:]
    z=box_probe['Z'][3:]
    t=box_probe[ux][:,1]
    X=np.swapaxes(x.reshape((2,10,20)),0,2)
    Y=np.swapaxes(y.reshape((2,10,20)),0,2)
    Z=np.swapaxes(z.reshape((2,10,20)),0,2)
    U=np.swapaxes(u.reshape((len(t),2,10,20)),1,3)
    
    # non-dimensionalize everything
    #   D_ref=0.00745
    #   U_ref=27.5

    #   X=X/D_ref
    #   Y=Y/D_ref
    #   Z=Z/D_ref
    #   U=U/U_ref
    #   t=t/(D_ref/U_ref)

    # check if reshape is same values
    #    print u.shape
    #    print t.shape
    #    print t
    #    print x.shape
    #    print U.shape
    #    loc=(9,7,0)
    #    num=loc[0]+20*loc[1]+20*10*loc[2]
    #    print loc,num
    #    print U[(8,)+loc] == u[8,num]
    #    print X[loc] == x[num]
    #    print Y[loc] == y[num]
    #    print Z[loc] == z[num]
    #    print 'if True, then reshape is correct and ordering is (t,x,y,z)'
     
    return t,X,Y,Z,U

def plot_csv(gridx,gridy,gridu,direction=0,colormap=True):
    fig=plt.figure()
    ax=plt.subplot(111)
    if colormap:
        if direction==0:
            vmin,vmid,vmax=0.,0.35,0.7
        elif direction==1:
            vmin,vmid,vmax=-0.05,0.,0.05
        elif direction==2:
            vmin,vmid,vmax=-0.05,0.,0.05
        colorvalues=np.linspace(vmin,vmax,100)
        ax.contour (gridx[:,:,0],gridy[:,:,0],gridu[:,:,0],colorvalues,cmap='jet',vmin=vmin,vmax=vmax)
        plt.colorbar(
        ax.contourf(gridx[:,:,0],gridy[:,:,0],gridu[:,:,0],colorvalues,cmap='jet',vmin=vmin,vmax=vmax)
        ,ticks=[vmin,vmid,vmax])
    else:
        ax.contour (gridx[:,:,0],gridy[:,:,0],gridu[:,:,0],30,cmap='jet')
        plt.colorbar(
        ax.contourf(gridx[:,:,0],gridy[:,:,0],gridu[:,:,0],30,cmap='jet')
        ,)
    ax.axis('equal')
    ax.set_title('instantaneous snapshot')
    labelfont={
            'fontsize':14,}
    plt.xlabel('X',**labelfont)
    plt.ylabel('Y',**labelfont)
    plt.tight_layout()
    #plt.savefig(output_csv+'.pdf',bbox_inches='tight')
    #plt.show()
def plot_rms(gridx,gridy,gridu,direction=0):
    homogeneous=(0,-1) # time and z
    fig=plt.figure()
    ax=plt.subplot(111)
    average=u.mean(axis=homogeneous)[np.newaxis,:,:,np.newaxis]
    ax.contour                  (gridx[:,:,0],gridy[:,:,0],average[0,:,:,0],30,cmap='jet')
    plt.colorbar( ax.contourf   (gridx[:,:,0],gridy[:,:,0],average[0,:,:,0],30,cmap='jet') ,)
    ax.axis('equal')
    ax.set_title('average')
    labelfont={ 'fontsize':14,}
    ax.set_xlabel('X',**labelfont)
    ax.set_ylabel('Y',**labelfont)
    plt.tight_layout()

    fig=plt.figure()
    ax=plt.subplot(111)
    u_prime=u-average
    urms=(u_prime**2).mean(axis=homogeneous)
    ax.contour                  (gridx[:,:,0],gridy[:,:,0],urms,30,cmap='jet')
    plt.colorbar( ax.contourf   (gridx[:,:,0],gridy[:,:,0],urms,30,cmap='jet') ,)
    ax.axis('equal')
    ax.set_title('urms')
    labelfont={ 'fontsize':14,}
    ax.set_xlabel('X',**labelfont)
    ax.set_ylabel('Y',**labelfont)
    plt.tight_layout()
    #plt.savefig(output_csv+'.pdf',bbox_inches='tight')
    #plt.show()

def calc_R11_tau_sim(t,gridu,direction=0,method=0,option=1,manufactured=False):
    homogeneous=(0,-1) # homogeneous axis are time=0 and z=-1
    
    # calc u' v' and w'  velocity prime values
    u_avg=gridu.mean(axis=(0,-1))
    u=gridu-u_avg
    
    # manufactured data (to check accuracy of method)
    if manufactured==True:
        # option 1
        if option==1:
            u=(np.cos(t)*np.exp(-1.*t)).reshape((len(t),1))*np.ones(u.shape)
        # option 2
        elif option==2:
            u=np.cos(t).reshape((len(t),1))*np.ones(u.shape)
    
    nt=int(gridu.shape[0]/2.) # half of time domain
    R11=np.zeros(nt)
    if method==0:
        tau=np.arange(0,nt)
        for taui in tau:
            R11[taui]=(u[:nt+1]*u[taui:nt+taui+1]).mean(axis=(homogeneous))
    elif method==1:
        tau=np.arange(1,nt)
        R11[0]= (u*u).mean(axis=(homogeneous))
        for taui in tau: 
            R11[taui]=((u[:-1*taui]*u[taui:]).mean(axis=(homogeneous)))
    elif method==2:
        tau=np.arange(0,nt)
        for tauispace in tau:
            for tauime in tau:
                R11[tauispace]=R11[tauispace]+(u[tauime]*u[(tauime+tauispace)]).mean()
        R11=R11/float(nt)
    else:
        print 'get right method...',method
    # save value
    #np.save('R11'+direction,R11)
    return R11
def plot_R11_tau_sim(ax,X,Y,t,gridu,method=0,direction=0,plot=True,manufactured=False,option=1):
    R11=calc_R11_tau_sim(t,gridu,direction=direction,method=method,manufactured=manufactured,option=option)
    if manufactured==True:
        title=          'cos wave in time,dir='+str(direction)+',and xyz=(%.4f,%.4f) %i'%(X,Y,method)
    else:
        title =         'dir='+str(direction)+', xy=(%.4f,%.4f) %i'%(X,Y,method)
    #tau=np.arange(0,0.24991,0.0001)
    #dt=0.5E-7 * 10000. # simulation time step between vtu files
    #tau=np.arange(0,dt*gridu.shape[0]/2.,dt)
    tau=t[:len(t)/2]-t[0]
    # normalize by variance
    R11_norm= R11/R11[0]
    if plot==True:
        #fig=plt.figure()
        #ax=plt.subplot(111)
        ax.plot(tau,R11_norm,'+',label=title)
        labelfont={
                'fontsize':14,
                'fontname':'Times New Roman'}
        ax.set_xlabel(r'$\tau[s]$',**labelfont)
        if direction==0:
            ax.set_ylabel(r'$R_{11}(\tau)$',**labelfont)
        elif direction==1:
            ax.set_ylabel(r'$R_{22}(\tau)$',**labelfont)
        elif direction==2:
            ax.set_ylabel(r'$R_{33}(\tau)$',**labelfont)
        ax.set_title(title,**labelfont)
        plt.tight_layout()
        print 'calc R11        = ',R11
        #print 'calc R11/R11[0] = ',R11_norm
        #np.savetxt('R11.txt',np.array([tau,R11_norm]).T)

    # calc tau integral scale
    T=np.trapz(R11_norm,x=tau)
    #print 'T= ',T
    return T
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
    print 'manufactured R11        =',R11
    R11_norm=R11/R11[0]
    print 'manufactured R11/R11[0] =',R11_norm
    # return normalized R11
    return R11_norm

###################################################
###################################################
###################################################

if __name__=="__main__":
    ###################################################
    # Experimental data set
    # user inputs
    option=2
    xi,yi=-3,48   # [0-169,0-126]
    manufactured=False
    fig=plt.figure()
    ax=plt.subplot(111)
    if manufactured:
        t=np.linspace(0,200*np.pi,5000)
        plot_R11_tau_exp(ax,xi,yi,t=t,direction='u',method=0,manufactured=manufactured,option=option,marker='+',label='method=0,option=%i'%option)
        plot_R11_tau_exp(ax,xi,yi,t=t,direction='u',method=1,manufactured=manufactured,option=option,marker='+',label='method=1,option=%i'%option)
        #plot_R11_tau_exp(ax,-1,0,t=t,direction='u',method=2,manufactured=manufactured,option=option,marker='+',label='method=1,option=%i'%option)
        ax.plot(t[:len(t)/2],Manufactured_R11(t,t[:len(t)/2],option=option),label='Manufactured')
    else:
        plot_R11_tau_exp(ax,xi,yi,direction='u',method=0,manufactured=manufactured,option=option,marker='+',label='method=0,option=%i'%option)
        plot_R11_tau_exp(ax,xi,yi,direction='u',method=1,manufactured=manufactured,option=option,marker='+',label='method=1,option=%i'%option)
        #plot_R11_tau_exp(ax,-1,0,direction='u',method=2,manufactured=manufactured,option=option,marker='x',label='method=1,option=%i'%option)
    ax.set_xlabel(r'$\tau[s]$',fontsize=14)
    ax.set_ylabel(r'$R_{11}(\tau)$',fontsize=14)
    plt.tight_layout()
    ax.legend(loc='best',numpoints=1)

###################################################
###################################################
###################################################

    # Simulation data
    os.chdir('../Create_Plots/DA_paper/OI/C_0/SOLUT/')

    # user options
    direction=1 # [0,1,2]
    option=2    # [1,2]
    manufactured=False   # [True,False]

    xi,yi=4,9   # [0-19,0-10]
    t,x,y,z,u=read_box_probe(direction=direction)

    # for manufactured solution only
    if manufactured:
        t=np.linspace(0,50*np.pi,len(t))

    # calculate and plot R11(tau)
    T,X=[],[]
    #for i in np.arange(20):
    for i in [xi,]:
        if i==xi:
            T.append(plot_R11_tau_sim(ax,x[i,yi,0],y[i,yi,0],t,u[:,i,yi,:],method=0,direction=direction,plot=True,manufactured=manufactured,option=option))
            T.append(plot_R11_tau_sim(ax,x[i,yi,0],y[i,yi,0],t,u[:,i,yi,:],method=1,direction=direction,plot=True,manufactured=manufactured,option=option))
            #T.append(plot_R11_tau_sim(ax,x[i,yi,0],y[i,yi,0],t,u[:,i,yi,:],method=2,direction=direction,plot=True,manufactured=manufactured,option=option))
        else:
            T.append(plot_R11_tau_sim(ax,x[i,yi,0],y[i,yi,0],t,u[:,i,yi,:],method=method,direction=direction,plot=False,manufactured=True,option=option))
        X.append(x[i,yi,0])
    #plot_R11_tau_sim(ax,x[0,0,0],y[0,0,0],t,u[:,0,0,:],method=1,direction=direction,manufactured=False)
    taus=t[:len(t)/2]-t[0]
    if manufactured:
        ax.plot(taus,Manufactured_R11(t-t[0],taus,option=option),label='Manufactured')
    ax.legend(loc='best',numpoints=1)

    plt.tight_layout()
    plt.show()
