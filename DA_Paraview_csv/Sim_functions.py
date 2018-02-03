import time
import os
import numpy as np
import matplotlib.pyplot as plt
import h5py
sim_dir= '../Create_Plots/DA_paper/OI/C_0/SOLUT/'
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
        print '  simulation reading '+fi.split(".")[-1]+' from '+sim_dir
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
    D_ref=0.00745
    U_ref=27.5
    X=X/D_ref
    Y=Y/D_ref
    Z=Z/D_ref
    U=U/U_ref
    t=t/(D_ref/U_ref)
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

def calc_R11_tau_sim(t,gridu,direction=0,method=0,option=1,manufactured=False):
    homogeneous=(0,-1) # homogeneous axis are time=0 and z=-1
    
    # calc u' v' and w'  velocity prime values
    u_avg=gridu.mean(axis=homogeneous)
    u=gridu-u_avg
    
    # manufactured data (to check accuracy of method)
    if manufactured==True:
        # option 1
        if option==1:
            u=(np.cos(t)*np.exp(-1.*t)).reshape((len(t),1)) + 0.08*(2.*np.random.random((len(t),)+u.shape[1:])-1.0)
        # option 2
        elif option==2:
            u=np.cos(t).reshape((len(t),1)) + 0.08*(2.*np.random.random((len(t),)+u.shape[1:])-1.0)
    
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
    return R11
def plot_R11_tau_sim(ax,X,Y,t,gridu,method=0,direction=0,plot=True,manufactured=False,option=1,**kwargs):
    R11=calc_R11_tau_sim(t,gridu,direction=direction,method=method,manufactured=manufactured,option=option)
    if manufactured==True:
        title=          'cos wave in time,dir='+str(direction)+',and xyz=(%.4f,%.4f) %i'%(X,Y,method)
    else:
        title =         'dir='+str(direction)+', xy=(%.4f,%.4f) %i'%(X,Y,method)
    tau=t[:len(t)/2]-t[0]
    # normalize by variance
    R11_norm= R11/R11[0]
    if plot==True:
        ax.plot(tau,R11_norm,**kwargs)
        labelfont={
                'fontsize':14,
                'fontname':'Times New Roman'}
        ax.set_xlabel(r'$\tau/\tau_\mathrm{ref}$',**labelfont)
        if direction==0:
            ax.set_ylabel(r'$R_{11}(\tau)$',**labelfont)
        elif direction==1:
            ax.set_ylabel(r'$R_{22}(\tau)$',**labelfont)
        elif direction==2:
            ax.set_ylabel(r'$R_{33}(\tau)$',**labelfont)
        plt.tight_layout()

    # calc tau integral scale
    T=np.trapz(R11_norm,x=tau)
    #print 'T= ',T
    return T,R11_norm
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
    R11_norm=R11/R11[0]
    return R11_norm

def R11_tau_calc_plot(ax):
    os.chdir(sim_dir)

    t,x,y,z,u=read_box_probe(direction=direction)

    # for manufactured solution only
    if manufactured:
        t=np.linspace(0,0.5/(0.00745/27.5),len(t))

    # calculate and plot R11(tau)
    T,X=[],[]
    for i in [xi,]:
        if i==xi:
            T1,R11sim1=plot_R11_tau_sim(ax,x[i,yi,0],y[i,yi,0],t,u[:,i,yi,:],method=method,direction=direction,plot=True,manufactured=manufactured,option=option,marker='+',linewidth=0,label='sim method=0,option=%i'%option)
            T.append(T1)
        else:
            T1,R11sim=(plot_R11_tau_sim(ax,x[i,yi,0],y[i,yi,0],t,u[:,i,yi,:],method=method,direction=direction,plot=False,manufactured=True,option=option))
            T.append(T1)
        X.append(x[i,yi,0])
    taus=t[:len(t)/2]-t[0]
    if manufactured:
        R11simMNF=Manufactured_R11(t-t[0],taus,option=option)
        ax.plot(taus,R11simMNF,marker='x',linewidth=0,label='sim Manufactured')
    else:
        R11simMNF=0.
    ax.legend(loc='best',numpoints=1)

    plt.tight_layout()
    #plt.show()
    #return R11sim1,R11sim2,R11simMNF
    return R11sim1,R11simMNF
if __name__=="__main__":
    # user options
    direction=1 # [0,1,2]
    option=2    # [1,2]
    method=0    # [0,1]
    manufactured=False
    xi,yi=4,9   # [0-19,0-10]

    fig=plt.figure()
    ax=plt.subplot(111)
    R11_tau_calc_plot(ax)
    plt.show()
