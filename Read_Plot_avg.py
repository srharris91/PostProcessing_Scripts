import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy.interpolate import griddata
import h5py
plt.style.use('seaborn-paper')

class sim_avg:
    def __init__(self,csv_filename,save_directory):
        self.csv_filename=csv_filename
        self.save_directory=save_directory

    def csv_slice_to_npy(self):
        # read in data
        d=np.genfromtxt(self.csv_filename,delimiter=',',skip_header=1)
        uexp,u,v,w,uavg,vavg,wavg,urms,vrms,wrms,x,y,z=[ d[:,i] for i in range(0,13) ]

        # non-dimensionalize everything
        D_ref = 0.00745
        U_ref = 27.5
        x = x/D_ref
        y = y/D_ref
        z = z/D_ref
        uexp=uexp/U_ref
        u=u/U_ref
        v=v/U_ref
        w=w/U_ref
        uavg=uavg/U_ref
        vavg=vavg/U_ref
        wavg=wavg/U_ref
        urms=urms/U_ref
        vrms=vrms/U_ref
        wrms=wrms/U_ref

        # interpolate to grid
        nx=3000j
        ny=260j
        self.Xall,self.Yall=np.mgrid[0:30:nx,-4:4:ny] # make grid of slice down the middle (if more data is given, can do slice of whole domain to get the whole experimental region)
        self.XCL,self.YCL=np.mgrid[0:30:nx,0:1:1j] # make grid of slice down the middle (if more data is given, can do slice of whole domain to get the whole experimental region)
        points=np.array([x,y]).T
        self.Uexp=griddata(points,uexp,(self.Xall,self.Yall),method='linear')
        self.Urms=griddata(points,urms,(self.Xall,self.Yall),method='linear')
        self.Uavg=griddata(points,uavg,(self.Xall,self.Yall),method='linear')
        self.Vavg=griddata(points,vavg,(self.Xall,self.Yall),method='linear')
        self.UCL=griddata(points,uavg,(self.XCL,self.YCL),method='linear')
        self.VCL=griddata(points,vavg,(self.XCL,self.YCL),method='linear')
        self.UrmsCL=griddata(points,urms,(self.XCL,self.YCL),method='linear')
        self.VrmsCL=griddata(points,vrms,(self.XCL,self.YCL),method='linear')

        # save interpolated data
        np.save(self.save_directory+'/uexp',self.Uexp)
        np.save(self.save_directory+'/uavg',self.Uavg)
        np.save(self.save_directory+'/vavg',self.Vavg)
        np.save(self.save_directory+'/urms',self.Urms)
        np.save(self.save_directory+'/ucl',self.UCL)
        np.save(self.save_directory+'/vcl',self.VCL)
        np.save(self.save_directory+'/urmscl',self.UrmsCL)
        np.save(self.save_directory+'/vrmscl',self.VrmsCL)
        np.save(self.save_directory+'/X',self.Xall)
        np.save(self.save_directory+'/Y',self.Yall)
        np.save(self.save_directory+'/Xcl',self.XCL)
        np.save(self.save_directory+'/Ycl',self.YCL)

    def load(self):
        # load data
        self.Uavg=np.load(self.save_directory+'/uavg.npy')
        self.Vavg=np.load(self.save_directory+'/vavg.npy')
        self.Urms=np.load(self.save_directory+'/urms.npy')
        self.Uexp=np.load(self.save_directory+'/uexp.npy')
        self.UCL =np.load(self.save_directory+'/ucl.npy')
        self.VCL =np.load(self.save_directory+'/vcl.npy')
        self.UrmsCL =np.load(self.save_directory+'/urmscl.npy')
        self.VrmsCL =np.load(self.save_directory+'/vrmscl.npy')
        self.Xall=np.load(self.save_directory+'/X.npy')
        self.Yall=np.load(self.save_directory+'/Y.npy')
        self.XCL=np.load(self.save_directory+'/Xcl.npy')
        self.YCL=np.load(self.save_directory+'/Ycl.npy')
        # change exp data to nan and 1000
        temp=self.Uexp.copy()
        temp[self.Uexp==0]=np.nan
        temp[np.isfinite(temp)]=1000
        # get bounding values for exp data
        self.exp_range=np.isfinite(temp)
        self.xmin=np.min(self.Xall[self.exp_range])
        self.xmax=np.max(self.Xall[self.exp_range])
        self.ymin=np.min(self.Yall[self.exp_range])
        self.ymax=np.max(self.Yall[self.exp_range])
        self.xwidth=self.xmax-self.xmin
        self.ywidth=self.ymax-self.ymin

class exp_avg:
    def __init__(self,ave_u_npy,ave_v_npy,rms_u_npy,rms_v_npy,mat_filename):
        self.ave_u_npy=ave_u_npy
        self.ave_v_npy=ave_v_npy
        self.rms_u_npy=rms_u_npy
        self.rms_v_npy=rms_v_npy
        self.mat_filename=mat_filename
    def load(self):
        expv=np.load(self.ave_v_npy)
        expu=np.load(self.ave_u_npy)
        rmsv=np.load(self.rms_v_npy)
        rmsu=np.load(self.rms_u_npy)
        # read in X and Y
        f=h5py.File(mat_filename)
        expX = np.array(f['X'][0,:])
        expY = np.array(f['Y'][0,:])
        D_ref=0.00745
        U_ref=27.5
        expY=expY/(1000.)+15.*D_ref # mm to m and 15 diam. downstream
        expX=expX/1000.-(1.15*D_ref) # centered between nozzles
        # non-dimensionalize everything
        expX = expX/D_ref
        expY = expY/D_ref
        expv = expv/U_ref
        expu = expu/U_ref
        rmsv = rmsv/U_ref
        rmsu = rmsu/U_ref
        # make mesh grid for plotting
        expX,expY=np.meshgrid(expX,expY,indexing='ij')

        # save non-dimensional data to self
        self.u=expu
        self.v=expv
        self.rmsu=rmsu
        self.rmsv=rmsv
        self.x=expX
        self.y=expY
    def center_line(self):
        self.XCL,self.YCL=np.mgrid[0:1:1j,self.y.min():self.y.max():100j] # make grid of slice down the middle (if more data is given, can do slice of whole domain to get the whole experimental region)
        points=np.array([self.x.flatten(),self.y.flatten()]).T
        self.VCL=griddata(points,self.v.flatten(),(self.XCL,self.YCL),method='linear')
        self.UCL=griddata(points,self.u.flatten(),(self.XCL,self.YCL),method='linear')
        self.VrmsCL=griddata(points,self.rmsv.flatten(),(self.XCL,self.YCL),method='linear')
        self.UrmsCL=griddata(points,self.rmsu.flatten(),(self.XCL,self.YCL),method='linear')
        
def plot_data(sim,exp,direction):
    #fig = plt.figure(figsize=(5.4,5))
    fig = plt.figure(figsize=(8.36,5))
    axavg = plt.subplot2grid((4,140),(0,0),rowspan=4,colspan=40,aspect='auto')
    if direction=='u':
        colormap=np.linspace(np.nanmin(sim.Uavg),np.nanmax(sim.Uavg),300)
        plt.colorbar(axavg.contourf(-sim.Yall,sim.Xall,sim.Uavg,colormap,cmap='jet'),ax=axavg,ticks=[0,1.25],)
        axexp1 = plt.subplot2grid((4,140),(0,40),rowspan=2,colspan=60)
        axexp2 = plt.subplot2grid((4,140),(2,40),rowspan=2,colspan=60,sharex=axexp1,sharey=axexp1)
        axexp3 = plt.subplot2grid((4,140),(0,100),rowspan=4,colspan=40,sharey=axavg)
        #colormap=np.linspace(np.nanmin(sim.Uavg[sim.exp_range]),np.nanmax(sim.Uavg[sim.exp_range]),300)
        colormap=np.linspace(0.,np.nanmax(sim.Uavg[sim.exp_range]),300)
        plt.colorbar(axexp1.tricontourf(-sim.Yall[sim.exp_range],sim.Xall[sim.exp_range],sim.Uavg[sim.exp_range],colormap,cmap='jet'),ax=axexp1,ticks=[0.,0.5])
        plt.colorbar(axexp2.contourf(-exp.x,exp.y,exp.v,colormap,cmap='jet',aspect='auto'),ax=axexp2,ticks=[0.,0.5])
    elif direction=='v':
        #colormap=np.linspace(np.nanmin(sim.Urms),np.nanmax(sim.Urms),300)
        colormap=np.linspace(0,np.nanmax(sim.Urms),300)
        plt.colorbar(axavg.contourf(-sim.Yall,sim.Xall,sim.Urms,colormap,cmap='jet'),ax=axavg,ticks=[0.,0.2])
        axexp1 = plt.subplot2grid((4,140),(0,40),rowspan=2,colspan=60)
        axexp2 = plt.subplot2grid((4,140),(2,40),rowspan=2,colspan=60,sharex=axexp1,sharey=axexp1)
        axexp3 = plt.subplot2grid((4,140),(0,100),rowspan=4,colspan=40,sharey=axavg)
        #colormap=np.linspace(np.nanmin(exp.u),np.nanmax(sim.Urms[sim.exp_range]),300)
        colormap=np.linspace(0,np.nanmax(sim.Urms[sim.exp_range]),300)
        plt.colorbar(axexp1.tricontourf(-sim.Yall[sim.exp_range],sim.Xall[sim.exp_range],sim.Urms[sim.exp_range],colormap,cmap='jet'),ax=axexp2,ticks=[0,0.1])
        plt.colorbar(axexp2.contourf(-exp.x,exp.y,exp.rmsv,colormap,cmap='jet',aspect='auto'),ax=axexp1,ticks=[0,0.1])
    # patch to axavg
    axavg.add_patch(patches.Rectangle(
        (-sim.ymin,sim.xmin),
        -sim.ywidth,
        sim.xwidth,
        fill=False,
        edgecolor='red',
        linewidth=2,
        ))

    axavg.set_xlabel(r'$y/D$')
    axavg.set_ylabel(r'$x/D$')
    #axavg.axis('equal')
    # stuff for exp axexp1
    axexp1.set_title('exp')
    axexp1.set_xlabel(r'$y/D$')
    axexp1.set_ylabel(r'$x/D$')
    #axexp1.axis('image')
    # stuff for sim axexp1
    axexp2.set_title('sim')
    axexp2.set_xlabel(r'$y/D$')
    axexp2.set_ylabel(r'$x/D$')
    #axexp2.axis('image')
    axexp2.axis([0,2.,14,15.5])
    # figure for centerline
    if direction=='u':
        lns1=axexp3.plot(sim.UCL,sim.XCL,'b-',label=r'sim $u$')
        lns2=axexp3.plot(exp.VCL.flatten(),exp.YCL.flatten(),'bo',label=r'exp $u$')
        axexp3.set_xlabel(r'$U_{CL}/U_{ref}$')
    elif direction=='v':
        lns1=axexp3.plot(sim.UrmsCL,sim.XCL,'b-',label=r'sim $u_{rms}$')
        lns2=axexp3.plot(exp.VrmsCL.flatten(),exp.YCL.flatten(),'bo',label=r'exp $u_{rms}$')
        axexp3.set_xlabel(r'$U_{rms}/U_{ref}$')
    lns=lns1+lns2#+lns3+lns4
    labs=[l.get_label() for l in lns]
    axexp3.set_ylabel(r'$x/D$')
    axexp3.legend(loc='upper right',framealpha=0.75,numpoints=1)
    fig.tight_layout()
    # arrow path for exp
    xyA=(-sim.ymax,sim.xmin)
    xyB=(-sim.ymin,sim.xmax)
    axexp1.add_artist(patches.ConnectionPatch(
        xyA=xyA,
        xyB=xyB,
        coordsA="data",
        coordsB="data",
        axesA=axexp1,
        axesB=axavg,
        linewidth=2,
        ))
    # arrow path for sim
    xyA=(-sim.ymax,sim.xmax)
    xyB=(-sim.ymin,sim.xmin)
    axexp2.add_artist(patches.ConnectionPatch(
        xyA=xyA,
        xyB=xyB,
        coordsA="data",
        coordsB="data",
        axesA=axexp2,
        axesB=axavg,
        linewidth=2,
        ))
    plt.savefig(sim.save_directory+'/ave__'+direction+'.png',bbox_inches='tight')

# average simulation data
#csv_filename='/home/shaun/Desktop/DA/Create_Plots/From_Jeff_new_BC/Data/data_slice_000_z.csv'
#save_directory='/home/shaun/Desktop/DA/Create_Plots/From_Jeff_new_BC/PostProcess'
csv_filename='../Create_Plots/From_Jeff_new_BC/Data2/Instant_and_last_ave/Last_ave.csv'
save_directory='../Create_Plots/From_Jeff_new_BC/Data2/Instant_and_last_ave'
s_avg=sim_avg(csv_filename,save_directory)

## read csv create npy files ( runs once )
#s_avg.csv_slice_to_npy()

## load npy
s_avg.load()

# average experimental data
ave_v_npy='../DA_ExperimentalData/From_JFrank_5000_cold_SandiaC/Ave_rms/ave_v.npy'
ave_u_npy='../DA_ExperimentalData/From_JFrank_5000_cold_SandiaC/Ave_rms/ave_u.npy'
rms_v_npy='../DA_ExperimentalData/From_JFrank_5000_cold_SandiaC/Ave_rms/RMSv.npy'
rms_u_npy='../DA_ExperimentalData/From_JFrank_5000_cold_SandiaC/Ave_rms/RMSu.npy'
mat_filename = '../DA_ExperimentalData/From_JFrank_5000_cold_SandiaC/T093013_01_Velocity_1-500.mat'
e_avg=exp_avg(ave_u_npy,ave_v_npy,rms_u_npy,rms_v_npy,mat_filename)
## load npy
e_avg.load()
e_avg.center_line()

# plot
#plt.style.use('classic')
#plt.style.use('ggplot')
#plt.style.use('seaborn-paper')
#plt.style.use('seaborn-notebook')
#plot_data_u(s_avg,e_avg)
plot_data(s_avg,e_avg,'u')
plot_data(s_avg,e_avg,'v')
#plt.savefig('Exp_vs_Sim.png',transparent=True,bbox_inches='tight')
#plot_data_v(s_avg,e_avg)

#plt.show()
