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
        d=np.genfromtxt(self.csv_filename,delimiter=',',skip_header=1)
        uexp,uexp1,u,v,w,x,y,z=[ d[:,i] for i in range(0,8) ]

        # non-dimensionalize everything
        D_ref = 0.00745
        U_ref = 27.5
        x = x/D_ref
        y = y/D_ref
        z = z/D_ref
        uexp=uexp/U_ref
        uexp1=uexp1/U_ref
        u=u/U_ref
        v=v/U_ref
        w=w/U_ref

        # interpolate to grid
        nx=3000j
        ny=260j
        self.Xall,self.Yall=np.mgrid[0:30:nx,-4:4:ny] # make grid of slice down the middle (if more data is given, can do slice of whole domain to get the whole experimental region)
        #self.XCL,self.YCL=np.mgrid[0:30:nx,0:1:1j] # make grid of slice down the middle (if more data is given, can do slice of whole domain to get the whole experimental region)
        points=np.array([x,y]).T
        self.Uexp=griddata(points,uexp,(self.Xall,self.Yall),method='linear')
        self.Uexp1=griddata(points,uexp1,(self.Xall,self.Yall),method='linear')
        self.U=griddata(points,u,(self.Xall,self.Yall),method='linear')
        #self.UCL=griddata(points,uavg,(self.XCL,self.YCL),method='linear')
        #self.VCL=griddata(points,vavg,(self.XCL,self.YCL),method='linear')

        # save interpolated data
        np.save(self.save_directory+'/uexp',self.Uexp)
        np.save(self.save_directory+'/uexp1',self.Uexp1)
        np.save(self.save_directory+'/u',self.U)
        #np.save(self.save_directory+'/ucl',self.UCL)
        #np.save(self.save_directory+'/vcl',self.VCL)
        np.save(self.save_directory+'/X',self.Xall)
        np.save(self.save_directory+'/Y',self.Yall)
        #np.save(self.save_directory+'/Xcl',self.XCL)
        #np.save(self.save_directory+'/Ycl',self.YCL)

    def load(self):
        # load data
        self.Uexp=np.load(self.save_directory+'/uexp.npy')
        self.Uexp1=np.load(self.save_directory+'/uexp1.npy')
        self.U=np.load(self.save_directory+'/u.npy')
        #self.UCL =np.load(self.save_directory+'/ucl.npy')
        #self.VCL =np.load(self.save_directory+'/vcl.npy')
        self.Xall=np.load(self.save_directory+'/X.npy')
        self.Yall=np.load(self.save_directory+'/Y.npy')
        #self.XCL=np.load(self.save_directory+'/Xcl.npy')
        #self.YCL=np.load(self.save_directory+'/Ycl.npy')
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

def plot_data(sim):
    fig = plt.figure(figsize=(3.0,5))
    #fig = plt.figure(figsize=(7.36,5))
    #axavg      =plt.subplot2grid((4,140),(0,0),rowspan=4,colspan=40)
    axavg=plt.subplot(111)
    colormap=np.linspace(np.nanmin(sim.U),np.nanmax(sim.U),300)
    plt.colorbar(axavg.contourf(-sim.Yall,sim.Xall,sim.U,colormap,cmap='jet'),ax=axavg)
    # patch to axavg
    axavg.add_patch(patches.Rectangle(
        (-sim.ymin,sim.xmin),
        -sim.ywidth,
        sim.xwidth,
        fill=False,
        edgecolor='red',
        linewidth=2,
        ))
    axavg.set_xlabel(r'$Y/D$')
    axavg.set_ylabel(r'$X/D$')
    axavg.axis('scaled')
    fig.tight_layout()
    plt.savefig(sim.save_directory+'/Instantaneous__u.png',bbox_inches='tight')

# average simulation data
csv_filename='/home/shaun/Desktop/DA/Create_Plots/DA_paper/Updated_BC/Reference/Instantaneous/Instant.csv'
save_directory='/home/shaun/Desktop/DA/Create_Plots/DA_paper/Updated_BC/Reference/Instantaneous'
s_avg=sim_avg(csv_filename,save_directory)

## read csv create npy files ( runs once )
#s_avg.csv_slice_to_npy()

## load npy
s_avg.load()

# plot
plot_data(s_avg)

plt.show()
