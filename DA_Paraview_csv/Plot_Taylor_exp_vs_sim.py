import numpy as np
import matplotlib.pyplot as plt

class base_class:
    def __init__(self,directory='./',fname='R11_output/R11.txt',cmap_range=100,cmap_color='jet'): # default
    # def __init__(self,directory='./',fname='R11_output/R11.txt',cmap_range=np.linspace(0,1.22,100),cmap_color='jet'): # Integral Length Scale
    # def __init__(self,directory='./',fname='R11_output/R11.txt',cmap_range=np.linspace(0,0.27,100),cmap_color='jet'): # Taylor length scale
        self.directory=directory
        self.fname=fname
        self.cmap_range=cmap_range
        self.cmap_color=cmap_color
    def plot_contourL(self,ax,title):
        cs1=ax.contour (self.x,self.y,self.L,self.cmap_range,cmap=self.cmap_color,extend='both')
        cs2=ax.contourf(self.x,self.y,self.L,self.cmap_range,cmap=self.cmap_color,extend='both')
        cs2.cmap.set_under('k') 
        plt.colorbar(cs2)
        ax.set_xlabel(r'$x/D$')
        ax.set_ylabel(r'$y/D$')
        ax.set_title(title+' Integral Length Scale')
        ax.axis('equal')
    def plot_contourTaylorL(self,ax,title):
        cs1=ax.contour (self.x,self.y,self.TaylorL,self.cmap_range,cmap=self.cmap_color,extend='both')
        cs2=ax.contourf(self.x,self.y,self.TaylorL,self.cmap_range,cmap=self.cmap_color,extend='both')
        cs2.cmap.set_under('k') 
        plt.colorbar(cs2)
        ax.set_xlabel(r'$x/D$')
        ax.set_ylabel(r'$y/D$')
        ax.set_title(title+' Taylor Length Scale')
        ax.axis('equal')
    def plot_points(self,points_list,list_ax_contours,ax_R11):
        for p in points_list:
            self.load_point_data(p)
            ax_R11.plot(self.y[self.xi,:]-self.y[self.xi,self.yi],self.R11,label=self.fnamept)
            for ax in list_ax_contours:
                ax.plot(self.x[self.xi,self.yi],self.y[self.xi,self.yi],'ko')

class sim(base_class):
    def read_data(self):
        print 'reading data from ',self.directory+self.fname
        XYLlam=np.genfromtxt(self.directory+self.fname)
        self.y=XYLlam[:,0].reshape(126,169).T
        self.x=XYLlam[:,1].reshape(126,169).T
        self.L=XYLlam[:,2].reshape(126,169).T
        self.TaylorL=XYLlam[:,3].reshape(126,169).T
    def load_point_data(self,p):
        self.xi,self.yi=p
        self.fnamept=self.fname[:-7]+'XYR'+self.fname[1:3]+'/ijXYR'+self.fname[1:3]+'_%i_%i_%.4f_%.4f'%(self.yi,self.xi,self.y[self.xi,self.yi],self.x[self.xi,self.yi])+'.npy'
        print 'loading ',self.fnamept
        self.R11=np.load(self.directory+self.fnamept)

class exp(base_class):
    def read_data(self):
        print 'reading data from ',self.directory+self.fname
        XYLlam=np.genfromtxt(self.directory+self.fname)
        self.x=XYLlam[:,0].reshape(169,126)
        self.y=XYLlam[:,1].reshape(169,126)
        self.L=XYLlam[:,2].reshape(169,126)
        self.TaylorL=XYLlam[:,3].reshape(169,126)
    def load_point_data(self,p):
        self.xi,self.yi=p
        self.fnamept=self.fname[:-7]+'XYR'+self.fname[1:3]+'/ijXYR'+self.fname[1:3]+'_%i_%i_%.4f_%.4f'%(self.xi,self.yi,self.x[self.xi,self.yi],self.y[self.xi,self.yi])+'.npy'
        print 'loading ',self.fnamept
        self.R11=np.load(self.directory+self.fnamept)

if __name__=="__main__":
    sim=sim(directory='../Create_Plots/OI_taylorLength/',fname='R11_output/R11.txt')
    sim.read_data()
    exp=exp(directory='../DA_ExperimentalData/From_JFrank_5000_cold_SandiaC/R11_output/LengthScales_Taylor/',fname='R22_output/R22.txt')
    exp.read_data()

    ##### Integral Length Scales
    # set colormap
    exp.cmap_range=np.linspace(0,1.22,100)
    sim.cmap_range=np.linspace(0,1.22,100)
    fig=plt.figure(figsize=(12,4))
    axexpL=plt.subplot(121)
    exp.plot_contourL(axexpL,'exp')
    axsimL=plt.subplot(122)
    sim.plot_contourL(axsimL,'sim')
    plt.tight_layout()

    ##### Taylor Length Scales
    # set colormap
    exp.cmap_range=np.linspace(0.1,0.4,100)
    sim.cmap_range=np.linspace(0.1,0.4,100)
    fig=plt.figure(figsize=(12,4))
    axexpTay=plt.subplot(121)
    exp.plot_contourTaylorL(axexpTay,'exp')
    axsimTay=plt.subplot(122)
    sim.plot_contourTaylorL(axsimTay,'sim')
    plt.tight_layout()

    # load pt
    p=[
            (47,48),        # outer middle
            (47,0),         # outer inlet
            (47,125),       # outer outlet
            #(106,48),      # middle middle
            #(166,48),      # inner middle
            #(166,0),       # inner inlet
            #(166,49),      # inner middle (nan taylor length)
            ]
    list_ax=[axexpL,axsimL,axexpTay,axsimTay]
    fig=plt.figure()
    ax=plt.subplot(111)
    exp.plot_points(p,list_ax,ax)
    sim.plot_points(p,list_ax,ax)
    ax.legend(loc='best',numpoints=1)

    # show plot
    plt.show()
