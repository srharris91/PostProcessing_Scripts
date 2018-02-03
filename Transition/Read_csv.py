import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

class read_csv:
    """read_csv class of functions
    
    Bunch of functions to be used with postprocessing simulation data from CFX output (csv exported data)
    """
    def __init__(self,calc_options=None):
        """Initialize class

        Parameters:
        calc_options -- Dictionary input containing several of input values [default = None]

        Info:
        If calc_options is None, then will call a bunch of default values listed in this function
        Otherwise, will use set_values(**calc_options) function to set the desired values
        """
        if calc_options==None:
            self.sim_dir='./'
            self.title='X'
            self.filename='OI_C'
            self.csv=self.filename+'.csv'
            self.sim_var='U'
            self.X_interp_pts=126j
            self.Y_interp_pts=169j
            self.Z_interp_pts=26j
            self.loaded=False
            print 'set default values'
        else:
            self.set_values(**calc_options)
    def set_values(self,sim_dir='./',title='X',filename='OI_C',sim_var='U',X_interp_pts=126j,Y_interp_pts=169j,Z_interp_pts=26j,loaded=False):
        """Initialize equation using calc_options dictionary input, if not specified, then the defaults are shown below

            self.sim_dir='./'
            self.title='X'
            self.filename='OI_C'
            self.csv=self.filename+'.csv'
            self.sim_var='U'
            self.X_interp_pts=126j
            self.Y_interp_pts=169j
            self.Z_interp_pts=26j
            self.loaded=False

        """
        self.sim_dir=sim_dir
        self.title=title
        self.filename=filename
        self.csv=self.filename+'.csv'
        self.sim_var=sim_var
        self.X_interp_pts=X_interp_pts
        self.Y_interp_pts=Y_interp_pts
        self.Z_interp_pts=Z_interp_pts
        self.loaded=loaded

    def read_files_list(self,startswith):
        """Return list of files that startswith

        Parameters:
        self -- sim class of functions
        startswith -- str to be matched with files in self.sim_dir location

        Returns:
        list of files in location of self.sim_dir
        """
        # read in files
        read_list=[]
        for file in os.listdir(self.sim_dir):
            if file.startswith(startswith):
                read_list.append(self.sim_dir+file)
        def atoi(text):
            return int(text) if text.isdigit() else text
        def natural_keys(text):
            '''
            alist.sort(key=natural_keys) sorts in human order
            http://nedbatchelder.com/blog/200712/human_sorting.html
            (See Toothy's implementation in the comments)
            '''
            return [ atoi(c) for c in re.split('(\d+)', text) ]
        read_list.sort(key=natural_keys)
        return read_list

    def read_csv(self):
        """ 
        Read in csv file using sim_dir and filename
        """
        
        d=np.genfromtxt(self.sim_dir+self.csv,delimiter=',',names=True)
        return d

    def read_csv_files(self):
        """
        Read list of files that matches self.csv pattern
        WARNING: creating time axis from user inputs from CharlesX.input file.  May need to edit later.
        """
        self.csv_files=self.read_files_list(self.csv)
        xl,yl,zl,ul=[],[],[],[]
        for filenamei in self.csv_files:
            print '  simulation reading '+filenamei
            xi,yi,zi,ui=self.read_csv(filenamei)
            xl.append(xi)
            yl.append(yi)
            zl.append(zi)
            ul.append(ui)
        self.xcsv=np.array(xl)
        self.ycsv=np.array(yl)
        self.zcsv=np.array(zl)
        self.ucsv=np.array(ul)
    def interpolate_csv_data(self):
        """
        Take xcsv,ycsv,zcsv,ucsv variables and interpolate them onto a grid using griddata
        and creates self.Xall, self.Yall, self.Zall, and self.Uall
        """
        x=self.xcsv
        y=self.ycsv
        z=self.zcsv
        u=self.ucsv
        #self.Xall,self.Yall,self.Zall=np.mgrid[x.min():x.max():126j,y.min():y.max():169j,z.min():z.max():26j] # make grid of slice down the middle (if more data is given, can do slice of whole domain to get the whole experimental region)
        print 'shapes',x.shape,y.shape,z.shape
        print 'mins',x[0].min(),x[0].max(),y[0].min(),y[0].max(),z[0].min(),z[0].max()
        self.Xall,self.Yall,self.Zall=np.mgrid[x[0].min():x[0].max():self.X_interp_pts,y[0].min():y[0].max():self.Y_interp_pts,z[0].min():z[0].max():self.Z_interp_pts] # make grid of slice down the middle (if more data is given, can do slice of whole domain to get the whole experimental region)
        Uall=[]
        for i,ui in enumerate(u):
            points=np.array([x[i],y[i],z[i]]).T
        #self.Uall=np.stack(np.array([ griddata(points,u[i],(self.Xall,self.Yall,self.Zall),method='linear') for i in np.arange(u.shape[0])]),axis=0)
            Uall.append(griddata(points,ui,(self.Xall,self.Yall,self.Zall),method='linear'))
        self.Uall=np.array(Uall)
        print 'Uall.shape=',self.Uall.shape
    def save_all_data_npy(self):
        """
        save all interpolated data to npy files
        {self.sim_dir}Xall.npy saved self.Xall
        {self.sim_dir}Yall.npy saved self.Yall
        {self.sim_dir}Zall.npy saved self.Zall
        {self.sim_dir+self.sim_var}all_{self.direction}.npy saved self.Uall
        """
        np.save(self.sim_dir+self.sim_var+'all_'+self.direction,self.Uall)
        np.save(self.sim_dir+'Xall',self.Xall)
        np.save(self.sim_dir+'Yall',self.Yall)
        np.save(self.sim_dir+'Zall',self.Zall)
    def load_all_data_npy(self):
        """
        Load all the interpolated data from npy files
        {self.sim_dir}Xall.npy saved self.Xall
        {self.sim_dir}Yall.npy saved self.Yall
        {self.sim_dir}Zall.npy saved self.Zall
        {self.sim_dir+self.sim_var}all_{self.direction}.npy saved self.Uall
        """
        self.Uall=np.load(self.sim_dir+self.sim_var+'all_'+self.direction+'.npy')
        self.Xall=np.load(self.sim_dir+'Xall.npy')
        self.Yall=np.load(self.sim_dir+'Yall.npy')
        self.Zall=np.load(self.sim_dir+'Zall.npy')
        self.t=np.arange(self.Uall.shape[0]) * 1.E-7 * 1000

    def plot_csv(self,colormap=True):
        """
        Plots the first slice of Uall in the Xall, Yall plane
        """
        fig = plt.figure()
        ax = plt.subplot(111)
        if colormap:
            if self.direction=='u':
                vmin,vmid,vmax = 0.,0.35,0.7
            elif self.direction=='v':
                vmin,vmid,vmax=-0.05,0.,0.05
            elif self.direction=='w':
                vmin,vmid,vmax = -0.05,0.,0.05
            colorvalues = np.linspace(vmin,vmax,100)
            ax.contour (self.Xall[:,:,0],self.Yall[:,:,0],self.Uall[0,:,:,0],colorvalues,cmap='jet',vmin=vmin,vmax=vmax)
            plt.colorbar(
            ax.contourf(self.Xall[:,:,0],self.Yall[:,:,0],self.Uall[0,:,:,0],colorvalues,cmap='jet',vmin=vmin,vmax=vmax)
            ,ticks=[vmin,vmid,vmax])
        else:
            ax.contour (self.Xall[:,:,0],self.Yall[:,:,0],self.Uall[0,:,:,0],30,cmap='jet')
            plt.colorbar(
            ax.contourf(self.Xall[:,:,0],self.Yall[:,:,0],self.Uall[0,:,:,0],30,cmap='jet')
            ,)
        ax.axis('equal')
        ax.set_title('instantaneous snapshot')
        plt.xlabel('X',**self.labelfont)
        plt.ylabel('Y',**self.labelfont)
        plt.tight_layout()

    def plot_rms(self):
        """
        Plots an average of simulation data along with the rms values of the first time z plane
        """
        fig = plt.figure()
        ax = plt.subplot(111)
        average = np.nanmean(self.Uall,axis=self.homogeneous)[np.newaxis,:,:,np.newaxis]
        ax.contour                  (self.Xall[:,:,0],self.Yall[:,:,0],average[0,:,:,0],30,cmap='jet')
        plt.colorbar( ax.contourf   (self.Xall[:,:,0],self.Yall[:,:,0],average[0,:,:,0],30,cmap='jet') ,)
        ax.axis('equal')
        ax.set_title('average')
        ax.set_xlabel('X',**self.labelfont)
        ax.set_ylabel('Y',**self.labelfont)
        plt.tight_layout()

        fig = plt.figure()
        ax = plt.subplot(111)
        u_prime = self.Uall-average
        urms = np.nanmean((u_prime**2),axis=self.homogeneous)
        ax.contour                  (self.Xall[:,:,0],self.Yall[:,:,0],urms,30,cmap='jet')
        plt.colorbar( ax.contourf   (self.Xall[:,:,0],self.Yall[:,:,0],urms,30,cmap='jet') ,)
        ax.axis('equal')
        ax.set_title('urms')
        ax.set_xlabel('X',**self.labelfont)
        ax.set_ylabel('Y',**self.labelfont)
        plt.tight_layout()




if __name__=="__main__":
    # user options
    calc_options1={
            'sim_dir':'/home/shaun/Documents/Winter2018/ME392_Research/FlatPlate/NASA_Turb/FlatPlate_NASA/MultipleGrids/Structured/Coarse_done/',
            'filename':'export_vel_tau_shortened_to_be_read',
            'sim_var':'Velocity u.Gradient Z [ s^-1 ]',
            'X_interp_pts':120j,
            'Y_interp_pts':160j,
            'Z_interp_pts':26j,
            }
    sim1 = read_csv(calc_options1)

    fig = plt.figure()
    ax = plt.subplot(111)
    data=sim1.read_csv()
    print data
    # for R11(tau)
    #sim1.R11_tau_calc_plot(ax)
    #sim1.plot_axis_labels_tau(ax)
    # for R11(r)
    plt.show()
