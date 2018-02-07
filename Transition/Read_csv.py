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
        # set flow density and u_inf freestream velocity
        self.rho=1.2 # kg/m^3
        self.uinf=50.1 #m/s
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
        self.x=d['X']
        self.y=d['Y']
        self.z=d['Z']
        self.mu=d['mu']
        self.dudz=d['dudz']
        try:
            d['tau_wall']
        except: 
            print 'tau_wall variable does not exist, calculating from mu and dudz'
        else: 
            print 'tau_wall is defined'
            self.tau_w=d['tau_wall']
    def read_exp(self,filename,names=None):
        d=np.genfromtxt(filename,delimiter=',',names=names)
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

    def plot_csv(self,ax):
        """
        Plots the first slice of u in the x, z plane
        """
        #ax.tricontour (    self.x,self.z,self.dudz,30,cmap='jet')
        plt.colorbar(
            ax.tricontourf(self.x,self.z,self.dudz,300,cmap='jet'))
        ax.axis('equal')
        ax.set_title('dudz')
        #plt.xlabel('X',**self.labelfont)
        #plt.ylabel('Y',**self.labelfont)
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.tight_layout()

    def plot_dudz(self,ax):
        """
        Plots the first slice of dudz in the x axis on the plate
        """
        #ax.tricontour (    self.x,self.z,self.dudz,30,cmap='jet')
        ax.set_title('dudz')
        x=self.z==0
        print x
        ax.plot(self.x[x],self.dudz[x],'.')
        #plt.xlabel('X',**self.labelfont)
        #plt.ylabel('Y',**self.labelfont)
        plt.xlabel('X')
        plt.ylabel('dudz')
        plt.tight_layout()

    def plot_Cf(self,ax,label='CFX'):
        """
        Plots skin Cf on the x axis on the plate
        """
        #ax.tricontour (    self.x,self.z,self.dudz,30,cmap='jet')
        ax.set_title(r'$C_f$')

        # plot CFX values
        #rho=1.1839 # kg/m^3

        x=self.z==0
        self.Rex=self.rho*self.uinf*self.x[x] / (np.mean(self.mu))
        try:
            self.tau_w
        except:
            print 'calculating tau_w from mu and dudz'
            self.tau_w = self.mu[x]*self.dudz[x]
        else:
            print 'extracting only points at wall for tau_wall'
            self.tau_w = self.tau_w[x]
        self.Cf = self.tau_w / (0.5 * self.rho * self.uinf**2)
        ax.plot(self.Rex,self.Cf,'*',label=label)



        plt.xlabel(r'$Re_x$')
        plt.ylabel(r'$C_f$')
        plt.tight_layout()

    def plot_turb_lam_Cf(self,ax):
        # plot turbulent (1/7 power law) and laminar equations
        x=np.linspace(0,2,1000)
        Rex=self.rho*self.uinf*x / (np.mean(self.mu))
        turb=0.0576 * Rex**(-1./5.)
        lam = 0.664/np.sqrt(Rex)
        ax.plot(Rex,turb,'--',label='turbulent')
        ax.plot(Rex,lam,'-',label='laminar')


if __name__=="__main__":
    # user options
    calc_options1={
            'sim_dir':'/home/shaun/Documents/Winter2018/ME392_Research/FlatPlate/NASA_Turb/FlatPlate_NASA/MultipleGrids/Structured/Coarse_done/',
            'filename':'export_vel_tau_shortened_to_be_read',
            #'sim_var':'dudz',
            #'X_interp_pts':120j,
            #'Y_interp_pts':160j,
            #'Z_interp_pts':26j,
            }
    # read in data from ANSYS sims
    sim1 = read_csv(calc_options1)
    sim1.read_csv()

    calc_options1['sim_dir']='/home/shaun/Documents/Winter2018/ME392_Research/FlatPlate/NASA_Turb/FlatPlate_NASA/MultipleGrids/Structured/Coarse/'
    calc_options1['filename']='export_shortened'
    sim2 = read_csv(calc_options1)
    sim2.read_csv()

    calc_options1['sim_dir']='/home/shaun/Documents/Winter2018/ME392_Research/FlatPlate/NASA_Turb/FlatPlate_NASA/MultipleGrids/Unstructured_Quad/Coarse_rightBCs/'
    calc_options1['filename']='export_shortened'
    sim3 = read_csv(calc_options1)
    sim3.read_csv()


    # open figure and set style
    plt.style.use('seaborn-paper')
    fig = plt.figure()
    ax = plt.subplot(111)

    # plot data from ANSYS
    #sim1.plot_csv(ax)
    #sim1.plot_dudz(ax)
    sim1.plot_Cf(ax,label='CFX structured wrong BCs')
    sim2.plot_Cf(ax,label='CFX structured right BCs')
    sim3.plot_Cf(ax,label='CFX unstructured right BCs')

    # plot turbulent and laminar equations
    sim2.plot_turb_lam_Cf(ax)

    # plot experimental and Menter2009 data
    menter2009SK = sim1.read_exp('/home/shaun/Documents/Winter2018/ME392_Research/FlatPlate/NASA_Turb/Exp_Previous_Paper_Data/SchubauerKlebanoff_FlatPlate.csv')
    #menter2009sim = sim1.read_exp('/home/shaun/Documents/Winter2018/ME392_Research/FlatPlate/Menter2009.csv')
    menter2009sim = sim1.read_exp('/home/shaun/Documents/Winter2018/ME392_Research/FlatPlate/NASA_Turb/Exp_Previous_Paper_Data/Menter2009_morepoints.csv')
    ax.plot(menter2009SK[:,0],menter2009SK[:,1],'s',label='S+K')
    ax.plot(menter2009sim[:,0],menter2009sim[:,1],'.',label='Menter 2009 sim')

    # add legend
    ax.legend(loc='best',numpoints=1,frameon=False)
    ax.axis([0,4500000,0,0.01])
    ax.set_title(r'$C_f$ for Schubauer and Klebanof plate (natural transition)')
    plt.tight_layout()

    plt.show()
