import time
import re
import os
import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy.interpolate import griddata
from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

class sim:
    def __init__(self,calc_options=None):
        if calc_options==None:
            self.sim_dir =  './'
            self.direction = 'u'
            self.manufactured = False
            self.xi = 0
            self.yi = 0
            self.method = 1
            self.option = 2
            self.t,self.X,self.Y,self.Z,self.U = 0,0,0,0,0
            self.homogeneous = (0,-1) # time and z
            self.labelfont = { 'fontsize':14, 'fontname':'Times New Roman'}
            self.RGB=[0., 0.0, 0.0, 1.0, 20., 1.0, 0.0, 0.0]
            self.title='X'
            self.vtu='OI_C'
            self.csv=self.vtu+'_csv'
            self.read_from='box_probe'
            self.sim_var='U'
            self.X_interp_pts=126j
            self.Y_interp_pts=169j
            self.Z_interp_pts=26j
            self.loaded=False
            print 'set default values'
        else:
            self.set_values(**calc_options)
    def set_values(self, sim_dir='./',direction='u',manufactured=False,xi=0,yi=0,method=1,option=2,t=0,X=0,Y=0,Z=0,U=0,homogeneous=(0,-1),labelfont={'fontsize':14,'fontname':'Times New Roman'},RGB=[0., 0.0, 0.0, 1.0, 20., 1.0, 0.0, 0.0],vtu='OI_C',read_from='box_probe',sim_var='U',X_interp_pts=6j,Y_interp_pts=9j,Z_interp_pts=1j):
            self.sim_dir = sim_dir
            self.direction = direction
            self.manufactured = manufactured
            self.xi = xi
            self.yi = yi
            self.method = method
            self.option = option
            self.t,self.X,self.Y,self.Z,self.U = t,X,Y,Z,U
            self.homogeneous=homogeneous
            self.labelfont=labelfont
            if direction=='u':
                self.title='X'
            elif direction=='v':
                self.title='Y'
            elif direction=='w':
                self.title='Z'
            self.RGB=RGB
            self.vtu=vtu
            self.csv=self.vtu+'_csv'
            self.read_from=read_from
            self.sim_var=sim_var
            self.X_interp_pts=X_interp_pts
            self.Y_interp_pts=Y_interp_pts
            self.Z_interp_pts=Z_interp_pts
            self.loaded=False

    def read_files_list(self,startswith):
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
    def slice_orient_save_vtu_to_csv(self):
        # read in data
        U=self.sim_var
        self.vtu_list=self.read_files_list(self.vtu)
        a1D_060000_vtu = XMLUnstructuredGridReader( FileName=self.vtu_list)
        a1D_060000_vtu.PointArrayStatus = [U]
        # CLIP data range
        #xmin = 0.1041  
        #xmax = 0.116
        xmin = 0.103713111302   # true exp region axial
        xmax = 0.11600928642    # true exp region axial
        #xmax = 0.14
        #ymin = -0.016
        ymin = -0.0163191665626 # true exp region radial
        #ymax = 0.0001
        ymax = 0.00020689279668 # true exp region radial
        zmin = -0.0012  # true exp region azimuthal
        zmax = 0.0012   # true exp region azimuthal
        # clip xmin
        Clip1 = Clip( ClipType="Plane" )
        Clip1.Scalars = ['POINTS', U]
        Clip1.ClipType.Origin = [xmin, 0.0, 0.0]
        Clip1.ClipType.Normal = [ 1.0, 0.0, 0.0]
        # clip xmax
        Clip2 = Clip( ClipType="Plane" )
        Clip2.Scalars = ['POINTS', U]
        Clip2.ClipType.Origin = [xmax, 0.0, 0.0]
        Clip2.ClipType.Normal = [-1.0, 0.0, 0.0]
        # clip ylow
        Clip3 = Clip( ClipType="Plane" )
        Clip3.Scalars = ['POINTS', U]
        Clip3.ClipType.Origin = [0.0, ymin, 0.0]
        Clip3.ClipType.Normal = [0.0,  1.0, 0.0]
        # clip yup
        Clip4 = Clip( ClipType="Plane" )
        Clip4.Scalars = ['POINTS', U]
        Clip4.ClipType.Origin = [0.0, ymax, 0.0]
        Clip4.ClipType.Normal = [0.0, -1.0, 0.0]
        # clip zlow
        Clip5 = Clip( ClipType="Plane" )
        Clip5.Scalars = ['POINTS', U]
        Clip5.ClipType.Origin = [0.0, 0.0, zmin]
        Clip5.ClipType.Normal = [0.0, 0.0, 1.0]
        # clip yup
        Clip6 = Clip( ClipType="Plane" )
        Clip6.Scalars = ['POINTS', U]
        Clip6.ClipType.Origin = [0.0, 0.0, zmax]
        Clip6.ClipType.Normal = [0.0, 0.0, -1.0]

        # slice
        #Slice1 = Slice( SliceType="Plane" )
        #Slice1.SliceOffsetValues = [0.0]
        #Slice1.SliceType.Origin = [0.0115, 0.0, 0.0]
        #Slice1.SliceType.Normal = [0.0, 0.0, 1.0]

        RenderView1 = GetRenderView()
        DataRepresentation1 = Show()

        writer=CreateWriter(self.sim_dir+self.csv+'.csv')
        writer.WriteAllTimeSteps=1
        writer.FieldAssociation="Points"
        writer.UpdatePipeline()
        del writer
        Delete(RenderView1)
        Delete(DataRepresentation1)
        del RenderView1
        del DataRepresentation1
    def read_csv(self,filename):
        d=np.genfromtxt(filename,delimiter=',',skip_header=1)
        if self.direction=='u':
            u,x,y,z=d[:,-6],d[:,-3],d[:,-2],d[:,-1]
        elif self.direction=='v':
            u,x,y,z=d[:,-5],d[:,-3],d[:,-2],d[:,-1]
        elif self.direction=='w':
            u,x,y,z=d[:,-4],d[:,-3],d[:,-2],d[:,-1]
        # non-dimensionalize everything
        D_ref = 0.00745
        U_ref = 27.5
        x = x/D_ref
        y = y/D_ref
        z = z/D_ref
        u = u/U_ref
        return [x,y,z,u]
    def read_csv_files(self):
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
        # create data for t ( need to know dt and output interval in input files to charlesX)
        # user input from CharlesX.input
        dt=0.5E-7
        interval=10000
        # calculate and create dt from user input
        print 'creating time axis from user inputs from CharlesX.input file.  May need to edit later.'
        D_ref = 0.00745
        U_ref = 27.5
        dt_sim=dt*interval
        t=np.arange(self.ucsv.shape[0])*dt_sim
        self.t = t/(D_ref/U_ref)
    def interpolate_csv_data(self):
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
        np.save(self.sim_dir+self.sim_var+'all_'+self.direction,self.Uall)
        np.save(self.sim_dir+'Xall',self.Xall)
        np.save(self.sim_dir+'Yall',self.Yall)
        np.save(self.sim_dir+'Zall',self.Zall)
    def load_all_data_npy(self):
        self.Uall=np.load(self.sim_dir+self.sim_var+'all_'+self.direction+'.npy')
        self.Xall=np.load(self.sim_dir+'Xall.npy')
        self.Yall=np.load(self.sim_dir+'Yall.npy')
        self.Zall=np.load(self.sim_dir+'Zall.npy')
        self.t=np.arange(self.Uall.shape[0]) * 1.E-7 * 1000

    def read_probe(self):
        if self.direction=='u':
            ux='U-X'
        elif self.direction=='v':
            ux='U-Y'
        elif self.direction=='w':
            ux='U-Z'
        else:
            print 'Problem with direction in self.read_probe, direction=',self.direction
        fn=[
        self.sim_dir+'box_probe.'+ux,
        self.sim_dir+'box_probe.X',
        self.sim_dir+'box_probe.Y',
        self.sim_dir+'box_probe.Z',
        ]
        box_probe = {}
        for fi in fn:
            print '  simulation reading '+fi.split(".")[-1]+' from '+self.sim_dir
            box_probe[fi.split(".")[-1]] = np.genfromtxt(fi,delimiter=' ')
        u = box_probe[ux][:,3:]
        x = box_probe['X'][3:]
        y = box_probe['Y'][3:]
        z = box_probe['Z'][3:]
        t = box_probe[ux][:,1]
        X = np.swapaxes(x.reshape((2,10,20)),0,2)
        Y = np.swapaxes(y.reshape((2,10,20)),0,2)
        Z = np.swapaxes(z.reshape((2,10,20)),0,2)
        U = np.swapaxes(u.reshape((len(t),2,10,20)),1,3)
        
        # non-dimensionalize everything
        D_ref = 0.00745
        U_ref = 27.5
        self.Xall = X/D_ref
        self.Yall = Y/D_ref
        self.Zall = Z/D_ref
        self.Uall = U/U_ref
        self.t = t/(D_ref/U_ref)
        # save specific pt to class variables

    def plot_csv(self,colormap=True):
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



    ###################################################################
    # Integral length scale functions
    def extract_row_from_all(self):
        self.X = self.Xall[:,self.yi,0]
        self.Y = self.Yall[0,self.yi,0]
        self.Z = self.Zall[0,self.yi,0]
        self.U = self.Uall[:,:,self.yi,:]

    def set_manufactured_r(self):
        self.r=self.X - self.X[self.xi]
        self.u=np.empty((30,)+(len(self.r),)+(20,)) # t,y,z
        if self.option==1:
            self.u = (np.cos(self.r)*np.exp(-1.*self.r))[np.newaxis,:,np.newaxis] + 0.08*(2.*np.random.random(self.u.shape)-1.0)
        elif self.option==2:
            self.u = np.cos(self.r).reshape((len(self.r),1)) + 0.08*(2.*np.random.random((len(self.r),)+self.u.shape[1:])-1.0)

    def set_uprime_r(self):
        # calc u' v' and w'  velocity prime values
        u_avg = np.nanmean(self.U,axis=self.homogeneous)[np.newaxis,:,np.newaxis]
        self.u = self.U-u_avg

    def calc_R11_r_sim(self):
        nx = int(self.u.shape[1]) # axial domain
        self.R11 = np.zeros(nx)
        if self.method==0:
            r = np.arange(0,nx)-self.xi
            for ri in r:
                self.R11[self.xi+ri] = np.nanmean((self.u[:,self.xi,:]*self.u[:,self.xi+ri,:]),axis=self.homogeneous)
        elif self.method==1:
            r = np.arange(0,nx)-self.xi
            #self.R11[0] =  (self.u*self.u).nanmean()
            for ri in r: 
                self.R11[self.xi+ri] = np.nanmean((self.u[:,self.xi,:]*self.u[:,self.xi+ri,:]),axis=self.homogeneous)
        elif self.method==2:
            print 'method 2 not implemented yet'
            r = np.arange(0,nx)
            for rispace in r:
                for rime in r:
                    self.R11[rispace] = np.nanmean(self.R11[rispace]+(self.u[rime]*self.u[(rime+rispace)]))
            self.R11 = self.R11/float(nx)
        else:
            print 'get right method...',self.method
        self.r = self.X-self.X[self.xi]
        #self.R11_norm =  self.R11/self.R11[0]
        self.R11_norm =  self.R11/self.R11[self.xi]

    def set_plot_label_r(self):
        print 'option=',self.option
        self.label='sim,%s,x=%.4f,y=%.4f,m=%i,o=%i'%(
                self.direction,
                self.X[self.xi] ,
                self.Y ,
                self.method ,
                self.option ,)
    def plot_R11_r_sim(self,ax,label='',**kwargs):
        self.set_plot_label_r()
        ax.plot(self.r,self.R11_norm,label=self.label+','+label,**kwargs)
    def plot_axis_labels_r(self,ax):
        ax.set_xlabel(r'$r/D$',**self.labelfont)
        if self.direction=='u':
            ax.set_ylabel(r'$R_{11}(r)/R_{11}(0)$',**self.labelfont)
        elif self.direction=='v':
            ax.set_ylabel(r'$R_{22}(r)/R_{22}(0)$',**self.labelfont)
        elif self.direction=='w':
            ax.set_ylabel(r'$R_{33}(r)/R_{33}(0)$',**self.labelfont)
        ax.legend(loc='best',numpoints=1)
        plt.tight_layout()

    def calc_integral_length_L(self):
        # calc L integral scale
        self.L = np.trapz(self.R11_norm,x=self.r)
        print 'L= ',self.L
    def calc_taylor_length_lambda(self,ax=None):
        r=self.r
        R11=self.R11_norm

        dR=np.gradient(R11,edge_order=self.option)
        d2R=np.gradient(dR,edge_order=self.option)
        dr=np.gradient(r,edge_order=self.option)
        d2Rdr2=d2R/dr**2
        f_prime_prime=d2Rdr2[self.xi]

        #print 'f_prime_prime = ',f_prime_prime,self.xi,r[self.xi]
        self.lambdaf=(-f_prime_prime/2.)**-0.5+r[self.xi] # r[self.xi] should equal 0

        #  self.set_plot_label_r()
        #  if ax!=None:
        #      P=1.+0.5*f_prime_prime*(r-r[self.xi])**2
        #      ax.plot(r,P,'s',label='taylor '+self.label)
        #      ax.plot(self.lambdaf,0.,'o',label=r'$\lambda_f$')
        print 'TaylorL=',self.lambdaf

    def Manufactured_R11_r(self):
        L = self.r[-1]-self.r[0]
        # option 1
        if self.option==1:
            R11 = (
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
            R11 = (1./L)*(
                    np.cos(self.r)
                    * (( np.sin(L)      *np.cos(L)+L )  /   2.)
                    - ((  np.sin(self.r)   *np.sin(L)**2)  /   2.))
        R11_norm = R11/R11[0]
        return R11_norm

    def R11_r_calc_plot(self,ax,label=''):
        # calculate and plot R11(r)
        # for manufactured solution only
        if self.manufactured:
            self.X = np.linspace(0,0.05/(0.00745/27.5),5000.)
            self.set_manufactured_r()
        else:
            if self.read_from=='box_probe':
                self.read_probe()
            elif self.read_from=='vtu':
                if os.path.isfile(self.sim_dir+'Uall_'+self.direction+'.npy')==False:
                    print 'vtu to csv'
                    self.slice_orient_save_vtu_to_csv()
                    print 'read csv'
                    self.read_csv_files()
                    print 'interpolate'
                    self.interpolate_csv_data()
                    print 'Saving Uall, Xall, Yall, Zall npy files'
                    self.save_all_data_npy()
                else:
                    self.load_all_data_npy()
            self.extract_row_from_all()
            self.set_uprime_r()

        self.calc_R11_r_sim()
        self.calc_integral_length_L()
        self.calc_taylor_length_lambda(ax)

        self.plot_R11_r_sim(ax,marker='+',linewidth=0,label=label)

        if self.manufactured:
            R11simMNF = self.Manufactured_R11_r()
            self.set_plot_label_r()
            ax.plot(self.r,R11simMNF,marker='x',linewidth=0,label=self.label+' analytic')
    def save_R11_r(self,f_handle,):
        np.savetxt(f_handle,[[self.X[self.xi],self.Y,self.L,self.lambdaf]],delimiter=' ')
        if self.direction=='u':
            np.save('R11_output/XYR11/ijXYR11_%i_%i_%.4f_%.4f'%(self.xi,self.yi,self.X[self.xi],self.Y),self.R11_norm)
        elif self.direction=='v':
            #np.save('R22_output/XYR22/ijXYR22_'+str(self.xi)+'_%.4f'%(self.X),self.R11_norm)
            np.save('R22_output/XYR22/ijXYR22_%i_%i_%.4f_%.4f'%(self.xi,self.yi,self.X[self.xi],self.Y),self.R11_norm)
        elif self.direction=='w':
            #np.save('R33_output/XYR33/ijXYR33_'+str(self.xi)+'_%.4f'%(self.X),self.R11_norm)
            np.save('R33_output/XYR33/ijXYR33_%i_%i_%.4f_%.4f'%(self.xi,self.yi,self.X[self.xi],self.Y),self.R11_norm)
        print 'saved ijxy=',self.xi,self.yi,self.X[self.xi],self.Y
    def R11_r_calc_save(self,f_handle,label=''):
        # calculate and save R11(r) and length scales
        # for manufactured solution only
        if self.manufactured:
            self.X = np.linspace(0,0.05/(0.00745/27.5),5000.)
            self.set_manufactured_r()
        else:
            if self.read_from=='box_probe':
                self.read_probe()
            elif self.read_from=='vtu':
                if os.path.isfile(self.sim_dir+'Uall_'+self.direction+'.npy')==False:
                    print 'vtu to csv'
                    self.slice_orient_save_vtu_to_csv()
                    print 'read csv'
                    self.read_csv_files()
                    print 'interpolate'
                    self.interpolate_csv_data()
                    print 'Saving Uall, Xall, Yall, Zall npy files'
                    self.save_all_data_npy()
                else:
                    self.load_all_data_npy()
            self.extract_row_from_all()
            self.set_uprime_r()

        #print 'Done! :) for i=%i j=%i'%(self.xi,self.yi)
        self.calc_R11_r_sim()
        self.calc_integral_length_L()
        self.calc_taylor_length_lambda()

        self.save_R11_r(f_handle,)
    ###################################################################
    ###################################################################
    # Time integral functions
    def extract_pt_from_all(self):
        self.X = self.Xall[self.xi,self.yi,0]
        self.Y = self.Yall[self.xi,self.yi,0]
        self.Z = self.Zall
        self.U = self.Uall[:,self.xi,self.yi,:]

    def set_manufactured_tau(self):
        self.u=np.empty((len(self.t),)+(20,))
        if self.option==1:
            self.u = (np.cos(self.t)*np.exp(-1.*self.t)).reshape((len(self.t),1)) + 0.08*(2.*np.random.random((len(self.t),)+self.u.shape[1:])-1.0)
        elif self.option==2:
            self.u = np.cos(self.t).reshape((len(self.t),1)) + 0.08*(2.*np.random.random((len(self.t),)+self.u.shape[1:])-1.0)
    def set_uprime_tau(self):
        # calc u' v' and w'  velocity prime values
        u_avg = np.nanmean(self.U,axis=self.homogeneous)
        self.u = self.U-u_avg

    def calc_R11_tau_sim(self):
        nt = int(self.u.shape[0]/2.) # half of time domain
        self.R11 = np.zeros(nt)
        if self.method==0:
            tau = np.arange(0,nt)
            for taui in tau:
                self.R11[taui] = np.nanmean(self.u[:nt+1]*self.u[taui:nt+taui+1],axis=(self.homogeneous))
        elif self.method==1:
            tau = np.arange(1,nt)
            self.R11[0] =  np.nanmean((self.u*self.u),axis=(self.homogeneous))
            for taui in tau: 
                self.R11[taui] = np.nanmean((self.u[:-1*taui]*self.u[taui:]),axis=(self.homogeneous))
        elif self.method==2:
            tau = np.arange(0,nt)
            for tauispace in tau:
                for tauime in tau:
                    self.R11[tauispace] = np.nanmean(self.R11[tauispace]+(self.u[tauime]*self.u[(tauime+tauispace)]))
            self.R11 = self.R11/float(nt)
        else:
            print 'get right method...',self.method
        self.tau = self.t[:len(self.t)/2]-self.t[0]
        self.R11_norm =  self.R11/self.R11[0]

    def plot_set_label_tau(self):
            self.label='sim,%s,x=%.4f,y=%.4f,m=%i,o=%i'%(
            self.direction,
            self.X ,
            self.Y ,
            self.method ,
            self.option ,)
    def plot_R11_tau_sim(self,ax,marker='-',label='',**kwargs):
        self.plot_set_label_tau()
        ax.plot(self.tau,self.R11_norm,marker,label=self.label+','+label,**kwargs)

        #return T,self.R11_norm
    def plot_axis_labels_tau(self,ax):
        ax.set_xlabel(r'$\tau/\tau_\mathrm{ref}$',**self.labelfont)
        if self.direction=='u':
            ax.set_ylabel(r'$R_{11}(\tau)/R_{11}(0)$',**self.labelfont)
        elif self.direction=='v':
            ax.set_ylabel(r'$R_{22}(\tau)/R_{22}(0)$',**self.labelfont)
        elif self.direction=='w':
            ax.set_ylabel(r'$R_{33}(\tau)/R_{33}(0)$',**self.labelfont)
        ax.legend(loc='best',numpoints=1)
        plt.tight_layout()
    def calc_integral_time_T(self):
        # calc tau integral scale
        self.T = np.trapz(self.R11_norm,x=self.tau)
        print 'T= ',self.T
    def calc_taylor_time_lambda(self,ax=None):
        t=self.tau
        R11=self.R11_norm

        dR=np.gradient(R11,edge_order=self.option)
        d2R=np.gradient(dR,edge_order=self.option)
        dt=np.gradient(t,edge_order=self.option)
        d2Rdr2=d2R/dt**2
        f_prime_prime=d2Rdr2[0]

        #print 'f_prime_prime = ',f_prime_prime,self.xi,r[self.xi]
        self.lambdaf=(-f_prime_prime/2.)**-0.5

        self.plot_set_label_tau()
        if ax!=None:
            P=1.+0.5*f_prime_prime*(t-t[0])**2
            ax.plot(t,P,'s',label='taylor '+self.label)
            ax.plot(self.lambdaf,0.,'o',label=r'$\lambda_f$')
        print 'TaylorT=',self.lambdaf

    def Manufactured_R11_tau(self):
        L = self.t[-1]-self.t[0]
        # option 1
        if self.option==1:
            R11 = (
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
            R11 = (1./L)*(
                    np.cos(self.tau)
                    * (( np.sin(L)      *np.cos(L)+L )  /   2.)
                    - ((  np.sin(self.tau)   *np.sin(L)**2)  /   2.))
        R11_norm = R11/R11[0]
        return R11_norm

    def R11_tau_calc_plot(self,ax,label=''):
        # calculate and plot R11(tau)
        # for manufactured solution only
        if self.manufactured:
            self.t = np.linspace(0,0.05/(0.00745/27.5),5000.)
            self.set_manufactured_tau()
        else:
            if self.read_from=='box_probe':
                self.read_probe()
            elif self.read_from=='vtu':
                if os.path.isfile(self.sim_dir+'Uall_'+self.direction+'.npy')==False:
                    print 'vtu to csv'
                    self.slice_orient_save_vtu_to_csv()
                    print 'read csv'
                    self.read_csv_files()
                    print 'interpolate'
                    self.interpolate_csv_data()
                    print 'Saving Uall, Xall, Yall, Zall npy files'
                    self.save_all_data_npy()
                else:
                    self.load_all_data_npy()
            self.extract_pt_from_all()
            self.set_uprime_tau()

        self.calc_R11_tau_sim()
        self.calc_integral_time_T()
        self.calc_taylor_time_lambda(ax=ax)

        self.plot_R11_tau_sim(ax,marker='+',linewidth=0,label=label)

        self.plot_axis_labels_tau(ax)

        if self.manufactured:
            R11simMNF = self.Manufactured_R11_tau()
            self.plot_set_label_tau()
            ax.plot(self.tau,R11simMNF,marker='x',linewidth=0,label=self.label+' analytic')
    def save_R11_tau(self,f_handle,):
        np.savetxt(f_handle,[[self.X,self.Y,self.T,self.lambdaf]],delimiter=' ')
        if self.direction=='u':
            np.save(self.sim_dir+'R11_output/XYR11/ijXYR11_%i_%i_%.4f_%.4f'%(self.xi,self.yi,self.X,self.Y),self.R11_norm)
        elif self.direction=='v':
            #np.save('R22_output/XYR22/ijXYR22_'+str(self.xi)+'_%.4f'%(self.X),self.R11_norm)
            np.save(self.sim_dir+'R22_output/XYR22/ijXYR22_%i_%i_%.4f_%.4f'%(self.xi,self.yi,self.X,self.Y),self.R11_norm)
        elif self.direction=='w':
            #np.save('R33_output/XYR33/ijXYR33_'+str(self.xi)+'_%.4f'%(self.X),self.R11_norm)
            np.save(self.sim_dir+'R33_output/XYR33/ijXYR33_%i_%i_%.4f_%.4f'%(self.xi,self.yi,self.X,self.Y),self.R11_norm)
        print 'saved ijxy=',self.xi,self.yi,self.X,self.Y
    def R11_tau_calc_save(self,f_handle,label=''):
        # calculate and plot R11(tau)
        if self.loaded==False:
            if self.read_from=='box_probe':
                self.read_probe()
                self.loaded=True
            elif self.read_from=='vtu':
                if os.path.isfile(self.sim_dir+'Uall_'+self.direction+'.npy')==False:
                    print 'vtu to csv'
                    self.slice_orient_save_vtu_to_csv()
                    print 'read csv'
                    self.read_csv_files()
                    print 'interpolate'
                    self.interpolate_csv_data()
                    print 'Saving Uall, Xall, Yall, Zall npy files'
                    self.save_all_data_npy()
                else:
                    self.load_all_data_npy()
                self.loaded=True
        self.extract_pt_from_all()
        self.set_uprime_tau()

        self.calc_R11_tau_sim()
        self.calc_integral_time_T()
        self.calc_taylor_time_lambda()
        self.save_R11_tau(f_handle,)

    ###################################################################

if __name__=="__main__":
    # user options
    calc_options1={
            'direction':'u',
            'sim_dir':'../Create_Plots/DA_paper/OI/C_0/SOLUT/',
            'manufactured':False,
            'method':0,
            'option':2,             # order of taylor derivative difference scheme (option 2 means second order accurate), or option for manufactured solution type
            'xi':4,                 # if box_probe [0-19] elif vtu [0-80]
            'yi':9,                 # if box_probe [0-10] elif vtu [0-90]
            'read_from':'vtu',      # ['box_probe','vtu']
            'vtu':'Box_Output',     # vtu file name
           #'read_from':'box_probe' # ['box_probe','vtu']
            'sim_var':'U',
            #'X_interp_pts':126j,
            #'Y_interp_pts':169j,
            'X_interp_pts':12j,
            'Y_interp_pts':16j,
            'Z_interp_pts':26j,
            # 'Z_interp_pts':1j,
            }
    sim1 = sim(calc_options1)
    sim1.xi=0

    #  fig = plt.figure()
    #  ax = plt.subplot(111)
    #  # for R11(tau)
    #  #sim1.R11_tau_calc_plot(ax)
    #  #sim1.plot_axis_labels_tau(ax)
    #  # for R11(r)
    #  sim1.R11_r_calc_plot(ax)
    #  sim1.plot_axis_labels_r(ax)
    #  # for R11(r) using vtu data

    #  plt.show()




    ############
    # save all length scales
    import sys
    direction=sys.argv[1]
    sim1.direction=direction
    sim1.xi=0
    sim1.yi=0
    if direction=='u':
        with open('R11_output/R11.txt','a') as f_handle:
            #for xi in range(126):
            for xi in range(int((sim1.X_interp_pts*-1j).real)):
                sim1.xi=xi
                #for yi in range(169):
                for yi in range(int((sim1.Y_interp_pts*-1j).real)):
                    sim1.yi=yi
                    sim1.R11_r_calc_save(f_handle)
    elif direction=='v':
        with open('R22_output/R22.txt','a') as f_handle:
            #for xi in range(126):
            for xi in range(int((sim1.X_interp_pts*-1j).real)):
                sim1.xi=xi
                #for yi in range(169):
                for yi in range(int((sim1.Y_interp_pts*-1j).real)):
                    sim1.yi=yi
                    sim1.R11_r_calc_save(f_handle)
    elif direction=='w':
        with open('R33_output/R33.txt','a') as f_handle:
            #for xi in range(126):
            for xi in range(int((sim1.X_interp_pts*-1j).real)):
                sim1.xi=xi
                #for yi in range(169):
                for yi in range(int((sim1.Y_interp_pts*-1j).real)):
                    sim1.yi=yi
                    sim1.R11_r_calc_save(f_handle)
