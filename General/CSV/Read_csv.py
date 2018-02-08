import numpy as np
import matplotlib.pyplot as plt

def read_csv(self,filename,delimiter=',',names=True,**genfromtxtkwargs):
    """ 
    Read in csv file using filename and given np.genfromtxt arguments
    """
    self.data=np.genfromtxt(filename,delimiter=',',names=True,**genfromtxtkwargs)
    return self.data

def contour_csv(self,ax,*args,**kwargs):
    """
    Plots the csv data in a contour
    *args for tricontourf 
    **kwargs for tricontourf
    """
    plt.colorbar( 
            ax.tricontourf(self.data[:,0],self.data[:,1],self.data[:,2],300,*args,cmap='jet',**kwargs))
    ax.axis('equal')
    #ax.set_title('dudz')
    ax.set_xlabel(self.data.dtype.names[0])
    ax.set_ylabel(self.data.dtype.names[1])
    plt.tight_layout()

if __name__=="__main__":
    # open figure and set style
    plt.style.use('seaborn-paper')
    fig = plt.figure()
    ax = plt.subplot(111)

    # plot using CSV

    # add legend
    ax.legend(loc='best',numpoints=1,frameon=False)
    ax.axis([0,4500000,0,0.01])
    ax.set_title(r'$C_f$ for Schubauer and Klebanof plate (natural transition)')
    plt.tight_layout()

    plt.show()
