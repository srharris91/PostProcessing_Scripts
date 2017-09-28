import time
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
plt.style.use('seaborn-paper')

calc_options1={
        'direction':'u',
        'sim_dir':'../Create_Plots/DA_paper/OI/C_0/SOLUT/',
        'read_from':'vtu',      # ['box_probe','vtu']
        'vtu':'Box_Output_spread',     # vtu file name
        'sim_var':'U',          # read in velocity
        }
sim1 = sim(calc_options1)

def plot_triang(ax,sim):
    x=sim.xcsv
    y=sim.ycsv
    u=sim.ucsv
    v=sim.vcsv
    w=sim.wcsv
    # create unstructured Delaunay triangulation
    triang=tri.Triangulation(x,y)

    # create plot
    #plt.figure(figsize=(2.5,3))
    if direction==0:
        ax.tricontour(triang,u,100,cmap='jet',vmin=0,vmax=0.7)
        #plt.colorbar(plt.tricontourf(triang,u,100,cmap='jet',vmin=0,vmax=20),ticks=[0,10,20])
        ax.tricontourf(triang,u,100,cmap='jet',vmin=0,vmax=0.7)
        ax.axis('scaled')
    elif direction==1:
        ax.tricontour(triang,v,100,cmap='jet',vmin=-0.05,vmax=0.05)
        #plt.colorbar(plt.tricontourf(triang,v,100,cmap='jet',vmin=-2,vmax=2),ticks=[-2,0,2])
        ax.tricontourf(triang,v,100,cmap='jet',vmin=-0.05,vmax=0.05)
        ax.axis('scaled')
    elif direction==2:
        ax.tricontour(triang,w,100,cmap='jet',vmin=-0.05,vmax=0.05)
        #plt.colorbar(plt.tricontourf(triang,w,100,cmap='jet',vmin=-2,vmax=2),ticks=[-2,0,2])
        ax.tricontourf(triang,w,100,cmap='jet',vmin=-0.05,vmax=0.05)
        ax.axis('scaled')
    #plt.gca().set_aspect('scaled')
    #plt.tricontour(triang,x,colors='k')
    ax.set_xlabel('X/D',)#fontsize=24)
    ax.set_ylabel('Y/D',)#fontsize=24)
    #plt.title('Contour plot')
    #plt.savefig('OI_'+title+'_'+cwd+'__0_4050000.png',bbox_inches='tight')

# import commandline arguments
import argparse

parser = argparse.ArgumentParser(description='Take OI_C_0.4050000.vtu file and save a screenshot of velocity in certain direction.')
    
parser.add_argument('direction', metavar='l', type=str, nargs=1,
                            help='which direction to save (x,y,z)')

args = parser.parse_args()
if args.direction==['x',]:
    print 'saving x'
    slice_orient_save(direction=0,title='X',RGB=[0., 0.0, 0.0, 1.0, 20., 1.0, 0.0, 0.0])
elif args.direction==['y',]:
    print 'saving y'
    slice_orient_save(direction=1,title='Y',RGB=[-2., 0.0, 0.0, 1.0, 2., 1.0, 0.0, 0.0])
elif args.direction==['z',]:
    print 'saving z'
    slice_orient_save(direction=2,title='Z',RGB=[-2., 0.0, 0.0, 1.0, 2., 1.0, 0.0, 0.0])
else:
    print 'something is wrong'

plt.tight_layout()
plt.show()

# usage
#       $ pvpython Render_Mesh2.py x
