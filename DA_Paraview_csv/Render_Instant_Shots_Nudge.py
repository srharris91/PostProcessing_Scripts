import time
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import Sim_functions2 as sim
plt.style.use('seaborn-paper')


def plot_triang(ax,sim,taui):
    x=sim.xcsv[0,:]
    y=sim.ycsv[0,:]
    u=sim.ucsv[0,:]
    # create unstructured Delaunay triangulation
    #triang=tri.Triangulation(x,y)
    triang=tri.Triangulation(-y,x)

    # create plot
    if sim.direction=='u':
        colormap=np.linspace(-0.1,0.8,300)
        ax.tricontour(triang,u,colormap,cmap='jet',vmin=0,vmax=0.7)
        plt.colorbar(ax.tricontourf(triang,u,colormap,cmap='jet',vmin=0,vmax=0.7),ax=ax,ticks=[0,0.35,0.7])
        ax.axis([0.,2.,14.,15.5])
    elif sim.direction=='v':
        colormap=np.linspace(-0.25,0.2,300)
        ax.tricontour(triang,u,colormap,cmap='jet',vmin=-0.05,vmax=0.05)
        plt.colorbar(ax.tricontourf(triang,u,colormap,cmap='jet',vmin=-0.05,vmax=0.05),ax=ax,ticks=[-0.25,-0.05,0.,0.05,0.2])
        ax.axis([0.,2.,14.,15.5])
    elif sim.direction=='w':
        colormap=np.linspace(-0.18,0.22,300)
        ax.tricontour(triang,u,colormap,cmap='jet',vmin=-0.05,vmax=0.05)
        plt.colorbar(ax.tricontourf(triang,u,colormap,cmap='jet',vmin=-0.05,vmax=0.05),ax=ax,ticks=[-0.18,-0.05,0.,0.05,0.22])
        ax.axis([0.,2.,14.,15.5])
    ax.set_xlabel(r'$y/D$',)
    ax.set_ylabel(r'$x/D$',)
    #ax.set_title(sim.sim_dir.split('/')[-2].split('_')[-1])
    if taui=='Exp':
        ax.set_title(taui)
    else:
        ax.set_title(r'$\tau = %s$'%taui)
    #plt.savefig('OI_'+title+'_'+cwd+'__0_4050000.png',bbox_inches='tight')

# import commandline arguments
import argparse

parser = argparse.ArgumentParser(description='Take OI_C_0.4050000.vtu file and save a screenshot of velocity in certain direction.')
    
parser.add_argument('direction', metavar='l', type=str, nargs=1,
                            help='which direction to save (u,v,w)')

args = parser.parse_args()
#if args.direction==['u',]:
# get list of subdirectories
sim_C_0 = '../Create_Plots/DA_paper/Updated_BC/Nudging/C_0/'
sim_C_0_000001 = '../Create_Plots/DA_paper/Updated_BC/Nudging/C_0.000001/'
sim_C_0_0001 = '../Create_Plots/DA_paper/Updated_BC/Nudging/C_0.0001/'
sim_C_0_1 = '../Create_Plots/DA_paper/Updated_BC/Nudging/C_0.1/'
sim_C_exp = '../Create_Plots/DA_paper/Updated_BC/Nudging/C_exp/'
calc_options_C_0={
    'direction':args.direction[0],
    'sim_dir':sim_C_0,#'../Create_Plots/DA_paper/Updated_BC/OI/C_0/',
    'read_from':'vtu',      # ['box_probe','vtu']
    'vtu':'OI_C_0.4050000',     # vtu file name
    'sim_var':'U',          # read in velocity
    }
calc_options_C_0_000001 = calc_options_C_0.copy()
calc_options_C_0_0001 = calc_options_C_0.copy()
calc_options_C_0_1 = calc_options_C_0.copy()
calc_options_C_exp = calc_options_C_0.copy()
calc_options_C_0_000001['sim_dir']=sim_C_0_000001
calc_options_C_0_0001['sim_dir']=sim_C_0_0001
calc_options_C_0_1['sim_dir']=sim_C_0_1
calc_options_C_exp['sim_dir']=sim_C_exp
# may need to uncomment this line below to get the experimental csv file data
#calc_options_C_exp['sim_var']='WTestExp'

# open plots and subplot axis
fig = plt.figure(figsize=(8.36,4))
#axsim  = plt.subplot(111,aspect='equal',adjustable='box-forced')
axsim0 = plt.subplot2grid((4,180),(0,15),rowspan=2,colspan=60,aspect='equal',adjustable='box-forced')
axsim1 = plt.subplot2grid((4,180),(2,0),rowspan=2,colspan=60,aspect='equal',adjustable='box-forced',sharex=axsim0,sharey=axsim0)
axsim2 = plt.subplot2grid((4,180),(2,60),rowspan=2,colspan=60,aspect='equal',adjustable='box-forced',sharex=axsim0,sharey=axsim0)
axsim3 = plt.subplot2grid((4,180),(2,120),rowspan=2,colspan=60,aspect='equal',adjustable='box-forced',sharex=axsim0,sharey=axsim0)
axexp0 = plt.subplot2grid((4,180),(0,105),rowspan=2,colspan=60,aspect='equal',adjustable='box-forced',sharex=axsim0,sharey=axsim0)
# for each simulation, plot in respective axis
for ax,calc_options in zip(
        [ axsim0,           axsim1,                     axsim2,                 axsim3,             axexp0,             ],
        [ calc_options_C_0, calc_options_C_0_0001,    calc_options_C_0_000001,  calc_options_C_0_1, calc_options_C_exp  ]
        ):

    sim1 = sim.sim(calc_options)
    #print 'saving x'
    if os.path.isfile(sim1.sim_dir+sim1.csv+'slice.0.csv'):
        print 'csv already exists... do not slice_and_save'
        pass
    else:
        sim1.slice_orient_save_slice_vtu_to_csv()
    sim1.read_csv_files()
    # read taui
    with open(sim1.sim_dir+'tau.txt') as tau_file:
        taui=tau_file.read().rstrip()

    # now plot
    plot_triang(ax,sim1,taui)
fig.tight_layout()
#plt.show()
plt.savefig('../Plots/'+args.direction[0] + '_velocity_Render_Instant_Shots_Nudge.pdf')
