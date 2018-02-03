import matplotlib.pyplot as plt
import numpy as np
import Sim_functions2 as sim
plt.style.use('seaborn-paper')

#user input
first_steps = 100

calc_options1={
        'direction':'u',
        'sim_dir':'/home/shaun/Desktop/DA/Create_Plots/DA_paper/Updated_BC/Lyp/No_per/',
        'manufactured':False,
        'method':0,
        'option':2,             # order of taylor derivative difference scheme (option 2 means second order accurate), or option for manufactured solution type
        'read_from':'box_probe',      # ['box_probe','vtu']
        #'vtu':'Box_Output',     # vtu file name
       #'read_from':'box_probe' # ['box_probe','vtu']
        'sim_var':'U',
        }
calc_options2=calc_options1.copy()
calc_options2['sim_dir']='/home/shaun/Desktop/DA/Create_Plots/DA_paper/Updated_BC/Lyp/per/'

sim1 = sim.sim(calc_options1)
sim1.read_probe()
print sim1.Uall.shape
u1=sim1.Uall[:first_steps]

sim2=sim.sim(calc_options2)
sim2.read_probe()
u2=sim2.Uall[:first_steps]
print sim2.Uall.shape

udiff=np.abs(u1-u2)

# plot
def plot(time=0):
    plt.figure(figsize=(3.0,5))
    colormap=np.linspace(np.min(udiff),np.max(udiff),300)
    plt.colorbar(plt.contourf(sim1.Yall[:,:,0],sim1.Xall[:,:,0],udiff[time,:,:,0],colormap,cmap='jet'))
    plt.axis('scaled')
    plt.xlabel(r'$Y/D$')
    plt.xlabel(r'$X/D$')
    plt.savefig('../Plots/Lyp_%i_.png'%i)
    plt.close()

[plot(i) for i in range(first_steps)]

#plt.show()
