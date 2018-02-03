import numpy as np
import matplotlib.pyplot as plt
import h5py
import Sim_functions2 as sim

# exp data
exp=np.load('../DA_ExperimentalData/From_JFrank_5000_cold_SandiaC/Ave_rms/ave_u.npy')
#exp=np.load('../DA_ExperimentalData/From_JFrank_5000_cold_SandiaC/Ave_rms/RMSu.npy')
# ../DA_ExperimentalData/From_JFrank_5000_cold_SandiaC/Ave_rms/ave_u.npy
# ../DA_ExperimentalData/From_JFrank_5000_cold_SandiaC/Ave_rms/ave_v.npy
# ../DA_ExperimentalData/From_JFrank_5000_cold_SandiaC/Ave_rms/ave_w.npy
# read in X and Y
f = h5py.File('../DA_ExperimentalData/From_JFrank_5000_cold_SandiaC/T093013_01_Velocity_1-500.mat')
expX = np.array(f['X'][0,:])
expY = np.array(f['Y'][0,:])
D_ref=0.00745
U_ref=27.5
expY=expY/(1000.)+15.*D_ref # mm to m and 15 diam. downstream
expX=expX/1000.-(1.15*D_ref) # centered between nozzles
# non-dimensionalize everything
expX = expX/D_ref
expY = expY/D_ref
exp = exp/U_ref
# make mesh grid for plotting
expX,expY=np.meshgrid(expX,expY,indexing='ij')

# simulation data
calc_options1={
        'direction':'v',
        'sim_dir':'../Create_Plots/for_S_AVG_vs_EXP/',
        'manufactured':False,
        'method':0,
        'option':2,             # order of taylor derivative difference scheme (option 2 means second order accurate), or option for manufactured solution type
        'xi':1000,                 # if box_probe [0-19] elif vtu [0-126]
        'yi':9999,                 # if box_probe [0-10] elif vtu [0-169]
        'read_from':'vtu',      # ['box_probe','vtu']
        #'read_from':'box_probe' # ['box_probe','vtu']
        'sim_var':'U_AVG',
        #'sim_var':'U_RMS',
        'X_interp_pts':126j,
        'Y_interp_pts':169j,
        'Z_interp_pts':1j,
        }
sim1 = sim.sim(calc_options1)
#sim1.slice_orient_save_vtu_to_csv()
#sim1.read_csv_files()
#print 'interpolating'
#sim1.interpolate_csv_data()
#print 'saving'
#sim1.save_all_data_npy()
sim1.load_all_data_npy()

# output a few values
#print np.abs(sim1.Yall[::-1,:,0].T - expX)
#print np.abs(sim1.Xall[::-1,:,0].T - expY)
#print np.abs(sim1.Uall[0,::-1,:,0].T - exp)
# creat plot
ax=plt.subplot(121)
fig=plt.gcf()
fig.set_size_inches(12,4)
# exp data
colors=np.linspace(-0.013,0.010,100)
#colors=np.linspace(0.02,0.09,100)
#colors=np.linspace(0.04,0.56,100) # axial avg
#colors=np.linspace(0.023,0.11,100) # axial RMS
#colors=100
plt.colorbar(ax.contourf(expX,expY,exp,colors,cmap='jet'))
ax.set_title('exp')
ax.axis('equal')
ax=plt.subplot(122)
plt.colorbar(ax.contourf(sim1.Yall[:,:,0].T,sim1.Xall[:,:,0].T,sim1.Uall[0,:,:,0].T,colors,cmap='jet'))
ax.set_title('sim')
ax.axis('equal')
plt.tight_layout()
plt.show()

