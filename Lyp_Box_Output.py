import matplotlib.pyplot as plt
import Sim_functions2 as sim
import numpy as np

calc_options1={
    'direction':'u',
    'sim_dir':'/home/shaun/Desktop/DA/Create_Plots/DA_paper/Updated_BC/Lyp/per/',
    'read_from':'vtu',      # ['box_probe','vtu']
    #'vtu':'Box_Output',     # vtu file name
    'vtu':'Test',
    'X_interp_pts':126j,
    'Y_interp_pts':169j,
    'Z_interp_pts':26j,
}

sim1=sim.sim(calc_options1)


# # output files to csv (Runs once)
#sim1.slice_orient_save_vtu_to_csv()
sim1.slice_orient_save_vtu_to_csv_with_volume()

