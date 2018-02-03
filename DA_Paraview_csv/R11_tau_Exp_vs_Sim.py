import numpy as np
import matplotlib.pyplot as plt
import Sim_functions as sim
import Exp_functions as exp

# open figure and subplot axes
fig=plt.figure()
ax=plt.subplot(111)
# user options for simulations
sim.sim_dir='Create_Plots/DA_paper/OI/C_0/SOLUT/'
sim.direction=1     # [0,1,2]
sim.option=2        # [1,2]
sim.method=1        # [0,1,2]
sim.manufactured=False
sim.xi,sim.yi=4,9   # [0-19,0-10]
# simulation calculations and add to subplot axes
R11sim1,R11simMNF=sim.R11_tau_calc_plot(ax)
# user inputs for experimental
exp.exp_dir='DA_ExperimentalData/From_JFrank_5000_cold_SandiaC/'
#exp.xi,exp.yi=-3,48 # [0-169,0-126]
exp.xi,exp.yi=-3,48 # [0-169,0-126]
calc_options={
        'direction':'u',        # ['u','v','w']
        'option':2,             # [1,2] for manufactured data
        'method':1,             # [0,1,2] half of domain, as much time domain as possible, and by for loops for all possible time domain
        'manufactured':False,   # use manufactured data instead of experimental data
        'skip_last':True,       # skip last time step
        'moving_mean':True,     # use running average definition when calculating the mean
        'average_z':True,       # average in almost azimuthal direction
        'nwindow':5,
            }
plot_options={ 
        'marker':'+',
        'linewidth':1,
        'label':'exp d=%s,m=%i,o=%i,sk=%i,mave=%i,avez=%i,nwin=%i'%(
            calc_options['direction'],
            calc_options['method'],
            calc_options['option'],
            calc_options['skip_last'],
            calc_options['moving_mean'],
            calc_options['average_z'],
            calc_options['nwindow'])}
# experimental  calculatons and add to subplot axes
nwin,T=[],[]
R11exp1,R11expMNF,T1=exp.R11_tau_calc_plot(ax,calc_options,plot_options)
nwin.append(calc_options['nwindow']),T.append(T1)

# change and replot
# for nwini in np.arange(10,1000,100):
#     calc_options['nwindow']=nwini
#     plot_options[ 'label']='exp d=%s,m=%i,o=%i,sk=%i,mave=%i,avez=%i,nwin=%i'%( calc_options['direction'], calc_options['method'], calc_options['option'], calc_options['skip_last'], calc_options['moving_mean'], calc_options['average_z'], calc_options['nwindow'])
#     R11exp1,R11expMNF,T1=exp.R11_tau_calc_plot(ax,calc_options,plot_options)
#     nwin.append(calc_options['nwindow']),T.append(T1)

calc_options['nwindow']=0
calc_options['moving_mean']=False
calc_options['method']=0
#for yii in np.arange(0,126):
#exp.yi=yii
plot_options[ 'label']='exp d=%s,m=%i,o=%i,sk=%i,mave=%i,avez=%i,nwin=%i'%( calc_options['direction'], calc_options['method'], calc_options['option'], calc_options['skip_last'], calc_options['moving_mean'], calc_options['average_z'], calc_options['nwindow'])
R11exp1,R11expMNF,T1=exp.R11_tau_calc_plot(ax,calc_options,plot_options)
nwin.append(calc_options['nwindow']),T.append(T1)

#
# change and replot
calc_options['method']=1
calc_options['average_z']=False
plot_options[ 'label']='exp d=%s,m=%i,o=%i,sk=%i,mave=%i,avez=%i,nwin=%i'%( calc_options['direction'], calc_options['method'], calc_options['option'], calc_options['skip_last'], calc_options['moving_mean'], calc_options['average_z'], calc_options['nwindow'])
R11exp1,R11expMNF,T1=exp.R11_tau_calc_plot(ax,calc_options,plot_options)

calc_options['moving_mean']=False
plot_options[ 'label']='exp d=%s,m=%i,o=%i,sk=%i,mave=%i,avez=%i,nwin=%i'%( calc_options['direction'], calc_options['method'], calc_options['option'], calc_options['skip_last'], calc_options['moving_mean'], calc_options['average_z'], calc_options['nwindow'])
R11exp1,R11expMNF,T1=exp.R11_tau_calc_plot(ax,calc_options,plot_options)

calc_options['average_z']=True
plot_options[ 'label']='exp d=%s,m=%i,o=%i,sk=%i,mave=%i,avez=%i,nwin=%i'%( calc_options['direction'], calc_options['method'], calc_options['option'], calc_options['skip_last'], calc_options['moving_mean'], calc_options['average_z'], calc_options['nwindow'])
R11exp1,R11expMNF,T1=exp.R11_tau_calc_plot(ax,calc_options,plot_options)
# calculate RMS error
RMSsim1=np.sqrt(((R11sim1-R11simMNF)**2).mean())
RMSexp1=np.sqrt(((R11exp1-R11expMNF)**2).mean())
print RMSsim1,RMSexp1


# plot nwin vs T
fig2=plt.figure()
ax2=plt.subplot(111)
ax2.plot(nwin,T,'k.')
# show plot
plt.show()
