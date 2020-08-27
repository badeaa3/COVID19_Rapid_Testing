'''
Authors: Anthony Badea and Beatrice Nash (June 2020)
Advisors: Bobby Brooke Herrera, Miguel Bosch, Irene Bosch, Ankita Reddy, Elena Naumova 
Construct and run the SIDARTH-Q model. 
'''

import os 
import sys
import copy
import numpy as np
import random
import modplot
import matplotlib.pyplot as plt
import argparse
import matplotlib.dates as mdates
from datetime import date, timedelta, datetime
# imports for MA data
sys.path.insert(1, '../data_processing')
from load_MA_data import *

# Input: 
#   - verbosose: to print out comments
# Output:
#   - dictionary of model and links
#   - dictionary of parameters
def make_model(verbose=False):
    # ## Constructing the SIDARTHE-Q Model
    # Constructing dictionaries to hold the state and parameters

    # Model dictionary
    # Each key is a state with value a list containing tuples of the connecting model
    # Each tuple contains (state ,parameters,coefficient) as a representation of a link
    model = {'S':[],
             'I':[],
             'D':[],
             'A':[],
             'R':[],
             'T':[],
             'H':[],
             'E':[],
             'Q':[],
             'UH':[],
             'HQ':[]
    }

    # Parameter dictionary
    # Each key is a parameter with value being a list of parameters
    # A list is used to enable documentation of the parameters as they vary in time
    # Alternatively, if constant parameters are chosen multiple constants can be looped over
    parameters = {'alpha':[],
                 'beta':[],
                 'gamma':[],
                  'delta':[],
                 'epsilon':[],
                 'theta':[],
                 'zeta':[],
                 'eta':[],
                 'mu':[],
                 'nu':[],
                 'tau':[],
                 'lambda':[],
                 'kappa':[],
                 'xi':[],
                 'rho':[],
                 'sigma':[],
                 'phi':[],
                 'omega':[],
                 'psi':[],
                 'pi':[]
                 }


    # ## Construct the links between the states of the SIDARTHE Model

    # In[3]:


    # Susceptible
    model['S'].append([['S','I'],['alpha'],-1])
    model['S'].append([['S','D'],['beta'],-1])
    model['S'].append([['S','A'],['gamma'],-1])
    model['S'].append([['S','R'],['delta'],-1])
    model['S'].append([['S'],['omega'],-1])
    model['S'].append([['Q'],['phi'],1])

    # Quarantined
    model['Q'].append([['S'],['omega'],1])
    model['Q'].append([['Q'],['phi'],-1])

    # Infected (asymptomatic, infected)
    model['I'].append([['S','I'],['alpha'],1])
    model['I'].append([['S','D'],['beta'],1])
    model['I'].append([['S','A'],['gamma'],1])
    model['I'].append([['S','R'],['delta'],1])
    model['I'].append([['I'],['epsilon'],-1])
    model['I'].append([['I'],['zeta'],-1])
    model['I'].append([['I'],['lambda'],-1])
    model['I'].append([['D'],['psi'],1])

    # Diagnosed (asymptomatic, detected)
    model['D'].append([['I'],['epsilon'],1])
    model['D'].append([['D'],['eta'],-1])
    model['D'].append([['D'],['rho'],-1])
    model['D'].append([['D'],['psi'],-1])

    # Ailing (developing symptoms, undetected)
    model['A'].append([['I'],['zeta'],1])
    model['A'].append([['A'],['theta'],-1])
    model['A'].append([['A'],['mu'],-1])
    model['A'].append([['A'],['kappa'],-1])
    model['A'].append([['R'],['pi'],1])

    # Recovered (developed symptoms, detected)
    model['R'].append([['D'],['eta'],1])
    model['R'].append([['A'],['theta'],1])
    model['R'].append([['R'],['nu'],-1])
    model['R'].append([['R'],['xi'],-1])
    model['R'].append([['R'],['pi'],-1])

    # Threatened (developed life threatening symptoms, undetected)
    model['T'].append([['A'],['mu'],1])
    model['T'].append([['R'],['nu'],1])
    model['T'].append([['T'],['sigma'],-1])
    model['T'].append([['T'],['tau'],-1])

    # (Recognized) Healed
    model['H'].append([['D'],['rho'],1])
    model['H'].append([['R'],['xi'],1])
    model['H'].append([['T'],['sigma'],1])

    # Unrecognized Healed
    model['UH'].append([['I'],['lambda'],1])
    model['UH'].append([['A'],['kappa'],1])
    model['UH'].append([['HQ'],['phi'],1])
    model['UH'].append([['UH'],['omega'],-1])

    # Healed Quarantined
    model['HQ'].append([['UH'],['omega'],1])
    model['HQ'].append([['HQ'],['phi'],-1])


    # Extinct (dead, due to COVID-19)
    model['E'].append([['T'],['tau'],1])


    # In[4]:


    if(verbose):
        # Print out the model
        print('<---- = Entering, ----> = Leaving')
        print('---------------------------------')
        for state in model:
            print('**** State %s links ****'%state)
            for link in model[state]:
                sts, params, coeff = link
                if len(sts) == 2:
                    if coeff < 0:
                        print("{}{} ---- {} ---->".format(sts[0],sts[1],params[0]))
                    if coeff > 0:
                        print("<---- {} ---- {}{}".format(params[0],sts[0],sts[1]))    
                if len(sts) == 1:
                    if coeff < 0:
                        print("{} ---- {} ---->".format(sts[0],params[0]))
                    if coeff > 0:
                        print("<---- {} ---- {}".format(params[0],sts[0]))

    return model,parameters 
  


# Input:
#   - model: dictionary of model
#   - states_init: initial state values
#   - parameters: dictionary of parameters 
#   - parameters_init: initial parameter values
#   - N: population size
#   - verbosose: to print out comments
# Output:
#   - initialized dictionary of state values
#   - initialized dictionary of parameter values
def initialize(model = {},
               states_init = {},
               parameters = {},
               parameters_init = {}, 
               N = 60*(10**6),
               verbose=False):

    # Initialize parameter values
    parameters.update(parameters_init)

    # State dictionary 
    # Copy of model dictionary but will contain the state values at each time step
    states = copy.deepcopy(model)

    # Initialize state values
    states.update(states_init)

    if(verbose):
        print(states)
        print(np.sum([states[s][-1] for s in states.keys()]))

    return states, parameters, N


# Update the parameter values
# Input: 
#   - parameter dictionary
#   - step index to update to
#   - list of testing days
#   - sensitivity: sensitivity of the test (probability of true/false negative/positive)
#   - freq: probability an individual will be tested
# Output: parameter update dictionary
def update_params(parameters = {}, 
                  i = 0, 
                  test_start_day = 15,
                  sensitivity = 1.0,
                  freq = 1.0/14.0,
                  specificity=0.9):
    if i < test_start_day:  
        parameters.update({'alpha':[0.285],
             'beta':[0.005],
             'gamma':[0.005],
             'delta':[0.005],
             'epsilon':[0],
             'theta':[sensitivity],
             'zeta':[0.034],
             'eta':[0.034],
             'mu':[0.008],
             'nu':[0.015],
             'tau':[0.01],
             'lambda':[0.08],
             'kappa':[0.017],
             'xi':[0.017],
             'rho':[0.017],
             'sigma':[0.017]})
    if i >= test_start_day:  
        parameters.update({'alpha':[0.285],
             'beta':[0.005],
             'gamma':[0.005],
             'delta':[0.005],
             'epsilon':[sensitivity*freq],
             'theta':[sensitivity],
             'zeta':[0.034],
             'eta':[0.034],
             'mu':[0.008],
             'nu':[0.015],
             'tau':[0.01],
             'lambda':[0.08],
             'kappa':[0.017],
             'xi':[0.017],
             'rho':[0.017],
             'sigma':[0.017]})
    # if i == 500:  
    #     parameters.update({'alpha':[0.210],
    #                  'beta':[0.005 ],
    #                  'gamma':[0.110],
    #                  'delta':[0.005 ],
    #                  'epsilon':[0.6], 
    #                  'theta':[0.9],
    #                  'zeta':[0.034],
    #                  'eta':[0.034],
    #                  'mu':[0.008],
    #                  'nu':[0.015],
    #                  'tau':[0.01],
    #                  'lambda':[0.08],
    #                  'kappa':[0.017],
    #                  'xi':[0.017],
    #                  'rho':[0.017],
    #                  'sigma':[0.017]})
    if i >= test_start_day:
        parameters.update({'epsilon':[freq*sensitivity]})
        parameters.update({'theta':[sensitivity]})
        parameters.update({'phi':[0.2]})
        parameters.update({'omega':[freq*(1-specificity)]})



# Input:
#   - time_step
#   - num_steps: dictionary of parameters 
#   - test_freq: frequency of testing
#   - sensitivity: sensitivity of the test (probability of true/false negative/positive)
#   - freq: probability an individual will be tested
#   - P_TruePositive: probability of a true positive    (NOT USED)
#   - P_TrueNegative: probability of a true negative    (NOT USED)
#   - P_FalsePositive: probability of a false postive   (NOT USED)
#   - P_FalseNegative: probability of a false negative  (NOT USED)
#   - verbosose: to print out comments
# Output:
#   - dictionary of state values
#   - dictionary of parameter values
def sim_model(time_step = 1,
            num_steps = 300,
            test_freq = 10, 
            sensitivity = 1.0,
            specificity = 0.90,
            freq = 1.0/14,
            P_TruePositive = 0.6,
            P_TrueNegative = 1.0,
            P_FalsePositive = 1.0,
            P_FalseNegative = 1.0,
            states_init = {},
            parameters_init = {},
            N = 60*(10**6),
            verbose=False,
            number_tests = 0,
            quarantine_length_in_days = 14
            ):
    
    model, parameters = make_model(verbose)
    states, parameters, N = initialize(model, states_init, parameters, parameters_init, N, verbose)

    prev_step = {key:states[key][-1] for key in states.keys()}

    test_start_day = 15

    if(verbose):
        print(prev_step)

    number_tests = []
    number_tests.append(0)

    # Finite Difference Numerical Method
    for i in range(num_steps):
        # Update the parameter values
        update_params_args = {'parameters':parameters,
                              'i':i,
                              'test_start_day':test_start_day,
                              'sensitivity': sensitivity,
                              'freq': freq,
                              'specificity':specificity
                            }
        update_params(**update_params_args)
                   
        # get last time step values
        prev_step = {key:states[key][-1] for key in states.keys()}
        for state in model:
            update = 0
            for link in model[state]:
                sts, params, coeff = link
                state_vals = [prev_step[i] for i in sts]
                param_vals = [parameters[i][-1] for i in params]
                coeff_vals = [coeff]
                # multiply the link values together
                update += np.array(state_vals + param_vals + coeff_vals).prod()
            # multiply by time step 
            update *= time_step
            # add the previous state value
            print(state)
            print(prev_step[state])
            update += prev_step[state]
            # append to the end of the state entry in the states dictionary
            if update < 0: 
                print("State {} is less than zero: {}".format(state,update))
                print("Setting equal to zero")
                update = 0
            states[state].append(update)
        if num_steps > 15:
            number_tests.append(freq*(states["S"][-1]+states["I"][-1])+(1/quarantine_length_in_days)*(states["D"][-1]+states["Q"][-1])+states["A"][-1])
        else:
            number_tests.append(0.4*states["A"][-1])
        if round(np.sum([states[s][-1] for s in states.keys()])) != 1:
            print("Total population doe not equal 1!")
            #break

    return states, parameters, number_tests



# Input:
#   - dictionary of state values
#   - dictionary of parameter values
#   - test_freq: frequency of testing
#   - sensitivity: sensitivity of the test           (NOT USED)
#   - outname: name base for all output files
#   - verbosose: to print out comments
# Output:
#   - 'outname'.pdf plot saved to SIDARTHE-Q_Plots directory
def plot_model(states,parameters,freq, sensitivity= 0.0,specificity=0.0, outName='State_Counts',verbose=False,number_tests=[], plot_MA_data=False):

    os.environ['PATH'] = os.environ['PATH'] + ':/Library/TeX/texbin/tex'
    plt.rcParams['figure.dpi'] = 100
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['text.usetex'] = False
    # Turn interactive plotting off
    # plt.ioff()

    # Number of steps in one of the states
    num_steps = len(states[list(states.keys())[0]])

    # a dictionary to hold plotting information about the observables
    obs = {}
    obs.setdefault('State_Counts',{}).update({
        'xlim': (0,num_steps), 'ylim': (0,1), 'ylim_ratio': (0,1.0),
        'xlabel': 'Time Step', 'ylabel': 'Fraction of Population', 'ylabel_ratio':"Ratio",
        'legend_loc': 'upper right', 'legend_ncol': 2, 'ytick_step': 0.1, 'stamp_xy': (0.05, 0.95)
    })
    color_cycle = ['#FA7200','#FAD400','#BBBBBB','black','#00CE0C','#C60000','#2B7BFE','#6C00B3','#FF1493','#EC7063']
    histStyle = {'ls': '-', 'marker': 'o', 'ms': 1, 'zorder': 1}
    legend_opts = {'handlelength': 2.0, 'loc': 'lower right', 'frameon': False, 'numpoints': 1}

    # Plot the results
    ob = obs['State_Counts']
    fig, [ax0, ax1] = modplot.axes(figsize=(12,6),**ob)

    if(plot_MA_data):
      # Format axis to be dates
      formatter = mdates.DateFormatter("%Y-%m-%d")
      locator = mdates.MonthLocator()
      locator_days = mdates.DayLocator()

      ax0.xaxis.set_major_formatter(formatter)
      ax0.xaxis.set_major_locator(locator)
      #ax0.xaxis.set_minor_locator(locator_days)
      ax1.xaxis.set_major_formatter(formatter)
      ax1.xaxis.set_major_locator(locator)
      #ax1.xaxis.set_minor_locator(locator_days)

      # Load the MA data and plot
      cases_by_date = load_cases_by_date()
      date_of_death = load_date_of_death()
      hosp_from_hosps = load_hospitilization_from_hospitals()

      scale_factor = 1000
      N_MA = 6.893*(10**6)/scale_factor
      case_dates = [datetime.strptime(d,"%m/%d/%Y").date() for d in cases_by_date['date']]
      ax0.plot(case_dates, np.array(cases_by_date['cases_pos_new'])/N_MA, label='MA Confirmed Daily Cases', color=color_cycle[0%len(color_cycle)], linestyle='dashed')

      hosp_dates = [datetime.strptime(d,"%m/%d/%Y").date() for d in hosp_from_hosps['date']]
      ax0.plot(hosp_dates, np.array(hosp_from_hosps['hosp_tot_today'])/N_MA, label="MA Total Daily Hospitalizations", color=color_cycle[1%len(color_cycle)], linestyle='dashed')

      death_dates = [datetime.strptime(d,"%m/%d/%Y").date() for d in date_of_death['date']]
      ax0.plot(death_dates, np.array(date_of_death['death_conf_tot'])/N_MA, label="MA Confirmed Total Deaths", color=color_cycle[2%len(color_cycle)], linestyle='dashed')

      # Plot the simulation 
      sdate = datetime.strptime(min(cases_by_date['date']),"%m/%d/%Y").date()
      model_dates = [sdate]
      for i in range(num_steps-1):
          day = sdate + timedelta(days=i)
          model_dates.append(day)
          
      ax0.set_xlim(min(model_dates), max(model_dates))
      ax1.set_xlim(min(model_dates), max(model_dates))
      fig.autofmt_xdate()

      fig.subplots_adjust(left=0.12, bottom=0.2, right=0.9, top=0.88, wspace=0.2,hspace=0.5)
      #print(model_dates)
      # Plot Infected, Threatened, Extinct states
      #ax0.plot(model_dates,np.array(states["D"]) + np.array(states["R"])*scale_factor, label='D+R', color=color_cycle[0%len(color_cycle)], **histStyle)
      #ax0.plot(model_dates,np.array(states["T"])*scale_factor, label='T', color=color_cycle[1%len(color_cycle)], **histStyle)
      #ax0.plot(model_dates,np.array(states["E"])*scale_factor, label='E', color=color_cycle[2%len(color_cycle)], **histStyle)
      i = 0
      for key in states.keys():
          # plot the distribution
          ax0.plot(model_dates,states[key], label='{}'.format(key),color=color_cycle[i%len(color_cycle)], **histStyle)
          i += 1 

      ax1.plot(model_dates,[(states["A"][x]+states["I"][x])/(states["A"][x]+states["I"][x]+states["D"][x]+states["T"][x]+states["R"][x]) for x in range(len(states["R"]))],color='grey')    
    
    else: 
      time_steps = range(num_steps)
      i = 0
      for key in states.keys():
          # plot the distribution
          ax0.plot(time_steps,states[key], label='{}'.format(key),color=color_cycle[i%len(color_cycle)], **histStyle)
          i += 1        

      ax1.plot(time_steps,[(states["A"][x]+states["I"][x])/(states["A"][x]+states["I"][x]+states["D"][x]+states["T"][x]+states["R"][x]) for x in range(len(states["R"]))],color='grey')    

    loc, ncol = ob['legend_loc'], ob['legend_ncol']
    order = range(len(states)+3) if plot_MA_data else range(len(states))
    modplot.legend(ax=ax0, frameon=False, order=order, loc=loc, ncol=ncol)

    # stamp to put on the plots
    
    modplot.stamp(*ob['stamp_xy'], delta_y=0.06, ax=ax0,
                  line_0='SIDARTHE-Q Model',
                  line_1='Test freq.: {} days'.format(freq),
                  line_2='Test sensitivity: {}'.format(sensitivity),
                  line_3='Test specificity: {}'.format(specificity))

    figname = outName
    path = os.path.join(os.getcwd(),'SIDARTHE-Q_Plots')
    modplot.save(fig, figname, add_watermark=True, tx=218.5, ty=261.5, plots_dir=path)
    plt.show()
    if len(number_tests) == range(num_steps):
        plt.plot(time_steps,number_tests)
        # plt.plot(time_steps,[sum(number_tests[:k]) for k in range(1,len(number_tests)+1)])
        plt.show()
    # plt.close(fig)

def setup_plot(num_steps=300,specificity=0.9):
    os.environ['PATH'] = os.environ['PATH'] + ':/Library/TeX/texbin/tex'
    plt.rcParams['figure.figsize'] = (10,10)
    plt.rcParams['figure.dpi'] = 240
    plt.rcParams['font.family'] = 'san-serif'
    plt.rcParams['text.usetex'] = False
    # Turn interactive plotting off
    # plt.ioff()

    # Number of steps in one of the states
    num_steps = num_steps

    # a dictionary to hold plotting information about the observables
    obs = {}
    obs.setdefault('State_Counts',{}).update({
        'xlim': (0,num_steps), 
        'xlabel': 'Time Step', 'ylabel': 'Fraction of Population','ylim':(0,1),
        'legend_loc': 'upper right', 'ytick_step': 0.1, 'stamp_xy': (0.05, 0.95),'ratio_plot':False,'legend_ncol':2,'fontsize':'small'
    })
    color_cycle = ['#FA7200','#FAD400','#BBBBBB','black','#00CE0C','#C60000','#2B7BFE','#6C00B3','#FF1493','#EC7063']
    histStyle = {'ls': '-', 'marker': 'o', 'ms': 1, 'zorder': 1}
    legend_opts = {'handlelength': 2.0, 'loc': 'lower right', 'frameon': False, 'numpoints': 1}

    # Plot the results
    ob = obs['State_Counts']
    fig, [ax0] = modplot.axes(**ob)
    
    # modplot.stamp(*ob['stamp_xy'], delta_y=0.06, ax=ax0,
    #               line_0='SIDARTHE-Q Model',
    #               line_1='Test specificity: {}'.format(specificity))
    return fig,ax0


def setup_freq_plot(specificity=0.9):
    os.environ['PATH'] = os.environ['PATH'] + ':/Library/TeX/texbin/tex'
    plt.rcParams['figure.figsize'] = (10,10)
    plt.rcParams['figure.dpi'] = 240
    plt.rcParams['font.family'] = 'san-serif'
    plt.rcParams['text.usetex'] = False
    # Turn interactive plotting off
    # plt.ioff()

    # Number of steps in one of the states

    # a dictionary to hold plotting information about the observables
    obs = {}
    obs.setdefault('State_Counts',{}).update({
        'xlim': (0,1), 
        'xlabel': 'Time Step', 'ylabel': 'Fraction of Population','ylim':(0,.5),
        'legend_loc': 'upper right', 'ytick_step': 0.05, 'stamp_xy': (0.05, 0.95),'ratio_plot':False,'legend_ncol':2,'fontsize':'small'
    })
    color_cycle = ['#FA7200','#FAD400','#BBBBBB','black','#00CE0C','#C60000','#2B7BFE','#6C00B3','#FF1493', '#EC7063']
    histStyle = {'ls': '-', 'marker': 'o', 'ms': 1, 'zorder': 1}
    legend_opts = {'handlelength': 2.0, 'loc': 'lower right', 'frameon': False, 'numpoints': 1}

    # Plot the results
    ob = obs['State_Counts']
    fig, [ax0] = modplot.axes(**ob)

    # ax0.set_yscale("log")
    # ax0.set_xscale("log")
    
    # modplot.stamp(*ob['stamp_xy'], delta_y=0.06, ax=ax0,
    #               line_0='SIDARTHE-Q Model',
    #               line_1='Test specificity: {}'.format(specificity))
    return fig,ax0




def plot_number_infected(states,parameters,freq, fig,ax0,sensitivity= 0.0,specificity=0.0, outName='State_Counts',verbose=False,color='red',linestyle='solid'):


    # additional axis settings
    num_steps = 501
    time_steps = range(num_steps)
    obs = {}
    obs.setdefault('State_Counts',{}).update({
        'xlim': (0,num_steps), 'ylim': (0,1), 
        'xlabel': 'Time Step', 'ylabel': 'Fraction of Population',
        'legend_loc': 'upper right', 'ytick_step': 0.1, 'stamp_xy': (0.05, 0.95),'ratio_plot':False,'legend_ncol':2
    })
    color_cycle = ['#FA7200','#FAD400','#BBBBBB','black','#00CE0C','#C60000','#2B7BFE','#6C00B3','#FF1493', '#EC7063']
    histStyle = {'ls': '-', 'marker': 'o', 'ms': 1, 'zorder': 1}
    legend_opts = {'handlelength': 2.0, 'loc': 'lower right', 'frameon': False, 'numpoints': 1}

    # Plot the results
    ob = obs['State_Counts']
    ax0.plot(time_steps,[states["Q"][x]+states["T"][x]+states["R"][x]+states["D"][x] for x in range(len(states["Q"]))],label='specificity = '+str(specificity)+',freq. = '+format(freq,'.2f'), color=color,linestyle=linestyle)
    ax0.plot(time_steps,[states["Q"][x] for x in range(len(states["Q"]))],label='specificity = '+str(specificity)+',freq. = '+format(freq,'.2f'), color=color,linestyle='dotted')

    # legend style and ordering
    # loc, ncol = ob.get('legend_loc', 'upper right'), ob.get('legend_ncol', 1)
    # order = range(len(states.keys()))

    # plot the legend
    



# Argument parser
def options():
    parser = argparse.ArgumentParser(usage=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-tstep", default=1, help="Length of a time step")
    parser.add_argument("-n", default=300, help="Number of time steps")
    parser.add_argument("-tf", default = 1.0, help="Testing frequency, i.e. number of days between a test")
    parser.add_argument("-tc",default=1.0,help="test specificity")
    parser.add_argument("-ts", default = 1.0, help="Test sensitivity (likely overlaps with probability of true/false negative/positive)")
    parser.add_argument("-ptest", default=1.0/14.0, help="Probability of an individual receiving a test")
    parser.add_argument("-pTP", default=0.0, help="Probability of a true positive")
    parser.add_argument("-pTN", default=0.0, help="Probability of a true negative")
    parser.add_argument("-pFP", default=0.0, help="Probability of a false positive")
    parser.add_argument("-pFN", default=0.0, help="Probability of a false negative")
    parser.add_argument("-name", default='State_Counts', help="Name of the output files")
    parser.add_argument("--v", action="store_true", help="Verbose mode to print out comments")
    parser.add_argument("--plot_MA", action="store_true", help="Plot the MA data on top of the model")
    return parser.parse_args()



# Get initial values for the states and parameters
def get_init_vals_Italy():
    N = 60*(10**6)

    states_init = {'S':[0],
             'I':[200/(N)],
             'D':[20/N],
             'A':[1/N],
             'R':[2/N],
             'T':[0],
             'H':[0],
             'E':[0],
             'Q':[0],
             'HQ':[0],
             'UH':[0]
    }
    states_init.update({'S':[1-np.sum([states_init[s][-1] for s in states_init.keys()])]})

    parameters_init = {'alpha':[0.570],
         'beta':[0.011],
         'gamma':[0.285],
         'delta':[0.011],
         'epsilon':[0.131],
         'theta':[0.371],
         'zeta':[0.034],
         'eta':[0.034],
         'mu':[0.008],
         'nu':[0.015],
         'tau':[0.01],
         'lambda':[0.034],
         'kappa':[0.017],
         'xi':[0.017],
         'rho':[0.034],
         'sigma':[0.017],
            'phi':[0.2], # arbitrarily set
            'omega':[0], # arbitrarily set
            'psi':[0], # arbitrarily set
            'pi':[0]} # arbitrarily set
                 
    

    return N, states_init, parameters_init


def get_init_vals_MA():
    N = 6.893*(10**6)

    cases_by_date = load_cases_by_date()
    date_of_death = load_date_of_death()
    hosp_from_hosps = load_hospitilization_from_hospitals()

    states_init = {'S':[0],
             'I':[np.array(cases_by_date['cases_pos_new'][0])/N],
             'D':[0],
             'A':[0],
             'R':[0],
             'T':[0],
             'H':[0],
             'E':[0],
             'Q':[0],
             'HQ':[0],
             'UH':[0]
    }
    states_init.update({'S':[1-np.sum([states_init[s][-1] for s in states_init.keys()])]})

    parameters_init = {'alpha':[1.3],
                 'beta':[0.011],
                 'gamma':[0.456],
                 'delta':[0.011],
                 'epsilon':[0.6],
                 'theta':[0.371],
                 'zeta':[0.125],
                 'eta':[0.125],
                 'mu':[0.017],
                 'nu':[0.027],
                 'tau':[0.01],
                 'lambda':[0.05],
                 'kappa':[0.017],
                 'xi':[0.017],
                 'rho':[0.034],
                 'sigma':[0.017],
                 'phi':[0.2], # arbitrarily set
                 'omega':[0.05], # arbitrarily set
                 'psi':[0], # arbitrarily set
                 'pi':[0]} # arbitrarily set

    return N, states_init, parameters_init

# Main method
def main():
    
    # get initial parameter values
    N, states_init, parameters_init = get_init_vals_Italy()

    # read in arguments
    ops = options()
    args = {'time_step': ops.tstep,
            'num_steps': ops.n,
            'freq': float(ops.tf),
            'sensitivity': float(ops.ts),
            'specificity':float(ops.tc),
            'P_TruePositive': ops.pTP,
            'P_TrueNegative': ops.pTN,
            'P_FalsePositive': ops.pFP,
            'P_FalseNegative': ops.pFN,
            'states_init': states_init,
            'parameters_init': parameters_init,
            'N': N,
            'verbose': ops.v}
    # run job
    states, parameters,number_tests = sim_model(**args)

    print(parameters,sum(number_tests))

    # write output to txt file
    path = os.path.join(os.getcwd(),'SIDARTHE-Q_Plots')
    state_file = open(os.path.join(path,'{}_states.txt').format(ops.name),'w')
    state_file.write( str(states) )
    state_file.close()
    parameters_file = open(os.path.join(path,'{}_parameters.txt').format(ops.name),'w')
    parameters_file.write( str(parameters) )
    parameters_file.close()

    # plot results
    args_plot = {'states': states,
             'parameters': parameters,
             'freq': float(ops.tf),
             'sensitivity': float(ops.ts),
             'specificity': float(ops.tc),
             'outName': ops.name,
             'verbose': ops.v,
             'number_tests': number_tests,
             'plot_MA_data': ops.plot_MA}
    plot_model(**args_plot)



if __name__ == '__main__':
    main()




