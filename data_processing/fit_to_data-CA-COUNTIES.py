#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy import integrate
import sys
import os
from load_CA_data import *
from scipy.optimize import fmin
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
import lmfit 
from IPython.core.display import display, HTML
display(HTML("<style>.container { width:90% !important; }</style>"))


# In[14]:


def series_interp(data = []):
    y = [k for k in data if k >= 0]
    x = [i for i in range(len(data)) if data[i] >= 0]
    f = interp1d(x,y,kind='linear',fill_value=0,bounds_error=False)
    return [data[i] if data[i] >= 0 else f(i) for i in range(len(data))]

# pick up CA county data
#=======================================================
cd = load_counties()
cases_tot,hosp_total,deaths_tot = {},{},{}
pops = {}
i = 0
for county, data in cd.items():
    i+=1
    county_pop = data['pop']
    pops[county] = county_pop
    sdate = '04/15/2020' 
    edate = '07/22/2020'
    
    cases_pos_tot  = np.array(data['cases']['cases_pos_tot'])
    cases_prob_tot = np.array(data['cases']['cases_prob_tot'])
    hosp_tot_today = np.array(data['hosp']['hosp_tot_today'])
    hosp_new_today = np.array(data['hosp']['hosp_new'])
    death_conf_tot = np.array(data['deaths']['death_conf_tot'])
    death_prob_tot = np.array(data['deaths']['death_prob_tot'])
    
    dates = data['cases']['date'][data['cases']['date'].index(sdate):data['cases']['date'].index(edate)]
    dates = [datetime.strptime(i,"%m/%d/%Y").date() for i in dates]
    
    cases = np.array(cases_pos_tot)[data['cases']['date'].index(sdate):data['cases']['date'].index(edate)] / county_pop
    cases_tot_c = series_interp(cases)

    
    hosp = np.array(hosp_tot_today)[data['hosp']['date'].index(sdate):data['hosp']['date'].index(edate)] / county_pop
    hosp_total_c = series_interp(hosp)

    
    deaths = np.array(death_conf_tot)[data['deaths']['date'].index(sdate):data['deaths']['date'].index(edate)] / county_pop
    deaths_tot_c = series_interp(deaths)
    cases_tot[county] = cases_tot_c
    hosp_total[county] = hosp_total_c
    deaths_tot[county] = deaths_tot_c
    


# In[15]:


#=======================================================
def eq(par,start_t,end_t,incr,county,sensitivity=0.8,specificity=0.8,freq=1/7,init_conds=[],taus=None):
    #-time-grid-----------------------------------
    start_t-=1
    end_t -=1
    t2  = np.linspace(start_t, end_t, incr)
    #differential-eq-system----------------------
    N = pops[county]
    ud_to_d = 11
    
    def funct2(y,t,taus=None,b=None):
        
        if t < 25:
            hosp_fat_rate = ((0.385-0.346)/25)*t + 0.385
        else:
            hosp_fat_rate = 0.346
        hosp_fat_rate= 0.34
        hosp_los_recov = 10.0
        hosp_los_fat = 7.0
        ud_time_recov = 9.0
        d_time_recov = 9.0
        eta_to_alpha = 0.5
        ''' 
        Note:
        beta = eta+alpha*phi
        With rapid test assume all positive go to quarantine hence eta/alpha appx= 0
        @Anthony Do we inc household transmission by setting eta/alpha = something small instead of 0?
        **For cost evaluation, do we use PCR confirmation of positive tests?  Then people will only quarantine until their result comes back (assume 5 days)
        '''
        
        q_time = 14
        freq = 0
        specificity = 0
        sensitivity = 0
        
       

        S,I,D,T,E,R,Qu,QRu,Ru,tot = y
        
        
        

        phi = ud_to_d

        rho = 1/d_time_recov
        lambd = 1/ud_time_recov
        tau = hosp_fat_rate/hosp_los_fat
        sigma = 1/hosp_los_recov
        gamma = freq*(1-specificity)
        psi = 1/q_time
        
        nu = 0

        eps = par['eps'].value
        beta = par['beta'].value
        
#         mu = (par['mu'].value)*(nu/(eps+nu*sensitivity))
        mu = par['mu'].value*(D+I)/D
        
        tau = par['tau'].value
        
        eta_to_alpha = 0.3
        

        alpha = beta
        

#         eta = alpha*eta_to_alpha

#         eta = 0.001
        if(round(np.sum(y[:9])) != N/N):
            print("Population not conserved")
            print(y,np.sum(y))
        dS = -0.01*S*D-alpha*S*I
        dI = - eps*I-lambd*I+0.01*D*S+alpha*S*I
        dD = eps*I - rho*D-mu*D
        dT = -tau*T-sigma*T+mu*D
        dE = tau*T
        dR = sigma*T+rho*D+lambd*I
        dRu = 0
        dQRu = 0
        dQu = 0
        dtot = eps*I

        return [dS, dI,dD, dT, dE, dR,dQu,dQRu,dRu,dtot]
    #integrate------------------------------------
    

#     phi = sum([1/(14*((1/((1-eps)**d))-1)) for d in range(1,14)])

    phi = ud_to_d
    if len(init_conds) ==0:
        D0 = cases_tot[county][start_t]
        I0 = 2*D0
        A0 = cases_tot[county][start_t]
        T0 = hosp_total[county][start_t]
        E0 = deaths_tot[county][start_t]
        tot0 = cases_tot[county][start_t]
        Qu0,QRu0,Ru0,R0 = 0,0,0,0
        S0 = 1-(I0+R0+D0+T0+E0+Ru0)
    else:
        S0,I0,D0,T0,E0,R0,Qu0,QRu0,Ru0,tot0 = init_conds
    
    initial_cond = [S0,I0,D0,T0,E0,R0,Qu0,QRu0,Ru0,tot0]
    ds = integrate.odeint(funct2,initial_cond,t2,args=(taus,None))
    return (ds[:,0],ds[:,1],ds[:,2],ds[:,3],ds[:,4],ds[:,5],ds[:,6],ds[:,7],ds[:,8],ds[:,9],t2)
#=======================================================


# In[16]:

number_tests = 0

#=======================================================
def eq_rt(par,start_t,end_t,incr,county,testing=True,sensitivity=0.8,specificity=0.8,freq=1/7,init_conds=[],taus=None):
    #-time-grid-----------------------------------
    start_t-=1
    end_t -=1
    testing = testing
    t2  = np.linspace(start_t, end_t, incr)
    #differential-eq-system----------------------
    N = 4.0*(10**6)
    ud_to_d = 11

    def funct2(y,t,taus=None,b=None,freq=0):

        global number_tests
        
        hosp_los_recov = 10.0
        hosp_los_fat = 7.0
        ud_time_recov = 9.0
        d_time_recov = 9.0
        eta_to_alpha = 0.5
        hosp_fat_rate = 0.0
        ''' 
        Note:
        beta = eta+alpha*phi
        With rapid test assume all positive go to quarantine hence eta/alpha appx= 0
        @Anthony Do we inc household transmission by setting eta/alpha = something small instead of 0?
        **For cost evaluation, do we use PCR confirmation of positive tests?  Then people will only quarantine until their result comes back (assume 5 days)
        '''
        
        q_time = 14
        
       

        S,I,D,T,E,R,Qu,QRu,Ru,tot = y

        # freq = 0
        if testing and cases_tot[county][min(int(np.floor(t)),len(cases_tot[county])-1)] >= 0.0005 and cases_tot[county][min(int(np.floor(t)),len(cases_tot[county])-1)] <= 0.003:
            print("LOW transmission")
            freq = 0.33
        elif testing and cases_tot[county][min(int(np.floor(t)),len(cases_tot[county])-1)] >= 0.003 and cases_tot[county][min(int(np.floor(t)),len(cases_tot[county])-1)] <= 0.01:
            print("MED transmission")
            freq = 0.5
        elif testing and cases_tot[county][min(int(np.floor(t)),len(cases_tot[county])-1)] > 0.01:
            print("  transmission")
            freq = 1
        elif testing:
            freq = 0
        else:
            freq = freq

        number_tests += freq*(S+I+R+Ru)*pops[county]
        

        phi = ud_to_d

        rho = 1/d_time_recov
        lambd = 1/ud_time_recov
        tau = par['tau'].value
        sigma = 1/hosp_los_recov
        gamma = freq*(1-specificity)
        psi = 1/q_time

        eps = freq*sensitivity
        
        nu = par['nu'].value
        

        beta = par['beta'].value
        
        '''
        mu is scaled by 1/(I_0 + D_0)
        Where I_0 + D_0 = phi
        Then multiplied by I+D
        And removed from population D
        (Assume that rapid test diagnoses individuals before they 
        develop symptoms severe enough to warrant a hospital visit)
        '''
        
        mu = par['mu'].value*(D+I)/D
        
        '''
        Eta to alpha is for the OLD ETA
        that is ETA WITHOUT RAPID TEST
        THIS IS BASED ON THE SIDARTHE PAPER,
        but is not really an important value choice
        given phi >> eta
        '''
        
        eta_to_alpha = 0.5
        
        alpha = beta/(phi+eta_to_alpha)
        
        '''Eta very small - this is the new eta in the rapid test scenario, in which quarantine is ideally highly effective'''
        
        eta = 0.0001

        if(round(np.sum(y[:9])) != N/N): 
            print("Population not conserved")
            print(y,np.sum(y))

        dS = -eta*S*D-alpha*S*I-gamma*S+psi*Qu
        dI = - eps*I-lambd*I+eta*D*S+alpha*S*I-nu*sensitivity*I
        dD = eps*I - rho*D -mu*D+nu*sensitivity*I
        dT = mu*D-tau*T-sigma*T
        dE = tau*T
        dR = sigma*T+rho*D+psi*QRu-gamma*R
        dRu = lambd*I-gamma*Ru+psi*QRu
        dQRu = gamma*Ru-2*psi*QRu+gamma*R
        dQu = gamma*S-psi*Qu
        dtot = eta*D*S+alpha*S*I

        

        return [dS, dI,dD, dT, dE, dR,dQu,dQRu,dRu,dtot]
    #integrate------------------------------------
    

#     phi = sum([1/(14*((1/((1-eps)**d))-1)) for d in range(1,14)])

    phi = ud_to_d
    if len(init_conds) ==0:
        D0 = cases_tot[county][start_t]
        I0 = 11*D0
        A0 = cases_tot[county][start_t]
        T0 = hosp_total[county][start_t]
        E0 = deaths_tot[county][start_t]
        tot0 = cases_tot[county][start_t]
        Qu0,QRu0,Ru0,R0 = 0,0,10*cases_tot[county][start_t],0
        S0 = 1-(I0+R0+D0+T0+E0+Ru0)
    else:
        S0,I0,D0,T0,E0,R0,Qu0,QRu0,Ru0,tot0 = init_conds
    
    initial_cond = [S0,I0,D0,T0,E0,R0,Qu0,QRu0,Ru0,tot0]

    ds = integrate.odeint(funct2,initial_cond,t2,args=(taus,None,freq))
    return (ds[:,0],ds[:,1],ds[:,2],ds[:,3],ds[:,4],ds[:,5],ds[:,6],ds[:,7],ds[:,8],ds[:,9],t2)


# In[18]:


# Score Fit of System
#=========================================================

def tot_resid(parms,start_t,end_t,n_t,county,init_conds=[],taus=None):
    # evaluate model and bar negative state values by strongly penalizing them
    tau = None
    if taus is not None:
        tau = taus
    S,I,D,T,E,R,Qu,QRu,Ru,tot,t=eq(parms,start_t,end_t,n_t,county,init_conds,tau)
    resid = []
    for k in range(start_t-1,end_t):
        resid.append(tot[k-start_t+1]-cases_tot[county][k])
        resid.append(T[k-start_t+1]-hosp_total[county][k])
        resid.append(E[k-start_t+1]-deaths_tot[county][k])
#     resid = list([tot-cases_tot[start_t-1:end_t]]) + list([T-hosp_total[start_t-1:end_t]]) + list([E-deaths_tot[start_t-1:end_t]])
    resid = np.array(resid)
    return resid.ravel()


# In[ ]:

def fit(county='Los Angeles'):
    spacing = len(cases_tot[county])-2
    end = len(cases_tot[county])

    n = end
    t2  = list(np.linspace(1, end, end))
    params_dict = {'beta':[],'eps':[],'mu':[],'tau':[],'phi':[]}

    init_conds = []

    states_dict = {'S':[],'I':[],'D':[],'T':[],'E':[],'R':[],'Qu':[],'QRu':[],'Ru':[],'tot':[]}

    mu = 0.3

    for i in range(0,len(t2)):
        if (t2[i]-1)%spacing == 0 and int(t2[i]) + spacing < end:
            k = t2[i]
            paramsD = lmfit.Parameters()
            if len(params_dict['beta']) == 0:
                paramsD.add('beta',0.2, min=0.1,max=1.0)
                paramsD.add('eps', 0.02, min=0.001,max=1.0)
                paramsD.add('mu', 0.04, min=0.0,max=0.2)
                paramsD.add('tau',0.03, min=0.0,max=0.05)
                paramsD.add('phi',11,vary=False)
            else:
                paramsD.add('beta', params_dict['beta'][-1], min=0.1,max=1.0)
                paramsD.add('eps', params_dict['eps'][-1],min=0.001,max=1.0)
                paramsD.add('mu', params_dict['mu'][-1],  min=0.0,max=0.2)
                paramsD.add('tau',params_dict['tau'][-1], min=0.0,max=0.05)
                paramsD.add('phi',params_dict['phi'][-1],vary=False)
            f = min(int(k)+spacing,end)
    #         tau = np.mean(taus[int(k)-1:f])
    #         print(tau)
            tau = 0
            resultD = lmfit.minimize(tot_resid, paramsD, args=(int(k),f,f-int(k)+1,county,init_conds,tau),method='leastsquares',**{'max_nfev':10000})
            lmfit.report_fit(resultD)
            for c in range(f-int(k)):
                params_dict['beta'].append(resultD.params['beta'])
                params_dict['eps'].append(resultD.params['eps'])
                params_dict['mu'].append(resultD.params['mu'])
                params_dict['tau'].append(resultD.params['tau'])
                params_dict['phi'].append(resultD.params['phi'])
            if np.abs(f-end) <= spacing:
                params_dict['beta'].append(resultD.params['beta'])
                params_dict['eps'].append(resultD.params['eps'])
                params_dict['mu'].append(resultD.params['mu'])
                params_dict['tau'].append(resultD.params['tau'])
                params_dict['phi'].append(resultD.params['phi'])
            S,I,D,T,E,R,Qu,QRu,Ru,tot,t=eq(resultD.params,int(k),f,f-int(k)+1,county,init_conds = init_conds,sensitivity=0,freq=0)
            if init_conds == []:
                states_dict['S'] += list(S)
                states_dict['I'] += list(I)
                states_dict['D'] += list(D)
                states_dict['T'] += list(T)
                states_dict['E'] += list(E)
                states_dict['R'] += list(R)
                states_dict['Qu'] += list(Qu)
                states_dict['QRu'] += list(QRu)
                states_dict['Ru'] += list(Ru)
                states_dict['tot'] += list(tot)
            else:
                states_dict['S'] += list(S[1:])
                states_dict['I'] += list(I[1:])
                states_dict['D'] += list(D[1:])
                states_dict['T'] += list(T[1:])
                states_dict['E'] += list(E[1:])
                states_dict['R'] += list(R[1:])
                states_dict['Qu'] += list(Qu[1:])
                states_dict['QRu'] += list(QRu[1:])
                states_dict['Ru'] += list(Ru[1:])
                states_dict['tot'] += list(tot[1:])
            init_conds = [S[-1],I[-1],D[-1],T[-1],E[-1],R[-1],Qu[-1],QRu[-1],Ru[-1],tot[-1]]
    return params_dict,states_dict


# In[7]:


def get_states(params_dict,end,county,spacing=13,freq=1/7,sensitivity=0.8,specificity=0.9,testing=False):
    global number_tests
    number_tests = 0
    spacing = len(cases_tot[county])-2
    n = end
    t2  = list(np.linspace(1, end, end))
    init_conds = []

    states_dict0 = {'S':[],'I':[],'D':[],'T':[],'E':[],'R':[],'Qu':[],'QRu':[],'Ru':[],'tot':[]}

    for i in range(0,len(t2)):
        if (t2[i]-1)%spacing == 0 and int(t2[i]) + spacing < end:
            k = t2[i]
            f = min(int(k)+spacing,end)
            params = {x:params_dict[x][int(k)] for x in params_dict}
            S,I,D,T,E,R,Qu,QRu,Ru,tot,t=eq_rt(params,int(k),f,f-int(k)+1,county,init_conds = init_conds,freq=freq,testing=testing,sensitivity=sensitivity,specificity=specificity)
            if init_conds == []:
                states_dict0['S'] += list(S)
                states_dict0['D'] += list(D)
                states_dict0['T'] += list(T)
                states_dict0['E'] += list(E)
                states_dict0['I'] += list(I)
                states_dict0['R'] += list(R)
                states_dict0['Qu'] += list(Qu)
                states_dict0['QRu'] += list(QRu)
                states_dict0['Ru'] += list(Ru)
                states_dict0['tot'] += list(tot)
            else:
                states_dict0['S'] += list(S[1:])
                states_dict0['D'] += list(D[1:])
                states_dict0['T'] += list(T[1:])
                states_dict0['E'] += list(E[1:])
                states_dict0['I'] += list(I[1:])
                states_dict0['R'] += list(R[1:])
                states_dict0['Qu'] += list(Qu[1:])
                states_dict0['QRu'] += list(QRu[1:])
                states_dict0['Ru'] += list(Ru[1:])
                states_dict0['tot'] += list(tot[1:])
            init_conds = [S[-1],I[-1],D[-1],T[-1],E[-1],R[-1],Qu[-1],QRu[-1],Ru[-1],tot[-1]]
    
    return states_dict0


# In[8]:




def make_fig_dates(tick_dates):
    f = plt.figure(figsize=(10,10))
    ax = f.add_subplot(111)
    formatter = mdates.DateFormatter("%Y-%m-%d")
    locator = mdates.MonthLocator()
    locator_days = mdates.DayLocator()
    ax.xaxis.set_major_formatter(formatter)
    ax.xaxis.set_major_locator(locator)
    ax.set_xlabel("Date",fontsize=18)
    ax.xaxis.set_ticks(tick_dates)
    ax.xaxis.set_tick_params(rotation=45)
    ax.set_ylabel("Percent of Population (%)",fontsize=18)
    ax.tick_params(labelsize=15, direction='out', length=6, width=2)
    #ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0),useMathText=True)
    ax.yaxis.offsetText.set_fontsize(20)
    ax.set_prop_cycle(color=[
    '#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c', '#98df8a',
    '#d62728', '#ff9896', '#9467bd', '#c5b0d5', '#8c564b', '#c49c94',
    '#e377c2', '#f7b6d2', '#7f7f7f', '#c7c7c7', '#bcbd22', '#dbdb8d',
    '#17becf', '#9edae5'])
    
    return f,ax

def make_fig():
    f = plt.figure(figsize=(10,10))
    ax = f.add_subplot(111)
    ax.set_xlabel("Frequency",fontsize=18)
    ax.set_ylabel("Fraction of Population",fontsize=18)
    ax.tick_params(labelsize=15, direction='out', length=6, width=2)
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0),useMathText=True)
    ax.yaxis.offsetText.set_fontsize(14)
    ax.set_prop_cycle(color=[
    '#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c', '#98df8a',
    '#d62728', '#ff9896', '#9467bd', '#c5b0d5', '#8c564b', '#c49c94',
    '#e377c2', '#f7b6d2', '#7f7f7f', '#c7c7c7', '#bcbd22', '#dbdb8d',
    '#17becf', '#9edae5'])
    return f,ax


# In[9]:


# test sensitivity
sensvs = [0.8]
# test specificity
specs = [0.9]
# plot output directory
outDir = '/Users/bnash/Desktop/'
# Geographical location of analysis 
place = 'CA Counties'
# convert all decimals to percents
toPerc = 100
# testing frequencies
freqs = [1]

# Used for saving all plots to one pdf
import matplotlib.backends.backend_pdf
pdf = matplotlib.backends.backend_pdf.PdfPages(outDir+"/covid_"+str(place)+".pdf")

# loop over sensitivities and specificities 
# Plot each state as a time series with each line corresponding to a different testing frequency
states_dict0 = {'S':[0],'I':[0],'D':[0],'T':[0],'E':[0],'R':[0],'Qu':[0],'QRu':[0],'Ru':[0],'tot':[0]}
states_dict = {'S':[0],'I':[0],'D':[0],'T':[0],'E':[0],'R':[0],'Qu':[0],'QRu':[0],'Ru':[0],'tot':[0]}
states_dict1 = {'S':[0],'I':[0],'D':[0],'T':[0],'E':[0],'R':[0],'Qu':[0],'QRu':[0],'Ru':[0],'tot':[0]}

cd = [np.mean([cases_tot[c][i] for c in cases_tot]) for i in range(len(cases_tot["Los Angeles"]))]
mx = [max([cases_tot[c][i] for c in cases_tot]) for i in range(len(cases_tot["Los Angeles"]))]
mn = [min([cases_tot[c][i] for c in cases_tot]) for i in range(len(cases_tot["Los Angeles"]))]
print(cd)
print(mx)
print(mn)

for sens in sensvs:
    for spec in specs:
        # Susceptible
        
        for freq in freqs:
            for c in cases_tot:
                axis = True
                params_dict,sd = fit(c)

                params_dict['nu'] = params_dict['eps']
                end = len(cases_tot[c])
                sd0 = get_states(params_dict,end,c,freq=0,sensitivity=sens,specificity=spec,testing=True)
                for k in states_dict0:
                    if len(sd0[k]) == len(states_dict0[k]):
                        states_dict0[k] = np.array(states_dict0[k]) + np.array(sd0[k])
                        states_dict[k] = np.array(states_dict[k]) + np.array(sd[k])
                    else:
                        states_dict0[k] = [0 for k in range(len(sd0[k]))]
                        states_dict[k] = [0 for k in range(len(sd[k]))]
                        states_dict0[k] = np.array(states_dict0[k]) + np.array(sd0[k])
                        states_dict[k] = np.array(states_dict[k]) + np.array(sd[k])
            freq = number_tests/(len(cases_tot[c])*sum([pops[c] for c in pops]))
            print(number_tests)
            print(1/freq,'FREQ')
            print(states_dict0['tot'])
            testpp = np.around(number_tests/sum([pops[c] for c in pops]),2)
            for c in cases_tot:
                print(freq)
                sd1 = get_states(params_dict,end,c,freq=0,sensitivity=sens,specificity=spec,testing=False)
                for k in states_dict1:
                    if len(sd1[k]) == len(states_dict1[k]):
                        states_dict1[k] = np.array(states_dict1[k]) + np.array(sd1[k]) 
                    else:
                        states_dict1[k] = [0 for k in range(len(sd1[k]))]
                        states_dict1[k] = np.array(states_dict1[k]) + np.array(sd1[k]) 
                end_t1 = len(states_dict0['tot'])
                end=len(cases_tot[c])
                sdate = '04/15/2020'
                sdate = datetime.strptime(sdate,"%m/%d/%Y").date()
            print(states_dict1['tot'])
            model_dates = [sdate]
            tick_dates = []
            for i in range(1,end_t1):
                day = sdate + timedelta(days=i)
                model_dates.append(day)
                if(day.day == 1):
                    tick_dates.append(day)
            if axis:
                fig,ax0 = make_fig_dates(tick_dates)
                axis = False
            label = c
            print(states_dict1['tot'])
            print(states_dict0['tot'])
            # ax0.plot(model_dates,np.array(states_dict['tot'])*toPerc,label='Testing Off') 
            ax0.plot(model_dates,np.array(states_dict0['tot'])*toPerc,label='Location-based \ntesting')
            ax0.plot(model_dates,np.array(states_dict1['tot'])*toPerc,label='Uniform \ntesting')
        ax0.text(0.715,0.85,place+' Total Infections',horizontalalignment='center',transform=ax0.transAxes,fontsize=18)
        ax0.text(0.715,0.9,"Sensitivity: {}, Specificity: {}".format(np.around(sens,1),np.around(spec,1)),horizontalalignment='center',transform=ax0.transAxes,fontsize=18)
        ax0.legend(title = "Intra-county transmission based testing:\n"+str(testpp)+" tests per person over "+str(len(states_dict1['tot']))+" days", title_fontsize=18,fontsize=18,frameon=False,ncol=2,loc="center left", bbox_to_anchor=(0.25,0.75))
        plt.savefig(outDir + "/ca-county-tot-vs-time-"+str(int(10*sens))+"sens-"+str(int(10*spec))+"spec"+".pdf", bbox_inches="tight")
        pdf.savefig(fig,bbox_inches="tight")
        
        # # Undetected Infected
        # fig,ax0 = make_fig_dates()
        # ax0.plot(model_dates,11*np.array(states_dict['D'])*toPerc,label='No Testing')
        # for freq in freqs:
        #     states_dict0 = get_states(params_dict,end,freq=freq,sensitivity=sens,specificity=spec)
        #     label='' #'Days Between Test:'
        #     if freq != 0:
        #         label += str(int(1/freq))
        #     else:
        #         label = 'No Rapid Testing'
        #     ax0.plot(model_dates,np.array(states_dict0['I'])*toPerc,label=label)
        # ax0.text(0.76,0.95,place+' Undetected Infected',horizontalalignment='center',transform=ax0.transAxes,fontsize=18)
        # ax0.text(0.715,0.9,"Sensitivity: {}, Specificity: {}".format(np.around(sens,1),np.around(spec,1)),horizontalalignment='center',transform=ax0.transAxes,fontsize=18)
        # ax0.legend(title = "Days Between Tests", title_fontsize=18,fontsize=18,frameon=False,ncol=2,loc="center left", bbox_to_anchor=(0.385,0.7))
        # plt.savefig(outDir + "/la-undetec-infec-vs-time-"+str(int(10*sens))+"sens-"+str(int(10*spec))+"spec"+".pdf", bbox_inches="tight")
        # pdf.savefig(fig,bbox_inches="tight") 
        
        # # Detected Infected
        # fig,ax0 = make_fig_dates()
        # ax0.plot(model_dates,np.array(states_dict['D'])*toPerc,label='No Testing')
        # for freq in freqs:
        #     states_dict0 = get_states(params_dict,end,freq=freq,sensitivity=sens,specificity=spec)
        #     label='' #'Days Between Test:'
        #     if freq != 0:
        #         label += str(int(1/freq))
        #     else:
        #         label = 'No Rapid Testing'
        #     ax0.plot(model_dates,np.array(states_dict0['D'])*toPerc,label=label)
        # ax0.text(0.78,0.95,place+' Detected Infected',horizontalalignment='center',transform=ax0.transAxes,fontsize=18)
        # ax0.text(0.715,0.9,"Sensitivity: {}, Specificity: {}".format(np.around(sens,1),np.around(spec,1)),horizontalalignment='center',transform=ax0.transAxes,fontsize=18)
        # ax0.legend(title = "Days Between Tests", title_fontsize=18,fontsize=18,frameon=False,ncol=2,loc="center left", bbox_to_anchor=(0.385,0.7))
        # plt.savefig(outDir + "/la-detec-infec-vs-time-"+str(int(10*sens))+"sens-"+str(int(10*spec))+"spec"+".pdf", bbox_inches="tight")  
        # pdf.savefig(fig,bbox_inches="tight")
        
        # # Hospitalized
        # fig,ax0 = make_fig_dates()
        # ax0.plot(model_dates,np.array(states_dict['T'])*toPerc,label='No Testing')
        # for freq in freqs:
        #     states_dict0 = get_states(params_dict,end,freq=freq,sensitivity=sens,specificity=spec)
        #     label='' #'Days Between Test:'
        #     if freq != 0:
        #         label += str(int(1/freq))
        #     else:
        #         label = 'No Rapid Testing'
        #     ax0.plot(model_dates,np.array(states_dict0['T'])*toPerc,label=label)
        # ax0.text(0.825,0.95,place+' Hospitalized',horizontalalignment='center',transform=ax0.transAxes,fontsize=18)
        # ax0.text(0.715,0.9,"Sensitivity: {}, Specificity: {}".format(np.around(sens,1),np.around(spec,1)),horizontalalignment='center',transform=ax0.transAxes,fontsize=18)
        # ax0.legend(title = "Days Between Tests", title_fontsize=18,fontsize=18,frameon=False,ncol=2,loc="center left", bbox_to_anchor=(0.385,0.74))
        # plt.savefig(outDir + "/la-hosp-vs-time-"+str(int(10*sens))+"sens-"+str(int(10*spec))+"spec"+".pdf", bbox_inches="tight")
        # pdf.savefig(fig,bbox_inches="tight")
        
        # # Recovered
        # fig,ax0 = make_fig_dates()
        # ax0.plot(model_dates,(np.array(states_dict['A']) - np.array(states_dict['D'])) *toPerc,label='No Testing')
        # for freq in freqs:
        #     states_dict0 = get_states(params_dict,end,freq=freq,sensitivity=sens,specificity=spec)
        #     label='' #'Days Between Test:'
        #     if freq != 0:
        #         label += str(int(1/freq))
        #     else:
        #         label = 'No Rapid Testing'
        #     ax0.plot(model_dates,np.array(states_dict0['R'])*toPerc,label=label)
        # ax0.text(0.15,0.96,place+' Recovered',horizontalalignment='center',transform=ax0.transAxes,fontsize=18)
        # ax0.text(0.275,0.91,"Sensitivity: {}, Specificity: {}".format(np.around(sens,1),np.around(spec,1)),horizontalalignment='center',transform=ax0.transAxes,fontsize=18)
        # ax0.legend(title = "Days Between Tests", title_fontsize=18,fontsize=18,frameon=False,ncol=2,loc="center left", bbox_to_anchor=(0.385,0.14))
        # plt.savefig(outDir + "/la-rec-vs-time-"+str(int(10*sens))+"sens-"+str(int(10*spec))+"spec"+".pdf", bbox_inches="tight")
        # pdf.savefig(fig,bbox_inches="tight")
        
        # # Extinct
        # fig,ax0 = make_fig_dates()
        # ax0.plot(model_dates,np.array(states_dict['E']) *toPerc,label='No Testing')
        # for freq in freqs:
        #     states_dict0 = get_states(params_dict,end,freq=freq,sensitivity=sens,specificity=spec)
        #     label='' #'Days Between Test:'
        #     if freq != 0:
        #         label += str(int(1/freq))
        #     else:
        #         label = 'No Rapid Testing'
        #     ax0.plot(model_dates,np.array(states_dict0['E'])*toPerc,label=label)
        # ax0.text(0.15,0.96,place+' Deceased',horizontalalignment='center',transform=ax0.transAxes,fontsize=18)
        # ax0.text(0.28,0.91,"Sensitivity: {}, Specificity: {}".format(np.around(sens,1),np.around(spec,1)),horizontalalignment='center',transform=ax0.transAxes,fontsize=18)
        # ax0.legend(title = "Days Between Tests", title_fontsize=18,fontsize=18,frameon=False,ncol=2,loc="center left", bbox_to_anchor=(0.385,0.48)) #(0.36,0.47)
        # plt.savefig(outDir + "/la-deceased-vs-time-"+str(int(10*sens))+"sens-"+str(int(10*spec))+"spec"+".pdf", bbox_inches="tight")
        # pdf.savefig(fig,bbox_inches="tight")
        
        # # Quarantined Uninfected
        # fig,ax0 = make_fig_dates()
        # ax0.plot(model_dates,np.array(states_dict['Qu'])*toPerc,label='No Testing')
        # for freq in freqs:
        #     states_dict0 = get_states(params_dict,end,freq=freq,sensitivity=sens,specificity=spec)
        #     label='' #'Days Between Test:'
        #     if freq != 0:
        #         label += str(int(1/freq))
        #     else:
        #         label = 'No Rapid Testing'
        #     ax0.plot(model_dates,np.array(states_dict0['Qu'])*toPerc,label=label)
        # ax0.text(0.73,0.15,place+' Quarantined Uninfected',horizontalalignment='center',transform=ax0.transAxes,fontsize=18)
        # ax0.text(0.715,0.1,"Sensitivity: {}, Specificity: {}".format(np.around(sens,1),np.around(spec,1)),horizontalalignment='center',transform=ax0.transAxes,fontsize=18)
        # ax0.legend(title = "Days Between Tests", title_fontsize=18,fontsize=18,frameon=False,ncol=2,loc="center left", bbox_to_anchor=(0.385,0.81)) #0.36,0.81
        # plt.savefig(outDir + "/la-qu-vs-time-"+str(int(10*sens))+"sens-"+str(int(10*spec))+"spec"+".pdf", bbox_inches="tight")
        # pdf.savefig(fig,bbox_inches="tight")
        
        # # Quarantined Recovered
        # fig,ax0 = make_fig_dates()
        # ax0.plot(model_dates,np.array(states_dict['QRu'])*toPerc,label='No Testing')
        # for freq in freqs:
        #     states_dict0 = get_states(params_dict,end,freq=freq,sensitivity=sens,specificity=spec)
        #     label='' #'Days Between Test:'
        #     if freq != 0:
        #         label += str(int(1/freq))
        #     else:
        #         label = 'No Rapid Testing'
        #     ax0.plot(model_dates,np.array(states_dict0['QRu'])*toPerc,label=label)
        # ax0.text(0.735,0.15,place+' Quarantined Recovered',horizontalalignment='center',transform=ax0.transAxes,fontsize=18)
        # ax0.text(0.715,0.1,"Sensitivity: {}, Specificity: {}".format(np.around(sens,1),np.around(spec,1)),horizontalalignment='center',transform=ax0.transAxes,fontsize=18)
        # ax0.legend(title = "Days Between Tests", title_fontsize=18,fontsize=18,frameon=False,ncol=2,loc="center left", bbox_to_anchor=(0.385,0.81))
        # plt.savefig(outDir + "/la-qr-vs-time-"+str(int(10*sens))+"sens-"+str(int(10*spec))+"spec"+".pdf", bbox_inches="tight")
        # pdf.savefig(fig,bbox_inches="tight")
        
pdf.close()



