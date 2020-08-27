#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy import integrate
import sys
sys.path.insert(1, '/Users/bnash/Downloads/Covid19Modeling-master3/data_processing')
print(sys.path)
import os
print(os.getcwd())
from load_NYC_data import *
from scipy.optimize import fmin
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
import lmfit 
from IPython.core.display import display, HTML
display(HTML("<style>.container { width:90% !important; }</style>"))


# In[2]:


def series_interp(data = []):
    y = [k for k in data if k >= 0]
    x = [i for i in range(len(data)) if data[i] >= 0]
    if 0 not in x:
        y = [0] + y
        x = [0] + x
    f = interp1d(x,y,kind='linear',bounds_error=False)
    return [data[i] if data[i] >= 0 else f(i) for i in range(len(data))]

def window_avg(data = [], window = 1):
    avg = []
    for i in range(len(data)):
        if i < window:
            avg.append(np.mean(data[:i+1]))
        else:
            avg.append(np.mean(data[i-window+1:i+1]))
    return avg


DATA = {'cases': load_cases_by_date('../Covid19Modeling-master3/data/NYC/raw'), 'hosp':load_hospitilization_from_hospitals('../Covid19Modeling-master3/data/NYC/raw'), 'deaths':load_date_of_death('../Covid19Modeling-master3/data/NYC/raw')}

pop = 8.39*(10**6)
sdate = '04/01/2020'
edate = '07/14/2020'

dates = DATA['cases']['date'][DATA['cases']['date'].index(sdate):DATA['cases']['date'].index(edate)]
dates = [datetime.strptime(i,"%m/%d/%Y").date() for i in dates]


cases_tot = np.array(DATA['cases']['cases_pos_tot'])[DATA['cases']['date'].index(sdate):DATA['cases']['date'].index(edate)] / pop
cases_tot = window_avg(series_interp(cases_tot))

hosp = np.array(DATA['hosp']['hosp_tot_today'])[DATA['hosp']['date'].index(sdate):DATA['hosp']['date'].index(edate)] / pop
hosp_total = hosp
deaths_tot = np.array(DATA['deaths']['death_conf_tot'])[DATA['deaths']['date'].index(sdate):DATA['deaths']['date'].index(edate)] / pop
deaths_tot = window_avg(series_interp(deaths_tot))



# In[3]:


# ratio undetected/detected infected (I/D)

# ratio undetected/detected infected (I/D)
phi = 11
# ratio of transmition rate for detected / undetected
eta_to_alpha = 0.01
# link Quarantine Undetected to S
q_time = 14.0
psi = 1/q_time
print(psi)
# Link Undetected Infected to Recovered
ud_time_recov = 14.0
lambd = 1/ud_time_recov
# Link Detected Infected to Hospitalized
d_time_recov = 14.0
rho = 1/d_time_recov
# Link Hospitalized to Recovered
hosp_los_recov = 14.0
sigma = 1/hosp_los_recov
# Used to Link Hospitalized to Dead
hosp_los_fat = 10.0
hosp_fat_rate = 0.346
tau = hosp_fat_rate/hosp_los_fat
# Ratio of infected detected to recovered
rRD = 10.0

def funct(y,par,dt=1.0):

    # Load parameter values
    S,Qu,I,D,H,R,QR,E,tot,total = y
    # Fit varying parameter values
    #beta = par['beta'].value
    alpha = par['alpha'].value
    eta   = eta_to_alpha*alpha
    eps   = par['eps'].value
    mu    = par['mu'].value
    tau   = par['tau'].value
    nu    = par['nu'].value
    # testing features
    freq = par['freq'].value
    sensitivity = par['sens'].value
    specificity = par['spec'].value
    if freq*sensitivity != 0:
        eps = freq*sensitivity
    if sensitivity != 0:
        nu *= sensitivity
    # link S to Quarantine Undetected
    gamma = freq*(1-specificity)
    if(round(S + Qu + I + D + H + R + QR + E) != 1.0): 
        print("Population not conserved")
        print([s for s in y])
        
    
    
    for k in range(100):
        dt = 1/100
        Sn   = S   + dt*(-eta*S*D            - alpha*S*I       - gamma*S         + psi*Qu)
        Qun  = Qu  + dt*(gamma*S             - psi*Qu)
        totn = tot + dt*((nu)*I)
        In   = I   + dt*(-(eps+nu+lambd)*I   + eta*S*D         + alpha*S*I)
        Dn   = D   + dt*((eps+nu)*I          - rho*D           - mu*(I+D))
        Hn   = H   + dt*(mu*(I+D)            - (sigma+tau)*H)
        Rn   = R   + dt*(rho*D               + lambd*I         + sigma*H         - gamma*R       + psi*QR)
        QRn  = QR  + dt*(gamma*R             - psi*QR)
        En   = E   + dt*(tau*H)
        totaln = total + dt*(eta*S*D         + alpha*S*I)
        S, Qu, I, D, H, R, QR, E, tot,total = Sn, Qun, In, Dn, Hn, Rn, QRn, En, totn,totaln

    return [S, Qu, I, D, H, R, QR, E, tot,total]

def resid(par,y=[],tot=[],H=[],E=[],Phi=[],window=1):
    Sn, Qun, In, Dn, Hn, Rn, QRn, En, totn, totaln = funct(y,par)
    residual = []
    for r in range(window):
        residual += list([10*(tot[r] - totn)/(tot[r]+10**-50), 10*(H[r] - Hn)/(H[r]+10**-50), 10*(E[r] - En)/(E[r]+10**-50), (0.015 - (Hn)/(In+Dn+Hn+10**-50))/0.015])
        y = [Sn, Qun, In, Dn, Hn, Rn, QRn, En, totn,totaln]
        Sn, Qun, In, Dn, Hn, Rn, QRn, En, totn, totaln = funct(y,par)
    residual = np.array(residual)
    #print(residual)
    return residual.ravel()

def get_PD(window=1):

    PD = {'beta':[],'alpha':[],'eps':[],'mu':[],'nu':[],'tau':[],'chisqr':[],'S':[],'Qu':[],'I':[],'D':[],'H':[],'R':[],'QR':[],'E':[],'tot':[],'total':[]}
    i = 0

    while i < len(cases_tot)-window:
        if i == 0:
            Qu0 = 0
    #         I0 = 30*cases_tot[i]
            D0 = cases_tot[i]
            H0 = hosp_total[i]
            R0 = rRD*cases_tot[i]
            I0 = (H0-0.015*(H0+D0))/(0.015)
            QR0 = 0
            E0 = deaths_tot[i]
            S0 = 1 - Qu0 - I0 - D0 - H0 - R0 - QR0 - E0
            tot0 = D0
            total0 = I0 + D0

            prev_step = [S0, Qu0, I0, D0, H0, R0, QR0, E0, tot0,total0]
            PD['S'].append(S0)
            PD['Qu'].append(Qu0)
            PD['I'].append(I0)
            PD['D'].append(D0)
            PD['H'].append(H0)
            PD['R'].append(R0)
            PD['QR'].append(QR0)
            PD['E'].append(E0)
            PD['tot'].append(tot0)
            PD['total'].append(total0)

                
        if i > 0:
            prev_step = [PD['S'][-1],PD['Qu'][-1],PD['I'][-1],PD['D'][-1],PD['H'][-1],PD['R'][-1],PD['QR'][-1],PD['E'][-1],PD['tot'][-1],PD['total'][-1]]
            
        params = lmfit.Parameters()
        #params.add('beta'  , 1.00  , min=0.5, max=2.5  )
        params.add('alpha' , 0.1         , min=0.025, max=1.0 )
        params.add('eps'   , 0.0          , vary=False        )
        params.add('nu'    , 0.17       , min=0.0,max=1.0 )
        params.add('mu'    , 0.008        , min=0.0,max=1.0 )
        params.add('tau'   , 0.05       , min=0.0, max=1.0  )
        params.add('freq'  , 0.0          , vary=False        )
        params.add('sens'  , 0.0          , vary=False        )
        params.add('spec'  , 0.0          , vary=False        )   
        
        
        result = lmfit.minimize(resid, params, args=(prev_step, cases_tot[i+1:i+1+window], hosp_total[i+1:i+1+window], deaths_tot[i+1:i+1+window], [phi for p in range(window)], window),method='leastsq',**{'max_nfev':10000,'epsfcn':1e-7})
    #     result = lmfit.minimize(resid, params, args=(prev_step, cases_tot[i+1], hosp_total[i+1], deaths_tot[i+1], phi, window),method='least_squares',**{'max_nfev':50000})
        #PD['beta'].append(result.params['beta'])
        
        # evaluate model
        for p in range(window):
            if len(PD['S']) < len(cases_tot):
                S,Qu,I,D,H,R,QR,E,tot,total = funct(prev_step,result.params)
                #print(cases_tot[i+1],tot,hosp_total[i+1],H,deaths_tot[i+1],E)
                PD['S'].append(S)
                PD['Qu'].append(Qu)
                PD['I'].append(I)
                PD['D'].append(D)
                PD['H'].append(H)
                PD['R'].append(R)
                PD['QR'].append(QR)
                PD['E'].append(E)
                PD['tot'].append(tot)
                PD['total'].append(total)
                PD['alpha'].append(result.params['alpha'])
                PD['nu'].append(result.params['nu'])
                PD['mu'].append(result.params['mu'])
                PD['tau'].append(result.params['tau'])
                PD['chisqr'].append(result.chisqr)
                prev_step = [S,Qu,I,D,H,R,QR,E,tot,total]
        i += window
    return PD


def run_avg(x, wind=14):
    avg = []
    for i in range(len(x)):
        avg.append(np.mean(x[i-wind+1:i+1]) if i >= wind else np.mean(x[:i+1]))
    return avg

# In[4]:


sdate = '03/05/2020'
sdate = datetime.strptime(sdate,"%m/%d/%Y").date()
model_dates = [sdate]
tick_dates = []
for i in range(1,len(cases_tot)):
    day = sdate + timedelta(days=i)
    model_dates.append(day)
    if(day.day == 1):
        tick_dates.append(day)

def make_fig_dates():
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


# In[5]:


def get_states(freq=1/7,sens=0.8,spec=0.8, wind=1,window=1,PD=None):

    if PD is None:
        PD = get_PD(window)
    else:
        PD = PD
    
    SD = {'S':[],'Qu':[],'I':[],'D':[],'H':[],'R':[],'QR':[],'E':[],'tot':[],'total':[]}
    for i in range(len(cases_tot)-window):
        prev_step = []
        if i == 0:
            S0   = PD['S'][0]
            Qu0  = PD['Qu'][0]
            I0   = PD['I'][0]
            D0   = PD['D'][0]
            H0   = PD['H'][0]
            R0   = PD['R'][0]
            QR0  = PD['QR'][0]
            E0   = PD['E'][0]
            tot0 = PD['tot'][0]
            total0 = PD['total'][0]
            prev_step = [S0, Qu0, I0, D0, H0, R0, QR0, E0, tot0,total0]
            
            SD['S'].append(S0)
            SD['Qu'].append(Qu0)
            SD['I'].append(I0)
            SD['D'].append(D0)
            SD['H'].append(H0)
            SD['R'].append(R0)
            SD['QR'].append(QR0)
            SD['E'].append(E0)
            SD['tot'].append(tot0)
            SD['total'].append(total0)
        if i > 0: 
            prev_step = [SD['S'][-1],SD['Qu'][-1],SD['I'][-1],SD['D'][-1],SD['H'][-1],SD['R'][-1],SD['QR'][-1],SD['E'][-1],SD['tot'][-1],SD['total'][-1]]
        params = lmfit.Parameters()
        params.add('alpha' ,   np.mean(PD['alpha'][i-wind+1:i+1]) if i >= wind else np.mean(PD['alpha'][:i+1] )                                            , vary=False)
        params.add('eps'   , 0      , vary=False)
        params.add('nu'    ,   np.mean(PD['nu'][i-wind+1:i+1])    if i >= wind else np.mean(PD['nu'][:i+1]    )                                            , vary=False)
        params.add('mu'    ,   np.mean(PD['mu'][i-wind+1:i+1])    if i >= wind else np.mean(PD['mu'][:i+1]    )                                            , vary=False)
        params.add('tau'   ,   np.mean(PD['tau'][i-wind+1:i+1])   if i >= wind else np.mean(PD['tau'][:i+1]   )                                            , vary=False)
        params.add('freq'  , freq            , vary=False)
        params.add('sens'  , sens            , vary=False)
        params.add('spec'  , spec            , vary=False)

        # evaluate model
        S,Qu,I,D,H,R,QR,E,tot,total = funct(prev_step,params)

        #print(cases_tot[i+1],tot,hosp_total[i+1],H,deaths_tot[i+1],E)
        SD['S'].append(S)
        SD['Qu'].append(Qu)
        SD['I'].append(I)
        SD['D'].append(D)
        SD['H'].append(H)
        SD['R'].append(R)
        SD['QR'].append(QR)
        SD['E'].append(E)
        SD['tot'].append(tot)
        SD['total'].append(total)
    return SD

