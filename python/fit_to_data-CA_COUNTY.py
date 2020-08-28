import sys, os, csv
sys.path.insert(1, '../data_processing')
import numpy as np
import lmfit 
# import matplotlib.pyplot as plt
from load_CA_data import *
from county_utilities import *
# import matplotlib.backends.backend_pdf
import time
import multiprocessing 
import pandas as pd 

# plot output directory
outDir = 'county_output'
csvOutDir = outDir
# Geographical location of analysis 
place = 'CA_COUNTIES'

# Set up global output PDF
#pdf = matplotlib.backends.backend_pdf.PdfPages(outDir+"/covid_" + place + ".pdf")

# Load the counties 
cd = load_counties( '../data/CA/raw/')
cases_tot,hosp_total,deaths_tot,pops = {},{},{},{}

sdate = '04/01/2020' 
edate = '07/22/2020'

region_pop = np.sum([cd[county]['pop'] for county in cd.keys()])
total_included_pop = 0

for county in list(cd.keys()):
    #if county not in ["Alameda"] : continue
    DATA = {'cases': cd[county]['cases'], 'hosp':cd[county]['hosp'], 'deaths':cd[county]['deaths']}

    pop = cd[county]['pop'] 
    dates = DATA['cases']['date'][DATA['cases']['date'].index(sdate):DATA['cases']['date'].index(edate)]
    dates = [datetime.strptime(i,"%m/%d/%Y").date() for i in dates]
    ct = np.array(DATA['cases']['cases_pos_tot'])[DATA['cases']['date'].index(sdate):DATA['cases']['date'].index(edate)] / pop
    ht = np.array(DATA['hosp']['hosp_tot_today'])[DATA['hosp']['date'].index(sdate):DATA['hosp']['date'].index(edate)] / pop
    dt = np.array(DATA['deaths']['death_conf_tot'])[DATA['deaths']['date'].index(sdate):DATA['deaths']['date'].index(edate)] / pop

    # Selections on data sets
    if pop/region_pop < 0.015: continue 
    if not any(c>0.0 for c in dt): continue
    total_included_pop += pop

    # Save the data sets passing selection
    pops[county] = pop
    cases_tot[county]  = window_avg(series_interp(ct))
    hosp_total[county] = window_avg(series_interp(ht))
    deaths_tot[county] = window_avg(series_interp(dt))

# Check how much of population is used
print("Percentage of California population used: {}".format(total_included_pop/region_pop))

# ratio undetected/detected infected (I/D)
phi = 11
# ratio of transmition rate for detected / undetected
eta_to_alpha = 0.01
# link Quarantine Undetected to S
q_time = 14.0
psi = 1/q_time
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

def funct(y,par,dt=1.0,thresh=0.0,p=2.0):

    # Load parameter values
    S,Qu,I,D,H,R,QR,E,tot,total = y
    
    # Fit varying parameter values
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
    if freq == -1.0:
        if D <= thresh:
            freq = 1/(p*np.log2(thresh/D)+1)
        else:
            freq = 1.0
        eps = freq*sensitivity
    if freq != -1.0 and freq*sensitivity != 0:
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
        residual += list([(tot[r] - totn)/(tot[r]+10**-50), (H[r] - Hn)/(H[r]+10**-50), (E[r] - En)/(E[r]+10**-50), (0.015 - (Hn)/(In+Dn+Hn+10**-50))/0.015])
        y = [Sn, Qun, In, Dn, Hn, Rn, QRn, En, totn,totaln]
        Sn, Qun, In, Dn, Hn, Rn, QRn, En, totn, totaln = funct(y,par)
    residual = np.array(residual)
    return residual.ravel()

def get_PD(window=1,county="Los Angeles"):

    PD = {'beta':[],'alpha':[],'eps':[],'mu':[],'nu':[],'tau':[],'chisqr':[],'S':[],'Qu':[],'I':[],'D':[],'H':[],'R':[],'QR':[],'E':[],'tot':[],'total':[]}
    i = 0

    while i < len(cases_tot[county])-window:
        if i == 0:
            Qu0 = 0
            D0 = cases_tot[county][i]
            H0 = hosp_total[county][i]
            R0 = rRD*cases_tot[county][i]
            QR0 = 0
            E0 = deaths_tot[county][i]
            I0 = (H0-0.015*(H0+D0))/(0.015)
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
        params.add('alpha' , 0.1         , min=0.01, max=1.0 )
        params.add('eps'   , 0.0          , vary=False        )
        params.add('nu'    , 0.17       , min=0.0,max=1.0 )
        params.add('mu'    , 0.008        , min=0.0,max=1.0 )
        params.add('tau'   , 0.05       , min=0.0, max=1.0  )
        params.add('freq'  , 0.0          , vary=False        )
        params.add('sens'  , 0.0          , vary=False        )
        params.add('spec'  , 0.0          , vary=False        )   
        
        
        result = lmfit.minimize(resid, params, args=(prev_step, cases_tot[county][i+1:i+1+window], hosp_total[county][i+1:i+1+window], deaths_tot[county][i+1:i+1+window], [phi for p in range(window)], window),method='leastsq',**{'max_nfev':10000,'epsfcn':1e-7})
        
        # evaluate model
        for p in range(window):
            if len(PD['S']) < len(cases_tot[county]):
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



def get_states(county='Los Angeles',freq=1/7,sens=0.8,spec=0.8,wind=1,thresh=0.0001,p=2.5,window=7,PD=None):
    SD = {'S':[],'Qu':[],'I':[],'D':[],'H':[],'R':[],'QR':[],'E':[],'tot':[],'total':[]}
    if PD is None:
        PD = get_PD(county=county,window=window)
    for i in range(len(PD['S'])-1):
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
            prev_step = [S0, Qu0, I0, D0, H0, R0, QR0, E0, tot0, total0]
            
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
        params.add('alpha' ,   np.mean(PD['alpha'][i-wind+1:i+1]) if i >= wind else np.mean(PD['alpha'][:i+1] ) , vary=False)
        params.add('eps'   ,   0               , vary=False)
        params.add('nu'    ,   np.mean(PD['nu'][i-wind+1:i+1])    if i >= wind else np.mean(PD['nu'][:i+1]    ) , vary=False)
        params.add('mu'    ,   np.mean(PD['mu'][i-wind+1:i+1])    if i >= wind else np.mean(PD['mu'][:i+1]    ) , vary=False)
        params.add('tau'   ,   np.mean(PD['tau'][i-wind+1:i+1])   if i >= wind else np.mean(PD['tau'][:i+1]   ) , vary=False)
        params.add('freq'  ,   freq            , vary=False)
        params.add('sens'  ,   sens            , vary=False)
        params.add('spec'  ,   spec            , vary=False)

        # evaluate model
        S,Qu,I,D,H,R,QR,E,tot,total = funct(prev_step,params,thresh=thresh,p=p)

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

print("Fitting to County Data...")
window=7
county_fits = {}
for county in cases_tot.keys():
    print(county)
    county_fits[county] = get_PD(county=county,window=window)
    df = pd.DataFrame.from_dict(county_fits[county], orient="index")
    df.to_csv(csvOutDir + "/" + "county_fit_time_series_{}.csv".format(county))
print("Done!")


'''
for county in county_fits.keys():
    # Save the plots of this 
    PD = county_fits[county]
    SD = get_states(county=county,freq=0,sens=0,spec=0)
    fig, (ax5,ax6,ax7) = plt.subplots(1, 3)
    fig.set_size_inches(28.0, 10.5)
    x = np.linspace(0,len(PD['tot']), len(PD['tot']))
    ax5.plot(cases_tot[county],color='black',label="Data")
    ax5.plot(PD['tot'],linestyle='dashdot',color='black',label="Simulation with parameters changing daily")
    ax5.set_xlabel("Time Step")
    ax5.set_ylabel("Fraction of Population")
    ax5.legend(title=county)
    ax6.plot(hosp_total[county],color='red',label="Data")
    ax6.plot(PD['H'],linestyle='dashdot',color='red',label="Simulation with parameters changing daily")
    ax6.set_xlabel("Time Step")
    ax6.set_ylabel("Fraction of Population")
    ax6.legend(title=county)
    ax7.plot(deaths_tot[county],color='green',label="Data")
    ax7.plot(PD['E'],linestyle='dashdot',color='green',label="Simulation with parameters changing daily")
    ax7.set_xlabel("Time Step")
    ax7.set_ylabel("Fraction of Population")
    ax7.legend(title=county)
    pdf.savefig(fig,bbox_inches="tight") 
    # plt.show()
    plt.close(fig)

    fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4)
    fig.set_size_inches(28.0, 10.5)
    x = np.linspace(0,len(PD['mu']), len(PD['mu']))
    ax1.plot(x,PD['alpha'],label="Time Varying alpha")
    ax1.set_xlabel(r"$\alpha$")
    ax1.set_ylabel("Value")
    ax1.legend(title=county)
    ax2.plot(x,PD['nu'],label="Time Varying nu")
    ax2.set_xlabel(r"$\nu$")
    ax2.set_ylabel("Value")
    ax2.legend(title=county)
    ax3.plot(x,PD['mu'],label="Time Varying mu")
    ax3.set_xlabel(r"$\mu$")
    ax3.set_ylabel("Value")
    ax3.legend(title=county)
    ax4.plot(x,PD['tau'],label="Time Varying tau")
    ax4.set_xlabel(r"$\tau$")
    ax4.set_ylabel("Value")
    ax4.legend(title=county)
    pdf.savefig(fig,bbox_inches="tight") 
    # plt.show()
    plt.close(fig)

    fig, (ax8,ax9,ax10) = plt.subplots(1, 3)
    fig.set_size_inches(28.0, 10.5)
    x = np.linspace(0,len(PD['I']), len(PD['I']))
    # z = np.linspace(0,len(cases_tot), len(cases_tot))
    # ax8.plot(x,PD['I'],linestyle='dashdot',color='black', label="Simulation state I with parameters changing daily")
    ax8.plot(PD['D'],linestyle='dashdot',color='red'  , label="Simulation state D with parameters changing daily")
    # ax8.plot(x,PD['R'],linestyle='dashdot',color='green', label="Simulation state R with parameters changing daily")
    # ax8.plot(x,PD['QR'],linestyle='dashdot',color='blue', label="Simulation state QR with parameters changing daily")
    ax8.set_xlabel("Time Step")
    ax8.set_ylabel("Fraction of Population")
    ax8.legend(title=county)
    ax9.plot(np.array(PD['I'])/np.array(PD['D']),linestyle='dashdot',color='red',label="Time Varying I/D")
    ax9.set_xlabel("Time Step")
    ax9.set_ylabel(r"\Phi")
    ax9.legend(title=county)
    ax10.plot(PD['S'],linestyle='dashdot',color='black', label="Simulation state S with parameters changing daily")
    ax10.plot(PD['Qu'],linestyle='dashdot',color='blue', label="Simulation state Qu with parameters changing daily")
    ax10.set_xlabel("Time Step")
    ax10.set_ylabel("Fraction of Population")
    ax10.legend(title=county)
    # plt.show()
    pdf.savefig(fig,bbox_inches="tight")  
    plt.close(fig)

pdf.close()
'''

# Convert to percent
toPerc = 100
# The fitting procedure cuts off some of the last elements if there is not a full window length at the end. This number is how much is cut off
len_diff = len(cases_tot['Alameda']) - len(county_fits['Alameda']['S'])
# zero array the length of the simulation 
nil = np.zeros(len(county_fits['Alameda']['S']))
# Total county population
POP = np.sum([pops[county] for county in pops.keys()])
# Total cases, hospitalized, dead from data
print("Computing total cases, hospitalizations, and cases from data...")
CASES, HOSP, DEATHS = nil, nil, nil
for county in cases_tot.keys():
    CASES  = CASES  + np.array(cases_tot[county][:-len_diff]  if len_diff > 0 else cases_tot[county])*pops[county] # ensure that the lengths are the same
    HOSP   = HOSP   + np.array(hosp_total[county][:-len_diff] if len_diff > 0 else hosp_total[county])*pops[county]
    DEATHS = DEATHS + np.array(deaths_tot[county][:-len_diff] if len_diff > 0 else deaths_tot[county])*pops[county]
CASES   =  CASES / POP * toPerc
HOSP    =  HOSP / POP * toPerc
DEATHS  =  DEATHS / POP * toPerc
print("Done!")


def get_map_states(series = {}, freqs=[], enable_protocol=0, sens=0.8, spec=0.9, p=1):
    for freq in freqs:
        print("Frequency: {}/{}".format(list(freqs).index(freq),len(freqs)))
        if freq not in series.keys():
            series[freq] = {}
        count = {}
        for county in county_fits.keys():
            if enable_protocol == -1: 
                SD = get_states(county=county,freq=enable_protocol,sens=sens,spec=spec,thresh=freq,p=p,PD=county_fits[county])
            else: 
                SD = get_states(county=county,freq=freq,sens=sens,spec=spec,PD=county_fits[county])
            count[county] = SD
        series[freq] = count

def get_protocol(scan={}, county_series={}, freqs=[], enable_protocol=0, sens=0.8, spec=0.9, p=1):
    # Sum up all of the county state values
    nil = np.zeros(len(county_fits['Alameda']['S']))
    for freq in freqs:
        print("Frequency: {}/{}".format(list(freqs).index(freq),len(freqs)))
        # Compute the total California states by summing and normalizing over counties
        if freq not in scan.keys():
            scan[freq] = {}
            county_series[freq] = {}

        s = {'S':nil.copy(),'Qu':nil.copy(),'I':nil.copy(),'D':nil.copy(),'I+D':nil.copy(),'H':nil.copy(),'R':nil.copy(),'QR':nil.copy(),'E':nil.copy(),'tot':nil.copy(),'total':nil.copy(),'cost':nil.copy()}
        ser = {}
        for county in county_fits.keys():
            if enable_protocol == -1:
                SD = get_states(county=county,freq=enable_protocol,sens=sens,spec=spec,thresh=freq,p=p,PD=county_fits[county])
            else:
                SD = get_states(county=county,freq=freq,sens=sens,spec=spec,PD=county_fits[county])
            s['S']     =  s['S']     + np.array(SD['S'])  * pops[county]    
            s['Qu']    =  s['Qu']    + np.array(SD['Qu']) * pops[county]    
            s['I']     =  s['I']     + np.array(SD['I'])  * pops[county]    
            s['D']     =  s['D']     + np.array(SD['D'])  * pops[county]    
            s['I+D']   =  s['I+D']   + (np.array(SD['D']) + np.array(SD['I'])) * pops[county] 
            s['H']     =  s['H']     + np.array(SD['H'])  * pops[county]    
            s['R']     =  s['R']     + np.array(SD['R'])  * pops[county]
            s['QR']    =  s['QR']    + np.array(SD['QR']) * pops[county]    
            s['E']     =  s['E']     + np.array(SD['E'])  * pops[county]    
            s['tot']   =  s['tot']   + np.array(SD['tot']) * pops[county]   
            s['total'] =  s['total'] + np.array(SD['total']) * pops[county] 
            ser[county] = SD
            
        # Convert to fraction of total CA population
        for key in s.keys():
            s[key] /= POP
            
        # Make sure total population is 1.0 as safekeeping
        for i in range(len(s['S'])):
            pop = 0
            for key in s.keys():
                pop += s[key][i]
            if round(pop) != 1.0: print(i,round(pop))
                
        # Convert to percent
        for key in s.keys():
            s[key] *= 100 

        scan[freq] = s
        county_series[freq] = ser

# sets range for frequency in testing protocol. See funct(...) in fitting area.
p = 1

# Run in parallel 
starttime = time.time()
processes = []
save = {}
manager = multiprocessing.Manager()


#print("Running Map Diagram for County Protocol...")
map_freqs_rapid_county = [5.00000000e-04, 3.53553391e-04, 2.50000000e-04, 1.76776695e-04, 1.25000000e-04, 8.83883476e-05, 6.25000000e-05, 4.41941738e-05, 3.12500000e-05, 2.20970869e-05, 1.56250000e-05, 1.10485435e-05, 7.81250000e-06]
enable_protocol = -1 # -1 turns on
sens=0.8
spec=0.9
map_series_rapid_county = manager.dict()
proc = multiprocessing.Process(target=get_map_states, args=(map_series_rapid_county, map_freqs_rapid_county, enable_protocol, sens, spec, p))
processes.append(proc)
proc.start()

#print("Running Map Diagram for PCR Uniform...")
map_freqs_pcr_uniform = [1./i for i in range(1,8)]
enable_protocol = 0 # -1 turns on
sens=1.0
spec=1.0
map_series_pcr_uniform = manager.dict()
proc = multiprocessing.Process(target=get_map_states, args=(map_series_pcr_uniform, map_freqs_pcr_uniform, enable_protocol, sens, spec, p))
processes.append(proc)
proc.start()


#print("Running County Based Testing Protocol...")
freqs = np.logspace(-17,-6.5,30,base=2) # user set threshold
enable_protocol = -1 # -1 turns on
sens=0.8
spec=0.9
county_scan, county_series = manager.dict(),manager.dict()
proc = multiprocessing.Process(target=get_protocol, args=(county_scan, county_series, freqs, enable_protocol, sens, spec, p))
processes.append(proc)
proc.start()


#print("Running PCR Nonuniform Testing Protocol...")
freqs_pcr = np.logspace(42,60,30,base=2) # thresholds for PCR county based
enable_protocol = -1 # -1 turns on
sens=1.0
spec=1.0
pcr_scan, pcr_series = manager.dict(),manager.dict()
proc = multiprocessing.Process(target=get_protocol, args=(pcr_scan, pcr_series, freqs_pcr, enable_protocol, sens, spec, p))
processes.append(proc)
proc.start()


#print("Running County Based Testing Uniform Protocol...")
freqs_county_based_uniform = np.linspace(0.02,1,30)
enable_protocol = 0
sens=0.8
spec=0.9
scan_uniform, county_series_uniform = manager.dict(),manager.dict()
proc = multiprocessing.Process(target=get_protocol, args=(scan_uniform, county_series_uniform, freqs_county_based_uniform, enable_protocol, sens, spec, p))
processes.append(proc)
proc.start()

#print("Running PCR Uniform Testing Protocol...")
freqs_pcr_based_uniform = np.linspace(0.02,1,30)
enable_protocol = 0
sens=1.0
spec=1.0
pcr_scan_uniform, pcr_series_uniform = manager.dict(),manager.dict()
proc = multiprocessing.Process(target=get_protocol, args=(pcr_scan_uniform, pcr_series_uniform, freqs_pcr_based_uniform, enable_protocol, sens, spec, p))
processes.append(proc)
proc.start()


print("Running {} Jobs...".format(len(processes)))
for process in processes:    
    process.join()       

save["map_series_rapid_county.csv"] = map_series_rapid_county
save["map_series_pcr_uniform.csv"]  = map_series_pcr_uniform
save["county_based_nonuniform_scan.csv"] = county_scan
save["county_based_nonuniform_series.csv"] = county_series
save["pcr_based_nonuniform_scan.csv"] = pcr_scan
save["pcr_based_nonuniform_series.csv"] = pcr_series
save["county_based_uniform_scan.csv"] = scan_uniform
save["county_based_uniform_series.csv"] = county_series_uniform
save["pcr_based_uniform_scan.csv"] = pcr_scan_uniform
save["pcr_based_uniform_series.csv"] = pcr_series_uniform
for file, dic in save.items():
    df = pd.DataFrame.from_dict(dic, orient="index")
    df.to_csv(csvOutDir + "/" + file)

print('That took {} seconds'.format(time.time() - starttime))
sys.exit()
