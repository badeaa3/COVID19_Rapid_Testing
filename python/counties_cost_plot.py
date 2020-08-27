import pandas as pd
import ast
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FuncFormatter
from matplotlib.lines import Line2D
from county_utilities import *
import sys, os, csv
sys.path.insert(1, '../data_processing')
from load_CA_data import *

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

# Total county population
POP = np.sum([pops[county] for county in pops.keys()])

# Check how much of population is used
print("Percentage of California population used: {}".format(total_included_pop/region_pop))

# Load the data
def load(scanName = "", seriesName=""):
    df = pd.read_csv("county_output/" + scanName, index_col=0)
    d = df.to_dict("split")
    scan = {}
    i = 0
    for key in d['index']:
        if key not in scan.keys():
            scan[key] = {}
        for state, data in zip(d['columns'],d['data'][i]):
            temp = data.split("[")[-1].split("]")[0].split(" ")
            l = []
            for tem in temp:
                try:
                    l.append(float(tem))
                except:
                    try:
                        l.append(float(tem.split("\n")[0]))
                    except:
                        continue
            scan[key][state] = l
        i += 1

    df = pd.read_csv("county_output/" + seriesName, index_col=0)
    d = df.to_dict("split")
    series = {}
    i = 0
    for key in d['index']:
        if key not in series.keys():
            series[key] = {}
        for county, data in zip(d['columns'],d['data'][i]):
            series[key][county] = ast.literal_eval(data) 
        i += 1

    freqs = d['index']
    return freqs, scan, series

county_freqs, county_scan, county_series = load(scanName="county_based_nonuniform_scan.csv", seriesName="county_based_nonuniform_series.csv")
pcr_freqs, pcr_scan, pcr_series = load(scanName="pcr_based_nonuniform_scan.csv", seriesName="pcr_based_nonuniform_series.csv")
county_freqs_uniform, county_scan_uniform, county_series_uniform = load(scanName="county_based_uniform_scan.csv", seriesName="county_based_uniform_series.csv")
pcr_freqs_uniform, pcr_scan_uniform, pcr_series_uniform = load(scanName="pcr_based_uniform_scan.csv", seriesName="pcr_based_uniform_series.csv")

# Print out some useful statistics 
print(np.array(county_scan[0.0005437278755471922]['I']) + np.array(county_scan[0.0005437278755471922]['D']))
print(np.array(county_scan_uniform[0.15517241379310345]['I']) + np.array(county_scan_uniform[0.15517241379310345]['D']))

outDir = 'county_output/Figure5'
sens = 0.8
spec = 0.9
toPerc = 100
p = 1
scale_pcr_total = 1.0
scale_pcr_H = 1.0
scale_pcr_E = 1.0
scale_pcr_cost = 40.
per_100k = 10**5 / toPerc # the state values are outputted as percents. CHANGE IF OTHER CODE IS CHANGED.
county_variable = 'County Based Variable Testing'
county_variable_pcr = 'County Based Variable Testing PCR\n(Cost divided by {})'#.format(scale_pcr)
uniform_pcr = 'PCR Based Uniform Testing (Cost divided by {})'.format(int(scale_pcr_cost))
uniform = 'County Based Uniform Testing'
xlabel = 'Cost per Person per Day ($)'
ylabel = 'People per 100k'
cost_min,cost_max = 0.5,2.5
cumulative_min, cumulative_max = 1.0, 1.75
hosp_min, hosp_max = 0.01448, 0.0148
dead_min, dead_max = 0.0045, 0.007
xlabelpad = 15
ylabelpad = 15


fig = plt.figure(figsize=(40,10))
ax0 = fig.add_subplot(131)
plt.suptitle('Sensitivity: {}%, Specificity: {}%'.format(str(np.around(100*sens,1)),np.around(100*spec,1)),fontsize=22,y=1.05)
fig,ax0 = make_fig(fig,ax0)

cp = []
price = 5
tp = []


for freq, sts in county_scan.items():
    length = len(sts['total'])
    tp.append(per_100k*sts['total'][-1])
    cost = 0
    for county in county_series[freq]:
        thresh = freq
        for i in range(len(county_series[freq][county]['S'])):
            D = county_series[freq][county]['D'][i]
            if D <= thresh:
                f = 1/(p*np.log2(thresh/(D))+1) #1/np.round(p*np.log2(thresh/D)+1,0) # D = thresh*(1/2)^((a-1)/2)
            else:
                f = 1.0
            cost += f*price*(county_series[freq][county]['S'][i]+county_series[freq][county]['I'][i]+county_series[freq][county]['R'][i])*(pops[county])
    cp.append(cost/(length*POP))
ax0.plot(cp,tp,color='black',linewidth=2.5,zorder=1,label=county_variable)
ax0.set_xlabel(xlabel)
ax0.set_ylabel(ylabel)
sc = ax0.scatter(cp,tp,marker='o',c = [np.log(f) for f in county_freqs],cmap = 'RdPu',zorder=2)


# cp = []
# price = 5
# tp = []

# for freq, sts in pcr_scan.items():
#     tp.append(per_100k*sts['total'][-1]*scale_pcr_total)
#     cost = 0
#     for county in pcr_series[freq]:
#         thresh = freq
#         for i in range(len(pcr_series[freq][county]['S'])):
#             D = pcr_series[freq][county]['D'][i]
#             if D <= thresh:
#                 f = 1/(p*np.log2(thresh/(D))+1) #1/np.round(p*np.log2(thresh/D)+1,0) # D = thresh*(1/2)^((a-1)/2)
#             else:
#                 f = 1.0
#             cost += f*100*(pcr_series[freq][county]['S'][i]+pcr_series[freq][county]['I'][i]+pcr_series[freq][county]['R'][i])*(pops[county])
#     cp.append(cost/(scale_pcr_cost*length*POP))
# print(cp)
# print(tp)
# ax0.plot(cp,tp,color='magenta',linewidth=2.5,label=county_variable_pcr.format(scale_pcr_total))



cp = []
price = 5
tp = []

for freq, sts in pcr_scan_uniform.items():
    tp.append(per_100k*sts['total'][-1]*scale_pcr_total)
    cost = 0
    for county in pcr_series_uniform[freq]:
        thresh = freq
        for i in range(len(pcr_series_uniform[freq][county]['S'])):
            D = pcr_series_uniform[freq][county]['D'][i]
            f = freq
            cost += f*100*(pcr_series_uniform[freq][county]['S'][i]+pcr_series_uniform[freq][county]['I'][i]+pcr_series_uniform[freq][county]['R'][i])*(pops[county])
    cp.append(cost/(scale_pcr_cost*length*POP))
ax0.plot(cp,tp,color='pink',linewidth=2.5,label=uniform_pcr.format(scale_pcr_total))

cp = []
price = 5
tp = []

for freq, sts in county_scan_uniform.items():
    tp.append(per_100k*sts['total'][-1])
    cost = 0
    for county in county_series_uniform[freq]:
        thresh = freq
        for i in range(len(county_series_uniform[freq][county]['S'])):
            D = county_series_uniform[freq][county]['D'][i]
            f= freq
            cost += f*price*(county_series_uniform[freq][county]['S'][i]+county_series_uniform[freq][county]['I'][i]+county_series_uniform[freq][county]['R'][i])*(pops[county])
    cp.append(cost/(length*POP))
ax0.plot(cp,tp,color='grey',linewidth=2.5,label=uniform)
ax0.set_xlim(cost_min,cost_max)
ax0.set_ylim(per_100k*cumulative_min,per_100k*cumulative_max)
ax0.get_yaxis().set_major_formatter(FuncFormatter(lambda x, p: human_format(x)))
ax0.xaxis.labelpad = xlabelpad
ax0.yaxis.labelpad = ylabelpad
ax0.set_title('Cumulative Infections',fontsize = 22,fontweight='bold',x=0.0,horizontalalignment='left')
#ax0.legend(title_fontsize=22,fontsize=15,ncol=1,loc="upper right",facecolor='white', framealpha=0.8,edgecolor='None')
plt.grid()

ax0 = fig.add_subplot(132)
fig,ax0 = make_fig(fig,ax0)
cp = []
price = 5
tp = []

for freq, sts in county_scan.items():
    tp.append(per_100k*np.max(sts['H']))
    cost = 0
    for county in county_series[freq]:
        thresh = freq
        for i in range(len(county_series[freq][county]['S'])):
            D = county_series[freq][county]['D'][i]
            if D <= thresh:
                f = 1/(p*np.log2(thresh/(D))+1) #1/np.round(p*np.log2(thresh/D)+1,0) # D = thresh*(1/2)^((a-1)/2)
            else:
                f = 1.0
            cost += f*price*(county_series[freq][county]['S'][i]+county_series[freq][county]['I'][i]+county_series[freq][county]['R'][i])*(pops[county])
    cp.append(cost/(POP*len(sts['H'])))

ax0.plot(cp,tp,color='black',linewidth=2.5,zorder=1,label=county_variable)
ax0.set_xlabel(xlabel)
ax0.set_ylabel(ylabel)
sc = ax0.scatter(cp,tp,marker='o',c = [np.log(f) for f in county_freqs],cmap = 'RdPu',zorder=2)


# cp = []
# price = 5
# tp = []

# for freq, sts in pcr_scan.items():
#     tp.append(per_100k*np.max(sts['H'])*scale_pcr_H)
#     cost = 0
#     for county in pcr_series[freq]:
#         thresh = freq
#         for i in range(len(pcr_series[freq][county]['S'])):
#             D = pcr_series[freq][county]['D'][i]
#             if D <= thresh:
#                 f = 1/(p*np.log2(thresh/(D))+1) #1/np.round(p*np.log2(thresh/D)+1,0) # D = thresh*(1/2)^((a-1)/2)
#             else:
#                 f = 1.0
#             cost += f*100*(pcr_series[freq][county]['S'][i]+pcr_series[freq][county]['I'][i]+pcr_series[freq][county]['R'][i])*(pops[county])
#     cp.append(cost/(scale_pcr_cost*length*POP))  
# ax0.plot(cp,tp,color='magenta',linewidth=2.5,label=county_variable_pcr.format(scale_pcr_H))


cp = []
price = 5
tp = []

for freq, sts in pcr_scan_uniform.items():
    tp.append(per_100k*np.max(sts['H'])*scale_pcr_H)
    cost = 0
    for county in pcr_series_uniform[freq]:
        thresh = freq
        for i in range(len(pcr_series_uniform[freq][county]['S'])):
            D = pcr_series_uniform[freq][county]['D'][i]
            f = freq
            cost += f*100*(pcr_series_uniform[freq][county]['S'][i]+pcr_series_uniform[freq][county]['I'][i]+pcr_series_uniform[freq][county]['R'][i])*(pops[county])
    cp.append(cost/(scale_pcr_cost*length*POP))    
ax0.plot(cp,tp,color='pink',linewidth=2.5,label=uniform_pcr.format(scale_pcr_H))

cp = []
price = 5
tp = []

for freq, sts in county_scan_uniform.items():
    tp.append(per_100k*np.max(sts['H']))
    cost = 0
    for county in county_series_uniform[freq]:
        thresh = freq
        for i in range(len(county_series_uniform[freq][county]['S'])):
            D = county_series_uniform[freq][county]['D'][i]
            f= freq
            cost += f*price*(county_series_uniform[freq][county]['S'][i]+county_series_uniform[freq][county]['I'][i]+county_series_uniform[freq][county]['R'][i])*(pops[county])
    cp.append(cost/(len(sts['H'])*POP))

ax0.plot(cp,tp,color='grey',linewidth=2.5,label=uniform)
ax0.set_xlim(cost_min,cost_max)
ax0.set_ylim(per_100k*hosp_min,per_100k*hosp_max)
ax0.xaxis.labelpad = xlabelpad
ax0.yaxis.labelpad = ylabelpad
ax0.set_title('Max Simultaneous Hospitalizations',fontsize = 22,fontweight='bold',x=0.0,horizontalalignment='left')
#ax0.legend(title_fontsize=22,fontsize=15,ncol=1,loc="upper right",facecolor='white', framealpha=0.8,edgecolor='None')
plt.grid()


ax0 = fig.add_subplot(133)
fig,ax0 = make_fig(fig,ax0)
cp = []
price = 5
tp = []
for freq, sts in county_scan.items():
    tp.append(per_100k*sts['E'][-1])
    cost = 0
    for county in county_series[freq]:
        thresh = freq
        for i in range(len(county_series[freq][county]['S'])):
            D = county_series[freq][county]['D'][i]
            if D <= thresh:
                f = 1/(p*np.log2(thresh/(D))+1) #1/np.round(p*np.log2(thresh/D)+1,0) # D = thresh*(1/2)^((a-1)/2)
            else:
                f = 1.0
            cost += f*price*(county_series[freq][county]['S'][i]+county_series[freq][county]['I'][i]+county_series[freq][county]['R'][i])*(pops[county])
    cp.append(cost/(length*POP))
     
l3 = ax0.plot(cp,tp,color='black',linewidth=2.5,zorder=1,label=county_variable.format(scale_pcr_total))
ax0.set_xlabel(xlabel)
ax0.set_ylabel(ylabel)
sc = ax0.scatter(cp,tp,marker='o',c = [np.log(f) for f in county_freqs],cmap = 'RdPu',zorder=2)


# cp = []
# price = 5
# tp = []

# for freq, sts in pcr_scan.items():
#     tp.append(per_100k*sts['E'][-1]*scale_pcr_E)
#     cost = 0
#     for county in pcr_series[freq]:
#         thresh = freq
#         for i in range(len(pcr_series[freq][county]['S'])):
#             D = pcr_series[freq][county]['D'][i]
#             if D <= thresh:
#                 f = 1/(p*np.log2(thresh/(D))+1) #1/np.round(p*np.log2(thresh/D)+1,0) # D = thresh*(1/2)^((a-1)/2)
#             else:
#                 f = 1.0
#             cost += f*100*(pcr_series[freq][county]['S'][i]+pcr_series[freq][county]['I'][i]+pcr_series[freq][county]['R'][i])*(pops[county])
#     cp.append(cost/(scale_pcr_cost*length*POP))
# ax0.plot(cp,tp,color='magenta',linewidth=2.5,label=county_variable_pcr.format(scale_pcr_E))


cp = []
price = 5
tp = []

for freq, sts in pcr_scan_uniform.items():
    tp.append(per_100k*sts['E'][-1]*scale_pcr_E)
    cost = 0
    for county in pcr_series_uniform[freq]:
        thresh = freq
        for i in range(len(pcr_series_uniform[freq][county]['S'])):
            D = pcr_series_uniform[freq][county]['D'][i]
            f = freq
            cost += f*100*(pcr_series_uniform[freq][county]['S'][i]+pcr_series_uniform[freq][county]['I'][i]+pcr_series_uniform[freq][county]['R'][i])*(pops[county])
    cp.append(cost/(scale_pcr_cost*length*POP))
l2 = ax0.plot(cp,tp,color='pink',linewidth=2.5,label=uniform_pcr.format(scale_pcr_E))

cp = []
price = 5
tp = []
for freq, sts in county_scan_uniform.items():
    tp.append(per_100k*np.max(sts['E']))
    cost = 0
    for county in county_series_uniform[freq]:
        thresh = freq
        for i in range(len(county_series_uniform[freq][county]['S'])):
            D = county_series_uniform[freq][county]['D'][i]
            f= freq
            cost += f*price*(county_series_uniform[freq][county]['S'][i]+county_series_uniform[freq][county]['I'][i]+county_series_uniform[freq][county]['R'][i])*(pops[county])
    cp.append(cost/(length*POP))

l1 = ax0.plot(cp,tp,color='grey',linewidth=2.5,label=uniform)
ax0.xaxis.labelpad = xlabelpad
ax0.yaxis.labelpad = ylabelpad
ax0.set_title('Deceased',fontsize = 22,fontweight='bold',x=0.0,horizontalalignment='left')
ax0.set_xlim(cost_min,cost_max)
ax0.set_ylim(per_100k*dead_min,per_100k*dead_max)
# ax0.legend(title_fontsize=22,fontsize=20,ncol=1,loc="upper right",facecolor='white', framealpha=0.8,edgecolor='None')

plt.grid()


c1 = l1[0].get_color()
c2 = l2[0].get_color()
c3 = l3[0].get_color()
l4 = Line2D([0],[0],color=c2,lw=4)
l5 = Line2D([0],[0],color=c1,lw=4)
l6 = Line2D([0],[0],color=c3,lw=4)
leg1 = fig.legend(handles=[l6,l5,l4],labels=[county_variable.format(scale_pcr_total),uniform,uniform_pcr.format("{}, {}, {}".format(scale_pcr_total,scale_pcr_H,scale_pcr_E))],fontsize=20,loc='upper center',bbox_to_anchor = (0.416,0.99),facecolor='white',ncol=3,edgecolor='None')

fig.add_artist(leg1)
        
# Formats the colobar 
def get_thresh(t,pos):
    j = per_100k*np.exp(t)
    thresh = ""
    if j < 1.0: thresh = "{:0.2f}".format(j)
    if j >= 1.0 and j < 10.0: thresh = "{:1.1f}".format(j)
    if j >= 10.0 and j < 100.0: thresh = "{:2.0f}".format(j)
    if j >= 100.0 and j < 1000.0: thresh = "{:3.0f}".format(j)
    return thresh

cbar_ax = fig.add_axes([0.93, 0.15, 0.02, 0.7])
cbar = plt.colorbar(sc,format=FuncFormatter(get_thresh),cax=cbar_ax)
cbar.ax.tick_params(labelsize=20) 
#cbar.ax.set_ylabel('Active Detected Infected Threshold (%)',fontsize=20,labelpad=-170)
cbar.ax.set_ylabel('Detected per 100k Threshold',fontsize=20,labelpad=-150)
plt.savefig(outDir+'/County_Cost_plot'+'.pdf',bbox_inches='tight')
#plt.show()


