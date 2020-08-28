import chart_studio.plotly as py
import plotly.figure_factory as ff
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd
import ast
import numpy as np
import sys, os, csv
sys.path.insert(1, '../data_processing')
from load_CA_data import *
from county_utilities import *
import time

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

# Load the map data
map_series_rapid_county = "map_series_rapid_county.csv"
map_series_pcr_uniform = "map_series_pcr_uniform.csv" 
inFileName = map_series_rapid_county
df = pd.read_csv("county_output/" + inFileName,index_col=0)
d = df.to_dict("split")
map_series = {}
i = 0
for key in d['index']:
    if key not in map_series.keys():
        map_series[key] = {}
    for county, data in zip(d['columns'],d['data'][i]):
        map_series[key][county] = ast.literal_eval(data) 
    i += 1

# Compute the infected populations
df_sample = pd.read_csv('https://raw.githubusercontent.com/plotly/datasets/master/minoritymajority.csv')
df_sample_r = df_sample[df_sample['STNAME'] == 'California']
ctyname = df_sample_r['CTYNAME'].tolist()
fips = df_sample_r['FIPS'].tolist()
# Issue: create_choropleth had issues displaying such small bins values. 
# Solution: scale all values up by this value and then create custom legend to display the corresponding percents
legend_hack_scaling = 10**5 
steps = {}
for freq, protocol in map_series.items():
    length = len(map_series[freq]['Los Angeles']['S'])
    if freq not in steps:
        steps[freq] = {}
    for i in range(length):
        infected = [0]*len(fips)
        for county, series in protocol.items():
            indx = [idx for idx, s in enumerate(ctyname) if county in s][0]
            infected[indx] = series['D'][i] * legend_hack_scaling #(series['D'][i]+series['I'][i])*toPerc 
        for county in list(cd.keys()):
            if county not in list(protocol.keys()):
                indx = [idx for idx, s in enumerate(ctyname) if county in s][0]
                infected[indx] = -999.0 # this indicates a county was not included in the study, i.e. that it recieves no rapid tests
        steps[freq][i] = infected


# In[ ]:


# Get the binning_endpoints and assign the color codes

# in increasing order - higher number of cases
# equals higher index color
not_included = ['#808080']
overflow = ['#000000']
colors = ['#4D006B','#660071','#830078','#98007B','#B20680','#c50f8a','#e23a9a','#eb5295','#f562a0','#f773a7','#f98fb0','#faa4b7','#fabfbf','#fdf0dc'] 
colorscale = not_included + list(reversed(colors)) # + overflow

buckets = {} # bucket[i] = upper bound on ith bucket
buckets_freqs = {} # freqs[i] = frequency in ith bucket
for thresh in steps.keys():
    buck = []
    bf = []
    D = thresh
    a = 1
    while D >= 0.00001:
        D = thresh*(1/2)**((a-1)/2)
        buck.append(D*legend_hack_scaling)
        a += 1
        bf.append((1/np.round((2*np.log2(thresh/D)+1),0)))
    if "pcr" in inFileName:
        buckets[thresh] = [50.0, 35.35533905932738, 25.0, 17.67766952966369, 12.5, 8.838834764831844, 6.25, 4.419417382415922, 3.125, 2.209708691207961, 1.5625, 1.1048543456039805, 0.78125]
        buckets_freqs[thresh] = [1.0, 0.5, 0.3333333333333333, 0.25, 0.2, 0.16666666666666666, 0.14285714285714285, 0.125, 0.1111111111111111, 0.1, 0.09090909090909091, 0.08333333333333333, 0.07692307692307693]
    else: 
        buckets[thresh] = buck
        buckets_freqs[thresh] = bf

    D = thresh
    a = 1
 
binning_endpoints = {}
legends = {}
days = {}
base_freq = list(buckets.keys())[0] #0.0005
toPerc = 100
per_100k = 10**5
for thresh in buckets.keys():
    binning_endpoints[thresh] = [0] + list(reversed(buckets[thresh]))
    '''
    # Legend displayed as percent
    binning_endpoints_legend = list( np.array(binning_endpoints[thresh]) / legend_hack_scaling * toPerc ) # the true binning endpoints in percent used for the legend
    legend = []
    for j in binning_endpoints_legend:
        if j < 1e-3: legend.append('{:0=1.2f}e-3'.format(j/1e-3))
        elif j < 1e-2 and j > 1e-3: legend.append('{:0=1.2f}e-2'.format(j/1e-2))
        elif j < 1e-1 and j > 1e-2: legend.append('{:0=1.2f}e-1'.format(j/1e-1))
        else: legend.append('{:0=1.3f}'.format(j))
    legends[thresh] = legend
    # Additional formatting
    legends[thresh] += legends[base_freq][len(legends[thresh]):-1]
    legend += ["1.00e-0"]
    '''
    # Legend displayed as infections per 100k
    binning_endpoints_legend = list( np.array(binning_endpoints[thresh]) / legend_hack_scaling * per_100k ) # the true binning endpoints in percent used for the legend
    legend = []
    for j in binning_endpoints_legend:
        if j < 1.0: legend.append("{:0.1f}".format(round(j,1)))
        if j >= 1.0 and j < 10.0: legend.append("{:1.1f}".format(j))
        if j >= 10.0 and j < 100.0: legend.append("{:2.0f}".format(j))# + ".")
        # if j < 1e-3: legend.append('{:0=1.2f}e-3'.format(j/1e-3))
        # elif j < 1e-2 and j > 1e-3: legend.append('{:0=1.2f}e-2'.format(j/1e-2))
        # elif j < 1e-1 and j > 1e-2: legend.append('{:0=1.2f}e-1'.format(j/1e-1))
        # else: legend.append('{:0=1.3f}'.format(j))
    legends[thresh] = legend
    # Additional formatting
    legends[thresh] += legends[base_freq][len(legends[thresh]):-1]
    hm =  human_format(0.01*per_100k)
    legend += [hm] #[hm[:1] + "." + hm[1:]]
    #print(legend)
    
    day = list(reversed((1./np.array(buckets_freqs[thresh])).astype(int)))
    day = ["No Testing Done"] + [k+1 for k in day] + [1]
    days[thresh] = day
    days[thresh] += [1] * (len(days[base_freq]) - len(days[thresh]))
    if "uniform" in inFileName:
        days[thresh] = ["No Testing Done"] + [1./thresh] * (len(days[thresh]) - 1)

    # print("Binning endpoints and corresponding days between tests:")
    # print(list(zip(legend,day)))
    # print(len(day),len(legend),len(colorscale))

# Calculate the costs per person 

map_cp = {}
price = 5
p = 1

for thresh in map_series:
    cost = 0
    length = len(map_series[thresh]['Alameda']['S'])
    for county in map_series[thresh]:
        if "uniform" in inFileName:
            for i in range(len(map_series[thresh][county]['S'])):
                cost += thresh*100*(map_series[thresh][county]['S'][i]+map_series[thresh][county]['I'][i]+map_series[thresh][county]['R'][i])*(pops[county])
        else:
            for i in range(len(map_series[thresh][county]['S'])):
                D = map_series[thresh][county]['D'][i]
                if D <= thresh:
                    f = 1/(p*np.log2(thresh/(D))+1) #1/np.round(p*np.log2(thresh/D)+1,0) # D = thresh*(1/2)^((a-1)/2)
                else:
                    f = 1.0
                cost += f*price*(map_series[thresh][county]['S'][i]+map_series[thresh][county]['I'][i]+map_series[thresh][county]['R'][i])*(pops[county])
    map_cp[thresh] = cost/(POP*length)

starttime = time.time()
# Output directory
outDir = 'county_output/Figure6/TimeSeries/' + inFileName.split(".csv")[0]
if not os.path.isdir(outDir):
    os.mkdir(outDir)
place = 'CA_COUNTIES'
panel_maps = {}
panel_times = [0, 30, 60, 90, 105]
for freq, series in steps.items():     
    figs = []
    print(freq)
    for t, infected in series.items():
        if t not in panel_times: continue 
        fig = ff.create_choropleth(
            fips=fips, 
            values=infected, 
            scope=['California'], 
            show_state_data=True,
            colorscale=colorscale,
            binning_endpoints=binning_endpoints[base_freq],
            round_legend_values=False, 
            plot_bgcolor='rgb(255,255,255)',
            paper_bgcolor='rgb(255,255,255)',
            county_outline={'color': 'rgb(0,0,0)', 'width': 0.75},
            exponent_format=False, 
            showlegend=False
        )
        
        # Custom Legend
        x0,y0 = 0.75,0.865
        y0_title_off = 0.1
        normal_txt_off = 0.004
        box_width = 0.02
        box_height = 0.02
        text_x_off = 0.01
        text_y_off = 0.042
        entry_spacing = 0.05
        
        fig.add_annotation(x = x0 - normal_txt_off,y = y0 + y0_title_off + entry_spacing,xanchor="left",yanchor="top",xref='paper',yref='paper',text="Time Step: {}".format(t),showarrow=False,font=dict(size=17),align='center')
        fig.add_annotation(x = x0 - normal_txt_off,y = y0 + y0_title_off,xanchor="left",yanchor="top",xref='paper',yref='paper',text="Detected per 100k | Days Between Tests",showarrow=False,font=dict(size=17),align='center')
        #fig.add_annotation(x = x0 - normal_txt_off,y = y0 + y0_title_off,xanchor="left",yanchor="top",xref='paper',yref='paper',text="Active Detected Infections | Days Between Tests",showarrow=False,font=dict(size=17),align='center')
        
        for i in range(len(legends[freq])):
            color = colorscale[i]
            fig.add_shape(type="rect",xref="paper",yref="paper",x0 = x0, y0 = y0 - i*entry_spacing, x1= x0 + box_width, y1= y0 - i*entry_spacing + box_width,line=dict(color=color,width=3),fillcolor=color)
            if i == 0: text = days[freq][i]
            else: text = "{} - {} | {}".format(legends[freq][i-1],legends[freq][i], days[freq][i])
            #else: text = "{}% - {}% | {}".format(legends[freq][i-1],legends[freq][i], days[freq][i])
            fig.add_annotation(x= x0 + box_width + text_x_off, y= y0 - i*entry_spacing + text_y_off, xanchor="left",yanchor="top",xref='paper',yref='paper',text=text,showarrow=False,font=dict(size=17),align='center')
        if "pcr" in inFileName and "uniform" in inFileName:
            i+=1
            fig.add_annotation(x = x0 - normal_txt_off,y = y0 - i*entry_spacing + text_y_off,xanchor="left",yanchor="top",xref='paper',yref='paper',text="PCR Based Uniform Testing",showarrow=False,font=dict(size=17),align='center')
        else:
            i+=1
            fig.add_annotation(x = x0 - normal_txt_off,y = y0 - i*entry_spacing + text_y_off,xanchor="left",yanchor="top",xref='paper',yref='paper',text="Test Sensitivity: 0.8",showarrow=False,font=dict(size=17),align='center')
            i+=1
            fig.add_annotation(x = x0 - normal_txt_off,y = y0 - i*entry_spacing + text_y_off,xanchor="left",yanchor="top",xref='paper',yref='paper',text="Test Specificity: 0.9",showarrow=False,font=dict(size=17),align='center')
        i+=1
        fig.add_annotation(x = x0 - normal_txt_off,y = y0 - i*entry_spacing + text_y_off,xanchor="left",yanchor="top",xref='paper',yref='paper',text="Percent of CA Under Protocol: {}%".format(round(total_included_pop/region_pop,3)*100),showarrow=False,font=dict(size=17),align='center')
        
        fig.update_layout(height=500,width=1200)
        fig.layout.margin.update({'t':15, 'b':15,'l':0,'r':200})
        fig.write_image(outDir+"/covid_" + place + "_map_thresh{}_iter_{}".format(freq,t) + ".pdf")
        fig.write_image(outDir+"/covid_" + place + "_map_thresh{}_iter_{}".format(freq,t) + ".png")
        #fig.show()
        
        figs.append(fig)
        
    panel_maps[freq] = {'figs':figs, 'binning_endpoints': binning_endpoints, 'freq':freq, 'cost': map_cp[freq]}
        
print('Creating individual maps took {} seconds'.format(time.time() - starttime))

# In[ ]:


def make_panel(figs, binning_endpoints, freq, cost, outDir = 'county_output/Figure6/Multipanels/' + inFileName.split(".csv")[0]):
    if not os.path.isdir(outDir):
        os.mkdir(outDir)
    ROWS = 2
    COLS = 3

    fig = make_subplots(
        rows=ROWS, cols=COLS,
        column_widths=[0.5]*COLS,
        row_heights=[0.5]*ROWS
    )

    k = 0
    annotations = []
    shapes = []
    for row in range(1,ROWS+1):
        for col in range(1,COLS+1):
            if k == len(figs): break
            for j in range(len(figs[k].data)):
                fig.add_trace(figs[k].data[j],row=row,col=col)

            fig['layout']['xaxis' if k == 0 else 'xaxis%i'%(k+1)]['showticklabels']=False
            fig['layout']['yaxis' if k == 0 else 'yaxis%i'%(k+1)]['showticklabels']=False

            # Custom Legend
            x0,y0 = -119.0,41.3
            y0_title_off = 0.8
            normal_txt_off = 0.14
            box_width = 0.5
            box_height = 0.1
            text_x_off = 0.5
            text_y_off = 0.3
            entry_spacing = 0.5
            xref = 'x' if k == 0 else 'x%i'%(k+1)
            yref = 'y' if k == 0 else'y%i'%(k+1)

            annotations.append(dict(x=x0-normal_txt_off,y=y0+y0_title_off,xanchor="left",yanchor="top",xref=xref,yref=yref,text="Time Step: {}".format(panel_times[k]),showarrow=False,font=dict(size=23),align='center'))
        
            k+=1

    # Information about the testing protocol
    xref = 'paper'
    yref = 'paper'   
    x0,y0 = 0.73,0.365
    y0_title_off = 0.035
    normal_txt_off = 0.004
    box_width = 0.015
    box_height = 0.006
    text_x_off = 0.01
    text_y_off = 0.0135
    entry_spacing = 0.02

    annotations.append(dict(x = x0 - normal_txt_off,y = y0 + y0_title_off,xanchor="left",yanchor="top",xref='paper',yref='paper',text="Detected per 100k | Days Between Tests",showarrow=False,font=dict(size=17),align='center'))
    #annotations.append(dict(x = x0 - normal_txt_off,y = y0 + y0_title_off,xanchor="left",yanchor="top",xref='paper',yref='paper',text="Active Detected Infections | Days Between Tests",showarrow=False,font=dict(size=17),align='center'))
        
    for i in range(len(legends[freq])):
        color = colorscale[i]
        shapes.append(dict(type="rect",xref="paper",yref="paper",x0 = x0, y0 = y0 - i*entry_spacing, x1= x0 + box_width, y1= y0 - i*entry_spacing + box_height,line=dict(color=color,width=3),fillcolor=color))
        if i == 0: text = days[freq][i]
        else: text = "{} - {} | {}".format(legends[freq][i-1],legends[freq][i], days[freq][i])
        #else: text = "{}% - {}% | {}".format(legends[freq][i-1],legends[freq][i], days[freq][i])
        annotations.append(dict(x= x0 + box_width + text_x_off, y= y0 - i*entry_spacing + text_y_off, xanchor="left",yanchor="top",xref='paper',yref='paper',text=text,showarrow=False,font=dict(size=17),align='center'))
    if "pcr" in inFileName and "uniform" in inFileName:
        i+=1
        annotations.append(dict(x = x0 - normal_txt_off,y = y0 - i*entry_spacing + text_y_off,xanchor="left",yanchor="top",xref='paper',yref='paper',text="PCR Based Uniform Testing",showarrow=False,font=dict(size=17),align='center'))
    else:
        i+=1
        annotations.append(dict(x = x0 - normal_txt_off,y = y0 - i*entry_spacing + text_y_off,xanchor="left",yanchor="top",xref='paper',yref='paper',text="Test Sensitivity: 0.8",showarrow=False,font=dict(size=17),align='center'))
        i+=1
        annotations.append(dict(x = x0 - normal_txt_off,y = y0 - i*entry_spacing + text_y_off,xanchor="left",yanchor="top",xref='paper',yref='paper',text="Test Specificity: 0.9",showarrow=False,font=dict(size=17),align='center'))
    i+=1
    annotations.append(dict(x = x0 - normal_txt_off,y = y0 - i*entry_spacing + text_y_off,xanchor="left",yanchor="top",xref='paper',yref='paper',text="Percent of CA Under Protocol: {}%".format(round(total_included_pop/region_pop,3)*100),showarrow=False,font=dict(size=17),align='center'))
    i += 1
    annotations.append(dict(x= x0-normal_txt_off, y= y0 - i*entry_spacing + text_y_off, xanchor="left",yanchor="top",xref=xref,yref=yref,text="Cost per Person per Day: ${}".format(round(cost,2)),showarrow=False,font=dict(size=17),align='center'))


    fig.layout.margin.update({'t':50, 'b':100})
    fig.update_layout(
        height=1500,
        width=1500,
        showlegend=False,
        plot_bgcolor='rgb(255,255,255)',
        paper_bgcolor='rgb(255,255,255)',
        annotations=annotations,
        shapes=shapes,
    )

    #fig.show()
    figname = outDir + "/" + place + "_map_progression" + "_thresh{}".format(freq) + ".pdf"
    fig.write_image(figname)
    figname1 = outDir + "/" + place + "_map_progression" + "_thresh{}".format(freq) + ".png"
    fig.write_image(figname1)
    
    return figname


starttime = time.time()
files = []
for key, entry in panel_maps.items():
    files.append(make_panel(**entry))
print('Creating multipanels took {} seconds'.format(time.time() - starttime))

combinepdfs = '/System/Library/Automator/Combine\ PDF\ Pages.action/Contents/Resources/join.py'
outDir = '../plots/Figure6/Multipanels'
cmd = combinepdfs + ' --output ' +  outDir + "/" + place + "_map_progression" + ".pdf"
for file in files: 
    cmd = cmd + " " + file

os.system('$combinepdfs $cmd')