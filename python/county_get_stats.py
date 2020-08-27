import pandas as pd
import sys, os, csv, ast
sys.path.insert(1, '../data_processing')
from load_CA_data import *
from county_utilities import *

# Load the counties 
cd = load_counties( '../data/CA/raw/')
pops = {}

sdate = '04/01/2020' 
edate = '07/22/2020'

region_pop = np.sum([cd[county]['pop'] for county in cd.keys()])
total_included_pop = 0

for county in list(cd.keys()):
    #if county not in ["Alameda"] : continue
    DATA = {'cases': cd[county]['cases'], 'hosp':cd[county]['hosp'], 'deaths':cd[county]['deaths']}

    pop = cd[county]['pop'] 
    dt = np.array(DATA['deaths']['death_conf_tot'])[DATA['deaths']['date'].index(sdate):DATA['deaths']['date'].index(edate)] / pop

    # Selections on data sets
    if pop/region_pop < 0.015: continue 
    if not any(c>0.0 for c in dt): continue
    total_included_pop += pop

    # Save the data sets passing selection
    pops[county] = pop

# Total county population
POP = np.sum([pops[county] for county in pops.keys()])

# Check how much of population is used
print("Percentage of California population used: {}".format(total_included_pop/region_pop))

# Load the dictionaries statistics of interest

map_series_pcr_uniform = "map_series_pcr_uniform.csv" 
df = pd.read_csv("county_output/" + map_series_pcr_uniform, index_col=0)
d = df.to_dict("split")
map_series_pcr_uniform  = {}
i = 0
for key in d['index']:
    if key not in map_series_pcr_uniform.keys():
        map_series_pcr_uniform[key] = {}
    for county, data in zip(d['columns'],d['data'][i]):
        map_series_pcr_uniform[key][county] = ast.literal_eval(data) 
    i += 1

map_series_rapid_county = "map_series_rapid_county.csv"
df = pd.read_csv("county_output/" + map_series_rapid_county, index_col=0)
d = df.to_dict("split")
map_series_rapid_county  = {}
i = 0
for key in d['index']:
    if key not in map_series_rapid_county.keys():
        map_series_rapid_county[key] = {}
    for county, data in zip(d['columns'],d['data'][i]):
        map_series_rapid_county[key][county] = ast.literal_eval(data) 
    i += 1


# Compute the statistics of interest 
total_pcr_uniform = {'start':0,'end':0,'means':{}}
D_pcr_uniform = {'start':0,'end':0,'means':{}}
I_pcr_uniform = {'start':0,'end':0,'means':{}}
ID_pcr_uniform = {'start':0,'end':0,'means':{}}
for key, val in map_series_pcr_uniform.items():
	if key != 1.0/7.0: continue
	for county, series in val.items():
		total_pcr_uniform['start'] += series['total'][0]  * pops[county] / POP
		total_pcr_uniform['end']   += series['total'][-1] * pops[county] / POP
		total_pcr_uniform['means'][county] = np.mean(np.array(series['total']))

		D_pcr_uniform['start'] += series['D'][0]  * pops[county] / POP
		D_pcr_uniform['end']   += series['D'][-1] * pops[county] / POP
		D_pcr_uniform['means'][county] = np.mean(np.array(series['D']))

		I_pcr_uniform['start'] += series['I'][0]  * pops[county] / POP
		I_pcr_uniform['end']   += series['I'][-1] * pops[county] / POP
		I_pcr_uniform['means'][county] = np.mean(np.array(series['I']))

		ID_pcr_uniform['start'] += ( series['D'][0] + series['I'][0]) * pops[county] / POP
		ID_pcr_uniform['end']   += ( series['D'][-1] + series['I'][-1] ) * pops[county] / POP
		ID_pcr_uniform['means'][county] = np.mean(np.array(( series['D'] + series['I'] )))


print(total_pcr_uniform)

total_rapid_county = {'start':0,'end':0,'means':{}}
D_rapid_county = {'start':0,'end':0,'means':{}}
I_rapid_county = {'start':0,'end':0,'means':{}}
ID_rapid_county = {'start':0,'end':0,'means':{}}
# Compute the statistics of interest 
for key, val in map_series_rapid_county.items():
	if key != 0.0005: continue
	for county, series in val.items():
		total_rapid_county['start'] += series['total'][0]  * pops[county] / POP
		total_rapid_county['end']   += series['total'][-1] * pops[county] / POP
		total_rapid_county['means'][county] = np.mean(np.array(series['total']))

		D_rapid_county['start'] += series['D'][0]  * pops[county] / POP
		D_rapid_county['end']   += series['D'][-1] * pops[county] / POP
		D_rapid_county['means'][county] = np.mean(np.array(series['D']))

		I_rapid_county['start'] += series['I'][0]  * pops[county] / POP
		I_rapid_county['end']   += series['I'][-1] * pops[county] / POP
		I_rapid_county['means'][county] = np.mean(np.array(series['I']))

		ID_rapid_county['start'] += ( series['D'][0] + series['I'][0] ) * pops[county] / POP
		ID_rapid_county['end']   += ( series['D'][-1] + series['I'][-1] ) * pops[county] / POP
		ID_rapid_county['means'][county] = np.mean(np.array(( series['D'] + series['I'] )))
print(total_rapid_county)

print("PCR Uniform | Rapid County")

print("Total")
for key in total_rapid_county.keys():
	if key == 'means': continue 
	print("{}: {} | {}".format(key, total_pcr_uniform[key],total_rapid_county[key]))

print("-------")

print("D")
for key in D_rapid_county.keys():
	if key == 'means': continue 
	print("{}: {} | {}".format(key, D_pcr_uniform[key],D_rapid_county[key]))

print("-------")

print("I")
for key in I_rapid_county.keys():
	if key == 'means': continue 
	print("{}: {} | {}".format(key, I_pcr_uniform[key],I_rapid_county[key]))

print("-------")

print("I+D")
for key in I_rapid_county.keys():
	if key == 'means': continue 
	print("{}: {} | {}".format(key, ID_pcr_uniform[key],ID_rapid_county[key]))

print("-------")


