import csv, glob, os
import matplotlib.pyplot as plt
import datetime
from datetime import date, timedelta, datetime
import numpy as np
import matplotlib.dates as mdates


def load_cases_by_date(in_dir = '../data/Brazil/raw'):
	casesCSV = 'brazil_confcases.csv'
	reader = csv.reader(open(os.path.join(in_dir,casesCSV)))

	data = {'date':[],'cases_pos_tot':[],'cases_pos_new':[],'cases_prob_tot':[],'cases_prob_new':[]}
	pos_tot = 0
	for row in reader:
		if 'DATE' in row[0] or '/' not in row[0]: continue
		pos_tot = row[3]
		data['date'].append(row[0]+'20')
		data['cases_pos_tot'].append(int(row[4]))
		data['cases_pos_new'].append(0)
		data['cases_prob_tot'].append(0)
		data['cases_prob_new'].append(0)
	return data

def load_hospitilization_from_hospitals(in_dir = '../data/Brazil/raw'):
	
	newHospitCSV = 'brazil_confcases.csv'
	outcomesCSV = 'outcomes2.csv'
	data_dates = {}

	data = {'date':[],'hosp_tot_today':[],'hosp_new':[],'icu':[], 'icu_new':[], 'intub_tot':[], 'intub_new':[]}
	sev = 0 
	exited = 0
	reader = csv.reader(open(os.path.join(in_dir,newHospitCSV)))
	for row in reader:
		if 'DATE' in row[0] or '/' not in row[0]: continue
		sev += int(row[1])
		data['date'].append(str(row[0])+'20')
		data['hosp_tot_today'].append(sev)
	reader = csv.reader(open(os.path.join(in_dir,outcomesCSV)))
	print(data['hosp_tot_today'],'hospitalized')
	print(data['date'])
	dates = []
	for row in reader:
		if 'DATE' in row[0] or '/' not in row[3]: continue
		date = row[3]
		mdy = date.split('/')
		if mdy[0][0] == '0': mdy[0] = mdy[0][1:]
		if mdy[1][0] == '0': mdy[1] = mdy[1][1:]
		mdy = '/'.join(mdy)
		index = data['date'].index(str(mdy))
		exited += int(row[2])
		exited += int(row[1])
		if mdy in dates:
			data['hosp_tot_today'][index] -= (int(row[2])+int(row[1]))
		else:
			data['hosp_tot_today'][index] -= exited
		dates.append(mdy)
	for d in range(len(data['hosp_tot_today'])):
		if data['date'][d] not in dates and d > 0:
			data['hosp_tot_today'][d] = data['hosp_tot_today'][d-1]
	print(exited)
	print(data['hosp_tot_today'],'htt')
	return data

def load_date_of_death(in_dir = '../data/Brazil/raw'):
	deathCSV = 'brazil_confcases.csv'
	outcomesCSV = 'outcomes2.csv'
	reader = csv.reader(open(os.path.join(os.path.join(in_dir),deathCSV)))
	data = {'date':[],'death_conf':[],'death_conf_tot':[],'death_prob':[],'death_prob_tot':[]}
	for row in reader:
		if 'DATE' in row[0] or '/' not in row[0]: continue
		data['date'].append(row[0]+'20')
		data['death_conf_tot'].append(0)
	reader = csv.reader(open(os.path.join(in_dir,outcomesCSV)))
	exited = 0
	for row in reader:
		if 'DATE' in row[0] or '/' not in row[3]: continue
		date = row[3]
		mdy = date.split('/')
		if mdy[0][0] == '0': mdy[0] = mdy[0][1:]
		if mdy[1][0] == '0': mdy[1] = mdy[1][1:]
		mdy = '/'.join(mdy)
		index = data['date'].index(str(mdy))
		exited += int(row[1])
		if data['death_conf_tot'][index] == 0:
			data['death_conf_tot'][index] = exited
		else:
			data['death_conf_tot'][index] += int(row[1])
	for d in range(len(data['death_conf_tot'])):
		if data['death_conf_tot'][d] == 0 and d > 0:
			data['death_conf_tot'][d] = data['death_conf_tot'][d-1]
	print(data['death_conf_tot'],'deaths')
	return data

'''
cases_by_date = load_cases_by_date()
date_of_death = load_date_of_death()
hosp_from_hosps = load_hospitilization_from_hospitals()

case_min, case_max = min(cases_by_date['date']), max(cases_by_date['date'])
death_min, death_max = min(date_of_death['date']), max(date_of_death['date'])
hosp_min, hosp_max = min(hosp_from_hosps['date']), max(hosp_from_hosps['date'])


ax = plt.gca()
formatter = mdates.DateFormatter("%Y-%m-%d")
ax.xaxis.set_major_formatter(formatter)
locator = mdates.MonthLocator()
locator_days = mdates.DayLocator()
ax.xaxis.set_major_locator(locator)
ax.xaxis.set_minor_locator(locator_days)


case_dates = [datetime.strptime(d,"%m/%d/%Y").date() for d in cases_by_date['date']]
plt.plot(case_dates, cases_by_date['cases_pos_new'], label='Confirmed Daily Cases')
#plt.plot(case_dates, cases_by_date['cases_prob_new'], label = 'Probable Cases')


death_dates = [datetime.strptime(d,"%m/%d/%Y").date() for d in date_of_death['date']]
plt.plot(death_dates, date_of_death['death_conf_tot'], label="Confirmed Total Deaths")
#plt.plot(death_dates, date_of_death['death_prob_tot'], label="Probable Total Deaths")
#plt.plot(death_dates, date_of_death['death_conf'], label="Confirmed Daily Deaths")
#plt.plot(death_dates, date_of_death['death_prob'], label="Probable Daily Deaths")

hosp_dates = [datetime.strptime(d,"%m/%d/%Y").date() for d in hosp_from_hosps['date']]
plt.plot(hosp_dates, hosp_from_hosps['hosp_tot_today'], label="Total Daily Hospitalizations")

plt.legend()
plt.show()
'''
