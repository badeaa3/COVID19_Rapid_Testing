import csv, glob, os
import matplotlib.pyplot as plt
import datetime
from datetime import date, timedelta, datetime
import numpy as np
import matplotlib.dates as mdates


def load_cases_by_date(in_dir = '../data/NYC/raw'):
	casesCSV = 'coronavirus-data-master/case-hosp-death.csv'
	reader = csv.reader(open(os.path.join(in_dir,casesCSV)))

	data = {'date':[],'cases_pos_tot':[],'cases_pos_new':[],'cases_prob_tot':[],'cases_prob_new':[]}
	pos_tot = 0
	for row in reader:
		if '/' not in row[0]: continue
		pos_tot += int(row[1])
		data['date'].append(row[0])
		data['cases_pos_tot'].append(pos_tot)
		data['cases_pos_new'].append(int(row[1]))
		data['cases_prob_tot'].append(0)
		data['cases_prob_new'].append(0)
	return data

def load_hospitilization_from_hospitals(in_dir = '../data/NYC/raw'):
	
	hospitCSV = 'Hospitalizations.csv'
	newHospitCSV = 'coronavirus-data-master/case-hosp-death.csv'
	reader = csv.reader(open(os.path.join(in_dir,hospitCSV),encoding='utf-16'))

	data_dates = {}

	data = {'date':[],'hosp_tot_today':[],'hosp_new':[],'icu':[], 'icu_new':[], 'intub_tot':[], 'intub_new':[]}
	for row in reader:
		row = row[0].split('\t')
		if '/' not in row[3]: continue
		else:
			sep = '/'
			date = row[3].split(sep)
			if len(date[0]) == 1:
				date[0] = '0' + date[0]
			if len(date[1]) == 1:
				date[1] = '0' + date[1]
			row[3] = sep.join(date)
			if str(row[3]) not in data_dates:
				data_dates[str(row[3])] = {'icu':0,'hosp_tot':0}
			if len(row[0]) > 0:
				data_dates[str(row[3])]['icu'] += int(row[0])
			if len(row[1]) > 0:
				data_dates[str(row[3])]['hosp_tot'] += int(row[1])

	reader = csv.reader(open(os.path.join(in_dir,newHospitCSV)))
	for row in reader:
		if '/' not in row[0]: continue
		data['date'].append(str(row[0]))
		try:
			data['hosp_tot_today'].append(data_dates[str(row[0])]['hosp_tot'])
		except:
			data['hosp_tot_today'].append(-1)
		data['hosp_new'].append(int(row[2]))
		try:
			data['icu'].append(data_dates[str(row[0])]['icu'])
		except:
			data['icu'].append(-1)
		try:
			data['intub_tot'].append(int(row[6]))
			data['intub_new'].append(int(row[7]))
		except:
			data['intub_tot'].append(-1)
			data['intub_new'].append(-1)
	return data

def load_date_of_death(in_dir = '../data/NYC/raw'):
	deathCSV = 'coronavirus-data-master/case-hosp-death.csv'
	excessDeathCSV = 'ExcessDeathsWeekly.csv'
	reader = csv.reader(open(os.path.join(in_dir,deathCSV)))
	data = {'date':[],'death_conf':[],'death_conf_tot':[],'death_prob':[],'death_prob_tot':[]}
	death_tot = 0
	for row in reader:
		if '/' not in row[0]: continue
		data['date'].append(row[0])
		data['death_conf'].append(int(row[3]))
		death_tot += int(row[3])
		data['death_conf_tot'].append(death_tot)

	reader = csv.reader(open(os.path.join(in_dir,excessDeathCSV)))

	sdate = datetime.strptime("2020-02-29","%Y-%m-%d").date()

	print(data['date'])

	dpt = 0
	
	for row in reader:
		if '-' not in row[0]: continue
		if 'New York City' not in row[1]: continue
		date = datetime.strptime(row[0],"%Y-%m-%d").date()
		date2 = '{:%m/%d/%Y}'.format(date)
		if date2 in data['date']:
			if len(data['death_prob']) == data['date'].index(date2):
				dpt += int(row[6])
				for k in range(7):
					if len(data['death_prob']) < len(data['date']):
						data['death_prob'].append(int(row[6])/7)
						data['death_prob_tot'].append(dpt/7)
	while len(data['death_prob']) < len(data['date']):
		data['death_prob'].append(0)
		data['death_prob_tot'].append(0)
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
