import csv, glob, os
import matplotlib.pyplot as plt
import datetime
from datetime import date, timedelta, datetime
import numpy as np
import matplotlib.dates as mdates

def load_cases_by_date(in_dir = '../data/MA/raw',
				 	   date = 'july-23-2020'):
	casesCSV = 'CasesByDate.csv'
	reader = csv.reader(open(os.path.join(os.path.join(in_dir,date),casesCSV)))

	data = {'date':[],'cases_pos_tot':[],'cases_pos_new':[],'cases_prob_tot':[],'cases_prob_new':[]}
	for row in reader:
		if '/' not in row[0]: continue
		data['date'].append(row[0])
		data['cases_pos_tot'].append(int(row[1]))
		data['cases_pos_new'].append(int(row[2]))
		data['cases_prob_tot'].append(int(row[3]))
		data['cases_prob_new'].append(int(row[4]))
	return data

def load_hospitilization_from_hospitals(in_dir = '../data/MA/raw',
				 	   date = 'july-23-2020'):
	hospitCSV = 'Hospitalization from Hospitals.csv'
	reader = csv.reader(open(os.path.join(os.path.join(in_dir,date),hospitCSV)))

	data = {'date':[],'hosp_tot_today':[],'hosp_new':[],'hosp_5day_avg':[],'icu':[], 'icu_new':[], 'intub_tot':[], 'intub_new':[]}
	for row in reader:
		if '/' not in row[0]: continue
		data['date'].append(row[0])
		data['hosp_tot_today'].append(int(row[1]))
		data['hosp_new'].append(int(row[2]))
		data['hosp_5day_avg'].append(float(row[3]))
		data['icu'].append(int(row[4]))
		data['icu_new'].append(int(row[5]))
		try:
			data['intub_tot'].append(int(row[6]))
			data['intub_new'].append(int(row[7]))
		except:
			data['intub_tot'].append(row[6])
			data['intub_new'].append(row[7])
	return data

def load_date_of_death(in_dir = '../data/MA/raw',
				 	   date = 'july-23-2020'):
	deathCSV = 'DateofDeath.csv'
	reader = csv.reader(open(os.path.join(os.path.join(in_dir,date),deathCSV)))
	data = {'date':[],'death_conf':[],'death_conf_tot':[],'death_prob':[],'death_prob_tot':[]}
	for row in reader:
		if '/' not in row[0]: continue
		data['date'].append(row[0])
		data['death_conf'].append(int(row[1]))
		data['death_conf_tot'].append(int(row[2]))
		data['death_prob'].append(int(row[3]))
		data['death_prob_tot'].append(int(row[4]))
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
