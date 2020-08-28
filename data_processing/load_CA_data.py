import csv, glob, os
import matplotlib.pyplot as plt
import datetime
from datetime import date, timedelta, datetime
import numpy as np
import matplotlib.dates as mdates


def load_cases_by_date(in_dir = '../data/CA/raw', county = 'Los Angeles'):
    casesCSV = 'CasesDeaths.csv'
    reader = csv.reader(open(os.path.join(in_dir,casesCSV),encoding='utf-8'))

    data = {'date':[],'cases_pos_tot':[],'cases_pos_new':[],'cases_prob_tot':[],'cases_prob_new':[]}
    pos_tot = 0
    for row in reader:
        if '-' not in row[5]: continue
        elif county not in row[0]: continue
        else:
            # row = row.split(',')
            sep = '-'
            date = row[5].split(sep)
            sep = '/'
            row[5] = [date[1],date[2],date[0]]
            row[5] = sep.join(row[5])
            data['date'].append(row[5])
            # For missing data put a -1. This value will be interpolated from other values at analysis level.
            try: 
                data['cases_pos_tot'].append(int(float(row[1])))
            except:
                data['cases_pos_tot'].append(-1)
    return data

def load_testing_data(in_dir = '../data/CA/raw', county = 'Los Angeles'):
    # url: http://dashboard.publichealth.lacounty.gov/covid19_surveillance_dashboard/
    casesCSV = 'LA_County_Covid19_tests_date_table.csv'
    reader = csv.reader(open(os.path.join(in_dir,casesCSV)))

    data = {'date':[],'tot_tests':[]}
    pos_tot = 0
    for row in reader:
        if '-' not in row[1]: continue
        else:
            # row = row.split(',')
            sep = '-'
            date = row[1].split(sep)
            sep = '/'
            row[1] = [date[1],date[2],date[0]]
            row[1] = sep.join(row[1])
            data['date'].append(row[1])
            # For missing data put a -1. This value will be interpolated from other values at analysis level.
            try: 
                data['tot_tests'].append(int(float(row[2])))
            except:
                data['tot_tests'].append(-1)
    return data

def load_hospitilization_from_hospitals(in_dir = '../data/CA/raw', county = 'Los Angeles'):
    
    newHospitCSV = 'Hospitalizations.csv'
    reader = csv.reader(open(os.path.join(in_dir,newHospitCSV)))

    data_dates = {}

    data = {'date':[],'hosp_tot_today':[],'hosp_new':[],'hosp_conf_today':[],'hosp_susp_today':[]}
    for row in reader:
        if '-' not in row[1]: continue
        elif county not in row[0]: continue
        else:
            # row = row.split(',')
            sep = '-'
            date = row[1].split(sep)
            sep = '/'
            row[1] = [date[1],date[2],date[0]]
            row[1] = sep.join(row[1])
            data['date'].append(row[1])
            if row[2] != '':
                data['hosp_conf_today'].append(int(float(row[2])))
            else:
                data['hosp_conf_today'].append(0)
            if row[3] != '':
                data['hosp_susp_today'].append(int(float(row[3])))
            else:
                data['hosp_susp_today'].append(0)
            data['hosp_tot_today'].append(data['hosp_susp_today'][-1]+data['hosp_conf_today'][-1])

    return data

def load_date_of_death(in_dir = '../data/CA/raw', county = 'Los Angeles'):
    deathCSV = 'CasesDeaths.csv'
    reader = csv.reader(open(os.path.join(in_dir,deathCSV)))
    data = {'date':[],'death_conf':[],'death_conf_tot':[],'death_prob':[],'death_prob_tot':[]}
    death_tot = 0
    for row in reader:
        if '-' not in row[5]: continue
        elif county not in row[0]: continue
        else:
            # row = row.split(',')
            sep = '-'
            date = row[5].split(sep)
            row[5] = [date[1],date[2],date[0]]
            sep = '/'
            row[5] = sep.join(row[5])
            data['date'].append(row[5])
            # For missing data put a -1. This value will be interpolated from other values at analysis level.
            try: 
                data['death_conf'].append(int(float(row[2])))
            except:
                data['death_conf'].append(-1)
            death_tot += int(float(row[4]))
            data['death_conf_tot'].append(death_tot)
    return data

def load_counties(in_dir = '../data/CA/raw'):
    counties = ['Alameda', 'Butte', 'Calaveras', 'Colusa', 'Contra Costa', 'Del Norte', 'El Dorado', 'Fresno', 'Glenn', 'Humboldt', 'Imperial', 'Inyo', 'Kern', 'Kings', 'Lake', 'Lassen', 'Los Angeles', 'Madera', 'Marin', 'Mariposa', 'Mendocino', 'Merced',  'Mono', 'Monterey', 'Napa', 'Nevada', 'Orange', 'Placer', 'Plumas', 'Riverside', 'Sacramento', 'San Benito', 'San Bernardino', 'San Diego', 'San Francisco', 'San Joaquin', 'San Luis Obispo', 'San Mateo', 'Santa Barbara', 'Santa Clara', 'Santa Cruz', 'Shasta', 'Siskiyou', 'Solano', 'Sonoma', 'Stanislaus', 'Sutter', 'Tehama', 'Trinity', 'Tulare','Ventura', 'Yolo', 'Yuba']
    county_data = {}
    pop_csv = 'CountyPops.csv'
    reader = csv.reader(open(os.path.join(in_dir,pop_csv),encoding='utf8'))

    for county in counties: 
        for row in reader:
            if county not in row[6]: continue
            else:
                if int(str(row[18])) > 100000:
                    county_data[county] = {'cases':load_cases_by_date(in_dir,county),'hosp':load_hospitilization_from_hospitals(in_dir,county),'deaths':load_date_of_death(in_dir,county)}
                    county_data[county]['pop'] = int(str(row[18]))
                break
    return county_data

