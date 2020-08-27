"""
Description: Program to download all available raw data from the Mass.gov website https://www.mass.gov/info-details/archive-of-covid-19-cases-in-massachusetts#related-. 

Author: Anthony Badea
Date: June 7, 2020

"""

import os
import sys
import glob
import csv
import requests, zipfile, io
from collections import OrderedDict

# Download files associated with date and return string to download directory
def download_data(out_dir = "../data/NYC/raw"):
	url = "https://github.com/nychealth/coronavirus-data/archive/master.zip"
	r = requests.get(url, allow_redirects=True)
	z = zipfile.ZipFile(io.BytesIO(r.content))
	# try:
	# z = zipfile.ZipFile(r.content)
	# except:
	# 	print("Date {} doesn't exist".format(date))
	# 	return None
	# out_dir = os.path.join(out_dir,date)
	# if os.path.isdir(out_dir) is not True:
	# 		os.mkdir(out_dir)
	z.extractall(out_dir)
	url = "https://covid19tracker.health.ny.gov/vizql/w/DailyHospitalizationSummary/v/Reopening-DailyHospitalization/vud/sessions/8FB12EBDC72945578208D367A3827F69-5:3/views/12302628778901485932_5034431464735168763?csv=true"
	r2 = requests.get(url, allow_redirects=True)
	f2 = open(out_dir+"/Hospitalizations.csv","wb+")
	f2.write(r2.content)
	# f = io.BytesIO(r.content)
	url = 'https://data.cdc.gov/api/views/xkkf-xrst/rows.csv?accessType=DOWNLOAD&bom=true&format=true%20target='
	r = requests.get(url, allow_redirects=True)
	f = open(out_dir+"/ExcessDeathsWeekly.csv","wb+")
	f.write(r.content)
	return out_dir
				 

# Return dictionary of available dates (key,val) = (date, "")
# The value will be updated by download_data into the local directory of the file
def make_available_data_dict():
	months = ['april','may','june']
	days = range(1,32)
	available_data = OrderedDict()
	for year in ['2020']:
		for month in ['april','may','june']:
			for day in range(1,32):
				available_data[month+'-{}'.format(day)+'-{}'.format(year)] = ""
	return available_data

# Download all days inside of dictionary outputted by make_available_data_dict
# Return a dictionary with (key,val) = (date, local directory to data)
# def download_data_batch(out_dir = "../data/NYC/raw"):
# 	data_dict = make_available_data_dict()
# 	if os.path.isdir(out_dir) is not True: os.mkdir(out_dir)
# 	for date in data_dict.keys():
# 		data_dict[date] = download_data(date,out_dir)
# 		if data_dict[date] is None: data_dict.pop(date)
# 	return data_dict

data_dict = download_data()
#print(data_dict)

# Process all data inside of in_dir and produce a single data file
# Anthony you are here processing the raw data!!!
def process_data(in_dir = "",
				 out_dir = ""):
	# Output data dictionary
	data = {}
	list_of_dirs = glob.glob(in_dir+'/*')
	if len(list_of_dirs) is 0:
		list_of_dirs = [in_dir]
	for lD in list_of_dirs:
		list_of_files = glob.glob(lD+'/*.csv')
		for file in list_of_files:
			csv_in = csv.DictReader(open(file))
			for row in csv_in:
				print(row)

# process_data('../data/NYC/raw/')

	








