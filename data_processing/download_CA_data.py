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
def download_data(out_dir = "../data/CA/raw"):
	# url = 'https://data.ca.gov/dataset/590188d5-8545-4c93-a9a0-e230f0db7290/resource/926fd08f-cc91-4828-af38-bd45de97f8c3/download/statewide_cases.csv'
	# r = requests.get(url, allow_redirects=True)
	# z = zipfile.ZipFile(io.BytesIO(r.content))
	# # try:
	# # z = zipfile.ZipFile(r.content)
	# # except:
	# # 	print("Date {} doesn't exist".format(date))
	# # 	return None
	# # out_dir = os.path.join(out_dir,date)
	# # if os.path.isdir(out_dir) is not True:
	# # 		os.mkdir(out_dir)
	# z.extractall(out_dir)
	url = 'https://data.ca.gov/dataset/590188d5-8545-4c93-a9a0-e230f0db7290/resource/926fd08f-cc91-4828-af38-bd45de97f8c3/download/statewide_cases.csv'
	r2 = requests.get(url, allow_redirects=True)
	f2 = open(out_dir+"/CasesDeaths.csv","wb+")
	f2.write(r2.content)
	# f = io.BytesIO(r.content)
	url = 'https://data.ca.gov/dataset/529ac907-6ba1-4cb7-9aae-8966fc96aeef/resource/42d33765-20fd-44b8-a978-b083b7542225/download/hospitals_by_county.csv'
	r = requests.get(url, allow_redirects=True)
	f = open(out_dir+"/Hospitalizations.csv","wb+")
	f.write(r.content)
	url = 'https://www2.census.gov/programs-surveys/popest/datasets/2010-2019/counties/totals/co-est2019-alldata.csv'
	r = requests.get(url, allow_redirects=True)
	f = open(out_dir+"/CountyPops.csv","wb+")
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

	








