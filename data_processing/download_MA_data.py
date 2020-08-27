"""
Description: Program to download all available raw data from the Mass.gov website https://www.mass.gov/info-details/archive-of-covid-19-cases-in-massachusetts#related-. 

Author: Anthony Badea
Date: June 7, 2020

"""

import os
import sys
import glob
import csv
import requests, zipfile
from collections import OrderedDict

# Download files associated with date and return string to download directory
def download_data(date = "july-23-2020",
				  out_dir = "../data/MA/raw"):
	url = "https://www.mass.gov/doc/covid-19-raw-data-{}/download".format(date)
	r = requests.get(url, allow_redirects=True)
	try: z = zipfile.ZipFile(r.content)
	except:
		print("Date {} doesn't exist".format(date))
		return None
	out_dir = os.path.join(out_dir,date)
	if os.path.isdir(out_dir) is not True:
		os.mkdir(out_dir)
	z.extractall(out_dir)
	return os.path.join(out_dir,date)
				 

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
def download_data_batch(out_dir = "../data/MA/raw"):
	data_dict = make_available_data_dict()
	if os.path.isdir(out_dir) is not True: os.mkdir(out_dir)
	for date in data_dict.keys():
		data_dict[date] = download_data(date,out_dir)
		if data_dict[date] is None: data_dict.pop(date)
	return data_dict

#data_dict = download_data_batch()
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
	print(list_of_dirs)
	for lD in list_of_dirs:
		list_of_files = glob.glob(lD+'/*.csv')
		for file in list_of_files:
			csv_in = csv.DictReader(open(file))
			for row in csv_in:
				print(row)

download_data()

process_data('../data/MA/raw/july-23-2020')







