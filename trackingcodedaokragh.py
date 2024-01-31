# Script: trackingcodedaokragh.py
# Author: Connor Kragh
# Organization: UMBC Observatory
# Description: A user inputs a series of images taken with the UMBC Observatory equipment,
# and the relative RA and Dec tracking rates are output into some file format

# NOTES: The camera used to capture the images being analyzed MUST be NORTH ALIGNED 
# or the RA and Dec rates will be inaccurate

# -*- coding: utf-8 -*-
import argparse
import os
from astropy.io import fits
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from photutils.detection import DAOStarFinder, IRAFStarFinder

# Sets up scopes array
scopes = [["Main Scope", 6500], ["Finderscope", 1140]]

# Adds and parses arguments from command line
parser = argparse.ArgumentParser(prog='Tracking Analysis', usage='%(prog)s Options')
parser.add_argument('-files', help='Folder containing input FITS files to be processed go here')
parser.add_argument('-interval', help='the length of time between captures (s)')
parser.add_argument('-telescope', help='input "MAIN" or "FINDER" based on scope used for data capture')
parser.add_argument('-object', help='input the name of the object being imaged')
parser.add_argument('-dateOfCapture', help='input date of capture in format YYYYMMDD')
parser.add_argument('-outpath', help='input the path to a folder that can accept the outputted data')
args = parser.parse_args()

# preps array of files
file = []

# Extracts files from computer and adds to file vector or takes input if not supplied
if args.files == None:
# Asks user how many frames were taken
  numOfCaps = int(input("How many frames were taken?\n" ))
  
  # PROVIDE INPUT FILES HERE
  for x in range(numOfCaps):
    file.append(input(f"What is file #{x+1}?\n"))
    print("")

# Adds group of files to the file vector
else:
  for x in os.listdir(args.files):
    file.append(args.files + '\\' + x)

  # Determines amount of captures in series
  numOfCaps = len(file)

# Extracts the time interval between captures or takes input if not supplied

if args.interval == None:
# timeInt is measured in Seconds and is the waiting time in between captures
  timeInt = int(input("What is the time interval in between captures (s)?\n" ))
  print(f"Taken {numOfCaps} captures with {timeInt} seconds between captures")

else:
  timeInt = int(args.interval)

# Extracts the telescope used from parser
# If either no telescope is entered or an invalid entry is provided, it asks for it
if args.telescope == "MAIN":
  scopeNum = 0
elif args.telescope == "FINDER":
  scopeNum = 1
else:
  # User selects scope used to take picture
  print(f"1. {scopes[0][0]}\n2. {scopes[1][0]}")
  scopeNum = int(input("Which scope was used to take the image? (#)\n")) - 1
  while scopeNum != 0 and scopeNum != 1:
    print(f"Invalid Telescope input. \n1. {scopes[0][0]}\n2. {scopes[1][0]}")
    scopeNum = int(input("Which scope was used to take the image? (#)\n")) - 1

# Extracts the object name from parser or takes input if not supplied
    
if args.object == None:
  # Uses name of object in naming convention
  name = input("What object was captured?\n")
else:
  name = args.object

#Extracts the date observed from parser or takes input if not supplied

if args.dateOfCapture == None:
  # Uses ndate observed in naming convention
  date = input("What was the date when these were captured? (YYYYMMDD)\n")
else:
  date = args.dateOfCapture

# Extracts the desired folder to output the data to or takes input if not supplied

if args.outpath == None:
  path = input("Where should the analyzed data be output? (Copy Path of desired folder)\n")

else:
  path = args.outpath

#===================================================================
# SET THE DESTINATION PATH FOR THE DATA OUTPUT
worksheetName = path + "\\" + name + "_" + date + ".csv"
#===================================================================


# Builds a list of captures
hdus = []
for x in file:
  hduList = fits.open(x)
  hdus.append(hduList[0])

# Converts mm per pixel to arcseconds per pixel for later conversion
xArcsecsPerPixel = (hduList[0].header["XPIXSZ"]/scopes[scopeNum][1])*3600
yArcsecsPerPixel = (hduList[0].header["YPIXSZ"]/scopes[scopeNum][1])*3600


captureArray = []
for x in hdus:
  # Plots the captures
  plt.figure()
  plt.imshow(x.data, cmap='gray')

  mean, median, std, max = np.mean(x.data), np.median(x.data), np.std(x.data), np.max(x.data)
  skyBrightness = 300
  # finding statistical values
  starFind = DAOStarFinder(threshold=median, fwhm=20.0, sky=skyBrightness, exclude_border=True, brightest=10, peakmax=70000)
  sources = starFind(x.data)

  # Finds number of sources pictured based on some threshhold
  numOfSources = 0
  for y in sources:
    if y["peak"] > 10*skyBrightness:
      numOfSources = numOfSources + 1

  # Builds array of source data
  captureArray.append(sources[:numOfSources])

  # Marks on plots where sources are
  for y in sources[:numOfSources]:
    plt.scatter(int(y["xcentroid"]), int(y["ycentroid"]), facecolors='none', edgecolors="r")

# Builds data structure of tracking data

# To access element: captureArray[capture number - 1][source number - 1]['column name']
# Example: captureArray[0][0]['xcentroid'] returns the xcentroid of the first source in the first capture
trackingRates = [['Source', 'Capture', 'X Centroid (pix)','X Movement (")', 'X Rate ("/s)', 'Y Centroid (pix)', 'Y Movement (")','Y Rate ("/s)']]
for i in range(numOfSources):
  for j in range(numOfCaps):

    if j == 0:
      trackingRates.append([i+1, j+1, captureArray[j][i]['xcentroid'], 'N/A', 'N/A', captureArray[j][i]['ycentroid'], 'N/A', 'N/A'])
    else:
      # finding x rate
      xMovement = (captureArray[j][i]['xcentroid'] - captureArray[j-1][i]['xcentroid']) * xArcsecsPerPixel
      xRate = xMovement/timeInt
      # finding y rate
      yMovement = (captureArray[j][i]['ycentroid'] - captureArray[j-1][i]['ycentroid']) * yArcsecsPerPixel
      yRate = yMovement/timeInt
      # adding to rates
      trackingRates.append([int(i+1), int(j+1), captureArray[j][i]['xcentroid'], xMovement, xRate, captureArray[j][i]['ycentroid'], yMovement, yRate])

#Loads Data into dataframe
df = pd.DataFrame(trackingRates[1:], columns=['Source', 'Capture', 'X Centroid (pix)','X Movement (")', 'X Rate ("/s)', 'Y Centroid (pix)', 'Y Movement (")','Y Rate ("/s)'])
print(df)

df_col_names = [df.columns.values.tolist()]
df_col_values = df.values.tolist()
#get the dimensions of the dataframe
df_rows = df.shape[0]
df_cols = df.shape[1]

df.to_csv(path_or_buf=worksheetName)
