# Script: trackingcodedaokragh.py
# Author: Connor Kragh
# Organization: UMBC Observatory
# Description: A user inputs a series of images taken with the UMBC Observatory equipment,
# and the relative RA and Dec tracking rates are output into some file format

# NOTES: The camera used to capture the images being analyzed MUST be NORTH ALIGNED 
# or the RA and Dec rates will be inaccurate

# -*- coding: utf-8 -*-
from astropy.io import fits
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from photutils.detection import DAOStarFinder, IRAFStarFinder

# Sets up scopes array
scopes = [["Main Scope", 6500], ["Finderscope", 1140]]

# Asks user how many frames were taken
numOfCaps = int(input("How many frames were taken?\n" ))

# timeInt is measured in Seconds and is the waiting time in between captures
timeInt = int(input("What is the time interval in between captures (s)?\n" ))
print(f"Taken {numOfCaps} captures with {timeInt} seconds between captures")

# preps file array
file = []
# PROVIDE INPUT FILES HERE
for x in range(numOfCaps):
  file.append(input(f"What is file #{x+1}?\n"))
  print("")

# User selects scope used to take picture
print(f"1. {scopes[0][0]}\n2. {scopes[1][0]}")
scopeNum = int(input("Which scope was used to take the image? (#)\n")) - 1


# Uses name of object and date observed as the naming convention
name = input("What object was captured?\n" )
date = input("What was the date when these were captured? (YYYYMMDD)\n")


#===================================================================
# SET THE DESTINATION PATH FOR THE DATA OUTPUT
worksheetName = "<DESIRED PATH HERE>" + "\\" + name + "_" + date + ".csv"
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

# THIS BLOCK OUTPUTS TO GOOGLE SHEETS DOC

# Morphs tracking data into structure that can be outputted
#sh = gc.open('TrackingCodeOutput')

#Loads Data into dataframe
df = pd.DataFrame(trackingRates[1:], columns=['Source', 'Capture', 'X Centroid (pix)','X Movement (")', 'X Rate ("/s)', 'Y Centroid (pix)', 'Y Movement (")','Y Rate ("/s)'])
print(df)

df_col_names = [df.columns.values.tolist()]
df_col_values = df.values.tolist()
#get the dimensions of the dataframe
df_rows = df.shape[0]
df_cols = df.shape[1]

df.to_csv(path_or_buf=worksheetName)
