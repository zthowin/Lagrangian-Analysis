#-------------------------------------------------------------------------------------------------
# This is a helper script to plot data to go alongside with the Lagrangian toolkit.
#
# Author:       Zachariah Irwin
# Institution:  University of Colorado, Boulder
# Last Edit:    October 2019
#-------------------------------------------------------------------------------------------------
import os
import numpy as np
import csv
import matplotlib.pyplot as plt
from scipy.integrate import simps

#------------
# MAIN SCRIPT
#------------
baseDir       = os.getcwd()

dataDir       = baseDir + '/output/'
if not os.path.exists(dataDir):
    os.mkdir(dataDir)

csvDir        = dataDir + 'CSV_Data/'
if not os.path.exists(csvDir):
    os.mkdir(csvDir)

pngDir        = dataDir + 'PNG_Data/'
if not os.path.exists(pngDir):
    os.mkdir(pngDir)

numFiles      = 91   # Number of files in the vtkDoubleArray dataset (91 for FDHS and FDLS)
numPlotPoints = 1000 # Change resolution of plot

BoundaryIn_Filename    = baseDir + '/inputs-HS/clot-boundary-release.vtk' # This is the filname for the vtkPolyData object
BoundaryOut_Filename   = dataDir + 'HS_Boundary_Points.csv'               # This is the filename for the boundary spline to load
                                                                          # into Paraview for probing
ProbeIn_X_Prefix       = 'HS_Boundary_X'                                  # Prefix for .csv you save with Paraview (x1 direction)
ProbeIn_Y_Prefix       = 'HS_Boundary_Y'                                  # Prefix for .csv you save with Paraview (x2 direction)

print("Reading vtkPolyData object")

Boundary_Points = extractPointData(BoundaryIn_Filename)        # Read in .vtk polydata

print("Finished reading vtkPolyData object\n")
print("Generating spline")

gridX, gridY       = getBoundaryCoordinates(Boundary_Points)      # Extract point coordinates
inds  = np.nonzero(gridY > 0.04)
gridY = gridY[gridY>0.04]
gridX = np.take(gridX, inds)
gridY, inds_unique = np.unique(gridY, return_index=True)
gridX = np.take(gridX, inds_unique)
inds_sorted        = np.argsort(gridX)                               # Get indices for small->large sort
gridX              = np.sort(gridX)                                  # Sort points from smallest to largest for splrep
gridY              = np.take(gridY, inds_sorted)  # Sort points according to gridX
x_Points, y_Points = getPlotPts(gridX, gridY, numPlotPoints)         # Get points along a b-spline interpolation

print("Finished generating spline\n")

np.savetxt(BoundaryOut_Filename, np.transpose([x_Points, y_Points]), header='X, Y', newline=os.linesep, delimiter=',')

print("Spline data saved to %s \n" %BoundaryOut_Filename)
input("Please generate the probe data in Paraview and save to the CSV_Data folder as %s and %s and press any key to continue..."%(ProbeIn_X_Prefix, ProbeIn_Y_Prefix))

bigXArr = np.zeros((numFiles, numPlotPoints), dtype=np.float64)
bigYArr = np.zeros((numFiles, numPlotPoints), dtype=np.float64)

for fileInd in range(numFiles):
  print("Generating plots for file index %i..."%fileInd)
  csvFileName_X = csvDir + ProbeIn_X_Prefix + '.%i.csv' %fileInd
  csvFileName_Y = csvDir + ProbeIn_Y_Prefix + '.%i.csv' %fileInd

  # Column header mapping
  with open(csvFileName_X) as f:
    reader  = csv.reader(f)
    columns = next(reader)
    colmap  = dict(zip(columns, range(len(columns))))

  # Load in data for stretch/strain calculated along x1
  arrayData_X = np.matrix(np.loadtxt(csvFileName_X, delimiter=",", skiprows=1))
  x           = arrayData_X[:,colmap['Points:0']]
  stretch     = arrayData_X[:,colmap['stretch']]
  bigXArr[fileInd, :] = stretch.ravel()

  """
  fig, ax1 = plt.subplots()

  # Plot x1 data along horizontal position
  color = 'tab:red'
  ax1.set_xlabel('Position (x)')
  ax1.set_ylabel('Stretch along x1', color=color)
  ax1.plot(x, stretch, color=color)
  ax1.tick_params(axis='y', labelcolor=color)

  # Plot clot boundary
  ax2 = ax1.twinx()
  color = 'tab:blue'
  ax2.set_ylabel('Clot boundary', color=color)
  ax2.plot(x, y_Points, color=color)
  ax2.tick_params(axis='y', labelcolor=color)
  fig.tight_layout()

  plt.savefig(pngDir + 'stretch_x_t-ind_%i' %fileInd)
  plt.close()
  """
  # Column header mapping
  with open(csvFileName_Y) as f:
    reader  = csv.reader(f)
    columns = next(reader)
    colmap  = dict(zip(columns, range(len(columns))))

  # Load in data for stretch/strain calculated along x2
  arrayData_Y = np.matrix(np.loadtxt(csvFileName_Y, delimiter=",", skiprows=1))
  y           = arrayData_Y[:,colmap['Points:0']]
  stretch     = arrayData_Y[:,colmap['stretch']]
  bigYArr[fileInd, :] = stretch.ravel()
  """
  fig2, ax3 = plt.subplots()

  # Plot x2 data along horizontal position
  color = 'tab:red'
  ax3.set_xlabel('Position (x)')
  ax3.set_ylabel('Stretch along x2', color=color)
  ax3.plot(y, stretch, color=color)
  ax3.tick_params(axis='y', labelcolor=color)

  # Plot clot boundary
  ax4 = ax3.twinx()  # instantiate a second axes that shares the same x-axis
  color = 'tab:blue'
  ax4.set_ylabel('Clot boundary', color=color)  # we already handled the x-label with ax1
  ax4.plot(y, y_Points, color=color)
  ax4.tick_params(axis='y', labelcolor=color)
  fig.tight_layout()

  plt.savefig(pngDir + 'stretch_y_t-ind_%i' %fileInd)
  plt.close()
  """

intXSimps = np.zeros((numPlotPoints,1), dtype=np.float64)
intYSimps = np.zeros((numPlotPoints,1), dtype=np.float64)
dt = 0.01
T  = 0.9

for point in range(numPlotPoints):

  intXSimps[point] = simps(bigXArr[:,point].ravel(), dx=dt)/T
  intYSimps[point] = simps(bigYArr[:,point].ravel(), dx=dt)/T

ratio = 0.3
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.figure()
plt.title(r'Time-averaged stretch along $x_1$')
ax1 = plt.subplot()

# Plot x1 data along horizontal position
color = 'tab:red'
ax1.set_xlabel(r'Position ($x$)')
ax1.set_ylabel(r'Stretch along $x_1$', color=color)
ax1.plot(x, intXSimps, color=color)
ax1.tick_params(axis='y', labelcolor=color)
ax1.set_ylim(0.0, 1400)


# Plot clot boundary
ax2 = ax1.twinx()
color = 'tab:blue'
ax2.set_ylabel('Clot boundary', color=color)
ax2.plot(x, y_Points, color=color)
ax2.tick_params(axis='y', labelcolor=color)
ax2.set_ylim(0.0, 10)
plt.savefig(pngDir + 'FDHS_Stretch_X_TimeAvg.png')
plt.close()

plt.figure()
plt.title(r'Time-averaged stretch along $x_2$')
ax1 = plt.subplot()

# Plot x1 data along horizontal position
color = 'tab:red'
ax1.set_xlabel(r'Position ($x$)')
ax1.set_ylabel(r'Stretch along $x_2$', color=color)
ax1.plot(x, intYSimps, color=color)
ax1.tick_params(axis='y', labelcolor=color)
ax1.set_ylim(0.0, 1400)


# Plot clot boundary
ax2 = ax1.twinx()
color = 'tab:blue'
ax2.set_ylabel('Clot boundary', color=color)
ax2.plot(x, y_Points, color=color)
ax2.tick_params(axis='y', labelcolor=color)
ax2.set_ylim(0.0, 10)
plt.savefig(pngDir + 'FDHS_Stretch_Y_TimeAvg.png')

np.savetxt('FDHS_X_StretchAvg.dat', np.transpose([x,intXSimps])[0,:,:], fmt='%.18g', delimiter=',')
np.savetxt('FDHS_Y_StretchAvg.dat', np.transpose([x,intYSimps])[0,:,:], fmt='%.18g', delimiter=',')

print("Done!")
