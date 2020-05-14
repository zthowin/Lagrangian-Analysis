#-------------------------------------------------------------------------------------------------
# This is a helper script to plot field data along a given line.
#
# Author:       Zachariah Irwin
# Institution:  University of Colorado, Boulder
# Last Edit:    October 2019
#-------------------------------------------------------------------------------------------------
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import simps

#------------
# MAIN SCRIPT
#------------
makePlots = False

no_T = True
w_T  = False

fdls = True
fdhs = False

fvcDataDir    = '/home/zach/Documents/Mukherjee/FlowVC/examples/clot_fdls/'

if fdls:
  horizontal_pos  = np.linspace(5., 27., 12, endpoint=True)
if fdhs:
  horizontal_pos  = np.linspace(8., 28., 6, endpoint=True)

if no_T:
  fvcCSVDir  = fvcDataDir + 'Line_Plots_noT/'
if w_T:
  fvcCSVDir  = fvcDataDir + 'Line_Plots_T/'

baseDir       = os.getcwd()

dataDir       = baseDir + '/Data_FD/FDLS/'
if not os.path.exists(dataDir):
    os.mkdir(dataDir)

csvDir        = dataDir + 'Line_Plots/'
if not os.path.exists(csvDir):
    os.mkdir(csvDir)

pngDir        = dataDir + 'PNG_Data/'
if not os.path.exists(pngDir):
    os.mkdir(pngDir)

numPositions  = 12
numFiles      = 90   # Number of files in the vtkDoubleArray dataset (91 for FDHS and FDLS)
arcLength     = 6.0
dataTiming    = 0.01

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

if makePlots:
  errors_noT = np.zeros((numPositions, numFiles), dtype=np.float64)
  errors_T   = np.zeros((numPositions, numFiles), dtype=np.float64)
  for posInd in range(numPositions):
    print('Making plots for position #%i...' %posInd)
    for fileInd in range(numFiles):
      time = float(fileInd)*dataTiming

      fvcName = fvcCSVDir + 'Pos%i.%i.csv' %(posInd, fileInd)
      dbjName = csvDir + 'Pos%i.%i.csv' %(posInd, fileInd)

      fvcDF = pd.read_csv(fvcName)
      dbjDF = pd.read_csv(dbjName)

      numDataPoints = fvcDF.shape[0] - 1

      s_Points = fvcDF['Points:1']
      s_Points = s_Points.to_numpy()

      fvcFTLE  = fvcDF['scalar_field']
      fvcFTLE  = fvcFTLE.to_numpy()

      if w_T:
        dbjFTLE_T    = dbjDF['ftle - T']
        dbjFTLE_T    = dbjFTLE_T.to_numpy()

        error = np.zeros(numDataPoints, dtype=np.float64)
        for point in range(numDataPoints):
          localError   = np.linalg.norm(fvcFTLE[point] - dbjFTLE_T[point])#/np.linalg.norm(fvcFTLE[point])
          if np.isnan(localError) and (np.linalg.norm(fvcFTLE[point] - dbjFTLE_T[point]) == 0. and np.linalg.norm(fvcFTLE[point]) == 0.):
            local_error = 0.
          error[point] = localError

        # Calculate total error between FTLE data sets along the line
        error_Tot                  = (1./arcLength) * simps(error, dx=arcLength/float(numDataPoints))
        errors_T[posInd, fileInd]  = error_Tot

        # Plot FTLE vs. vertical position for discrete horizontal location in the domain
        plt.figure()
        plt.plot(dbjFTLE_T, s_Points)
        plt.ylim(0., 6.)
        plt.xlim(0., 80.)
        plt.ylabel(r'Arc Length')
        plt.xlabel(r'FTLE')
        plt.title(r'FTLE vs. Vertical Position - With Time Scaling at t = %.2f' %time)
        plt.savefig(pngDir + 'FTLE_wT-Pos%i-TimeInd:%i.png' %(posInd, fileInd))
        #plt.show()
        plt.close()

      if no_T:
        dbjFTLE_noT  = dbjDF['ftle - No T']
        dbjFTLE_noT  = dbjFTLE_noT.to_numpy()

        error = np.zeros(numDataPoints, dtype=np.float64)
        for point in range(numDataPoints):
          localError   = np.linalg.norm(fvcFTLE[point] - dbjFTLE_noT[point])#/np.linalg.norm(fvcFTLE[point])
          #if np.isnan(localError) and (np.linalg.norm(fvcFTLE[point] - dbjFTLE_noT[point]) == 0. and np.linalg.norm(fvcFTLE[point]) == 0.):
            #local_error = 0.
          error[point] = localError

        # Calculate total error between FTLE data sets along the line
        error_Tot                    = (1./arcLength) * simps(error, dx=arcLength/float(numDataPoints))
        errors_noT[posInd, fileInd]  = error_Tot

        # Plot FTLE vs. vertical position for discrete horizontal location in the domain
        plt.figure()
        plt.plot(dbjFTLE_noT, s_Points)
        plt.ylim(0., 6.)
        plt.xlim(0., 10.)
        plt.ylabel(r'Arc Length')
        plt.xlabel(r'FTLE')
        plt.title(r'FTLE vs. Vertical Position - Without Time Scaling at t = %.2f' %time)
        plt.savefig(pngDir + 'FTLE_noT_old-Pos%i-TimeInd:%i.png' %(posInd, fileInd))
        #plt.show()
        plt.close()

  # Save error data
  if no_T:
    np.savetxt(dataDir + 'Error_NoT_old.dat', errors_noT.T, newline=os.linesep, delimiter=',')
  if w_T:
    np.savetxt(dataDir + 'Error_wT.dat', errors_T.T, delimiter=',')

# Load error data
if not makePlots:
  if no_T:
    errors_noT = np.loadtxt(dataDir + 'Error_NoT.dat', delimiter=',')
    errors_noT = errors_noT.T
  if w_T:
    errors_T   = np.loadtxt(dataDir + 'Error_wT.dat', delimiter=',')
    errors_T   = errors_T.T

if fdls:
  # Load boundary data
  boundaryPts = np.loadtxt('/home/zach/Documents/Mukherjee/Lagrangian-Toolkit/Code/Data_FD/LS_Boundary_Points.csv', delimiter=',')
  # Make boxplots for errors during the whole integration window
  if w_T:
    fig, ax = plt.subplots()
    line1 = ax.boxplot(errors_T.T, positions=horizontal_pos, meanline=True, vert=True, showfliers=False)
    ax.set_xlabel(r'Horizontal positions of line plots')
    ax.set_ylabel(r'Error')
    ax.set_title(r'Comparison of FVC with LT - Error in FTLE (with time scaling)')
    ax.set_ylim(-1.25)

    ax2 = ax.twinx()
    line2 = ax2.plot(boundaryPts[:,0], boundaryPts[:,1] - 12.)
    ax2.set_ylim(-12., -2.)
    ax2.yaxis.set_visible(False)
    plt.savefig(pngDir + 'FDLS_wT_Errors.png')
    plt.show()
    plt.close()

  if no_T:
    fig, ax = plt.subplots()
    line1 = ax.boxplot(errors_noT.T, positions=horizontal_pos, meanline=True, vert=True, showfliers=True)
    ax.set_xlabel(r'Horizontal positions of line plots')
    ax.set_ylabel(r'Error')
    ax.set_title(r'Comparison of FVC without LT - Error in FTLE (without time scaling)')
    ax.set_ylim(-0.25, 0.5)

    ax2 = ax.twinx()
    line2 = ax2.plot(boundaryPts[:,0], boundaryPts[:,1] - 12.)
    ax2.set_ylim(-12., -2.)
    ax2.yaxis.set_visible(False)
    plt.savefig(pngDir + 'FDLS_NoT_Errors_zoomed.png')
    plt.show()
    plt.close()


print("Done!")
