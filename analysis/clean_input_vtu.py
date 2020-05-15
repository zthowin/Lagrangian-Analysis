#-----------------------------------------------------------
# Helper script to remove trailing 0's from .vtu files
#
# Author:       Zachariah Irwin
# Institution:  University of Colroado Boulder
# Last Edits:   May 2020
#-----------------------------------------------------------
import os

numFiles = 181
structNames = ['M1', 'M2', 'M3', 'M4', 'M7', 'M8', 'M9']
leakNames   = ['0PL', '10PL', '20PL', '40PL']
baseDir     = '/dir/'

for clotGeo in ['HS', 'LS']:
  for struct in structNames:
    for leak in leakNames:
      leakNo = leak.split('P')[0]
      for fileInd in range(0, numFiles+1):

        if len(str(fileInd)) == 1:
          add = '00'
        elif len(str(fileInd)) == 2:
          add = '0'
        else:
          add = ''
        velocityDir = '/Velocity-%s-%s-%spL/' %(clotGeo, struct, leakNo)
        filename    = baseDir + velocityDir + 'Velocity-%s-%s-%spL000%s.vtu' %(clotGeo, struct, leakNo, add + str(fileInd))

        # Remove first cardiac cycle, extra file
        if (fileInd < 90) or (fileInd >= numFiles - 1):
          os.remove(filename)
          print("Removed %s" %filename)

        else:
          newfilename = baseDir + velocityDir + 'Velocity-%s-%s-%spL000%s.vtu' %(clotGeo, struct, leakNo, str(fileInd))
          os.rename(filename, newfilename)
          print("Renamed %s" %filename)
