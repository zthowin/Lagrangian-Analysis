#------------------------------------------------------------------------
# Helper script to generate input files for batch runs (.dat and .sbatch)
#
# Author:       Zachariah Irwin
# Institution:  University of Colorado Boulder
# Last Edit:    May 2020
#------------------------------------------------------------------------
import os
import shutil
import linecache

sbatch = True
dat    = False
tracer = False
ftle   = True

structNames = ['M1', 'M2', 'M3', 'M4', 'M7', 'M8']
leakNames   = ['0PL', '10PL', '20PL', '40PL']

if tracer:
  locNames  = ['Clot', 'Wall', 'Drug']
elif ftle:
  locNames  = ['']

if dat:
  if tracer:
    linesOI   = [6, 8, 18]
  elif ftle:
    linesOI   = [6, 8, 16]
if sbatch:
  linesOI   = [5, 20]

for clotGeo in ['HS', 'LS']:

  for loc in locNames:
    if dat and tracer:
      templateFileName = os.getcwd() + '/input-%s-%s-Temp.dat' %(clotGeo, loc)

    elif dat and ftle:
      templateFileName = os.getcwd() + '/input-%s-Temp.dat' %(clotGeo)

    else:
      templateFileName = os.getcwd() + '/run_LT.sbatch'

    for struct in structNames:
      for leak in leakNames:

        leakNo = leak.split('P')[0]

        if dat and tracer:
          fileName = os.getcwd() + '/%s-Tracer-Input-Files/input-%s-%s-%s-%s.dat' %(clotGeo, clotGeo, loc, struct, leak)
        elif dat and ftle:
          fileName = os.getcwd() + '/%s-Input-Files/input-%s-%s-%s.dat' %(clotGeo, clotGeo, struct, leak)
        elif sbatch and tracer:
          fileName = os.getcwd() + '/Batch-%s-Tracer/Batch-%s-%s-%s-%s.sbatch' %(clotGeo, clotGeo, loc, struct, leak)
        elif sbatch and ftle:
          fileName = os.getcwd() + '/Batch-%s/Batch-%s-%s-%s.sbatch' %(clotGeo, clotGeo, struct, leak)

        shutil.copy(templateFileName, fileName)

        with open(templateFileName, 'r+') as tempFile:
          with open(fileName, 'w+') as file:

            lineContents   = tempFile.readlines()

            if dat:
              dataDir        = lineContents[linesOI[0]]
              newDataDir     = dataDir.split('/')
              newDataDir[3]  = struct
              newDataDir[4]  = leak
              newDataDir     = '/'.join(newDataDir)

              lineContents[linesOI[0]] = newDataDir

              if struct == 'M1':
                  add = '_0_'
              else:
                  add = '-'
              dataTag        = lineContents[linesOI[1]]
              newDataTag     = dataTag.split('-')
              newDataTag[2]  = struct
              newDataTag[3]  = leakNo + 'pL' + add + '\n'
              newDataTag     = '-'.join(newDataTag)

              lineContents[linesOI[1]] = newDataTag

              if struct == 'M1':
                  dataname = 'Flow data field name = f_22-0' + '\n'
                  lineContents[10] = dataname

              outDir         = lineContents[linesOI[2]]
              newOutDir      = outDir.split('/')
              newOutDir[3]   = struct
              newOutDir[4]   = leak
              newOutDir      = '/'.join(newOutDir)
              newOutDir      = newOutDir.split('-')
              if tracer:
                newOutDir[3]   = struct
                newOutDir[4]   = leakNo + 'pL' + '\n'
              elif ftle:
                newOutDir[2]   = struct
                newOutDir[3]   = leakNo + 'pL'
              newOutDir      = '-'.join(newOutDir)

              lineContents[linesOI[2]] = newOutDir

            if sbatch:
              outDir         = lineContents[linesOI[0]]
              newOutDir      = outDir.split('-')

              if tracer:
                newOutDir[2]   = 'output=/projects/zair9172/LT/sbatch_out/%s-%s-Tracer-%s-%s.out' %(clotGeo, loc, struct, leak)
              elif ftle:
                newOutDir[2]   = 'output=/projects/zair9172/LT/sbatch_out/%s-%s-%s.out' %(clotGeo, struct, leak)

              newOutDir      = '-'.join(newOutDir) + '\n'

              lineContents[linesOI[0]] = newOutDir

              inputFile      = lineContents[linesOI[1]]
              dirs           = inputFile.split('/')

              if tracer:
                dirs[10]     = '%s-Tracer-Input-Files' %(clotGeo)
              elif ftle:
                dirs[10]     = '%s-Input-Files' %(clotGeo)

              inputFile      = '/'.join(dirs)

              if tracer:
                newInputFile   = inputFile.split('\n')[0] + 'input-%s-%s-%s-%s.dat' %(clotGeo, loc, struct, leak) + '\n'
              elif ftle:
                newInputFile   = inputFile.split('\n')[0] + 'input-%s-%s-%s.dat' %(clotGeo, struct, leak) + '\n'

              lineContents[linesOI[1]] = newInputFile

            file.writelines(lineContents)
            print("Wrote file %s" %fileName)

          file.close()
        tempFile.close()
