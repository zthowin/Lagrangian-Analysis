import os
import shutil
import linecache

structNames = ['M1', 'M2', 'M3', 'M4', 'M7', 'M8', 'M9']
leakNames   = ['0PL', '10PL', '20PL', '40PL']
linesOI     = [6, 8, 16]
# linesOI     = [5, 20]

#templateFileName = os.getcwd() + '/run_LT.sbatch'
templateFileName = os.getcwd() + '/input-LS-Temp.dat'

for struct in structNames:
  for leak in leakNames:

    leakNo = leak.split('P')[0]

    #fileName = os.getcwd() + '/Batch-HS/Batch-HS-%s-%s.sbatch' %(struct, leak)
    fileName = os.getcwd() + '/LS-Input-Files/input-LS-%s-%s.dat' %(struct, leak)
    shutil.copy(templateFileName, fileName)

    with open(templateFileName, 'r+') as tempFile:
      with open(fileName, 'w+') as file:

        lineContents   = tempFile.readlines()

        # outDir         = lineContents[linesOI[0]]
        # newOutDir      = outDir.split('-')
        # newOutDir[2]   = 'output=/projects/zair9172/LT/sbatch_out/HS-%s-%s.out' %(struct, leak)
        # newOutDir      = '-'.join(newOutDir) + '\n'

        # lineContents[linesOI[0]] = newOutDir

        # inputFile      = lineContents[linesOI[1]]
        # dirs           = inputFile.split('/')
        # dirs[10]       = 'HS-Input-Files'
        # inputFile      = '/'.join(dirs)
        # newInputFile   = inputFile.split('\n')[0] + 'input-HS-%s-%s.dat' %(struct, leak) + '\n'

        # lineContents[linesOI[1]] = newInputFile

        dataDir        = lineContents[linesOI[0]]
        newDataDir     = dataDir.split('/')
        newDataDir[3]  = struct
        newDataDir[4]  = leak
        newDataDir     = '/'.join(newDataDir)

        lineContents[linesOI[0]] = newDataDir

        dataTag        = lineContents[linesOI[1]]
        newDataTag     = dataTag.split('-')
        newDataTag[2]  = struct
        newDataTag[3]  = leakNo + 'pL' + '-' + '\n'
        newDataTag     = '-'.join(newDataTag)

        lineContents[linesOI[1]] = newDataTag

        outDir         = lineContents[linesOI[2]]
        newOutDir      = outDir.split('/')
        newOutDir[3]   = struct
        newOutDir[4]   = leak
        newOutDir      = '/'.join(newOutDir)
        newOutDir      = newOutDir.split('-')
        newOutDir[2]   = struct
        newOutDir[3]   = leakNo + 'pL'
        newOutDir      = '-'.join(newOutDir)

        lineContents[linesOI[2]] = newOutDir

        file.writelines(lineContents)

      file.close()
    tempFile.close()


