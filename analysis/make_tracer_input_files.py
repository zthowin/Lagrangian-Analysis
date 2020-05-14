import os
import shutil
import linecache

structNames = ['M1', 'M2', 'M3', 'M4', 'M7', 'M8']
leakNames   = ['0PL', '10PL', '20PL', '40PL']
locNames    = ['Clot', 'Wall', 'Drug']
linesOI     = [6, 8, 18]
# linesOI     = [5, 20]

# templateFileName = os.getcwd() + '/run_LT.sbatch'

for loc in locNames:
    templateFileName = os.getcwd() + '/input-HS-%s-Temp.dat' %loc

    for struct in structNames:
      for leak in leakNames:

        leakNo = leak.split('P')[0]

        #fileName = os.getcwd() + '/Batch-HS-Tracer/Batch-HS-%s-%s-%s.sbatch' %(loc, struct, leak)
        fileName = os.getcwd() + '/HS-Tracer-Input-Files/input-HS-%s-%s-%s.dat' %(loc, struct, leak)
        shutil.copy(templateFileName, fileName)

        with open(templateFileName, 'r+') as tempFile:
          with open(fileName, 'w+') as file:

            lineContents   = tempFile.readlines()

            # outDir         = lineContents[linesOI[0]]
            # newOutDir      = outDir.split('-')
            # newOutDir[2]   = 'output=/projects/zair9172/LT/sbatch_out/HS-%s-Tracer-%s-%s.out' %(loc, struct, leak)
            # newOutDir      = '-'.join(newOutDir) + '\n'

            # lineContents[linesOI[0]] = newOutDir

            # inputFile      = lineContents[linesOI[1]]
            # dirs           = inputFile.split('/')
            # dirs[10]       = 'HS-Tracer-Input-Files'
            # inputFile      = '/'.join(dirs)
            # newInputFile   = inputFile.split('\n')[0] + 'input-HS-%s-%s-%s.dat' %(loc, struct, leak) + '\n'

            # lineContents[linesOI[1]] = newInputFile

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
            newOutDir[3]   = struct
            newOutDir[4]   = leakNo + 'pL'
            newOutDir      = '-'.join(newOutDir)

            lineContents[linesOI[2]] = newOutDir

            file.writelines(lineContents)

          file.close()
        tempFile.close()


