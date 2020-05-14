import vtk
import numpy as np
import sys

def makeReaderObject(a_FileName):

  if a_FileName.endswith('vtk'):
    reader = vtk.vtkUnstructuredGridReader()
  elif a_FileName.endswith('vtu'):
    reader = vtk.vtkXMLUnstructuredGridReader()

  return reader

def extractFlowDataFromFile(a_FileName, a_ReaderObj, a_DataName=None):

    a_ReaderObj.SetFileName(a_FileName)
    a_ReaderObj.Update()

    if a_DataName is not None:
        return a_ReaderObj.GetOutput().GetPointData().GetArray(a_DataName)
    else:
        return a_ReaderObj.GetOutput()

def getMeshVelocity(a_FileName, a_ReaderObj, a_DataName=None):

  data       = extractFlowDataFromFile(a_FileName, a_ReaderObj, a_DataName)
  numNodes   = data.GetNumberOfTuples()
  velocities = np.zeros((numNodes, 3), dtype=np.float64, order='F')

  for p in range(numNodes):
    velocities[p,:]    = data.GetTuple(p)[0], data.GetTuple(p)[1], data.GetTuple(p)[2]

  return velocities


numFiles = 91
prefix   = '/home/zach/Documents/Mukherjee/Lagrangian-Toolkit/Code/inputs-LS/FDLS-Velocity-._'

for fileInd in range(numFiles):

  print(fileInd)

  saveFile = prefix + '%i.npy'%fileInd
  vtuFile  = prefix + '%i.vtu'%fileInd
  print(vtuFile)

  if fileInd == 0:
    readerObj = makeReaderObject(vtuFile)

  velocityArr = getMeshVelocity(vtuFile, readerObj, a_DataName='f_11-0')
  np.save(saveFile, velocityArr)
