import os

# numFiles = 180
# structNames = ['M9']#, 'M3', 'M4', 'M7', 'M8']
# leakNames   = ['0PL', '10PL', '20PL', '40PL']

# for struct in structNames:
#   for leak in leakNames:
#     leakNo = leak.split('P')[0]
#     for fileInd in range(90, 180):

#       if len(str(fileInd)) == 1:
#         add = '00'
#       elif len(str(fileInd)) == 2:
#         add = '0'
#       else:
#         add = ''
#       fname = '/%s/Velocity-HS-%s-%spL/Velocity-HS-%s-%spL000%s.vtu' %(struct, struct, leakNo, struct, leakNo, add + str(fileInd))
#       #fname = '/inputs-LS/Leakage/%s/%s/Velocity-LS-%s-%spL000%s.vtu' %(struct, leak, struct, leakNo, add + str(fileInd))
#       # fname = '/inputs-HS/Leakage/%s/%s/Velocity-HS-%s-%spL%s.vtu' %(struct, leak, struct, leakNo,str(fileInd))
#       #newfname = '/inputs-LS/Leakage/%s/%s/Velocity-LS-%s-%spL-%s.vtu' %(struct, leak, struct, leakNo,str(fileInd))
#       newfname = '/%s/Velocity-HS-%s-%spL/Velocity-HS-%s-%spL-%s.vtu' %(struct, struct, leakNo, struct, leakNo, str(fileInd))
#       os.rename(os.getcwd() + fname, os.getcwd() + newfname)
#       #os.remove(os.getcwd() + fname)

#     print("Leakage %s completed successfully.\n" %leak)
#   print("Structure %s completed successfully.\n" %struct)

for fileInd in range(90, 180):

  if len(str(fileInd)) == 1:
    add = '00'
  elif len(str(fileInd)) == 2:
    add = '0'
  else:
    add = ''

  fname = '/Velocity-HS-IM-V2/Velocity-HS-IM000%s.vtu' %(add + str(fileInd))
  newfname = '/Velocity-HS-IM-V2/Velocity-HS-IM-%s.vtu' %(str(fileInd))

  os.rename(os.getcwd() + fname, os.getcwd() + newfname)