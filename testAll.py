#! /usr/bin/env python3
import os, time
from testNameList import nameList
# This script will perform tests on test cases defined in testNameList.

# No need for MPI for now.
# MPIRUN='mpirun.mpich' # please modify MPIRUN to fit your system accordingly. 
# print('testAll: MPIRUN is ', MPIRUN)
# print('testAll: please modify MPIRUN to fit your system accordingly.')

os.system('rm -rf test')
os.system('rm -rf bin/eqdyna.2dcycle')
os.system('mkdir test')

os.system('./install.eqdyna.2dcycle.sh -m ubuntu')
os.chdir('test')

startTime = time.time()
def runTest(testDir, compSet):
    cmd = 'create.newcase '+testDir+' '+compSet
    os.system(cmd)
    os.chdir(testDir)
    os.system('./case.setup')
    os.system('bash FE_run.sh')
    os.chdir('..')
    
for testName in nameList:
    runTest(testName, testName)

os.chdir('..')
os.system('python3 check.test.py')

endTime = time.time()

print('Time consumed for all the tests are ', endTime-startTime, ' s')
