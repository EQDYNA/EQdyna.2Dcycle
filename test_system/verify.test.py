#! /usr/bin/env python3
import os, sys, time
import numpy as np
import xarray as xr
try:
    from .testNameList import nameList
except ImportError:
    from testNameList import nameList

n_fail = 0
#testIDList   = ['tpv8', 'tpv104','tpv1053d','tpv1053d.6c','meng2023a','meng2023cb']
fileNameList = ['cyclelog.txt1','interval.txt1','totalop.txt1', 'totalop.txt2',
    'totalop.txt3', 'totalop.txt4', 'totalop.txt5', 'totalop.txt6', 'totalop.txt7',
    'totalop.txt8', 'totalop.txt9', 'totalop.txt10']
    
refRoot  = 'test_system/reference.results'
testRoot = 'work/test_results'

def compare_nc_files(fn1,fn2,threshold=1e-3):
    isTheSame = 'SUCCESS '+fn1+' '+fn2
    f1 = xr.open_dataset(fn1)
    f2 = xr.open_dataset(fn2)
    metadata_equal = f1.identical(f2)
    data_equal = (f1==f2).all().items()
    
    # Compare variables
    for var in f1.variables:
        var1 = f1[var]
        var2 = f2[var]
        if var1.dims != var2.dims:
            isTheSame = 'FAIL var dim '+fn1+' '+fn2
        if not np.allclose(var1,var2,rtol=threshold, atol=threshold):
            isTheSame = 'FAIL var numbers '+fn1+' '+fn2
        
    if not metadata_equal:# and data_equal:
        isTheSame = 'FAIL metadata '+fn1+' '+fn2
    
    try:
        xr.testing.assert_allclose(f1,f2)
        print('SUCCESS by xarray.testing.assert_allclose')
        print('SUCCESS '+fn1 +' '+fn2)
    except AssertionError as e:
        print(e)
    print(isTheSame)

    return isTheSame

def compare_txt_files(fn1,fn2,threshold=1e-3):
    isTheSame = 'SUCCESS '+fn1+' '+fn2
    with open(fn1,'r') as f1, open(fn2,'r') as f2:
        result1 = f1.read().split()
        result2 = f2.read().split()
    if len(result1) != len(result2):
        isTheSame = 'FAIL '+fn1+' '+fn2

    for num1,num2 in zip(result1,result2):
        fnum1, fnum2 = float(num1),float(num2)
        if abs(fnum1-fnum2) > threshold:
            isTheSame = 'FAIL '+fn1+' '+fn2
    print(isTheSame)
    return isTheSame

    
for testid in nameList:
    print(' ')
    for filename in fileNameList:
        refPath  = refRoot+'/'+testid+'/'+filename
        testPath = testRoot+'/'+testid+'/aRawSimuData/'+filename
        print(f'Checking {testid}: {filename}')
        if os.path.exists(refPath):
            if os.path.exists(testPath):
                if 'nc' in filename:
                    result = compare_nc_files(refPath, testPath, 1e-3)
                elif 'frt' in filename:
                    result = compare_txt_files(refPath, testPath, 1e-3)
                else:
                    result = compare_txt_files(refPath, testPath, 1e-3)
                if result and result.startswith('FAIL'):
                    n_fail += 1
            else:
                print(f'WARNING: Test output not found: {testPath}')
                n_fail += 1
        else:
            print(f'WARNING: Reference not found: {refPath}')

if n_fail > 0:
    print(f'\n❌ {n_fail} verification failure(s)')
    sys.exit(1)
print('\n✅ all verifications passed')

