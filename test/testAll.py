#! /usr/bin/env python3
"""
EQdyna.2Dcycle Testing System
Automated test runner for organized case structure
"""
import os, sys, time
try:
    from .testNameList import nameList
except ImportError:
    from testNameList import nameList

print("🧪 EQdyna.2Dcycle Testing System")
print("="*50)

# Change to project root if we're in test subdirectory
if os.path.basename(os.getcwd()) == 'test':
    os.chdir('..')
    print(f"📁 Changed to project root: {os.getcwd()}")

# Clean previous test results
print("🧹 Cleaning previous test results...")
os.system('rm -rf work/test_results')
os.system('rm -rf bin/eqdyna.2dcycle')

# Build fresh binary
print("🔨 Building EQdyna.2Dcycle...")
build_result = os.system('./install.eqdyna.2dcycle.sh -m macos')  # Changed to macos
if build_result != 0:
    print("❌ Build failed! Exiting.")
    sys.exit(1)

# Set environment variables  
print("🌍 Setting environment variables...")
os.environ['EQDYNA2DCYCLEROOT'] = os.getcwd()
os.environ['PATH'] = f"{os.getcwd()}/bin:{os.getcwd()}/scripts:{os.environ['PATH']}"

# Create test results directory
os.makedirs('work/test_results', exist_ok=True)
os.chdir('work/test_results')

startTime = time.time()

def runTest(testDir, compSet):
    """Run a complete test case with organized workflow"""
    print(f"\n📋 Running test: {testDir}")
    print("-" * 30)
    
    # Create new case using organized structure
    print(f"  🏗️  Creating case from {compSet}...")
    cmd = f'python3 ../../scripts/create.newcase {testDir} {compSet}'
    result = os.system(cmd)
    if result != 0:
        print(f"  ❌ Case creation failed for {testDir}")
        return False
        
    os.chdir(testDir)
    
    # Generate mesh (if needed - organized structure may have pre-generated mesh)
    print("  🌐 Generating mesh...")
    mesh_result = os.system('python3 meshgen.py')
    if mesh_result != 0:
        print(f"  ❌ Mesh generation failed for {testDir}")
        os.chdir('..')
        return False
    
    # Setup case configuration
    print("  ⚙️  Setting up case...")
    setup_result = os.system('python3 case.setup')
    if setup_result != 0:
        print(f"  ❌ Case setup failed for {testDir}")
        os.chdir('..')
        return False
    
    # Run simulation with updated script
    print("  🚀 Running simulation...")
    sim_result = os.system('bash run.sh')  # Updated from FE_run.sh to run.sh
    if sim_result != 0:
        print(f"  ❌ Simulation failed for {testDir}")
        os.chdir('..')
        return False
    
    print(f"  ✅ Test {testDir} completed successfully")
    os.chdir('..')
    return True

# Run all tests
print(f"\n🔄 Running {len(nameList)} test case(s):")
success_count = 0
for testName in nameList:
    if runTest(testName, testName):
        success_count += 1

# Return to project root for validation
os.chdir('../..')

# Run validation
print(f"\n🔍 Running result validation...")
validation_result = os.system('python3 test/verify.test.py')

endTime = time.time()
elapsed = endTime - startTime

# Summary
print(f"\n📊 Test Summary:")
print(f"   Tests run: {len(nameList)}")
print(f"   Successful: {success_count}")  
print(f"   Failed: {len(nameList) - success_count}")
print(f"   Time elapsed: {elapsed:.2f} seconds")

if success_count == len(nameList):
    print("✅ All tests passed!")
else:
    print("❌ Some tests failed!")
    sys.exit(1)
