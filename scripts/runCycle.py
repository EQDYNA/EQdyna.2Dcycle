import os
import subprocess

# Set environment variables
os.environ["OMP_NUM_THREADS"] = "20"

# Create necessary directories
workdir = "."
os.makedirs(workdir, exist_ok=True)

# Initialize variables
numcycle = 1000
icycle = 1

# Create res.tar
subprocess.run(["tar", "-cvf", "res.tar", "eqdyna.2dcycle"])

# Loop over cycles
while icycle <= numcycle:
    print(f"Cycle #:{icycle}")
    with open("cycles.log", "a") as log_file:
        log_file.write(f"Cycle #:{icycle}\n")
        log_file.write(f"The {icycle} cycle started at\n")
        log_file.write(subprocess.check_output(["date"]).decode("utf-8"))

    # Build up stress
    if icycle == 1:
        subprocess.run(["./ini"])
    else:
        subprocess.run(["./bld"])

    os.remove("output4plot_dy.txt")
    os.makedirs(f"{workdir}/{icycle}", exist_ok=True)
    subprocess.run(["cp", "finalstress.txt", f"{workdir}/{icycle}"])
    subprocess.run(["mv", "stress4plot.txt", f"{workdir}/{icycle}"])

    # Simulate rupture
    subprocess.run(["./eqdyna.2dcycle"])
    subprocess.run(["cp", "output4plot_dy.txt", f"{workdir}/{icycle}"])
    subprocess.run(["mv", "nucloc.txt", f"{workdir}/{icycle}"])
    os.remove("finalstress.txt")

    # Update res.tar
    subprocess.run(["tar", "-rf", "res.tar", f"{workdir}/{icycle}"])
    subprocess.run(["tar", "-rf", "res.tar", "interval.time"])
    subprocess.run(["tar", "-rf", "res.tar", "output4plot_dy.txt"])
    subprocess.run(["tar", "-rf", "res.tar", "numofnucpoints.txt"])
    subprocess.run(["tar", "-rf", "res.tar", "cycles.log"])

    # Clean up
    subprocess.run(["rm", "-r", f"{workdir}/{icycle}"])

    with open("cycles.log", "a") as log_file:
        log_file.write(f"The {icycle} cycle ended at\n")
        log_file.write(subprocess.check_output(["date"]).decode("utf-8"))
        log_file.write("=====================================\n")

    icycle += 1

# Move log files
subprocess.run(["mv", "cycles.log", workdir])
subprocess.run(["mv", "interval.time", workdir])

print("Execution completed successfully.")