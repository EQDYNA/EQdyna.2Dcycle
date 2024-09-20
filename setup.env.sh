#! /bin/bash 

# The shell script is to set up computing environments for EQdyna.2Dcycle.

# Currently, the machines supported are:
#	ls6:	Lonestar6 at TACC
#	ubuntu: Ubuntu 22.04

# Usage: setup.env.sh [Machine_name]


echo "Usage: ./setup.env.sh [Machine_name]  "
echo "                                                                     "
echo "Examples:                                                            "
echo "                                                                     "
echo "setup.env.sh                                                         "
echo " -----Display this help message                                      "
echo "                                                                     "
echo "setup.env.sh macos                                                   "
echo " -----Install EQdyna.2Dcycle on Lonestar6 at TACC                    "
echo "                                                                     "
echo "setup.env.sh ubuntu                                  "
echo " -----Simply set up envs for EQdyna.2Dcycle without installation     "
echo " -----on ubuntu                                                      "
echo "                                                                     "
echo "Currently supported machines include:                                "
echo " ls6/ubuntu/macos                                                    "



export MACHINE=$1
if [ $MACHINE == "ls6" ]; then 
    echo "Installing EQdyna.2Dcycle on Lonestar6 at TACC ... ..."
    
elif [ $MACHINE == "ubuntu" ]; then 
    echo "Installing EQdyna.2Dcycle on Ubuntu 22.04 ... ..."

elif [ $MACHINE == "macos" ]; then 
    echo "Installing EQdyna.2Dcycle on MacOS with Homebrew ... ..."
    echo "It is assumed Homebrew is installed either vis .pkg or cmd"
    brew install gcc python 
    pip install numpy matplotlib xarray
fi 

export EQDYNA2DCYCLEROOT=$(pwd)
export PATH=$(pwd)/bin:$PATH
export PATH=$(pwd)/scripts:$PATH

chmod -R 755 scripts

export EQDYNA2DCYCLEROOT=$(pwd)
export PATH=$(pwd)/bin:$PATH
export PATH=$(pwd)/scripts:$PATH

echo EQDYNA2DCYCLEROOT
echo PATH 
