#! /bin/bash 

# The shell script is to set up environments for EQdyna.2Dcycle and 
#	install it. It will call the makefile inside src/ and generate 
#	an executable eqdyna.2dcycle and move it to bin/.

# Currently, the machines supported are:
#	ls6:    Lonestar6 at TACC
#	ubuntu: Ubuntu 22.04
#   macos:  MacOS Sonoma 14.7

# Usage: install.eqdyna.2dcycle.sh [-h] [-m Machine_name] [-e Machine_name] [-c Machine_name]

while getopts "m:e:c:h" OPTION; do
    case $OPTION in 
        m)
            MACH=$OPTARG
            ;;
        e) 
            MACH=$OPTARG
            ENV="True"
            ;;
        c)
            MACH=$OPTARG
            CONFIG="True"
            ;;
        h)
            echo "Usage: ./install-eqdyna.sh [-h] [-m Machine_name] [-c Machine_name]  "
            echo "                                                                     "
            echo "Examples:                                                            "
            echo "                                                                     "
            echo "install.eqdyna.2dcycle.sh -h                                         "
            echo " -----Display this help message                                      "
            echo "                                                                     "
            echo "install.eqdyna.2dcycle.sh -m ls6                                     "
            echo " -----Install EQdyna.2Dcycle on Lonestar6 at TACC                    "
            echo "                                                                     "
            echo "install.eqdyna.2dcycle.sh -e macos                                   "
            echo " -----Setup env and install EQdyna.2Dcycle on MacOS                  "
            echo "                                                                     "
            echo "source install.eqdyna.2dcycle.sh                                     "
            echo " -----Activate ENV VAR EQDYNA2DCYCLEROOT and add exes to PATH        "
            echo "                                                                     "
            echo "Currently supported machines include:                                "
            echo " ls6/ubuntu/grace                                                    "
            ;;
    esac
done 

if [ -n "$MACH" ]; then 
    export MACHINE=$MACH
    if [ $MACHINE == "ls6" ]; then 
        echo "Installing EQdyna.2Dcycle on Lonestar6 at TACC ... ..."
        if [ -n "$ENV" ]; then
            echo "Set up ENV on LS6 ... ..."
        fi 
    elif [ $MACHINE == "ubuntu" ]; then 
        echo "Installing EQdyna.2Dcycle on Ubuntu 22.04 ... ..."
    elif [ $MACHINE == "macos" ]; then 
        echo "Installing EQdyna.2Dcycle on MacOS ... ..."
        if [ -n "$ENV" ]; then
            echo "Setting up Env for EQdyna.2Dcycle on MacOS with Homebrew ... ..."
            echo "It is assumed Homebrew is installed either vis .pkg or cmd"
            brew install gcc python 
            pip3 install --break-system-packages numpy matplotlib xarray gmsh meshio nbconvert
        fi 
    elif [ $MACHINE == "grace" ]; then 
        echo "Installing EQdyna.2Dcycle on Grace at TAMU ... ..."
    fi 
    
    if [ -n "$CONFIG" ]; then 
        echo "Simply configure EQdyna without installation ... ..."
    else
        cd src
        make
        cd ..
        mkdir bin
        mv src/eqdyna.2dcycle bin
    fi

    export EQDYNA2DCYCLEROOT=$(pwd)
    export PATH=$(pwd)/bin:$PATH
    export PATH=$(pwd)/scripts:$PATH
    
    chmod -R 755 scripts
fi

export EQDYNA2DCYCLEROOT=$(pwd)
export PATH=$(pwd)/bin:$PATH
export PATH=$(pwd)/scripts:$PATH

echo EQDYNA2DCYCLEROOT
echo PATH 
