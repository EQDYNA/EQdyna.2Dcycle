#! /bin/bash 

# The shell script is to set up environments for EQdyna.2Dcycle and 
#	install it. It will call the makefile inside src/ and generate 
#	an executable eqdyna.2dcycle and move it to bin/.

# Currently, the machines supported are:
#	ls6:	Lonestar6 at TACC
#	ubuntu: Ubuntu 22.04

# Usage: install.eqdyna.2dcycle.sh [-h] [-m Machine_name] [-c Machine_name]

while getopts "m:c:h" OPTION; do
    case $OPTION in 
        m)
            MACH=$OPTARG
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
            echo "install.eqdyna.2dcycle.sh -c ubuntu                                  "
            echo " -----Simply set up envs for EQdyna.2Dcycle without installation     "
            echo " -----on ubuntu                                                      "
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
        
    elif [ $MACHINE == "ubuntu" ]; then 
        echo "Installing EQdyna.2Dcycle on Ubuntu 22.04 ... ..."
    
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
