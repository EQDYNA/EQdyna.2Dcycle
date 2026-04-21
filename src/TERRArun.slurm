#!/bin/bash
##SBATCH --export=NONE
##SBATCH --get-user-env=L


##
#SBATCH --job-name=3c-EQdyna2D-omp-P16ModD.slurm
###SBATCH -p mpi
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=3
#SBATCH --mem=30G
#SBATCH --account=122792052368
##SBATCH --partition=short
#SBATCH --output=job.%j
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dunyuliu@tamu.edu
###
echo "srun -l /bin/hostname"
srun -l /bin/hostname


#First Executable Line
tcsh runcycle.tcsh
#
#for ipatch in {39..39}
#do 
#  echo "Patch #:" $ipatch
#  echo "The x"$ipatch" patch started at" >>log.log
#  date >> log.log
#  rm patchnumber.txt
#  rm direction.txt
#  echo $ipatch >>patchnumber.txt
#  for idirection in {1..1}
#	do 
#  	echo $idirection >>direction.txt
#	done
#  mkdir "./x$((ipatch))"
#  mpirun -np 100 ./eqdyna-hyb
#  mv body* frt* fslip* "./x$((ipatch))"
#  
#  echo "y"$ipatch" started at" >>log.log
#  date >> log.log
#  rm direction.txt
#  for idirection in {2..2}
#	do
#  	echo $idirection >>direction.txt
#	done
#  mkdir "./y$((ipatch))"
#  mpirun -np 100 ./eqdyna-hyb
#  mv body* frt* fslip* "./y$((ipatch))"
#  
#  echo "z"$ipatch" patch started at" >>log.log
#  date >> log.log
#  rm direction.txt
#  for idirection in {3..3}
#	do
#  	echo $idirection >>direction.txt
#	done
#  mkdir "./z$((ipatch))"
#  mpirun -np 100 ./eqdyna-hyb
#  mv body* frt* fslip* "./z$((ipatch))"
#done
