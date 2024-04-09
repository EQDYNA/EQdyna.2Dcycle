set numcycle = 1000
set icycle = 1
set workdir = .
#setenv PATH /usr/local/mpich1/bin:$PATH
setenv OMP_NUM_THREADS 20
#echo "Interval time between events in years" >> interval.time
tar -cvf res.tar eqdyna2drun-omp
while ($icycle <= $numcycle)
  echo "Cycle #:" $icycle
  echo "The" $icycle "cycle started at" >> cycles.log
  date >> cycles.log
  # build up stress
  # echo "Cycle #:" $icycle >> interval.time
  if ($icycle == 1) then
    ./ini
  else
    ./bld
  endif 
  rm output4plot_dy.txt

  mkdir $workdir/$icycle
  \cp finalstress.txt $workdir/$icycle
  \mv stress4plot.txt $workdir/$icycle
  # simulate rupture
  #foreach size (30)
  #  setenv OMP_NUM_THREADS ${size}
    ./eqdyna2drun-omp
  #end
  \cp output4plot_dy.txt $workdir/$icycle
  mv nucloc.txt $workdir/$icycle
  rm finalstress.txt
  
  tar -rf res.tar $workdir/$icycle
  tar -rf res.tar interval.time
  tar -rf res.tar output4plot_dy.txt
  tar -rf res.tar numofnucpoints.txt
  tar -rf res.tar cycles.log

  rm -r $workdir/$icycle 

  echo "The" $icycle "cycle ended at" >> cycles.log
  date >> cycles.log
  echo "=====================================" >> cycles.log
  @ icycle++
end
mv cycles.log $workdir
mv interval.time $workdir
