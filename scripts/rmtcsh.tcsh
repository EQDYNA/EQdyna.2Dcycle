set numcycle = 1000
set icycle = 1
while ($icycle <= $numcycle)
  rm -r $icycle
  @ icycle++
end
rm interval.time
rm cycles.log
rm numofnucpoints.txt
rm nucloc.txt
rm finalstress.txt
