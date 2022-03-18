set outfile [lindex $argv 0]
append outfile ".mol2"
set frames [lindex $argv 1]

animate write mol2 $outfile sel [atomselect 0 "name C CA N"] beg 1 end $frames skip 1 0
exit
