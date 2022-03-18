#! bin/bash

IFS=

# This linux shell script requires VMD to be installed and is optimized to extract the coordinates of Desmond generated trajectories

base=$(pwd)
MD_simulations=./MD_simulations # Folder of the MD simulations
traj_dir=$(cd "$MD_simulations"; pwd)

frames=1000   ### Set the number of frames ###

  if [ \! -d Extracted_coordinates ]; then
  mkdir Extracted_coordinates
  fi

# Creating the evaluation folders for different systems (beginning of the outer cycle). Trajectory folders should be ./System_name_MD/System_name_MD_simulation_length/System_name_MD_simulation_length-out.cms and ./System_name_MD/System_name_MD_simulation_length/System_name_MD_simulation_length_trj/clickme.dtr

  for system in SH2_wt; do ### Set the system names ###

  if [ \! -d $system\_MD ]; then
  mkdir $system\_MD
  fi

  cd $system\_MD

  # Creating the evaluation folders of different runs (beginning of the inner cycle start)

  for ns in 10; do ### Set the simulation lengths ###

    if [ \! -d $system\_MD_$ns ]; then
    mkdir $system\_MD_$ns
    fi

    cd $system\_MD_$ns

      #Setting the starting frame
      start=1

      # Extractin the backbone coordinates in mol2 format using vmd

	echo "Extracting coordinates of System: $system Length: $ns"

      bash /home/user/Program_files/VMD/vmd.sh $traj_dir/$system\_MD/$system\_MD_$ns/$system\_MD_$ns\-out.cms $traj_dir/$system\_MD/$system\_MD_$ns/$system\_MD_$ns\_trj/clickme.dtr -dispdev none -e $base/Desmondtraj_to_mol2.tcl -args $traj_dir/$system\_MD/$system\_MD_$ns/$system\_MD_$ns\_backbone $frames > vmd.log ### Set here the path to vmd.sh ###

      # Creating list of the atoms and cordinates

      grep ' C \|CA \| N ' $traj_dir/$system\_MD/$system\_MD_$ns/$system\_MD_$ns\_backbone.mol2 | awk '{print $6"-"$7"-"$8" "$3","$4","$5}' > List_of_coordinates.dat

      grep ' C \|CA \| N ' $traj_dir/$system\_MD/$system\_MD_$ns/$system\_MD_$ns\_backbone.mol2 | awk '{print $6"-"$7"-"$8}' | sort | uniq > List_of_atoms.dat

      test -f Data_all_numbered.dat && rm Data_all_numbered.dat

      #Data sorting and numbering
      
      echo "Creating coordinate file of System: $system Length: $ns"

      while read p; do

         grep "$p" List_of_coordinates.dat | grep -n "$p" | awk '{print " "$1" "$2" "$3}' >> Data_all_numbered.dat

         done <List_of_atoms.dat

      sed -i -r 's/:/ /g' Data_all_numbered.dat

      #Transposing the data into table form

      test -f $system\_MD_$ns\_final_formatted.dat && rm $system\_MD_$ns\_final_formatted.dat

      cat List_of_atoms.dat | awk '{print $0"-X,"$0"-Y,"$0"-Z,"}' | rs -T | awk '{print "Frame "$0}' > $base/Extracted_coordinates/$system\_MD_$ns\_final_formatted.dat

      for (( i=$start; i<=$frames; i++ )); do

         grep " $i " Data_all_numbered.dat | awk '{print $3}' | rs -T | awk -v Frame=$i '{print 'Frame'" "$0}' >> $base/Extracted_coordinates/$system\_MD_$ns\_final_formatted.dat

         done

      sed -i 's/,/ /g' $base/Extracted_coordinates/$system\_MD_$ns\_final_formatted.dat

      # Cleaning up (optional)

      rm List_of_coordinates.dat

      rm List_of_atoms.dat

      rm Data_all_numbered.dat

      rm vmd.log

      cd ..

      done # End of the inner cycle

  rm -r $system\_MD_$ns

  cd ..

rm -r $system\_MD

done # End of the outer cycle
