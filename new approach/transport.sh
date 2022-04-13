sshpass -f pass scp -r src kriti:./Documents/axial_symmetry/
wait
printf "Source files transfered\n"
sshpass -f pass scp CMakeLists.txt kriti:./Documents/axial_symmetry/
wait
printf "Cmake list transfered\n"
sshpass -f pass ssh kriti "rm Documents/data/*; cd Documents/axial_symmetry/build; cmake .. ; make ; ./axialSym; wait"
printf "Execution complete, transfering results back \n"
sshpass -f pass scp kriti:./Documents/data/* ./data
printf "Complete\n"
