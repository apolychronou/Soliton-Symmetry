sshpass -f pass scp -r src kriti:./Documents/axial_symmetry/src
wait
printf "Source files transfered\n"
sshpass -f pass ssh kriti "rm Documents/data/*; cd Documents/axial_symmetry/build; cmake ..; wait ; make ; wait; ./axialSym; wait"
printf "Execution complete, transfering results back \n"
sshpass -f pass scp kriti:./Documents/data/* ./data
printf "Complete\n"
