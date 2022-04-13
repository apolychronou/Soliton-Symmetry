gcc  -Wall -O3 -x c axial_sym.c rk4.c functions.c -o axial -lm
wait
sshpass -f pass scp  axial kriti:./Documents
printf "Connected\n"
sshpass -f pass ssh kriti "rm Documents/data/*; cd Documents; ./axial; wait"
printf "Transfering Files\n"
sshpass -f pass scp kriti:./Documents/data/* ./data
printf "Complete\n"
