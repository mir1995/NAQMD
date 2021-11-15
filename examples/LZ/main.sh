for s in 1 2 3 4 5 6 7 8 9 10
do 
    for npart in 10000000 100000000 1000000000
    do 
        sed -n '11p' mass_transitioned_seed${s}_npart$npart.txt | awk -v npart=$npart '{print npart, $2}' >> sa123_seed$s.txt
    done
done
