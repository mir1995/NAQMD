for s in 1 2 3 4 5 6 7 8 9 10
do
    for rate in lzdia lzadia sa3multid
    do
        while [ $(jobs | wc -l) -ge 10 ] 
        do
            sleep 1
        done 
        ./jahn_teller_dynamics $rate $s &
    done
done
