for s in 1 2 3 4 5 6 7 8 9 10
do
    for rate in sa3multid lzdia lzadia
    do 
        ./jahn_teller_dynamics $rate $s &
    done
done
