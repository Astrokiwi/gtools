for i in 64 32 16 8 4 2
do
    export OMP_NUM_THREADS=$i
    #echo $OMP_NUM_THREADS
    echo "nprocs=" $i
    ./rebeamcomp
done
