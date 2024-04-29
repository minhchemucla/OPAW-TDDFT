#!/bin/bash


for nvtot in 32 64 128
do
    for nxch in 3 9 27
    do
	rm dim.inp
	cp dim1.inp tmp.inp
	sed "22s/32/$nvtot/g" tmp.inp > dim.inp 
	cp dim.inp tmp.inp
	sed "34s/27/$nxch/g" tmp.inp > dim.inp 
	grep nxch dim.inp
	grep 'nvtot ' dim.inp

	for i in `seq 1 5`;
	do
	    cp fort.4_$i fort.4
	    cat fort.4
	    # cp fort.4 tmp_x$nxch.I$nvtot.r$i  
	    /home/samuelh/app-src/openmpi-1.6/build-lahey/bin/mpirun -n 32 ./a.out > logq
	    mv logq log_x$nxch.I$nvtot.r$i  
	done
    done
done





