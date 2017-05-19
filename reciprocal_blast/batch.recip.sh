#!/bin/bash
while IFS='$\n' read line; do
	pair=($line)
	a=${pair[0]}
	b=${pair[1]}
	echo "$a"
	echo "$b"
	qsub -q smallq -v "a=$a,b=$b" reciprocal.pbs
done
