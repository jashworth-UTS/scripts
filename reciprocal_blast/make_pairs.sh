#!/bin/bash
for p1 in $*; do
	for p2 in $*; do
		if [ "$p1" == "$p2" ]; then
			continue
		fi
		echo "$p1	$p2"
	done
done
