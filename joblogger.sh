#!/bin/bash

# gets current PBSPro qstat -f info every [N] seconds, pipes to a Python parser, concats to a running stats file for later analysis

rm jobstats

while true; do
	qstat -f | ./jobstats.py >> jobstats
	sleep 10
done
