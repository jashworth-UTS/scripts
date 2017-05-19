#!/bin/bash

queue=c3b
queue=smallq

./make_pairs.sh $* > pairs
j1=`wc -l < pairs`
echo "$j1 pairs"

./unique_pairs.py $* > unique_pairs
j2=`wc -l < unique_pairs`
echo "$j2 unique pairs"

j1="1-$j1"
j2="1-$j2"

# indepedent one-directional blasts
blastp_best=$(qsub -q $queue -J $j1 -v pairs=pairs blastp_best.pbs.array)


# now figure out unique pairs and process reciprocally as such
qsub -W depend=afterok:$blastp_best -q $queue -J $j2 -v pairs=unique_pairs reciprocal.pbs
