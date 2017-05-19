for f in `cat blastps`; do qsub -q smallq -v "f=$f" best.pbs ; done
