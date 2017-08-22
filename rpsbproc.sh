q=workq
pbs=rpsbproc.pbs

echo "pbs script:"
echo "[BEGIN SCRIPT]"
cat $pbs
echo "[END SCRIPT]"
echo ""

for f in $*; do
	echo "qsub $pbs for $f"
	qsub -q $q -v "xml=$f" $pbs
done

# array job style
#rm cmds
#for x in $*; do
#	echo "rpsbproc -i $x -m full -t -e 1e-2 -o $x.proc" >> cmds
#done
#
#j=$(wc -l < cmds)
#echo "$j cmds"
#cat $cmds
#qsub -q workq -J 1-$j -v "cmds=cmds" rpsbproc.pbs
