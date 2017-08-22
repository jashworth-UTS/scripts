q=smallq
pbs=rpsblast.pbs
db=Cdd

echo "pbs script:"
echo "[BEGIN SCRIPT]"
cat $pbs
echo "[END SCRIPT]"
echo ""

for f in $*; do
	echo "qsub $pbs on $q for $f"
	qsub -q $q -v "fa=$f,db=$db" $pbs
done

# array job style
#in=($@)
#fas="fas"
#db="Cdd"
#printf "%s\n" "${in[@]}" > $fas
##echo $* > fas
#j=$(wc -l < $fas)
#echo "$j fasta queries $fas being searched against db $db"
#cat $fas
#qsub -q workq -J 1-$j -v "fas=$fas,db=$db" rpsblast.pbs
