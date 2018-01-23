for filename in ../GRHL2_filtered/*.bam 
do
 echo $filename
 output=`basename $filename`
 bamCoverage -b $filename -o $filename.bw
done
