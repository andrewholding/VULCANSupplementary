cd mergedReads
mkdir sorted

for bam in *bam
do
	samtools sort $bam -o sorted/$bam
done

cd sorted

for bam in *bam
do
	samtools index $bam
done

