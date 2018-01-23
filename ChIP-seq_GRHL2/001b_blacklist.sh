#mkdir blacklists
#cd blacklists
#wget http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg38-human/hg38.blacklist.bed.gz
#gunzip hg38.blacklist.bed.gz
#cd ..


bl=../../blacklists/hg38.blacklist.bed

mkdir ./GRHL2_filtered
cd mergedReads/sorted/
for f in *.bam
do
	echo $f
	bedtools intersect -v -abam $f -b $bl  > ../../GRHL2_filtered/$f
done

### Re-index

cd ../../GRHL2_filtered
for f in *.bam
do
samtools index $f
done

cd ..

