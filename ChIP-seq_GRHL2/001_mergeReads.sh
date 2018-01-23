mkdir ./mergedReads

cd ./mergedReads


for bam in $(find ../SLX-14333/SLX-14333.D???_D???.HM3WHBBXX.s_2.*.bam -type f);
do
	
	bamOut=${bam/.HM3WHBBXX.s_2.bwa.homo_sapiens/}	
	bam2=${bam/HM3WHBBXX/HM3WTBBXX}
	bam2=${bam2/s_2/s_4}
	samtools merge `basename $bamOut` $bam  $bam2  
done




