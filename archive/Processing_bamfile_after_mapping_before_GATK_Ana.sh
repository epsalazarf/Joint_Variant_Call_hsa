### Written by Federico Sanchez starting from previous script by MCAA
## Modifications by Ana Hernandez

#Your job name
#$-N Proccss_bamfile
#Use current working directory
#$-cwd

#Join stdout and stderr
#$-j y
# Run job through bash shell
#$-S /bin/bash

#Send an email after the job has finished
#$-m e
#$-M annhled@gmail.com

#If modules are needed, source modules environment (Do not delete the next line):
 . /etc/profile.d/modules.sh

#change email to receive results there 
EMAIL='annhled@gmail.com'

####Preparing data: Cleaning and reporting stats
module load htslib/1.2.1
module load gcc/5.1.0
module load bwa/0.7.13
module load samtools/1.9
module load picard/2.6.0
module load gatk/3.4-46
module load fastqc/0.11.3
module load bcftools/1.9
module load mosdepth/0.3.3

## Start analysis
export base=$1; ## This is sample name
echo $base
echo -e "### This is the report for sample "${base}"  ####" > ${base}_REPORT.txt

## Sort mapped reads with samtools
#samtools sort -n -o ${base}.sort.bam  ${base}.bam   
#echo "${base}.sort.bam"
#countR=`samtools view -F 4 -q 30 -c ${base}.sort.bam`
#echo $countR
#echo "Counts after sort  $countR"
#echo -e "##Unprocessed bamfile \t ${base} \t #reads \t "${countR}""  >> ${base}_REPORT.txt

## Remove PCR clones ### se adapto a markdup por version usada de samtools
#samtools fixmate -m -O BAM ${base}.sort.bam ${base}.sortf.bam ##add ms and MC tags for markdup to use later 
#samtools sort -o ${base}.sorted.bam ${base}.sortf.bam ## Markdup needs position order
#samtools markdup -r ${base}.sorted.bam ${base}.sort.rmdup.bam 
#echo "${base}.sorted.rmdup.bam"
#countR=`samtools view -F 4 -q 30 -c ${base}.sort.rmdup.bam`
#echo "Counts after markdup $countR"
#echo -e "##After markdup removing \t ${base} \t #reads \t "${countR}"" >> $base\_REPORT.txt

## Filter reads according to mapping quality (30)
#samtools view -q 30 ${base}.sort.rmdup.bam -o ${base}.sort.rmdup.qc.bam
#countR=`samtools view -c ${base}.sort.rmdup.qc.bam`
#echo "Counts after mapping quality filtering $countR"
#echo -e "##After mapping quality filtering \t ${base} \t #reads\t "${countR}"" >> ${base}_REPORT.txt

## Create index for bamfile
samtools index ${base}.sort.rmdup.qc.bam ${base}.sort.rmdup.qc.bam.bai

## Removing files for saving memory space
#rm -f ${base}.sort.rmdup.bam
#rm -f ${base}.sortf.bam
#rm -f ${base}.sorted.bam 
#rm -f ${base}.sort.bam
#echo "files removed" 

#Run fastqc of mapped reads and Compute de average lenght of reads
#fastqc -f bam ${base}.sort.rmdup.qc.bam &
#echo "fastqc ran on mapped reads"

## Calculation of genotype likelihoods at each genomic position with coverage 
export hg38="../../../alhernandez/Lupus/WGS/hg38.fa"
echo $hg38
echo ${base}.sort.rmdup.qc.bam

cov_st=`samtools mpileup -f $hg38 -C 50 ${base}.sort.rmdup.qc.bam | tee ${base}.sort.rmdup.qc.pileup | awk '{if($3!="*"){sum+=$4}}END{print  sum/NR "\n" NR}'`    
echo -e "##Coverage mpileup samtools \t ${base} \t "${cov_st}"" >> ${base}_REPORT.txt
echo "mpile up ran $cov_st"

#cov=`bcftools mpileup -f $hg38 -C 50 ${base}.sort.rmdup.qc.bam | tee ${base}.sort.rmdup.qc.pileup | awk '{if($3!="*"){sum+=$4}}END{print sum/NR "\n" NR}'`
#echo $cov
#echo -e "##Coverage mpileup bcftools \t ${base} \t "${cov}"" >> ${base}_REPORT.txt
#echo -e "##Mapped reads average length" >> ${base}\_REPORT.txt

## Calculation of coverage with mosdepht from samtools 
mosdepth -t 4 --fast-mode -n -f $hg38 $base ${base}.sort.rmdup.qc.bam     
cov_st2=`awk '$1=="total"{print $4;}' ${base}.mosdepth.summary.txt`;
echo -e  "##MosDepth coverage samtools \t ${base} \t "${cov_st2}"" >> ${base}_REPORT.txt
echo "mosdepth ran " 

## Calculation of coverage with mosdepht from bcftools
#mosdepth -t 4 --fast-mode -n -f $hg38 $base $base.sort.rmdup.qc.bam
#cov=`awk '$1=="total"{print $4;}' $base.mosdepth.summary.txt`;
#echo -e  "##MosDepth coverage bcftools \t ${base} \t "${cov}"" >> ${base}_REPORT.txt

## Calculation of the real lenght of the reads
#RL=`samtools view ${base}.sort.rmdup.qc.bam |awk '{SUM+=length($10);DIV++}END{print SUM/DIV}'`
#echo $RL
#echo -e "##Real Lenght \t ${base} \t "${RL}"" >> ${base}_REPORT.txt

## Make index
#samtools index ${base}.sort.rmdup.qc.bam

## Replace readgroups with ID from sequencing
#java -jar $PICARDHOME/picard.jar AddOrReplaceReadGroups \
 #     I=${base}.sort.rmdup.qc.bam\
 #     O=${base}.sort.rmdup.qc.RG.bam\
  #    RGLB=${base} 
