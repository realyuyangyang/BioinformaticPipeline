#!/bin/bash
raw_fq_dir=$1
sra_list=$2

trim_galore=/data/nas2/software/miniconda3/envs/ref_trans/bin//trim_galore
cutadapt=/data/nas2/software/miniconda3/envs/ref_trans/bin/cutadapt
fastqc=/data/nas2/software/miniconda3/envs/snakemake/bin/fastqc
multiqc=/data/nas2/software/miniconda3/bin/multiqc
samtools=/data/nas2/software/miniconda3/envs/ref_trans/bin//samtools
hisat2=/data/nas2/software/miniconda3/envs/ref_trans/bin//hisat2
ss_file=/data/nas1/zhiyq/project/02_BJTC116/00_rawdata/04_Mouse_genome/mm39.ss
featurecounts=/data/nas2/software/miniconda3/envs/ref_trans/bin/featureCounts
hisat_index_dir=/data/nas1/zhiyq/project/02_BJTC116/00_rawdata/04_Mouse_genome/mm39
ref_gtf=/data/nas1/zhiyq/project/02_BJTC116/00_rawdata/04_Mouse_genome/mm39.ncbiRefSeq.gtf

# 01 trim galore
trim_dir="01_trim"
if [[ ! -d $trim_dir ]]
    then
        mkdir $trim_dir
    fi
find $raw_fq_dir -name "*gz" | while read f
    do
        echo "trim_galore work started at [$(date)] for ${f}"
        $trim_galore --quality 20 \
                     -j 4 \
                     --fastqc \
                     --path_to_cutadapt $cutadapt \
                     -o $trim_dir \
                     ${f}
        echo "trim_galore work finished at [$(date)] for ${f}"
        echo "=========================================================="
    done

# 02 fastqc
fastqc_dir="02_fastqc"
if [[ ! -d $fastqc_dir ]]
    then
        mkdir $fastqc_dir
    fi
find $trim_dir -name "*gz" | while read f
    do
        echo "fastqc work started at [$(date)] for ${f}"
        $fastqc --outdir $fastqc_dir \
                --threads 16 \
                ${f}
        $multiqc $fastqc_dir/*zip -o ${fastqc_dir}/multiqc
        echo "fastqc work finished at [$(date)] for ${f}"
        echo "=========================================================="
    done

# 03 hisat align
bam_dir="03_bam"
if [[ ! -d $bam_dir ]]
    then
        mkdir $bam_dir
    fi
cat $sra_list | while read f
    do
        echo "[hisat2 align] sample ${f} start at [$(date)]"
        $hisat2 -p 16 \
                -q \
                --known-splicesite-infile ${ss_file} \
                -x ${hisat_index_dir} \
		        -U ${trim_dir}/${f}_trimmed.fq.gz \
                -S ${bam_dir}/${f}.sam
	    ${samtools} sort -@ 16 -o ${bam_dir}/${f}.sort.bam ${bam_dir}/${f}.sam
	    ${samtools} index -@ 16 ${bam_dir}/${f}.sort.bam
        echo "[hisat2 align] sample ${f} finished at [$(date)]"
        echo "=========================================================="
    done

# 04 read count
count_dir="04_count"
if [[ ! -d $count_dir ]]
    then
        mkdir $count_dir
    fi
${featurecounts} -T 16 \
                -g gene_id \
                -a ${ref_gtf} \
                -o ${count_dir}/final_counts.txt \
                ${bam_dir}/*sort.bam

sed -i -e 's![/]*\S*/!!g' -e 's/\.sort\.bam//g' ${count_dir}/final_counts.txt
