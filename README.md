# Exploring differences in sequencing quality and variant-detectionconcordance across DNA sources

### 22126 Next Generation Sequencing Analysis

*Mariana Chichkova (s205694), Mette Christoffersen (s192033), Laura Sans-Comerma (s192437),Natasia Thornval (s143493), and Huijiao Yang (s202360*



### What can I find here?

| What           | Where                              |
| -------------- | ---------------------------------- |
| Plots          | [R](https://github.com/laurasansc/NGS_group_7/blob/main/plots/plots.R)  [Python](https://github.com/laurasansc/NGS_group_7/blob/main/plots/venn.py)                  |
| FastQC         | [Pre-trimming Forward](https://github.com/laurasansc/NGS_group_7/blob/main/FastQC/SRR8595497_1_fastqc.png) [Post-trimming]() |
| Samtools STATS | [Stats](https://github.com/laurasansc/NGS_group_7/blob/main/STATS/SRR8595497_merged_markdup.txt)                                   |
| BAMstats       | [Coverage](https://github.com/laurasansc/NGS_group_7/blob/main/STATS/SRR8595497_markup_plot-coverage.png) , [GC content]( https://github.com/laurasansc/NGS_group_7/blob/main/STATS/SRR8595497_markup_plot-gc-content.png) , [More](https://github.com/laurasansc/NGS_group_7/tree/main/STATS)                                |
| Poster         | [PDF]()                            |

* To see the whole FastQC analysis, download the `.html` file and open in browser. 

## Commands and plots

**Downloading the data set**

```bash
wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR859/008/SRR85954XX/SRR85954XX_Y.fastq.gz
```

### SEQUENCING DATA PRE-PROCESSING AND MAPPING

**Quality check - FastQC**

```bash
fastqc -o ./fastqc SRR85954XX_X.fastq.gzTo open the fastqc files -> firefox SRR85954XX_X_fastqc.htmld
```

**Trimming (adapter sequences)** 
*leeHom used to infer adapter sequences*

```bash
leeHom --auto -fq1 SRR85954XX_1.fastq.gz -fq2 SRR85954XX_2.fastq.gz -fqo SRR85954XX_trimmed
Adapter sequences inferred:Inferred fwd adapter: AGATCGGAAGAGCACACGTCTGAACTCCAGTInferred rev adapter: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGT
```

*fastp used for adapter trimming - to speed up the process*

```bash
fastp --thread 2 --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGT --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGT -i /data/shared/groups/group_7/data/SRR85954XX_1.fastq.gz -I /data/shared/groups/group_7/data/SRR85954XX_2.fastq.gz -o /data/shared/groups/group_7/data/trimmed/SRR85954XX_1_trim.fastq.gz -O /data/shared/groups/group_7/data/trimmed/SRR85954XX_2_trim.fastq.gz
```

**Quality check - FastQC after trimming of adapters**

```bash
fastqc -o ./fastqc SRR85954XX_X_trim.fastq.gzTo open the fastqc files -> firefox SRR85954XX_X_trim_fastqc.html
```

**Alignment - BWA MEM**

Human genome reference : `/path/in/cluster/data/references/human/human_g1k_v37.fasta`

```bash
bwa mem -t 2 ref.fasta read1.fq read2.fq | samtools view -uS - | samtools sort /dev/stdin > SRR85954XXXX_aln.bam
```

**NB!** Since the mapping was taking > 48 hours, 210 temporary files were merged and the BAM subsampling step was skipped. 

**Temporary BAM mapping files merging**

Use of AddRG for running GATK smoothly

```bash
samtools merge --threads 4 -f  /dev/stdout trimmed/samtools.1301.4578.tmp.0{000..210}.bam | /path/to/Software/libbam/addRG /dev/stdin SRR85954XX_merged.bam SRR85954XX
```

**Index BAM files**

```bash
samtools index SRR85954XX.bam
```

**Alignment statistics**

```bash
samtools stat SRR85954XX.bamplot-bamstats -p SRR85954XX.stat
```


For average coverage per sample: 

```bash
 mosdepth SRR85954XX SRR85954XX.bam
```

**BAM file subsampling - not done**

Files should be subsampled to the level of the sample with the lowest original mean depth of coverage. The proportion of aligned reads to retainfrom each sample was calculated as 𝐷(M)/𝐷(X), where 𝐷(M) is the minimum original mean read depth among all the samples and 𝐷(X) is the original mean read depth of sample 𝑋.

```bash
samtools view -s 22.[fraction] -b SRR85954XX.bam > SRR85954XX_sub.bam
```

where 22 is a seed used for randomness of the selection.

**Mark duplicate reads - Picard **

```bash
java -Xmx10g -jar /usr/local/bin/picard.jar MarkDuplicates -I SRR85954XX_merged.bam -M SRR85954XX_sub_markdup.metrics.txt -O SRR85954XX_sub_markdup.bam
```

**Remove duplicate reads - samtools rmdup**

```bash
samtools rmdup SRR95854XX_merged.bam SRR95854XX_merged_rmdup.bam
```

**Index BAM files**

```bash
samtools index SRR85954XX_merged.bam
```

**Alignment statistics - when the duplicates are removed or marked**

```bash
samtools stat SRR85954XX.bamplot-bamstats -p SRR85954XX.stat
```

### VARIANT CALLING AND OTHERS

**Variant calling - GATK**

HaplotypeCaller is used to identify and annotate SNPs and indels. The output from this command is a gVCF-file (genomic VCF). The gVCF file format is a special VCF-file, which contains genotype likelihoods for all positions in the genome - opposed to regular VCF files, which only include positions with SNPs and indels.

```bash
gatk --java-options "-Xmx10g" HaplotypeCaller  -R /data/shared/data/references/human/human_g1k_v37.fasta -I SRR85954XX_sub_markdup.bam -L 1 -O SRR85954XX.vcf.gz --dbsnp /data/shared/groups/group_7/data/databases/00-All.vcf.gz 
```

**Concatenate chromosomal VCFs**

The VCFs should be sorted before combining them.

```bash
bcftools concat -o SRR85954XX.vcf.gz --threads 2 SRR85954XX_chr_1.vcf.gz SRR85954XX_chr_2.vcf.gz SRR85954XX_chr_3.vcf.gz SRR85954XX_chr_4.vcf.gz SRR85954XX_chr_5.vcf.gz SRR85954XX_chr_6.vcf.gz SRR85954XX_chr_7.vcf.gz SRR85954XX_chr_8.vcf.gz SRR85954XX_chr_9.vcf.gz SRR85954XX_chr_10.vcf.gz SRR85954XX_chr_11.vcf.gz SRR85954XX_chr_12.vcf.gz SRR85954XX_chr_13.vcf.gz SRR85954XX_chr_14.vcf.gz SRR85954XX_chr_15.vcf.gz SRR85954XX_chr_16.vcf.gz SRR85954XX_chr_17.vcf.gz SRR85954XX_chr_18.vcf.gz SRR85954XX_chr_19.vcf.gz SRR85954XX_chr_20.vcf.gz SRR85954XX_chr_21.vcf.gz SRR85954XX_chr_22.vcf.gz
```



​	













**Check concordance of variant calling across sample type**

Concordance statistics can be obtained from

```bash
bcftools isec -p concordance SRR8595490.vcf.gz SRR8595494.vcf.gz
```

where `-p concordance` defines the directory to which the output is saved.

The command saves 4 vcf-files :

* `0000.vcf `- records private to vcf

* `10001.vcf` - records private to vcf

* `20002.vcf` - shared records

* `0003.vcf` - shared records

  

  Afterwards relevant numbers can be extracted from the output files using

  ```
  bcftools view -H -v [snps or indels] 0000.vcf | wc -l
  ```

  

**Further unmapped reads analysis**

The following was prepared to do, but not done due to time restrains. The BLAST database `nt`  takes some time to build. Since the goal of our project was to reproduce the paper and we had some time restrains, we did not carry out the unmapped reads analysis. Even though we are aware of the existence of other tools to process the unmapped reads analysis. 

*Get data ready for the BLAST*

Specifically, we selected 0.2% of the read pairs from each sample for which both reads in the pair were unmapped, and used the first read in each pair as a BLAST query. 

```bash
for i in `seq X X `; do echo "nice -n 19 samtools view -b -f12 SRR859549"$i"_merged_markdup.bam |  samtools view -b -s 20.0002 /dev/stdin > SRR859549"$i"_sub.bam"; done | parallel -j 20
```

Default parameters were used to blastn except `tblastn` and `-evalue 1e-5` , and we used the default nucleotide db, `nt`.  First, extract the unmapped reads and transform them into FASTQ format and then feed it to `tblastn`. 

```bash
samtools view --threads 4 --write-index -f SRR859549XX_merged_markdup.bam | awk '{printf(">%s/%s\n%s\n",$1,(and(int($2),0x40)?1:2),$10);}' |  tblastn -db nt -evalue 1e-5 -out SRR859549XX_blast.txt
```

An option to run BLAST if it is not build in the cluster that you are using is using BLAST webserver. However this has some restrains on the input size.

[BLAST online](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome)

For each match, the esummary program from NCBI was used to determine the taxonomic ID of the organism from which the database sequence was derived. The domain (e.g., bacteria or eukaryota) associated with that taxonomic ID was determined using the `fullnamelineage.dmp` file, which can be downloaded from: https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz. 



