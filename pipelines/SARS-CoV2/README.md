# SARS-CoV2 pipeline README #
This pipeline was written with nextflow to trim, align and call variants from SARS-CoV-2 sequence output. The pipeline was specifically designed to work with data generated from the [ARTIC consortium protocol](https://artic.network/ncov-2019) and adapted versions. We specifically recommend the protocol developed by the Ghedin Lab which is available upon request.

**Amplicons should be generated in two pools, A and B, from separate preparations from amplification through library generation and sequencing.** The pipeline takes single or paired-end fastq files and trims both adapters and amplification primers and aligns them to the provided reference file, then merges the A and B files to created merged BAM files. The pipeline additionally calls minority variants using haplotypeCaller from GATK and has an option for uploading the sequences to jbrowse for visualization. The pipeline will work on any number of samples, provided that each sample has 2 (single end) or 4 (paired end) input fastq files (read files from both A and B amplicon preparations).

Files in fasta format should be generated separately for A and B primers to be referenced by the pipeline for primer trimming from sequence data. We use the Wuhan-Hu-1 sequence as the reference sequence for SARS-CoV-2 data (GenBank: MN908947.3).

![Alt text](https://github.niaid.nih.gov/rapleeid/nephele2/blob/next_release/pipelines/SARS-CoV2/images/SARSCOV2Artic.png?raw=true "Pipeline")
Solid lines represent the core functions of the pipeline.
Dashed lines represent additionally functionality for viewing on jbrowse. (*coming soon*)

# Example Usage #
SARS-CoV2 pipeline is for Variant calling of SARS-CoV2 sequences.

## Unique sample pairing requirements ##
Input sample files need to be matched to their respective primer pair. Sample1+n should have two input files for SE and 4 input files for PE.
SE Run: Sample1_primerA.fastq and Sample1_primerB.fastq
PE Run: Sample1_primerA_Forward.fastq Sample1_primerA_Reverse.fastq Sample1_primerB_Forward.fastq Sample1_primerB_Reverse.fastq.

Each sample needs an option to upload Single End sampleN primerA and sampleN primerB or Paired End sampleN_forward primerA, sampleN_reverse primerA, sampleN_forward primerB, sampleN_reverse primerB

## Unique sample naming requirements ##
Currently this pipeline requires samples to be named in the following format:

*_L001_R{1,2}_001.fastq.gz

The * in the sample name needs to be your sample identifier and primer identifier. Example N85-A.

The full sample name would like N85-A_S7_L001_R1_001.fastq.gz.
In this example, QC and alignment outputs are identified by N85-A until A and B portions of the sample are merged.
Further downstream outputs would be identified by N85.
## Dependencies ##

- [TRIMMOMATIC](http://www.usadellab.org/cms/?page=trimmomatic)
- [BWA 0.7.17](http://bio-bwa.sourceforge.net/bwa.shtml)
- [PICARD 2.17.11](https://broadinstitute.github.io/picard/)
- [GATK 4.1.3.0](https://gatk.broadinstitute.org/hc/en-us)
- [SAMTOOLS 1.9](http://www.htslib.org)
- [HTSLIB 1.4.1](http://www.htslib.org)
- [DEEPTOOLS 3.4.1](https://deeptools.readthedocs.io/en/develop/)
- [JVARKIT](https://lindenb.github.io/jvarkit/)
- [PILON 1.23](https://github.com/broadinstitute/pilon)
- [BCFTOOLS 1.9](http://www.htslib.org/download/)
- [BEDTOOLS 2.27.1](https://bedtools.readthedocs.io/en/latest/)
- [R 3.4.3](https://www.r-project.org)
- [PYSAM](https://pysam.readthedocs.io/en/latest/)
- [PYPAIRIX](https://pypi.org/project/pypairix/)
- ref = "SARS-CoV2.fa"

## Files/Directories Necessary for Nextflow to run on AWS ##
- /samples/
- /refs/
- /snpeff_data/
- main.nf
- nextflow.config
- filter_variants.py
- parse_metrics.sh
- snpEff.config
- update_trackList.py
- copy-data-to-jbrowse.sh

**Dependent files for testing pipeline on Nephele with sample data from the Ghedin Lab**

*Primer files for TRIMMOMATIC*
- new_A.fa
- new_B.fa

*Sample files for the entire run*
- N85-A_S7_L001_R1_001.fastq.gz
- N85-A_S7_L001_R2_001.fastq.gz
- N85-B_S8_L001_R1_001.fastq.gz
- N85-B_S8_L001_R2_001.fastq.gz

## Trimmomatic PE primerA ##
The following code represents a paired end run with primerA
```
java -jar /directory/to/trimmomatic.jar PE -phred33 -threads <threads> N85-A_S7_L001_R1_001.fastq.gz N85-A_S7_L001_R2_001.fastq.gz N85-A_trimmed_1.fq.gz N85-A_unpaired_trimmed_1.fq.gz N85-A_trimmed_2.fq.gz N85-A_unpaired_trimmed_2.fq.gz ILLUMINACLIP:/directory/to/adapterfile.fasta:2:30:10:8:true ILLUMINACLIP:/directory/to/primerA.fasta:2:30:10:8:true LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:20
```
## Trimmomatic PE primerB ##
The following code represents a paired end run with primerB
```
java -jar /directory/to/trimmomatic.jar PE -phred33 -threads <threads> N85-B_S7_L001_R1_001.fastq.gz N85-B_S7_L001_R2_001.fastq.gz N85-B_trimmed_1.fq.gz N85-B_unpaired_trimmed_1.fq.gz N85-B_trimmed_2.fq.gz N85-B_unpaired_trimmed_2.fq.gz ILLUMINACLIP:/directory/to/adapterfile.fasta:2:30:10:8:true ILLUMINACLIP:/directory/to/primerB.fasta:2:30:10:8:true LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:20
```

### Trimmomatic parameters ###
- ``SE``: Single end file input
- ``PE``: Paired end file inputs
- `-phred33``-phred64`: Specifies the base quality encoding. If no quality encoding is specified it will be determined automatically (>=version 0.32).
- `-threads`: Number of threads to use.
- `-trimlog`: Specifying a trimlog file which logs read trimmings, indicating read name, surviving sequence length, location of first surviving base, and amount trimmed from end.
- ``ILLUMINACLIP:<fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clip threshold>:<simple clip threshold><minAdapterLength><keepBothReads>``: Find and remove adapters: file location of adapters: a number which specifies the maximum mismatch count which will still allow a full match to be performed: a number which specifies how accurate the match between the two'adapter ligated' reads must be for PE palindrome read alignment: a number which specifies how accurate the match between any adapter sequence must be against a read: a number which is the minimum length of the adapter to be detected: True or False to keep both reads.
- ``LEADING:``: Cut bases off the start of the read, if below a threshold quality
- ``TRAILING:``: Cut bases off the end of a read, if below a threshold quality
- ``SLIDINGWINDOW:<windowSize>:<requiredQuality>``: Perform a sliding window trimming, cutting once the average quality within the window falls below a threshold.
- ``MINLEN:``: Minimal read length.

## BWA Alignment ##
Example of BWA alignment with sample N85 ran with primer set A.
```
bwa mem -K 100000000 -v 3 -t ${task.cpus} -Y -R \RG\\tID:N85-A\\tLB:N85-A\\tPL:illumina\\tPM:miseq\\tSM:N85-A $ref N85-A_trimmed_1.fq.gz N85-A_trimmed_2.fq.gz > N85-A_aligned_reads.sam
```

## PICARD SortSam ##
```
java -jar ${params.picardjar} SortSam \
I=N85-A_aligned_reads.sam \
O=N85-A_aligned_reads.bam \
SORT_ORDER=coordinate \
CREATE_INDEX=true
```

## PICARD MergeSamFiles ##
```
java -jar ${params.picardjar} MergeSamFiles \
I=N85-A_aligned_reads.bam \
I=N85-B_aligned_reads.bam \
O=N85_merged.bam
```

## PICARD AddOrReplaceReadGroups ##
```
java -jar ${params.picardjar} AddOrReplaceReadGroups \
I=N85_merged.bam \
O=N85_fixed.bam \
RGID=N85 \
RGLB=N85 \
RGPL=illumina \
RGPU=000000000-CP3RN \
RGSM=N85
```

## GATK MarkDuplicates ##
```
java -jar ${params.gatkjar} \
MarkDuplicates \
-I N85_fixed.bam \
-M N85_dedup_metrics.txt \
-O N85_sorted_dedup.bam \
--CREATE_INDEX true
mv N85_sorted_dedup.bai N85_sorted_dedup.bam.bai
```

## PICARD Metrics ##
```
java -jar ${params.picardjar} \
    CollectAlignmentSummaryMetrics \
R=${params.ref} \
I=N85_sorted_dedup.bam \
O=N85_alignment_metrics.txt
java -jar ${params.picardjar} \
    CollectInsertSizeMetrics \
INPUT=N85_sorted_dedup.bam \
OUTPUT=N85_insert_metrics.txt \
HISTOGRAM_FILE=N85_insert_size_histogram.pdf
samtools depth -a N85_sorted_dedup.bam > N85_depth_out.txt
```

## GATK HaplotypeCaller ##
```
java -jar ${params.gatkjar} HaplotypeCaller \
-R $ref \
-I N85_sorted_dedup.bam \
-O N85_raw_variants.vcf \
-bamout N85_haplotypecaller_bamout.bam \
-ploidy 1
```

## GATK SelectVariants SNPS and Indels ##
```
java -jar ${params.gatkjar} SelectVariants \
-R $ref \
-V N85_raw_variants.vcf \
-select-type SNP \
-O N85_raw_snps.vcf
java -jar ${params.gatkjar} SelectVariants \
    -R $ref \
    -V N85_raw_variants.vcf \
    -select-type INDEL \
    -O N85_raw_indels.vcf
```

## GATK VariantFiltration SNPS ##
```  java -jar ${params.gatkjar} VariantFiltration \
-R $ref \
-V N85_raw_snps.vcf \
-O N85_filtered_snps.vcf \
      -filter-name "QD_filter" -filter "QD < 2.0" \
      -filter-name "FS_filter" -filter "FS > 60.0" \
      -filter-name "MQ_filter" -filter "MQ < 40.0" \
      -filter-name "SOR_filter" -filter "SOR > 4.0" \
      -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"
```

## Filter Variants for consensus ##
Run filter_variants.py after VariantFiltration on SNPS
```
filter_variants.py N85
```

## GATK VariantFiltration Indels ##
```
java -jar ${params.gatkjar} VariantFiltration \
    -R $ref \
    -V N85_raw_indels.vcf \
    -O N85_filtered_indels.vcf \
-filter-name "DP_filter" -filter "DP < 20.0" \
-filter-name "QD_filter" -filter "QD < 2.0" \
-filter-name "FS_filter" -filter "FS > 200.0" \
-filter-name "SOR_filter" -filter "SOR > 10.0"
```

## Consensus Sequence Generation ##
```
java -jar ${params.gatkjar} IndexFeatureFile \
-I N85_filtered_snps.vcf
java -jar ${params.gatkjar} FastaAlternateReferenceMaker \
-R $ref \
-O N85.fasta \
-V N85_filtered_snps.vcf

# chromosome ID needs to match ID in bam for bedtools (maskfasta)
sed -i 's/1 SARS-CoV2:1-29903/SARS-CoV2/g' N85.fasta
for x in {6,10,20}
do
# make bedfile with regions below x coverage
    # genomecov generates bedgraph file
# genomecov input is filtered for min MAPQ (20)
# and to remove dups and non-primary alignments
# first awk filters bedgraph for coverage <= x
# second awk converts bedgraph to 3-col bedfile
samtools view \
-bq 20 \
-F 1284 \
N85_sorted_dedup.bam | \
bedtools genomecov \
-ibam stdin \
-bga | \
awk -v threshold="\$x" '\$4<threshold' | \
awk '{print \$1 "\t" \$2 "\t" \$3}' \
> N85_below_\${x}_cov.bed

# mask all regions in bedfile produced above
bedtools maskfasta \
-fi N85.fasta \
-bed N85_below_\${x}_cov.bed \
-fo N85_below_\${x}_masked.fasta

# rename the fasta header from ref name to sample id
sed -i 's/SARS-CoV2/N85/g' N85_below_\${x}_masked.fasta
done
```

## snpEFF ##
```
java -jar ${params.snpeffjar} -v \
    -c snpeff_config \
    SARS-CoV2_NC_045512.2 \
    N85_filtered_snps.vcf > N85_filtered_snps.ann.vcf
```

## QC ##
This script uses the follwing inputs to generate a QC report
- N85_alignment_metrics.txt
- N85_insert_metrics.txt
- N85_insert_size_histogram.pdf
- N85_depth_out.txt
- N85_dedup_metrics.txt
- N85_raw_snps.vcf
- N85_filtered_snps.vcf
```
parse_metrics.sh N85 > N85_report.csv
```
