/*  COVID-19 Pipeline
 *  Usage: nextflow run /path/to/main.nf
 *
 *  Author: Mohammed Khalfan < mkhalfan@nyu.edu >
 *  NYU Center for Genetics and System Biology 2020
 */

// Setting some defaults here,
// can be overridden in config or via command line
params.out = "${params.outdir}/out"
params.tmpdir = "${params.outdir}/gatk_temp"
params.snpEff_config = "${params.snpeff}"

// Define modules here

// println "reads: $params.reads"
println "inputs dir: $params.inputs_dir"
println "mapping file: $params.map_file"
println "ref: $params.ref"
println "outdir: $params.out"

ref = file(params.ref)
snpeff_config = file(params.snpEff_config)
//primers = file(params.primers)
//primers_a = file(params.primers_a)
//primers_b = file(params.primers_b)

// Prepare the fastq read pairs for input.
// Use the size parameter to not auto-group, and instead
// use the mapping through getBaseName() and subtract
// two regexs to get the ID.
// This enables support for CGSB sequence data file naming format
// Channel
//     .fromFilePairs( params.reads, size: -1)
//     { file -> file.getBaseName() - ~/${params.read_pair_regex}/ - ~/.fastq/ }
//     .set { read_pairs_ch }

// Reading mapping --> generate pairs of: sample ID, sequences and primer files
map_file = file(params.map_file)
pairs = []
lines  = map_file.readLines()
for( line : lines ) {
    if(line == null || line.isEmpty()) {
        continue
    }
    row = line.split("\t")
    if(row[0] == "#SampleID") {
        continue
    } else {
        pair_A = [row[0], "A", [params.inputs_dir+row[1], params.inputs_dir+row[2]], params.inputs_dir+row[5]]
        pair_B = [row[0], "B", [params.inputs_dir+row[3], params.inputs_dir+row[4]], params.inputs_dir+row[6]]
        pairs.add(pair_A)
        pairs.add(pair_B)
    }
}
// pairs = [
//     ["N87", "A", ["/home/admin/samples/N87-A_S9_L001_R1_001.fastq", "/home/admin/samples/N87-A_S9_L001_R2_001.fastq"], "/usr/local/src/nephele2/pipelines/SARS-CoV2/new_A.fa"],
//     ["N87", "B", ["/home/admin/samples/N87-B_S10_L001_R1_001.fastq", "/home/admin/samples/N87-B_S10_L001_R2_001.fastq"], "/usr/local/src/nephele2/pipelines/SARS-CoV2/new_B.fa"]
// ]
// pairs = [
//     ["N85", "A", ["/home/admin/samples/N85-A_S7_L001_R1_001.fastq", "/home/admin/samples/N85-A_S7_L001_R2_001.fastq"], "/usr/local/src/nephele2/pipelines/SARS-CoV2/new_A.fa"],
//     ["N85", "B", ["/home/admin/samples/N85-B_S8_L001_R1_001.fastq", "/home/admin/samples/N85-B_S8_L001_R2_001.fastq"], "/usr/local/src/nephele2/pipelines/SARS-CoV2/new_B.fa"]
// ]
println "pairs: $pairs"

// Channel
//     .fromList(pairs)
//     .set { read_pairs_ch }

Channel
    .fromFilePairs( params.bams, size: 1)
    { file -> file.getBaseName() - ~/.bam/ }
    .set { bams_in_ch }

process trim {
    publishDir "${params.out}/trimmed", mode:'copy'
    errorStrategy 'finish'

    input:
    // set pair_id, reads, primer from pairs
    set sample_id, pool, reads, primer from pairs


    output:
    set val(sample_id),
    val(pool),
	file("${sample_id}_${pool}_trimmed_1.fq.gz"),
	file("${sample_id}_${pool}_trimmed_2.fq.gz") \
	into trimmed_ch

    script:
    // Set the A or B primer file according to the sample
    // trim_primer_cmd = null
    // if(pair_id.contains("-A_") || pair_id.contains("_A_")){
	// trim_primer_cmd = "ILLUMINACLIP:${primers_a}:2:30:10:8:true"
    // }
    // else if(pair_id.contains("-B_") || pair_id.contains("_B_")){
    // 	trim_primer_cmd = "ILLUMINACLIP:${primers_b}:2:30:10:8:true"
    // }
    // else{
	// trim_primer_cmd = ""
    // }
    trim_primer_cmd = "ILLUMINACLIP:${file(primer)}:2:30:10:8:true"
    """
    java -jar ${params.trimjar} \
	PE \
	-phred33 \
	-threads ${task.cpus} \
	${reads[0]} \
	${reads[1]} \
	${sample_id}_${pool}_trimmed_1.fq.gz \
	${sample_id}_${pool}.unpair_trimmed_1.fq.gz \
	${sample_id}_${pool}_trimmed_2.fq.gz \
	${sample_id}_${pool}.unpair_trimmed_2.fq.gz \
	ILLUMINACLIP:${params.adapters}:2:30:10:8:true \
	${trim_primer_cmd} \
	LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:20
    """
}

process align {
    publishDir "${params.out}/aligned_reads", mode:'copy'
    errorStrategy 'finish'

    input:
    set sample_id,
    pool,
	file(read_1),
	file(read_2) from trimmed_ch

    output:
    val(pair_id_unmerged) into jbrowse_pair_id_ch
    set val(sample_id),
	file("${sample_id}_${pool}_aligned_reads.bam") \
	into aligned_reads_ch
    set val(pair_id_unmerged),
        file("${sample_id}_${pool}_aligned_reads.bam"),
	file("${sample_id}_${pool}_aligned_reads.bai") \
	into individual_bw_ch

    script:
    // sample_id=pair_id - ~/${params.grouping_regex}/
    pair_id_unmerged=sample_id + "_${pool}_unmerged"
    readGroup = "@RG\\tID:${sample_id}_${pool}\\tLB:${sample_id}_${pool}\\tPL:${params.pl}\\tPM:${params.pm}\\tSM:${sample_id}_${pool}"
    """
    bwa mem \
	-K 100000000 \
	-v 3 -t ${task.cpus} \
	-Y \
	-R \"${readGroup}\" \
	$ref \
	$read_1 \
	$read_2 \
	> ${sample_id}_${pool}_aligned_reads.sam

    java -jar ${params.picardjar} SortSam \
	I=${sample_id}_${pool}_aligned_reads.sam \
	O=${sample_id}_${pool}_aligned_reads.bam \
	SORT_ORDER=coordinate \
	CREATE_INDEX=true
    """
}

process mergeBam{
    publishDir "${params.out}/merged", mode:'copy'
    errorStrategy 'finish'

    input:
    set sample_id, file(sams) \
	from aligned_reads_ch
	.groupTuple(size: 2, remainder: true)

    output:
    set val(sample_id),
	file("{*fixed.bam,*unmerged.bam}") \
	into merged_bam_ch

    script:
    if( sams.size() == 2 )
    """
    java -jar ${params.picardjar} MergeSamFiles \
	I=${sams[0]} \
	I=${sams[1]} \
	O=${sample_id}_merged.bam
    java -jar ${params.picardjar} AddOrReplaceReadGroups \
	I=${sample_id}_merged.bam \
        O=${sample_id}_fixed.bam \
        RGID=${sample_id} \
        RGLB=${sample_id} \
        RGPL=${params.pl} \
        RGPU=${params.fcid} \
        RGSM=${sample_id}
    """

    else
    """
    # need to 'output' the file so it goes in the channel
    mv ${sams[0]} ${sample_id}_unmerged.bam
    """

}

process markDuplicatesSpark  {
    publishDir "${params.out}/sorted", mode:'copy'
    errorStrategy 'finish'

    input:
    set val(sample_id),
	file(bam) from merged_bam_ch
	.mix(bams_in_ch)

    output:
    val(sample_id) into jbrowse_sample_id_ch
    set val(sample_id),
	file("${sample_id}_sorted_dedup.bam") \
	into sorted_dedup_bam_ch, sorted_dedup_ch_for_metrics, downsample_bam_ch, pilon_ch, bcftools_ch, consensus_bam_ch
    set val(sample_id),
        file("${sample_id}_sorted_dedup.bam"),
	file("${sample_id}_sorted_dedup.bam.bai") \
	into merged_bw_ch
    set val(sample_id),
	file("${sample_id}_dedup_metrics.txt") into dedup_qc_ch

    script:
    """
    java -jar ${params.gatkjar} \
	 MarkDuplicates \
	-I ${bam} \
	-M ${sample_id}_dedup_metrics.txt \
	-O ${sample_id}_sorted_dedup.bam \
	--CREATE_INDEX true
    mv ${sample_id}_sorted_dedup.bai ${sample_id}_sorted_dedup.bam.bai
    """
}

process getMetrics {
    publishDir "${params.out}/metrics", mode:'copy'
    errorStrategy 'finish'

    input:
    set val(sample_id),
	file(sorted_dedup_reads) from sorted_dedup_ch_for_metrics

    output:
    set val(sample_id),
            file("${sample_id}_alignment_metrics.txt"),
            file("${sample_id}_insert_metrics.txt"),
            file("${sample_id}_insert_size_histogram.pdf"),
            file("${sample_id}_depth_out.txt") \
            into metrics_output

    script:
    """
    java -jar ${params.picardjar} \
        CollectAlignmentSummaryMetrics \
	R=${params.ref} \
        I=${sorted_dedup_reads} \
	O=${sample_id}_alignment_metrics.txt
    java -jar ${params.picardjar} \
        CollectInsertSizeMetrics \
        INPUT=${sorted_dedup_reads} \
	OUTPUT=${sample_id}_insert_metrics.txt \
        HISTOGRAM_FILE=${sample_id}_insert_size_histogram.pdf
    samtools depth -a ${sorted_dedup_reads} > ${sample_id}_depth_out.txt
    """
}

process pilon{
    publishDir "${params.out}/pilon", mode:'copy'
    errorStrategy 'finish'

    input:
    set val(sample_id),
	file(preprocessed_bam) from pilon_ch

    output:
    file("${sample_id}_pilon.vcf") into pilon_bzip_tabix_vcf_ch
    file '*' into pilon_out_ch

    script:
    """
    java -Xmx16G -jar ${params.pilon} \
	--genome $ref \
	--bam $preprocessed_bam \
	--fix bases \
	--changes \
	--vcf \
	--threads ${task.cpus} \
	--mindepth 10 \
	--output ${sample_id}_pilon_g

    java -jar ${params.gatkjar} SelectVariants \
	-V ${sample_id}_pilon_g.vcf \
	-O ${sample_id}_pilon.vcf \
	--exclude-non-variants \
	--exclude-filtered
    """
}

process bcftools{
    publishDir "${params.out}/bcftools", mode:'copy'
    errorStrategy 'finish'

    input:
    set val(sample_id),
        file(preprocessed_bam) from bcftools_ch

    output:
    file("${sample_id}_bcftools.vcf") into bcftools_bzip_tabix_vcf_ch

    script:
    """
    bcftools mpileup \
	--redo-BAQ \
	--adjust-MQ 50 \
	--gap-frac 0.05 \
	--max-depth 10000 \
	--max-idepth 200000 \
	--fasta-ref $ref \
	$preprocessed_bam | bcftools call \
	--ploidy 1 \
	--keep-alts \
	--multiallelic-caller \
	--variants-only \
	--output ${sample_id}_bcftools.vcf
    """
}

process haplotypeCaller {
    errorStrategy 'finish'
    input:
    set val(sample_id),
	file(preprocessed_bam) from sorted_dedup_bam_ch

    output:
    set val(sample_id),
	file("${sample_id}_raw_variants.vcf") into hc_output_ch
    set val(hc_bamout_sample_id),
	file("${sample_id}_haplotypecaller_bamout.bam"),
	file("${sample_id}_haplotypecaller_bamout.bai") \
	into hc_bam_bw_ch

    script:
    hc_bamout_sample_id = sample_id + "-hc_bamout"
    """
    java -jar ${params.gatkjar} HaplotypeCaller \
	-R $ref \
	-I $preprocessed_bam \
	-O ${sample_id}_raw_variants.vcf \
	-bamout ${sample_id}_haplotypecaller_bamout.bam \
	-ploidy 1
    """
}

process selectVariants {
    errorStrategy 'finish'
    input:
    set val(sample_id),
	file(raw_variants) from hc_output_ch

    output:
    set val(sample_id),
	file("${sample_id}_raw_snps.vcf") \
	into raw_snps_ch, raw_snps_qc_ch
    set val(sample_id),
	file("${sample_id}_raw_indels.vcf") into raw_indels_ch

    script:
    """
    java -jar ${params.gatkjar} SelectVariants \
	-R $ref \
	-V $raw_variants \
	-select-type SNP \
	-O ${sample_id}_raw_snps.vcf
    java -jar ${params.gatkjar} SelectVariants \
        -R $ref \
        -V $raw_variants \
        -select-type INDEL \
        -O ${sample_id}_raw_indels.vcf
    """
}

process filterSnps {
    publishDir "${params.out}/filtered_snps", mode:'copy'
    errorStrategy 'finish'

    input:
    set val(sample_id),
	file(raw_snps) from raw_snps_ch

    output:
    set val(sample_id),
        file("${sample_id}_filtered_snps.vcf") \
        into filtered_snps_qc_ch
    set val(sample_id),
	file("${sample_id}_filtered_snps_eaf.vcf") \
	into snpeff_ch
    set val(sample_id),
        file("${sample_id}_consensus_snps.vcf") \
        into consensus_snps_ch
    file("${sample_id}_consensus_snps.vcf") \
	into cons_bzip_tabix_vcf_ch

    script:
    """
    java -jar ${params.gatkjar} VariantFiltration \
	-R $ref \
	-V $raw_snps \
	-O ${sample_id}_filtered_snps.vcf \
        -filter-name "QD_filter" -filter "QD < 2.0" \
        -filter-name "FS_filter" -filter "FS > 60.0" \
        -filter-name "MQ_filter" -filter "MQ < 40.0" \
        -filter-name "SOR_filter" -filter "SOR > 4.0" \
        -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"

    # This script generates the _consensus_snps.vcf
    # and _eaf.vcf (empirical AF) files
    python3 ${params.filter_variants_script} ${sample_id}
    """
}

process filterIndels {
    publishDir "${params.out}/filtered_indels", mode:'copy'
    errorStrategy 'finish'
    input:
    set val(sample_id),
	file(raw_indels) from raw_indels_ch

    output:
    file("${sample_id}_filtered_indels.vcf") into indel_bzip_tabix_vcf_ch

    script:
    """
    java -jar ${params.gatkjar} VariantFiltration \
        -R $ref \
        -V $raw_indels \
        -O ${sample_id}_filtered_indels.vcf \
	-filter-name "DP_filter" -filter "DP < 20.0" \
	-filter-name "QD_filter" -filter "QD < 2.0" \
	-filter-name "FS_filter" -filter "FS > 200.0" \
	-filter-name "SOR_filter" -filter "SOR > 10.0"
    """
}

process consensus {
    publishDir "${params.out}/consensus", mode:'copy'
    errorStrategy 'finish'

    input:
    set val(sample_id),
	file(filtered_snps),
	file(bam) \
	from consensus_snps_ch
	.join(consensus_bam_ch)

    output:
    file("${sample_id}*.fasta") into consensus_ch

    script:
    """
    java -jar ${params.gatkjar} IndexFeatureFile \
	-I $filtered_snps
    java -jar ${params.gatkjar} FastaAlternateReferenceMaker \
	-R $ref \
	-O ${sample_id}.fasta \
	-V $filtered_snps

    # chromosome ID needs to match ID in bam for bedtools (maskfasta)
    sed -i 's/1 NC_045512.2:1-29903/NC_045512.2/g' ${sample_id}.fasta
    for x in {1,6,10,20}
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
		$bam | \
		bedtools genomecov \
		-ibam stdin \
		-bga | \
		awk -v threshold="\$x" '\$4<threshold' | \
		awk '{print \$1 "\t" \$2 "\t" \$3}' \
		> ${sample_id}_below_\${x}_cov.bed

	# mask all regions in bedfile produced above
	bedtools maskfasta \
		-fi ${sample_id}.fasta \
		-bed ${sample_id}_below_\${x}_cov.bed \
		-fo ${sample_id}_below_\${x}_masked.fasta

	# rename the fasta header from ref name to sample id
	sed -i 's/NC_045512.2/${sample_id}/g' ${sample_id}_below_\${x}_masked.fasta
    done

    rm -f ${sample_id}.fasta
    """
}

process snpEff{
    publishDir "${params.out}/snpEff", mode:'copy'
    errorStrategy 'finish'

    input:
    set val(sample_id),
	file(snps) \
	from snpeff_ch

    output:
    file '*' into snpeff_out
    file("${sample_id}_filtered_snps.ann.vcf") into snpeff_bzip_tabix_vcf_ch

    script:
    """
    java -jar ${params.snpeffjar} -v \
        -c $snpeff_config \
        SARS-CoV2_NC_045512.2 \
        $snps > ${sample_id}_filtered_snps.ann.vcf
    """
}

process make_bw{
    publishDir "${params.out}/bigwig", mode:'copy'
    errorStrategy 'finish'

    input:
    /* id can be sample_id or pair_id */
    set val(id),
	file(bam),
	file(bam_index) \
	from individual_bw_ch
	.mix(merged_bw_ch)
	.mix(hc_bam_bw_ch)

    output:
    file("${id}_coverage.bam.bw") into jbrowse_bw_ch

    when:
    id != "CV-40-hc_bamout" && id != "CV-70-hc_bamout" && id != "CV-62-hc_bamout"

    script:
    """
    bamCoverage \
        -p max  \
        --bam $bam \
	--binSize 1 \
	--ignoreDuplicates \
	--minMappingQuality 20 \
        -o ${id}_coverage.bam.bw
    """
}

process downsample_bam{
    errorStrategy 'finish'
    input:
    set val(sample_id), file(bam) from downsample_bam_ch

    output:
    set file("${sample_id}_downsampled.bam"),
        file("${sample_id}_downsampled.bam.bai") into jbrowse_bam_ch

    script:
    """
    java -jar ${params.jvarkitsortjar} \
        --samoutputformat BAM \
        $bam |\
        java -jar ${params.jvarkitbiojar} \
        -n 75 \
        --samoutputformat BAM |\
        samtools sort -o ${sample_id}_downsampled.bam
    samtools index ${sample_id}_downsampled.bam
    """
}

process bzip_tabix_vcf{
    errorStrategy 'finish'
    input:
    file(vcf) from pilon_bzip_tabix_vcf_ch
	.mix(cons_bzip_tabix_vcf_ch)
	.mix(indel_bzip_tabix_vcf_ch)
	.mix(snpeff_bzip_tabix_vcf_ch)
	.mix(bcftools_bzip_tabix_vcf_ch)

    output:
    file("*.vcf.gz*") into jbrowse_vcf_ch

    script:
    """
    bgzip -c ${vcf} > ${vcf}.gz
    tabix -p vcf ${vcf}.gz
    """
}

process jbrowse{
    publishDir "${params.out}/trackList", mode:'copy'
    errorStrategy 'finish'

    input:
    val pair_ids from jbrowse_pair_id_ch.collect().ifEmpty("")
    val sample_ids from jbrowse_sample_id_ch.collect()
    file '*' from jbrowse_bw_ch.collect()
    file '*' from jbrowse_bam_ch.collect()
    file '*' from jbrowse_vcf_ch.collect()

    output:
    file '*.json' into trackList_ch

    when:
    false

    script:
    """
    copy-data-to-jbrowse.sh "${pair_ids}" "${sample_ids}" $params.fcid "$params.grouping_regex"
    """
}

process qc {
    errorStrategy 'finish'
    input:
    set val(sample_id),
	file("${sample_id}_alignment_metrics.txt"),
	file("${sample_id}_insert_metrics.txt"),
	file("${sample_id}_insert_size_histogram.pdf"),
	file("${sample_id}_depth_out.txt"),
	file("${sample_id}_dedup_metrics.txt"),
	file("${sample_id}_raw_snps.vcf"),
        file("${sample_id}_filtered_snps.vcf") \
	from metrics_output
	.join(dedup_qc_ch)
	.join(raw_snps_qc_ch)
	.join(filtered_snps_qc_ch)

    output:
    file("${sample_id}_report.csv") into parse_metrics_output

    script:
    """
    /bin/bash ${params.parse_metrics_script} ${sample_id} > ${sample_id}_report.csv
    """
}

/* Process qc above creates a report for each sample.
 * Below we compile these into a single report.
 */
parse_metrics_output.collectFile(name: "${workflow.runName}_report.csv", keepHeader: true, storeDir: "${params.out}/reports")
