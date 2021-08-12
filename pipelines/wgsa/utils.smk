import csv
import json
import pandas as pd
from functools import reduce
import biom
from biom.cli.util import write_biom_table



def read_mapping_file():
    samples = {}
    with open(config["map_file"]) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if row['#SampleID'] == '':
                continue
            samples[row['#SampleID']] = {}
            if row.get('ForwardFastqFile'):
                samples[row['#SampleID']]["f"] = config["input_dir"] + row["ForwardFastqFile"]
            if row.get('ReverseFastqFile'):
                samples[row['#SampleID']]["r"] = config["input_dir"] + row["ReverseFastqFile"]
    return samples

def get_fastq_gz(wildcards):
    return samples[wildcards.sample]

def get_num_threads_per_job(num_samples):
    num_cpus = config["num_cpus"]
    default_threads_per_job = 16
    return min(max((num_cpus / num_samples), default_threads_per_job), num_cpus)

def run_fastp_cmd():
    base_cmd = "fastp -i {input.f} -I {input.r} -o {output.f} -O {output.r} -h {output.s_l_html} -j {output.s_l_json} -y -c --trim_poly_x -w {threads}"
    if config["trim_filter"]:
        trim_filter = "-W 4 -e {config[average_read_quality]} -l {config[min_read_length]} -5 --cut_front_mean_quality {config[trim_of_5]} -3 --cut_tail_mean_quality {config[trim_of_3]}"
        return base_cmd + " " + trim_filter
    else:
        return base_cmd

def get_read_length(fastp_json):
    """
    Given output json file of fastp, parse and return
    average R1 read length after filtering.
    """
    with open(fastp_json, 'r') as fastpfile:
        fastp = json.load(fastpfile)
    return(fastp["summary"]["after_filtering"]["read1_mean_length"])



def make_biom_file(profiles, taxfiles, samples, map_file, biom_file_name, tax_col, prof_col):
    """
    Make biom file from checkm profile and checkm tree_qa taxonomy across
    all samples.
    
    Args:
        profiles (str): profile filename template with ``{sample}`` substring to format.
        taxfiles (str): bin_TAX filename template with ``{sample}`` substring to format.
        samples (list): list of sample names
        map_file (str): mapping filename
        biom_file_name (str): output biom filename
        tax_col (str): column from bin_TAX to use
        prof_col (str): column from profiles to use

    Returns:
        nothing, but produces biom file named *biom_file_name*.
    """
    def merge_single_profile_tax(profile, taxfile, sample, tax_col=tax_col, prof_col=prof_col):
        """
        Merge profile.txt and bin_TAX.txt file for single sample, and sum ``prof_col`` by ``tax_col``.
        Returns `pandas.DataFrame` with only 2 columns with tax_col column name replaced
        by sample name.
        """
        prof = pd.read_csv(profile, sep='\t')
        tax = pd.read_csv(taxfile, sep='\t')
        tax = tax.astype({tax_col: str})
        tab = pd.merge(prof, tax, how="left", on="Bin Id")
        tab = tab[[tax_col, prof_col]]
        tab = tab.groupby(tax_col).sum()
        tab.rename(columns={prof_col : sample}, inplace=True)
        return(tab)

    ## merge tables over all samples; format string to ensure files match
    tablist = [merge_single_profile_tax(profiles.format(sample=s), taxfiles.format(sample=s), s) for s in samples]
    if len(samples) > 1:
        mergedtable = reduce(lambda l,r: pd.merge(l,r, how="outer", on=tax_col), tablist)
    else:
        mergedtable = tablist[0]
    mergedtable = mergedtable.rename(index={'nan':'unclassified'}).fillna(0)

    ## make biom https://biom-format.org/documentation/table_objects.html
    observ_ids = [f"O{i}" for i in range(mergedtable.shape[0])]
    observ_metadata = [{ 'taxonomy' : x.split(';')} for x in mergedtable.index]
    biomtab = biom.table.Table(data=mergedtable.to_numpy(),
                               observation_ids=observ_ids,
                               sample_ids=list(mergedtable.columns),
                               observation_metadata=observ_metadata,
                               type="Taxon table")
    
    ## read in mapping file
    sample_met = biom.parse.MetadataMap.from_file(map_file)
    biomtab.add_metadata(sample_met, axis='sample')
    write_biom_table(biomtab, fmt="hdf5", filepath=biom_file_name)

def clean_output():
    asm_bam=" ".join(expand(ASM + "*.bam", sample=samples))
    asm_bai=" ".join(expand(ASM + "*.bai", sample=samples))
    asm_bt2=" ".join(expand(ASM + "*.bt2", sample=samples))
    asm_fastg=" ".join(expand(ASM + "*.fastg", sample=samples))
    asm_spades_log=" ".join(expand(ASM + "spades.log", sample=samples))
    annotations_log=" ".join(expand(ANNOTATION_PATH + ".log", sample=samples))
    binqc_storage=" ".join(expand(BINQC + "storage", sample=samples))
    qc_plots_done=" ".join(expand(QC_PLOTS + "*.done", sample=samples))

    command = f"""
    set +o pipefail
    echo "Cleaning output..."
    rm -rf {asm_bam} {asm_bai} {asm_bt2} {asm_fastg}\
    {asm_spades_log} {annotations_log} {binqc_storage} {qc_plots_done}
    for f in $(find {METASPADES} -type f -name hmmer*); do
    rm -f $f
    done
    """
    return command

def handle_TED():
    if config["include_ted_files"]:
        command = f"""
        if [ -d {TED_TMP} ] && [ "$(ls -A {TED_TMP})" ]
        then
        tar -cf {TED_FILE_NAME} -C {config["output_dir"]} {TED_TMP_NAME}
        rm -rf {TED_TMP}
        fi
        """
    else:
        command = f"""
        rm -rf {TED_TMP}
        """
    return command
