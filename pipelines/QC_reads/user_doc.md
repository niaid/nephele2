# Pre-processing QC Pipeline


## Pipeline Steps and Tools

  - You can choose which of the steps in the pipeline you would like to
    run, with the exception of FastQC which is always the first step.
  - If you do not know very much about your data quality, we recommend
    that you run the pipeline without choosing any of the steps. This
    will run only FastQC and give you a summary report of your data
    quality.

### 1\. FastQC

  - FastQC is *always run* first in the pipeline with default
    parameters.
  - FastQC analyzes the input FASTQ files and reports summary statistics
    about each file in both tabular and graphical format, including
    number of reads, average per base quality score, etc.
  - Nephele uses
    [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
    v0.11.7.

### 2\. Adapter/Primer trimming

  - Adapter and/or primer trimming is *optional*.
  - Nephele uses the [cutadapt
    plugin](https://docs.qiime2.org/2018.6/plugins/available/cutadapt/)
    from QIIME 2 version 2018.6 for primer and adapter trimming.  
  - cutadapt trims the sequences specified by the user from either the
    5’ or 3’ ends of reads.
  - To trim primers for amplicon sequence data, the primer sequence
    should be specified as the Forward and/or Reverse front 5’ adapter.
  - More information on adapter and primer trimming can be found in the
    [QIIME 2
    docs](https://docs.qiime2.org/2018.6/plugins/available/cutadapt/),
    and on [cutadapt’s help
    page](https://cutadapt.readthedocs.io/en/stable/).

### 3\. Quality trimming

  - Quality trimming is *optional*.
  - [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) v0.38
    is used for quality trimming.
  - Trimmomatic uses a sliding window from the 5’ to the 3’ end of the
    read. When the average quality in the window drops below the
    required quality score, the read is trimmed to remove the 3’ low
    quality end of the read.
  - Poor quality sequence can also be trimmed at the start of a read
    using the [`Leading quality`](#quality-trimming-1) trim option.
  - Further help can be found in the [Trimmomatic
    manual](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf).

### 4\. Paired-end read merging

  - Read merging is only for paired-end data sets and is *optional*.
  - For merging read pairs, Nephele uses
    [FLASH2](https://github.com/dstreett/FLASH2) v2.2.00 which is based
    on the [FLASH](https://ccb.jhu.edu/software/FLASH/) read merger.
  - FLASH is designed to merge pairs of reads when the original fragment
    length is shorter than twice the length of reads. It merges read
    pairs if the rate of mismatches in the overlapping region is less
    than the user-specified [`Maximum mismatch
    density`](#merge-read-pairs). See user options below for more
    information.

### 5\. Summary graphs

  - Summary graphs are *always made* as the final step of the pipeline.
  - Nephele uses [MultiQC](http://multiqc.info/) v1.8dev.  
  - MultiQC summarizes the output of FastQC, cutadapt, Trimmomatic, and
    FLASh (if any of those steps are chosen by the user) and produces a
    single html report with graphs.

## User Options

### Adapter/Primer trimming

  - At least one of the following options for the adapter/primer
    sequences must be specified:
      - **3’ adapter**: Sequence of an adapter ligated to the 3’ end.
        The adapter and any subsequent bases are trimmed. If a `$` is
        appended, the adapter is only found if it is at the end of the
        read. Maximum 250bp.
      - **Front 5’ adapter**: Sequence of an adapter or primer ligated
        to the 5’ end. The adapter and any preceding bases are trimmed.
        Partial matches at the 5’ end are allowed. If a `^` character is
        prepended, the adapter is only found if it is at the beginning
        of the read. Maximum 250bp.
      - **Anywhere adapter**: Sequence of an adapter that may be ligated
        to the 5’ or 3’ end. Both types of matches as described under
        `3' adapter` and `Front 5' adapter` are allowed. If the first
        base of the read is part of the match, the behavior is as with
        `Front 5'`, otherwise as with `3' adapter`. This option is
        mostly for rescuing failed library preparations - do not use if
        you know which end your adapter was ligated to. Maximum 250bp.
  - **Error rate**: Maximum allowed error rate. (default: 0.1)
  - **Indels**: Allow insertions or deletions of bases when matching
    adapters. (default: True)
  - **Adapter overlap**: Require at least this many bases of overlap
    between read and adapter for an adapter to be found. (default: 3)
  - **Match read wildcards**: Interpret IUPAC wildcards (e.g., N) in
    reads. (default: False)
  - **Match adapter wildcards**: Interpret IUPAC wildcards (e.g., N) in
    adapters. (default: True)

### Quality trimming

  - **Window size**: Sliding window size for quality trimming, cutting
    once the average quality within the window falls below the `required
    quality`. Specifies the number of bases to average across. (default:
    4)
  - **Required quality**: Specifies the average quality required in the
    sliding window. (default: 12)
  - **Leading quality**: Remove low quality bases from the beginning (5’
    end) of the read. As long as a base has a value below this threshold
    the base is removed and the next base will be investigated.
    (default: 0)
  - **Trailing quality**: Remove low quality bases from the end of the
    read. As long as a base has a value below this threshold the base is
    removed and the next base (which as Trimmomatic is starting from the
    3’ end would be base preceding the just removed base) will be
    investigated. (default: 0)
  - **Minimum length**: Remove reads that fall below the specified
    minimal length. (default: 60)
  - **Average quality**: Drop the read if the average quality across the
    entire read is below the specified level. (default: 0)

### Merge Read Pairs

  - **Minimum overlap**: The minimum required overlap length between two
    reads to provide a confident overlap. (default: 10)
  - **Maximum overlap**: Maximum overlap length expected in
    approximately 90% of read pairs. Overlaps longer than the maximum
    overlap parameter are still considered as good overlaps, but the
    `mismatch density` (explained below) is *only calculated over the
    first `max_overlap` bases in the overlapped region* rather than the
    entire overlap. (default: 315)
  - **Minimum outie overlap**: The minimum required overlap length
    between two reads to provide a confident overlap in an “outie”
    scenario. (default: 35)

<img src="../../nephele/static/images/outie.svg" width="30%" />

  - **Maximum mismatch density**: Maximum allowed ratio between the
    number of mismatched base pairs and the overlap length. (default:
    0.25)

## Output Files

  -   - *logfile.txt*  
        log file containing messages from running the pipeline.

  -   - **cutadapt\_trimmed**  
        adapter/primer trimmed FASTQ files output by q2-cutadapt.

  -   - **merged**  
        merged FASTQ files and other output from FLASH2
          - *\_merged.extendedFrags.fastq.gz*: FASTQ file of merged
            reads.
          - *\_merged.notCombined\_1/2.fastq.gz*: FASTQ files of reads
            that were not able to be merged.
          - *\_merged.hist*: numeric histogram of merged read lengths.
          - *\_merged.histogram*: visual histogram of merged read
            lengths.

  -   - **multiqc\_input**  
        contains log files and FastQC output used by MultiQC.
          - *\*fastqc.html/.zip*: FastQC reports. See the [FastQC
            Analysis module
            docs](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/)
            for more information.

  -   - *multiqc\_report.html*  
        [MultiQC HTML
        report](http://multiqc.info/docs/#using-multiqc-reports) with
        summary graphs.

  -   - **multiqc\_data**  
        files containing [parsed
        data](http://multiqc.info/docs/#parsed-data-directory) output by
        MultiQC.

  -   - **qtrimmed\_seqs**  
        quality trimmed FASTQ files output by Trimmomatic. For
        paired-end data, there will be 4 files:
          - *.trimmed\_1/2P.fastq.gz*: paired-end output where both
            reads survived processing.
          - *.trimmed\_1/2U.fastq.gz*: corresponding unpaired output
            where a read survived, but the partner read did not.

## References

1.  Andrews, S. Babraham Bioinformatics - FastQC A Quality Control tool
    for High Throughput Sequence Data. (n.d.). Retrieved September 11,
    2018, from
    <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>.
2.  Martin, M. (2011). Cutadapt removes adapter sequences from
    high-throughput sequencing reads. *EMBnet.Journal*, *17*(1), 10-12.
    doi: [10.14806/ej.17.1.200](https://doi.org/10.14806/ej.17.1.200).
3.  Bolyen, Evan, et al. “Reproducible, Interactive, Scalable and
    Extensible Microbiome Data Science Using QIIME 2.” *Nature
    Biotechnology*, July 2019,
    [doi:\[10.1038/s41587-019-0209-9](doi:%5B10.1038/s41587-019-0209-9)\](<https://doi.org/10.1038/s41587-019-0209-9>).
4.  Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: a
    flexible trimmer for Illumina sequence data. *Bioinformatics*,
    *30*(15), 2114. doi:
    [10.1093/bioinformatics/btu170](https://doi.org/10.1093/bioinformatics/btu170).
5.  Magoc, T., & Salzberg, S. L. (2011). FLASH: fast length adjustment
    of short reads to improve genome assemblies. *Bioinformatics*,
    *27*(21), 2957-2963. doi:
    [10.1093/bioinformatics/btr507](https://doi.org/10.1093/bioinformatics/btr507).
6.  Streett, D. (2018). *Flash2 has some improvements from flash\_1
    including new logic from innie and outie overlaps as well as some
    initial steps for flash for amplicons: dstreett/FLASH2*. C.
    Retrieved from <https://github.com/dstreett/FLASH2> (Original work
    published 2015).
7.  Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC:
    summarize analysis results for multiple tools and samples in a
    single report. *Bioinformatics*, *32*(19), 3047-3048. doi:
    [10.1093/bioinformatics/btw354](https://doi.org/10.1093/bioinformatics/btw354).
