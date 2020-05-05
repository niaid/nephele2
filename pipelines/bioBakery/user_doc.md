bioBakery WGS Pipeline
================

-   [bioBakery Workflows](#biobakery-workflows)
-   [User Options](#user-options)
-   [Output Files](#output-files)
-   [Tools and References](#tools-and-references)

bioBakery Workflows
-------------------

-   [bioBakery](https://bitbucket.org/biobakery/biobakery/wiki/Home) is developed by the [Huttenhower Lab](http://huttenhower.sph.harvard.edu).
-   Nephele runs the [Whole Metagenome Shotgun (wmgx)](https://bitbucket.org/biobakery/biobakery_workflows/wiki/Home#!whole-metagenome-shotgun-wmgx) and [Visualization for Whole Metagenome Shotgun (wmgx\_vis)](https://bitbucket.org/biobakery/biobakery_workflows/wiki/Home#!visualization-for-whole-metagenome-shotgun-wmgx_vis) bioBakery workflows. The workflows are run using the default parameters.
-   The visualization pipeline is only run for datasets with **at least 3 samples**.
-   More information about the individual tools that make up the pipelines can be found on the [bioBakery wiki](https://bitbucket.org/biobakery/biobakery/wiki/Home).

User Options
------------

-   **Strainphlan:** Should strain profiling with [StrainPhlAn](http://segatalab.cibio.unitn.it/tools/strainphlan/) be run? Strain profiling can greatly increase the runtime of your job depending on the size and diversity of your samples. (Logical. Default: False)
-   **Project name:** A project name to go at the top of the html graphical output.

Output Files
------------

-   The [workflows tutorial](https://bitbucket.org/biobakery/biobakery/wiki/biobakery_workflows) goes through the pipelines step-by-step including information about [all intermediate and final output files](https://bitbucket.org/biobakery/biobakery/wiki/biobakery_workflows#rst-header-output-files). We list some of the output files that may be of interest to our users here, as well as any output files made or removed by Nephele.
-   **log files:**
    -   *logfile.txt*: contains the messages associated with the Nephele backend, such as transferring files
    -   [*anadama.log*](https://bitbucket.org/biobakery/biobakery/wiki/biobakery_workflows#rst-header-log-file): produced by the bioBakery wmgx workflow and contains [all the associated information and error messages](https://bitbucket.org/biobakery/biobakery/wiki/biobakery_workflows#rst-header-standard-output) from the analysis
    -   *wgmx\_vis/anadama.log*: produced by the bioBakery wmgx\_vis workflow
-   **renamed\_inputs:** If you submit paired-end data, Nephele makes renamed links to the sequence files to suit the workflow's convention in this directory.
-   [**kneaddata:**](https://bitbucket.org/biobakery/biobakery/wiki/biobakery_workflows#rst-header-quality-control-data)
    -   *main*: Nephele removes all FASTQ files produced by Kneaddata, so this folder will only contain log files for each sample.
    -   *merged/kneaddata\_read\_count\_table.tsv*: merged data file containing read counts for each step of the QC process for each input file
-   [**metaphlan2:**](https://bitbucket.org/biobakery/biobakery/wiki/biobakery_workflows#rst-header-taxonomic-profiling-data)
    -   *main*: Nephele removes all bowtie2 sam files produced by MetaPhlAn2. So, this folder only contains *sample\_name\_taxonomic\_profile.tsv*, a taxonomic profile for each sample
    -   *merged/metaphlan2\_taxonomic\_profiles.tsv*: merged taxonomic profiles for all samples
    -   *merged/metaphlan2\_species\_counts\_table.tsv*: total number of species identified for each sample
-   [**humann2:**](https://bitbucket.org/biobakery/biobakery/wiki/biobakery_workflows#rst-header-functional-profiling-data)
    -   *main*: for each sample, a file of gene family and pathway abundances, pathway coverage, and a log
    -   *merged/\*.tsv*: gene families, ecs, and pathways files for all samples merged into single files
    -   *merged/\*\_relab.tsv*: data sets normalized to relative abundance
    -   *counts/humann2\_feature\_counts.tsv*: feature counts of gene families, ecs, and pathways for all samples
    -   *counts/humann2\_read\_and\_species\_count\_table.tsv*: total species identified after filtering and total reads aligning (for nucleotide and translated search) for each sample
-   [**strainphlan:**](https://bitbucket.org/biobakery/biobakery/wiki/biobakery_workflows#rst-header-strain-profiling-data) if the Strainphlan option is chosen, contains [core output](https://bitbucket.org/biobakery/metaphlan2/overview#markdown-header-some-other-useful-output-files) for profiling up to 10 species found in the sample
    -   *RAxML.\**: trees generated for each species, may not exist if species are not found (more information can be found in the [StrainPhlAn manual](https://bitbucket.org/biobakery/metaphlan2/overview#markdown-header-metagenomic-strain-level-population-genomics))
    -   *clade\_name.fasta*: the alignment file of all metagenomic strains
    -   *\*.info*: general information like the total length of the concatenated markers (full sequence length), number of used markers, etc.
    -   *\*.polymorphic*: statistics on the polymorphic site, [details here](https://bitbucket.org/biobakery/metaphlan2/overview#markdown-header-some-other-useful-output-files)
    -   *\*.marker\_pos*: this file shows the starting position of each marker in the strains.
-   **wmgx\_vis:** If you submit at least 3 samples, the [html report from the visualization pipeline](https://bitbucket.org/biobakery/biobakery/wiki/biobakery_workflows#rst-header-id10) will be created here. It includes the software versions as well as the individual commands used.

Tools and References
--------------------

-   McIver, L. J., Abu-Ali, G., Franzosa, E. A., Schwager, R., Morgan, X. C., Waldron, L., ... Huttenhower, C. (n.d.). "bioBakery: a meta'omic analysis environment." Bioinformatics. <https://doi.org/10.1093/bioinformatics/btx754>
-   Nephele runs the [biobakery/workflows docker image](https://hub.docker.com/r/biobakery/workflows/) from September 2017, which lists the following software versions:
    -   kneaddata v0.6.1
    -   MetaPhlAn version 2.6.0 (19 August 2016)
    -   humann2 v0.11.1
