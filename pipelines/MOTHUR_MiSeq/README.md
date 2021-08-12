# Synopsis #
``mothur`` is an open source software package for bioinformatics data processing initiated by
[Patrick Schloss](https://medicine.umich.edu/dept/microbiology-immunology/patrick-schloss-phd) and his development team in
the [Department of Microbiology and Immunology at the University of Michigan](https://medicine.umich.edu/dept/microbiology-immunology).
The ``mothur`` project aims to develop a single piece of open source, expandable software to fill the bioinformatics needs of
the microbial ecology community. ``mothur`` incorporates the functionality of [DOTUR](https://aem.asm.org/content/71/3/1501),
[SONS](https://aem.asm.org/content/72/10/6773), [TreeClimber](https://aem.asm.org/content/72/4/2379),
[s-libshuff](https://aem.asm.org/content/70/9/5485.long) and others. The current version of ``mothur`` on Nephele is 1.40.5.
The ``mothur`` pipeline follows the standard operating procedure (SOP) published by the
[Schloss lab](https://www.mothur.org/wiki/MiSeq_SOP) to process 16S rRNA gene sequences. These sequences are assumed
generated using [Illumina MiSeq platform](https://www.illumina.com/systems/sequencing-platforms/miseq.html) with either
paired or single ends. The SOP is largely the product of a series of manuscripts published by the Schloss lab. The MiSeq
SOP is described in the manuscript:

Kozich JJ, Westcott SL, Baxter NT, Highlander SK, Schloss PD (2013) Development of a dual-index sequencing strategy and
curation pipeline for analyzing amplicon sequence data on the MiSeq Illumina sequencing platform. *Appl Environ Microbiol*,
79(17): 5112-51120. [doi: 10.1128/AEM.01043-13](https://aem.asm.org/content/79/17/5112.long).

# Example Usage #

The following shows an example usage of the ``mothur`` pipeline script [m2.py](m2.py):
```bash
m2.py --input_dir input --output_dir output --data_type PE --map_file map_file.csv --maxlength 0 --picrust --optimize start-end --criteria 80 --remove_lineage Archaea-Chloroplast-Mitochondria-Eukaryote-Unknown --sampling_depth 1500 --ref_db sv99
```

# Pipeline Options #
- ``--job_id``: job id
- ``--inputs_dir``: input folder
- ``--outputs_dir``: output folder
- ``--data_type``: sequence data type; the default is paired-end (PE)
- ``--map_file``: sample mapping file

**User options**
- ``--maxlength``: ``maxlength`` largely works as quality assurance on the sequences. This option removes sequences that are
either too long or too short. Excessive long sequences are likely due to incorrect pairing while short sequences indicate
low quality reads. The algorithm constructs a histogram based on the sequence lengths and determine the most appropriate
cutoff for analysis. See the [<b>Max Length</b>](#maxlength) section below for more information.
- ``--picrust``: run [PICRUSt](https://github.com/picrust/picrust2/wiki). PICRUSt (Phylogenetic Investigation of Communities
by Reconstruction of Unobserved States) is a software for predicting functional abundances based only on marker gene
sequences.
- ``--optimize ``: remove sequences that fail certain criteria. Options are ``start-end``, ``start``, and ``end``.
- ``--criteria``: optimize aligned sequences. This option will remove any sequence that starts after the given position, or
ends before the position. This parameter ensures that sequences align to the same region.
- ``--remove_lineage``: remove taxa that belong to the lineages such as archaea, chloroplast, mitochondria, eukaryote, or
unknown. The default lineages are chloroplast, mitochondria, eukaryote, and unknown.
- ``--sampling_depth``: sampling depth for downstream analysis. The default value is ``None``.
- ``--ref_db``: reference database; either SILVA v128 99 (sv99) or HOMD (homd) databases.

**Internal/testing options**
- ``--keep``: keep intermediate files
- ``--maxee``: maxee param for ``[make.contigs](https://mothur.org/wiki/Make.contigs)``; if specified, param will be passed,
otherwise (default) ``make.contigs`` runs without this param.

# Analysis Steps and Commands
The [original pipeline specification](mothur_spec.md) (:doc:`mothur_spec`) has been modified to contend with the realities of our user's data.
The current steps are below:
- ``make.contigs(file=combo.txt)``: reads forward and reverse FASTQ files and output new FASTA and report files. This command
assumes the file ``combo.txt`` lists the sample IDs and corresponding forward and reverse FASTQ files.
- ``rename.seqs(file=current, group=current)``: rename the sequences appending the group name to the sequence name. <b>Note</b>:
This command should reduce the size of distance matrix for the downstream analysis.
- ``summary.seqs(fasta=current)``: summarize the quality of sequences in an unaligned or aligned fasts formatted sequence
file.
- ``screen.seqs(fasta=current, maxambig=0, maxlength=262, group=current)``: remove sequences that exceed a certain length.
The ``maxlength`` largely works as quality assurance on the sequences. This option removes sequences that are either too long
or too short. Excessive long sequences are likely due to incorrect pairing while short sequences indicate low quality reads.
The algorithm constructs a histogram based on the sequence lengths and determine the most appropriate cutoff for analysis.
- ``unique.seqs()``: deduplicate sequences.
- ``count.seqs(group=current)``: generate descriptive statistics of the current collection of sequences.
- ``align.seqs(fasta=current, reference=silva.nr_v128.align, flip=T)``: taxonomically assign contigs using the
[custom SILVA database](https://mothur.org/wiki/Silva_reference_files). The
[HOMD (Human Oral Microbiome Database)](http://www.homd.org/) is also available.
- ``screen.seqs(optimize=start-end, criteria=90, fasta=current, count=current)``: remove sequences that fail certain criteria.
The ``start`` and ``end`` arguments remove sequences that do not start and end at given positions based on the alignments. These
arguments are determined based on the median start and end positions of sequences. The ``criteria`` argument will remove any
sequences that starts after the position that 90% of the sequences do, or ends before the position that 90% of the sequences
do, or whose length is shorter than 90% of the sequences.
- ``filter.seqs(fasta=current, vertical=T, trump=.)``: remove columns from alignments with either a ``-`` or ``.`` characters.
These columns are not included in calculating distances because they have no information in them.
- ``unique.seqs(fasta=current, count=current)``: remove duplicate sequences and reduce the comping time in the downstream
analysis.
- ``pre.cluster(diffs=2, fasta=current, count=current)``: remove sequences that are likely due to pyrosequencing errors. This
command implements a pseudo-single linkage algorithm that aims to reduce sequencing errors. The original algorithm was
published by [Sue Huse](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2909393/). The basic idea is that abundant sequences
are more likely to generate erroneous sequences than rare sequences. The algorithm proceeds by ranking sequences in order of
their abundance. Sequences with abundance of a certain threshold are merged with the larger collection of sequences.
- ``chimera.vsearch(fasta=current, count=current, dereplicate=t)``: reads fasta and reference file, and output potentially
chimeric sequences using [vsearch](https://github.com/torognes/vsearch). ``vsearch`` uses a fast heuristic based on words
shared by the query and target sequences in order to quickly identify similar sequences. The algorithm performs optimal
global sequence alignment of the query against potential target sequences, using full dynamic programming instead of the
seed-and-extend heuristic approach. Pairwise alignments are computed in parallel using vectorization and multiple threads.
- ``remove.seqs(fasta=current, accnos=current)``: generate a new file that does not contain the sequences specified in the
argument ``accnos``.
- ``classify.seqs(count=current, reference=ref, taxonomy=tax_ref, cutoff=80, probs=f)``: assign contigs to the taxonomy
outline. This command implements the *k*-nearest neighbor consensus and Wang approach. Taxonomy outlines and reference
sequences can be obtained from the [taxonomy outline](https://www.mothur.org/wiki/Taxonomy_outline) page.
- ``remove.lineage(count=current, taxonomy=current, taxon=lineage)``: reads the taxonomy and taxon files, and generates a new
file that contains only the sequences not containing the taxon.
- ``summary.tax(taxonomy=current, count=current)``: reads a taxonomy file and an optional name and or group file, and
summarizes the taxonomy information.
- ``cluster.split(fasta={fasta}, count={counts}, taxonomy={taxonomy}, splitmethod=classify, taxlevel=4, large=T, cluster=F, cutoff=0.03)``:
split sequences into groupings based on taxonomy and compute distance matrices. Each matrix contains pairwise distances
between sequences in grouping. Outputs ``*.pick.file`` with distance matrices and count tables for each grouping.
- We check the sizes of the distance matrices, and if too large, exit the pipeline (See
:py:func:`~nephele2.pipelines.MOTHUR_MiSeq.m2.MothurNeph.check_dist_matrix_size`).
- ``cluster.split(file={pick_file}, cutoff=0.03, runsensspec=f)``: cluster sequences from each grouping in into OTUs using
default opticlust method, and merge result . outputs merged distance matrix and list file of sequences in each OTU.
See [cluster.split](https://mothur.org/wiki/Cluster.split).
- ``make.shared(list=current, count=current, label=0.03)``: reads a list and group file or BIOME file, and create a shared
file as well as a ``rabund`` file for each group.
- ``classify.otu(list=current, count=current, taxonomy=current, label=0.03)``: get a consensus taxonomy for an OTU.
- ``make.biom(shared=current, constaxonomy=current, reftaxonomy=tax_ref)``: convert the shared file to a
[BIOM file](http://biom-format.org/documentation/biom_format.html). ``mothur`` only support format version 1.0.
- ``dist.seqs(fasta=current, output=lt)``: calculate uncorrected pairwise distances between aligned sequences. The ``lt``
option will output a phylip formatted lower triangle matrix, or to <q>square</q> for a phylip formatted square matrix. The
phylip output is required for the subsequent ``clearcut`` command.
- ``clearcut(phylip=current)``: run the [clearcut algorithm](https://github.com/ibest/clearcut/). The algorithm was
developed by the [Institute of Bioinformatics and Evolutionary Studies](https://www.ibest.uidaho.edu/) at the University of
Idaho. The clearcut program implements the relaxed neighbor joining algorithm and achieves speedup with a significant
reduction in the quality of the inferred trees. Relaxed neighbor joining algorithm was developed by Jason Evans, Luke
Sheneman, and James Foster at the University of Idaho.

# Max Length
The ``maxlength`` parameter largely works as quality assurance on the contigs. The ``maxlength`` option removes contigs that
are either too long or too short. Excessive long contigs are likely due to incorrect pairing while short contigs indicate
low quality reads. In general, the length a PCR fragment is determined by the locations where the forward and reverse primers
bind on the template sequences. However, the length of PCR fragments is also limited by the sequencing instrument. As an
example, if the instrument generate reads of 250bp, then the joined contig length cannot possibly exceed 450bp with an
overlap region of 25bp. The ``mothur`` pipeline implements a heuristic algorithm to determine the most probable maximum length
of contigs. The algorithm constructs a histogram based on the sequence lengths and determine the most appropriate cutoff for
downstream analysis.

![Histogram of contigs](histogram.png "The length of the most abundant contigs is used to determine the cutoff length (`maxlength`)")
The length of the most abundant contigs is used to determine the cutoff length (`maxlength`)

In many cases, the actual maximum length of contigs is difficult to determine. If samples are heavily contaminated, then
sequencing instrument might not be able to generate reads with sufficient quality. In this case, many nucleotides in the
reads will get removed, which inevitably shortens the configs prematurely. The algorithm, although rudimentary, seems to
perform satisfactory.

# References #
* Kozich JJ, Westcott SL, Baxter NT, Highlander SK, Schloss PD (2013) Development of a dual-index sequencing strategy and
curation pipeline for analyzing amplicon sequence data on the MiSeq Illumina sequencing platform. *Appl Environ Microbiol*.
[doi: 10.1128/AEM.01043-13](https://aem.asm.org/content/79/17/5112.long).
* PD Schloss, SL Westcott, T Ryabin, JR Hall, M Hartmann, EB Hollister, RA Lesniewski, BB Oakley, DH Parks, CJ Robinson,
JW Sahl, B Stres, GG Thallinger11, DJ Van Horn and CF Weber (2009) Introducing mothur: open-source, platform-independent,
community-supported software for describing and comparing microbial communities, *Appl Envion Microbiol*.
[doi: 10.1128/AEM.01541-09](https://aem.asm.org/content/75/23/7537).
* MGI Langille, J Zaneveld, JG Caporaso, D McDonald, D Knights, JA Reyes, JC Clemente, DE Burkepile, RL Vega Thurber,
R Knight, RG Beiko and C Huttenhower (2013) Predictive functional profiling of microbial communities using 16S rRNA marker
gene sequences. *Nature Biotechnology*. [doi: 10.1038/nbt.2676](https://www.nature.com/articles/nbt.2676).
* J Evans, L Sheneman, and JA Foster. Relaxed neighbor-joining: a fast distance-based phylogenetic tree construction method,
*Journal of Molecular Evolution*, 62(6): 785-792.
[doi: 10.1007/s00239-005-0176-2](https://link.springer.com/article/10.1007/s00239-005-0176-2).


*Last updated on February 6, 2020*.
