#### QIIME2 workflow

***The new QIIME2 pipeline will be based on the QIIME2 tutorials available on docs.qiime2.org***
***It will allow users to select between a method of denoising (deblur) or the revised clustering methods similar to the previous QIIME1 (open, closed and denovo)***

The code is on the provided jupyter notebook file: ***[qiime2_nephele_paired_pipe_nodada.ipynb](https://github.niaid.nih.gov/quinonesm/nephele2/blob/next_release/pipelines/QIIME2_16S/qiime2_nephele_paired_pipe_nodada.ipynb)*** 

For reference, see QIIME2 pipeline diagrams: https://docs.qiime2.org/2020.11/tutorials/overview/#let-s-get-oriented-flowcharts 
in particular https://docs.qiime2.org/2020.11/tutorials/overview/#derep-denoise

Also see basic tutorial: https://docs.qiime2.org/2020.11/tutorials/moving-pictures/ 

#### Reference files that need to be obtained
1. Greengenes 13_8 OTU files (99, 97, 85) - FTP URL|[ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_8_otus.tar.gz]
2. Tree files for insertion - download both files (greengenes and silva) from https://docs.qiime2.org/2020.11/data-resources/#sepp-reference-databases (sepp-refs-silva-128.qza or sepp-refs-gg-13-8.qza)
3. Taxonomy sklearn classifier for full 16S using OTUs 99%
```
wget https://data.qiime2.org/2020.11/common/silva-138-99-nb-classifier.qza
wget https://data.qiime2.org/2020.11/common/gg-13-8-99-nb-classifier.qza
```
4. Taxonomy skleran classifier for partial (V4 region only) 16S using OTUs 99%
```
wget https://data.qiime2.org/2020.11/common/gg-13-8-99-515-806-nb-classifier.qza
wget https://data.qiime2.org/2020.11/common/silva-138-99-515-806-nb-classifier.qza
```
5. Taxonomy reads and taxonomy assignment to be used with vsearch classify
``` 
wget https://data.qiime2.org/2020.11/common/silva-138-99-seqs.qza
wget https://data.qiime2.org/2020.11/common/silva-138-99-tax.qza
```
___________
***We will limit the file inputs to demultiplexed single, prejoined paires or paired end fastq files.  The filenames will be listed on the mapping files.  The validation of the mapping file will be different than what was established for QIIME version 1.0 .  For paired-end libraries, the mapping file needs one column with header "sample-id", a second column with header "forward-absolute-filepath" and a third columns called "reverse-absolute-filepath".  For single-end, the file only needs columns "sample-id" and "absolute-filepath".  Additional metadata files can be provided as long as the headers are with Unicode characters***

#### Part 1: File Import
##### Part 1a: Import of demultiplexed PairedEnd reads using indicated filepaths
```
--p-min-quality (default:20)***
```

##### Part 1b: Import of demultiplexed SingleEnd reads using indicated filepaths 
```
--p-min-quality (default:20)***
```
##### Part 1c: Import of demultiplexed prejoined reads using indicated filepaths
```
--p-min-quality (default:20)***
```
#### Part 2: Clustering or Denoising
***User will need to have the option of selecting Clustering (vsearch) or Denoising (Deblur). For clustering, users will first have the option of running chimera removal *Yes or no*, then they need select from  denovo, open or closed reference with a default selection of closed-reference.  Users will then need to select between Greengenes or SILVA 85, 97 or 99% identity (default 97).***

##### Part 2a: Denovo clustering
```
--p-perc-identity (default:0.97)***
```
##### Part 2b: Open-Reference clustering
```
--p-perc-identity (default:0.97)
--i-reference-sequences (default:97_otus.qza)***
```
##### Part 2c: Closed-Reference clustering
```
--p-perc-identity (default:0.97)***
--i-reference-sequences (default:97_otus.qza)***
```
#### Part 3: Denoising with Deblur.  
***It only works with single end data (or previously merged pairs and treated as single end)***
```
--p-trim-length (default 200). Sequence trim length, specify -1 to disable trimming. (default 200)
```
#### Part 4: Create a phylogenetic tree to be used in downstream diversity steps
***This step will take input from the Denoising steps or clustering steps above.  There are two options of methods: align to tree using mafft or fragment-insertion using fasttree***

#### Part 5: Alpha and Beta Diversity - This is similar as the Downstream Analysis pipeline except that it uses the phylogenetic tree
```
--p-sampling-depth*** (allow user to enter but when empty calculate automatically as done in other Nephele pipelines)
```
#### Part 6: Alpha rarefaction
```
--p-max-depth*** (allow user to input but when empty use same value as --p-sampling-depth but multiplied by 2)
```
#### Part 7: Taxonomy classification and barplots
***This step will allow for use of premade classifiers trained on 99% OTUs from Greengenes and SILVA databases.  In addition, by default we will classify using the classifier trained on the entire 16S but a faster and more appropriate option will be made available for V4 region (primers 515-806)***

##### Part 7a:Classify method sklearn
```
--p-reads-per-batch 100 (default is 100 to speed up the analysis somewhat)
--i-classifier (provide V4 options: silva-138-99-515-806-nb-classifier.qza, gg-13-8-99-515-806-nb-classifier.qza and full length options silva-138-99-nb-classifier.qza and gg-13-8-99-nb-classifier.qza)
```
##### Part 7b:Classify method vsearch with SILVA database reference
***(maybe we can prepare the greengenes files but this is for future)***

#### Part 8:Barplots
```
--p-min-frequency*** (allow user input or use same value as --p-sampling-depth)
--i-taxonomy*** Select greengenes or SILVA database
```
#### Part 9: Optional analyses 
***(I wonder if these can be done here or in the downstream analysis pipeline).  It will need metadata file columns headers to be presented to the user for the user to select from.  For each of the metadata colums, the beta-group-significance function will be run.***
