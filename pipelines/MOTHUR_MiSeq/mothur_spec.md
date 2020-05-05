:orphan:
# mothur pipeline spec

.. sectionauthor:: Mariam Quiñones <quinonesm@niaid.nih.gov>


```bash
mothur "#make.contigs(file=combo.txt);
 summary.seqs(fasta=current);
 screen.seqs(fasta=current, maxambig=0, maxlength=262, group=current);
 unique.seqs();
 count.seqs(group=current);
 align.seqs(fasta=current, reference=silva.nr_v128.align, flip=T);
 summary.seqs();
 screen.seqs(optimize=start-end, criteria=90, fasta=current, count=current);
 filter.seqs(fasta=current, vertical=T, trump=.);
 unique.seqs(fasta=current, count=current);
 pre.cluster(diffs=2, fasta=current, count=current);
 chimera.vsearch(fasta=current, count=current, dereplicate=t);
 remove.seqs(fasta=current, accnos=current);
 summary.seqs(count=current);
 classify.seqs(count=current, reference={ref}, taxonomy={tax_ref}, cutoff=80, probs=f);
 remove.lineage(count=current, taxonomy=current, taxon={remove_lineage});
 summary.tax(taxonomy=current, count=current);
 dist.seqs(fasta=current, cutoff=0.03);
 cluster.split(fasta=current, count=current, taxonomy=current, splitmethod=classify, taxlevel=4, cutoff=0.03);
 make.shared(list=current, count=current, label=0.03);
 classify.otu(list=current, count=current, taxonomy=current, label=0.03);
 count.groups(count=current);
 make.biom(shared=current, constaxonomy=current, reftaxonomy={tax_ref});
 dist.seqs(fasta=current, output=lt);
 clearcut(phylip=current);"
```

# step by step (of above)

```bash
make.contigs(file=combo.txt)
summary.seqs(fasta=combo.trim.contigs.fasta)
screen.seqs(maxambig=0, maxlength=262, group=combo.contigs.groups)
unique.seqs()
```

Output File Names:

   - combo.trim.contigs.good.names

   - combo.trim.contigs.good.unique.fasta

  ```bash
  count.seqs(group=combo.contigs.good.groups)
  ```

  Output File Names:

   - combo.trim.contigs.good.count_table

```bash
align.seqs(fasta=combo.trim.contigs.good.unique.fasta, reference=silva.nr_v128.pcr.align, flip=T)
```

  Output File Names:

- combo.trim.contigs.good.unique.align
- combo.trim.contigs.good.unique.align.report

```bash
summary.seqs()
screen.seqs(optimize=start-end, criteria=90, fasta=combo.trim.contigs.good.unique.align, count=combo.trim.contigs.good.count_table)
```

Output File Names:

- combo.trim.contigs.good.unique.bad.accnos
- combo.trim.contigs.good.unique.good.align
- combo.trim.contigs.good.good.count_table

```bash
filter.seqs(fasta=combo.trim.contigs.good.unique.good.align, vertical=T, trump=.)
```

Output File Names:

- combo.filter

- combo.trim.contigs.good.unique.good.filter.fasta



```bash
unique.seqs(fasta=combo.trim.contigs.good.unique.good.filter.fasta, count=combo.trim.contigs.good.good.count_table)
```

Output File Names:

- combo.trim.contigs.good.unique.good.filter.count_table

- combo.trim.contigs.good.unique.good.filter.unique.fasta

```bash
pre.cluster(diffs=2, fasta=combo.trim.contigs.good.unique.good.filter.unique.fasta, count=combo.trim.contigs.good.unique.good.filter.count_table)
```

Output File Names:

- combo.trim.contigs.good.unique.good.filter.unique.precluster.fasta

- combo.trim.contigs.good.unique.good.filter.unique.precluster.count_table

- combo.trim.contigs.good.unique.good.filter.unique.precluster.*.map (one for each sample)



```bash
chimera.vsearch(fasta=combo.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=combo.trim.contigs.good.unique.good.filter.unique.precluster.count_table, dereplicate=t)
```

Output File Names:

- combo.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table

- combo.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.chimeras

- combo.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.accnos



```bash
remove.seqs(fasta=combo.trim.contigs.good.unique.good.filter.unique.precluster.fasta, accnos=combo.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.accnos)
```

Output File Names:

- combo.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta



```bash
summary.seqs(count=combo.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, fasta=combo.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta)
  classify.seqs(fasta=combo.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=combo.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, reference=silva.nr_v128.align, taxonomy=silva.nr_v128.tax, cutoff=80, probs=f)
```

Output File Names:

- combo.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v128.wang.taxonomy

- combo.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v128.wang.tax.summary



```bash
remove.lineage(count=combo.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, taxonomy=combo.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v128.wang.taxonomy, taxon=Chloroplast-Mitochondria-unknown-Eukaryota)
```

Output File Names:

- combo.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v128.wang.pick.taxonomy

- combo.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table

```bash
dist.seqs(fasta=combo.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, cutoff=0.03)
```

Output File Names:

- combo.trim.contigs.good.unique.good.filter.unique.precluster.pick.dist

```bash
cluster.split(fasta=combo.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=combo.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy=combo.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v128.wang.pick.taxonomy, splitmethod=classify, taxlevel=4, cutoff=0.03)
```

Output File Names:

- combo.trim.contigs.good.unique.good.filter.unique.precluster.pick.dist

- combo.trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.list

- combo.trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.sensspec

```bash
make.shared(list=combo.trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.list, count=combo.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=0.03)
```

Output File Names:

- combo.trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.shared

```bash
classify.otu(list=combo.trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.list, count=combo.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy=combo.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v128.wang.pick.taxonomy, label=0.03)
```

Output File Names:

- combo.trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.0.03.cons.taxonomy

- combo.trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.0.03.cons.tax.summary

```bash
count.groups(count=combo.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table)
make.biom(shared=combo.trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.shared, constaxonomy=combo.trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.0.03.cons.taxonomy, reftaxonomy=silva.nr_v128.tax)
```

Output File Names:

- combo.trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.0.03.biom

```bash
dist.seqs(fasta=combo.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, output=lt)
```

Output File Names:

- combo.trim.contigs.good.unique.good.filter.unique.precluster.pick.phylip.dist

```bash
clearcut(phylip=combo.trim.contigs.good.unique.good.filter.unique.precluster.pick.phylip.dist)
```

Output File Names:

- combo.trim.contigs.good.unique.good.filter.unique.precluster.pick.phylip.tre
