#!/usr/bin/env python3

import pysam
import sys

# give it an ID
sample_id = sys.argv[1]

# and a threshold
#threshold = sys.argv[2]
threshold = .5

# read vcf
vcf = pysam.VariantFile(sample_id + "_filtered_snps.vcf", "r")

# create 2 vcfs
major = pysam.VariantFile(sample_id + "_consensus_snps.vcf", 'w', header=vcf.header)
minor =  pysam.VariantFile(sample_id + "_filtered_snps_eaf.vcf", 'w', header=vcf.header)

for var in vcf:
	for sample in var.samples:
		ad = var.samples[sample]['AD']
		dp = var.samples[sample]['DP']
	
	af = var.info['AF']

	new_afs = []
	# first ad value is for ref
	ref_af = ad[0] / dp
	for i in range(len(ad))[1:]:
		eaf = ad[i] / dp
		new_afs.append(eaf)
	var.info['AF'] = tuple(new_afs)
			
	# add ALL variants (with updated AF) to minor file, regardless
	minor.write(var)
			
	# add only variants which PASS filter and have the highest frequency to major
	if not 'PASS' in var.filter: continue
	# biallelic
	if len(ad) == 2:
		if new_afs[0] > threshold:
			major.write(var)
	
	# multiallelic
	elif len(ad) > 2:
		if not max(new_afs) > ref_af: continue
		# if we're here, one of the alt alleles has an AF greater than the ref, which one?
		af_dict = {}
		af_dict['ref'] = ref_af
		for i in range(len(new_afs)):
			af_dict[var.alleles[i+1]] = new_afs[i]
		max_allele = max(af_dict, key=af_dict.get)
		# var.alleles returns (ref, alt_1, alt_2, ...)
		# set it to (ref, alt) where alt is max_allele
		var.alleles = (var.alleles[0], max_allele)
		new = major.new_record(contig = var.contig, start = var.start, stop = var.stop, alleles = var.alleles, id = var.id, qual = var.qual, filter = var.filter, info = {"AF": af_dict[max_allele]})
		major.write(new)

	else:
		print("len(ad) < 2 (?). ad = ", ad)
		sys.exit(1)

