#!/usr/bin/env python

from __future__ import print_function

import bed
import dataset_utils

outbed_dir = "project/barakatt_projects/HRS/results/SCCS_notphased_20140923/outbed"

testfile = outbed_dir + "TGWAS_0001_A_cM.bed"

i = bed.Individual.from_dir_and_name(outbed_dir, "GWAS_0001", dataset_utils.sccs_name_parser)

print("Haplotype A")
print(i.to_debugstr("A"))

print("Haplotype B")
print(i.to_debugstr("B"))
