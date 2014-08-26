#!/usr/bin/env python

import sys
from sys import argv as args

from os import path

import ibd
import bed
import match

import jerrington_tools as j

PROJECT_DIR = "/lb/project/gravel"
IBD_DIR = path.join(PROJECT_DIR,
        "baharian_projects/HRS/data/dbGaP/AfrAm/phased/3_GERMLINE")

BED_DIR = path.join(PROJECT_DIR,
        "barakatt_projects/HRS/results/HRS_AFRAM_20140609/outbed")

make_ibd_path = lambda i: path.join(IBD_DIR,
        "".join(["AfrAm.chr", str(i), ".IBD.match.gz"]))

def main(args):
    ibd_chromosomes = map(
            j.compose(ibd.IBDEntry.from_GERMLINE, make_ibd_path),
            xrange(1, 2)) # just chromosome 1 for testing

    flipcurry2 = j.compose(j.curry2, j.flip)

    match_from_ibd_segment__ = match.IBDAncestryMatch.from_ibd_segment
    my_from_ibd_segment = j.supply(match_from_ibd_segment, {"generate":True})

    match_from_ibd_segment = flipcurry2(my_from_ibd_segment)(BED_DIR)

    # now match_from_ibd_segment is a unary function that takes one IBDEntry
    # instance and produces from it

    matches = map(j.map_c(match_from_ibd_segment), ibd_chromosomes)

if __name__ == "__main__":
    main(args)
