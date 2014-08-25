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
            xrange(1, 2)) # 1 to 22

    matches = map(
            match.IBDAncestryMatch


if __name__ == "__main__":
    main(args)
