#!/usr/bin/env python

from __future__ import print_function

import sys
from sys import argv as args

from itertools import imap

import jerrington_tools as jt
import match
import bed
import ibd
import dataset_utils

def show_usage():
    errprint = jt.supply(print, {"file":sys.stderr})
    map(errprint, [
        "calculate_ancestry_proportions_ibd.py -- calculate the proportion of",
        "  each ancestry in a given population within IBD regions.",
        "",
        "usage: ./calculate_ancestry_proportions_ibd.py",
        "  (-d|--dataset) <dataset filter> (-b|--beddir) <rfmix outbed directory>",
        "  ((-i|--ibd) <germline output file for one chromosome>)+",
        ""])

    map(errprint, ["Arbitrarily many IBD chromosome files can be given, in which ",
        "case statistics are calculated for each chromosome and then averaged."])

def main(ibd_paths, outbed_path, ibd_filter):
    if ibd_filter is dataset_utils.is_afram_hrs:
        INDIVS = set(jt.with_file(jt.map_c(lambda x: x[:-1]), indiv_list_path))
        ibd_filterf = jt.curry2(ibd_filter)(INDIVS)
    else:
        ibd_filterf = ibd_filter

    handles = imap(jt.maybe_gzip_open, ibd_paths)

    matches = match.IBDAncestryMatch.from_ibds_and_bedpath(
            handles, outbed_path, ibd_filterf)

    # Initialize the code counts to zero
    total_sizes = {}
    for code in bed.AncestryCode.CODENAMES: total_sizes[code] = 0

    def elemwise_append(d1, d2):
        """ Append the values in d2 to the corresponding ones in d1.
            d1 is returned for chaining.
            """
        for (k, v) in d2.items():
            d1[k] += v
        return d1

    (total_ibd_length, total_sizes) = reduce(
            lambda (ibd_len, tot_size), m: (
                ibd_len + len(m.ibd_segment),
                elemwise_append(tot_size, m.calculate_ibd_ancestry_sizes())),
            matches, (0, total_sizes))

    print(total_sizes)

    for (code, size) in total_sizes.items():
        print(code, size / float(total_ibd_length), sep='\t')

    for x in handles: x.close()

if __name__ == "__main__":
    if len(args) == 1:
        print("No command-line arguments given.", file=sys.stderr)
        show_usage()
        sys.exit(1)

    if jt.any_eq(args[1], ["-h", "--help"]):
        show_usage()
        sys.exit(0)

    ibds = []
    beddir = None
    dataset = "id"

    i = 1
    while i < len(args):
        arg = args[i]
        nextarg = lambda: args[i + 1]
        if jt.any_eq(arg, ["-i", "--ibd"]):
            ibds.append(nextarg())
            i += 1
        elif jt.any_eq(arg, ["-b", "--bed-dir"]):
            beddir = nextarg()
            i += 1
        elif jt.any_eq(arg, ["-d", "--dataset"]):
            dataset = nextarg()
            i += 1
        else:
            print("Unrecognized command-line argument:", arg, file=sys.stderr)
            show_usage()
            sys.exit(1)
        i += 1

    if not (ibds and beddir and dataset):
        print("One or more missing command-line arguments.", file=sys.stderr)
        show_usage()
        sys.exit(1)

    filters = {"sccs":dataset_utils.is_sccs, "hrs":dataset_utils.is_afram_hrs,
            "id":jt.const(True)}

    if not jt.any_eq(dataset, ["sccs", "hrs", "id"]):
        print("Unknown dataset filter. Known dataset filters are:\n" +
                ", ".join(filters.keys()) +
                "\nThe 'id' filter is a filter that simply"
                "does not filter elements at all.", file=sys.stderr)
        show_usage()
        sys.exit(1)

    main(ibds, beddir, dataset_utils.is_sccs)
