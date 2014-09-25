#!/usr/bin/env python

import jerrington_tools as jt
import match
import bed
import ibd
import dataset_utils

def main(ibd_paths, outbed_path, ibd_filterf):
    matches = match.from_ibds_and_bedpath(ibd_paths, outbed_path, ibd_filterf)

    ancestry_sizes = []

    map(lambda m: ancestry_sizes.extend(m.calculate_ibd_ancestry_sizes()),
            matches)

    # determine the sum of all the sizes of the ibd segments
    total_ibd_length = float(sum(len(match.ibd_segment) for match in matches))

    # Initialize the code counts to zero
    total_sizes = {}
    for code in bed.AncestryCode.CODENAMES:
        total_sizes[code] = 0

    for ancestry_size in ancestry_sizes:
        for (code, size) in ancestry_sizes.items():
            total_sizes[code] += size

    for (code, size) in total_sizes.items():
        print(code, size, sep='\t')

if __name__ == "__main__":
    main(["project/baharian_projects/MergedData/phased/3_GERMLINE/cMcorrected/MERGED_chr1.cM.IBD.match.gz"],
        "project/barakatt_projects/HRS/results/SCCS_notphased_20140923/outbed",
        dataset_utils.is_sccs)
