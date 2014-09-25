#!/usr/bin/env python

from __future__ import print_function

import jerrington_tools as jt
import match
import dataset_utils
import ibd
import bed

from itertools import imap, ifilter, chain

ibd_paths=["project/baharian_projects/MergedData/phased/3_GERMLINE/cMcorrected/MERGED_chr1.cM.IBD.match.gz"]
outbed_path="project/barakatt_projects/HRS/results/SCCS_notphased_20140923/outbed"

def get_chromosome_data(handles, filterf):
    """ Construct a generator to yield all the IBD entries for the
        African-American HRS dataset. We are pulling Soheil's MergedData
        dataset, so the GERMLINE output files need to be filtered to remove
        entries from the other datasets.

        Arguments:
            handles (list of file handles):
                Handles on all the files to load.
        """
    return ifilter(filterf, chain(*map(ibd.IBDEntry.ifrom_GERMLINE,
        handles)))

# a utility function
flipcurry2 = jt.compose(jt.curry2, jt.flip)

# construct a function that takes an IBDEntry and generates the match
# object from it.
match_from_ibd_segment__ = match.IBDAncestryMatch.from_ibd_segment
my_from_ibd_segment = jt.supply(match_from_ibd_segment__,
        {"generate":True, "cache":True,
            "filename_parserf":dataset_utils.sccs_name_parser})
match_from_ibd_segment = flipcurry2(my_from_ibd_segment)(outbed_path)

# Open the relevant files
handles = map(jt.maybe_gzip_open, ibd_paths)

# compute the ancestry matches for those individuals
matches = filter(lambda x: len(x) > 0, imap(
    match_from_ibd_segment,
    get_chromosome_data(handles, dataset_utils.is_sccs)))

# close the files now that they've been read and parsed
map(lambda x: x.close(), handles)

print("Found", len(matches), "matches")

# Run the actual check.
for match in matches:
    for individual in match.individuals:
        for hcode in bed.Individual.HAPLOTYPE_CODES:
            debug_bed = individual.to_debugstr(hcode).split('\n')
            bed_file_path = path.join(outbed_path,
                    id_to_bedfile(individual.name, hcode))
            h = open(bed_file_path)
            file_bed = [line for line in h]
            h.close()
            for i in xrange(len(file_bed)):
                if not file_bed[i].startswith(debug_bed[i]):
                   print("Local ancestry mismatch!\n",
                       "\tINTERNAL: " + debug_bed[i] + "\n",
                       "\tFILE:" + file_bed[i] + "\n",
                       file=sys.stderr, sep="")
                   raise Exception()
