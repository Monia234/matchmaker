#!/usr/bin/env python

from __future__ import print_function

import sys
from sys import argv as args

from os import path

import ibd
import bed
import match

import jerrington_tools as j

import ibd_anc_plot_config as conf

from itertools import islice, imap
import operator as op

from PIL import Image, ImageDraw

def plot_matches(matches):
    """ Create an image configured as in ibd_anc_plt_config.py displaying the
        given list of matches. The list of matches must be ordered, but this is
        not checked!
        """
    def get_code_color(anc_code, inside_ibd):
        """ Get the color for the given ancestry code, with the opacity
            adjusted for whether we're inside an IBD region. """
        return anc_code.color + ((1.0,) if inside_ibd else (0.5,))

    im = Image.new("RGBA", conf.FIGURE_SIZE, "white")
    draw = ImageDraw.Draw(im)

    WIDTH_SCALE = conf.FIGURE_WIDTH / float(len(matches[0]))
    def scale_x(x):
        return x * WIDTH_SCALE

    # take the height of the figure, subtract by the number of margins (one
    # less than the number of matches) multiplied by the margin size, and
    # divide by the number of entries to plot.
    ENTRY_HEIGHT = (conf.FIGURE_HEIGHT -
            (len(matches) - 1) * conf.INTER_ENTRY_MARGIN) / float(len(matches))
    # simply two individuals in each match
    INDIVIDUAL_HEIGHT = ENTRY_HEIGHT / 2

    for (i, entry) in enumerate(matches):
        def rebase_x(x):
            return scale_x(x - entry.shared_segment.start)
        # the scaled width of this entry
        my_width = scale_x(len(entry))
        # the x origin for this entry
        my_x0    = conf.FIGURE_WIDTH / 2 - my_width / 2
        # the chromosome of this entry (rebound to save space)
        my_chr   = entry.chromosome
        for (indiv_i, individual) in enumerate(entry.individuals):
            # the y origin for this individual
            my_y0    = (i * (ENTRY_HEIGHT + conf.INTER_ENTRY_MARGIN)
                     + indiv_i * INDIVIDUAL_HEIGHT)
            # the haplotype code to use for this individual
            my_haplotype = bed.Individual.bed_code_from_IBD(
                    entry.ibd_segment.haplotype[indiv_i])
            # the ancestry segments on the IBD haplotype, on the IBD
            # chromosome, for this individual
            segments = individual.ancestries[my_haplotype][my_chr]
            # control variable: whether we are inside the IBD segment
            inside_ibd = False
            for seg in segments:
                def draw_seg(lower_bound, upper_bound, inside_ibd):
                    draw.rectangle([my_x0 + lower_bound, my_y0,
                        my_x0 + upper_bound, my_y0 + INDIVIDUAL_HEIGHT],
                        fill=get_code_color(seg.code, inside_ibd))

                lower_bound = rebase_x(seg.interval_bp.start)
                if entry.ibd_segment.start in seg:
                    upper_bound = rebase_x(entry.ibd_segment.start)
                elif entry.ibd_segment.end in seg:
                    upper_bound = rebase_x(entry.ibd_segment.end)
                else:
                    upperBound = rebase_x(seg.interval_bp.end)
                draw_seg(lower_bound, upper_bound, inside_ibd)
                if upper_bound < seg.interval_bp.end:
                    inside_ibd = not inside_ibd # invert the state
                    # set out bounds to the second part of the segment
                    lower_bound = upper_bound
                    upper_bound = seg.interval_bp.end
                    draw_seg(lower_bound, upper_bound, inside_ibd)
    return im

def n_most(seq, n, comp=op.lt):
    """ Return the n smallest elements in the seq. If n is greater than or
        equal to len(seq), this is just selection sort, and better sorting
        algorithms should be used instead.
        """
    outseq = list(seq) # copy the input sequence
    def swap(a, b):
        t = outseq[b]
        outseq[b] = outseq[a]
        outseq[a] = t

    s = outseq[0]
    for i in xrange(min(n, len(seq))):
        v = seq[i]
        for j in xrange(i + 1, len(outseq)):
            if comp(outseq[j], v):
                swap(i, j)
                break
    return outseq if n >= len(seq) else outseq[:n]

def main(project_dir):
    global PROJECT_DIR
    PROJECT_DIR = project_dir

    IBD_DIR = path.join(PROJECT_DIR,
            "baharian_projects/HRS/data/dbGaP/AfrAm/phased/3_GERMLINE")

    BED_DIR = path.join(PROJECT_DIR,
            "barakatt_projects/HRS/results/HRS_AFRAM_20140609/outbed")

    make_ibd_path = lambda i: path.join(IBD_DIR,
            "".join(["AfrAm.chr", str(i), ".IBD.match.gz"]))

    # a utility function
    flipcurry2 = j.compose(j.curry2, j.flip)

    # construct a function that takes an IBDEntry and generates the match
    # object from it.
    match_from_ibd_segment__ = match.IBDAncestryMatch.from_ibd_segment
    my_from_ibd_segment = j.supply(match_from_ibd_segment__, {"generate":True})
    match_from_ibd_segment = flipcurry2(my_from_ibd_segment)(BED_DIR)

    for chromosome_number in xrange(1, 2):
        ibd_chromosome_data = ibd.IBDEntry.from_GERMLINE(
                make_ibd_path(chromosome_number))
        matches = filter(lambda x: len(x) > 0, imap(
            match_from_ibd_segment, ibd_chromosome_data))
        print("Found", len(matches), "matches")

        longest_matches = n_most(matches, 120, op.gt)
        fig = plot_matches(longest_matches)

if __name__ == "__main__":
    if len(args) != 2:
        raise ValueError("incorrect number of arguments")
    main(args[1])
