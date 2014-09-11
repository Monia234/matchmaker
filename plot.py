#!/usr/bin/env python

from __future__ import print_function

import sys
from sys import argv as args

from os import path

import ibd
import bed
import match

import jerrington_tools as j
jt = j

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
        return anc_code.color + ((255,) if inside_ibd else (128,))

    im = Image.new("RGBA", conf.FIGURE_SIZE, "white")
    draw = ImageDraw.Draw(im)

    ibd_center = lambda m: (m.ibd_segment.interval.start +
            m.ibd_segment.interval.end) / 2
    shared_center = lambda m: (m.shared_segment.start +
            m.shared_segment.end) / 2
    delta = lambda m: shared_center(m) - ibd_center(m)

    min_bp = min(imap(lambda m: m.shared_segment.start + delta(m), matches))
    max_bp = max(imap(lambda m: m.shared_segment.end + delta(m), matches))
    total_width = max_bp - min_bp
    print("WIDTH: ", total_width)

    WIDTH_SCALE = conf.FIGURE_WIDTH / float(total_width)
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
        jt.errprint("ENTRY:", i)
        # the scaled width of this entry
        ibd_width_true = len(entry.ibd_segment.interval)
        ibd_width = scale_x(len(entry.ibd_segment.interval))
        jt.errprint("\tIBD WIDTH:", ibd_width_true, "->", ibd_width)

        # where the starting point of the IBD segment needs to be drawn for it
        # to be centered
        ibd_x0    = conf.FIGURE_WIDTH / 2 - ibd_width / 2
        ibd_start_x = scale_x(entry.ibd_segment.interval.start)
        ibd_center_offset = ibd_x0 - ibd_start_x

        jt.errprint("\tSEGMENT OFFSET:", ibd_center_offset)

        # the ibd center offset needs to be added to all the positions, so that
        # the ibd segment is properly centered in the figure.

        # the chromosome of this entry (rebound to save space)
        my_chr   = entry.chromosome
        for (indiv_i, individual) in enumerate(entry.individuals):
            jt.errprint("\tINDIVIDUAL:", indiv_i)

            # the y origin for this individual
            my_y0    = (i * (ENTRY_HEIGHT + conf.INTER_ENTRY_MARGIN)
                     + indiv_i * INDIVIDUAL_HEIGHT)

            # the haplotype code to use for this individual
            my_haplotype = bed.Individual.bed_code_from_IBD(
                    entry.ibd_segment.haplotype[indiv_i])

            # the ancestry segments on the IBD haplotype, on the IBD
            # chromosome, for this individual
            my_chromosome = individual.ancestries[my_haplotype][my_chr]
            segments = my_chromosome.segments

            # control variable: whether we are inside the IBD segment
            inside_ibd = False
            for (j, seg) in enumerate(segments):
                # construct the drawing function for this segment
                def draw_rect(upper_bound, lower_bound, inside_ibd):
                    segment_width_true = upper_bound - lower_bound
                    segment_width = scale_x(segment_width_true)
                    rect = Image.new("RGBA",
                            (int(segment_width), int(INDIVIDUAL_HEIGHT)))
                    rect_draw = ImageDraw.Draw(rect)
                    rect_draw.rectangle([0, 0, segment_width, INDIVIDUAL_HEIGHT],
                        fill=get_code_color(seg.code, inside_ibd))
                    xstart = int(scale_x(lower_bound) + ibd_center_offset)
                    xend   = xstart + segment_width
                    im.paste(rect,
                            (xstart, int(my_y0)),
                            mask=rect)
                    return xstart, xend

                # identify the width of the segment, accounting for possible
                # IBD segment beginnings or ends within it.
                lower_bound = seg.interval_bp.start
                if entry.ibd_segment.interval.start in seg.interval_bp:
                    jt.errprint("\t\tIBD segment starting in this segment.")
                    upper_bound = entry.ibd_segment.interval.start
                elif entry.ibd_segment.interval.end in seg.interval_bp:
                    jt.errprint("\t\tIBD segment ending in this segment.")
                    upper_bound = entry.ibd_segment.interval.end
                else:
                    # the fact that upper_bound /= seg.interval_bp.end will be
                    # used to determine that a cut was made in this segment.
                    # The other portion of the segment will be drawn a bit
                    # later.
                    upper_bound = seg.interval_bp.end

                # draw the rectangle
                (xstart, xend) = draw_rect(upper_bound, lower_bound, inside_ibd)

                jt.errprint("\t\tSEGMENT: #", j, ": (", lower_bound, ", ",
                        upper_bound, ") -> (", xstart, ", ", xend, ")", sep='')

                # if the upper bound does not reach to the end of the segment,
                # then we conclude that the segment was split due to the
                # presence of an IBD segment boundary within it.
                if upper_bound < seg.interval_bp.end:
                    # this represents a change in whether we're inside the IBD
                    inside_ibd = not inside_ibd # so we invert the state
                    # set out bounds to the second part of the segment
                    lower_bound = upper_bound

                    if entry.ibd_segment.interval.end in jt.Interval(
                            lower_bound, seg.interval_bp.end):
                        upper_bound = entry.ibd_segment.interval.end
                    else:
                        upper_bound = seg.interval_bp.end

                    (xstart, xend) = draw_rect(
                            upper_bound, lower_bound, inside_ibd)

                    if upper_bound < seg.interval_bp.end:
                        inside_ibd = not inside_ibd
                        lower_bound = upper_bound
                        upper_bound = seg.interval_bp.end
                        draw_rect(upper_bound, lower_bound, inside_ibd)
    return im

def n_most(seq, n, comp=op.lt):
    """ Return the n smallest elements in the seq. If n is greater than or
        equal to len(seq), this is just selection sort, and better sorting
        algorithms should be used instead.
        """
    outseq = list(seq) # copy the input sequence
    def swap(s, a, b):
        t = s[b]
        s[b] = s[a]
        s[a] = t

    for i in xrange(min(n, len(seq))):
        v = outseq[i]
        for j in xrange(i + 1, len(outseq)):
            if comp(outseq[j], v):
                swap(outseq, i, j)
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

        # set to 10 for testing
        longest_matches = sorted(
                matches, key=lambda m: len(m.ibd_segment.interval))[-50:][::-1]

        map(lambda m: print(len(m.ibd_segment.interval)), longest_matches)

        # a sanity check
        for m in longest_matches:
            for individual in m.individuals:
                for (hcode, chromosomes) in individual.ancestries.items():
                    for chromosome in chromosomes:
                        for seg in chromosome.segments:
                            assert(len(seg.interval_bp) > 0)

        im = plot_matches(longest_matches)
        im.save("%d.jpg" % chromosome_number, "PNG")

if __name__ == "__main__":
    if len(args) != 2:
        raise ValueError("incorrect number of arguments")
    main(args[1])
