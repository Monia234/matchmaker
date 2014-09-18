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

    # the length of the longest match scaled to image coordinates
    max_size = scale_x(len(matches[0]))

    # the x values before and after which nothing must be drawn
    min_x = conf.FIGURE_WIDTH / 2 - max_size / 2 - 1
    max_x = conf.FIGURE_WIDTH / 2 + max_size / 2 + 1

    # the interval (in image coordinates) outside which nothing must be drawn
    drawable_interval = jt.Interval(min_x, max_x)

    print("DRAWABLE INTERVAL: ", drawable_interval)

    # take the height of the figure, subtract by the number of margins (one
    # less than the number of matches) multiplied by the margin size, and
    # divide by the number of entries to plot.
    ENTRY_HEIGHT = (conf.FIGURE_HEIGHT -
            (len(matches) - 1) * conf.INTER_ENTRY_MARGIN) / float(len(matches))
    # simply two individuals in each match
    INDIVIDUAL_HEIGHT = jt.intround(ENTRY_HEIGHT / 2)

    for (i, entry) in enumerate(matches):
        print("ENTRY:", i)
        # the scaled width of this entry
        ibd_width_true = len(entry.ibd_segment.interval)
        ibd_width = scale_x(len(entry.ibd_segment.interval))
        (ibd_start_true, ibd_end_true) = entry.ibd_segment.interval.to_tuple()
        (ibd_start, ibd_end) = (scale_x(ibd_start_true), scale_x(ibd_end_true))

        print("\tIBD: ", "(", ibd_start_true, ", ", ibd_end_true, ") -> (",
                ibd_start, ", ", ibd_end, ") \n",
                "\tIBD WIDTH: ", ibd_width_true, " -> ", ibd_width, sep='')

        # where the starting point of the IBD segment needs to be drawn for it
        # to be centered
        ibd_x0    = conf.FIGURE_WIDTH / 2 - ibd_width / 2
        ibd_start_x = scale_x(entry.ibd_segment.interval.start)
        ibd_center_offset = ibd_x0 - ibd_start_x

        # used to properly align individuals
        last_y = None

        print("\tSEGMENT OFFSET:", ibd_center_offset)

        # the ibd center offset needs to be added to all the positions, so that
        # the ibd segment is properly centered in the figure.

        # the chromosome of this entry (rebound to save space)
        my_chr   = entry.chromosome
        for (indiv_i, individual) in enumerate(entry.individuals):
            print("\tINDIVIDUAL: ", indiv_i, " (", individual.name, ")",
                    sep='')

            # the y origin for this individual
            my_y0    = jt.intround(last_y or
                       i * (ENTRY_HEIGHT + conf.INTER_ENTRY_MARGIN)
                     + indiv_i * INDIVIDUAL_HEIGHT)

            # reset the last_x counter for each individual, unlike last_y which
            # only needs to be reset per entry.
            last_x = None

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
                # the drawing function also verifies that the region to draw is
                # within the drawable bounds. If part or all of it isn't, then
                # it will draw only the portion that is contained.
                def draw_rect(upper_bound, lower_bound, inside_ibd, last_x):
                    segment_width_true = upper_bound - lower_bound
                    segment_width = scale_x(segment_width_true)

                    if last_x is not None:
                        xstart = last_x
                    else:
                        xstart = jt.intround(scale_x(lower_bound) + ibd_center_offset)
                    xend   = jt.intround(xstart + segment_width)
                    drawn_interval = jt.Interval(xstart, xend)
                    intersection = drawn_interval.intersection(
                            drawable_interval)
                    print("\t\tDRAWING INTERSECTION:", intersection)

                    # Perform the draw, if there is something to draw
                    if intersection:
                        rect = Image.new("RGBA",
                                (jt.intround(len(intersection)),
                                    jt.intround(INDIVIDUAL_HEIGHT)))
                        rect_draw = ImageDraw.Draw(rect)
                        rect_draw.rectangle([0, 0, len(intersection),
                            INDIVIDUAL_HEIGHT],
                            fill=get_code_color(seg.code, inside_ibd))
                        im.paste(rect,
                                (jt.intround(intersection.start), my_y0),
                                mask=rect)
                    return (intersection.start, my_y0,
                            xend, my_y0 + INDIVIDUAL_HEIGHT)

                splitting = True # whether the segment will be split

                # identify the width of the segment, accounting for possible
                # IBD segment beginnings or ends within it.
                lower_bound = seg.interval_bp.start
                if entry.ibd_segment.interval.start in seg.interval_bp:
                    print("\t\tIBD START (SEGMENT ", j, ")", sep='')
                    upper_bound = entry.ibd_segment.interval.start
                elif entry.ibd_segment.interval.end in seg.interval_bp:
                    print("\t\tIBD END (SEGMENT ", j, ")", sep='')
                    upper_bound = entry.ibd_segment.interval.end
                else:
                    # the fact that upper_bound /= seg.interval_bp.end will be
                    # used to determine that a cut was made in this segment.
                    # The other portion of the segment will be drawn a bit
                    # later.
                    upper_bound = seg.interval_bp.end
                    splitting = False

                # draw the rectangle
                rect = draw_rect(upper_bound, lower_bound, inside_ibd, last_x)
                last_x = rect[2] # the x-coordinate of the right side
                last_y = rect[3] # set the last y to the bottom of this rect

                print("\t\tSEGMENT #", j, ": (", lower_bound, ", ",
                        upper_bound, ") -> DRAW [",", ".join(map(str, rect)),
                        "] ", seg.code.name, sep='')

                # splitting is False only if the upper bound goes to the end of
                # the segment. If splitting is True, then we are splitting the
                # segment.
                if splitting:
                    # this represents a change in whether we're inside the IBD
                    inside_ibd = not inside_ibd # so we invert the state
                    # save the value of upper_bound as the split point
                    split_point1 = upper_bound
                    # set out bounds to the second part of the segment
                    lower_bound = upper_bound + 1

                    # if the ibd segment is ending in the between the split
                    # point and the end of the segment, then we need to split
                    # again
                    if entry.ibd_segment.interval.end in jt.Interval(
                            lower_bound, seg.interval_bp.end):
                        upper_bound = entry.ibd_segment.interval.end
                        splitting = True
                    else:
                        upper_bound = seg.interval_bp.end
                        splitting = False

                    rect2 = draw_rect(upper_bound, lower_bound, inside_ibd,
                                      last_x)
                    last_x = rect2[2]

                    print("\t\tSEGMENT #", j, ".b: (", lower_bound, ", ",
                            upper_bound, ") -> DRAW [",", ".join(map(str, rect2)),
                            "] ", seg.code.name, sep='')

                    if splitting:
                        inside_ibd = not inside_ibd
                        lower_bound = upper_bound
                        upper_bound = seg.interval_bp.end
                        rect3 = draw_rect(upper_bound, lower_bound, inside_ibd,
                                          last_x)
                        last_x = rect3[2]

                        print("\t\tSEGMENT #", j, ".c: (", lower_bound, ", ",
                                upper_bound, ") -> DRAW [",", ".join(map(str, rect3)),
                                "] ", seg.code.name, sep='')
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

def main(project_dir, plot_name):
    global PROJECT_DIR
    PROJECT_DIR = project_dir

    IBD_DIR = path.join(PROJECT_DIR,
            "baharian_projects/HRS/data/dbGaP/fromRFMix/AfrAm/3_GERMLINE")

    BED_DIR = path.join(PROJECT_DIR,
            "barakatt_projects/HRS/results/HRS_AFRAM_20140609/outbed")

    make_ibd_path = lambda i: path.join(IBD_DIR,
            "".join(["chr", str(i), ".IBD.match.gz"]))

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

        # a sanity check
        for m in longest_matches:
            assert(m.individuals)
            for individual in m.individuals:
                assert(individual.ancestries)
                for (hcode, chromosomes) in individual.ancestries.items():
                    assert(chromosomes)
                    for chromosome in chromosomes:
                        assert(chromosome.segments)
                        for seg in chromosome.segments:
                            assert(len(seg.interval_bp) > 0)

        im = plot_matches(longest_matches)
        im.save("%s-%d.png" % (plot_name, chromosome_number), "PNG")

if __name__ == "__main__":
    if len(args) != 3:
        raise ValueError("incorrect number of arguments")
    main(args[1], args[2])
