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

import dataset_utils

import ibd_anc_plot_config as conf

from itertools import islice, imap, chain, ifilter
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
    shared_center = lambda m: (m.shared_segments[0].interval_bp.start +
            m.shared_segments[-1].interval_bp.end) / 2
    delta = lambda m: shared_center(m) - ibd_center(m)

    min_bp = min(imap(
        lambda m: m.shared_segments[0].interval_bp.start + delta(m), matches))
    max_bp = max(imap(
        lambda m: m.shared_segments[0].interval_bp.end + delta(m), matches))
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
                def draw_rect(upper_bound, lower_bound, inside_ibd, last_x):
                    segment_width_true = upper_bound - lower_bound
                    segment_width = scale_x(segment_width_true)
                    rect = Image.new("RGBA",
                            (jt.intround(segment_width),
                                jt.intround(INDIVIDUAL_HEIGHT)))
                    rect_draw = ImageDraw.Draw(rect)
                    # Draw the rectangle, getting the color form the anc code
                    rect_draw.rectangle([0, 0, segment_width,
                        INDIVIDUAL_HEIGHT],
                        fill=get_code_color(seg.code, inside_ibd))

                    if last_x is not None:
                        xstart = last_x
                    else:
                        xstart = jt.intround(scale_x(lower_bound) + ibd_center_offset)
                    xend   = jt.intround(xstart + segment_width)
                    im.paste(rect,
                            (xstart, my_y0),
                            mask=rect)
                    return xstart, my_y0, xend, my_y0 + INDIVIDUAL_HEIGHT

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
                                upper_bound, ") "
                                "-> DRAW [", ", ".join(map(str, rect3)), "] ",
                                seg.code.name, sep='')
        for shared_segment in entry.shared_segments:
            rect_height = ENTRY_HEIGHT

            segment_width_true = len(shared_segment.interval_bp)
            segment_width = scale_x(segment_width_true)

            xstart = jt.intround(
                    scale_x(shared_segment.interval_bp.start)
                    + ibd_center_offset)
            xend = jt.intround(xstart + segment_width)

            shared_rect = Image.new("RGBA",
                    (jt.intround(segment_width),
                        jt.intround(rect_height)),
                    color=(128, 128, 128, 128))

            sa_y0 = jt.intround(
                    i * (ENTRY_HEIGHT + conf.INTER_ENTRY_MARGIN))

            im.paste(shared_rect, (xstart, sa_y0), mask=shared_rect)
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

def main(ibd_paths, outbed_path, indiv_list_path, plot_name):

    # the set of individuals in the HRS AfrAm
    #INDIVS = set(jt.with_file(jt.map_c(lambda x: x[:-1]), indiv_list_path))

    # Define some filter functions for the merged data.

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
    flipcurry2 = j.compose(j.curry2, j.flip)

    # construct a function that takes an IBDEntry and generates the match
    # object from it.
    match_from_ibd_segment__ = match.IBDAncestryMatch.from_ibd_segment
    my_from_ibd_segment = j.supply(match_from_ibd_segment__,
            {"generate":True,
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

    # take the N longest matches
    # this is unnecessarily slow since it sorts the list in full prior to
    # selecting the N longest. Better would be to use a partial selection sort.
    # TODO: implement a partial selection sort for this.
    longest_matches = sorted(
            matches, key=lambda m: len(m.ibd_segment.interval))[-50:][::-1]


    # sanity check
    for m in longest_matches:
        # each match relates exactly two individuals
        assert(len(m.individuals) == 2)
        for individual in m.individuals:
            # each individual has at least one ancestry
            assert(individual.ancestries)
            for (hcode, chromosomes) in individual.ancestries.items():
                # each haplotype has some chromosomes
                assert(chromosomes)
                for chromosome in chromosomes:
                    # each chromosome has some segments
                    assert(chromosome.segments)
                    for seg in chromosome.segments:
                        # each segment has a nonzero length
                        assert(len(seg.interval_bp) > 0)

    # plot the matches
    im = plot_matches(longest_matches)
    im.save(plot_name, "PNG")

USAGE = "This is a usage message." # TODO: write real usage message

def die_with_usage():
    print(USAGE, file=sys.stderr)
    sys.exit(1)

if __name__ == "__main__":
    i = 1
    ibd_paths = []
    outbed_path = None
    indiv_list_path = None
    output_file_path = None
    while i < len(args):
        arg = args[i]
        nextarg = lambda: args[i+1]
        if arg == "--ibd":
            ibd_paths.append(nextarg())
            i += 1
        elif arg == "--bed":
            outbed_path = nextarg()
            i += 1
        elif arg == "--ids":
            indiv_list_path = nextarg()
            i += 1
        elif arg == "-o":
            output_file_path = nextarg()
            i += 1
        else:
            print("Unrecognized command line option ``%s''." % arg,
                    file=sys.stderr)
            die_with_usage()
        i += 1

    if not (ibd_paths and outbed_path and indiv_list_path and
            output_file_path):
        print("One or more necessary arguments have no value.", file=sys.stderr)
        die_with_usage()

    main(ibd_paths, outbed_path, indiv_list_path, output_file_path)
