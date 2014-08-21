#!/usr/bin/env python

from __future__ import print_function
from functools import total_ordering
from itertools import izip, repeat
import re

import jerrington_tools as je

# The colors used to generate the ancestry codes
class UnknownAncestryError(Exception):
    """ Represents a parse failure regarding ancestry names. """
    pass

class AncestryCode(object):
    COLOR_RED    = (1.0, 0,   0)
    COLOR_BLUE   = (0,   0,   1.0)
    COLOR_YELLOW = (0,   1.0, 1.0)

    CODENAMES = ["AFR", "EUR", "NAT"]

    """ This class represents an identifier for a given ancestry.
        Ancestries have an associated color used to render themselves when
        generating plots.
        Creating custom instances of this class is discouraged, in favour of
        using the factory methods for the ancestries that we are interested in.
        """
    def __init__(self, name, color):
        self.name = name
        self.color = color

    def __eq__(self, other):
        """ Compare two ancestry codes for equality.
            The sufficient condition for the equality of ancestry codes is the
            equality of their names. In other words, the associated color is
            irrelevant.
            """
        return self.name == other.name

    @staticmethod
    def make_EUR():
        """ Generate a European ancestry code. """
        return AncestryCode("EUR", AncestryCode.COLOR_RED)

    @staticmethod
    def make_AFR():
        """ Generate an African ancestry code. """
        return AncestryCode("AFR", AncestryCode.COLOR_BLUE)

    @staticmethod
    def make_NAT():
        """ Generate a Native American ancestry code. """
        return AncestryCode("NAT", AncestryCode.COLOR_YELLOW)

    @staticmethod
    def make_from_name(name):
        """ Generate an ancestry code from the codename of the ancestry.

            Law: AC.make_EUR() == AC.make_from_name(AC.make_EUR().name)
                 where AC = AncestryCode

            An UnknownAncestryError is raised if an invalid codename is given.
            """
        if name == "EUR":
            return make_EUR()
        elif name == "NAT":
            return make_NAT()
        elif name == "AFR":
            return make_AFR()
        else:
            raise UnknownAncestryError()

    @staticmethod
    def is_valid_codename(name):
        return any(name == x for x in AncestryCode.CODENAMES)

class IntervalOperationError(Exception):
    pass

class Interval(object):
    """ Represents a (closed) interval, providing the interface for a
        non-iterable container. The type of the underlying values is
        irrelevant, provided that values of that type are totally ordered and
        have a sensible implementation of subtraction.

        Two Interval objects are considered equal if they have the same start
        and end points. An Interval A is considered less than another Interval
        B if A's end point is less than B's start point. The converse applies
        to a greater-than comparison. If the Intervals overlap, then the
        ordering comparisons will raise IntervalOperationError.

        Intervals are said to be _adjacent_ if the start point of one is equal
        to the end point of the other, or vice-versa. In the case of adjacent
        intervals, they can be straightforwardly joined into a larger one.
        """

    def __init__(self, start, end):
        """ Construct an Interval with the given start and end points. """
        self.start = start
        self.end = end

    @staticmethod
    def zero():
        """ Return an interval with zero length. """
        return Interval(0, 0)

    @staticmethod
    def from_tuple(tup):
        """ Construct an interval object from a tuple. """
        start, end = tup
        return Interval(start, end)

    def is_empty(self):
        return len(self) == 0

    def to_tuple(self):
        """ Convert this interval into a tuple. """
        return (self.start, self.end)

    def is_adjacent_to(self, other):
        """ Decide whether this interval is adjacent to another one.
            Because intervals are closed, adjacency of intervals implies that
            they overlap. """
        return self.start == other.end or other.start == self.end

    def joined_to(self, other):
        """ Return a new Interval that is the union of this one and another
            one, if the two intervals are adjacent. If they are not adjacent,
            an IntervalOperationError is raised.
            """
        if self.is_adjacent_to(other):
            return Interval(
                    min(self.start, other.start), max(self.end, other.end))
        else:
            raise IntervalOperationError("Trying to join disjoint intervals.")

    def overlaps(self, other):
        return (self.start in other or self.end in other
             or other.start in self or other.end in self)

    def is_disjoint_with(self, other):
        """ The negation of `overlaps`. """
        return not self.overlaps(other)

    def gap_to(self, other):
        """ Construct a new interval that represents the space between this
            interval and another one. If the two intervals are overlapping, a
            zero-length interval is constructed.
            """
        return (Interval.zero() if self.overlaps(other)
                else (Interval(self.end, other.start) if self < other
                      else Interval(other.start, self.end)
                     )
               )

    def __contains__(self, value):
        """ Determine whether the given value is contained within the interval.
            """
        return self.start <= value and value <= self.end

    def __len__(self):
        """ Calculate the length of the interval. """
        return self.end - self.start

class AncestrySegment(object):
    """ A simple product type from AncestryCode and Interval, which represents
        the ancestry of an individual at a particular location in its genome.

        Two AncestrySegment objects are considered equal if they range over the
        same interval, have the same AncestryCode, and are on the same
        chromosome.

        An AncestrySegment A is considered less than an AncestrySegment B if:
            * A's chromosome number is less than B's chromosome number; or
            * A's Interval is less than B's interval.

        See the Inverval class's documentation on how Interval's are compared.

        The serialization and deserialization methods to_string and from_string
        obey the following law:

            AncestrySegment.from_string(as.to_string()) == as
        """

    def __init__(self, code, chromosome, interval_bp, interval_cm):
        """ Constructor.

            Arguments:
                code (AncestryCode):
                    the code to use for this segment. one of the factory
                    methods from the AncestryCode class.
                chromosome (int):
                    the number of the chromosome on which the interval is
                    located.
                interval_bp (Interval):
                    the interval along which this ancestry holds, in base
                    pairs.
                interval_bp (Interval):
                    the interval along which this ancestry holds, in
                    centimorgan.

            Notes:
                Two intervals are required such that a valid .bed file can
                potentially be generated from a large list of AncestrySegment
                objects.
            """
        self.interval_bp = interval_bp
        self.interval_cm = interval_cm
        self.code = code
        self.chromosome = int(chromosome)

    def to_string(self):
        """ Convert this AncestrySegment into a string that could stand be a
            line in a .bed file.
            """
        return ("%d\t%d\t%d\t%s\t%f\t%f"
                % (self.chromosome,
                   self.interval_bp.start, self.interval_bp.end,
                   self.code.name,
                   self.interval_cm.start, self.interval_cm.end))

    def __contains__(self, value):
        if isinstance(value, int):
            return value in self.interval_bp
        elif isinstance(value, float):
            return value in self.interval_cm
        else:
            raise TypeError("cannot determine whether value of type ``%s'' "
                    "is a member of an interval." % str(type(value)))

    def __lt__(self, other):
        """ Compare this segment to another. If this segment is on a lower
            chromosome, it is automatically considered less than the other.
            If it is ona higher chromosome, it is automatically considered not
            less than the other one. If the chromosome numbers are the same,
            then the intervals are compared according to the rules of interval
            comparison described in the Interval class.
            """
        return (self.chromosome < other.chromosome or
                (self.chromosome == other.chromosome and
                 self.interval_bp < other.interval_bp and
                 self.interval_cm < other.interval_cm))

    def __eq__(self, other):
        return (self.interval_cm == other.interval_cm and
                self.interval_bp == other.interval_bp and # for consistency
                self.code == other.code and
                self.chromosome == other.chromosome)

    @staticmethod
    def from_string(input_string):
        """ Parse a line of a .bed file into an AncestrySegment.
            If the parsing of fields to ints or floats fails, then a ValueError
            will be raised by those conversion functions. If the values for the
            intervals are not in the correct order, then a ValueError will be
            raised. If the ancestry codename read from the string is not
            recognized, an UnknownAncestryError will be raised. If the given
            string does not contain the correct number of fields (six), then a
            ValueError will be raised.
            """
        words = input_string.split()
        id = lambda x: x # identity function
        chromosome, start_bp, end_bp, ancestry, start_cm, end_cm = [
                f(x) for (f, x)
                     in izip([int, int, int, id, float, float], words)]
        return AncestrySegment(
                AncestryCode.make_from_name(ancestry),
                chromosome,
                Interval(start_bp, end_bp),
                Interval(start_cm, end_cm))

class Individual(object):
    """ The complete ancestry information of one individual.
        This class simply wraps around a dictionary relating haplotype codes
        ("A" or "B") to a list of chromosomes. Each chromosome is just a list
        of AncestrySegment objects. AncestrySegment is a simple product type
        over Interval, which describes over what portion of the chromosome the
        ancestry is valid, and AncestryCode, which determines what ancestry is
        present there.
        Static methods are provided by this class to facilitate parsing the
        bed files that describe the local ancestry of the individual.
        Methods are also provided for writing out bed files. """

    # Example bed file name: T161269150_A_cM.bed
    # capture #1 is the name of the individual
    # capture #2 is the haplotype that this bed file concerns itself with
    _BED_FILE_REGEX_STR = "T(\d+)_(A|B)_cM\.bed"
    BED_FILE_REGEX = re.compile(Individual._BED_FILE_REGEX_STR)

    @staticmethod
    def _ancestry_pre_to_string(ancestry_pre):
        """ Convert an ancestry_pre object, i.e. a list of 22 lists of
            AncestrySegment objects, into a string, spanning multiple lines,
            that would constitute a valid bed file.

            Law
            """

    @staticmethod
    def _decorate_name(individual_name):
        """ Decorate an individual name such that it would successfully be
            parsed by the default regex. First, we verify the name to check
            that it consists solely of digits. Then, we prepend a T, the
            haplotype codes, the centimorgan identifier, and the bed extension.

            e.g. 1236195 -> T1236195_A_cM.bed, T1236191_B_cM.bed

            The two decorated names are returned as a tuple.
            """
        underscore = "_".join
        haplotype_codes = ["A", "B"]

        return tuple(map(
                lambda code: underscore("T" + individual_name, code, "cM.bed"),
                haplotype_codes))

    @staticmethod
    def _ancestry_pre_from_lines(lines):
        ancestry_pre = list(repeat([], 22)) # create the chromosome buckets
        for segment in (AncestrySegment.from_string(line) for line in lines):
            # add that segment to the appropriate bucket.
            ancestry_pre[segment.chromosome].append(segment)
        return ancestry_pre

    @staticmethod
    def _id_data_from_filename(filename, regex=None):
        """ Parse the given filename with the given regular expression. This
            method may seem useless, since all it really does is match the
            string against the regex, but it goes a bit further than that.
            First, if no regex is given (or the value None is provided), then a
            default regex is used. The provided regex must specify two capture
            groups: the first must yield the name of the individual, and the
            second must yield the haplotype identifier.
            """
        if regex is None:
            regex = Individual._BED_FILE_REGEX_STR
            prog = Individual.BED_FILE_REGEX
        else:
            prog = re.compile(regex)

        result = prog.match(filename)

        if result is None:
            raise ValueError("the given filename could not be matched against"
                    " the regular expression for parsing bed file names,"
                    " namely ``%s''." % regex)

        name, haplotype = result.groups()

        return (name, haplotype)

    @staticmethod
    def from_files(path_to_bed_a, path_to_bed_b):
        """ Parse the two files, one for each haplotype, as well as their
            filenames for identity information. The names parsed from the
            filenames must match, and the haplotype codes must be different.
            """
        # check that the given paths aren't the same
        if path_to_bed_a == path_to_bed_b:
            raise ValueError("different files must be given, one for each"
                    " haplotype")

        I = Individual # to save space on other static method calls

        # we can aggregate the paths, since we don't need them separately
        # anymore
        paths = [path_to_bed_a, path_to_bed_b]
        # create a helper function since we'll be mapping over the paths a lot
        for_each_path = for_each_c(paths)

        # Individual.id_data_from_filename takes just the filename, but we have
        # full paths, so we construct a function that will parse a full path.
        id_data_from_path = je.compose(
                I._id_data_from_filename, path.basename)

        # Parse the paths to extract the name and haplotype information
        namehaps = for_each_path(id_data_from_path) # :: [(name, hap)]
        # transpose so that [(name, hap)] -> ([name], [hap])
        (name_a, name_b), haps = zip(*namehaps) # pull apart the names list

        # check that the names are the same
        if name_a != name_b:
            raise ValueError("the ancestry information to parse is from "
                    "different individuals.")

        # check that the haplotypes are different
        if haps[0] == haps[1]:
            raise ValueError("the haplotype information is the same.")
        # technically, the paths being different, but the names being the same
        # sort of implies that it's the haplotypes that are different, but to
        # avoid weird cases, we might as well check anyway.


        # ancestry_pre_from_lines needs a list of lines, but we have paths.
        # je.file_as_lines represents a file handle as a list of lines
        # je.with_file lets us run a handle-needing function on a path, by opening
        # the file for us:
        # je.with_file_c(file_as_lines) :: path -> lines
        # I._ancestry_pre_from_lines :: lines -> ancestry_pre
        # V :: path -> ancestry_pre
        load_ancestry_pre_from_path = je.compose(
                I._ancestry_pre_from_lines, je.with_file_c(file_as_lines))
        ancestry_pres = for_each_path(load_ancestry_pre_from_path)

        # the following will throw if there're any issues with either of the
        # ancestry_pre objects generated just above
        je.for_each(ancestry_pres, I.check_ancestry_pre)

        # construct the dict linking the haplotype codes to the relevant
        # chromosome list. Since we maintained the order throughout, by using
        # map, we can simply zip together these two lists. Furthermore, we're
        # guaranteed that they're the same size.
        ancestry_dict = dict(zip(haps, ancestry_pres))
        return Individual(name_a, ancestry_dict)

    @staticmethod
    def from_dir_and_name(bed_dir, name):
        return Individual.from_files(
                *map(lambda dname: os.path.join(bed_dir, dname),
                    Individual._decorate_name(name)))

    @staticmethod
    def check_ancestry_pre(ancestry_pre):
        """ Perform a simple sanity check on an ancestry_pre object.
            An ancestry_pre object is simply a list of 22 lists of
            AncestrySegment objects. They must all be from the same haplotype,
            although it is impossible to verify this, simply given the object.
            This check ensures that the following criteria are met:
                * for each chromosome bucket N, segments in that bucket have
                  chromosome number N.
                * for each segment S except the last one in each bucket, S is
                  less than segment S + 1.
            If a criterion is not met, then a ValueError is raised with an
            appropriate message.
            """
        for (i, segments_for_chromosome_i) in enumerate(ancestry):
            for (j, segment_j) in enumerate(segments_for_chromosome_i):
                if segment_j.chromosome - 1 != i:
                    raise ValueError("segment with chromosome number %d found"
                            " in bucket for chromosome %d"
                            % (segment_j.chromosome, i + 1))
                if (j + 1 < len(segments_for_chromosome_i)
                        and not segment_j < segments_for_chromosome_i[j + 1]):
                    raise ValueError("segments not sorted correctly.")
        return True

    def __init__(self, individual_name, ancestries):
        """ Construct an object that represents an individual, complete with
            their ancestry.

            An _ancestry object_ is a dictionary relating a haplotype code
            (usually "A" or "B") to an ancestry_pre object, which is a list of
            22 lists of AncestrySegments.
            """
        self.ancestries = ancestries

    def __getitem__(self, i):
        """ Simply, a wrapper around the inner dict's __getitem__ function. """
        return self.ancestries[i]
