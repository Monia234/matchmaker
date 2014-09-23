#!/usr/bin/env python

from __future__ import print_function
from functools  import total_ordering
from itertools  import izip, repeat, imap
from operator   import eq

import re
import os

import jerrington_tools as je

# The colors used to generate the ancestry codes
class UnknownAncestryError(Exception):
    """ Represents a parse failure regarding ancestry names. """
    pass

class AncestryCode(object):
    COLOR_RED    = (255, 0,   0)
    COLOR_BLUE   = (0,   0,   255)
    COLOR_YELLOW = (0,   255, 255)

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

    def __repr__(self):
        return "AncestryCode.from_name('%s')" % (self.name,)

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

    @classmethod
    def from_name(cls, name):
        """ Generate an ancestry code from the codename of the ancestry.

            Law: AC.make_EUR() == AC.from_name(AC.make_EUR().name)
                 where AC = AncestryCode

            An UnknownAncestryError is raised if an invalid codename is given.
            """
        if name == "EUR":
            return cls.make_EUR()
        elif name == "NAT":
            return cls.make_NAT()
        elif name == "AFR":
            return cls.make_AFR()
        else:
            raise UnknownAncestryError()

    @staticmethod
    def is_valid_codename(name):
        return any(name == x for x in AncestryCode.CODENAMES)

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

        i.e. they are inverses of each other.
        """

    # How far apart must common-ancestry segments be to be considered distinct
    # this is used in Individual.shared_ancestry_with
    DISTANCE_CUTOFF = 100000
    # What fraction of the ancestry must be shared to warrant merging the
    # segments into one
    SMOOTH_CUTOFF = 0.8

    def __init__(self, code, chromosome, interval_bp, interval_cm):
        """ Construct an AncestrySegment object.

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
        """ Convert this AncestrySegment into a string that could stand as a
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

    def __len__(self):
        """ Return the length of this segment, in basepairs. To get the length
            in centimorgan, it is necessary to use
                len(seg.interval_cm)
            where seg is the ancestry segment to measure.
            """
        return len(self.interval_bp)

    def __eq__(self, other):
        return (self.interval_cm == other.interval_cm and
                self.interval_bp == other.interval_bp and # for consistency
                self.code        == other.code        and
                self.chromosome  == other.chromosome)

    def __repr__(self):
        return "AncestrySegment(%s, %s, %s, %s)" % tuple(
                map(repr, [self.code, self.chromosome, self.interval_bp,
                           self.interval_cm]))

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
                AncestryCode.from_name(ancestry),
                chromosome,
                je.Interval(start_bp, end_bp),
                je.Interval(start_cm, end_cm))

class SwitchPoint(object):
    """ Represents a change in ancestry at a given location.
        SwitchPoint instances are equal if all their fields are equal, and
        SwitchPoint instances are compared on the basis of their positions.
        """

    @staticmethod
    def from_segment_pair(left, right):
        if left.chromosome != right.chromosome:
            raise ValueError("cannot construct SwitchPoint across "
                    "chromosomes.")
        position_bp = (left.interval_bp.end + right.interval_bp.start) / 2
        position_cm = (left.interval_cm.end + right.interval_cm.start) / 2
        return SwitchPoint(left.chromosome, position_bp, position_cm,
                left.code, right.code)

    def __init__(self, chromosome, position_bp, position_cm,
            source, destination):
        """ Constructor.

            Arguments:
                chromosome (int):
                    the number of the chromosome on which the switch occurs.
                position_bp (bp):
                    the position of the switch point, in base pairs.
                position_cm (float):
                    the position of the switch point, in centimorgan.
                source (AncestryCode):
                    the ancestry on the left of this switch point.
                destination (AncestryCode):
                    the ancestry on the right of htis switch point.

            Notes:
                It is possible to represent the local ancestry on a chromosome
                as a list of switch points together with the start and end
                positions on that chromosome. Chromosome.switch_points creates
                such a representation.
            """
        if source == destination:
            raise ValueError("invalid SwitchPoint: the source and destination "
                    "ancestries must be different.")

        self.chromosome  = chromosome
        self.position    = position
        self.source      = source
        self.desitnation = destination

class Chromosome(object):
    """ Represents the ancestry information along one chromosome.
        Essentially wraps around a list of AncestrySegment objects.
        """

    @staticmethod
    def check_ancestry_segments(chromosome, segments):
        """ Verify that segments are ordered and are all on the same
            chromosome.
            """
        for (i, segment) in enumerate(segments[:-1]):
            next_segment = segments[i + 1]
            if not segment < next_segment:
                raise ValueError("segments are not ordered")
            if segment.chromosome != chromosome:
                raise ValueError("not all segments are on the same chromosome")
        return True

    def __init__(self, number, segments):
        """ Constructor.

            Arguments:
                number (int):
                    the number of this chromosome.
                segments (list of AncestrySegment):
                    the list of ancestry segments along this chromosome.

            Note:
                the segments are checked for consistency: they must be
                non-overlapping, ordered, and all indexed to the same
                chromosome number, which must also match the number passed in.
            """
        if not segments:
            raise ValueError("a non-empty list of segments must be given to "
                    "construct a Chromosome.")
        if number != segments[0].chromosome:
            raise ValueError("one or more segments given do not match the "
                    "given chromosome number.")
        Chromosome.check_ancestry_segments(number, segments),
        self.segments = segments
        self.number   = number
        self.switches = None # the memoized result of as_switch_points
        self.len      = None # memoized

    def segment_index_of(self, position):
        for (i, segment) in enumerate(self.segments):
            if position in segment:
                return i
        return -1

    def __len__(self):
        """ Return the length of this chromosome, in basepairs, by adding up
            the lengths of all the segments. The result of this is cached.
            """
        if self.len is None:
            self.len = sum(map(len, self.segments))
        return self.len

    def __getitem__(self, index):
        """ Get the AncestrySegment object associated with the given position.
            Using a floating-point value will perform a check in centimorgan,
            and using an integer value will perform a check in basepairs.
            If you want to access the nth AncestrySegment from the underlying
            list, then you must do so directly, via a.segments[n], where a is a
            Chromosome instance.
            """
        hit = je.flip(filter)(self.segments,
                lambda s: index in s)
        if len(hit) > 1:
            raise ValueError("inconsistency: more than one AncestrySegment "
            "contains the given position.")
        elif not hit:
            raise IndexError("cannot find the given position in this "
                    "chromosome.")
        else:
            return hit[0]

    def __repr__(self):
        return "Chromosome(%s, %s)" % tuple(
                map(repr, [self.number, self.segments]))

    def as_switch_points(self):
        """ Represent the internal list of AncestrySegment objects as a list of
            ``switch-points'', which describes how the ancestry fluctuates
            along the chromosome.

            Since Chromosome objects should not be mutated, this function
            performs memoization of its return value: subsequent calls will run
            in constant time.
            """
        if self.switches is not None:
            return self.switches

        switches = []
        for (left, right) in je.ipairs(self.segments):
            switches.append(SwitchPoint(left, right))
        self.switches = switches # memoize
        return switches

_BED_FILE_REGEX_STR = "T(\d+)_(A|B)_cM\.bed"

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
        Methods are also provided for writing out bed files.
        """

    # Example bed file name: T161269150_A_cM.bed
    # capture #1 is the name of the individual
    # capture #2 is the haplotype that this bed file concerns itself with
    BED_FILE_REGEX = re.compile(_BED_FILE_REGEX_STR)

    # Codes used to refer to haplotypes in a bed file
    # The index corresponding to the code is the value used in an IBD file
    # to refer to that same haplotype. In other words, "A" is at index 0, and
    # "0" is the haplotype code used in the IBD match file.
    HAPLOTYPE_CODES = ["A", "B"]

    @staticmethod
    def _ancestry_pre_to_string(ancestry_pre):
        """ Convert an ancestry_pre object, i.e. a list of 22 lists of
            AncestrySegment objects, into a string, spanning multiple lines,
            that would constitute a valid bed file.

            Law: something
            """
        return NotImplemented # TODO

    @staticmethod
    def _decorate_name(individual_name):
        """ Decorate an individual name such that it would successfully be
            parsed by the default regex. First, we verify the name to check
            that it consists solely of digits. Then, we prepend a T, the
            haplotype codes, the centimorgan identifier, and the bed extension.

            e.g. 1236195 -> T1236195_A_cM.bed, T1236191_B_cM.bed

            The two decorated names are returned as a tuple.
            """
        underscore = je.unsplat("_".join)

        return tuple(map(
                lambda code: underscore("T" + str(individual_name), code,
                        "cM.bed"),
                Individual.HAPLOTYPE_CODES))

    @staticmethod
    def _ancestry_pre_from_lines(lines):
        ancestry_pre = map(lambda _: [], xrange(22))
        # each line is one ancestry segment
        for segment in imap(AncestrySegment.from_string, lines):
            # add that segment to the appropriate bucket.
            ancestry_pre[segment.chromosome - 1].append(segment)
        # each bucket is one chromosome, so we construct the Chromosome
        # instances

        chromosomes = map(
                lambda (i, v): Chromosome(i, v),
                enumerate(ancestry_pre, 1))
        return chromosomes

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
            regex = _BED_FILE_REGEX_STR
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
        for_each_path = je.for_each_c(paths)

        # Individual.id_data_from_filename takes just the filename, but we have
        # full paths, so we construct a function that will parse a full path.
        id_data_from_path = je.compose(
                I._id_data_from_filename, os.path.basename)


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
                I._ancestry_pre_from_lines, je.with_file_c(je.file_as_lines))
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
        for (i, segments_for_chromosome_i) in enumerate(ancestry_pre):
            for (j, segment_j) in enumerate(segments_for_chromosome_i):
                if segment_j.chromosome - 1 != i:
                    raise ValueError("segment with chromosome number %d found"
                            " in bucket for chromosome %d"
                            % (segment_j.chromosome, i + 1))
                if (j + 1 < len(segments_for_chromosome_i)
                        and not segment_j < segments_for_chromosome_i[j + 1]):
                    raise ValueError("segments not sorted correctly.")
        return True

    @staticmethod
    def bed_code_from_IBD(ibd_haplotype_code):
        """ bed files and IBD match files use different conventions for
            haplotypes. This method will convert the IBD match notation into
            the bed notation.
            Specifically, IBD match files use the strings "0" and "1" (dotted
            onto the individual name, e.g. 16128895.0), whereas bed files use
            the strings "A" and "B".

            Argument:
                ibd_haplotype_code (single-character string):
                    The code from an IBD match file, either "0" or "1", to
                    convert to the bed file notation.
            """
        # Since the ibd haplotype code is just a number, we can parse it and
        # use it to index a list of bed-style haplotype codes.
        return Individual.HAPLOTYPE_CODES[int(ibd_haplotype_code)]

    def __init__(self, individual_name, ancestries):
        """ Construct an object that represents an individual, complete with
            their ancestry.

            An _ancestry object_ is a dictionary relating a haplotype code
            (usually "A" or "B") to an ancestry_pre object, which is a list of
            22 lists of AncestrySegments.
            """
        self.name = individual_name
        self.ancestries = ancestries

    def __getitem__(self, i):
        """ Simply, a wrapper around the inner dict's __getitem__ function. """
        return self.ancestries[i]

    def __repr__(self):
        return "Individual(%s, %s)" % tuple(
                map(repr, (self.name, self.ancestries)))

    def shared_ancestry_with(self, other, haplo_self, haplo_other, chromosome):
        """ For each haplotype, determine an interval along which this
            Individual has the same ancestry as another Individual.

            Arguments:
                other (Individual):
                    The other Individual with whom the ancestry is compared.
                haplo_self (haplotype code):
                    Which haplotype, for this individual, should be examined.
                haplo_other (haplotype code):
                    Which haplotype, for the other individual, should be
                    examined.

            Note:
                Recall that haplotype codes are the strings 'A' or 'B', used as
                indices for the underlying ancestry information dict.

            Law: this computation is commutative.
                 If A and B are instances of Individual, then
                 A.shared_ancestry_with(B) == B.shared_ancestry_with(A)
            """
        # pairs will be a list of two tuples, each tuple is
        # (self's ancestry, other's ancestry) on the same haplotype.
        pairs = zip(self.ancestries.values(), other.ancestries.values())

        # To represent whether we are in a shared region.
        shared = False

        haplos = (self[haplo_self][chromosome],
                  other[haplo_other][chromosome])

        haplo_sizes = map(lambda hap: len(hap.segments), haplos)

        # haplo has the chromosome on the same haplotype for each individual.
        # haplo[0] refers to set of chromosomes in haplotype N of individual 0
        # (self), where N is the iteration number of this for loop.

        # take the greatest of the two starting points, since we can't
        # really compare where there's nothing there!
        start_position = max(haplos[0].segments[0].interval_bp.start,
                             haplos[1].segments[0].interval_bp.start)

        # take the smallest of the two end points, since there's no point
        # in trying to check past there.
        last_position = min(haplos[0].segments[-1].interval_bp.end,
                            haplos[1].segments[-1].interval_bp.end)

        position = start_position # we start at the beginning, of course

        # but which segments does "the beginning" correspond to?
        segment_counter = [haplos[0].segment_index_of(position),
                           haplos[1].segment_index_of(position)]

        # verify that the segments were found properly
        if any(map(lambda x: x == -1, segment_counter)):
            raise ValueError("could not match start position in one or "
                    "more of the haplotypes")

        regions = [] # we'll collect Interval instances here

        region_start = -1 # some dummy initial values for these
        region_end   = -1

        # until one of the two chromosomes ends
        while position < last_position:
            # determine what the ancestries are for this iteration
            my_anc    = haplos[0].segments[segment_counter[0]]
            other_anc = haplos[1].segments[segment_counter[1]]

            if shared: # if we are in a shared region
                if my_anc.code == other_anc.code: #if the codes match
                    pass # then we remain in the shared region
                else: #the codes don't match
                    shared = False # we exit the shared region
                    region_end = position # we mark the end position
                    regions.append(je.Interval(region_start, region_end))
            else: # we are not in a shared region
                if my_anc.code == other_anc.code: # if the codes match
                    shared = True # we enter the shared region
                    region_start = position # we mark the start position
                    if regions: # if there is at least one region so far
                        # this condition is the same as the one above, so
                        # if it fails, there is an inconsistency
                        if region_end == -1:
                            raise ValueError("inconsistency: there are "
                                    "shared regions, but the last end "
                                    "position is the dummy initial value")
                        # we check that the distance between the two
                        # regions meets the cutoff requirement. If it does,
                        # then we merge the last region with the one
                        # beginning at this position.
                        if (region_start - region_end <
                                AncestrySegment.DISTANCE_CUTOFF):
                            # we set the start to that of the last region
                            region_start = regions[-1].start
                            del regions[-1] # remove the last region
                        else: #i.e. the regions are distinct
                            pass #no big deal.
                    else: #there are no past regions
                        pass #doesn't matter.
                else: #the ancestry codes don't match
                    # doesn't matter since we aren't in a shared region now
                    pass # we remain unchanged

            # Now we need to figure out what to move the position to
            # We're going to look at the next segment for each chromosome
            # we're traversing, to see which next one begins the earliest.

            #first, we get the next indices for each haplo
            next_indices = map(je.succ, segment_counter)

            # We check whether we've reached the end of a chromosome.
            if any(imap(eq,
                    haplo_sizes,
                    next_indices)):
                # If so, then we need to exit the loop
                # We need to make sure that if we're in a shared region, then
                # we close it off.
                if shared:
                    shared = False
                    region_end = last_position
                    regions.append(je.Interval(region_start, region_end))
                break

            # then, we associate each of these indices with the relevant index
            # via ``enumerate''. We then sort according to the lowest start
            # position.
            bests = sorted(list(enumerate(next_indices)),
                    key=
                    lambda (i, j): haplos[i].segments[j].interval_bp.start)
            best = bests[0] # we take the smallest value

            # remember this is a tuple s.t. (the index of the haplo, the
            # index for the segment)
            haplo_index, segment_index = best # unpack the tuple

            # assign the value we got to the relevant haplo
            segment_counter[haplo_index] = segment_index

            # actually assign the new position now
            position = (haplos[haplo_index].
                        segments[segment_index].
                        interval_bp.start)
        # end of the while loop

        # now, ``regions'' has been populated. We need to smooth this list
        # so that there remains only one segment.
        def total_length(segments):
            # Interval defines __len__ to give its size, so we can write
            # the following fold.
            return sum(map(len, segments))

        shared_fraction = total_length(regions) / float(len(haplos[0]))

        #if shared_fraction >= AncestrySegment.SMOOTH_CUTOFF:
        #    # join all the regions together into the final region
        #    final_region = je.Interval(regions[0].start, regions[-1].end)
        #else:
        #    final_region = je.Interval.zero()

        if regions:
            final_region = je.Interval(regions[0].start, regions[-1].end)
        else:
            final_region = je.Interval.zero()

        # TODO think about making it such that it's the caller of this
        # method who must smooth the intervals, that way, maybe, it can do
        # something more sophisticated, or do something that involves all
        # the segments separately.
        return final_region
