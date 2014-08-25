#!/usr/bin/env python

import ibd
import bed
import jerrington_tools as je

from itertools import imap
from operator import eq

ANCESTRY_DATA_CACHE = {}
def CLEAR_ANCESTRY_CACHE():
    ANCESTRY_DATA_CACHE = {}

class IBDAncestryMatch:
    """ For some pair of individuals, this class represents the product of
        their shared ancestry and their IBD segment. This class supports rich
        comparison, where instances are compared on the basis of the length of
        the IBD segment.
        """

    @staticmethod
    def from_ibd_segment(ibd_segment, bed_dir, generate=False, robust=False):
        """ Construct an IBDAncestryMatch from only an IBD segment and a
            repository of bed files, optionally generating the shared ancestry
            segment, and optionally checking for robustness.
            In other words, this factory method will do the work of loading the
            necessary bed files, using the default name decoration procedure
            (see bed.Individual) from the given directory.

            Arguments:
                ibd_segment (ibd.IBDEntry):
                    The segment from which the names of the related individuals
                    are used to load the relevant local ancestry data.
                bed_dir (path):
                    The path to the directory in which reside all the bed
                    files.
                generate (boolean) (default: False):
                    Whether to call the `compute` method on the newly-created
                    IBDAncestryMatch instance, which would generate the shared
                    ancestry segment.
                robust (boolean) (default: False):
                    Whether to perform robustness checking during the
                    computation of the shared ancestry segment. The value of
                    this parameter is irrelevant if `generate` is False.

            Returns (IBDAncestryMatch):
                The constructed value, with the shared ancestry segment having
                been computed if `generate` was set to True.

            Note:
                If in the dataset being used, the same individual is IBD with
                more than one other individual, then this factory method will
                load the bed file more than once, which will result in poor
                memory performance. For better performance, it is advisable to
                set the `cache` parameter to True, which will cause
                newly-loaded ancestry data to be cached, and recalled later if
                necessary.
            """
        individuals = map(je.curry2(bed.Individual.from_dir_and_name)(bed_dir),
                ibd_segment.name)
        match = IBDAncestryMatch(ibd_segment, individuals)
        if generate:
            match.compute(robust=robust)
        return match

    @staticmethod
    def generate(ibd_segment, individuals, robust=False):
        """ Factory method to construct an IBDAncestryMatch and immediately
            perform the shared ancestry determination, with optional
            robustness checking.
            The arguments are as given to `__init__` and `compute`, so it is
            advised to consult the documentation on those methods.
            """
        match = IBDAncestryMatch(ibd_segment, individuals)
        match.compute(robust=robust)
        return match

    def __init__(self, ibd_segment, individuals):
        """ Constructor.

            Arguments:
                ibd_segment (ibd.IBDEntry):
                    the identity by descent (IBD) segment, as loaded from
                    (optionally filtered) GERMLINE output.
                individuals (list of bed.Individual):
                    the individuals related in this IBD segment.
                    Warning: the order in which the individuals is given in the
                    argument list must match the order in which their names
                    appear in the IBD entry.

            Note:
                The constructor does not automatically generate the
                shared_segment member. To generate it, call the `compute`
                method on an instance, or use the static factory method
                `generate`, which will construct the instance and immediately
                call compute on it.
            """

        self.ibd_segment = ibd_segment
        self.individuals = individuals

        # TODO increase robustness by automatically switching the individuals
        # if it is determined that the IBD order is the opposite of the given
        # order
        if not all(imap(eq,
                je.for_each(self.individuals, je.project_c("name"))
                self.ibd_segment.name)):
            raise ValueError("the given IBD segment does not relate the given "
                    "individuals.")

    def compute(self, robust=False):
        """ Perform shared ancestry determination, with optional robustness
            checking.

            Argument:
                robust (boolean) (default: False):
                    Verify that the shared segment determination satisfies the
                    commutativity law by performing the determination both
                    ways and checking for the equality of the determined
                    segments. This effectively doubles the runtime of this
                    method, however.

            Returns (boolean):
                Whether the determined segment is nonempty.
            """
        self.shared_segment = self.individuals[0].shared_ancestry_with(individual2,
                Individual.bed_code_from_ibd(ibd_segment.haplotype[0]),
                Individual.bed_code_from_ibd(ibd_segment.haplotype[1]))

        return not self.is_empty()

    def is_empty(self):
        """ Wrapper for the inner Interval's is_empty method. """
        return self.shared_segment.is_empty()

    def __lt__(self, other):
        """ Compare this IBDAncestryMatch with another on the basis of their
            IBD segments. IBD segments are compared on the basis of their
            lengths.
            Therefore, the default ordering for IBDAncestryMatch instances is
            from smallest to largest.
            """
        return self.ibd_segment < other.ibd_segment