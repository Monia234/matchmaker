#!/usr/bin/env python

import ibd
import bed

class IBDAncestryMatch:
    """ For some pair of individuals, this class represents the product of
        their shared ancestry and their IBD segment. This class supports rich
        comparison, where instances are compared on the basis of the length of
        the IBD segment.
        """

    def __init__(self, ibd_segment, individual_1, individual_2):
        """ Constructor.

            Arguments:
                ibd_segment (ibd.IBDEntry):
                    the identity by descent (IBD) segment, as loaded from
                    (optionally filtered) GERMLINE output.
                individual_1, individual_2 (bed.Individual):
                    the individuals related in this IBD segment.
            """

        self.ibd_segment = ibd_segment
        self.
