#!/usr/bin/env python

import operator as op
import jerrington_tools as je

_EQUALS = je.curry2(op.eq)

class IBDEntry:
    """ Represents one entry of IDB between two individuals. Represents exactly
        one line of GERMLINE data.  This class defines __str__ to produce
        GERMLINE data.
        IBDEntry instances support rich comparison, on the basis of the length
        of the inner interval.
        """

    @staticmethod
    def from_GERMLINE(path):
        """ Parse a GERMLINE output file (possibly gzipped) into a list of
            IDBEntry objects.
            """
        L = []
        with je.maybe_gzip_open(path) as fileH:
            for line in fileH:
                L.append(IBDEntry.from_string(line))
        return L

    @staticmethod
    def from_string(input_str):
        # Split up the input string according to the GERMLINE output format
        (fam1, id1hap, fam2, id2hap,
                chr, start, end, dat) = input_str.split(None, 7)
        id1, hap1 = id1hap.split('.')
        id2, hap2 = id2hap.split('.')

        return IBDEntry(chr, id1, id2, hap1, hap2, fam1, fam2, start, end, dat)

    def __init__(self, chr, name1, name2, hap1, hap2, fam1, fam2,
            start, end, dat):
        def maybeint(d):
            return d if type(d) == int else int(d)

        self.chromosome = maybeint(chr)
        self.name = (name1, name2)
        self.haplotype = map(maybeint, [hap1, hap2])
        self.family = (fam1, fam2)
        self.interval = je.Interval(maybeint(start), maybeint(end))
        self.dat = dat

    def is_involved(self, individual_name):
        return any(map(_EQUALS(individual_name), self.name))

    def to_string(self):
        """ Convenience function for readability. """
        return self.__str__()

    def __str__(self):
        """ Convert this IBD entry into a valid line of GERMLINE output. """
        return "{}\t{}.{}\t{}\t{}.{}\t{}\t{}\t{}\t{}" \
               .format(self.family[0], self.name[0], self.haplotype[0],
                        self.family[1], self.name[1], self.haplotype[1],
                        self.chromosome,
                        self.interval.start, self.interval.end,
                        self.dat)

    def __repr__(self):
        return "IBDEntry({}, '{}', '{}', '{}', '{}', '{}', '{}', {}, {}, '{}')" \
               .format(self.chromosome, self.name[0], self.name[1],
                       self.haplotype[0], self.haplotype[1], self.family[0],
                       self.family[1], self.interval.start, self.interval.end,
                       self.dat)

    def __lt__(self, other):
        """ Compare this IBDEntry with another, on the basis of the length of
            the inner interval.
            """
        return len(self.interval) < len(other.interval)
