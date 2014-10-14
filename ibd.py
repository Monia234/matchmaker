#!/usr/bin/env python

import operator as op
import jerrington_tools as jt

_EQUALS = jt.curry2(op.eq)

BASEPAIR    = "MB"
CENTIMORGAN = "cM"

class IBDEntry:
    """ Represents one entry of IDB between two individuals. Represents exactly
        one line of GERMLINE data.  This class defines __str__ to produce
        GERMLINE data.
        IBDEntry instances support rich comparison, on the basis of the length
        of the inner interval.
        """

    @staticmethod
    def from_GERMLINE(path_or_handle):
        """ Parse a GERMLINE output file (possibly gzipped) into a list of
            IDBEntry objects.

            Arguments:
                path_or_handle (string or file handle):
                    If the argument is a string, then a handle will be opened
                    and closed by this method. Otherwise, the caller is
                    expected to manage the resource.
            """
        needs_close = False
        if isinstance(path_or_handle, str):
            handle = jt.maybe_gzip_open(path_or_handle)
            needs_close = True
        else:
            handle = path_or_handle

        L = map(IBDEntry.from_string, handle)

        if needs_close:
            handle.close()
        return L

    @staticmethod
    def ifrom_GERMLINE(handle):
        """ Lazily parse GERMLINE output (possible gzipped) as a generator.
            Since this function is a generator, the caller must handle opening
            and closing the file"""
        for line in handle:
            yield IBDEntry.from_string(line)

    @staticmethod
    def from_string(input_str):
        """ Parse a string of one line of GERMLINE output as an IBDEntry. """
        (fam1, id1hap, fam2, id2hap, chr, start, end,
                a1, a2, a3, a4, type, dat) = input_str.split(None, 12)
        id1, hap1 = id1hap.split('.')
        id2, hap2 = id2hap.split('.')

        dat = '\t'.join([a1, a2, a3, a4, type, dat])

        return IBDEntry(
                chr, id1, id2, hap1, hap2, fam1, fam2, start, end, dat, type)

    def __init__(self, chr, name1, name2, hap1, hap2, fam1, fam2,
            start, end, dat, type=BASEPAIR):
        self.chromosome = numparse(chr)
        self.name = (name1, name2)
        self.haplotype = map(int, [hap1, hap2])
        self.family = (fam1, fam2)

        # parse as int or float if string, else don't parse.
        numparse = lambda x: (float if type == "cM" else int)(x) if isinstance(x, str) else x
        self.interval = jt.Interval(numparse(start), numparse(end))
        self.dat = dat

    def complement(self):
        """ Construct an IBD segment ranging over the same region of the
            genome, but with the haplotype identifiers switched for both
            individuals.
            """
        return IBDEntry(self.chromosome, self.name[0], self.name[1],
                1 - self.haplotype[0], 1 - self.haplotype[1],
                self.family[0], self.family[1],
                self.interval.start, self.interval.end, self.dat, self.type)

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

    def __len__(self):
        return len(self.interval)
