#!/usr/bin/env python

""" Dataset-specific utilities. """

# HRS

def is_afram_hrs(individual_list, entry):
    """ We are interested only in those IBD entries that relate two
        African American individuals from the HRS dataset. Those
        individuals are listed in INDIVS, so we check for membership in
        that set.
        A good way to use this function is to construct a function from it with
        the individual list already supplied by currying.

            import jerrington_tools as jt
            my_is_afram_hrs = jt.curry2(is_afram_hrs)(my_individual_list)

        Now my_is_afram_hrs is a unary function taking IBDEntry instances,
        determining if they relate HRS African Americans.
        """
    return all(map(lambda x: x in individual_list, entry.name))

# SCCS

def sccs_name_parser(filename):
    """ Parse an SCCS local ancestry (.bed) filename to extract the subject ID
        and haplotype identifier.
        """
    return (filename[1:10], filename[11])

def is_sccs(entry):
    """ We are interested only in those IBD entries that relate
        individuals from the SCCS dataset. These individuals have
        subject IDs starting with "GWAS_".
        """
    return all(map(lambda x: x.startswith("GWAS_"), entry.name))

