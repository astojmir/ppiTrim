#
# ===========================================================================
#
#                            PUBLIC DOMAIN NOTICE
#               National Center for Biotechnology Information
#
#  This software/database is a "United States Government Work" under the
#  terms of the United States Copyright Act.  It was written as part of
#  the author's official duties as a United States Government employee and
#  thus cannot be copyrighted.  This software/database is freely available
#  to the public for use. The National Library of Medicine and the U.S.
#  Government have not placed any restriction on its use or reproduction.
#
#  Although all reasonable efforts have been taken to ensure the accuracy
#  and reliability of the software and data, the NLM and the U.S.
#  Government do not and cannot warrant the performance or results that
#  may be obtained by using this software or data. The NLM and the U.S.
#  Government disclaim all warranties, express or implied, including
#  warranties of performance, merchantability or fitness for any particular
#  purpose.
#
#  Please cite the author in any work or product based on this material.
#
# ===========================================================================
#
# Code author:  Aleksandar Stojmirovic
#
"""
Routines for parsing iRefIndex PSI-MI TAB 2.6 files as well as ppiTrim
modifications.
"""

from interaction import Interaction


def full_mitab_iterator(input_fp):
    """
    Runs over all lines from the file.
    """
    input_fp.seek(0)
    line = input_fp.next()  # header
    j = 1          # line number
    i = len(line)  # line offset
    for line in input_fp:
        yield line, (i, j)
        i += len(line)
        j += 1


def partial_mitab_iterator(input_fp, line_offsets, line_numbers=None):
    """
    Runs only over lines with offsets in line_offsets.
    """
    if line_numbers is None:
        line_iterator = enumerate(line_offsets)
    else:
        line_iterator = zip(line_numbers, line_offsets)

    for j, i in line_iterator:
        input_fp.seek(i)
        yield input_fp.next(), (i, j)


def full_interaction_consumer(scanner):
    """
    Interaction scanner parse_modified_mitab_file yields a pair (interaction,
    is_complex). This consumer collects all interactions in two lists, pairs
    and complexes and returns them as a result.
    """
    pairs = []
    complexes = []
    for interaction, _ in scanner:
        if interaction.is_complex():
            complexes.append(interaction)
        else:
            pairs.append(interaction)
    return pairs, complexes


def offsets_by_pmid_consumer(scanner):
    """
    Constructs a map pmid -> list-of-offsets for subsequent processing
    """

    offsets_map = {}
    for interaction, lines in scanner:
        line_offsets = lines[0]
        key = interaction.pmid
        if key is not None:
            if key not in offsets_map:
                offsets_map[key] = []
            offsets_map[key].extend(line_offsets)
    return offsets_map


def parse_mitab_file(input_fp, line_iterator_func, line_iterator_args=None,
                     interaction_class=Interaction):
    """
    Parse raw lines from iRefIndex MITAB 2.6.
    Yield interactions and their line offsets and line numbers.
    """

    if line_iterator_args is None:
        line_iterator_args = tuple()

    # First pass - parse and yield only binary interactions, save locations of
    # complexes.
    complexes_map = {}
    line_iterator = line_iterator_func(input_fp, *line_iterator_args)
    for line, line_locs in line_iterator:
        i, j = line_locs  # line offset and line number
        cols = line.strip().split('\t')
        idata = interaction_class.get_complex_data(cols)
        if idata is None:
            # binary interaction
            interaction = interaction_class()
            interaction.from_cols(cols)
            if len(interaction.interactors):
                yield interaction, ([i], [j])
        else:
            # complex saved for later
            if idata not in complexes_map:
                complexes_map[idata] = [[], []]
            complexes_map[idata][0].append(i)
            complexes_map[idata][1].append(j)

    # Second pass - parse each complex separately
    for line_offsets, line_numbers in complexes_map.itervalues():
        line_iterator = partial_mitab_iterator(input_fp, line_offsets,
                                               line_numbers)
        interaction = interaction_class()
        for line, _ in line_iterator:
            cols = line.strip().split('\t')
            interaction.from_cols(cols)
        if len(interaction.interactors):
            yield interaction, (line_offsets, line_numbers)
