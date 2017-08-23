#! /usr/bin/env python
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
A helper script to summarize complexes from ppiTrim.
"""

import sys
import getopt
import random
from parser import parse_mitab_file
from parser import full_mitab_iterator

_HELP = \
"""
Summarize complexes from a ppiTrim file.

SYNOPSIS:

    %s ppiTrim_file output_file_A output_file_N
    %s -h|--help

ARGUMENTS:

    ppiTrim_file        A PSI-MI TAB 2.6 file output by ppiTrim
    output_file_A
    output_file_N

OPTIONS:

    -h, --help          Print this message and exit

""" % (__file__, __file__)

_EMSG = "Insufficient arguments. Use '%s -h' to print help.\n" % __file__


def summarize_complexes(ppiTrim_file, output_file):
    """
    Output a summary of generated complexes of 'A' and 'N' type.
    These are randomly shuffled for sampling.
    """

    lines = []
    input_fp = open(ppiTrim_file, 'rU')
    output_fp = open(output_file, 'w')

    scanner = parse_mitab_file(input_fp, full_mitab_iterator)

    for intr, _ in scanner:
        edgetype = intr.edgetype
        n = len(intr.interactors)
        if not intr.is_complex() or n > 20:
            continue
        if ('C' in edgetype or 'G' in edgetype or 'R' in edgetype):
            continue

        if not ('A' in edgetype):
            continue

        geneids = ', '.join(p.uid.acc for p in intr.interactors)
        symbols = ', '.join(p.alias.ids[0].acc for p in intr.interactors)

        for item in intr.confidence.ids:
            if item.db == 'sourcedbs':
                sourcedbs = item.acc
                break

        line = '\t'.join([intr.complex.uid.acc,
                          '%d' % len(intr.interactors),
                          edgetype,
                          sourcedbs,
                          '%d' % intr.pmid,
                          symbols,
                          geneids,
                          '*****'
                          ])
        lines.append(line)

    random.shuffle(lines)
    for line in lines:
        output_fp.write(line)
        output_fp.write('\n')

    input_fp.close()
    output_fp.close()


if __name__ == "__main__":

    opts, args = getopt.getopt(sys.argv[1:], 'h', ['help'])
    for o, a in opts:
        if o in ("-h", "--help"):
            sys.stdout.write(_HELP)
            sys.exit()

    if len(args) < 2:
        sys.stderr.write(_EMSG)
        sys.exit()

    summarize_complexes(*args)
