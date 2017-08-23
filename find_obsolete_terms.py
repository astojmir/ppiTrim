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
List all obsolete PSI-MI terms from a column from a tab-delimited file.
"""

import sys
import getopt
import re
import obo

_HELP = \
"""
List all obsolete PSI-MI terms from a column from a tab-delimited file.

SYNOPSIS:

    %s input_file obo_file column0
    %s -h|--help

ARGUMENT:

    input_file          A tab-delimited file
    obo_file            A file in OBO format holding the PSI-MI ontology
    column0             The column with PSI-MI terms (indexed from 0)

OPTIONS:

    -h, --help          Print this message and exit

OUTPUT:

    Outputs to stdout a list of terms that are obsolete according to the
    obo_file.
""" % (__file__, __file__)

_EMSG = "Insufficient arguments. Use '%s -h' to print help.\n" % __file__


def find_obsolete_terms(input_file, obo_file, column0):

    i = int(column0)
    cvterm_pattern = re.compile('(MI:[0-9]+)')

    def _get_term_id(s):
        mtch = cvterm_pattern.search(s)
        if mtch is not None:
            return mtch.group(1)
        return None

    obsolete_terms = set()
    obo_fp = open(obo_file, 'rU')
    ontology = obo.OBOntology(obo_fp)
    input_fp = open(input_file, 'rU')
    for line in input_fp:
        if line[0] == '#':
            continue
        cols = line.strip().split('\t')
        term_id = _get_term_id(cols[i])
        if term_id is not None:
            term = ontology.get_term(term_id)
            if hasattr(term, 'is_obsolete'):
                obsolete_terms.add('%s(%s)' % (term.term_id, term.name))
    input_fp.close()
    obo_fp.close()

    for s in sorted(obsolete_terms):
        sys.stdout.write("%s\n" % s)


if __name__ == "__main__":

    opts, args = getopt.getopt(sys.argv[1:], 'h', ['help'])
    for o, a in opts:
        if o in ("-h", "--help"):
            sys.stdout.write(_HELP)
            sys.exit()

    if len(args) < 3:
        sys.stderr.write(_EMSG)
        sys.exit()

    find_obsolete_terms(*args)
