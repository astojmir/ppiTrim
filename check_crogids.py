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
Check correspondence between iRefIndex cannonical rogids and canonical groups
obtained by ppiTrim.
"""

import sys
import getopt
import os
from collections import defaultdict

def get_mappings(id_mapping_file):

    cg2crogid = defaultdict(set)
    crogid2cg = defaultdict(set)
    with open(id_mapping_file, 'rU') as fp:

        for i in xrange(3):
            line = fp.next()
            assert line[0] == '#'
        line = fp.next()
        while line[0] != '#':
            fields = line.strip().split('\t')
            cg = int(fields[0])
            crogid = int(fields[5])
            cg2crogid[cg].add(crogid)
            crogid2cg[crogid].add(cg)
            line = fp.next()
    return cg2crogid, crogid2cg


def write_mismatches(cg2crogid, crogid2cg, fp=sys.stdout):

    items = [('ONE-TO-MANY ppiTrim GROUPS', cg2crogid),
             ('ONE-TO-MANY iRefIndex CROGIDs', crogid2cg),
             ]

    for msg, mapping in items:
        fp.write('#\n')
        fp.write('# %s\n' % msg)
        fp.write('#\n')

        for id1 in sorted(mapping.iterkeys()):
            ids2 = mapping[id1]
            if len(ids2) > 1:
                fp.write('%d\t-\t' % id1)
                fp.write('%s\n' % '\t'.join(map(str, sorted(ids2))))



if __name__ == "__main__":

    opts, args = getopt.getopt(sys.argv[1:], 'h', ['help'])
    for o, a in opts:
        if o in ("-h", "--help"):
            sys.exit()

    id_mapping_file = args[0]
    cg2crogid, crogid2cg = get_mappings(id_mapping_file)
    write_mismatches(cg2crogid, crogid2cg)
