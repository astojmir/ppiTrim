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
Assign PSI-MI interaction type terms to interactions from the BioGRID with the
type 'Biochemical Activity'.
"""

import sys
import os
import urllib
import urlparse
import zipfile
from download import download_file

ptm2itype = {
    'Acetylation': ('MI:0192', 'acetylation reaction'),
    'Deacetylation': ('MI:0197', 'deacetylation reaction'),
    'Demethylation': ('MI:0871', 'demethylation reaction'),
    'Dephosphorylation': ('MI:0203', 'dephosphorylation reaction'),
    'Desumoylation': ('MI:0568', 'desumoylation reaction'),
    'Deubiquitination': ('MI:0204', 'deubiquitination reaction'),
    'Glycosylation': ('MI:0559', 'glycosylation reaction'),
    'Methylation': ('MI:0213', 'methylation reaction'),
    'Nedd(Rub1)ylation': ('MI:0567', 'neddylation reaction'),
    'No Modification': ('MI:0414', 'enzymatic reaction'),
    'Phosphorylation': ('MI:0217', 'phosphorylation reaction'),
    'Prenylation': ('MI:0211', 'lipid addition'),
    'Proteolytic Processing': ('MI:0570', 'protein cleavage'),
    'Ribosylation': ('MI:0557', 'adp ribosylation reaction'),
    'Sumoylation': ('MI:0566', 'sumoylation reaction'),
    'Ubiquitination': ('MI:0220', 'ubiquitination reaction'),
}

# BioGRID TAB 2.0 FORMAT (October 2010)

# COLUMN  DESCRIPTION
# ------------------------------------------------
# 1       BioGRID Interaction ID.
# 2       Entrez Gene ID for Interactor A.
# 3       Entrez Gene ID for Interactor B.
# 4       BioGRID ID for Interactor A.
# 5       BioGRID ID for Interactor B.
# 6       Systematic name for Interactor A.
# 7       Systematic name for Interactor B.
# 8       Official symbol for Interactor A.
# 9       Official symbol for Interactor B.
# 10      Synonyms/Aliases for Interactor A.
# 11      Synonyms/Aliases for Interactor B.
# 12      Experimental System Name.
# 13      Experimental System Type.
# 14      First author surname
# 15      Pubmed ID
# 16      Organism ID for Interactor A
# 17      Organism ID for Interactor B.
# 18      Interaction Throughput.
# 19      Quantitative Score.
# 20      Post Translational Modification.
# 21      Phenotypes.
# 22      Qualifications.
# 23      Tags.
# 24      Source Database.


def _get_archived_filename(zf):
    # Here we assume that the archive contains exactly one file - the one that
    # we need, in TAB2 format. Later, if this is no longer the case, we will
    # need to look for it more carefully.
    return zf.namelist()[0]


def download_and_extract_biogrid(src, download_dir):
    """
    Download and extract a BioGRID dataset (zipped).
    """
    pathname = urllib.url2pathname(src)
    biogrid_zipfile = os.path.basename(urlparse.urlsplit(pathname)[2])
    zipfile_path = os.path.join(download_dir, biogrid_zipfile)

    if not os.path.exists(zipfile_path):
        print "Downloading latest BioGRID."
        download_file(src, zipfile_path)
    zf = zipfile.ZipFile(zipfile_path, 'r', allowZip64=True)
    biogrid_file = _get_archived_filename(zf)
    extracted_path = os.path.join(download_dir, biogrid_file)

    if not os.path.exists(extracted_path):
        print "Extracting latest BioGRID."
        out_fp = open(extracted_path, 'w')
        out_fp.write(zf.read(biogrid_file))
        out_fp.close()

    zf.close()
    return extracted_path


def biogrid_ptm_codes(tab2_file, output_file):
    """
    Retrieve the post-translational modifications associated with BioGRID's
    'biochemical activity' interactions.
    """
    if output_file == '-':
        out_fp = sys.stdout
    else:
        out_fp = open(output_file, 'w')

    fp = open(tab2_file, 'rU')
    for line in fp:
        # Skip any line starting with #
        if line[0] == '#':
            continue
        cols = line.split('\t')
        if cols[11] == 'Biochemical Activity':
            out_fp.write('%s\t' % cols[0])
            out_fp.write('%s\t%s\n' % ptm2itype[cols[19]])
    fp.close()
    if output_file != '-':
        out_fp.close()


if __name__ == "__main__":

    args = sys.argv[1:]
    biogrid_ptm_codes(args[0], args[1])
