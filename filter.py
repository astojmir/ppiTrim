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
Filtering of iRefIndex interactions.
"""


class NullFile(object):
    """
    Fake file object - does nothing.
    """
    def write(self, s):
        pass

    def close(self):
        pass


def _print_lines(logfile_fp, lines):
    _lines = ['%d' % k for k in lines]
    logfile_fp.write('   Removed line(s): %s.\n' % ', '.join(_lines))


def is_filtered(interaction, filtered_pmids, logfile_fp, lines,
                accepted_taxids=None):
    """
    Returns True if the raw interaction is to be removed.
    """

    # Remove anything from BIND (original) and OPHID. We remove BIND because
    # BIND_Translation is more up-to-date and OPHID because it containes
    # computationally derived interactions
    if interaction.source_db.name in ('bind', 'ophid'):
        logfile_fp.write('Removing old BIND and OPHID\n')
        _print_lines(logfile_fp, lines)
        return True

    # Filter genetic interactions - this may have to be adjusted later
    if interaction.detection_method.term_id == 'MI:0254':
        logfile_fp.write('Genetic interaction\n')
        _print_lines(logfile_fp, lines)
        return True

    # Parse pubmed id and check against pmids to be filtered
    if interaction.pmid in filtered_pmids:
        logfile_fp.write('Filtered pubmed id (%d)\n' % interaction.pmid)
        _print_lines(logfile_fp, lines)
        return True

    # Filtering by taxonomy
    taxids = set(p.organism.taxid for p in interaction.interactors)
    # Either explicitly accept only specified IDs ...
    if accepted_taxids is not None:
        if not (taxids <= accepted_taxids):
            logfile_fp.write('Forbidden taxonomy ID found in: %s\n' %
                             ', '.join(sorted(map(str, taxids))))
            _print_lines(logfile_fp, lines)
            return True

    return False
