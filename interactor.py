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
Classes for storing interactor data from iRefIndex (and ppiTrim).
"""

from fields import SingleID
from fields import MultipleIDs
from fields import TaxonomyID
from fields import PSIMITerm
from fields import EmptyField
from fields import MultiplePSIMITerms


class Interactor(object):
    """
    General interactor parsed from PSI-MI 2.6 file.
    """

    _fields = [('uid', 0, SingleID),
               ('alt', 2, MultipleIDs),
               ('alias', 4, MultipleIDs),
               ('organism', 9, TaxonomyID),
               ('biological_role', 16, MultiplePSIMITerms),
               ('experimental_role', 18, MultiplePSIMITerms),
               ('interactor_type', 20, MultiplePSIMITerms),
               ('xrefs', 22, EmptyField),
               ('annotations', 25, EmptyField),
               ('checksum', 32, SingleID)]

    def __init__(self):

        for field_name, _, field_class in self._fields:
            setattr(self, field_name, field_class())

    def __str__(self):

        data = ['%s\t%s\n' % (field_name, str(getattr(self, field_name)))
                for field_name, _, _ in self._fields]
        return ''.join(data)

    def from_cols(self, cols, offset):
        """
        Sets all relevant field names from a list of column strings and
        interactor offset (0 or 1)
        """

        for field_name, i, _ in self._fields:
            if i < len(cols):
                field = getattr(self, field_name)
                field.from_string(cols[i + offset])

    def to_cols(self, cols, offset):
        """
        Writes relevant field names to a list of column strings, given an
        interactor offset (0 or 1). cols must be a list large enough so that
        its interactor items can be replaced with data from this instance.
        """

        for field_name, i, _ in self._fields:
            if i + offset < len(cols):
                field = getattr(self, field_name)
                cols[i + offset] = field.__str__()

    def is_complex(self):
        """
        Return True if the interactor represents a complex.
        """

        if self.uid.db == 'complex':
            return True
        type_terms = set(item.term_id for item in self.interactor_type)
        if 'MI:0315' in type_terms:
            return True
        return False

    def copy_ids(self, other):
        """
        Copy (by reference) ID fields from other.
        """
        assert False
        for field_name in ['uid', 'alt', 'alias', 'organism']:
            setattr(self, field_name, getattr(other, field_name))


class iRefIndexInteractor(Interactor):
    """
    Interactor parsed from iRefIndex.
    """

    _fields = [('uid', 0, SingleID),
               ('alt', 2, MultipleIDs),
               ('alias', 4, MultipleIDs),
               ('organism', 9, TaxonomyID),
               ('biological_role', 16, PSIMITerm),
               ('experimental_role', 18, PSIMITerm),
               ('interactor_type', 20, PSIMITerm),
               ('xrefs', 22, EmptyField),
               ('annotations', 25, EmptyField),
               ('checksum', 32, SingleID),
               ('final_reference', 38, SingleID)]

    def __init__(self):

        super(iRefIndexInteractor, self).__init__()
        self.id = None

    def from_cols(self, cols, offset):

        super(iRefIndexInteractor, self).from_cols(cols, offset)
        assert all((self.checksum.db is not None,
                    self.checksum.acc is not None))
        self.id = self.checksum.acc  # rogid

    def is_complex(self):

        if self.uid.db == 'complex' or \
           self.interactor_type.term_id == 'MI:0315':  # (protein complex):
            return True
        return False

    def normalize_ids_with_gene(self, gene_ids, gene_symbols):
        """
        Set primary ID (uid) to NCBI Gene. Also keep only Gene IDs and rogid
        as secondary IDs.
        """

        # Reset the ID fields
        if self.is_complex():
            # We do not touch uid for a complex
            reset_fields = ['alt', 'alias']
        else:
            reset_fields = ['uid', 'alt', 'alias']

        for field_name, _, field_class in self._fields:
            if field_name in reset_fields:
                setattr(self, field_name, field_class())

        # Reassign IDs - only if not complex
        if not self.is_complex():
            self.uid.set('entrezgene/locuslink', '%d' % min(gene_ids))
            for gene_id, gene_symbol in zip(gene_ids, gene_symbols):
                new_id_field = SingleID()
                new_symbol_field = SingleID()
                new_id_field.set('entrezgene/locuslink', '%d' % gene_id)
                new_symbol_field.set('entrezgene/locuslink', gene_symbol)
                self.alt.ids.append(new_id_field)
                self.alias.ids.append(new_symbol_field)

        # Add rogid for all interactors (including complexes) - if it exists
        rogid_field = SingleID()
        rogid_field.set('rogid', self.id)
        self.alias.ids.append(rogid_field)
