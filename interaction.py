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
Classes for storing interaction data from iRefIndex (and ppiTrim).
"""

import hashlib
import base64
from fields import SingleID
from fields import MultipleIDs
from fields import TaxonomyID
from fields import PSIMITerm
from fields import MultiplePSIMITerms
from fields import EmptyField
from fields import PlainText
from fields import ConfidenceScores
from fields import DateField
from interactor import Interactor
from interactor import iRefIndexInteractor

HEADER_FIELDS = \
['uidA',
 'uidB',
 'altA',
 'altB',
 'aliasA',
 'aliasB',
 'method',
 'author',
 'pmids',
 'taxa',
 'taxb',
 'interactionType',
 'sourcedb',
 'interactionIdentifier',
 'confidence',
 'expansion',
 'biological_role_A',
 'biological_role_B',
 'experimental_role_A',
 'experimental_role_B',
 'interactor_type_A',
 'interactor_type_B',
 'xrefs_A',
 'xrefs_B',
 'xrefs_Interaction',
 'Annotations_A',
 'Annotations_B',
 'Annotations_Interaction',
 'Host_organism_taxid',
 'parameters_Interaction',
 'Creation_date',
 'Update_date',
 'Checksum_A',
 'Checksum_B',
 'Checksum_Interaction'
 ]


class Interaction(object):
    """
    General interaction parsed from PSI-MI 2.6 TAB file
    """

    _term_cache = {}

    _fields = [('detection_method', 6, PSIMITerm),
               ('author', 7, PlainText),
               ('publications', 8, MultipleIDs),
               ('interaction_type', 11, PSIMITerm),
               ('source_db', 12, MultiplePSIMITerms),
               ('ids', 13, MultipleIDs),
               ('confidence', 14, ConfidenceScores),
               ('expansion', 15, PlainText),
               ('xrefs', 24, EmptyField),
               ('annotations', 27, EmptyField),
               ('host_taxid', 28, TaxonomyID),
               ('parameters', 29, EmptyField),
               ('creation_date', 30, DateField),
               ('update_date', 31, DateField),
               ('checksum', 34, SingleID),
               ('negative', 35, PlainText),
               ]
    _interactor_class = Interactor

    @classmethod
    def get_complex_data(cls, cols):
        """
        Check if the cols come from a line describing a complex. If yes, return
        all identifying fields, otherwise return None.
        """
        p = cls._interactor_class()
        p.from_cols(cols, 0)
        if p.is_complex():
            idata = tuple([p.uid.acc] + [cols[i] for _, i, _ in cls._fields])
        else:
            idata = None
        return idata

    @staticmethod
    def write_header(fp):
        """
        Write the header with field names (36 columns)
        """
        fp.write('#')
        fp.write('\t'.join(HEADER_FIELDS))
        fp.write('\n')

    def __init__(self):

        self.complex = None
        self.idata = None
        self.num_lines = 0

        self.interactors = []
        for field_name, _, field_class in self._fields:
            setattr(self, field_name, field_class())

        # Extra fields related to iRefIndex
        self._edgetype = None
        self.pmid = None
        self.rigid = None

        # Used only for ppiTrim complexes. Bait is a distinguished point so two
        # complexes with same members, properties etc. are not the same if
        # their baits are different. Just in case, we allow more than one bait,
        # although this should not happen for complexes.
        self.baits = None

    def __str__(self):

        data = ['%s\t%s\n' % (field_name, str(getattr(self, field_name)))
                for field_name, _, _ in self._fields]
        return ''.join(data)

    def find_baits(self):

        baits = []
        # Template-deflated complexes have the bait of their template annotated
        bait_ids = None
        for item in self.confidence.ids:
            if item.db == 'templatebaits':
                bait_ids = item.acc.split(',')
                break

        if bait_ids is not None:
            for p in self.interactors:
                if p.uid.acc in bait_ids:
                    baits.append(p)
            assert len(baits) > 0
        else:
            # Normal case - just look for the bait
            for p in self.interactors:
                for erole in p.experimental_role:
                    if erole.term_id == 'MI:0496':
                        baits.append(p)
                        # Here we force just one bait but this could be removed
                        # later
                        break
        return baits

    def set_template_baits(self, bait_ids):

        templatebaits = SingleID()
        bait_ids = list(bait_ids)
        templatebaits.set('templatebaits', ','.join(bait_ids))
        self.confidence.ids += [templatebaits]

    def get_ppiTrim_id(self):
        """
        Compute ppiTrim_id (hash) for the interaction.
        Also assign the bait attribute.
        """

        # To assign the bait(s) we either need to do it when each interactor is
        # added (not possible now) or wait until all interactors are input but
        # before this information is used anywhere. This function is where the
        # bait is used first (to compute SHA1 digest) so it makes sense to set
        # the bait right here. Also, this function is one of the last steps in
        # construction of ppiTrim consolidated interactions.
        self.baits = self.find_baits()

        hash_factory = hashlib.sha1()
        uids = [p.uid.acc for p in self.interactors]
        baits = [p.uid.acc for p in self.baits]
        pubs = [pub.__str__() for pub in self.publications.ids]
        props = map(str, [item.name for item in self.source_db] +
                         [self.detection_method.term_id,
                          self.interaction_type.term_id])
        msg = '.'.join(pubs + uids + props + [self.edgetype] + baits)
        hash_factory.update(msg)
        return base64.b64encode(hash_factory.digest())

    def is_complex(self):
        """
        Return True if interaction is considered a 'complex'
        """
        if self.complex is None:
            return False
        return True

    def _get_source_ids(self):
        srcids = [item.get() for item in self.ids.ids
                  if item.db not in ('rigid', 'irigid', 'edgetype', 'ppiTrim')]
        return srcids

    def _get_irigids(self):
        irigids = [item.get() for item in self.ids.ids
                   if item.db == 'irigid']
        return irigids

    def _set_extra_fields(self):

        for item in self.publications.ids:
            if item.db == 'pubmed':
                self.pmid = int(item.acc)
                break

        for item in self.ids.ids:
            if item.db == 'edgetype':
                self._edgetype = item
            elif item.db == 'rigid' and self.rigid is None:
                self.rigid = item.acc

    def from_cols(self, cols):
        """
        Sets all relevant field names and interactors from a list of column
        strings.
        """

        if not self.is_complex():
            # First time - add everything
            p1 = self._interactor_class()
            p1.from_cols(cols, 0)
            if p1.is_complex():
                self.complex = p1
                self.idata = tuple([cols[i] for _, i, _ in self._fields])
            else:
                self.interactors.append(p1)

            for field_name, i, _ in self._fields:
                field = getattr(self, field_name)
                field.from_string(cols[i])
            self._set_extra_fields()

        # Subsequently - only if complex. The calling routine should ensure
        # that the line being added is the part of the same complex as the
        # first one.
        p2 = self._interactor_class()
        p2.from_cols(cols, 1)
        self.interactors.append(p2)
        self.num_lines += 1

    def to_file(self, fp):
        """
        Write entire interaction in PSI-MI 2.6 TAB format (one or more lines)
        """

        cols = [None] * 36
        for field_name, i, _ in self._fields:
            field = getattr(self, field_name)
            cols[i] = field.__str__()

        if self.is_complex():
            self.complex.to_cols(cols, 0)
            for p in self.interactors:
                p.to_cols(cols, 1)
                fp.write('\t'.join(cols))
                fp.write('\n')
        else:
            # There are exactly two interactors
            for j, p in enumerate(self.interactors):
                p.to_cols(cols, j)
            fp.write('\t'.join(cols))
            fp.write('\n')

    def binary_from_complex(self, p1, p2, code):
        """
        Extract a binary interaction given by p1 and p2 from a complex
        """

        assert self.is_complex()
        interaction = Interaction()
        for field_name, _, _ in self._fields:
            prop1 = getattr(self, field_name)
            prop2 = getattr(interaction, field_name)
            prop2.set(*prop1.get())

        interaction.expansion.txt = 'none'
        interaction.creation_date.set_today()
        interaction.update_date.set_today()
        interaction.edgetype = code

        interaction.interactors = [p1, p2]
        return interaction

    def _get_edgetype(self):
        if self._edgetype is not None:
            return self._edgetype.acc
        return None

    def _set_edgetype(self, code):
        if self._edgetype is not None:
            self._edgetype.acc = code
        else:
            self._edgetype = SingleID()
            self._edgetype.set('edgetype', code)
            self.ids.ids.append(self._edgetype)

    edgetype = property(_get_edgetype, _set_edgetype)


class iRefIndexInteraction(Interaction):
    """
    Interaction parsed from iRefIndex PSI-MI 2.6 TAB file
    """
    _fields = [('detection_method', 6, PSIMITerm),
               ('author', 7, PlainText),
               ('publications', 8, MultipleIDs),
               ('interaction_type', 11, PSIMITerm),
               ('source_db', 12, PSIMITerm),
               ('ids', 13, MultipleIDs),
               ('confidence', 14, ConfidenceScores),
               ('expansion', 15, PlainText),
               ('xrefs', 24, EmptyField),
               ('annotations', 27, EmptyField),
               ('host_taxid', 28, TaxonomyID),
               ('parameters', 29, EmptyField),
               ('creation_date', 30, DateField),
               ('update_date', 31, DateField),
               ('checksum', 34, SingleID),
               ('negative', 35, PlainText),
               ]
    _interactor_class = iRefIndexInteractor
