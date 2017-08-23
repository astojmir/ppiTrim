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
ppiTrim's consolidated interaction.
"""

from fields import SingleID
from fields import PlainText
from interaction import Interactor
from interaction import Interaction


class ConsolidatedInteraction(Interaction):
    """
    Consolidated interaction from ppiTrim.
    """

    def __init__(self, group, ontology, edgetype):

        Interaction.__init__(self)
        self._conflicts = []

        self._consolidate_all(group, ontology, edgetype)

    def _consolidate_interactors(self, group):
        # Order of interactors cannot be important - biological/experimental
        # role should be annotated. Interactors are grouped by uid, so this
        # should be well-characterized.
        old_interactors = {}
        for intr in group:
            for p in intr.interactors:
                if p.uid.get() not in old_interactors:
                    old_interactors[p.uid.get()] = []
                old_interactors[p.uid.get()].append(p)

        self_edge = (len(old_interactors) == 1)
        new_interactors = []
        for pgroup in old_interactors.itervalues():
            q = Interactor()

            # uid is the same by construction (smallest geneid)
            q.uid.set(*pgroup[0].uid.get())

            # alt consists of one or more gene ids. If more than one,
            # order is important because the symbols should be in the same
            # order in the alias field. So, to do this properly, we need to
            # join gene ids and symbols, and do icrogids (in alias)
            # separately.
            icrogids = set(item.acc for p in pgroup for item in p.alias.ids
                           if item.db == 'icrogid')
            genes = set()
            for p in pgroup:
                for geneid_, symbol_ in zip(p.alt.ids, p.alias.ids):
                    genes.add((geneid_.acc, symbol_.acc))

            for geneid, symbol in sorted(genes):
                geneid_ = SingleID()
                geneid_.set(q.uid.db, geneid)
                symbol_ = SingleID()
                symbol_.set(q.uid.db, symbol)
                q.alt.ids.append(geneid_)
                q.alias.ids.append(symbol_)

            for icrogid in sorted(icrogids):
                icrogid_ = SingleID()
                icrogid_.set('icrogid', icrogid)
                q.alias.ids.append(icrogid_)

            # Organism ought to be the same everywhere due to filtering in
            # Phase 1 - NOT TRUE ANY MORE
            q.organism.set(*pgroup[0].organism.get())

            # biological role
            bio_roles = [item for p in pgroup for item in p.biological_role
                         if item.term_id is not None]
            if len(bio_roles):
                q.biological_role.from_items(bio_roles)
            else:
                q.biological_role.add('MI:0499', 'unspecified role')

            # experimental role - set if consistent, otherwise unspecified
            exp_roles = [item for p in pgroup for item in p.experimental_role
                         if item.term_id is not None]
            if len(exp_roles):
                q.experimental_role.from_items(exp_roles)
            else:
                q.experimental_role.add('MI:0499', 'unspecified role')

            # interactor type is always MI:0326(protein) for PPIs
            q.interactor_type.add('MI:0326', 'protein')

            # xrefs, annotations and checksum are left as null
            new_interactors.append(q)

        new_interactors.sort(key=lambda p: p.uid.acc)

        if self_edge:
            # Self-edge
            q = new_interactors[0]
            new_interactors.append(q)
        self.interactors = new_interactors

    def _consolidate_complex_node(self, group):
        complex_nodes = [intr.complex for intr in group]
        if complex_nodes[0] is not None:
            self.complex = Interactor()
            self.complex.uid.set('complex', self.ppiTrim_id)
            self.complex.organism = complex_nodes[0].organism
            self.complex.biological_role.add('MI:0499', 'unspecified role')
            self.complex.experimental_role.add('MI:0499', 'unspecified role')
            self.complex.interactor_type.add('MI:0315', 'protein complex')

    def _consolidate_detection_method(self, group, ontology):
        dterms = [intr.detection_method.term_id for intr in group
                  if intr.detection_method.term_id is not None]
        if len(dterms):
            term_id = ontology.finest_consistent_term(dterms, self._term_cache)
        else:
            term_id = 'MI:0686'  # unspecified method
        self.detection_method.term_id = term_id
        self.detection_method.normalize_name(ontology)

        self._get_dm_consistency(ontology, dterms)

        # Counting resolvable conflicts
        if len(set(dterms)) > 1:
            items = frozenset((intr.detection_method.term_id,
                               intr.source_db[0].name)
                              for intr in group
                              if intr.detection_method.term_id is not None)
            self.resolvable_conflicts = items
        else:
            self.resolvable_conflicts = None

    def _collect_author(self, group):
        authors = set()
        for intr in group:
            if intr.author.txt is not None:
                authors.add(intr.author.txt)
        new_author = PlainText()
        new_author.txt = '|'.join(sorted(authors))
        self.author = new_author

    def _consolidate_publications(self, group):
        pubs = set(pub.get()
                   for intr in group
                   for pub in intr.publications.ids)
        for db, acc in pubs:
            pub = SingleID()
            pub.set(db, acc)
            self.publications.ids.append(pub)
        for db, acc in pubs:
            if db == 'pubmed':
                self.pmid = int(acc)
                break

    def _consolidate_interaction_type(self, group, ontology):
        iterms = [intr.interaction_type.term_id for intr in group
                  if intr.interaction_type.term_id is not None]
        if len(iterms):
            term_id = ontology.finest_consistent_term(iterms, self._term_cache)
            self.interaction_type.term_id = term_id
            self.interaction_type.normalize_name(ontology)
        else:
            self.interaction_type.term_id = None

    def _consolidate_source_db(self, group):
        sourcedbs = [item for intr in group for item in intr.source_db]
        self.source_db.from_items(sourcedbs)
        self._sourcedbs = [item.name for item in self.source_db]

    def _consolidate_ids(self, group):
        # We keep original IDs plus irigids plus our own ppiTrimID. The latter
        # will be repeated under interaction checksum.
        srcids = set()
        irigids = set()
        for intr in group:
            srcids.update(intr._get_source_ids())
            irigids.update(intr._get_irigids())

        all_ids = sorted(srcids) + sorted(irigids)
        for db, acc in all_ids:
            item = SingleID()
            item.set(db, acc)
            self.ids.ids.append(item)

    def _consolidate_host_taxid(self, group):
        # Here the result may be inconsistent, that is, grouped interactions
        # may carry different tags. For now we will chose the smallest (as
        # string)
        taxids = sorted(set(intr.host_taxid.get() for intr in group
                            if intr.host_taxid.taxid is not None))
        if len(taxids):
            self.host_taxid.set(*taxids[0])

    def _get_maxsources(self):
        # Count the maximum number of independent sources.

        if self.is_complex():
            # For complexes, we must assume only one experiment per
            # publication, although it is possible that original papers have a
            # much better measure of confidence
            self._maxsources = 1
        else:
            # For binary interactions, this is the maximum number of
            # independent IDs over all source databases.
            src_count = {}
            for db, acc in self._get_source_ids():
                if db not in src_count:
                    src_count[db] = 0
                src_count[db] += 1
            self._maxsources = max(src_count.itervalues())

    def _get_dm_consistency(self, ontology, term_ids):
        if len(term_ids) == 0:
            self._dmconsistency = 'invalid'
            return
        elif len(term_ids) == 1:
            self._dmconsistency = 'single'
            return
        has_max = ontology.maximum(term_ids, self._term_cache) is not None
        has_min = ontology.minimum(term_ids, self._term_cache) is not None
        if has_max and has_min:
            self._dmconsistency = 'full'
            return
        elif has_max and not has_min:
            self._dmconsistency = 'min'
            return
        elif has_min and not has_max:
            self._dmconsistency = 'max'
            return
        else:
            self._dmconsistency = 'invalid'

    def _set_confidence(self):

        maxsources = SingleID()
        maxsources.set('maxsources', self._maxsources)
        dmconsistency = SingleID()
        dmconsistency.set('dmconsistency', self._dmconsistency)
        self.confidence.ids += [maxsources, dmconsistency]

    def _consolidate_all(self, group, ontology, edgetype):
        # Note that the order of consolidation here is important in some cases

        # The fields for: xrefs, annotations, and parameters are left null.
        self._consolidate_interactors(group)
        self._consolidate_detection_method(group, ontology)
        self._collect_author(group)
        self._consolidate_publications(group)
        self._consolidate_interaction_type(group, ontology)
        self._consolidate_ids(group)
        self._consolidate_source_db(group)
        self._consolidate_host_taxid(group)
        self.creation_date.set_today()
        self.update_date.set_today()
        self.negative.txt = 'false'

        self.edgetype = edgetype

        self.ppiTrim_id = self.get_ppiTrim_id()
        self.checksum.set('ppiTrim', self.ppiTrim_id)
        self.ids.ids.insert(0, self.checksum)

        self._consolidate_complex_node(group)

        if self.is_complex():
            self.expansion.txt = 'bipartite'
        else:
            self.expansion.txt = 'none'

        self._get_maxsources()
        self._set_confidence()

    def insert_conflict(self, ppiTrim_id):
        """
        Insert ppiTrim ID of the conflicting interaction
        """

        if len(self._conflicts) == 0:
            cnflct = SingleID()
            self.confidence.ids.append(cnflct)
        else:
            cnflct = self.confidence.ids[-1]

        self._conflicts.append(ppiTrim_id)
        cnflct.set('conflicts', ','.join(self._conflicts))


class DeflatedComplex(ConsolidatedInteraction):
    """
    A complex deflated from a set of binary interactions.
    """

    def __init__(self, group, edgetype, bait_ids=None):

        Interaction.__init__(self)
        # We first set the putative bait so it can be used later for the
        # ppiTrim_id
        if bait_ids is not None:
            self.set_template_baits(bait_ids)
        self._create_complex(group, edgetype)

    def _set_joint_property(self, group, field_name):
        prop_strs = set()
        new_prop = None
        for intr in group:
            new_prop = getattr(intr, field_name)
            prop_strs.add(new_prop.__str__())
        assert len(prop_strs) <= 1
        setattr(self, field_name, new_prop)

    def _create_complex(self, group, edgetype):
        # The fields for: confidence, xrefs, annotations, and parameters are
        # left null.

        # * Set interactors
        self._consolidate_interactors(group)

        # * Set interaction properties that should be unique
        for field_name in ['detection_method',
                           'interaction_type',
                           'host_taxid']:
            self._set_joint_property(group, field_name)

        # Set all others
        self._collect_author(group)
        self._consolidate_publications(group)
        self._consolidate_ids(group)
        self._consolidate_source_db(group)

        self.expansion.txt = 'bipartite'
        self.creation_date.set_today()
        self.update_date.set_today()
        self.negative.txt = 'false'

        self.edgetype = edgetype

        self.ppiTrim_id = self.get_ppiTrim_id()
        self.complex = Interactor()
        self.complex.uid.set('complex', self.ppiTrim_id)
        self.complex.organism = self.interactors[0].organism
        self.complex.biological_role.add('MI:0499', 'unspecified role')
        self.complex.experimental_role.add('MI:0499', 'unspecified role')
        self.complex.interactor_type.add('MI:0315', 'protein complex')
