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
Process candidates for complexes.
"""

import sys
from operator import itemgetter
from idmapping import read_filtered_pmids
from parser import parse_mitab_file
from parser import partial_mitab_iterator
from parser import full_mitab_iterator
from parser import full_interaction_consumer
from parser import offsets_by_pmid_consumer
from counts import Phase2Counts
from interaction import Interaction
from consolidated import DeflatedComplex


class ComplexDeflator(object):

    def __init__(self, logfile_fp, max_complex_size, min_complex_size):

        self.logfile_fp = logfile_fp
        self.max_complex_size = max_complex_size
        self.min_complex_size = min_complex_size

        self.new_complexes = None
        self.used_pairs = None
        self.templates = None
        self.pmid = None

    def __call__(self, pmid, pairs, complexes):

        self.new_complexes = []
        self.used_pairs = set()
        self.templates = {}
        self.pmid = pmid

        grppairs = self._group_binary_interactions(pairs)
        # *** 1. Collapse spoke-expanded associations (must have a bait)
        self.collapse_spokes(grppairs)

        # *** 2. Match already annotated complexes (templates) - only sets
        #        without a bait
        self.match_templates(grppairs, complexes)

        unused_pairs = [intr for intr in pairs if intr not in self.used_pairs]
        retval = (self.new_complexes, unused_pairs)
        self._log_summary(pairs, complexes, unused_pairs)

        self.new_complexes = None
        self.used_pairs = None
        self.templates = None
        self.pmid = None

        return retval

    def collapse_spokes(self, grppairs):
        """
        Collapse spoke-expanded complexes based on bait-prey relationships.
        """
        for key, pair_set in grppairs.iteritems():
            author, db_name, bait, dterm, iterm, host_taxid = key
            if bait is not None:
                if db_name == 'biogrid' and \
                   (dterm, iterm) in [('MI:0401', 'MI:0914'),
                                      ('MI:0401', 'MI:0403')]:
                    # (a) 'Co-purification' and 'Co-fractionation'
                    # interactions. Here they picked a center at random and
                    # then expanded. Hence it is safe to simply go through the
                    # adjacency list and collapse.
                    # Code: G
                    self._add_complex(pair_set, 'G')
                else:
                    # (b) Everything else (Affinity chromatography and
                    # similar).
                    # Code: A
                    self._add_complex(pair_set, 'A')

    def match_templates(self, grppairs, complexes):
        """
        Match sets of interactions without bait/prey labels to existing
        complexes.
        """
        self._insert_templates(complexes, 'R')  # (original complexes)
        self._log_templates()
        self._insert_templates(self.new_complexes, 'N')  # (deflated spokes)
        self._clean_templates()
        for key, pair_set in grppairs.iteritems():
            author, db_name, bait, dterm, iterm, host_taxid = key
            if bait is None:
                self._match_templates(pair_set)

    def _add_complex(self, pairs, code, bait_ids=None):
        if bait_ids is not None:
            # Ensure that bait ids are strings
            bait_ids = map(str, bait_ids)
        cmplx = DeflatedComplex(pairs, code, bait_ids)
        self.new_complexes.append(cmplx)
        self.used_pairs.update(pairs)

    def _group_binary_interactions(self, pairs):
        # ** Group pairs by experiment (for each database)
        # All interactions from an experiment must have the same:
        #
        # - author tag (IntAct and MINT use them to distinguish experiments)
        #
        # - source database (We use source_db_name because several databases
        #                    are assigned MI:0000 as the PSI-MI term. We assume
        #                    exactly one term for source_db)
        #
        # - bait protein (if annotated, otherwise this field is set to None)
        #
        # - detection method
        #
        # - interaction type (if present)
        #
        # - host taxonomy id (if present)
        grppairs = {}
        for intr in pairs:

            baits = intr.find_baits()
            if len(baits) > 0:
                bait_acc = baits[0].uid.acc
            else:
                bait_acc = None

            key = (intr.author.txt,
                   intr.source_db[0].name,
                   bait_acc,
                   intr.detection_method.term_id,
                   intr.interaction_type.term_id,
                   intr.host_taxid.taxid)

            if key not in grppairs:
                grppairs[key] = []
            grppairs[key].append(intr)

        protein_sets = {}
        for key, pair_set in grppairs.items():
            geneids = frozenset(int(p.uid.acc)
                                for intr in pair_set
                                for p in intr.interactors)
            if self.min_complex_size <= len(geneids) <= self.max_complex_size:
                protein_sets[key] = geneids
            else:
                grppairs.pop(key)

        self._log_grouped_pairs(grppairs, protein_sets)
        return grppairs

    def _insert_templates(self, complexes, template_code='R'):
        # ** Reduce existing complexes to unique representatives to be used as
        #    templates (code R). We only require a list of proteins (geneids)
        #    but we also keep track of original interactions for debugging.
        for intr in complexes:
            geneids = frozenset(int(p.uid.acc) for p in intr.interactors)
            if len(geneids) >= self.min_complex_size:
                if geneids not in self.templates:
                    self.templates[geneids] = []
                baits = intr.find_baits()
                if baits:
                    bait_ids = tuple(sorted(p.uid.acc for p in baits))
                else:
                    bait_ids = None
                self.templates[geneids].append((template_code, intr, bait_ids))

    def _clean_templates(self):
        # It is only necessary to have the templates with baits. If none could
        # be found, we will only need a single template

        # Also, if there are complexes with baits in the group, insert putative
        # baits in the original complexes without baits so they could be
        # brought together at Phase 3
        for geneids in self.templates:
            with_baits = [item for item in self.templates[geneids]
                          if item[2] is not None]
            without_baits = [item for item in self.templates[geneids]
                             if item[2] is None]
            if with_baits:
                self.templates[geneids] = with_baits
                bait_ids = with_baits[0][2]
                for _, cmplx, _ in without_baits:
                    cmplx.set_template_baits(bait_ids)
            else:
                self.templates[geneids] = [self.templates[geneids][0]]

    def _get_nbhd_sets(self, pairs):
        # Assume that all interactions in pairs are undirected (i.e. no
        # bait/prey distinction)
        nbhd_sets = {}
        for intr in pairs:
            geneids = [int(p.uid.acc) for p in intr.interactors]
            for gene_id in geneids:
                if gene_id not in nbhd_sets:
                    nbhd_sets[gene_id] = [set(), []]
                nbhd_sets[gene_id][0].update(geneids)
                nbhd_sets[gene_id][1].append(intr)
        return nbhd_sets

    def _get_template_items(self):
        for template, items in self.templates.iteritems():
            for code, cmplx, bait_ids in items:
                baitids = None
                if bait_ids is not None:
                    baitids = map(int, bait_ids)  # must convert to ints
                yield template, code, cmplx, baitids

    def _match_templates(self, pairs):

        for center, val in self._get_nbhd_sets(pairs).iteritems():
            geneids, nbhd_pairs = val
            assert center in geneids

            deflations = []
            for template, code, cmplx, baitids in self._get_template_items():
                if center in template and \
                   (baitids is None or center in baitids) and \
                   geneids.issuperset(template):

                    matched_pairs = []
                    for intr in nbhd_pairs:
                        match = [int(p.uid.acc) in template
                                 for p in intr.interactors]
                        if all(match):
                            matched_pairs.append(intr)

                    jaccard = 1.0 * len(template) / len(geneids)
                    mdata = (jaccard, matched_pairs, template, code, cmplx,
                             baitids)
                    deflations.append(mdata)

            # Insert all matched templates, removing any that are redundant
            coverage = set()
            deflations.sort(reverse=True, key=itemgetter(0))
            for mdata in deflations:
                jaccard, matched_pairs, template, code, cmplx, baitids = mdata
                if geneids - coverage:
                    coverage |= geneids
                    self._add_complex(matched_pairs, code, baitids)
                    self._log_template_match(code, cmplx, matched_pairs,
                                             jaccard, template, geneids,
                                             center, baitids)

    def _log_grouped_pairs(self, grppairs, protein_sets):
        self.logfile_fp.write('------------- PMID: %d' % self.pmid)
        self.logfile_fp.write(' -------------\n')
        sorted_keys = sorted(grppairs.keys(), key=itemgetter(2))
        for key in sorted_keys:
            bait = key[2]
            if bait is not None:
                code = 'A'
            else:
                code = 'T'
            pair_set = grppairs[key]
            keystr = '\t'.join(map(str, key))
            self.logfile_fp.write('@ %s  %s\t%d\t%d\n' %
              (code, keystr, len(pair_set), len(protein_sets[key])))

    def _log_templates(self):
        for template, code, cmplx, baitids in self._get_template_items():
            data = [str(cmplx.author.txt),
                    cmplx.source_db[0].name,
                    str(baitids),
                    cmplx.detection_method.term_id,
                    cmplx.interaction_type.term_id,
                    str(len(template))]
            self.logfile_fp.write('& C  %s\n' % '\t'.join(data))

    def _log_template_match(self, code, cmplx, matched_pairs, jaccard,
                            template, geneids, center, baitids):
        intr = matched_pairs[0]
        merge_size = '%.2f (%d/%d)' % (jaccard, len(template), len(geneids))
        data = [code,
                merge_size,
                cmplx.source_db[0].name,
                str(baitids),
                cmplx.detection_method.term_id,
                intr.source_db[0].name,
                str(center),
                intr.detection_method.term_id]
        self.logfile_fp.write('# %s\n' % '  '.join(data))

    def _log_summary(self, pairs, complexes, unused_pairs):
        reduction = 100.0 - 100.0 * len(unused_pairs) / len(pairs)
        mncs = 0
        if len(self.new_complexes):
            mncs = max(len(intr.interactors) for intr in self.new_complexes)
        self.logfile_fp.write('pairs=%d ' % len(pairs))
        self.logfile_fp.write('complexes=%d ' % len(complexes))
        self.logfile_fp.write('unused_pairs=%d ' % len(unused_pairs))
        self.logfile_fp.write('new_complexes=%d ' % len(self.new_complexes))
        self.logfile_fp.write('reduction=%.1f%% ' % reduction)
        self.logfile_fp.write('max_new_complex_size=%d\n\n' % mncs)


def _parse_by_pmid(input_fp):

    # Will go smallest to largest pubmed id, in two passes
    # First pass: collect pubmed ids with pointers to their interactions
    scanner = parse_mitab_file(input_fp, full_mitab_iterator)
    pmids_map = offsets_by_pmid_consumer(scanner)
    sorted_pmids = sorted((len(v), k) for k, v in pmids_map.iteritems())

    # Second pass
    for _, pmid in sorted_pmids:
        scanner = parse_mitab_file(input_fp, partial_mitab_iterator,
                                   (pmids_map[pmid],))
        pairs, complexes = full_interaction_consumer(scanner)
        yield pmid, pairs, complexes


def _write_unused_pairs(output_fp, unused_pairs):
    for interaction in unused_pairs:
        interaction.to_file(output_fp)


def _write_existing_complexes(output_fp, complexes):
    for interaction in complexes:
        ppiTrim_id = interaction.get_ppiTrim_id()
        interaction.complex.uid.set('complex', ppiTrim_id)
        interaction.to_file(output_fp)


def _write_new_complexes(output_fp, new_complexes):

    for cmplx in new_complexes:
        cmplx.to_file(output_fp)


def process_ppi2(input_file, output_file, output_logfile,
                 skipped_pmids_file=None, max_complex_size=120,
                 min_complex_size=3):

    counts = Phase2Counts()
    logfile_fp = open(output_logfile, 'w')
    skipped_pmids = read_filtered_pmids(skipped_pmids_file)
    skipped_pmids.add(0)  # this is invalid pmid - not really a paper
    input_fp = open(input_file, 'rU')

    output_fp = open(output_file, 'w')
    Interaction.write_header(output_fp)

    deflator = ComplexDeflator(logfile_fp, max_complex_size, min_complex_size)

    for pmid, pairs, complexes in _parse_by_pmid(input_fp):
        counts.initial_pairs += len(pairs)
        counts.C  += len(complexes)
        counts.pmids += 1

        if pmid in skipped_pmids or len(pairs) < (min_complex_size - 1):
            new_complexes = []
            unused_pairs = pairs
        else:
            new_complexes, unused_pairs = deflator(pmid, pairs, complexes)

        _write_unused_pairs(output_fp, unused_pairs)
        _write_existing_complexes(output_fp, complexes)
        _write_new_complexes(output_fp, new_complexes)

        counts.unused_pairs += len(unused_pairs)
        for intr in new_complexes:
            code = intr.edgetype
            counter = getattr(counts, code)
            counter += 1
            setattr(counts, code, counter)

    counts.to_file(logfile_fp)
    input_fp.close()
    output_fp.close()
    logfile_fp.close()


if __name__ == "__main__":

    args = sys.argv[1:]
    process_ppi2(*args)
