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
Produce final consolidated interactions.
"""
import sys
from parser import parse_mitab_file
from parser import full_mitab_iterator
from parser import partial_mitab_iterator
from parser import full_interaction_consumer
from counts import Phase3Counts
from consolidated import ConsolidatedInteraction
import obo


def _collect_by_pmid(pairs, key_func, pmid_links=None, unassigned=None):
    if pmid_links is None:
        pmid_links = {}
    if unassigned is None:
        unassigned = []
    for intr in pairs:
        key = key_func(intr)
        if key is None:
            unassigned.append(intr)
        else:
            if key not in pmid_links:
                pmid_links[key] = []
            pmid_links[key].append(intr)
    return pmid_links, unassigned


def _get_biochem_key(intr):
    # Assumes intr is truly directed
    p0 = intr.interactors[0]
    p1 = intr.interactors[1]

    if intr.source_db[0].term_id == 'MI:0463':
        # All BioGRID's biochemical interactions are in the form bait -> prey
        # and hence they are confirmed. We ignore 'biological role' since it is
        # likely in error.
        key = (p0.uid.acc, p1.uid.acc, intr.pmid)
    elif p0.uid.acc == p1.uid.acc:
        # Self-targetted actions are also OK
        key = (p0.uid.acc, p1.uid.acc, intr.pmid)
    # For the rest, check biological roles
    elif not len(p0.biological_role) or not len(p1.biological_role):
        key = None
    elif p0.biological_role[0].term_id == 'MI:0501' and \
         p1.biological_role[0].term_id == 'MI:0502':
        key = (p0.uid.acc, p1.uid.acc, intr.pmid)
    elif p1.biological_role[0].term_id == 'MI:0501' and \
         p0.biological_role[0].term_id == 'MI:0502':
        key = (p1.uid.acc, p0.uid.acc, intr.pmid)
    else:
        # Cannot get the direction
        key = None
    return key


def _get_undirected_key(intr):
    ids = sorted(p.uid.acc for p in intr.interactors)
    key = tuple(ids) + (intr.pmid,)
    return key


def _get_complex_key(intr):

    baits = intr.find_baits()
    bait_ids = sorted(p.uid.acc for p in baits)
    has_baits = len(bait_ids) > 0
    if has_baits:
        key = (None, ) + (has_baits,) + tuple(bait_ids) + (intr.pmid,)
    else:
        ids = sorted(p.uid.acc for p in intr.interactors)
        key = tuple(ids) + (has_baits,) + (None, ) + (intr.pmid,)
    return key


def _consume_complex_offsets_only(scanner):
    offsets_map = {}
    for interaction, lines in scanner:
        line_offsets = lines[0]
        if interaction.is_complex():
            key = _get_complex_key(interaction)
            if key not in offsets_map:
                offsets_map[key] = []
            offsets_map[key].extend(line_offsets)
    return offsets_map


def _consume_undirected_offsets_only(scanner):
    offsets_map = {}
    for interaction, lines in scanner:
        line_offsets = lines[0]
        if not interaction.is_complex():
            key = _get_undirected_key(interaction)
            if key not in offsets_map:
                offsets_map[key] = []
            offsets_map[key].extend(line_offsets)
    return offsets_map


def _write_directional_conflict(key, src1, src2, logfile_fp):
    logfile_fp.write('* DIRECTED INTERACTION METHOD CONFLICT *\n')
    logfile_fp.write('  pmid=%d, interactors=[%s, %s]\n' %
                     (key[2], key[0], key[1]))
    logfile_fp.write('  1. %s\n' % ', '.join(sorted(src1)))
    logfile_fp.write('  2. %s\n\n' % ', '.join(sorted(src2)))


def _get_biochem_key_conflicts(pmid_links, logfile_fp):
    # Find the cases (keys) where there is a directed interaction A -> B
    # as well as B -> A reported by the same pmid. The function returns the
    # list of all such conflicts (key_conflicts) as pairs of keys, plus the
    # list of conflicts between different source databases  (e.g. IntAct and
    # BioGRID  directly contradict each other about which protein is an enzyme
    # and which is its target). The second list (called real_conflicts) is most
    # likely to point to curation errors by one or the other database.
    _key_conflicts = set()
    for accA_0, accA_1, pmidA in pmid_links:
        if accA_0 == accA_1:
            continue
        for accB_0, accB_1, pmidB in pmid_links:
            if pmidA == pmidB and (accA_0, accA_1) == (accB_1, accB_0):
                _key_conflicts.add((frozenset((accA_0, accA_1)), pmidA))
    key_conflicts = []
    for proteins, pmid in _key_conflicts:
        acc0, acc1 = tuple(proteins)
        key_conflicts.append(((acc0, acc1, pmid), (acc1, acc0, pmid)))
    real_conflicts = []
    for key, key2 in key_conflicts:
        src1 = set('%s(%s)' % (intr.source_db[0].term_id,
                               intr.source_db[0].name)
                   for intr in pmid_links[key])
        src2 = set('%s(%s)' % (intr.source_db[0].term_id,
                               intr.source_db[0].name)
                   for intr in pmid_links[key2])
        if src1 != src2:
            real_conflicts.append((key, key2))
            _write_directional_conflict(key, src1, src2, logfile_fp)
    return key_conflicts, real_conflicts


def _collect_biochem_unassigned(unassigned, pmid_links):
    # Attempt to match the biochemical interactions without specified
    # biological role (so far only from DIP) to the directed interactions
    # already collected from other databases.
    unresolved = []
    for intr in unassigned:
        p0 = intr.interactors[0]
        p1 = intr.interactors[1]
        if not len(p0.biological_role) or not len(p1.biological_role) or \
           (p0.biological_role[0].term_id in (None, 'MI:0499') and
            p1.biological_role[0].term_id in (None, 'MI:0499')):
            key = (p0.uid.acc, p1.uid.acc, intr.pmid)
            key2 = (p1.uid.acc, p0.uid.acc, intr.pmid)
            # Here we allow for the possibility that this 'undirected'
            # interaction matches interactions directed both ways
            if key in pmid_links:
                pmid_links[key].append(intr)
            if key2 in pmid_links:
                pmid_links[key2].append(intr)
            if key not in pmid_links and key2 not in pmid_links:
                unresolved.append(intr)
        else:
            unresolved.append(intr)
    return unresolved


def _extract_biochem_from_complexes(complexes):
    # Attempt to extract directed biochemical reactions from 'complexes'
    # This will lose additional information since only enzymes and their
    # targets are kept

    extracted_pairs = []
    remaining_complexes = []
    for cmplx in complexes:
        enzymes = [p for p in cmplx.interactors
                   if len(p.biological_role) and
                      p.biological_role[0].term_id == 'MI:0501']
        targets = [p for p in cmplx.interactors
                   if len(p.biological_role) and
                      p.biological_role[0].term_id == 'MI:0502']

        if len(enzymes) == 1 and len(targets) == 1:
            intr = cmplx.binary_from_complex(enzymes[0], targets[0], 'H')
            extracted_pairs.append(intr)
        else:
            remaining_complexes.append(cmplx)

    return extracted_pairs, remaining_complexes


def _correct_biogrid_biochem(interactions):
    # Must correct biological roles because BioGRID's assignment is
    # wrong. A bait is always enzyme, prey is the target.
    for intr in interactions:
        if intr.source_db[0].name == 'biogrid':
            baits = intr.find_baits()
            assert len(baits) == 1
            bait = baits[0]
            prey = next(p for p in intr.interactors if p is not bait)

            # BioGRID entries are assumed to always have experimental role set
            # as bait/prey
            assert bait.experimental_role[0].term_id == 'MI:0496'
            assert prey.experimental_role[0].term_id == 'MI:0498'

            bait.biological_role[0].set('MI:0501', 'enzyme')
            prey.biological_role[0].set('MI:0502', 'enzyme target')


def _group_sources(sources, term_func, ontology, cache):

    # This ought to speed up the most common case
    if len(sources) == 1:
        return [sources]

    # For each element, get all elements 'smaller' than it ('None' is smaller
    # than any true term)
    groups_ = set()
    obj_dict = {}
    for a in sources:
        grp = set()
        term_id1 = term_func(a)
        if term_id1 is None:
            for b in sources:
                term_id2 = term_func(b)
                if term_id2 is None:
                    obj_dict[id(b)] = b
                    grp.add((term_id2, id(b)))
        else:
            term1 = ontology.get_term(term_id1)
            for b in sources:
                term_id2 = term_func(b)
                if term_id2 is None:
                    is_smaller = True
                else:
                    term2 = ontology.get_term(term_id2)
                    is_smaller = term2.compare_to(term1, cache) >= 0
                if is_smaller:
                    obj_dict[id(b)] = b
                    grp.add((term_id2, id(b)))
        groups_.add(frozenset(grp))
    # Remove redundant covers
    for grp in groups_.copy():
        for grp2 in groups_.copy():
            if grp < grp2:
                groups_.remove(grp)
                break
    groups = []
    for grp in groups_:
        groups.append([obj_dict[oid] for _, oid in grp])
    return groups


def _log_method_conflict(intr1, intr2, logfile_fp, counts):

    logfile_fp.write('* DETECTION METHOD CONFLICT *\n')
    logfile_fp.write('  pmid=%s, interactors=%s\n' % (str(intr1.pmid),
                     ', '.join('%s' % p.uid.acc for p in intr1.interactors)))
    term1 = '%s (%s)' % intr1.detection_method.get()
    term2 = '%s (%s)' % intr2.detection_method.get()
    int_ids1 = ['%s:%s' % (db, ac) for db, ac in intr1._get_source_ids()]
    int_ids2 = ['%s:%s' % (db, ac) for db, ac in intr2._get_source_ids()]
    logfile_fp.write('  1. %s -- %s [%s]\n' %
                     (intr1.ppiTrim_id, term1, ', '.join(int_ids1)))
    logfile_fp.write('  2. %s -- %s [%s]\n\n' %
                     (intr2.ppiTrim_id, term2, ', '.join(int_ids2)))
    dbs1 = set(intr1._sourcedbs)
    dbs2 = set(intr2._sourcedbs)
    cnflct = sorted([(term1, dbs1), (term2, dbs2)])
    key = (cnflct[0][0], cnflct[1][0])
    counts.add_unresolvable(key, cnflct[0][1], cnflct[1][1])


def _check_method_conflicts(consolidated, logfile_fp, counts):

    for i, intr1 in enumerate(consolidated):
        # Resolvable conflicts are recorded for each consolidated interaction
        # at its creation time
        if intr1.resolvable_conflicts is not None:
            counts.add_resolvable(intr1.resolvable_conflicts)

        # Check for unresolvable conflicts
        dbs1 = set(intr1._sourcedbs)
        for intr2 in consolidated[i + 1:]:
            dbs2 = set(intr2._sourcedbs)
            if len(dbs1 & dbs2) == 0:
                # Annotate interactions with their conflicts
                intr1.insert_conflict(intr2.ppiTrim_id)
                intr2.insert_conflict(intr1.ppiTrim_id)
                # Write conflict to logfile
                _log_method_conflict(intr1, intr2, logfile_fp, counts)


def _get_detection_method_term(intr):
    term_id = intr.detection_method.term_id
    # Labels for 'unspecified method', 'in-vitro' and 'in-vivo' are turned
    # into 'None' to be compatible with all other terms. The latter two
    # come from HPRD only
    if term_id == 'MI:0686' or \
       term_id == 'MI:0492' or \
       term_id == 'MI:0493':
        term_id = None
    return term_id


def _log_missing_itype(intr1, logfile_fp):

    logfile_fp.write('* MISSING INTERACTION TYPE LABEL *\n')
    logfile_fp.write('  pmid=%d ' % intr1.pmid)
    logfile_fp.write('interactors=%s\n' %
                      ', '.join('%s' % p.uid.acc for p in intr1.interactors))
    term1 = '%s (%s)' % intr1.detection_method.get()
    int_ids1 = ['%s:%s' % (db, ac) for db, ac in intr1._get_source_ids()]
    logfile_fp.write('  %s [%s]\n' % (term1, ', '.join(int_ids1)))


def _consolidate_links(sources, ontology, code, logfile_fp, counts):

    consolidated = []
    groups = _group_sources(sources, _get_detection_method_term, ontology,
                            ConsolidatedInteraction._term_cache)
    for grp in groups:
        new_intr = ConsolidatedInteraction(grp, ontology, code)
        consolidated.append(new_intr)

        # If interaction type is unknown, set it to physical association
        iterm_id = new_intr.interaction_type.term_id
        if iterm_id is None:
            new_intr.interaction_type.set('MI:0915', 'physical association')
            _log_missing_itype(new_intr, logfile_fp)

    _check_method_conflicts(consolidated, logfile_fp, counts)
    return consolidated


def _log_complex_conflict(conflicting, logfile_fp):

    logfile_fp.write('* COMPLEX CONFLICT (pmid=%d) (%d complexes)*\n' %
                     (conflicting[0].pmid, len(conflicting)))

    for i, intr in enumerate(conflicting):
        prots = [p.uid.acc for p in intr.interactors]
        baits = '[%s]' % ', '.join(p.uid.acc for p in intr.baits)
        dbs = '[%s]' % ', '.join(term.name for term in intr.source_db)
        dmterm = '%s (%s)' % intr.detection_method.get()
        # int_ids = ['%s:%s' % (db, ac) for db, ac in intr._get_source_ids()]

        logfile_fp.write(' %d. bait=%s (%d members) %s\n' %
                         (i+1, baits, len(prots), intr.ppiTrim_id))
        logfile_fp.write('     %s  %s\n' % (dbs, dmterm))
        logfile_fp.write('     interactors=%s\n\n' %
                         ', '.join('%s' % p.uid.acc for p in intr.interactors))


def _consolidate_complexes(sources, ontology, logfile_fp):
    # This is used to consolidate together protein complexes.
    # We assume that interactions (complexes) in sources list either all have
    # the same bait(s) or they don't have any bait but have the same sets of
    # interactors.

    # Groups based on compatibility of detection methods
    dm_groups = _group_sources(sources, _get_detection_method_term,
                               ontology, ConsolidatedInteraction._term_cache)

    # If the sources are grouped by bait, we need to split each group further
    # based on members. If they are grouped by members (no bait), they will all
    # form one cluster anyway.
    dm_consolidated = []
    for grp in dm_groups:
        # Here, we split based on members
        members_map = {}
        for intr in grp:
            # This is relaxed - it allows that a protein can be
            # bait and prey in one case but only bait in another
            key = frozenset(p.uid.acc for p in intr.interactors)
            if key not in members_map:
                members_map[key] = []
            members_map[key].append(intr)

        # NOTE: we might be computing baits too many times:
        # - once to get the key for grouping sources
        # - once here, and
        # - once when we construct a ConsolidatedInteraction
        # This is wasteful but looks unavoidable because we don't know when
        # a complex is fully finished after parsing. At least this way, this
        # computation is not hidden away somewhere.
        baits = sources[0].find_baits()
        have_baits = len(baits) > 0
        assert have_baits or len(members_map) == 1

        dm_subgroup = []
        for subgrp in members_map.values():
            # Assign a concatenated edgetype to consolidated complexes
            all_codes = set()
            for intr in subgrp:
                all_codes.update(list(intr.edgetype))
            code = ''.join(sorted(all_codes))

            # Construct a consolidated interaction
            new_intr = ConsolidatedInteraction(subgrp, ontology, code)

            # Interaction type of complexes cannot be finer than 'association'
            iterm_id = new_intr.interaction_type.term_id
            if iterm_id is None:
                new_intr.interaction_type.set('MI:0914', 'association')

            # Record all subgroups - if there is more than one, they are in
            # conflict (same pubmed, same bait, same detection method,
            # different sets of preys)
            dm_subgroup.append(new_intr)
        dm_consolidated.append(dm_subgroup)

    # There are two types of conflicts:
    # - within each dm_group, if complex members are different
    # - between dm groups, based on detection method (as for binary
    # interactions)
    # Both types are merged (treated as 'complex conflict')

    # First label within dm_subgroup conflicts and collect source databases
    source_dbs = []
    for dm_subgroup in dm_consolidated:
        dbs = set()
        members_conflict = False
        for i, intr1 in enumerate(dm_subgroup):
            dbs |= set(intr1._sourcedbs)

            # Members conflict: this loop will run only if dm_subgroup has
            # length more than 1
            for intr2 in dm_subgroup[i+1:]:
                members_conflict = True
                intr1.insert_conflict(intr2.ppiTrim_id)
                intr2.insert_conflict(intr1.ppiTrim_id)
        source_dbs.append((members_conflict, dbs))

    # Now compare dm_subgroups and construct conflict groups. Each conflict
    # group should represent what we beleive ought to be a single evidence. So,
    # conflict here means that something is stopping these records to be merged
    # into a single one. Note that dm_conflicts are not transitive (it is
    # theoretically possible to have A -> B and B -> C but not A -> C) so this
    # code could possibly not detect all the conflicts. This issue is too
    # complicated to be properly solved right now and this case is very
    # unlikely anyway.
    conflict_groups = []
    unvisited = set(range(len(source_dbs)))

    while unvisited:
        i = unvisited.pop()
        dm_subgroup1 = dm_consolidated[i]
        srcdbs1 = source_dbs[i][1]
        sim = set()
        for j in unvisited:
            dm_subgroup2 = dm_consolidated[j]
            srcdbs2 = source_dbs[j][1]
            if len(srcdbs1 & srcdbs2) == 0:
                # dm_conflict
                for intr1 in dm_subgroup1:
                    for intr2 in dm_subgroup2:
                        intr1.insert_conflict(intr2.ppiTrim_id)
                        intr2.insert_conflict(intr1.ppiTrim_id)
                sim.add(j)
        for j in sim:
            unvisited.remove(j)
        if len(sim) or source_dbs[i][0]:
            sim.add(i)
            conflict_groups.append(sim)

    # Now write to log each conflict group
    for cgrp_ix in conflict_groups:
        cgrp = [intr for i in cgrp_ix for intr in dm_consolidated[i]]
        _log_complex_conflict(cgrp, logfile_fp)

    # Finally, collect all consolidated interactions
    consolidated = [intr for dm_sbgrp in dm_consolidated for intr in dm_sbgrp]
    return consolidated


def process_biochem(biochem_file, output_fp, ontology, logfile_fp, counts):
    """
    Process all candidate 'biochemical' interactions in biochem_file and write
    consolidated interactions to output_fp.
    """
    input_fp = open(biochem_file, 'rU')
    scanner = parse_mitab_file(input_fp, full_mitab_iterator)
    pairs, complexes = full_interaction_consumer(scanner)

    counts.pmids.update(intr.pmid for intr in pairs)
    counts.pmids.update(intr.pmid for intr in complexes)

    # Extract additional interactions from 'complexes', which are in fact
    # biochemical reactions described in more detail
    pairs2, complexes2 = _extract_biochem_from_complexes(complexes)
    # Collect source interactions by (directed) pair of proteins plus pubmed id
    pmid_links, unassigned = _collect_by_pmid(pairs + pairs2, _get_biochem_key)
    # Find conflicts and inconsistencies
    _get_biochem_key_conflicts(pmid_links, logfile_fp)
    # Attempt to assign interactions from DIP, which are not directional
    unresolved = _collect_biochem_unassigned(unassigned, pmid_links)

    # Consolidate each interaction and write to file
    sorted_keys = sorted(pmid_links.keys())
    for key in sorted_keys:
        interactions = pmid_links[key]
        counts.biochem += len(interactions)

        _correct_biogrid_biochem(interactions)

        for new_intr in _consolidate_links(interactions, ontology, 'D',
                                           logfile_fp, counts):
            new_intr.to_file(output_fp)
            counts.directed += 1

    # Write all 'unused' interactions
    for item in unresolved:
        key = _get_undirected_key(item)
        for new_intr in _consolidate_links([item], ontology, 'B',
                                           logfile_fp, counts):
            new_intr.to_file(output_fp)
            counts.undirected += 1

    # Write unused complexes
    for cmplx in complexes2:
        key = _get_undirected_key(cmplx)
        for new_intr in _consolidate_links([cmplx], ontology, 'C',
                                           logfile_fp, counts):
            new_intr.to_file(output_fp)
            counts.complexes += 1

    input_fp.close()


def process_complexes(complexes_file, output_fp, ontology, logfile_fp, counts):
    """
    Process all candidate 'complexes' in complexes_file and write
    consolidated interactions to output_fp.
    """
    input_fp = open(complexes_file, 'rU')
    scanner = parse_mitab_file(input_fp, full_mitab_iterator)
    offsets_map = _consume_complex_offsets_only(scanner)

    # Consolidate each interaction and write to file
    sorted_keys = sorted(offsets_map.iterkeys())
    for key in sorted_keys:
        scanner = parse_mitab_file(input_fp, partial_mitab_iterator,
                                   (offsets_map[key], ))
        _, complexes = full_interaction_consumer(scanner)
        counts.pmids.update(intr.pmid for intr in complexes)

        for new_intr in _consolidate_complexes(complexes, ontology,
                                               logfile_fp):
            new_intr.to_file(output_fp)
            counts.complexes += 1
    input_fp.close()


def process_undirected(binary_file, output_fp, ontology, logfile_fp, counts):
    """
    Process all 'other' interactions - undirected physical interactions. Write
    consolidated interactions to output_fp.
    """
    input_fp = open(binary_file, 'rU')
    scanner = parse_mitab_file(input_fp, full_mitab_iterator)
    offsets_map = _consume_undirected_offsets_only(scanner)

    # Consolidate each interaction and write to file
    sorted_keys = sorted(offsets_map.iterkeys())
    for key in sorted_keys:
        scanner = parse_mitab_file(input_fp, partial_mitab_iterator,
                                   (offsets_map[key], ))
        pairs, _ = full_interaction_consumer(scanner)
        counts.pmids.update(intr.pmid for intr in pairs)
        counts.other += len(pairs)

        for new_intr in _consolidate_links(pairs, ontology, 'X', logfile_fp,
                                           counts):
            new_intr.to_file(output_fp)
            counts.undirected += 1
    input_fp.close()


def process_ppi3(biochem_file, complexes_file, binary_file, obo_file,
                 output_file, output_logfile):

    counts = Phase3Counts()
    output_fp = open(output_file, 'w')
    ConsolidatedInteraction.write_header(output_fp)
    logfile_fp = open(output_logfile, 'w')
    obo_fp = open(obo_file, 'rU')
    ontology = obo.OBOntology(obo_fp)

    process_biochem(biochem_file, output_fp, ontology, logfile_fp, counts)
    process_complexes(complexes_file, output_fp, ontology, logfile_fp, counts)
    process_undirected(complexes_file, output_fp, ontology, logfile_fp, counts)
    process_undirected(binary_file, output_fp, ontology, logfile_fp, counts)

    counts.to_file(logfile_fp, ontology)

    obo_fp.close()
    output_fp.close()
    logfile_fp.close()


if __name__ == "__main__":

    args = sys.argv[1:]
    process_ppi3(*args)
