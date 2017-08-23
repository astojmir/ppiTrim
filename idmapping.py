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
Map ROGIDs to NCBI Gene IDs.
"""

import time
import sys
import re
import os
import hashlib
import base64
import urllib
import cPickle as pickle
from download import download_file
from download import download_and_extract_file
from parser import parse_mitab_file
from parser import full_mitab_iterator
from filter import is_filtered
from filter import NullFile
from interaction import iRefIndexInteraction
from counts import Phase1Counts

uniprot_url = 'http://www.uniprot.org/uniprot/%s.txt'
gn_ptrn = re.compile(r'Name=(.+?);')   # canonical gene name
gn_ptrn2 = re.compile(r'OrderedLocusNames=(.+?);')   # canonical gene name
pe_ptrn = re.compile(r'([1-5]:\s.+);')  # protein existence
geneid_ptrn = re.compile(r'GeneID; ([0-9].+?);')
acc_ptrn = re.compile(r'([A-Z][0-9][A-Z,0-9][A-Z,0-9][A-Z, 0-9][0-9]);')
ox_ptrn = re.compile(r'NCBI_TaxID=([0-9].+?);')   # NCBI Taxonomy ID

gene_history_url = 'ftp://ftp.ncbi.nih.gov/gene/DATA/gene_history.gz'
gene_info_url = 'ftp://ftp.ncbi.nih.gov/gene/DATA/gene_info.gz'


uniprot_ftp_url = 'ftp://ftp.uniprot.org/pub/databases/uniprot/' \
                  'current_release/knowledgebase/complete/'
db_fetch_url = 'http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=uniprotkb' \
               '&id=%s&format=uniprot&style=raw'


def read_filtered_pmids(filtered_pmids_file):
    """
    Read a list of IDs to be filtered from a file
    """
    filtered_pmids = set()
    if filtered_pmids_file is not None:
        fp = open(filtered_pmids_file, 'rU')
        for line in fp:
            cols = line.strip().split('\t')
            filtered_pmids.add(int(cols[0]))
        fp.close()
    return filtered_pmids


def extract_protein_ids(mitab_file, filtered_pmids, accepted_taxids):
    """
    Extract protein IDs from iRefIndex PSI-MI TAB 2.6 file
    """

    protein_ids = set()
    nullfile = NullFile()

    input_fp = open(mitab_file, 'rU')
    scanner = parse_mitab_file(input_fp, full_mitab_iterator, None,
                               iRefIndexInteraction)
    for interaction, lines in scanner:
        line_numbers = lines[1]
        if is_filtered(interaction, filtered_pmids, nullfile, line_numbers,
                       accepted_taxids):
            continue

        for p in interaction.interactors:
            # p stands for 'protein interactor'
            protein_ids.add(p.id)

    input_fp.close()
    return protein_ids


def _add_GR_link(GRgraph, geneid, rogid, code):

    if rogid not in GRgraph:
        GRgraph[rogid] = dict()
    if geneid not in GRgraph:
        GRgraph[geneid] = dict()
    GRgraph[rogid][geneid] = code
    GRgraph[geneid][rogid] = code


def gene_ids_from_mapping_txt(protein_ids, mapping_file):
    """
    Extract mappings of ROGIDs to Gene ID from mapping.txt file provided by
    iRefIndex. Also records mapping of ROGIDS to canonical ROGIDs for
    reference.
    """

    GRgraph = dict()
    resolved = set()
    orphan_ids = dict()
    icrogids = dict()
    fp = open(mapping_file, 'rU')
    for line in fp:
        fields = line.split('\t')
        db, acc, geneid_, _, rogid, icrogid = fields[:6]
        geneid = int(geneid_)
        if rogid in protein_ids:
            icrogids[rogid] = int(icrogid)
            if geneid != -1:
                # There could be multiple Gene IDs associated with a rogid. We
                # want to get all of them.
                _add_GR_link(GRgraph, geneid, rogid, 'M')
                resolved.add(rogid)
                orphan_ids.pop(rogid, None)
            elif rogid not in resolved:
                if rogid not in orphan_ids:
                    orphan_ids[rogid] = set()
                orphan_ids[rogid].add((db, acc))
    fp.close()
    return GRgraph, orphan_ids, icrogids


def _retrieve_uniprot_file(uniprot_acc, download_dir):
    """
    Download a Uniprot file by accession.
    """
    src = uniprot_url % uniprot_acc
    src = db_fetch_url % uniprot_acc
    dst = os.path.join(download_dir, '%s.txt' % uniprot_acc)
    dst = download_file(src, dst)
    return dst


def _get_ids_hash(ids):

    hash_factory = hashlib.sha1()
    hash_factory.update(ids)
    return base64.b32encode(hash_factory.digest())


def _batch_retrieve_uniprot_files(uniprot_accessions, download_dir,
                                  dataset_prefix, initial_batch_size=120,
                                  sleep_interval=15):
    """ Download Uniprot files in batches """

    N = len(uniprot_accessions)
    batch_size = initial_batch_size
    i = 0   # downloaded records counter
    j = 0   # number of downloaded files

    # expected number of total files to be downloaded
    M = N // batch_size + (1 if N % batch_size else 0)

    entries = []
    while i < N:
        while True:
            repeat_download = False
            ids = ','.join(uniprot_accessions[i: i + batch_size])
            src = db_fetch_url % ids
            fn = '%s_uniprot_batch_%s.txt' % (dataset_prefix, _get_ids_hash(ids))
            dst = os.path.join(download_dir, fn)

            if not os.path.exists(dst):
                print "    Downloading file: %s (%d/%d)" % (fn, j+1, M)
                try:
                    dst = download_file(src, dst)
                    time.sleep(sleep_interval)  # wait a little
                except Exception:
                    print "      FAILURE"
                    if batch_size > 1:
                        batch_size = batch_size // 2
                        print "      Halving batch size to %d." % batch_size
                        repeat_download = True
                    else:
                        print "      Skipping %s." % ids

            if not repeat_download:
                entries += _process_uniprot_file(dst)
                i += batch_size
                batch_size = initial_batch_size
                r = 1 if (N-i) % batch_size else 0
                j += 1
                M = j + (N-i) // batch_size + r
                break

    acc_map = {}
    for ac, gn, orgn, taxid, pe, geneids in entries:
        for acc in ac:
            if acc in uniprot_accessions:
                acc_map[acc] = (gn, orgn, taxid, pe, geneids)
    return acc_map


def _process_uniprot_file(uniprot_file):
    """
    Extract gene name, organism and protein existence code from a Uniprot
    file.
    """
    uniprot_fp = open(uniprot_file, 'rU')
    old_tag = None
    joined_lines = None
    geneids = []
    ac = []
    gn = []
    orgn = '-'
    taxid = -1
    pe = '-'
    entries = []
    for line in uniprot_fp:
        tag = line[:2]
        if tag == '//':
            entries.append((ac, gn, orgn, taxid, pe, geneids))
            old_tag = joined_lines = None
            geneids = []
            ac = []
            gn = []
            orgn = pe = '-'
            continue

        if tag == 'DR':
            mtch = geneid_ptrn.match(line[5:-1])
            if mtch is not None:
                geneids.append(int(mtch.group(1)))
            continue

        if old_tag is not None:
            if tag != old_tag:
                if old_tag == 'AC':
                    ac = []
                    for mtch in acc_ptrn.finditer(joined_lines):
                        ac.append(mtch.group(1))
                elif old_tag == 'GN':
                    gn = []
                    for mtch in gn_ptrn.finditer(joined_lines):
                        gn.append(mtch.group(1))
                    for mtch in gn_ptrn2.finditer(joined_lines):
                        gn.extend(mtch.group(1).split(', '))
                elif old_tag == 'OS':
                    orgn = joined_lines.split('.')[0].split('(')[0]
                elif old_tag == 'OX':
                    mtch = ox_ptrn.match(joined_lines)
                    assert mtch is not None
                    taxid = int(mtch.group(1))
                elif old_tag == 'PE':
                    mtch = pe_ptrn.match(joined_lines)
                    if mtch is not None:
                        pe = mtch.group(1)
            else:
                joined_lines += line[5:-1]
                continue
        if tag not in ('AC', 'GN', 'OS', 'OX', 'PE'):
            tag = None
            joined_lines = None
        else:
            joined_lines = line[5:-1]
        old_tag = tag
    uniprot_fp.close()

    return entries


def geneids_from_symbols(gns_records, download_dir):
    """
    Retrieve Gene IDs from gene_info using a potential gene
    symbol
    """

    # each item in gns_records |-> (rogid, acc, gns, orgn, taxid, pe)

    gene_info = os.path.join(download_dir, 'gene_info')
    download_and_extract_file(gene_info_url, gene_info)
    gene_history = os.path.join(download_dir, 'gene_history')
    download_and_extract_file(gene_history_url, gene_history)

    # Create an 'index' of query symbols, then traverse gene_info once to try
    # and match canonical symbol to query symbols from Uniprot.
    # Index here is a simple hash of (orgn, gene_name) but the implementation
    # can be changed later if necessary.
    hits = dict()
    gn_index = dict()
    rogid2gns = dict()
    for rogid, _, gns, _, taxid, _ in gns_records:
        rogid2gns[rogid] = gns
        for gene_name in gns:
            assert gene_name != '-'
            if (taxid, gene_name) not in gn_index:
                gn_index[(taxid, gene_name)] = []
            gn_index[(taxid, gene_name)].append(rogid)


    gene_info_fp = open(gene_info, 'rU')
    gene_info_fp.next()  # skip first line
    for line in gene_info_fp:
        fields = line.strip().split('\t')
        tax_id, gene_id = map(int, fields[0:2])
        symbol, locus_tag = fields[2:4]

        gene_name = None
        if (tax_id, symbol) in gn_index:
            gene_name = symbol
            code = 'S'
        elif (tax_id, locus_tag) in gn_index:
            gene_name = locus_tag
            code = 'L'
        rogids = gn_index.pop((tax_id, gene_name), [])
        for rogid in rogids:
            hits[rogid] = (gene_id, code)
            # Match only one potential symbol
            for gene_name in rogid2gns[rogid]:
                gn_index.pop((tax_id, gene_name), None)

    gene_info_fp.close()

    # Also check for obsolete Gene IDs
    history_fp = open(gene_history, 'rU')
    history_fp.next()  # skip first line
    for line in history_fp:
        fields = line.strip().split('\t')
        tax_id = int(fields[0])
        new_gene_id = fields[1]
        old_gene_id = int(fields[2])
        symbol = fields[3]

        if (tax_id, symbol) in gn_index:
            rogids = gn_index.pop((tax_id, symbol))
            for rogid in rogids:
                if new_gene_id == '-':
                    hits[rogid] = (old_gene_id, 'B')
                else:
                    hits[rogid] = (int(new_gene_id), 'S')
                for gene_name in rogid2gns[rogid]:
                    gn_index.pop((tax_id, gene_name), None)
    history_fp.close()

    gns_mapping = []
    for rogid, acc, gns, orgn, taxid, pe in gns_records:
        if rogid in hits:
            gene_id, code = hits[rogid]
            geneids = [gene_id]
        else:
            code = 'U'
            geneids = None
        gns_mapping.append((rogid, acc, gns, orgn, pe, geneids, code))
    return gns_mapping


def gene_ids_from_uniprot_files(GRgraph, orphan_ids, download_dir,
                                mapping_logfile, dataset_prefix):
    """
    Retrieve Gene IDs using values of various field within a Uniprot file.

    Outputs tab-delimited logfile with the following columns:
    1. Uniprot ID
    2. Species
    3. Uniprot protein existence code (PE line)
    4. Uniprot canonical gene name (from GN line)
    5. NCBI Gene ID(s) assigned - '|' separated if multiple
    6. Processing code:
         G - matched to Gene ID using dbxref field in Uniprot record
         S - matched to Gene ID by using gene name from Uniprot record
               to match canonical gene symbol
         L - matched to Gene ID by using gene name from Uniprot record
               to match locus tag
         B - matched to an obsolete gene record - later removed
         U - could not match using gene name
         N - Uniprot record does not contain gene name
         F - Could not retrieve Uniprot record from www.uniprot.org
    """

    mapping_details = []
    resolved_orphans = set()

    uniprot_accessions = sorted(acc for rogid in orphan_ids
                                for db, acc in orphan_ids[rogid]
                                if db == 'uniprotkb')
    acc_map = _batch_retrieve_uniprot_files(uniprot_accessions, download_dir,
                                            dataset_prefix)

    print "    Mapping Uniprot IDs to Gene IDs using dbxref annotation."
    gns_records = []
    for rogid in orphan_ids:
        for db, acc in orphan_ids[rogid]:
            if db != 'uniprotkb':
                continue

            if acc not in acc_map:
                # Code F: failed to retrieve Uniprot file
                mapping_details.append((acc, [], '-', '-', '-', 'F'))
                continue

            gns, orgn, taxid, pe, geneids = acc_map[acc]
            if len(geneids):
                # Code G: direct match to gene id from Uniprot db xref
                mapping_details.append((acc, geneids, gns, orgn, pe, 'G'))
                for geneid in geneids:
                    _add_GR_link(GRgraph, geneid, rogid, 'G')
                resolved_orphans.add(rogid)

            elif len(gns):
                gns_records.append((rogid, acc, gns, orgn, taxid, pe))

            else:
                # Code N: no gene name found
                mapping_details.append((acc, [], [], orgn, pe, 'N'))


    print "    Mapping Uniprot IDs to Gene IDs using symbols."
    gns_mapping = geneids_from_symbols(gns_records, download_dir)
    for rogid, acc, gns, orgn, pe, geneids, code in gns_mapping:
        if geneids is not None:
            # Code S or L: succesfully matched gene id using gene name
            # Code B: matched obsolete record
            mapping_details.append((acc, geneids, gns, orgn, pe, code))
            if code in ('S', 'L', 'B'):
                for geneid in geneids:
                    _add_GR_link(GRgraph, geneid, rogid, code)
            resolved_orphans.add(rogid)
        else:
            # Code U: could not match gene id using gene name
            mapping_details.append((acc, [], gns, orgn, pe, code))


    # Logfile is written at the end
    print "    Writing Uniprot map."
    output_fp = open(mapping_logfile, 'w')
    for acc, geneids, gns, orgn, pe, code in mapping_details:
        if len(gns) == 0:
            gns = ['-']
        if len(geneids) == 0:
            geneids = ['-']
        else:
            geneids = map(str, geneids)

        output_fp.write('%s\n' % '\t'.join([acc, orgn, pe, '|'.join(gns),
                                            '|'.join(geneids), code]))
    output_fp.close()
    return resolved_orphans


def validate_gene_ids(used_geneids, download_dir):
    """
    Validate all mapped Gene IDs using gene_history file.
    """

    dst = os.path.join(download_dir, 'gene_history')
    download_and_extract_file(gene_history_url, dst)

    # Traverse the history file multiple times to follow through replacement
    # ids. This is to ensure that we catch the cases where the replacement id
    # is also obsoleted.

    obsolete_geneids = {}
    unverified_geneids = used_geneids
    while len(unverified_geneids):
        replacement_ids = set()
        history_fp = open(dst, 'rU')
        history_fp.next()  # skip first line
        for line in history_fp:
            cols = line.split('\t')
            obs_geneid = int(cols[2])
            if obs_geneid in unverified_geneids:
                if cols[1] == '-':
                    new_geneid = None
                else:
                    new_geneid = int(cols[1])
                    replacement_ids.add(new_geneid)
                obsolete_geneids[obs_geneid] = new_geneid
        history_fp.close()
        unverified_geneids = set(geneid for geneid in replacement_ids
                                 if geneid not in obsolete_geneids)

    # Follow through the chains of obsolete ids to obtain the final
    # replacement id
    geneid_replacement = {}
    for obs_geneid, new_geneid in obsolete_geneids.iteritems():
        while new_geneid in obsolete_geneids:
            new_geneid = obsolete_geneids[new_geneid]
        geneid_replacement[obs_geneid] = new_geneid
    return geneid_replacement


def get_gene_symbols(used_geneids, download_dir):
    """
    Assign symbols to mapped Gene IDs without them.
    """

    geneid2symbol = dict()
    gene_info = os.path.join(download_dir, 'gene_info')
    download_and_extract_file(gene_history_url, gene_info)
    gene_info_fp = open(gene_info, 'rU')
    gene_info_fp.next()  # skip first line
    for line in gene_info_fp:
        fields = line.strip().split('\t')
        geneid = int(fields[1])
        symbol = fields[2]
        if geneid in used_geneids:
            geneid2symbol[geneid] = str(symbol)
    gene_info_fp.close()
    return geneid2symbol


def construct_final_map(protein_ids, GRgraph, geneid_replacement, geneid2symbol):
    """
    Construct the final map protein_id -> (geneids, code) where all
    geneids must be valid.
    """

    # Replace invalid Gene IDs
    invalidated = []
    geneid_nodes = [id_node for id_node in GRgraph
                    if id_node not in protein_ids]
    for geneid in geneid_nodes:
        if geneid in geneid_replacement:
            new_geneid = geneid_replacement[geneid]
            if new_geneid is not None:
                for rogid, code in GRgraph[geneid].iteritems():
                    _add_GR_link(GRgraph, new_geneid, rogid, 'V')
                    del GRgraph[rogid][geneid]
            else:
                for rogid, code in GRgraph[geneid].iteritems():
                    invalidated.append((rogid, geneid, code))
                    del GRgraph[rogid][geneid]
            del GRgraph[geneid]

    # (Re-)Construct canonical groups
    canonical_groups = list()
    mapped_rogids = [u for u in GRgraph if u in protein_ids]
    grouped_rogids = set()
    for rogid in mapped_rogids:
        if rogid in grouped_rogids:
            continue
        unvisited = set([rogid])
        cgroup = set()
        while unvisited:
            u = unvisited.pop()
            cgroup.add(u)
            unvisited.update(v for v in GRgraph[u].iterkeys()
                             if v not in cgroup)

        grp_rogids = sorted(v for v in cgroup if v in protein_ids)
        grp_geneids = sorted(v for v in cgroup if v not in protein_ids)
        grp_symbols = [geneid2symbol[v] for v in grp_geneids]
        grp_codes = list()
        for rogid in grp_rogids:
            codes = sorted(GRgraph[rogid].itervalues())
            assert 'B' not in codes
            if not codes:
                # Here we could check invalidated and assign 'O' if any links
                # were 'I'. However, in this version, there are no such links.
                assert not grp_geneids
                codes = ['B']
            grp_codes.append(codes)
        canonical_groups.append((grp_rogids, grp_codes, grp_geneids,
                                 grp_symbols))
        grouped_rogids.update(grp_rogids)

    canonical_groups.sort(key=lambda grp: grp[0][0])
    return canonical_groups


def write_final_mappings(canonical_groups, id_logfile, icrogids):
    """
    Write final mappings to logfile.
    """

    id_lines = []
    mapping_counts = {}.fromkeys('IMGSLVOB', 0)

    for i, grp in enumerate(canonical_groups):
        rogids, map_codes, geneids, symbols = grp
        for rogid, codes in zip(rogids, map_codes):
            code = codes[0]
            id_lines.append((i+1, code, rogid, geneids, symbols,
                             icrogids[rogid]))
            mapping_counts[code] += 1

    fp = open(id_logfile, 'w')
    fp.write('#\n')
    fp.write('# ROGID MAPPINGS\n')
    fp.write('#\n')
    for line in id_lines:
        i, code, rogid, geneids, symbols, icrogid = line
        if len(geneids):
            written_geneids = ', '.join(map(str, geneids))
            written_symbols = ', '.join(symbols)
        else:
            written_geneids = '-'
            written_symbols = '-'
        fp.write('%d\t%s\t%s\t%s\t%s\t%d\n' % (i, rogid, code, written_geneids,
                                               written_symbols, icrogid))
    fp.close()

    mapping_table_line = [mapping_counts[k] for k in 'IVOMGSB']
    return mapping_table_line


def map_protein_ids_to_gene_ids(mitab_file, filtered_pmids_file, download_dir,
                                mapping_logfile, id_logfile, geneid_logfile,
                                mapping_file=None, accepted_taxids=None):
    """
    Master routine for mappings of rogids to gene_ids
    """

    dataset_prefix = os.path.split(mitab_file)[1].split('.')[0]
    filtered_pmids = read_filtered_pmids(filtered_pmids_file)
    counts = Phase1Counts()

    print "   Extracting protein IDs from MITAB file."
    protein_ids = extract_protein_ids(mitab_file, filtered_pmids,
                                      accepted_taxids)

    print "   Obtaining Gene IDs from mapping file."
    # mapping.txt could be large so we will store parsing results so we can
    # quickly get them if we need to debug
    mapping_data_path = os.path.join(download_dir,
                                     dataset_prefix + 'mapping.pkl')
    if not os.path.exists(mapping_data_path):
        mapping_data = gene_ids_from_mapping_txt(protein_ids, mapping_file)
        with open(mapping_data_path, 'wb') as fp:
            pickle.dump(mapping_data, fp, 2)
    else:
        with open(mapping_data_path, 'rb') as fp:
            mapping_data = pickle.load(fp)
    GRgraph, orphan_ids, icrogids = mapping_data

    print "   Obtaining Gene IDs from Uniprot records."
    uniprot_resolved = gene_ids_from_uniprot_files(GRgraph, orphan_ids,
                                                   download_dir,
                                                   mapping_logfile,
                                                   dataset_prefix)

    print "   Validating Gene IDs."
    used_geneids = set(id_node for id_node in GRgraph
                       if id_node not in protein_ids)
    geneid_replacement = validate_gene_ids(used_geneids, download_dir)

    print "   Retrieving Gene symbols."
    geneid2symbol = get_gene_symbols(used_geneids, download_dir)

    print "   Constructing canonical groups."
    canonical_groups = construct_final_map(protein_ids, GRgraph,
                                           geneid_replacement, geneid2symbol)

    print "   Writing final mappings."
    mapping_table_line = write_final_mappings(canonical_groups, id_logfile,
                                              icrogids)

    counts.total_rogids = len(protein_ids)
    counts.initial_mappings = len(protein_ids) - len(orphan_ids)
    counts.orphans = len(orphan_ids)
    for _oid in orphan_ids.iterkeys():
        _dbs = set(db for db, _ in orphan_ids[_oid])
        if 'uniprotkb' in _dbs:
            counts.uniprot_orphans += 1
        elif 'pdb' in _dbs:
            counts.pdb_orphans += 1
    counts.new_mappings = len(uniprot_resolved)

    initial_invalid = 0  # code 'O' does not occur now
    new_invalid = mapping_table_line[-1]  # code 'B'

    counts.initial_valid = counts.initial_mappings - initial_invalid
    counts.new_valid = counts.new_mappings - new_invalid
    counts.final_geneids = sum(1 for geneid in used_geneids
       if geneid_replacement.get(geneid, geneid) is not None)

    counts.final_rogids = sum(len(grp[0]) for grp in canonical_groups
                               if grp[2])

    lost = 100 * (1 - 1.0 * counts.final_rogids / counts.total_rogids)
    counts.lost_rogids = lost

    mapping_table_line[3:3] = [counts.orphans,
                               counts.pdb_orphans,
                               counts.uniprot_orphans]
    counts.mapping_table_line = mapping_table_line

    fp = open(id_logfile, 'a')
    counts.write_mappings(fp)
    fp.close()


def read_id_mapping(id_logfile):

    fp = open(id_logfile, 'rU')

    # Section header
    for i in xrange(3):
        line = fp.next()
        assert line[0] == '#'

    # First section is our mapping
    id_map = {}
    line = fp.next()
    while line[0] != '#':
        fields = line.strip().split('\t')
        protein_id = fields[1]
        code = fields[2]
        if code not in ('O', 'B'):
            geneids = map(int, fields[3].split(', '))
            symbols = fields[4].split(', ')
            mapped = (code != 'I')
            id_map[protein_id] = (geneids, symbols, mapped)
        line = fp.next()

    fp.close()
    return id_map
