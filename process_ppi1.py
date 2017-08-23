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
Divide interactions into disjoint groups for further processing.
"""

import sys
import obo
from parser import parse_mitab_file
from parser import full_mitab_iterator
from idmapping import read_filtered_pmids
from idmapping import read_id_mapping
from filter import is_filtered
from filter import NullFile
from counts import Phase1Counts
from interaction import Interaction
from interaction import iRefIndexInteraction


def _fix_detection_method(interaction, ontology):

    obsolete_detection = {'MI:0021': ('MI:0428', 'MI:0403'),
                          'MI:0022': ('MI:0428', 'MI:0403'),
                          'MI:0023': ('MI:0428', 'MI:0403'),
                          'MI:0025': ('MI:0401', 'MI:0914'),
                          'MI:0044': ('MI:0411', 'MI:0407'),  # elisa
                          'MI:0059': ('MI:0096', 'MI:0915'),
                          'MI:0061': ('MI:0096', 'MI:0915'),
                          # This one is not obsolete but here MI:0218 should
                          # not map to MI:0915 but to MI:0407
                          'MI:0065': ('MI:0065', 'MI:0407'),
                          'MI:0079': ('MI:0401', 'MI:0914'),
                          'MI:0109': ('MI:0676', 'MI:0915'),
                          'MI:0229': ('MI:0809', 'MI:0915'),  # bifc
                          # this should not exist apart from HPRD
                          'MI:0493': ('MI:0045', None),
                          }
    mi0000_detection = {
        'affinMI 0030- cross-linking studies': 'MI:0030',
        'affinity chromatography technologies': 'MI:0004',
        'alfa-screen': 'MI:0905',
        'cence resonance energy transfer': 'MI:0012',
        'cent resonance energy transfer': 'MI:0055',
        'chromatography technologies': 'MI:0091',
        'cosedimentation through density gradients': 'MI:0029',
        'coip  coimmunoprecipitation': 'MI:0019',
        'elisa  enzyme-linked immunosorbent assay': 'MI:0411',
        'enzymatic studies': 'MI:0415',
        'fluorescence technologies': 'MI:0051',
        'gdp_gtp exchange': 'MI:0949',
        'gst pull down': 'MI:0096',
        'hat': 'MI:0887',
        'in gel phosphatase assay': 'MI:0514',
        'linked immunosorbent assay': 'MI:0411',
        'OPHID Predicted Protein Interaction': 'MI:0063',
        'p elisa': 'MI:0813',
        'prediction': 'MI:0063',
        'predictive text mining': 'MI:0087',
        'protease access': 'MI:0814',
        'text mining': 'MI:0110',
        }

    old_dm_id = interaction.detection_method.term_id
    new_dm_id = None
    new_itype_id = None

    if old_dm_id is None or old_dm_id == 'MI:0001':
        # Totally unspecified - be nice and map to MI:0686
        # Note: if itype is missing, it will be assumed to be
        # MI:0914, which is probably correct
        new_dm_id = 'MI:0686'

    elif interaction.source_db.term_id == 'MI:0468':
        # Treat HPRD as a special case - transform MI:0045 back to (obsolete)
        # MI:0492 (in vitro). Also, do not transform MI:0493 if it comes from
        # HPRD.
        if interaction.detection_method.term_id == 'MI:0045':
            new_dm_id = 'MI:0492'

    elif old_dm_id == 'MI:0000':
        if interaction.detection_method.name in mi0000_detection:
            new_dm_id = mi0000_detection[interaction.detection_method.name]
        else:
            new_dm_id = 'MI:0686'  # unknown

    elif old_dm_id in obsolete_detection:
        new_dm_id, new_itype_id = obsolete_detection[old_dm_id]

    # Set new ids
    if new_dm_id is not None:
        term = ontology.get_term(new_dm_id)
        interaction.detection_method.set(term.term_id, term.name)
    if new_itype_id is not None:
        term = ontology.get_term(new_itype_id)
        interaction.interaction_type.set(term.term_id, term.name)


def _fix_itype(interaction, ontology):

    obsolete_itype = {'MI:0191': 'MI:0915',
                      'MI:0218': 'MI:0915',
                      }

    # We try to infer interaction type from detection method labels here. This
    # association is not always unique so there could be possible problems. We
    # try to be reasonably generous based on what other databases commonly use.
    dm_match = {'MI:0001': 'MI:0190',  # TOTALLY GENERIC - SHOULD NOT HAPPEN
                'MI:0004': 'MI:0915',  # affinity chrom
                'MI:0006': 'MI:0915',  # anti bait coip
                'MI:0007': 'MI:0915',  # anti tag coip
                'MI:0008': 'MI:0915',  # array technology
                'MI:0012': 'MI:0915',  # bret
                'MI:0013': 'MI:0915',  # biophysical
                'MI:0014': 'MI:0915',  # adenylate cyclase complementation
                'MI:0016': 'MI:0407',  # circular dichroism
                'MI:0017': 'MI:0407',  # classical fluorescence spectroscopy
                'MI:0018': 'MI:0407',  # two hybrid
                'MI:0019': 'MI:0915',  # coip
                'MI:0020': 'MI:0403',  # transmission electron microscopy
                'MI:0024': 'MI:0915',  # confirmational text mining
                'MI:0027': 'MI:0403',  # cosedimentation
                'MI:0028': 'MI:0403',  # solution sedimentation
                'MI:0029': 'MI:0403',  # density sedimentation
                'MI:0030': 'MI:0915',  # cross-linking
                'MI:0031': 'MI:0915',  # protein crosslink
                'MI:0038': 'MI:0407',  # dynamic light scattering
                'MI:0040': 'MI:0915',  # electron microscopy
                'MI:0042': 'MI:0407',  # epr
                'MI:0043': 'MI:0407',  # electron resonance
                'MI:0045': 'MI:0914',  # experimental interaction detection
                'MI:0046': 'MI:0914',  # experimental knowledge based
                'MI:0047': 'MI:0407',  # far western blotting
                'MI:0049': 'MI:0915',  # filter binding
                'MI:0051': 'MI:0407',  # fluorescence
                'MI:0053': 'MI:0407',  # fps
                'MI:0054': 'MI:0915',  # fluorescence-activated cell sorting
                'MI:0055': 'MI:0407',  # fluorescent resonance energy transfer
                'MI:0065': 'MI:0407',  # itc
                'MI:0063': 'MI:0915',  # interaction prediction
                'MI:0067': 'MI:0915',  # light scattering
                'MI:0069': 'MI:0915',  # ms of complexes
                'MI:0071': 'MI:0915',  # molecular sieving
                'MI:0077': 'MI:0407',  # nmr
                'MI:0081': 'MI:0407',  # peptide array
                'MI:0084': 'MI:0915',  # phage display
                'MI:0087': 'MI:0915',  # predictive text mining
                'MI:0089': 'MI:0915',  # protein array
                'MI:0090': 'MI:0407',  # protein complementation assay
                'MI:0091': 'MI:0915',  # chromatography technology
                'MI:0095': 'MI:0914',  # seldi chip
                'MI:0096': 'MI:0915',  # pull down
                'MI:0099': 'MI:0915',  # spa
                'MI:0104': 'MI:0915',  # static light scattering
                'MI:0107': 'MI:0407',  # spr
                'MI:0108': 'MI:0915',  # t7 phage display
                'MI:0110': 'MI:0915',  # text mining
                'MI:0112': 'MI:0915',  # ubiquitin reconstruction
                'MI:0113': 'MI:0915',  # immunoblotting
                'MI:0114': 'MI:0407',  # x-ray diffraction
                'MI:0225': 'MI:0914',  # chromatin immunoprecipitation array
                'MI:0226': 'MI:0915',  # ion exchange chrom
                'MI:0227': 'MI:0915',  # reverse phase chromatography
                'MI:0232': 'MI:0915',  # transcriptional complementation assay
                'MI:0254': 'MI:0931',  # genetic interference
                'MI:0276': 'MI:0915',  # bn-page
                'MI:0363': 'MI:0915',  # inferred by author
                'MI:0370': 'MI:0915',  # toxcat
                'MI:0397': 'MI:0915',  # two hybrid array
                'MI:0398': 'MI:0915',  # two hybrid pooling
                'MI:0399': 'MI:0915',  # two hybrid fragment pooling
                'MI:0400': 'MI:0915',  # affinity technology
                'MI:0401': 'MI:0914',  # biochemical
                'MI:0402': 'MI:0914',  # chromatin immunoprecipitation assays
                'MI:0404': 'MI:0915',  # comig non denat gel
                'MI:0405': 'MI:0915',  # competition binding
                'MI:0406': 'MI:0197',  # deacetylase assay
                'MI:0410': 'MI:0915',  # electron tomography
                'MI:0411': 'MI:0407',  # elisa
                'MI:0412': 'MI:0915',  # emsa supershift
                'MI:0413': 'MI:0914',  # emsa
                'MI:0415': 'MI:0414',  # enzymatic study
                'MI:0416': 'MI:0403',  # fluorescence imaging
                'MI:0417': 'MI:0915',  # footprinting
                'MI:0423': 'MI:0217',  # in-gel kinase assay
                'MI:0424': 'MI:0217',  # phosphorylation
                'MI:0426': 'MI:0403',  # light microscopy
                'MI:0428': 'MI:0403',  # imaging technique
                'MI:0430': 'MI:0915',  # nucleic acid uv cross-linking
                'MI:0434': 'MI:0203',  # dephosphorylation
                'MI:0435': 'MI:0570',  # protease assay
                'MI:0437': 'MI:0914',  # protein tri hybrid
                'MI:0440': 'MI:0915',  # saturation binding
                'MI:0492': 'MI:0915',  # HPRD in-vivo
                'MI:0493': 'MI:0915',  # HPRD in-vitro
                'MI:0514': 'MI:0203',  # in gel phosphatase assay
                'MI:0515': 'MI:0213',  # methyltransferase assay
                'MI:0516': 'MI:0213',  # methyltransferase radiometric assay
                'MI:0655': 'MI:0407',  # lambda repressor two hybrid
                'MI:0663': 'MI:0403',  # confocal microscopy
                'MI:0676': 'MI:0915',  # tap
                'MI:0678': 'MI:0915',  # antibody array
                'MI:0686': 'MI:0914',  # UNSPECIFIED DETECTION METHOD
                'MI:0729': 'MI:0915',  # lumier
                'MI:0807': 'MI:0915',  # comigration in gel electrophoresis
                'MI:0808': 'MI:0915',  # comigration in sds page
                'MI:0809': 'MI:0915',  # bifc
                'MI:0813': 'MI:0407',  # proximity elisa
                'MI:0814': 'MI:0914',  # protease accessibility laddering
                'MI:0825': 'MI:0407',  # x-ray fiber diffraction
                'MI:0826': 'MI:0915',  # x ray scattering
                'MI:0833': 'MI:0914',  # autoradiography
                'MI:0841': 'MI:0844',  # phosphotransfer assay
                'MI:0859': 'MI:0407',  # intermolecular force
                'MI:0870': 'MI:0871',  # demethylase assay
                'MI:0872': 'MI:0407',  # atomic force microscopy
                'MI:0880': 'MI:0414',  # atpase assay
                'MI:0888': 'MI:0407',  # small angle neutron scattering
                'MI:0889': 'MI:0192',  # acetylation assay
                'MI:0892': 'MI:0915',  # solid phase assay
                'MI:0894': 'MI:0915',  # electron diffraction
                'MI:0899': 'MI:0915',  # p3 filamentous phage display
                'MI:0905': 'MI:0407',  # alfa-screen
                'MI:0921': 'MI:0407',  # surface-plasmon-resonance-chip
                'MI:0943': 'MI:0915',  # detection by mass spectrometry
                'MI:0947': 'MI:0915',  # bead aggregation assay
                'MI:0949': 'MI:0414',  # gdp/gtp exchange assay
                }

    old_itype_id = interaction.interaction_type.term_id
    new_itype_id = None
    if old_itype_id is None:
        # Try to fix missing interaction type based on detection method
        dm_id = interaction.detection_method.term_id
        new_itype_id = dm_match.get(dm_id, 'MI:0190')

    elif old_itype_id in obsolete_itype:
        # Replace obsolete
        new_itype_id = obsolete_itype[old_itype_id]
    else:
        return

    assert new_itype_id is not None
    term = ontology.get_term(new_itype_id)
    interaction.interaction_type.set(term.term_id, term.name)


def _fix_biological_role(interaction, ontology):
    mi0000_map = {'phosphate acceptor': 'MI:0843',
                  'phosphate donor': 'MI:0842',
                  'unspecified': 'MI:0499',
                  'competitor': 'MI:0941',
                  }

    for p in interaction.interactors:
        if p.biological_role.term_id == 'MI:0000':
            new_term_id = None
            if p.biological_role.name in mi0000_map:
                new_term_id = mi0000_map[p.biological_role.name]
                term = ontology.get_term(new_term_id)
                p.biological_role.set(term.term_id, term.name)
            assert new_term_id is not None


def _fix_experimental_role(interaction, ontology):

    mi0000_map = {'unspecified': 'MI:0499',
                  'fret pair': 'MI:0865',
                  'fluorescence accept': 'MI:0584',
                  'unspecified role': 'MI:0499',
                  'bait': 'MI:0496',
                  'prey': 'MI:0497',
                  'experimental role': 'MI:0499'  # assumed unspecified
                  }

    for p in interaction.interactors:
        if p.experimental_role.term_id == 'MI:0000':
            new_term_id = None
            names = p.experimental_role.name.split('|')
            for name in names:
                if name in mi0000_map:
                    new_term_id = mi0000_map[name]
                    term = ontology.get_term(new_term_id)
                    p.experimental_role.set(term.term_id, term.name)
                    break
            assert new_term_id is not None


class _BiochemFilter(object):
    """
    A routine for selecting biochemical reaction (with static local variables)
    """

    def __init__(self, ontology, biogrid_ptm_codes_file):
        self.ontology = ontology
        self.biogrid_ptms = self._read_ptm_codes(biogrid_ptm_codes_file)
        self.cache = {}
        self.enz_reaction_term = ontology.get_term('MI:0414')
        self.enz_study_term = ontology.get_term('MI:0415')

    def _get_biogrid_interaction_id(self, interaction):

        biogrid_ids = [int(item.acc) for item in interaction.ids.ids
                       if item.db == 'biogrid']
        if len(biogrid_ids):
            return biogrid_ids[0]
        return None

    def _read_ptm_codes(self, biogrid_ptm_codes_file):
        fp = open(biogrid_ptm_codes_file, 'rb')
        biogrid_ptms = {}
        for line in fp:
            cols = line.strip().split('\t')
            biogrid_ptms[int(cols[0])] = (cols[1], cols[2])
        fp.close()
        return biogrid_ptms

    def __call__(self, interaction):

        # Check if detection method is subterm of MI:0415(enzymatic study)
        tmp_id = interaction.detection_method.term_id
        if tmp_id is not None:
            term = self.ontology.get_term(tmp_id)
            val = term.compare_to(self.enz_study_term, self.cache)
            if val >= 0:
                if interaction.source_db.term_id == 'MI:0463':
                    # BioGRID - replace codes
                    i = self._get_biogrid_interaction_id(interaction)
                    if i in self.biogrid_ptms:
                        interaction.interaction_type.set(*self.biogrid_ptms[i])
                return True

        # Check if interaction type is subterm of MI:0414(enzymatic reaction)
        tmp_id = interaction.interaction_type.term_id
        if tmp_id is not None:
            term = self.ontology.get_term(tmp_id)
            val = term.compare_to(self.enz_reaction_term, self.cache)
            if val >= 0:
                return True

        # Check if biological role for any interactor is enzyme or enzyme
        # target
        for p in interaction.interactors:
            if p.biological_role.term_id in ('MI:0501', 'MI:0502'):
                return True

        return False


class _ComplexFilter(object):
    """
    A routine for selecting candidate complex parts (with static local
    variables)
    """

    def __init__(self, ontology):
        self.ontology = ontology
        self.cache = {}
        self.term0004 = ontology.get_term('MI:0004')
        self.term0403 = ontology.get_term('MI:0403')
        self.term0407 = ontology.get_term('MI:0407')

    def __call__(self, interaction):

        # Check if interaction labelled as true complex
        if interaction.is_complex():
            return True

        # Get terms for detection method and interaction type
        dterm_id = interaction.detection_method.term_id
        if dterm_id is not None:
            dterm = self.ontology.get_term(dterm_id)
        else:
            dterm = None
        iterm_id = interaction.interaction_type.term_id
        if iterm_id is not None:
            iterm = self.ontology.get_term(iterm_id)
        else:
            iterm = None

        if iterm is not None:
            # Don't accept anything labelled with 'direct interaction'
            # (MI:0407)
            if iterm.compare_to(self.term0407, self.cache) >= 0:
                return False
            # Colocalization (MI:0403) is accepted
            if iterm.compare_to(self.term0403, self.cache) >= 0:
                return True

        if dterm is not None:
            # Affinity chromatogrphy (MI:004) is accepted
            if dterm.compare_to(self.term0004, self.cache) >= 0:
                return True
            # BioGRID's 'Co-purification' is also accepted
            if iterm is not None and \
               dterm.term_id == 'MI:0401' and \
               iterm.term_id == 'MI:0914':
                return True

        return False


# Filter results
REMOVE = 0
BINARY = 1
COMPLEX = 2
BIOCHEM = 3


def _process_interaction(interaction, id_map, filtered_pmids, logfile_fp,
                         counts, lines, ontology, biochem_filter,
                         complex_filter, accepted_taxids):

    assert interaction.rigid is not None
    counts.add_raw(interaction)

    if is_filtered(interaction, filtered_pmids, logfile_fp, lines,
                   accepted_taxids):
        counts.add_filtered(interaction)
        return REMOVE

    # Get Gene IDs
    geneids_lst = [id_map.get(p.id, ([], [], False))
                   for p in interaction.interactors]
    saved = False  # interaction 'saved' (not filtered) due to Uniprot mapping
    for i, tmp in enumerate(geneids_lst):
        gene_ids, gene_symbols, mapped = tmp
        if len(gene_ids) == 0:
            logfile_fp.write('At least one interactor cannot be mapped into'
                             ' a valid NCBI Gene ID\n')

            non_mapped = interaction.interactors[i].id
            if non_mapped is None:
                non_mapped = -1
            logfile_fp.write('   Non-mapped protein id: %s\n' % non_mapped)
            logfile_fp.write('   Removed line(s): %s.\n' %
                             ', '.join('%d' % k for k in lines))
            counts.add_noid(interaction)
            return REMOVE
        if mapped:
            saved = True
    if saved:
        counts.add_saved(interaction)

    for p, tmp in zip(interaction.interactors, geneids_lst):
        gene_ids, gene_symbols = tmp[:2]
        p.normalize_ids_with_gene(gene_ids, gene_symbols)
    if interaction.is_complex():
        interaction.complex.normalize_ids_with_gene(None, None)

    counts.add_new(interaction)

    # Replace invalid or obsolete terms for detection method, interaction type,
    # biological role and experimental role by the closest valid terms.
    _fix_detection_method(interaction, ontology)
    _fix_itype(interaction, ontology)
    _fix_biological_role(interaction, ontology)
    _fix_experimental_role(interaction, ontology)

    # Check if this is potentially a directed interaction
    if biochem_filter(interaction):
        return BIOCHEM

    # Check if this is a candidate for clustering into complex
    if complex_filter(interaction):
        return COMPLEX

    return BINARY


def process_ppi1(irefindex_file, id_logfile, output_logfile, biochem_file,
                 binary_file, complexes_file, obo_file, biogrid_ptm_codes_file,
                 filtered_pmids_file=None, accepted_taxids=None):

    counts = Phase1Counts()
    filtered_pmids = read_filtered_pmids(filtered_pmids_file)

    id_map = read_id_mapping(id_logfile)

    obo_fp = open(obo_file, 'rU')
    ontology = obo.OBOntology(obo_fp)

    biochem_filter = _BiochemFilter(ontology, biogrid_ptm_codes_file)
    complex_filter = _ComplexFilter(ontology)

    input_fp = open(irefindex_file, 'rU')
    removed_fp = NullFile()
    biochem_fp = open(biochem_file, 'w')
    binary_fp = open(binary_file, 'w')
    complex_fp = open(complexes_file, 'w')
    logfile_fp = open(output_logfile, 'w')

    output_fps = (removed_fp, binary_fp, complex_fp, biochem_fp)
    for fp in output_fps:
        Interaction.write_header(fp)

    scanner = parse_mitab_file(input_fp, full_mitab_iterator, None,
                               iRefIndexInteraction)
    for interaction, lines in scanner:
        line_numbers = lines[1]
        res = _process_interaction(interaction, id_map, filtered_pmids,
                                   logfile_fp, counts, line_numbers, ontology,
                                   biochem_filter, complex_filter,
                                   accepted_taxids)
        interaction.to_file(output_fps[res])

    counts.to_file(logfile_fp)
    input_fp.close()
    obo_fp.close()
    logfile_fp.close()

    for fp in output_fps:
        fp.close()


if __name__ == "__main__":

    args = sys.argv[1:]
    process_ppi1(*args)
