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
Various statistics collected during each phase.
"""


class Phase1Counts(object):

    def __init__(self):

        self.total_rogids = 0
        self.initial_mappings = 0
        self.new_mappings = 0
        self.initial_valid = 0
        self.new_valid = 0
        self.orphans = 0
        self.uniprot_orphans = 0
        self.pdb_orphans = 0
        self.final_rogids = 0
        self.final_geneids = 0
        self.lost_rogids = 0.0

        self.raw_lines = 0
        self.noid_lines = 0
        self.filtered_lines = 0
        self.new_lines = 0
        self.saved_lines = 0

        self.raw_i = 0
        self.noid_i = 0
        self.filtered_i = 0
        self.new_i = 0
        self.saved_i = 0

        self.raw_rigids = set()
        self.noid_rigids = set()
        self.new_rigids = set()
        self.saved_rigids = set()

        self.not_mapped = set()

        self.mapping_table_line = []

    def write_mappings(self, fp):

        fp.write('#\n')
        fp.write('# INTERACTOR MAPPING STATISTICS\n')
        fp.write('#\n')
        fp.write('Total rogids: %d\n' % self.total_rogids)
        fp.write('Initially mapped to Gene IDs (total): %d\n'
                 % self.initial_mappings)
        fp.write('Initially mapped to Gene IDs (valid): %d\n'
                 % self.initial_valid)
        fp.write('Total Orphans: %d\n' % self.orphans)
        fp.write('Orphans with PDB IDs: %d\n' % self.pdb_orphans)
        fp.write('Orphans with Uniprot IDs: %d\n' % self.uniprot_orphans)

        fp.write('Additionally mapped to Gene IDs (total): %d\n' %
                 self.new_mappings)
        fp.write('Additionally mapped to Gene IDs (valid): %d\n' %
                 self.new_valid)
        fp.write('Final accepted rogids: %d\n' % self.final_rogids)
        fp.write('Final valid Gene IDs: %d\n' % self.final_geneids)
        fp.write('Lost rogids: %.1f %%\n' % self.lost_rogids)

        fp.write('\n****** TABLE 3 Lines (LaTeX)  ******\n\n')
        fp.write('& %d & %d & %d & %d & %d & %d & %d \\\\ \n'
                 % (self.total_rogids, self.initial_mappings,
                    self.orphans, self.new_mappings, self.new_valid,
                    self.final_rogids, self.final_geneids))

        fp.write('\n****** TABLE S-id_mappings (LaTeX)  ******\n\n')
        fp.write('%s \\\\ \n' % ' & '.join(map(str, self.mapping_table_line)))

    def to_file(self, fp):

        # This was only used for debugging
        if False:
            fp.write('****** SAVED RIGIDS  ******\n')
            for rigid in self.saved_rigids:
                fp.write(str(rigid))
                fp.write('\n')

        fp.write('****** PHASE 1 STATISTICS  ******\n')
        fp.write('**** Interactions ****\n')

        fp.write('Total lines scanned: %d\n' % self.raw_lines)
        fp.write('Filtered lines: %d\n' % self.filtered_lines)
        fp.write('Lines removed because of missing Gene IDs: %d\n'
                 % self.noid_lines)
        fp.write('Total retained lines: %d\n' % self.new_lines)
        fp.write('Lines retained due to mapping: %d\n' % self.saved_lines)

        fp.write('\n')
        fp.write('Total raw interactions scanned: %d\n' % self.raw_i)
        fp.write('Filtered raw interactions: %d\n' % self.filtered_i)
        fp.write('Raw interactions removed because of missing Gene IDs: %d\n'
                 % self.noid_i)
        fp.write('Total retained raw interactions: %d\n' % self.new_i)
        fp.write('Raw Interactions retained due to mapping: %d\n'
                 % self.saved_i)

        fp.write('\n')
        fp.write('Total rigids scanned: %d\n' % len(self.raw_rigids))
        fp.write('Rigids removed because of missing Gene IDs: %d\n'
                 % len(self.noid_rigids))
        fp.write('Total retained rigids: %d\n' % len(self.new_rigids))
        fp.write('Rigids retained due to mapping: %d\n'
                 % len(self.saved_rigids))

        fp.write('\n****** TABLE 2 Lines (LaTeX)  ******\n\n')
        fp.write('& %d & %d & %d & %d & %d \\\\ \n' %
                 (self.raw_i, self.filtered_i, self.noid_i, self.new_i,
                  self.saved_i))

        fp.write('\n****** TABLE S-rigids Lines (LaTeX)  ******\n\n')
        fp.write('& %d & %d & %d & %d \\\\ \n' %
                 (len(self.raw_rigids), len(self.noid_rigids),
                  len(self.new_rigids), len(self.saved_rigids)))

    def add_raw(self, interaction):
        self.raw_lines += interaction.num_lines
        self.raw_i += 1
        self.raw_rigids.add(interaction.rigid)

    def add_filtered(self, interaction):
        self.filtered_lines += interaction.num_lines
        self.filtered_i += 1

    def add_noid(self, interaction):
        self.noid_lines += interaction.num_lines
        self.noid_i += 1
        self.noid_rigids.add(interaction.rigid)

    def add_saved(self, interaction):
        self.saved_lines += interaction.num_lines
        self.saved_i += 1
        self.saved_rigids.add(interaction.rigid)

    def add_new(self, interaction):
        self.new_lines += interaction.num_lines
        self.new_i += 1
        self.new_rigids.add(interaction.rigid)


class Phase2Counts(object):

    _counters = ['initial_pairs',
                 'unused_pairs',
                 'C',
                 'G',
                 'R',
                 'A',
                 'N',
                 'GR',
                 'AR',
                 'NR',
                 'pmids']

    def __init__(self):

        for attr in self._counters:
            setattr(self, attr, 0)

    def to_file(self, fp):
        fp.write('****** SUMMARY ******\n')
        for attr in self._counters:
            fp.write('%s: %d\n' % (attr, getattr(self, attr)))

        table_line = [self.pmids,
                      self.initial_pairs,
                      self.unused_pairs,
                      self.C,
                      self.G + self.GR,
                      self.R + self.AR + self.NR,
                      self.A,
                      self.N]
        fp.write('\n****** TABLE 3 Line (LaTeX)  ******\n\n')
        fp.write('& %s \\\\ \n' % ' & '.join(map(str, table_line)))


short_db_names = {'bind_translation': 'N',
                  'corum': 'C',
                  'mpact': 'A',
                  'mppi': 'P',
                  'ophid': 'O',
                  'bind': 'U',
                  'biogrid': 'B',
                  'dip': 'D',
                  'hprd': 'H',
                  'intact': 'I',
                  'mint': 'M',
                  'innatedb': 'E',
                  'mpi-imex': 'J',
                  'mpi-lit': 'L',
                  'matrixdb': 'X',
                  }


class Phase3Counts(object):

    _counters = ['biochem', 'other', 'complexes', 'directed', 'undirected']

    def __init__(self):

        for attr in self._counters:
            setattr(self, attr, 0)
        self.pmids = set()
        self.resolvable = {}
        self.unresolvable = {}

    def add_resolvable(self, key):
        if key not in self.resolvable:
            self.resolvable[key] = 0
        self.resolvable[key] += 1

    def add_unresolvable(self, key, dbs1, dbs2):
        if key not in self.unresolvable:
            self.unresolvable[key] = [set(), set(), 0]
        self.unresolvable[key][0].update(dbs1)
        self.unresolvable[key][1].update(dbs2)
        self.unresolvable[key][2] += 1

    @staticmethod
    def _expand_term(term_id, ontology):
        term = ontology.get_term(term_id)
        fmt = '%s (%s)'
        return fmt % (term.term_id, term.name)

    @staticmethod
    def _abbrv_db(dbs):
        short = sorted(short_db_names[db] for db in dbs)
        return ''.join(short)

    def _res_cnflct_key(self, old_key, ontology):
        term2dbs = {}
        for term_id, db in old_key:
            if term_id not in term2dbs:
                term2dbs[term_id] = set()
            term2dbs[term_id].add(db)

        key = [(self._expand_term(term_id, ontology), self._abbrv_db(dbs))
               for term_id, dbs in term2dbs.iteritems()]
        return key

    def to_file(self, fp, ontology):

        fp.write('\n\n****** RESOLVABLE CONFLICTS SUMMARY ******\n')
        # For easier reference it might be better if we don't carry the
        # information about databases involved. So we will combine the stats
        # here before writing the table
        cmb_resolvable = {}
        for old_key, val in self.resolvable.iteritems():
            new_key = frozenset(term_id for term_id, _ in old_key)
            if new_key not in cmb_resolvable:
                cmb_resolvable[new_key] = 0
            cmb_resolvable[new_key] += val

        items = sorted((-val, new_key)
                       for new_key, val in cmb_resolvable.iteritems())
        for val, new_key in items:
            expanded_terms = [self._expand_term(term_id, ontology)
                              for term_id in sorted(new_key)]
            fp.write('%s & %d \\tabularnewline \n' %
                     (', '.join(expanded_terms), -val))

        fp.write('\n\n****** UNRESOLVABLE CONFLICTS SUMMARY ******\n')
        items = ((-val[2],
                  key[0], ''.join(self._abbrv_db(val[0])),
                  key[1], ''.join(self._abbrv_db(val[1])),
                  val[2])
                 for key, val in self.unresolvable.iteritems())
        for item in sorted(items):
            new_item = item[1:]
            fp.write('%s & %s & %s & %s & %d \\\\ \n' % new_item)

        num_conflicts = sum(v[2] for v in self.unresolvable.itervalues())
        res_conflicts = sum(self.resolvable.itervalues())

        fp.write('\n\n****** MAIN SUMMARY ******\n')
        fp.write('pmids: %d\n' % len(self.pmids))
        for attr in self._counters:
            fp.write('%s: %d\n' % (attr, getattr(self, attr)))
        fp.write('conflicts: %d\n' % num_conflicts)

        table_line = [len(self.pmids),
                      self.biochem,
                      self.other,
                      self.complexes,
                      self.directed,
                      self.undirected,
                      res_conflicts,
                      num_conflicts]
        fp.write('\n****** TABLE 4 Line (LaTeX)  ******\n\n')
        fp.write('& %s \\\\ \n' % ' & '.join(map(str, table_line)))
