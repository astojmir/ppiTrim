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
Code for dumping ppiTrim's final results into a normalized relational database
"""

import sqlite3
from parser import parse_mitab_file
from parser import full_mitab_iterator


class SqliteLoader(object):

    sql_tables = \
"""
CREATE TABLE psimi(
  id          INTEGER PRIMARY KEY,
  accession   TEXT,
  name        TEXT
);

CREATE TABLE taxonomy(
  id          INTEGER PRIMARY KEY,
  species     TEXT
);

CREATE TABLE entities(
  id          INTEGER PRIMARY KEY,
  symbol      TEXT,
  organism    INTEGER,
  FOREIGN KEY(organism) REFERENCES taxonomy(id)
);

CREATE TABLE aliases(
  entityid    INTEGER,
  database    TEXT,
  accession   TEXT,
  FOREIGN KEY(entityid) REFERENCES entities(id)
);

CREATE TABLE interactions(
  id                    INTEGER PRIMARY KEY,
  checksum              TEXT,
  method                INTEGER,
  author                TEXT,
  pubmed                INTEGER,
  itype                 INTEGER,
  edgetype              TEXT,
  maxsources            INTEGER,
  dmconsistency         TEXT,
  hostorganism          INTEGER,
  creationdate          TEXT,
  updatedate            TEXT,
  negative              INTEGER,
  FOREIGN KEY(method) REFERENCES psimi(id),
  FOREIGN KEY(itype) REFERENCES psimi(id),
  FOREIGN KEY(hostorganism) REFERENCES taxonomy(id)
);

CREATE TABLE interactors(
  entityid              INTEGER,
  interactionid         INTEGER,
  interactor_type       INTEGER,
  FOREIGN KEY(entityid) REFERENCES entities(id),
  FOREIGN KEY(interactionid) REFERENCES interactions(id),
  FOREIGN KEY(interactor_type) REFERENCES psimi(id)
);

CREATE TABLE conflicts(
  interactionid1        INTEGER,
  interactionid2        INTEGER,
  FOREIGN KEY(interactionid1) REFERENCES interactions(id),
  FOREIGN KEY(interactionid2) REFERENCES interactions(id)
);

CREATE TABLE interaction_xrefs(
  interactionid         INTEGER,
  database              TEXT,
  accession             TEXT,
  FOREIGN KEY(interactionid) REFERENCES interactions(id)
);

CREATE TABLE interaction_sources(
  interactionid         INTEGER,
  termid                INTEGER,
  FOREIGN KEY(interactionid) REFERENCES interactions(id),
  FOREIGN KEY(termid) REFERENCES psimi(id)
);

CREATE TABLE biological_roles(
  entityid              INTEGER,
  interactionid         INTEGER,
  role                  INTEGER,
  FOREIGN KEY(entityid) REFERENCES entities(id),
  FOREIGN KEY(interactionid) REFERENCES interactions(id),
  FOREIGN KEY(role) REFERENCES psimi(id)
);

CREATE TABLE experimental_roles(
  entityid              INTEGER,
  interactionid         INTEGER,
  role                  INTEGER,
  FOREIGN KEY(entityid) REFERENCES entities(id),
  FOREIGN KEY(interactionid) REFERENCES interactions(id),
  FOREIGN KEY(role) REFERENCES psimi(id)
);


"""
    sql_placeholder = '?'

    def __init__(self):

        self.conn = None
        self.cur = None

        self.intr_counter = 0
        self.intrid_map = {}
        self.used_entities = set()
        self.used_terms = {}
        self.term_counter = -32768
        self.used_taxons = set()
        self.conflicts = []

        self.sql_psimi_insert = self.insert_stmt('psimi', 3)
        self.sql_taxonomy_insert = self.insert_stmt('taxonomy', 2)
        self.sql_entities_insert = self.insert_stmt('entities', 3)
        self.sql_aliases_insert = self.insert_stmt('aliases', 3)
        self.sql_interactions_insert = self.insert_stmt('interactions', 13)
        self.sql_interactors_insert = self.insert_stmt('interactors', 3)
        self.sql_conflicts_insert = self.insert_stmt('conflicts', 2)
        self.sql_ixrefs_insert = self.insert_stmt('interaction_xrefs', 3)
        self.sql_isources_insert = self.insert_stmt('interaction_sources', 2)
        self.sql_broles_insert = self.insert_stmt('biological_roles', 3)
        self.sql_eroles_insert = self.insert_stmt('experimental_roles', 3)

    def insert_stmt(self, table, num_vals, placeholder='?'):
        stmt = 'INSERT INTO %s  VALUES (%s)' % \
               (table,
                ', '.join([self.sql_placeholder] * num_vals))
        return stmt

    def connect_to_db(self, sqlite_file):

        self.conn = sqlite3.connect(sqlite_file)
        self.cur = self.conn.cursor()
        self.cur.executescript(self.sql_tables)

    def __call__(self, ppiTrim_file, *dbargs, **dbkwargs):

        self.connect_to_db(*dbargs, **dbkwargs)
        input_fp = open(ppiTrim_file, 'rU')
        scanner = parse_mitab_file(input_fp, full_mitab_iterator)
        for intr, _ in scanner:
            self._insert_interaction(intr)

        self._insert_all_conflicts()

        input_fp.close()
        self.cur.close()
        self.conn.commit()
        self.conn.close()

    def _insert_interaction(self, intr):

        self.intr_counter += 1
        self.intrid_map[intr.checksum.acc] = self.intr_counter

        # Store conflicts - inserted only at the end
        _ms, _dc, conflicts = self._confidence_scores(intr)
        for ppiTrim_id in conflicts:
            self.conflicts.append((intr.checksum.acc, ppiTrim_id))

        # Insert interaction data
        data = (self.intr_counter,
                intr.checksum.acc,
                self._insert_psimi_term(intr.detection_method),
                intr.author.txt,
                int(intr.publications.ids[0].acc),
                self._insert_psimi_term(intr.interaction_type),
                intr.edgetype,
                _ms,
                _dc,
                self._insert_taxon(intr.host_taxid),
                str(intr.creation_date),
                str(intr.update_date),
                0  # always false
                )
        self.cur.execute(self.sql_interactions_insert, data)

        # Insert source databases
        for term in intr.source_db:
            tid = self._insert_psimi_term(term)
            self.cur.execute(self.sql_isources_insert,
                             (self.intr_counter, tid))

        # Insert xrefs
        for xid in intr.ids.ids:
            db, acc = xid.get()
            if db == 'ppiTrim':
                continue
            self.cur.execute(self.sql_ixrefs_insert,
                             (self.intr_counter, db, acc))

        # Insert interactors
        for p in intr.interactors:
            self._insert_interactor(self.intr_counter, p)

    def _confidence_scores(self, intr):
        maxsources = None
        dmconsistency = None
        conflicts = []
        for item in intr.confidence.ids:
            if item.db == 'maxsources':
                maxsources = int(item.acc)
            elif item.db == 'dmconsistency':
                dmconsistency = item.acc
            elif item.db == 'conflicts':
                conflicts.extend(item.acc.split(','))
        return maxsources, dmconsistency, conflicts

    def _insert_all_conflicts(self):

        for id1, id2 in self.conflicts:
            self.cur.execute(self.sql_conflicts_insert,
                             (self.intrid_map[id1], self.intrid_map[id2]))

    def _insert_interactor(self, iid, p):

        # Principal ID is always Gene ID
        geneid = int(p.uid.acc)

        if geneid not in self.used_entities:
            # Symbol is the first entry in the alias field
            symbol = p.alias.ids[0].acc
            organism = self._insert_taxon(p.organism)
            self.cur.execute(self.sql_entities_insert,
                             (geneid, symbol, organism))

            for item in p.alt.ids[1:] + p.alias.ids[1:]:
                db, acc = item.get()
                self.cur.execute(self.sql_aliases_insert, (geneid, db, acc))

            self.used_entities.add(geneid)

        for term in p.biological_role:
            tid = self._insert_psimi_term(term)
            self.cur.execute(self.sql_broles_insert, (geneid, iid, tid))

        for term in p.experimental_role:
            tid = self._insert_psimi_term(term)
            self.cur.execute(self.sql_eroles_insert, (geneid, iid, tid))

        self.cur.execute(self.sql_interactors_insert,
                         (geneid, iid,
                          self._insert_psimi_term(p.interactor_type[0])))

    def _insert_psimi_term(self, term):

        if term.is_null():
            return None

        if str(term) not in self.used_terms:
            num_id = int(term.term_id.split(':')[1])
            if num_id == 0:
                num_id = self.term_counter
                self.term_counter -= 1
            self.cur.execute(self.sql_psimi_insert,
                             (num_id, term.term_id, term.name))
            self.used_terms[str(term)] = num_id
        else:
            num_id = self.used_terms[str(term)]
        return num_id

    def _insert_taxon(self, taxon):

        if taxon.taxid not in self.used_taxons:
            self.cur.execute(self.sql_taxonomy_insert,
                             (taxon.taxid, taxon.species))
            self.used_taxons.add(taxon.taxid)
        return taxon.taxid
