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
Classes for various field types from PSI-MI TAB 2.6 format.
"""

import re
import datetime


class EmptyField(object):
    def __str__(self):
        return '-'

    def set(self):
        pass

    def get(self):
        return tuple()

    def from_string(self, data):
        assert data == '-'


class PlainText(object):

    def __init__(self):
        self.txt = None

    def __str__(self):
        if self.txt is None:
            return '-'
        return "%s" % self.txt

    def set(self, txt):
        self.txt = txt

    def get(self):
        return self.txt,

    def from_string(self, data):
        if data == '-':
            self.txt = None
        else:
            self.txt = data


class SingleID(object):

    _id_pattern = re.compile('(.+):(.+)')

    def __init__(self):
        self.db = None
        self.acc = None

    def __str__(self):
        if self.db is None and self.acc is None:
            return '-'
        return "%s:%s" % (self.db, self.acc)

    def set(self, db, acc):
        self.db = db
        self.acc = acc

    def get(self):
        return self.db, self.acc

    def from_string(self, data):
        if data == '-':
            self.db = None
            self.acc = None
        else:
            mtch = self._id_pattern.match(data)
            if mtch is not None:
                self.db = mtch.group(1)
                self.acc = mtch.group(2)


class MultipleIDs(object):

    def __init__(self):
        self.ids = []

    def __str__(self):
        if len(self.ids) == 0:
            return '-'
        return '|'.join(item.__str__() for item in self.ids)

    def set(self, ids):
        self.ids = ids

    def get(self):
        return self.ids,

    def from_string(self, data):
        if data == '-':
            self.ids = []
        else:
            for t in data.split('|'):
                item = SingleID()
                item.from_string(t)
                self.ids.append(item)


class TaxonomyID(object):

    _taxid_pattern = re.compile('taxid:([0-9]+)\((.+)\)')

    def __init__(self):
        self.taxid = None
        self.species = None

    def __str__(self):
        if self.taxid is None:
            return '-'
        return "taxid:%d(%s)" % (self.taxid, self.species)

    def set(self, taxid, species):
        self.taxid = taxid
        self.species = species

    def get(self):
        return self.taxid, self.species

    def from_string(self, data):
        if data == '-':
            self.taxid = None
            self.species = None
        else:
            mtch = self._taxid_pattern.match(data)
            if mtch is not None:
                self.taxid = int(mtch.group(1))
                self.species = mtch.group(2)


class DateField(object):

    _date_pattern = re.compile('([0-9]+)/([0-9]+)/([0-9]+)')

    def __init__(self):
        self.date = None

    def __str__(self):
        if self.date is None:
            return '-'
        return self.date.strftime('%Y/%m/%d')

    def set(self, date):
        self.date = date

    def set_today(self):
        self.date = datetime.datetime.today()

    def get(self):
        return self.date,

    def from_string(self, data):
        if data == '-':
            self.date = None
        else:
            mtch = self._date_pattern.match(data)
            if mtch is not None:
                year = int(mtch.group(1))
                month = int(mtch.group(2))
                day = int(mtch.group(3))
                self.date = datetime.datetime(year, month, day)


class ConfidenceScores(MultipleIDs):

    def from_string(self, data):
        if data == '-':
            self.ids = []
        else:
            for t in data.split('|'):
                if ':' in t:
                    item = SingleID()
                    item.from_string(t)
                    self.ids.append(item)


class PSIMITerm(object):

    _cvterm_pattern = re.compile('(MI:[0-9]+)\((.+)\)')

    def __init__(self):
        self.term_id = None
        self.name = None

    def __str__(self):
        if self.is_null():
            return '-'
        return "%s(%s)" % (self.term_id, self.name)

    def is_null(self):
        return self.term_id is None and self.name is None

    def set(self, term_id, name):
        self.term_id = term_id
        self.name = name

    def get(self):
        return self.term_id, self.name

    def from_string(self, data):
        if data == '-':
            self.term_id = None
            self.name = None
        else:
            mtch = self._cvterm_pattern.match(data)
            if mtch is not None:
                self.term_id = mtch.group(1)
                self.name = mtch.group(2)

    def normalize_name(self, ontology):
        if self.term_id is not None:
            term = ontology.get_term(self.term_id)
            self.name = term.name


class MultipleItems(list):
    """
    A list of various items.
    """

    _item_class = None

    def __init__(self):
        list.__init__(self)

    def __str__(self):
        if len(self) == 0:
            return '-'
        return '|'.join(item.__str__() for item in self)

    def set(self, terms):
        self[:] = terms

    def get(self):
        return self[:],

    def from_string(self, data):
        if data == '-':
            self[:] = []
        else:
            for t in data.split('|'):
                item = self._item_class()
                item.from_string(t)
                self.append(item)

    def from_items(self, items):
        """
        Builds a list of items from a collection of items. Eliminates
        duplicates (by both term_id and name).
        """
        items_set = set(term.get() for term in items)
        for item_data in sorted(items_set):
            item = self._item_class()
            item.set(*item_data)
            self.append(item)

    def add(self, *item_data):

        item = self._item_class()
        item.set(*item_data)
        self.append(item)


class MultiplePSIMITerms(MultipleItems):
    """
    A list of PSI-MI Terms.
    """
    _item_class = PSIMITerm
