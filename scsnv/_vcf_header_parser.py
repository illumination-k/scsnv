import codecs
import collections
import csv
import gzip
import itertools
import os
import re
import sys

try:
    from collections import OrderedDict
except ImportError:
    from ordereddict import OrderedDict

try:
    import cparse
except ImportError:
    cparse = None

# Metadata parsers/constants
RESERVED_INFO = {
    "AA": "String",
    "AC": "Integer",
    "AF": "Float",
    "AN": "Integer",
    "BQ": "Float",
    "CIGAR": "String",
    "DB": "Flag",
    "DP": "Integer",
    "END": "Integer",
    "H2": "Flag",
    "H3": "Flag",
    "MQ": "Float",
    "MQ0": "Integer",
    "NS": "Integer",
    "SB": "String",
    "SOMATIC": "Flag",
    "VALIDATED": "Flag",
    "1000G": "Flag",
    # Keys used for structural variants
    "IMPRECISE": "Flag",
    "NOVEL": "Flag",
    "SVTYPE": "String",
    "SVLEN": "Integer",
    "CIPOS": "Integer",
    "CIEND": "Integer",
    "HOMLEN": "Integer",
    "HOMSEQ": "String",
    "BKPTID": "String",
    "MEINFO": "String",
    "METRANS": "String",
    "DGVID": "String",
    "DBVARID": "String",
    "DBRIPID": "String",
    "MATEID": "String",
    "PARID": "String",
    "EVENT": "String",
    "CILEN": "Integer",
    "DPADJ": "Integer",
    "CN": "Integer",
    "CNADJ": "Integer",
    "CICN": "Integer",
    "CICNADJ": "Integer",
}

RESERVED_FORMAT = {
    "GT": "String",
    "DP": "Integer",
    "FT": "String",
    "GL": "Float",
    "GLE": "String",
    "PL": "Integer",
    "GP": "Float",
    "GQ": "Integer",
    "HQ": "Integer",
    "PS": "Integer",
    "PQ": "Integer",
    "EC": "Integer",
    "MQ": "Integer",
    # Keys used for structural variants
    "CN": "Integer",
    "CNQ": "Float",
    "CNL": "Float",
    "NQ": "Integer",
    "HAP": "Integer",
    "AHAP": "Integer",
}

# Spec is a bit weak on which metadata lines are singular, like fileformat
# and which can have repeats, like contig
SINGULAR_METADATA = ["fileformat", "fileDate", "reference"]

# Conversion between value in file and Python value
field_counts = {
    ".": None,  # Unknown number of values
    "A": -1,  # Equal to the number of alternate alleles in a given record
    "G": -2,  # Equal to the number of genotypes in a given record
    "R": -3,  # Equal to the number of alleles including reference in a given record
}


_Info = collections.namedtuple(
    "Info", ["id", "num", "type", "desc", "source", "version"]
)
_Filter = collections.namedtuple("Filter", ["id", "desc"])
_Alt = collections.namedtuple("Alt", ["id", "desc"])
_Format = collections.namedtuple("Format", ["id", "num", "type", "desc"])
_SampleInfo = collections.namedtuple(
    "SampleInfo", ["samples", "gt_bases", "gt_types", "gt_phases"]
)
_Contig = collections.namedtuple("Contig", ["id", "length"])


class _vcf_metadata_parser(object):
    """Parse the metadata in the header of a VCF file."""

    def __init__(self):
        super(_vcf_metadata_parser, self).__init__()
        self.info_pattern = re.compile(
            r"""\#\#INFO=<
            ID=(?P<id>[^,]+),\s*
            Number=(?P<number>-?\d+|\.|[AGR])?,\s*
            Type=(?P<type>Integer|Float|Flag|Character|String),\s*
            Description="(?P<desc>[^"]*)"
            (?:,\s*Source="(?P<source>[^"]*)")?
            (?:,\s*Version="?(?P<version>[^"]*)"?)?
            >""",
            re.VERBOSE,
        )
        self.filter_pattern = re.compile(
            r"""\#\#FILTER=<
            ID=(?P<id>[^,]+),\s*
            Description="(?P<desc>[^"]*)"
            >""",
            re.VERBOSE,
        )
        self.alt_pattern = re.compile(
            r"""\#\#ALT=<
            ID=(?P<id>[^,]+),\s*
            Description="(?P<desc>[^"]*)"
            >""",
            re.VERBOSE,
        )
        self.format_pattern = re.compile(
            r"""\#\#FORMAT=<
            ID=(?P<id>.+),\s*
            Number=(?P<number>-?\d+|\.|[AGR]),\s*
            Type=(?P<type>.+),\s*
            Description="(?P<desc>.*)"
            >""",
            re.VERBOSE,
        )
        self.contig_pattern = re.compile(
            r"""\#\#contig=<
            ID=(?P<id>[^>,]+)
            (,.*length=(?P<length>-?\d+))?
            .*
            >""",
            re.VERBOSE,
        )
        self.meta_pattern = re.compile(r"""##(?P<key>.+?)=(?P<val>.+)""")

    def vcf_field_count(self, num_str):
        """Cast vcf header numbers to integer or None"""
        if num_str is None:
            return None
        elif num_str not in field_counts:
            # Fixed, specified number
            return int(num_str)
        else:
            return field_counts[num_str]

    def read_info(self, info_string):
        """Read a meta-information INFO line."""
        match = self.info_pattern.match(info_string)
        if not match:
            raise SyntaxError("One of the INFO lines is malformed: %s" % info_string)

        num = self.vcf_field_count(match.group("number"))

        info = _Info(
            match.group("id"),
            num,
            match.group("type"),
            match.group("desc"),
            match.group("source"),
            match.group("version"),
        )

        return (match.group("id"), info)

    def read_filter(self, filter_string):
        """Read a meta-information FILTER line."""
        match = self.filter_pattern.match(filter_string)
        if not match:
            raise SyntaxError(
                "One of the FILTER lines is malformed: %s" % filter_string
            )

        filt = _Filter(match.group("id"), match.group("desc"))

        return (match.group("id"), filt)

    def read_alt(self, alt_string):
        """Read a meta-information ALTline."""
        match = self.alt_pattern.match(alt_string)
        if not match:
            raise SyntaxError("One of the ALT lines is malformed: %s" % alt_string)

        alt = _Alt(match.group("id"), match.group("desc"))

        return (match.group("id"), alt)

    def read_format(self, format_string):
        """Read a meta-information FORMAT line."""
        match = self.format_pattern.match(format_string)
        if not match:
            raise SyntaxError(
                "One of the FORMAT lines is malformed: %s" % format_string
            )

        num = self.vcf_field_count(match.group("number"))

        form = _Format(match.group("id"), num, match.group("type"), match.group("desc"))

        return (match.group("id"), form)

    def read_contig(self, contig_string):
        """Read a meta-contigrmation INFO line."""
        match = self.contig_pattern.match(contig_string)
        if not match:
            raise SyntaxError(
                "One of the contig lines is malformed: %s" % contig_string
            )
        length = self.vcf_field_count(match.group("length"))
        contig = _Contig(match.group("id"), length)
        return (match.group("id"), contig)

    def read_meta_hash(self, meta_string):
        # assert re.match("##.+=<", meta_string)
        items = meta_string.split("=", 1)
        # Removing initial hash marks
        key = items[0].lstrip("#")
        # N.B., items can have quoted values, so cannot just split on comma
        val = OrderedDict()
        state = 0
        k = ""
        v = ""
        for c in items[1].strip("[<>]"):

            if state == 0:  # reading item key
                if c == "=":
                    state = 1  # end of key, start reading value
                else:
                    k += c  # extend key
            elif state == 1:  # reading item value
                if v == "" and c == '"':
                    v += c  # include quote mark in value
                    state = 2  # start reading quoted value
                elif c == ",":
                    val[k] = v  # store parsed item
                    state = 0  # read next key
                    k = ""
                    v = ""
                else:
                    v += c
            elif state == 2:  # reading quoted item value
                if c == '"':
                    v += c  # include quote mark in value
                    state = 1  # end quoting
                else:
                    v += c
        if k != "":
            val[k] = v
        return key, val

    def read_meta(self, meta_string):
        if re.match("##.+=<", meta_string):
            return self.read_meta_hash(meta_string)
        match = self.meta_pattern.match(meta_string)
        if not match:
            # Spec only allows key=value, but we try to be liberal and
            # interpret anything else as key=none (and all values are parsed
            # as strings).
            return meta_string.lstrip("#"), "none"
        return match.group("key"), match.group("val")
