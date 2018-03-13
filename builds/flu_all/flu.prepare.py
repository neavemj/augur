from __future__ import print_function
import os, sys
# we assume (and assert) that this script is running from the virus directory, i.e. inside H7N9 or zika
sys.path.append(os.path.join(os.path.dirname(__file__), '../..'))
sys.path.append('reference_segments')
from reference_info import references
import base.prepare
from base.prepare import prepare
from datetime import datetime
from base.utils import fix_names
import argparse


def collect_args():
    parser = base.prepare.collect_args()
    segments = ['pb2', 'pb1', 'pa', 'ha', 'np', 'na', 'mp', 'ns']
    parser.add_argument('-s', '--segments', choices=segments, default=segments, nargs='+', type = str,  help = "segments to prepare")
    parser.set_defaults(
        viruses_per_month=0
    )
    return parser

dropped_strains = [
]

config = {
    "dir": "flu_all", # the current directory. You must be inside this to run the script.
    "title": "testing",
    "maintainer": ["James Hadfield", "http://bedford.io/team/james-hadfield/"],
    "input_paths": [
        "../../../flora/downloaded_data/flua.PB2.fasta",
        "../../../flora/downloaded_data/flua.PB1.fasta",
        "../../../flora/downloaded_data/flua.PA.fasta",
        "../../../flora/downloaded_data/flua.HA.fasta",
        "../../../flora/downloaded_data/flua.NP.fasta",
        "../../../flora/downloaded_data/flua.NA.fasta",
        "../../../flora/downloaded_data/flua.MP.fasta",
        "../../../flora/downloaded_data/flua.NS.fasta"
    ],
    "header_fields": {
        0: 'strain',
        1: 'date',
        2: 'country',
        3: 'region',
        4: 'segment',
        5: 'type',
        6: 'ha_type',
        7: 'na_type',
        8: 'host',
        9: 'accession'
    },
    "filters": (
        ("Dropped Strains", lambda s: s.id not in [fix_names(x) for x in dropped_strains]),
        ("Missing Month & Year", lambda s: "XX-XX" not in s.attributes['raw_date']),
        ("Exclude bad host", lambda s: s.attributes["host"] not in ["laboratoryderived", "other", "unknown", 'bat']),
        ("Missing Country", lambda s: s.attributes["country"] not in ["", "?"]),
        ("Missing Segment", lambda s: s.attributes["segment"] not in ["", "?", None]),
        ("Sequence Length", {
            "pb2": lambda s: len(s.seq)>=2200,
            "pb1": lambda s: len(s.seq)>=2200,
            "pa": lambda s: len(s.seq)>=2100,
            "ha": lambda s: len(s.seq)>=1500,
            "np": lambda s: len(s.seq)>=1400,
            "na": lambda s: len(s.seq)>=1200,
            "mp": lambda s: len(s.seq)>=900,
            "ns": lambda s: len(s.seq)>=800
        })
    ),
    "subsample": {
        "category": lambda x: (x.attributes['ha_type'], x.attributes['host'], x.attributes['region'], x.attributes['date'].year),
        "priority": None,
    },
    # see the docs for what's going on with colours (sic) & lat/longs
    "colors": ["type", "ha_type", "na_type", "host", "lineage", "country", "region"],
    "color_defs": ["./colors.tsv"],
    "auspice_filters": ["country", "host", "region", "ha_type", "na_type"],
    "lat_longs": ["country", "region"],
    "references": references,
}


if __name__=="__main__":
    parser = collect_args()
    params = parser.parse_args()

    if params.viruses_per_month == 0:
        config["subsample"] = False
    else:
        config["subsample"]["threshold"] = params.viruses_per_month

    if params.sequences is not None:
        config["input_paths"] = params.sequences

    # if params.file_prefix is not None:
    #     config["file_prefix"] = params.file_prefix
    #     segment_addendum = False
    # else:

    config["segments"] = params.segments
    config["file_prefix"] = "flu_all"

    runner = prepare(config)
    runner.load_references()
    runner.applyFilters()
    runner.ensure_all_segments()
    runner.subsample(segment_idx = params.segments.index('ha'))
    runner.colors()
    runner.latlongs()
    runner.write_to_json(segment_addendum=True)
