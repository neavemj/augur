from __future__ import print_function
import os, sys
# we assume (and assert) that this script is running from the virus directory, i.e. inside rhdv
sys.path.append(os.path.join(os.path.dirname(__file__), '../..'))
import base.prepare
from base.prepare import prepare
from datetime import datetime
from base.utils import fix_names
import argparse

def collect_args():
    """Returns a RHDV-specific argument parser."""
    parser = base.prepare.collect_args()
    parser.set_defaults(
        viruses_per_month=1,
    )
    return parser

dropped_strains = [
]

config = {
    "dir": "rhdv",
    "file_prefix": "rhdv",
    "title": "Monitoring the initial spread and genomic variability of rabbit hemorrhagic disease virus 2 (GI.2) in "
             "the Australian landscape: 2015-2016",
    "maintainer": ["Robyn Hall and Matthew Neave", ""],
    "input_paths": ["./data/mahar_RHDV.fasta"],
    # >AUS/ACT/BLMT-3/2015|RHDV|MF421563.1|GI.1c_(RHDV1_G2)|2015-06-18|Australia|ACT|Mahar et al|Monitoring the ini...
    "header_fields": {0:'strain', 3:'rhdv_strain', 4:'date', 5:'country', 6:'state', 7: 'authors', 9:'journal',
                      8:'title'},
    "filters": (
        ("Dropped Strains", lambda s: s.id not in [fix_names(x) for x in dropped_strains]),
        # ("Restrict Date Range", lambda s: s.attributes['date'] >= datetime(2012,01,1).date()),
        # ("Restrict Date Range", lambda s: s.attributes['date'] <= datetime(2018,01,1).date()),
        #("Sequence Length", lambda s: len(s.seq)>=10000),
    ),
    "subsample": {
        "category": lambda x:(x.attributes['date'].year, x.attributes['date'].month),
    },
    "colors": ["country", "state", "authors", "strain"],
    # "color_defs": ["./colors.tsv"],
    "lat_longs": ["country"],
    "lat_longs": ["state"],
    "lat_long_defs": "./aus_state_lat_longs.tsv",
    "auspice_filters": ["country", "state", "authors", "strain"],
    "reference": {
        "path": "sequence.gb",
        "metadata": {
            'strain': "GI.1c_(RHDV1_G2)", "date": "2015-06-18", "country": "AUSTRALIA", "state": "NSW"
        },
        "include": 0,
        "genes": ['ORF1', 'ORF2']
    }
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

    if params.file_prefix is not None:
        config["file_prefix"] = params.file_prefix

    runner = prepare(config)
    runner.load_references()
    runner.applyFilters()
    runner.ensure_all_segments()
    #runner.subsample()
    runner.colors()
    runner.latlongs()
    runner.write_to_json()
