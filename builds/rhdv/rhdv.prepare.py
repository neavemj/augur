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

dropped_strains = ["AUS/NSW/YAR-1/2015", "AUS/WA/GER-2/2016", "AUS/WA/GER-1/2016", "AUS/WA/GER-5/2016",
                   "AUS/WA/GER-4/2016", "AUS/WA/GER-3/2016", "AUS/ACT/BLMT-3/2015", "AUS/ACT/AIN-3/2015",
                   "AUS/TAS/GEE-2/2016", "AUS/VIC/CWS-1/2016", "AUS/TAS/BRK-1/2016", "AUS/TAS/GRA-1/2016"
] # these get dropped because they are too different. If I don't remove them here the filter goes the other way and
# basically removes all the sequences - only when blmt-1 is included though.

config = {
    "dir": "rhdv",
    "file_prefix": "rhdv",
    "title": "Monitoring the initial spread and genomic variability of rabbit hemorrhagic disease virus 2 (GI.2) in "
             "the Australian landscape: 2015-2016",
    "maintainer": ["Robyn Hall and Matthew Neave", "mailto:robyn.hall@csiro.au"],
    "input_paths": ["./data/mahar_RHDV.edit.fasta"],
    #>AUS/ACT/BLMT-3/2015|blmt-3|RHDV|MF421563.1|RHDV1_G2|2015-06-18|Australia|ACT|Mahar et al|Monitoring the init
    "header_fields": {0:'strain', 1:'isolate', 4:'rhdv_strain', 5:'date', 6:'country', 7:'state', 8: 'authors',
                      10:'journal', 9:'title'},
    "filters": (
        ("Dropped Strains", lambda s: s.id not in [fix_names(x) for x in dropped_strains]),
        # ("Restrict Date Range", lambda s: s.attributes['date'] >= datetime(2012,01,1).date()),
        # ("Restrict Date Range", lambda s: s.attributes['date'] <= datetime(2018,01,1).date()),
        #("Sequence Length", lambda s: len(s.seq)>=10000),
    ),
    "subsample": {
        "category": lambda x:(x.attributes['date'].year, x.attributes['date'].month),
    },
    "colors": ["state", "authors", "rhdv_strain"],
    # "color_defs": ["./colors.tsv"],
    "lat_longs": ["isolate"],
    "lat_long_defs": "./RHDV2_coords.edit.txt",
    "auspice_filters": ["state", "authors", "rhdv_strain"],
    "reference": {
        "path": "sequence2.gb",
        "metadata": {
            'strain': "RHDV1_G2", "date": "2015-05-13", "country": "AUSTRALIA", "state": "ACT"
        },
        "include": 0,
        "genes": ['polyprotein', 'VP10']
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
