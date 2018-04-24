from __future__ import print_function
import os, sys
# we assume (and assert) that this script is running from the virus directory, i.e. inside H7N9 or zika
sys.path.append(os.path.join(os.path.dirname(__file__), '../..'))
import base.prepare
from base.prepare import prepare
from datetime import datetime
from base.utils import fix_names
import argparse

dropped_strains = [
]
forced_strains = [
]

def collect_args():
    """Returns an Pangea-specific argument parser.
    """
    parser = base.prepare.collect_args()
    segments = ["env", "gag", "pol"]
    parser.add_argument('-s', '--segments', choices=segments, default=segments, nargs='+', type = str,  help = "segments to prepare")
    parser.set_defaults(
        viruses_per_month=0
    )
    return parser

def make_config(params):
    return {
        "dir": "pangea",
        "file_prefix": "pangea",
        "title": "Phylogenetics And Networks for Generalised HIV Epidemics in Africa",
        "maintainer": ["PANGEA Consortium", "http://www.pangea-hiv.org/"],
        "input_paths": [
            "../../../pangea-ii/sequences/pangea_ref_env.fasta",
            "../../../pangea-ii/sequences/pangea_ref_gag.fasta",
            "../../../pangea-ii/sequences/pangea_ref_pol.fasta"
        ],
        # >PG16-BW003077|env|BW-PRI-INF|Botswana|2007-06-11
        "header_fields": {0:'strain', 1: 'accession', 2: 'subtype', 3: 'segment', 4:'authors', 5:'country', 6:'date'},
        "filters": (
            ("Dropped Strains", lambda s: s.id not in [fix_names(x) for x in dropped_strains]),
            ("Sequence Length", {
                "env": lambda s: len(s.seq)>=1000,
                "gag": lambda s: len(s.seq)>=1000,
                "pol": lambda s: len(s.seq)>=1000
            })
        ),
        "subsample": {
            "category": lambda x:(x.attributes['country'], x.attributes['date'].year),
            "threshold": params.viruses_per_month,
            "priority": lambda x:x.id in forced_strains
        },
        "colors": ["country"],
        "color_defs": ["./colors.tsv"],
        "lat_longs": ["country"],
        "auspice_filters": ["country", "authors", "subtype"],
        "references": {
            # references are pinneo strain. Same as Kristian's Cell paper.
            # Pinneo paper: http://jvi.asm.org/content/74/15/6992.long
            # Cell paper: http://www.cell.com/cell/pdfExtended/S0092-8674(15)00897-1
            "env": {
                "path": "metadata/pangea.gb",
                "metadata": {
                    'strain': "R880_MPL_C11", "accession": "KP223844", "date": "2007-01-12",
                    'country': "rwanda", 'segment': 'env'
                },
                "include": 0,
                "genes": ['env']
            },
            "gag": {
                "path": "metadata/pangea.gb",
                "metadata": {
                    'strain': "R880_MPL_C11", "accession": "KP223844", "date": "2007-01-12",
                    'country': "rwanda", 'segment': 'gag'
                },
                "include": 0,
                "genes": ['gag']
            },
            "pol": {
                "path": "metadata/pangea.gb",
                "metadata": {
                    'strain': "R880_MPL_C11", "accession": "KP223844", "date": "2007-01-12",
                    'country': "rwanda", 'segment': 'pol'
                },
                "include": 0,
                "genes": ['pol']
            },
        }
    }


if __name__=="__main__":
    parser = collect_args()
    params = parser.parse_args()
    config = make_config(params)
    config["segments"] = params.segments

    if params.viruses_per_month == 0:
        config["subsample"] = False

    runner = prepare(config)
    runner.load_references()
    runner.applyFilters()
    # runner.ensure_all_segments()
    runner.subsample()
    runner.colors()
    runner.latlongs()
    runner.write_to_json(segment_addendum=True)
