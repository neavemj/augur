from __future__ import print_function
import os, sys
# we assume (and assert) that this script is running from the virus directory, i.e. inside H7N9 or zika
sys.path.append(os.path.join(os.path.dirname(__file__), '../..'))
import base.process
from base.process import process


def collect_args():
    """Returns a Zika-specific argument parser.
    """
    parser = base.process.collect_args()
    parser.set_defaults(
        json="prepared/rhdv.json"
    )
    return parser


config = {
    "dir": "rhdv",
    "in": "prepared/rhdv.json",
    "newick_tree_options": {"nthreads": 4},
    "clock_filter": {
        "n_iqd": 4,
    },
    "geo_inference": False,
    # "geo_inference": ['country', 'state'], # what traits to perform this on
    "geo_inference_options": {
        "root_state": {
            "state": "nsw",
            "country": "australia",
        },
    },
    "auspice": { ## settings for auspice JSON export
        "color_options": {
            "authors": {"key":"authors", "legendTitle":"Authors", "menuItem":"authors", "type":"discrete"},
            "state": {"key":"state", "legendTitle":"State", "menuItem":"state", "type":"discrete"},
            "country": {"key":"country", "legendTitle":"Country", "menuItem":"country", "type":"discrete"},
            "rhdv_strain": {"key":"rhdv_strain", "legendTitle":"Strain", "menuItem":"strain", "type":"rhdv_strain"}
        },
        "controls": {},
        "defaults": {
            'mapTriplicate': False,
            'colorBy': "state"
        }
    },
    "timetree_options": {
        "Tc": 'opt',
        "confidence":True
    }
}

if __name__=="__main__":
    parser = collect_args()
    params = parser.parse_args()

    if params.clean:
        config["clean"] = True

    if params.json:
        config["in"] = params.json

    config["newick_tree_options"]["method"] = params.tree_method

    runner = process(config)
    runner.align(fill_gaps=True)
    runner.build_tree()
    runner.timetree_setup_filter_run()
    runner.run_geo_inference()
    runner.save_as_nexus()
    runner.auspice_export()
