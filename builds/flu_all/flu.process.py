from __future__ import print_function
import os, sys
# we assume (and assert) that this script is running from the virus directory, i.e. inside H7N9 or zika
sys.path.append(os.path.join(os.path.dirname(__file__), '../..'))
import base.process
from base.process import process
import argparse


def collect_args():
    """Returns an avian flu-specific argument parser.
    """
    parser = base.process.collect_args()
    # parser.set_defaults(
    #     json="prepared/avian_h7n9_ha.json"
    # )
    return parser


config = {
    "dir": "flu_all",
    "geo_inference": ["segment", "type", "lineage"], # what traits to perform this on
    "auspice": { ## settings for auspice JSON export
        "panels": ['tree', 'map', 'entropy'],
        "defaults": {
            "colorBy": "type",
            "geoResolution": "country"
        },
        "extra_attr": [],
        "date_range": {},
        "color_options": {
            "country":{"key":"country", "legendTitle":"Country", "menuItem":"country", "type":"discrete", "color_map": []},
            "region":{"key":"region", "legendTitle":"Region", "menuItem":"region", "type":"discrete"},
            "lineage":{"key":"lineage", "legendTitle":"lineage", "menuItem":"lineage", "type":"discrete", "color_map": []},
            "type":{"key":"type", "legendTitle":"type", "menuItem":"type", "type":"discrete", "color_map": []},
            "ha_type":{"key":"ha_type", "legendTitle":"ha_type", "menuItem":"ha_type", "type":"discrete"},
            "na_type":{"key":"na_type", "legendTitle":"na_type", "menuItem":"na_type", "type":"discrete"},
            "host":{"key":"host", "legendTitle":"Host", "menuItem":"host", "type":"discrete", "color_map": []},
            "num_date":{"key":"num_date", "legendTitle":"Sampling date", "menuItem":"date", "type":"continuous"},
            "gt":{"key":"genotype", "legendTitle":"Genotype", "menuItem":"genotype", "type":"discrete"}
        },
        "controls": {'geographic location':['country']}
    },
    "newick_tree_options": {}
}

if __name__=="__main__":
    parser = collect_args()
    params = parser.parse_args()

    config["in"] = params.json
    config["clean"] = params.clean
    config["newick_tree_options"]["raxml"] = not params.no_raxml

    runner = process(config)
    runner.align()
    runner.build_tree()
    runner.timetree_setup_filter_run()
    #runner.run_geo_inference()
    runner.save_as_nexus()
    runner.auspice_export()
