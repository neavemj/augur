"""
all of this information can one day make it to the DB I would think
It lives in a seperate file simply to make flu.prepare.py less cluttered
"""

segments = ["pb2", "pb1", "pa", "ha", "np", "na", "mp", "ns"]

# regions is list of tuples (region, acronym, popsize)
# acronym = "" means ignore for frequency calcs
regions = [
    ('africa',            "", 1.216e9),
    ('europe',            "EU", 0.739e9),
    ('north_america',     "NA", 0.58e9),
    ('china',             "AS", 1.4e9),
    ('south_asia',        "AS", 1.75e9),
    ('japan_korea',       "AS", 0.2e9),
    ('south_pacific',     "OC", 0.002e9),
    ('oceania',           "OC", 0.038e9),
    ('south_america',     "", 0.422e9),
    ('southeast_asia',    "AS", 0.618e9),
    ('west_asia',         "AS", 0.245e9)
]

outliers = {
    'fluB':[
        "A/Malaysia/438/2016", "B/Togo/LNG/419/2013", "B/Brisbine/33/2008", "B/Kol/2024/2008",
        "B/Kolkata/1373/2008", "B/Kolkata/2024/2008", "B/Kolkata/372/2010", "B/Cambodia/26/2011",
        "B/Cambodia/30/2011", "B/Cambodia/62/2011", "B/Cambodia/89/2011", "B/Cambodia/V1005378/2011",
        "B/Stockholm/7/2011", "B/Bangkok/SI17/2012", "B/Bangkok/SI58/2012", "B/SouthAustralia/81/2012",
        "B/Netherlands/76/2014", "B/NewCaledonia/119/2015", "B/Thailand/CU-B11637/2015",
        "B/Brisbane/14/2016", "B/Sydney/6/2016", "B/Netherlands/883/2016",
        "B/Kisumu/7/2005", "B/Nairobi/351/2005", "B/Kolkata/2546/2009", "B/Kolkata/N-1272/2009",
        "B/Kolkata/N-2047/2009", "B/Riyadh/3/2010", "B/Riyadh/4/2010", "B/England/581/2012",
        "B/Thailand/CU-B10303/2014", "B/Catalonia/NSVH100562319/2017", "B/Norway/2155/2017"
    ]
}

reference_maps = {
####### Consider using yam outgroup instead of vic outgroup
    # "yam": {
    #     "ha": {
    #         "path": "metadata/yam_outgroup.gb",
    #         "metadata": {
    #             'strain': "B/Singapore/11/1994",
    #             'isolate_id': "CY019707",
    #             'date': "1994-05-10",
    #             'region': "southeast_asia",
    #             'country': "Singapore",
    #             "city": "Singapore",
    #             "passage": "unknown",
    #             'lab': "unknown",
    #             'age': "unknown",
    #             'gender': "M"
    #         },
    #         "include": 0,
    #         "genes": ["SigPep", "HA1", "HA2"]
    #     },
    #     "na": {
    #         "path": "metadata/yam_na_outgroup.gb",
    #         "metadata": {
    #             'strain': "B/Singapore/11/1994",
    #             'isolate_id': "CY019707",
    #             'date': "1994-05-10",
    #             'region': "southeast_asia",
    #             'country': "Singapore",
    #             "city": "Singapore",
    #             "passage": "unknown",
    #             'lab': "unknown",
    #             'age': "unknown",
    #             'gender': "M"
    #         },
    #         "include": 0,
    #         "genes": ["NA", "NB"]
    #     }
    # },
    "fluB": { # previously vic outgroup
        "ha": {
            "path": "metadata/vic_outgroup.gb",
            "metadata": {
                'strain': "B/Hong Kong/02/1993",
                'isolate_id': "CY018813",
                'date': "1993-02-15",
                'region': "china",
                'country': "Hong Kong",
                "city": "Hong Kong",
                "passage": "4",
                'lab': "unknown",
                'age': "unknown",
                'gender': "unknown"
            },
            "include": 0,
            "genes": ["SigPep", "HA1", "HA2"]
        },
        "na": {
            "path": "metadata/vic_na_outgroup.gb",
            "metadata": {
                'strain': "B/Hong Kong/02/1993",
                'isolate_id': "CY018813",
                'date': "1993-02-15",
                'region': "china",
                'country': "Hong Kong",
                "city": "Hong Kong",
                "passage": "4",
                'lab': "unknown",
                'age': "unknown",
                'gender': "unknown"
            },
            "include": 0,
            "genes": ["NB", "NA"]
        }
    }
}

## lots of the references share data
reference_maps["yam"]["na"]["metadata"] = reference_maps["yam"]["ha"]["metadata"]

# these will be included in narrower builds, for example 'A/Michigan/15/2014' will appear in a 2015-2017
# place vaccine strains here - this ensures they'll be sampled
reference_viruses = {
    'fluB':[
        'B/Shangdong/7/1997', 'B/HongKong/330/2001', 'B/Malaysia/2506/2004', 'B/Ohio/1/2005',
        'B/Brisbane/60/2008', 'B/Utah/8/2012', 'B/Montana/5/2012', 'B/Texas/2/2013', 'B/Florida/33/2014',
        'B/Indiana/25/2015', 'B/Florida/78/2015', 'B/Colorado/6/2017',
        'B/Beijing/184/1993', 'B/Sichuan/379/1999', 'B/Shanghai/361/2002', 'B/Florida/4/2006',
        'B/Wisconsin/1/2010', 'B/Massachusetts/2/2012', 'B/Phuket/3073/2013', 'B/Utah/9/2014',
        'B/Brisbane/9/2014', 'B/Guangdong-Liwan/1133/2014', 'B/California/12/2015', 'B/Arizona/10/2015',
        'B/NewHampshire/1/2016'
    ]
}

# vaccine choices copied from https://github.com/blab/nextflu/tree/de1c8323f75b0fbac9cf26451380d2758288290e/auspice/_includes feb 1 2018
vaccine_choices = {
    "fluB": {
        'B/Shangdong/7/1997': "1999-09-25",
        'B/HongKong/330/2001': "2002-09-25",
        'B/Malaysia/2506/2004': "2006-09-25",
        'B/Brisbane/60/2008': "2009-09-25",
        'B/Colorado/6/2017': "2018-02-22",
        'B/Beijing/184/1993': "1998-11-01",
        'B/Sichuan/379/1999': "2001-09-25",
        'B/Shanghai/361/2002': "2004-09-25",
        'B/Florida/4/2006': "2008-09-25",
        'B/Wisconsin/1/2010': "2012-02-25",
        'B/Massachusetts/2/2012': "2013-02-25",
        'B/Phuket/3073/2013': "2014-09-25"
    },
}

# Local Branching Index (LBI) params
LBI_params = {
    '3y': {"tau": 0.4, "time_window": 0.6},
    '12y': {"tau": 0.0005, "time_window": 0.5}
}

# Frequency Params
frequency_params = {
    '3y': {"dfreq_dn": 6},
    '12y': {"dfreq_dn": 6}
}

# Map resolution to pivot spacing
resolution_to_pivot_spacing = {
    "3y": 1. / 12.,
    "12y": 3. / 12.
}

clade_designations = {
    "fluB":{
        '1A': [('nuc',206,'G'), ('nuc',644,'C'), ('nuc',1340,'T'), ('nuc',1821,'T'),
            ('HA1',165,'K'), ('HA1',172,'P')],
        '1B': [('nuc',1034,'G'), ('nuc',1172,'G'), ('HA1',165,'K'), ('HA1',172,'P')],
        '117V': [('HA1', 75,'K'), ('HA1', 58, 'L'), ('HA1', 165, 'K'), ('HA1', 129, 'D'), ('HA1', 117, 'V')],
        'DV': [('HA1', 75,'K'), ('HA1', 58, 'L'), ('HA1', 165, 'K'), ('HA1', 117, 'V'), ('HA1', 180, 'V'), ('HA2', 152, 'K'), ('HA1', 162, 'X')],
        '2':  [('HA1', 48,'K'), ('HA1', 108, 'A'), ('HA1', 150, 'S')],
        '3':  [('HA1', 48,'R'), ('HA1', 108, 'P'), ('HA1', 150, 'I')],
        '3a': [('HA1', 37,'A'), ('HA1', 298, 'E'), ('HA1', 48,'R'), ('HA1', 105, 'P'), ('HA1', 150, 'I')],
        '172Q': [('HA1', 48,'R'), ('HA1', 108, 'P'), ('HA1', 150, 'I'), ('HA1', 116, 'K'), ('HA1', 172, 'Q'), ('HA1', 298, 'E'), ('HA1', 312, 'K')]
    }
}
