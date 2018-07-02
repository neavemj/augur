# currently the RHDV co-ordinates are grouped by state
# want to add their actual co-ordinates
# Matthew J. Neave 02.07.2018

import sys

fasta_fl = sys.argv[1]
coords_fl = sys.argv[2]
fasta_output = open(sys.argv[3], "w")
coords_output = open(sys.argv[4], "w")

# make a dict of the sample name and its coords

coords_dict = {}

with open(coords_fl) as fl:
    next(fl)
    for line in fl:
        line = line.strip()
        cols = line.split("\t")
        isolate = cols[0].lower()
        if len(cols) < 4:
            gps = None
        else:
            gps = cols[3]
        if isolate in coords_dict:
            print("isolate name {} is not unique!".format(isolate))
            sys.exit(1)
        else:
            coords_dict[isolate] = gps

# read through fasta file, get isolate name and print out nicely

with open(fasta_fl) as fl:
    coords_output.write("location\tcountry_code\tlatitude\tlongitude\n")
    for line in fl:
        if line.startswith(">"):
            line = line.strip()
            cols = line.split("|")
            isolate = cols[0].split("/")[2].lower()
            if isolate not in coords_dict:
                print("{} does not contain GPS co-ordinates!".format(isolate))
            else:
                gps_edit = [gps.strip().strip('"') for gps in coords_dict[isolate].split(",")]
                coords_output.write(isolate + "\tAUS\t" + "\t".join(gps_edit) + "\n")
            cols.insert(1, isolate)
            fasta_output.write("|".join(cols) + "\n")
        else:
            fasta_output.write(line)

