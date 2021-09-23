#!/usr/bin/env python

import analyse_phylome as ap
import argparse
from ete3 import Tree
import json
import sys

# Parse data
parser = argparse.ArgumentParser(
    description="Script to concatenate PhylomeDB alignments."
)
parser.add_argument(
    "-p",
    "--phylomeID",
    dest="phyID",
    action="store",
    default="None",
    help="Current number with which the phylome was registered (PhylomeID). Consult the phylomeDB website if you're not sure about which number was assigned to your phylome.",
    required=True,
)
parser.add_argument(
    "-o",
    "--outdir",
    action="store",
    default="out_cat",
    help="Directory where the data will be stored",
)
parser.add_argument(
    "-r",
    "--readalPath",
    dest="readal",
    action="store",
    default="None",
    required=True,
    help="Directory of readal executable",
)
parser.add_argument(
    "-s",
    "--species_tree",
    dest="species_tree",
    action="store",
    help="Provide species tree to build spe2age dictionry only use it with seed anc method 2. Otherwise provide json root_phy file.",
)
parser.add_argument(
    "--seed",
    dest="seed",
    action="store",
    help="Seed species of phylome. This is only needed if a species tree is given and concat method 2 is chosen.",
)
parser.add_argument(
    "--spe2age_file",
    dest="spe2age_file",
    action="store",
    help="If no species tree is present use root_phy.json to build spe2age. Only use it if concat method 2.",
)
parser.add_argument(
    "-a",
    "--alndir",
    dest="alndir",
    action="store",
    default="None",
    help="Directory containing aln to concatenate",
    required=True,
)
parser.add_argument(
    "-t",
    "--besttrees",
    action="store",
    help="File containing a list of trees. Needed if taking data from a file",
    required=True,
)
parser.add_argument(
    "-m",
    "--method",
    dest="method",
    choices=["std", "nodup", "ortho"],
    action="store",
    default="std",
    help="Concatenate alignment with more than 'prop' species",
)
parser.add_argument(
    "--prop",
    default=0.9,
    type=float,
    action="store",
    help="Proportion of species present in alignment in order to concatenate. Only use it with method 1.",
)
parser.add_argument(
    "--at_least",
    default=500,
    type=int,
    action="store",
    help="Max Number of genes to concatenate before finishing",
)
parser.add_argument(
    "--max_length",
    default=500000,
    type=int,
    action="store",
    help="Max length of concatenated aln",
)
parser.add_argument(
    "--partition", action="store_true", help="Get partition file.",
)
# parser.add_argument("--image_off",dest="image_off",action="store_true",help="Will switch off image printing.")

args = parser.parse_args()

ap.create_folder(args.outdir)

phy_id = str(args.phyID)
pathsFile = args.outdir + "/paths.txt"
trees121File = args.outdir + "/trees121.txt"
tag = "alg_aa"
path = args.alndir
treeFile = args.besttrees

if args.method == "nodup":
    if args.species_tree is not None and args.seed is not None:
        sp = Tree(args.species_tree)
        spe2age = ap.build_sp2age(sp, str(args.seed))
    elif args.spe2age_file is not None:
        with open(args.spe2age_file) as f:
            all_dicts = json.load(f)
        spe2age = all_dicts[phy_id]
    else:
        sys.exit(
            "If concat method nodup give either species tree + seed or json root file"
        )
else:
    spe2age = ap.get_all_species(treeFile)

ap.create_pathfile(pathsFile, path, tag)

if args.method == "std":
    ap.obtain_121_trees(treeFile, trees121File)
    ap.build_concatenated_alg(
        trees121File,
        spe2age,
        pathsFile,
        args.outdir,
        args.readal,
        prop=args.prop,
        at_least=args.at_least,
        max_length=args.max_length,
        partition=args.partition,
        treeFile=treeFile,
    )

if args.method == "nodup":
    ap.build_extra_concatenated_alg2(
        treeFile,
        pathsFile,
        args.outdir,
        args.readal,
        spe2age,
        prop=args.prop,
        at_least=args.at_least,
        max_length=args.max_length,
        partition=args.partition,
    )

if args.method == "ortho":
    ap.build_extra_concatenated_alg3(
        treeFile,
        pathsFile,
        args.outdir,
        args.readal,
        spe2age,
        prop=args.prop,
        at_least=args.at_least,
        max_length=args.max_length,
        partition=args.partition,
    )
