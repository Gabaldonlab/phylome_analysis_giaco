#!/usr/bin/env python

import analyse_phylome as ap
import argparse
from ete3 import Tree
import os.path

# Parse data
parser = argparse.ArgumentParser(
    description="Script to run discrdance analysis on a phylome."
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
    default="out_disc",
    help="Directory where the data will be stored",
)
parser.add_argument(
    "-e",
    "--exe",
    dest="exe",
    action="store",
    help="Directory of chosen executable",
    required=True,
)
parser.add_argument(
    "-t",
    "--besttrees",
    action="store",
    help="File containing a list of trees. Needed if taking data from a file",
)
parser.add_argument(
    "-m",
    "--method",
    dest="method",
    choices=["ecce", "grax", "abac_phylo", "abac_taxa"],
    action="store",
    required=True,
    help="Method to use",
)
parser.add_argument(
    "--threads", type=int, action="store", default=1, help="Number of threads to use.",
)
parser.add_argument(
    "-s",
    "--species_tree",
    dest="species_tree",
    action="store",
    help="Species tree to be reconcilde",
)
parser.add_argument(
    "-a",
    "--alndir",
    dest="alndir",
    action="store",
    help="Directory containing alns. Only use it with GeneRax",
)
parser.add_argument(
    "--phy_info", action="store", help="Phylome info file. Only use it with abac",
)
parser.add_argument(
    "--arguments",
    dest="arguments",
    default="",
    action="store",
    help="Extra arguments for the programs",
)
args = parser.parse_args()

ap.create_folder(args.outdir)

trees = args.besttrees
phy_id = str(args.phyID)

# print_all = "print.reconciliations=1 print.info=1  orthology.output=1 print.newick=1 print.graph=1"

if args.method == "ecce":
    ecce_trees = args.outdir + "/" + "ecce_best_trees" + phy_id + ".txt"
    ap.get_ecce_data(trees, args.species_tree, args.outdir)
    log_file = args.outdir + "/phylome_" + phy_id + "_ecce.out"
    cmd = (
        args.exe
        + " species.file="
        + args.species_tree
        + " gene.file="
        + ecce_trees
        + " dated=0 verbose=1 print.reconciliations=1 recPhyloXML.reconciliation=true output.dir="
        + args.outdir
        + "/ecce_results"
        + args.arguments
        + " > "
        + log_file
    )
    ap.run_command(cmd)

# maybe option to clean up all output files
if args.method == "grax":
    # ap.scan_for_Us(args.alndir, replace=True)
    out_grax_dir = args.outdir + "/grax_phylome_" + phy_id
    ap.get_generax_data(trees, out_grax_dir, args.alndir)
    # If you want oother options directly modify this part!
    # It's bad but otherwise too many params.
    cmd = (
        "mpiexec -np "
        + str(args.threads)
        + " "
        + args.exe
        + " -s "
        + args.species_tree
        + " --families "
        + out_grax_dir
        + "/family.txt"
        + " --strategy EVAL --rec-model UndatedDTL --per-family-rates --prefix "
        + out_grax_dir
    )
    ap.run_command(cmd)


if "abac" in args.method:
    sptree_tree = Tree(args.species_tree)
    tax_dict_i = ap.get_tax_dict_info(args.phy_info)
    whole_tax_dict = ap.get_taxonomic_df(tax_dict_i)
    if args.method == "abac_taxa":
        out_tax = args.outdir + "/phylome_" + phy_id + "_taxon_euka.csv"
        ap.get_abaccus_taxo(tax_dict_i, out_tax)
        ap.run_abaccus(args.exe, trees, args.outdir, taxonomy=out_tax)
    if args.method == "abac_phylo":
        ap.get_proka_file(whole_tax_dict, out_dir=args.outdir)
        proka_file = args.outdir + "/proka_mnemo.txt"
        if os.path.isfile(proka_file):
            proka = proka_file
        else:
            proka = None
            ap.run_abaccus(
                args.exe, trees, args.outdir, proka=proka, taxonomy=args.species_tree
            )
    abac_out = args.outdir + "/abaccus_output.txt"
    rename_abac = args.outdir + "/phylome_" + phy_id + "_" + args.method + ".txt"
    os.rename(abac_out, rename_abac)
