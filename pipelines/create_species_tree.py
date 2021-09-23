#!/usr/bin/env python

import analyse_phylome as ap
import argparse
import sys

# Parse data
parser = argparse.ArgumentParser(
    description="Script to Create a species tree from PhylomeDB data."
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
    default="out_sp",
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
    choices=["cat", "sprax", "dup", "apro"],
    action="store",
    default="dup",
    help="Method to create species tree",
)
parser.add_argument(
    "--threads", type=int, action="store", default=1, help="Number of threads to use.",
)
parser.add_argument(
    "--starting_sprax",
    action="store",
    default="MiniNJ",
    help="Starting species tree for generax. The default is MININJ method. other options are filepath with another species tree or random",
)
parser.add_argument(
    "--cat_alignment",
    action="store",
    help="Concatenated alignment. Use it with cat method",
)
# parser.add_argument(
#     "--arguments",
#     dest="arguments",
#     default="",
#     action="store",
#     help="Extra arguments for the programs",
# )


args = parser.parse_args()

ap.create_folder(args.outdir)

trees = args.besttrees
phy_id = str(args.phyID)

if args.method == "dup":
    duptreeFile = args.outdir + "/duptree_in.txt"
    tagName = ""
    out_file = args.outdir + "/phylome_" + phy_id + "_sp_duptree.nwk"
    ap.obtain_duptree_file(trees, duptreeFile)
    ap.launch_duptree(duptreeFile, out_file, args.exe, args.outdir, tagName)

if args.method == "apro":
    out_file = args.outdir + "/phylome_" + phy_id + "_astral_trees.txt"
    out_tree = args.outdir + "/phylome_" + phy_id + "_apro.nwk"
    out_log = args.outdir + "/phylome_" + phy_id + "_apro.log"
    ap.get_astral_pro_data(trees, out_file)
    cmd = (
        args.exe
        + " -o "
        + out_tree
        + " -u 1 -t "
        + str(args.threads)
        + " "
        + out_file
        + " 2>&1 | tee "
        + out_log
    )
    ap.run_command(cmd)

# maybe option to clean up all output files
if args.method == "sprax":
    out_sprax_dir = args.outdir + "/sprax_phylome_" + phy_id
    ap.get_generax_data(trees, out_sprax_dir)
    # If you want oother options directly modify this part!
    # It's bad but otherwise too many params.
    cmd = (
        "mpiexec -np "
        + str(args.threads)
        + " "
        + args.exe
        + " --families "
        + out_sprax_dir
        + "/family.txt"
        + " --strategy SKIP --si-strategy HYBRID --species-tree "
        + args.starting_sprax
        + " --rec-model UndatedDTL --prefix "
        + out_sprax_dir
    )
    ap.run_command(cmd)


if args.method == "cat":
    if args.cat_alignment is None:
        sys.exit("To run cat method please provide an alignment")
    cmd = (
        "python "
        + args.exe
        + " -s t -i "
        + args.cat_alignment
        + " -o "
        + args.outdir
        + " --treeThreads "
        + str(args.threads)
    )
    ap.run_command(cmd)
