import ete3
import sys
from .utils import load_species_name, load_species_name_whole, run_command


def obtain_duptree_file(treeFile, duptreeFile):
    print("Creating duptree file")
    outfile = open(duptreeFile, "w")
    for line in open(treeFile):
        line = line.strip()
        dades = line.split()
        t = ete3.PhyloTree(dades[-1], sp_naming_function=load_species_name)
        for leaf in t.iter_leaves():
            leaf.name = leaf.species
        # this is very important! before each tree was rooted with midpoint rooting I don't know why.
        t.resolve_polytomy()
        # #Root tree for new duptree
        # try:
        #     t.set_outgroup(t.get_midpoint_outgroup())
        # except:
        #     node = t.get_leaves()[0]
        #     t.set_outgroup(node)
        outfile.write(t.write(format=9) + "\n")
    outfile.close()


# Loads species tree and gets the spe2age dictionary
def load_species_tree(spTreeFile):
    try:
        t = ete3.PhyloTree(spTreeFile, sp_naming_function=load_species_name_whole)
    except:
        print(spTreeFile)
        exit("Unable to load species tree")
    if len(t.get_children()) != 2:
        sys.exit("Species tree needs to be rooted")
    return t


def launch_duptree(fileName, spTreeFile, duptreePath, outDir, tagName):
    # Launches duptree program to calculate the species tree and creates a species tree file
    cmd = (
        duptreePath
        + " -i "
        + fileName
        + " -o "
        + outDir
        + "/"
        + tagName
        + "results_duptree.txt 1>duptree.log 2>duptree.error"
    )
    run_command(cmd)
    cmd = "rm duptree.log"
    run_command(cmd)
    cmd = "rm duptree.error"
    run_command(cmd)
    cmd = "tail -n1 " + outDir + "/" + tagName + "results_duptree.txt > " + spTreeFile
    run_command(cmd)


# ASTRAL

# ./astral_c/ASTER/astral-pro -o all_astral/apro_phylome_$i.nwk -u 1 -t 4 all_best_trees/best_trees$i\_renamed.txt 2>&1 | tee all_astral/apro_phylome_$i.log
