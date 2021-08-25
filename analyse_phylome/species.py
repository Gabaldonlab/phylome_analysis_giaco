import ete3
import sys
from .utils import load_species_name, load_species_name_whole, run_command


def obtain_duptree_file(treeFile, duptreeFile, midpoint=False):
    """Prepare files for duptree.

    Parameters
    ----------
    treeFile : str
        best trees file from PhylomeDB
    duptreeFile : str
        File where the tress will be written
    midpoint: bool
        Midpoint root gene trees?

    Returns
    -------
    type
        File that will be the input of duptree

    """
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
        #Root tree for new duptree
        # if root:
        if midpoint:
            t.set_outgroup(t.get_midpoint_outgroup())
            #if spe2age:
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
    """Run duptree.

    Parameters
    ----------
    fileName : type
        file produced by obtain_duptree_file()
    spTreeFile : type
        Output species tree file
    duptreePath : str
        Path of duptree exe
    outDir : str
        output directory
    tagName : str
        tag

    Returns
    -------
    type
        Duptree species tree

    """
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
