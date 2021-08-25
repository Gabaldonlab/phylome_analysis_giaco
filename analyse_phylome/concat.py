import glob
import ete3
import os
from .utils import (
    run_command,
    load_sequences,
    load_species_name,
    load_tree,
    root_tree,
    create_folder,
    loadPaths,
    print_sequence,
)


def obtain_121_trees(treeFile, trees121File):
    """Get 1-to-1 trees from a phylome.

    Parameters
    ----------
    treeFile : str
        Best_trees file of a phylome from PhylomeDB.
    trees121File : str
        Output file where the 121 trees will be written.

    Returns
    -------
    type
        121 trees file.

    """
    outfile121 = open(trees121File, "w")
    num = 0
    for line in open(treeFile):
        line = line.strip()
        dades = line.split()
        t = ete3.PhyloTree(dades[-1], sp_naming_function=load_species_name)
        if len(t.get_leaves()) == len(t.get_species()):
            outfile121.write(
                dades[0] + "\t" + str(len(t.get_species())) + "\t" + dades[-1] + "\n"
            )
            num += 1
    outfile121.close()
    print("there are ", num, " 121 trees")


def concatenate_alignments_from_path(path, outFileName, spe2age, readalPath):
    """Concatenate alignment given the path of a directory where alns are stored.

    Parameters
    ----------
    path : str
        Dir of the alns.
    outFileName : str
        Output file name.
    spe2age : dict
        species2age dictionary obtained with build_spe2age()
    readalPath : str
        path of readal exe

    Returns
    -------
    type
        concatenated fasta.
    """
    # species_codes = []
    concatenatedAlg = {}
    for fileName in glob.glob(path + "/*"):
        # Ensuring the sequences are in fasta format
        cmd = readalPath + " -in " + fileName + " -out t -fasta"
        run_command(cmd, True)
        seqs = load_sequences("t", "")
        added = []
        # Add sequences to concatenated alignment
        for code in seqs:
            spe = code.split("_")[1].split("|")[0]
            if spe not in concatenatedAlg:
                concatenatedAlg[spe] = ""
            concatenatedAlg[spe] += seqs[code]
            added.append(spe)
        # Adding gaps to species that don't have the gene
        for spe in spe2age:
            if spe not in added:
                if spe not in concatenatedAlg:
                    concatenatedAlg[spe] = ""
                s = "-" * len(seqs[code])
                concatenatedAlg[spe] += s
    # Print concatenated alignment
    outfile = open("t", "w")
    for code in concatenatedAlg:
        # This worked in py2
        # print >> outfile, ">" + code + "\n" + concatenatedAlg[code]
        outfile.write(">" + code + "\n" + concatenatedAlg[code] + "\n")
    outfile.close()
    cmd = readalPath + " -in t -out " + outFileName + " -fasta"
    run_command(cmd, True)
    cmd = (
        readalPath
        + " -in "
        + outFileName
        + " -out "
        + outFileName.replace("fasta", "phy")
        + " -phylip"
    )
    run_command(cmd, True)


def collapse_lineage_specific_expansions(t, seed):
    events = t.get_descendant_evol_events()
    taken = set([])
    for ev in events:
        if ev.etype == "D":
            if len(ev.node.get_species()) == 1:
                leaves = set(ev.node.get_leaf_names())
                common = leaves.intersection(taken)
                if len(common) == 0:
                    for leaf in leaves:
                        taken.add(leaf)
                    if seed in leaves:
                        chosen = seed
                    else:
                        chosen = list(leaves)[0]
                    ev.node.add_feature("collapse", chosen)
    for node in t.traverse():
        if "collapse" in node.features:
            if node.is_root():
                t = None
            else:
                parent = node.up
                newName = node.collapse
                node.detach()
                parent.add_child(name=newName)
    return t


def collapse_all_lineage_specific_expansions(t, code):
    t = load_tree(t.write())  # , spe2age, midpoint)
    t = collapse_lineage_specific_expansions(t, code)
    if t:
        t = load_tree(t.write())  # , spe2age, midpoint)
        try:
            t.set_outgroup(t.get_leaves_by_name(code)[0])
        except:
            pass
        t = collapse_lineage_specific_expansions(t, code)
    return t


def create_parsed_alg_file(code, listLeaves, pathsFile, workingDir, readalPath):
    # Takes an alignment and deletes sequences from it
    pathsFile = loadPaths(pathsFile, "alg_aa")
    cmd = readalPath + " -in " + pathsFile[code] + " -out t -fasta"
    run_command(cmd, False)
    seqs = load_sequences("t", "")
    outfile = open(workingDir + "/" + code + ".alg.clean", "w")
    for code in seqs:
        if code in listLeaves:
            print_sequence(code, seqs[code], outfile)
    outfile.close()


def build_extra_concatenated_alg2(
    treeFile,
    pathsFile,
    outDir,
    readalPath,
    spe2age,#=None,
    min=5,
    # midpoint=False,
    at_least=100,
    max_length=50000
):
    """Concatenate alignments collapsing clade specific duplications.

    Parameters
    ----------
    treeFile : str
        best_trees file from phylomeDB
    pathsFile : str
        Pathfile created with create_pathfile()
    outDir : str
        output directory
    readalPath : str
        path of readal exe
    spe2age : dict
        Species 2 age dictionary or species list
    min : num
        Consider only trees with at least min species
    at_least : num
        Stop the concatenation when there are at least `at_least` gene families concatenated.
    max_length : num
        Stop the concatenation when the aln is `max_length` positions long.

    Returns
    -------
    type
        A outDir/concatenated_extra2 dir where the alignments, the concatenated fasta and phy and some stats are stored.

    """
    # Create output folders
    concatDir = outDir + "/concatenated_extra2/"
    create_folder(concatDir)
    # concatDir = concatDir + "/collapsed_leaves/"
    # create_folder(concatDir)
    cAlgsDir = concatDir + "/algs/"
    create_folder(cAlgsDir)
    # Collapse species specific duplications and obtain list of one-to-one trees
    counter = {}
    for line in open(treeFile):
        line = line.strip()
        dades = line.split("\t")
        code = dades[0]
        t = load_tree(dades[-1])  # , spe2age, midpoint)
        t = root_tree(t, spe2age)
        if code in t.get_leaf_names():
            t = collapse_all_lineage_specific_expansions(t, code)
            if t:
                t = load_tree(t.write())  # , spe2age, midpoint)
                num = len(t.get_leaves())
                if num == len(t.get_species()) and num > min:
                    if num not in counter:
                        counter[num] = {}
                    for i in range(min, num + 1):
                        if i not in counter:
                            counter[i] = {}
                        counter[i][code] = set(t.get_leaf_names())
    statsFile = concatDir + "/stats.txt"
    outfile = open(statsFile, "w")
    outfile.write("Number of species\tNumber of single gene trees\tLength\n")
    limit = False
    keys = counter.keys()
    keys = sorted(keys)
    keys.reverse()
    for i in keys:
        if i >= min and not limit:
            workingDir = cAlgsDir + "/" + str(i)
            create_folder(workingDir)
            for code in counter[i]:
                create_parsed_alg_file(
                    code, counter[i][code], pathsFile, workingDir, readalPath
                )
            concatFileName = concatDir + "/concatenated_" + str(i) + ".fasta"
            concatenate_alignments_from_path(
                workingDir, concatFileName, spe2age, readalPath
            )
            phy_file = concatDir + "/concatenated_" + str(i) + ".phy"
            with open(phy_file) as c:
                first = c.readline()
            length = first.strip().split()[1]
            outfile.write(str(i) + "\t" + str(len(counter[i])) + "\t" + length + "\n")
            if len(counter[i]) > at_least or int(length) > max_length:
                limit = True
    outfile.close()


# Script used to generate concatenated alignments
def build_concatenated_alg(
    trees121File, spe2age, pathsFile, outDir, readalPath, prop=0.9, at_least=100, max_length = 50000
):
    """Align those 121 gene families where there at least `prop` species are represented.

    Parameters
    ----------
    trees121File : str
        File of 121 trees created with ap.obtain_121_trees().
    spe2age : dict
        Species 2 age dictionary
    pathsFile : str
        Pathfile created with create_pathfile()
    outDir : str
        output directory
    readalPath : str
        path of readal exe
    prop : num
        a decimal from 0 to 1. This represent
    at_least : num
        Stop the concatenation when there are at least `at_least` gene families concatenated.
    max_length : num
        Stop the concatenation when the aln is `max_length` positions long.

    Returns
    -------
    type
        A outDir/concatenated dir where the alignments, the concatnated fasta and phy and some stats are stored.

    """
    if not 0 <= prop <=1:
        sys.exit("Prop must be between 0 and 1")

    concatDir = outDir + "/concatenated/"
    if os.path.exists(concatDir):
        print("Concatenation directory already exists")
        fastaFiles = [x for x in glob.glob(concatDir + "/*fasta")]
        if len(fastaFiles) == 1:
            fastaFile = fastaFiles[0]
        else:
            fastaFiles = sorted(
                fastaFiles,
                key=lambda x: int(x.split("/")[-1].split("_")[1].split(".")[0]),
            )
            fastaFile = fastaFiles[0]
    else:
        # Create paths where data will be stored
        cAlgsDir = concatDir + "/algs/"
        create_folder(concatDir)
        create_folder(cAlgsDir)
        # codes = []
        counter = {}
        # Calculate how many species can be missing
        # if not extra:
        threshold = int(len(spe2age) - len(spe2age) * (1 - prop) + 1)
        # else:
        #     # If any extra concatenation method is called, just keep values over 5 is there are more than 5 species
        #     if len(spe2age) < 5:
        #         threshold = len(spe2age)
        #     else:
        #         threshold = 5
        # Obtain the trees with 121 relationships with a number of leaves above the threshold, assign them to the leafNumber and lower leafNumbers
        for line in open(trees121File):
            line = line.strip()
            dades = line.split()
            num = int(dades[1])
            if num >= threshold:
                if num not in counter:
                    counter[num] = set([])
                for i in range(threshold, num + 1):
                    if i not in counter:
                        counter[i] = set([])
                    counter[i].add(dades[0])
        # Open stats file
        statsFile = concatDir + "/stats.txt"
        outfile = open(statsFile, "w")
        outfile.write("Number of species\tNumber of single gene trees\tLength\n")
        limit = False
        keys = counter.keys()
        keys = sorted(keys)
        keys.reverse()
        # Load paths
        pathsList = loadPaths(pathsFile, "alg_aa")
        # For each of the trees for the maximum number of leaves
        for i in keys:
            # outfile.write(str(i) + "\t" + str(len(counter[i])) + "\n")
            if limit:
                pass
            else:
                workingDir = cAlgsDir + "/" + str(i)
                create_folder(workingDir)
                # For each alignment copy the alignment in the working folder
                for code in counter[i]:
                    if code in pathsList.keys():
                        cmd = "cp " + pathsList[code] + " " + workingDir
                        run_command(cmd, True)
                concatFileName = concatDir + "/concatenated_" + str(i) + ".fasta"
                concatenate_alignments_from_path(
                    workingDir, concatFileName, spe2age, readalPath
                )
                fastaFile = concatDir + "/concatenated_" + str(i) + ".fasta"
                phy_file = concatDir + "/concatenated_" + str(i) + ".phy"
                with open(phy_file) as c:
                    first = c.readline()
                length = first.strip().split()[1]
            # If enough alignments have been concatenated then the lower leaf numbers will not be processed
                outfile.write(str(i) + "\t" + str(len(counter[i])) + "\t" + length + "\n")
            if len(counter[i]) > at_least or int(length) > max_length:
                limit = True
        outfile.close()
    # return fastaFile


def get_orthologous_tree(t, code):
    subTrees = t.get_speciation_trees(map_features=["dist"])
    valid_subtree = set([])
    for subT in subTrees[2]:
        if code in subT.get_leaf_names() and len(subT.get_leaf_names()) > 1:
            subT.set_species_naming_function(load_species_name)
            subT = collapse_all_lineage_specific_expansions(subT, code)
            if subT:
                if len(subT.get_species()) == len(subT.get_leaves()):
                    if len(valid_subtree) == 0:
                        valid_subtree = subT.get_leaf_names()
                    else:
                        if len(valid_subtree) < len(subT.get_leaves()):
                            valid_subtree = subT.get_leaf_names()
    return valid_subtree



def build_extra_concatenated_alg3(
    treeFile, pathsFile, outDir, readalPath, spe2age, min=5, at_least=500, max_length=50000
):
    """Concatenate alignments from orthologous subtrees.

    Parameters
    ----------
    treeFile : str
        best_trees file from phylomeDB
    pathsFile : str
        Pathfile created with create_pathfile()
    outDir : str
        output directory
    readalPath : str
        path of readal exe
    spe2age : dict
        Speceis 2 age dictionary
    min : num
        Consider only trees with at least min species
    at_least : num
        Stop the concatenation when there are at least `at_least` gene families concatenated.
    max_length : num
        Stop the concatenation when the aln is `max_length` positions long.

    Returns
    -------
    type
        A outDir/concatenated_extra3 dir where the alignments, the concatnated fasta and phy and some stats are stored.

    """
    # Create output folders
    concatDir = outDir + "/concatenated_extra3/"
    create_folder(concatDir)
    # concatDir = concatDir + "/orthologous_trees/"
    # create_folder(concatDir)
    cAlgsDir = concatDir + "/algs/"
    create_folder(cAlgsDir)
    # Collapse species specific duplications and obtain list of one-to-one trees
    counter = {}
    for line in open(treeFile):
        line = line.strip()
        dades = line.split()
        code = dades[0]
        t = load_tree(dades[-1])
        if code in t.get_leaf_names():
            listLeaves = get_orthologous_tree(t, code)
            if len(listLeaves) >= min:
                num = len(listLeaves)
                if num not in counter:
                    counter[num] = {}
                for i in range(min, num + 1):
                    if i not in counter:
                        counter[i] = {}
                    counter[i][code] = listLeaves
    statsFile = concatDir + "/stats.txt"
    outfile = open(statsFile, "w")
    outfile.write("Number of species\tNumber of single gene trees\tLength\n")
    limit = False
    keys = counter.keys()
    keys = sorted(keys)
    keys.reverse()
    for i in keys:
        if i >= min and not limit:
            workingDir = cAlgsDir + "/" + str(i)
            create_folder(workingDir)
            for code in counter[i]:
                create_parsed_alg_file(
                    code, counter[i][code], pathsFile, workingDir, readalPath
                )
            concatFileName = concatDir + "/concatenated_" + str(i) + ".fasta"
            concatenate_alignments_from_path(
                workingDir, concatFileName, spe2age, readalPath
            )
            phy_file = concatDir + "/concatenated_" + str(i) + ".phy"
            with open(phy_file) as c:
                first = c.readline()
            length = first.strip().split()[1]
            outfile.write(str(i) + "\t" + str(len(counter[i])) + "\t" + length + "\n")
            if len(counter[i]) > at_least or int(length) > max_length:
                limit = True
    outfile.close()
