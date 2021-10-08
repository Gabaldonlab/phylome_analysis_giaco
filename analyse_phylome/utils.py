# -*- coding: utf-8 -*-
import subprocess as sp
import sys
import os
import ete3
import glob
import pandas as pd
from ftplib import FTP

# Run bash commands
def run_command(cmd, ommit=False):
    if ommit:
        try:
            process = sp.Popen(cmd, shell=True)
        except:
            pass
        process.communicate("Y\n")
        if process.wait() != 0:
            print("Error occurred, but you chose to ommit it")
    else:
        try:
            process = sp.Popen(cmd, shell=True)
        except OSError:
            sys.exit("Error: Execution cmd failed")
        process.communicate("Y\n")
        if process.wait() != 0:
            sys.exit("ERROR: Execution cmd failed")


def run_command_with_return(cmd):
    try:
        process = sp.Popen(cmd, shell=True, stdout=sp.PIPE).stdout
    except:
        sys.exit()
    return process.readlines()


# Create a folder
def create_folder(name):
    if not os.path.exists(name):
        try:
            os.mkdir(name)
        except:
            print("Unable to create directory " + name)


def loadPaths(pathsFile, tag):
    """Loads paths from a pathfile based on the tag parameter.

    Parameters
    ----------
    pathsFile : str
        File created with create_pathfile().
    tag : str
        Tag of the files of interest

    Returns
    -------
    type
        Dictionary of gene:corresponding_file

    """
    pathsList = {}
    for line in open(pathsFile):
        line = line.strip()
        dades = line.split("\t")
        if dades[0] == tag:
            pathsList[dades[1]] = dades[2]
    return pathsList


def build_sp2age(sptree, seed):
    """Build a species to age dictionary given a tree and its main seq.

    Parameters
    ----------
    sptree : tree
        Phylogenetic tree of interest
    seed : str
        Main seq where the dictionary will start.

    Returns
    -------
    type
        Returns a dict with species:age_to_main as key:value

    """
    spe2age = {}
    node = sptree.get_leaves_by_name(seed)[0]
    age = 1
    while not node.is_root():
        for leaf in node.get_leaves():
            if leaf.name not in spe2age:
                spe2age[leaf.name] = age
        age += 1
        node = node.up
    for leaf in sptree.get_leaves():
        if leaf.name not in spe2age:
            spe2age[leaf.name] = age
    return spe2age


def load_sequences(contigFile, delimiter):
    """Load sequences into a variable.

    Parameters
    ----------
    contigFile : str
        File where the seqs are stored.
    delimiter : str
        Delimiter of sequences.

    Returns
    -------
    type
        Variable with sequences stored.

    """
    seqs = {}
    name = ""
    s = []
    for line in open(contigFile):
        line = line.strip()
        if ">" in line:
            if name != "":
                seqs[name] = "".join(s)
            if delimiter == "":
                name = line.replace(">", "")
            else:
                name = line.replace(">", "").split(delimiter)[0]
            s = []
        else:
            s.append(line.upper())
    if len(s) != 0:
        seqs[name] = "".join(s)
    return seqs


def load_species_name(node):
    return node.split("_")[1]


def load_species_name_whole(node):
    return node


def load_tree(tree):
    t = ete3.PhyloTree(tree, sp_naming_function=load_species_name)
    return t


def root_tree(t, spe2age=None, midpoint=False):
    """Root a Phylogenetic tree.

    Parameters
    ----------
    t : Tree
        Tree to root.
    spe2age : dict
        Species2age dictionary if available.
    midpoint : bool
        Use midpoint rooting.

    Returns
    -------
    tree
        Rooted phylogenetic tree.

    """
    if spe2age is not None:
        try:
            t.set_outgroup(t.get_farthest_oldest_leaf(spe2age))
        except:
            sys.exit("Something went wrong with sp2age dict!")
    elif midpoint:
        t.set_outgroup(t.get_midpoint_outgroup())
    else:
        sys.exit("Something went wrong with rooting!")
    return t


def root_species_tree(t, spe2age=None, midpoint=False, out_list=None, force=False):
    """Root a species tree.

    Parameters
    ----------
    t : tree
        Species tree to root.
    spe2age : dict
        Species 2 age dictionary if available.
    midpoint : bool
        Do midpoint rooting.
    out_list : list
        Lisst of outgroup species.
    force:
        if furthest species from species 2 age dictionary remove one and try again until a monophyletic outgroup can be set. The furthest species are removed in order of topological distance from seed (closest are removed first!). May not have much sense, therefore use it with caution!
    Returns
    -------
    tree
        Rooted species tree.

    """
    p = ete3.PhyloTree(t.write(), sp_naming_function=None)
    for n in p.get_leaves():
        n.species = n.name
    if spe2age is not None:
        furthest = [k for k, v in spe2age.items() if v == max(spe2age.values())]
        if len(furthest) > 1:
            print(
                "there are more than one furthest sequence, the root will be the ancestor node that comprise all of them."
            )
            print(furthest)
            if p.get_common_ancestor(furthest) == p.get_tree_root():
                if force:
                    red = sort_furthest_by_dist(spe2age, p)
                    num = len(furthest)
                    while num > 0:
                        red = red[1:]
                        if p.get_common_ancestor(red) == p.get_tree_root():
                            num = num - 1
                            continue
                        else:
                            print("Found monophyletic outgroup:")
                            print(red)
                            p.set_outgroup(p.get_common_ancestor(red))
                            break
                else:
                    raise ValueError(
                        "Could not find a monophyletic clade containing all furthest species, try midppoint rooting or give a list of outgroup"
                    )
            # p.set_outgroup(p.get_midpoint_outgroup())
            else:
                p.set_outgroup(p.get_common_ancestor(furthest))
        else:
            print("only one furthest sequence, this will be the outgroup")
            print(furthest)
            p.set_outgroup(p.get_farthest_oldest_leaf(spe2age))
    if midpoint:
        p.set_outgroup(p.get_midpoint_outgroup())
    if out_list is not None:
        if len(out_list) == 1:
            p.set_outgroup(out_list[0])
        else:
            p.set_outgroup(p.get_common_ancestor(out_list))
    return p


def sort_furthest_by_dist(s2a, tree):
    furthest = [k for k, v in s2a.items() if v == max(s2a.values())]
    seed = [k for k, v in s2a.items() if v == min(s2a.values())]
    dist_furthest = {}
    for el in furthest:
        dest = tree & el
        arr = tree & seed[0]
        dist_furthest[el] = arr.get_distance(dest)
    furthest = sorted(furthest, key=dist_furthest.get)
    return furthest


def print_sequence(code, sequence, outfile):
    outfile.write(">" + code + "\n")
    i = 0
    if sequence[-1] == "*":
        sequence = sequence[:-1]
    while i < len(sequence):
        outfile.write(sequence[i : i + 60] + "\n")
        i += 60


def create_pathfile(pathsFile, dir, tag):
    """Create a file with paths for a file in a specific dir.

    Parameters
    ----------
    pathsFile : str
        outputfile
    dir : str
        Directory containing the files of interest
    tag : str
        Tag for the type of file. For now only "alg_aa" is supported though.

    Returns
    -------
    type
        Returns a file with paths for each object of interest.

    """
    outfilePath = open(pathsFile, "w")
    for firstDir in glob.glob(dir + "/*"):
        code = firstDir.split("/")[-1].split(".")[0]
        outfilePath.write(tag + "\t" + code + "\t" + firstDir + "\n")
    outfilePath.close()


# get speciesrax data

# set alignment args as optional way better than this!
# mapping and family file not imporant


def get_generax_data(gene_trees, out_dir, aln_dir=None, keep_model=True):
    """Prepare data for generax or speciesrax.

    Parameters
    ----------
    gene_trees : str
        Best trees file from phylomeDB.
    out_dir : str
        output directory.
    aln_dir : str
        directory where the alns are stored. It can be ignored if the alns are not useful.
    keep_model : bool
        Keep the best model found in phylomeDB

    Returns
    -------
    type
        A directory with the structure and files needed for GeneRax to run.

    """

    gene_file = os.path.basename(gene_trees)
    gene_noex = os.path.splitext(gene_file)[0]
    gene_dir = os.path.dirname(gene_trees)

    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
    if not os.path.isdir(out_dir + "/trees/"):
        os.mkdir(out_dir + "/trees/")
    if not os.path.isdir(out_dir + "/mapping/"):
        os.mkdir(out_dir + "/mapping/")

    out_map = out_dir + "/mapping/mapping"
    out_family = out_dir + "/family.txt"

    all_leaves = []

    with open(out_family, "w") as f:
        f.write("[FAMILIES]\n")
        with open(gene_trees, "r") as g:
            for line in g:

                line = line.strip().split()
                id = line[0]
                tree = ete3.Tree(line[3])

                leaves = [leaf.name for leaf in tree.iter_leaves()]
                species = list(set(leaf.split("_")[1] for leaf in leaves))

                sp_gene_dict = {}

                for s in species:
                    sp_genes = []
                    for gene in leaves:
                        nm = gene.split("_")
                        if nm[1] == s:
                            sp_genes.append(gene)
                    sp_gene_dict[s] = sp_genes

                out_mapfile = out_map + "_" + id + ".txt"
                with open(out_mapfile, "w") as m:
                    for k in sp_gene_dict.keys():
                        if len(sp_gene_dict[k]) > 0:
                            m.write(str(k) + ":" + ";".join(sp_gene_dict[k]) + "\n")

                f.write("- " + id + "\n")
                id_file = out_dir + "/trees/" + id + ".gene.newick"
                tree.write(outfile=id_file)
                f.write("starting_gene_tree = " + id_file + "\n")
                if keep_model:
                    model = line[1]
                    f.write("subst_model = " + model + "\n")
                if aln_dir is not None:
                    aln_file = aln_dir + "/" + id + ".clean.fasta"
                    aln_file_u = aln_dir + "/" + id + ".clean_noU.fasta"
                    if os.path.exists(aln_file_u):
                        f.write("alignment = " + aln_file_u + "\n")
                    elif os.path.exists(aln_file):
                        f.write("alignment = " + aln_file + "\n")
                    else:
                        sys.exit("could not find file: " + aln_file)
                f.write("mapping = " + out_mapfile + "\n")


def scan_for_Us(aln_dir, what="U", replace=False, with_what="C"):
    """Scan for U (Selenocysteine) characters in alignment. useful
    if using generax as alignments containing Us won't be parsed.

    Parameters
    ----------
    aln_dir : str
        Directory where .
    what : str
        Character to find (U is default).
    replace : bool
        If true create a file with no Us named with _noU added.

    Returns
    -------
    type
        Returns if the alignments contain non canonical aminoacid and eventually a file with those aminoacid replaced.

    """
    if os.path.exists(aln_dir):
        is_u_ever = False
        for file in os.listdir(aln_dir):
            if file.endswith("clean.fasta"):
                toread = os.path.join(aln_dir, file)
                is_U = False
                cmd = "grep -h -v '>' " + toread + " | grep " + what
                files_u = run_command_with_return(cmd)
                if len(files_u) > 0:
                    is_U = True
                    is_u_ever = True
                    if replace:
                        outname = "clean_no" + what
                        towrite = toread.replace("clean", outname)
                        cmd = (
                            "sed '/^>/! {s/"
                            + what
                            + "/"
                            + with_what
                            + "/g}' "
                            + toread
                            + " > "
                            + towrite
                        )
                        run_command(cmd)
                    if is_U:
                        print(file + " Contains " + what)
        if not is_u_ever:
            print("No files contain " + what)


def get_astral_pro_data(gene_trees, out_file):
    """Obtain data to run astral-pro.

    Parameters
    ----------
    gene_trees : str
        Best trees file from PhylomeDB.
    out_file : str
        outputfile.

    Returns
    -------
    file
        Write outputfile in specified directory.

    """

    with open(out_file, "w") as o:
        with open(gene_trees) as t:
            for line in t:
                line = line.split()
                tree = ete3.Tree(line[3])
                for leaf in tree.iter_leaves():
                    leaf.name = leaf.name.split("_")[1]
                string = tree.write()
                o.write(string + "\n")


def get_all_species(fileTree):
    """Get all mmnemonics code in best trees file in phylomedb.

    Parameters
    ----------
    fileTree : str
        best tree file from PhylomeDB

    Returns
    -------
    list
        List of mnemonics codes found in file.

    """
    set_sp = set()
    with open(fileTree) as t:
        for line in t:
            line = line.split()
            tree = ete3.Tree(line[3])
            for leaf in tree.iter_leaves():
                leaf.name = leaf.name.split("_")[1]
                set_sp.add(leaf.name)
    return list(set_sp)


def get_all_species_counts(fileTree):
    """Get all mmnemonics code  and number of trees in which the species is present in best trees file in phylomedb.

    Parameters
    ----------
    fileTree : str
        best tree file from PhylomeDB

    Returns
    -------
    list
        dictionary of mnemonics codes found in file + number of trees.

    """
    dict_occurence = {}
    with open(fileTree) as t:
        for line in t:
            line = line.split()
            tree = ete3.Tree(line[3])
            set_sp = set([name.split("_")[1] for name in tree.get_leaf_names()])
            for leaf in set_sp:
                if dict_occurence.get(leaf) is None:
                    dict_occurence[leaf] = 0
                dict_occurence[leaf] += 1
    return dict_occurence


def normalize_counts(counts_df, norm_dict):
    df = counts_df[counts_df.index.isin(norm_dict.keys())]
    # df = df[~df.index.str.contains("node_")]
    df = df.apply(lambda value: value / norm_dict[value.name], axis=1)
    return df


def get_all_models_counts(treefile):
    with open(treefile) as t:
        models = [line.split()[1] for line in t]
        dict_counts = dict((x, models.count(x)) for x in set(models))
        return dict_counts


def get_species_name(node_name_string):
    spcode = node_name_string.split("_")[1]
    return spcode


def get_ecce_data(
    gene_trees, species_tree, out_dir, root=False, midpoint=False, ultrametric=False
):
    """Build data useful for ecceTERA.

    Parameters
    ----------
    gene_trees : str
        Best trees file from PhylomeDB.
    species_tree : str
        Species tree
    out_dir : str
        output directory.
    root : bool
        Root the gene trees before comparing.
    midpoint : bool
        If root=True root by midpoint? If false the trees will be rooted with spe2age dictionary
    ultrametric : bool
        Convert the species tree to ultrametric (to run dated analysis)

    Returns
    -------
    type
    Data needed by ecceTERA in outdir ecce_best_trees file.

    """

    gene_file = os.path.basename(gene_trees)
    # gene_noex = os.path.splitext(gene_file)[0]
    # gene_dir = os.path.dirname(gene_trees)

    sp_file = os.path.basename(species_tree)
    # sp_noex = os.path.splitext(sp_file)[0]
    # sp_dir = os.path.dirname(species_tree)
    sp = ete3.Tree(species_tree)

    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    if ultrametric:
        sp.convert_to_ultrametric()
        sp.write(outfile=out_dir + "/ecce_ultra_" + sp_file, format=5)

    sp_names = sp.get_leaf_names()
    out_gene = out_dir + "/ecce_" + gene_file

    if root and not midpoint:
        with open(gene_trees) as g:
            first_line = g.readline().strip().split()
            seed = first_line[0].split("_")[1]
            sp2age = build_sp2age(sp, seed)

    with open(out_gene, "a") as o:
        with open(gene_trees) as g:
            for line in g:
                line = line.strip().split()
                tree = line[3]
                if root:
                    t = ete3.PhyloTree(tree, sp_naming_function=load_species_name)
                    if midpoint:
                        t.set_outgroup(t.get_midpoint_outgroup())
                    else:
                        t.set_outgroup(t.get_farthest_oldest_leaf(sp2age))
                    # t.resolve_polytomy()
                else:
                    t = ete3.Tree(tree)
                    # t.resolve_polytomy()
                for spec in sp_names:
                    num = 1
                    for leaf in t.iter_leaves():
                        nm = leaf.name.split("_")
                        if spec == nm[1]:
                            leaf.name = nm[1] + "_" + str(num) + "_" + nm[0]
                            num = num + 1
                string = t.write()
                o.write(string + "\n")


def get_tax_dict_uniprot(tree, uniprot_df):  # add solve species??
    """Get taxonomic code from mnemonic starting from uniprot speclist file.

    Parameters
    ----------
    tree : tree
        Tree with mnemo code as leaf names.
    uniprot_df : str
        file of uniprot speclist. It can be retireved with: "curl https://www.uniprot.org/docs/speclist.txt -o data/speclist_$(date +'%d_%m_%y').txt"
    Returns
    -------
    dict
        Dictionary with mnemos as key and ncbi taxid as values.

    """
    ncbi = ete3.NCBITaxa()
    sp = [nm for nm in tree.get_leaf_names()]
    sp = list(set([s.split("_")[-1] if "_" in s else s for s in sp]))
    sp_str = [el for el in sp if not el.isdecimal()]
    sp_num = [el for el in sp if el.isdecimal()]

    with open(uniprot_df) as f:
        sp_line = [line.split() for line in f.readlines()]
        sp_line = [line for line in sp_line if line != []]
        taxo_dict = {
            line[0]: line[2].replace(":", "") for line in sp_line if line[0] in sp_str
        }

    for el in sp_num:
        if len(ncbi.get_taxid_translator([el])) > 0:
            taxo_dict[el] = el

    absent = set([abs for abs in sp if abs not in taxo_dict.keys()])

    if absent:
        print("Warning: " + " ".join(absent) + " was not found in the dictionary.")
        print("You may want to add it manually")

    return taxo_dict


def get_tax_dict_info(info_txt):
    """Get taxonomic code from phylomedb info file.

    Parameters
    ----------
    info_txt : str
        info.txt file from phylomeDB ftp server.
    Returns
    -------
    dict
        Dictionary with mnemos as key and ncbi taxid as values.

    """
    with open(info_txt) as i:
        flag = False
        lines = []
        for line in i:
            if "TaxaID" in line:
                flag = True
            if flag:
                lines.append(line)

    tax_dict = {el.split()[1].split(".")[0]: el.split()[0] for el in lines[2:]}
    return tax_dict


def get_taxonomic_df(tax_dict, set_cols=False, fill=False):
    ncbi = ete3.NCBITaxa()
    tax_resolved = ncbi.get_taxid_translator([sp for sp in tax_dict.values()])
    tax_dict = {k: v for k, v in tax_dict.items() if int(v) in tax_resolved.keys()}
    mnemo_sp = {}
    for mnemo in tax_dict:
        id = int(tax_dict[mnemo])
        try:
            species = tax_resolved[id]
        except KeyError:
            species = mnemo
        mnemo_sp[mnemo] = species

    whole_tax_dict = {}

    for key, value in tax_resolved.items():
        # sp_name = "".join(ncbi.get_taxid_translator([key]).values())
        lineage = ncbi.get_lineage(key)
        names = ncbi.get_taxid_translator(lineage)
        rank = ncbi.get_rank(lineage)
        ordered_names = [names[taxid] for taxid in lineage]
        ordered_clades = [rank[taxid] for taxid in lineage]
        taxonomy = list(zip(ordered_names, ordered_clades))
        # d = {k:{rank[k]:lineage[k]} for k in lineage.keys()}
        whole_tax_dict[value] = taxonomy

    if set_cols:
        common_clades = [
            "species",
            "genus",
            "family",
            "order",
            "class",
            "phylum",
            "kingdom",
            "superkingdom",
        ]
    else:
        seen = []

        for sp in whole_tax_dict:
            set_single = set()
            for id in whole_tax_dict[sp]:
                set_single.add(id[1])
                #     for clade in whole_tax_dict[sp][id]:
                #         set_single.add(clade)
                seen.append(set_single)

        common_clades = list(set.intersection(*seen) - set(["no rank", "clade"]))

    for key in whole_tax_dict:
        if not set_cols:
            new_tuple = [el for el in whole_tax_dict[key] if el[1] in common_clades]
        else:
            clades_present = [
                el[1] for el in whole_tax_dict[key] if el[1] in common_clades
            ]
            memory = ()
            new_tuple = []
            for clade in common_clades[::-1]:
                if clade in clades_present:
                    toadd = [el for el in whole_tax_dict[key] if el[1] == clade][0]
                    memory = toadd
                elif fill and len(memory) > 0:
                    toadd = (memory[0], clade)
                else:
                    toadd = ("", clade)
                new_tuple.append(toadd)
        whole_tax_dict[key] = new_tuple
    return whole_tax_dict


def get_dataframe(whole_tax_dict):
    data = []
    for key in whole_tax_dict:
        row = [el[0] for el in whole_tax_dict[key]]
        col_names = [el[1] for el in whole_tax_dict[key]][::-1]
        data.append(row[::-1])

    df = pd.DataFrame(data, columns=col_names)
    return df


def annotate_taxonomy(sptree, whole_tax_dict, tax_dict):
    out_sptree = sptree.copy("cpickle")
    ncbi = ete3.NCBITaxa()
    tax_resolved = ncbi.get_taxid_translator([sp for sp in tax_dict.values()])
    tax_dict = {k: v for k, v in tax_dict.items() if int(v) in tax_resolved.keys()}
    mnemo_sp = {}
    for mnemo in tax_dict:
        id = int(tax_dict[mnemo])
        try:
            species = tax_resolved[id]
        except KeyError:
            species = mnemo
        mnemo_sp[mnemo] = species
    for node in out_sptree.iter_leaves():
        if node.name in mnemo_sp.keys():
            node.species = mnemo_sp[node.name]
        else:
            node.species = node.name

    all_phyla = set()

    no_species_dict = {}

    for key in whole_tax_dict:
        no_sp = [el for el in whole_tax_dict[key] if el[1] != "species"]
        no_species_dict[key] = no_sp

    for key in whole_tax_dict:
        for el in no_species_dict[key]:
            all_phyla.add(el[0])

    color_dict = {}
    colors = ete3.random_color(num=len(set(all_phyla)))
    color_dict = {el[1]: el[0] for el in list(zip(colors, list(all_phyla)))}
    color_dict[""] = "white"

    for key in no_species_dict:
        new_tuple = [el + (color_dict[el[0]],) for el in no_species_dict[key]]
        no_species_dict[key] = new_tuple

    num_clades = [len(no_species_dict[key]) for key in no_species_dict][0]

    for node in out_sptree.iter_leaves():
        if node.name in tax_dict.keys():
            sp_df = no_species_dict[node.species][::-1]
        else:
            sp_df = [("", "", "")] * num_clades
        node.add_feature("col_df", sp_df)
    return out_sptree


def annotate_boxes(taxo_sptree, whole_tax_dict, target=5):

    clade_dict = {}
    for tax in whole_tax_dict.values():
        for el in tax:
            if el[1] not in clade_dict:
                clade_dict[el[1]] = set()
            clade_dict[el[1]].add(el[0])

    num_dict = {el: len(clade_dict[el]) for el in clade_dict}
    key, value = min(num_dict.items(), key=lambda kv: abs(kv[1] - target))

    out_sptree = taxo_sptree.copy("cpickle")

    classes = set()
    for node in out_sptree:
        for row in node.col_df:
            if key in row:
                node.tax_class = row[0]
                classes.add(node.tax_class)

    color_dict = {}
    colors = ete3.random_color(num=len(classes))
    color_dict = {el[1]: el[0] for el in list(zip(colors, list(classes)))}

    for tax_class in list(classes):
        nodes = []
        for node in out_sptree:
            if node.species not in whole_tax_dict:
                node.tax_class = ""
            elif node.tax_class == tax_class:
                nodes.append(node.name)
        # nodes = [
        #     node.name for node in out_sptree node.tax_class == tax_class
        # ]
        nst = ete3.NodeStyle()
        nst["bgcolor"] = color_dict[tax_class]
        if len(nodes) > 1:
            mrca = out_sptree.get_common_ancestor(nodes)
            mrca.set_style(nst)
        else:
            node = out_sptree & nodes[0]
            node.set_style(nst)
    return out_sptree, color_dict


def layout_species_taxon(node):
    width = 100  # try to find a mulitplicator or something
    # If node is a leaf, add the nodes name and a its scientific
    # name
    node.img_style["size"] = 0
    if node.is_leaf():
        name_face = ete3.faces.AttrFace("species")
        ete3.faces.add_face_to_node(name_face, node, column=0, position="branch-right")
        col_idx = 1
        for clade in node.col_df:
            rect = ete3.faces.RectFace(
                width,
                20,
                bgcolor=clade[2],
                fgcolor=clade[2],
                label={"text": clade[0], "color": "white", "fontsize": 6},
            )
            ete3.faces.add_face_to_node(rect, node, column=col_idx, aligned=True)
            col_idx += 1


def layout_species_std(node):
    node.img_style["size"] = 0
    if node.is_leaf():
        name_face = ete3.faces.AttrFace("species")
        ete3.faces.add_face_to_node(name_face, node, column=0, position="branch-right")


def viz_species_tree(
    sptree,
    legend=None,
    taxonomy=False,
    circular=False,
    show=True,
    render=None,
    bs=False,
):
    ts = ete3.TreeStyle()
    ts.show_leaf_name = False
    # ts.allow_face_overlap = True
    ts.draw_aligned_faces_as_table = True
    if legend is not None:
        ts.legend_position = 4
        for el in legend.items():
            ts.legend.add_face(
                ete3.faces.RectFace(width=10, height=10, fgcolor=el[1], bgcolor=el[1]),
                column=0,
            )
            ts.legend.add_face(ete3.faces.TextFace(" " + el[0]), column=1)
    if taxonomy:
        ts.layout_fn = layout_species_taxon
    else:
        ts.layout_fn = layout_species_std
    if circular:
        ts.mode = "c"
    if show:
        sptree.show(tree_style=ts)
    if bs:
        ts.show_branch_support = True
    if render is not None:
        sptree.render(render, tree_style=ts)


def annotate_genetree(genetree, taxo_dict):
    """Annotate genetree in order to visualize the events.

    Parameters
    ----------
    genetree : tree
        Gene tree.
    taxo_dict : dict
        Taxonomic dict obtained with get_tax_dict_*

    Returns
    -------
    tree
        tree where each leaf has a color based on taxonomic id

    """

    colors = ete3.random_color(num=len(set(taxo_dict)))

    mnemo_sp = {}
    num = 0
    for mnemo in taxo_dict:
        # id = int(tax_dict[mnemo])
        # species = tax_resolved[id]
        col = colors[num]
        num += 1
        mnemo_sp[mnemo] = col

    for node in genetree.iter_leaves():
        node.mnemo = node.name.split("_")[1]
        # node.species = mnemo_sp[node.mnemo][0]
        if node.mnemo in taxo_dict.keys():
            node.col = mnemo_sp[node.mnemo]
        else:
            node.col = "black"
    return genetree


def layout_grax_nhx(node):
    # If node is a leaf, add the nodes name and a its scientific
    # name
    if node.is_leaf():
        nameFace = ete3.faces.TextFace(node.name, fgcolor=node.col)
        ete3.faces.add_face_to_node(nameFace, node, column=0)
        node.img_style["size"] = 7
        node.img_style["shape"] = "circle"
        node.img_style["fgcolor"] = node.col
    if node.D != "N":
        node.img_style["size"] = 15
        node.img_style["shape"] = "circle"
        node.img_style["fgcolor"] = "darkred"
    elif node.H != "N":
        node.img_style["size"] = 15
        node.img_style["shape"] = "circle"
        node.img_style["fgcolor"] = "darkgreen"
        HGTFace = ete3.faces.TextFace("-".join(node.H.split("@")[1:]), fsize=5)
        ete3.faces.add_face_to_node(HGTFace, node, column=0, position="branch-bottom")
    else:
        node.img_style["size"] = 0


def viz_grax_tree(genetree, show=True, render=None):
    ts = ete3.TreeStyle()
    ts.show_leaf_name = False
    ts.show_branch_length = True
    ts.show_branch_support = True
    ts.layout_fn = layout_grax_nhx
    if show:
        genetree.show(tree_style=ts)
    if render is not None:
        genetree.render(render, tree_style=ts)


def read_per_sp_counts(events_file):
    df = pd.read_csv(events_file, sep=" ")
    df["#S"] = df["#S"].str.replace(r"S=", "")
    df["#SL"] = df["#SL"].str.replace(r"SL=", "")
    df["#D"] = df["#D"].str.replace(r"D=", "")
    df["#T"] = df["#T"].str.replace(r"T=", "")
    df["#TL"] = df["#TL"].str.replace(r"TL=", "")
    df[df.columns] = df[df.columns].apply(pd.to_numeric, errors="coerce")
    df.columns = [col.replace("#", "") for col in df.columns]
    return df


def annotate_sptree_counts(sptree, counts_df):
    out_sptree = sptree.copy("cpickle")
    for node in out_sptree.traverse():
        values = [int(el) for el in counts_df.loc[node.name]]
        node.values = values
    return out_sptree


def annotate_sptree_counts_both(sptree, counts_ecce, counts_grax, matches):
    out_sptree = sptree.copy("cpickle")
    for node in out_sptree.traverse():
        ecce_values = [int(el) for el in counts_ecce.loc[node.name]]
        node.ecce_values = ecce_values
        grax_key = matches[node.name]
        grax_values = [int(el) for el in counts_grax.loc[grax_key]]
        node.grax_values = grax_values
    return out_sptree


def layout_both_bars(node):
    palette = ["#ed4c67", "#f69e1f", "#1289a7"] + ["#ed4c67", "#f69e1f", "#1289a7"]
    values = node.ecce_values + node.grax_values[1:-1]
    # values = [j - i for i, j in zip(node.grax_values[1:-1], node.ecce_values)]
    B = ete3.faces.BarChartFace(values, colors=palette)  # SL,D,T!
    # Add node name to laef nodes
    ete3.faces.add_face_to_node(B, node, 0, position="branch-right")


def layout_pies(node):
    # if node.is_leaf():
    # values = [int(el) for el in df.loc[node.name][2:-1]]
    # norm_vals = [(val / sum(values))*100 for val in values]
    # labels = ["D", "T"]
    # else:
    palette = ["#ed4c67", "#f69e1f", "#1289a7"]
    if len(node.values) > 3:
        values = node.values[1:-1]
    else:
        values = node.values
    if sum(values) > 0:
        norm_vals = [(val / sum(values)) * 100 for val in values]
    else:
        norm_vals = [0, 0, 0]
    # B = faces.BarChartFace(values=values, labels=labels)
    P = ete3.faces.PieChartFace(norm_vals, 100, 100, colors=palette)  # SL,D,T!
    # Add node name to laef nodes
    # faces.add_face_to_node(B, node, 0)#, position="branch-right")
    ete3.faces.add_face_to_node(P, node, 0, position="branch-right")


def viz_pies_tree(sptree, circular=False, show=True, render=None):
    palette = ["#ed4c67", "#f69e1f", "#1289a7"]
    labels = ["SL", "D", "T"]
    legend = list(zip(labels, palette))

    ts = ete3.TreeStyle()
    ts.legend_position = 4
    ts.layout_fn = layout_pies
    if circular:
        ts.mode = "c"

    for el in legend:
        ts.legend.add_face(
            ete3.faces.RectFace(width=50, height=50, fgcolor=el[1], bgcolor=el[1]),
            column=0,
        )
        ts.legend.add_face(ete3.faces.TextFace(" " + el[0], fsize=20), column=1)
    if show:
        sptree.show(tree_style=ts)
    if render is not None:
        sptree.render(render, tree_style=ts)


def viz_both_bars_tree(sptree, circular=False, show=True, render=None):
    palette = ["#ed4c67", "#f69e1f", "#1289a7"]
    labels = ["L", "D", "T"]
    legend = list(zip(labels, palette))

    ts = ete3.TreeStyle()
    ts.legend_position = 4
    ts.layout_fn = layout_both_bars
    if circular:
        ts.mode = "c"

    for el in legend:
        ts.legend.add_face(
            ete3.faces.RectFace(width=50, height=50, fgcolor=el[1], bgcolor=el[1]),
            column=0,
        )
        ts.legend.add_face(ete3.faces.TextFace(" " + el[0], fsize=20), column=1)
    if show:
        sptree.show(tree_style=ts)
    if render is not None:
        sptree.render(render, tree_style=ts)


def annotate_sptree_ancsize(sptree, counts_df):
    out_sptree = sptree.copy("cpickle")
    # very ugly code!
    scale = 0
    for node in out_sptree.traverse():
        if not node.is_leaf():
            anc_size = sum([int(el) for el in counts_df.loc[node.name][0:1]])
            if anc_size > scale:
                scale = anc_size
    for node in out_sptree.traverse():
        if not node.is_leaf():
            anc_size = (
                sum([int(el) for el in counts_df.loc[node.name][0:1]]) / scale * 10
            )
            node.anc_size = anc_size
    return out_sptree, scale


# Comment 2: the best way to estimate the size of an ancestral genome is to count the number of S and SL events associated with this genome. Do not forget the SL events.


def layout_anc(node):
    node.img_style["size"] = 0
    if not node.is_leaf():
        # norm_vals = [(val / sum(values)) * 100 for val in values]
        C = ete3.faces.CircleFace(radius=node.anc_size, color="teal", style="sphere")
        # C.opacity = 0.1
        # Add node name to laef nodes
        ete3.faces.add_face_to_node(
            C, node, 0, position="float"
        )  # , position="branch-right")


def viz_ancsize_tree(sptree, scale, circular=False, show=True, render=None):
    ts_anc = ete3.TreeStyle()
    ts_anc.layout_fn = layout_anc
    if circular:
        ts_anc.mode = "c"
    ts_anc.legend_position = 4
    ts_anc.legend.add_face(
        ete3.faces.CircleFace(10, color="teal", style="sphere"), column=0
    )
    ts_anc.legend.add_face(ete3.faces.TextFace(" " + "S+SL: " + str(scale)), column=1)
    if show:
        sptree.show(tree_style=ts_anc)
    if render is not None:
        sptree.render(render, tree_style=ts_anc)


def rename_tree_tax(tree, taxo_dict):

    for leaf in tree.iter_leaves():
        mnemo = leaf.name.split("_")[-1]
        if mnemo in taxo_dict.keys():
            leaf.name = taxo_dict[mnemo]
        else:
            leaf.name = mnemo
    return tree


def analyze_rooting(spe2age):
    # fracs = {}
    num = len([k for k, v in spe2age.items() if v == max(spe2age.values())])
    den = len(spe2age)
    frac = str(num) + "/" + str(den)
    dec = num / den
    fracs = [frac, dec, num, den]

    return fracs


def get_ftp_stats(phylome_ids=None, all=False):
    """Get files and dimensions for a phylome stored in FTP server.

    Parameters
    ----------
    phylome_ids : list
        List of phylome ids (if not given all must be set to True).
    all : bool
        If set to True it will return the name-size dictionary for all phylomes in ftp.

    Returns
    -------
    dict
        A dicitonary where filenames are key and byte size are values.

    """
    ftp = FTP("ftp.phylomedb.org")
    ftp.login()
    ftp.cwd("phylomedb")
    ftp.cwd("phylomes")

    if not all:
        subdir = ["phylome_" + str(format(int(num), "04d")) for num in phylome_ids]
    else:
        subdir = [id for id in ftp.nlst() if "phylome_" in id]

    ftp_dict = {}

    for dir in subdir:
        ftp.cwd(dir)
        dim_id = {}
        for el in ftp.nlst():
            ftp.sendcmd("TYPE i")  # Switch to Binary mode
            dim = ftp.size(el)
            dim_id[el] = dim
        ftp_dict[dir] = dim_id
        ftp.cwd("../")
    ftp.close()

    return ftp_dict


def get_ftp_files(
    phylome_id, outdir=None, to_get=["best_trees", "all_algs", "phylome_info"]
):
    """Get files stored in ftp directory for specific phylome.

    Parameters
    ----------
    phylome_id : str or int
        Phylome id.
    outdir : str
        Output directory where files will be stored. If not given it will be the same as in ftp server.
    to_get : list
        List of file names (without extension) to get from ftp. Default is to_get=["best_trees", "all_algs", "phylome_info"].

    Returns
    -------
    type
        An outdir with the requested files if present in ftp.

    """
    ftp = FTP("ftp.phylomedb.org")
    ftp.login()
    ftp.cwd("phylomedb")
    ftp.cwd("phylomes")

    subdir = "phylome_" + str(format(int(phylome_id), "04d"))

    ftp.cwd(subdir)

    files_in_ftp = ftp.nlst()

    if outdir is None:
        outdir = subdir

    create_folder(outdir)

    for file in files_in_ftp:
        # you could add warning if a requested file is not in nlst
        file_nm = file.split(".")[0]
        if file_nm in to_get:
            local_filename = os.path.join(outdir, file)
            print("file: " + file + " found, will be written in " + local_filename)
            if os.path.isfile(local_filename):
                print(local_filename + " already present!")
            else:
                with open(local_filename, "wb") as f:
                    ftp.retrbinary("RETR " + file, f.write)

    ftp.close()
