import subprocess as sp
import sys
import os
import ete3
import glob
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


def root_species_tree(t, spe2age=None, midpoint=False, out_list=None):
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
                print(
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
                    aln_file_u = aln_dir + "/" + id + ".clean_noUs.fasta"
                    if os.path.exists(aln_file_u):
                        f.write("alignment = " + aln_file_u + "\n")
                    elif os.path.exists(aln_file):
                        f.write("alignment = " + aln_file + "\n")
                    else:
                        sys.exit("could not find file: " + aln_file)
                f.write("mapping = " + out_mapfile + "\n")


def scan_for_Us(aln_dir, replace=False, with_what="C"):  # , replace=False):
    if os.path.exists(aln_dir):
        is_u_ever = False
        for file in os.listdir(aln_dir):
            if file.endswith("clean.fasta"):
                toread = os.path.join(aln_dir, file)
                is_U = False
                cmd = "grep -v '>' " + toread + " | grep U"
                files_u = run_command_with_return(cmd)
                if len(files_u) > 0:
                    is_U = True
                    is_u_ever = True
                    if replace:
                        towrite = toread.replace("clean", "clean_noUs")
                        cmd = (
                            "sed '/^>/! {s/U/"
                            + with_what
                            + "/g}' "
                            + toread
                            + " > "
                            + towrite
                        )
                        run_command(cmd)
                    if is_U:
                        print(file + " Contains Us")
        if not is_u_ever:
            print("No files contain Us.")


# # parser.add_argument("-s", "--species_tree", action="store", default="None", help="Species tree, must have internal names and unique names")
# parser.add_argument("-g", "--gene_trees", required=True, action="store", help="Best trees file from phylomeDB")
# # change g param to a directory if you'll use empress
# parser.add_argument("-o", "--outdir", action="store", help="Dir where the resulting files will be stored", default="./generax_data")
# parser.add_argument("-a", "--alignment", action="store", help="alignemnt dir where the alignments are stored")
# parser.add_argument("-km", "--keepmodel", action="store_true", help="Use the same model as in best_trees")
# parser.add_argument("-m", "--mapping", action="store", help="Name of mapping file prefix", default="mapping")
# parser.add_argument("-f", "--family", action="store", help="name of family file", default="family.txt")


def get_astral_pro_data(gene_trees, out_file):

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
    set_sp = set()
    with open(fileTree) as t:
        for line in t:
            line = line.split()
            tree = ete3.Tree(line[3])
            for leaf in tree.iter_leaves():
                leaf.name = leaf.name.split("_")[1]
                set_sp.add(leaf.name)
    return list(set_sp)


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


def get_taxonomic_df(tax_dict, sptree):
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
    for node in sptree.iter_leaves():
        if node.name in mnemo_sp.keys():
            node.species = mnemo_sp[node.name]
        else:
            node.species = node.name

    # tax_resolved_ids = [str(k) for k in tax_resolved]
    whole_tax_dict = {}

    for key, value in tax_resolved.items():
        sp_name = "".join(ncbi.get_taxid_translator([key]).values())
        lineage = ncbi.get_lineage(key)
        names = ncbi.get_taxid_translator(lineage)
        rank = ncbi.get_rank(lineage)
        ordered_names = [names[taxid] for taxid in lineage]
        ordered_clades = [rank[taxid] for taxid in lineage]
        taxonomy = list(zip(ordered_names, ordered_clades))
        # d = {k:{rank[k]:lineage[k]} for k in lineage.keys()}
        whole_tax_dict[value] = taxonomy

    seen = []

    for sp in whole_tax_dict:
        set_single = set()
        for id in whole_tax_dict[sp]:
            set_single.add(id[1])
            #     for clade in whole_tax_dict[sp][id]:
            #         set_single.add(clade)
            seen.append(set_single)

    common_clades = list(set.intersection(*seen) - set(["no rank", "clade", "species"]))
    num_clades = len(common_clades)
    all_phyla = set()

    for key in whole_tax_dict:
        new_tuple = [el for el in whole_tax_dict[key] if el[1] in common_clades][::-1]
        whole_tax_dict[key] = new_tuple
        for el in new_tuple:
            all_phyla.add(el[0])

    color_dict = {}
    colors = ete3.random_color(num=len(set(all_phyla)))
    color_dict = {el[1]: el[0] for el in list(zip(colors, list(all_phyla)))}

    for key in whole_tax_dict:
        new_tuple = [el + (color_dict[el[0]],) for el in whole_tax_dict[key]]
        whole_tax_dict[key] = new_tuple

    # get order and order keys
    for node in sptree.iter_leaves():
        if node.name in tax_dict.keys():
            sp_df = whole_tax_dict[node.species]
        else:
            sp_df = [("", "", "")] * num_clades
        node.add_feature("col_df", sp_df)

    return sptree
    # return whole_tax_dict


def layout_species(node):
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


def layout_species_circular(node):
    node.img_style["size"] = 0
    if node.is_leaf():
        name_face = ete3.faces.AttrFace("species")
        ete3.faces.add_face_to_node(name_face, node, column=0, position="branch-right")


def viz_species_tree(sptree, show=True, render=None):
    ts = ete3.TreeStyle()
    ts.show_leaf_name = False
    ts.allow_face_overlap = True
    ts.draw_aligned_faces_as_table = True
    ts.layout_fn = layout_species
    # sptree.render("prova.png",tree_style = ts)
    if show:
        sptree.show(tree_style=ts)
    if render is not None:
        sptree.render(render, tree_style=ts)


def viz_species_tree_circular(sptree, show=True, render=None):
    ts = ete3.TreeStyle()
    ts.show_leaf_name = False
    ts.layout_fn = layout_species_circular
    ts.mode = "c"
    if show:
        sptree.show(tree_style=ts)
    if render is not None:
        sptree.render(render, tree_style=ts)


def annotate_genetree(genetree, taxo_dict):

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
    a = [frac, dec, num, den]
    fracs = [frac, dec, num, den]

    return fracs


def get_ftp_stats(phylome_ids=None, all=False):
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
