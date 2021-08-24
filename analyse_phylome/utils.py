import subprocess as sp
import sys
import os
import ete3
import glob
from ftplib import FTP
from Bio import SeqIO, Seq

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
        cmd = "mkdir " + name
        try:
            run_command(cmd, False)
        except:
            print("Unable to create directory " + name)


def loadPaths(pathsFile, tag):
    pathsList = {}
    for line in open(pathsFile):
        line = line.strip()
        dades = line.split("\t")
        if dades[0] == tag:
            pathsList[dades[1]] = dades[2]
    return pathsList


def build_sp2age(sptree, seed):
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
    # if spe2age:
    #     try:
    #         t.set_outgroup(t.get_farthest_oldest_leaf(spe2age))
    #     except:
    #         sys.exit("Something went wrong with sp2age dict!")
    # elif midpoint:
    #     t.set_outgroup(t.get_midpoint_outgroup())
    # else:
    #     sys.exit("Something went wrong with rooting!")
    return t


def root_tree(t, spe2age=None, midpoint=False):
    if spe2age:
        try:
            t.set_outgroup(t.get_farthest_oldest_leaf(spe2age))
        except:
            sys.exit("Something went wrong with sp2age dict!")
    elif midpoint:
        t.set_outgroup(t.get_midpoint_outgroup())
    else:
        sys.exit("Something went wrong with rooting!")
    return t


def print_sequence(code, sequence, outfile):
    outfile.write(">" + code + "\n")
    i = 0
    if sequence[-1] == "*":
        sequence = sequence[:-1]
    while i < len(sequence):
        outfile.write(sequence[i : i + 60] + "\n")
        i += 60


def create_pathfile(pathsFile, algs):
    outfilePath = open(pathsFile, "w")
    for firstDir in glob.glob(algs + "/*"):
        code = firstDir.split("/")[-1].split(".")[0]
        outfilePath.write("alg_aa\t" + code + "\t" + firstDir + "\n")
    outfilePath.close()


# get speciesrax data

# set alignment args as optional way better than this!
# mapping and family file not imporant


def get_generax_data(gene_trees, out_dir, aln_dir=None, keep_model=True):

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
                    f.write("alignment = " + aln_file + "\n")
                f.write("mapping = " + out_mapfile + "\n")


def scan_for_Us(aln_dir):#, replace=False):
    if os.path.exists(aln_dir):
        for file in os.listdir(aln_dir):
            if file.endswith("clean.fasta"):
                toread = os.path.join(aln_dir, file)
                for record in SeqIO.parse(toread, 'fasta'):
                    # if replace:
                    if 'U' in record.seq:
                    #         record.seq = str(record.seq).replace('U','X')
                    #     else:
                    #         pass
                    #
                        # else:
                        print(file  + " contains Us, Raxml won't parse this alignment, replace them!")


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


def get_species_name(node_name_string):
    spcode = node_name_string.split("_")[1]
    return spcode


def get_ecce_data(
    gene_trees, species_tree, out_dir, root=False, ultrametric=False, midpoint=False
):

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


def get_tax_dict(tree, uniprot_df):

    sp = [nm for nm in tree.get_leaf_names() if not nm.isdecimal()]
    sp = [s.split("_")[-1] if "_" in s else s for s in sp]

    with open(uniprot_df) as f:
        sp_line = [line.split() for line in f.readlines()]
        sp_line = [line for line in sp_line if line != []]
        taxo_dict = {
            line[0]: line[2].replace(":", "") for line in sp_line if line[0] in sp
        }

    absent = [abs for abs in sp if abs not in taxo_dict.keys()]
    if absent:
        print("Warning: " + " ".join(absent) + " was not found in the dictionary.")
        print("You may want to add it manually")

    return taxo_dict


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
