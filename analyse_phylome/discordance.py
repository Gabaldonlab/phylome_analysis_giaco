import pandas as pd
import os
from ete3 import Tree, NCBITaxa
import plotly.express as px

# import plotly.figure_factory as ff
# from plotly.subplots import make_subplots
import glob
from .utils import run_command, get_taxonomic_df, get_dataframe


def get_generax_df(res_dir):
    pattern = res_dir + "/reconciliations/*eventCounts*"
    list_files = glob.glob(pattern)

    main_df = pd.DataFrame()

    for file in list_files:
        id = os.path.basename(file).replace("_eventCounts.txt", "")
        df = pd.read_csv(file, delimiter=":", header=None).transpose()
        df = df.rename(columns=df.iloc[0]).drop(df.index[0])
        df["Index"] = id
        main_df = main_df.append(df, ignore_index=True)
    main_df["Leaf"] = pd.to_numeric(main_df["Leaf"])
    for i in ["S", "SL", "D", "T", "TL", "L"]:
        new_nm = i + "_norm"
        main_df[new_nm] = main_df[i] / main_df["Leaf"]
    # main_df.to_csv(args.outfile, index=False)
    return main_df


def get_ecce_df(treeFile, ecce_out):  # , out_dir):

    # outfile_nm = os.path.basename(ecce_out)
    # outfile_nm_noex = os.path.splitext(outfile_nm)[0]

    # output = out_dir + "/" + outfile_nm_noex + "_results.csv"

    list_lens = []
    with open(treeFile) as t:
        for line in t:
            line = line.split()
            seed = line[0]
            single_t = Tree(line[3])
            num_sp = len(single_t.get_leaves())
            num_all = len([n for n in single_t.traverse()])
            list_lens.append([seed, num_sp, num_all])

    vals = []
    with open(ecce_out) as o:
        for line in o:
            if "nAm" in line:
                line = line.split(",")
                dups = int(line[1].replace("duplications=", "").strip())
                transf = int(line[2].replace("transfers=", "").strip())
                losses = int(line[3].replace("losses=", "").strip())
                vals.append([dups, transf, losses])

    # # why????????? THIS IS AN ISSUE
    # mult_factor = int(len(vals) / len(list_lens))
    # list_lens = list_lens * mult_factor

    vals = [m + n for m, n in zip(list_lens, vals)]

    df = pd.DataFrame(vals, columns=["Index", "Leaf", "Leaves_int", "D", "T", "L"])
    for i in ["D", "T", "L"]:
        new_nm = i + "_norm"
        df[new_nm] = df[i] / df["Leaves_int"]
    # df.reset_index(inplace=True)
    # df.to_csv(output, index=False)
    # if mult_factor > 1:
    #     print("Warning: genes were repeated " + str(mult_factor) + " times!")
    return df


def ternary_ecce_plot(df_ecce, renderer="notebook_connected"):
    fig = px.scatter_ternary(
        df_ecce,
        a="D_norm",
        b="T_norm",
        c="L_norm",
        color="Leaf",
        hover_name="Index",
        template="simple_white",
    )  # , title="Ternary plots of each gene family ecceTera DTL normalized rates colored by number of leaf in the family")
    fig.show(renderer=renderer)


def ternary_grax_plot(
    df_grax, axes=["D_norm", "T_norm", "SL_norm"], renderer="notebook_connected"
):
    fig = px.scatter_ternary(
        df_grax,
        a=axes[0],
        b=axes[1],
        c=axes[2],
        color="Leaf",
        hover_name="Index",
        template="simple_white",
    )  # , title="Ternary plots of each gene family ecceTera DTL normalized rates colored by number of leaf in the family")
    fig.show(renderer=renderer)


def get_proka_file(whole_tax_dict, out_dir="./"):
    proka_list = []
    for key in whole_tax_dict:
        for el in whole_tax_dict[key]:
            if el[0] == "Bacteria":
                proka_list.append(key)

    if len(proka_list) > 0:
        proka_file = out_dir + "proka_mnemo.txt"

        with open(proka_file, "w") as f:
            for item in proka_list:
                f.write("%s\n" % item)
    else:
        print("No prokaryotic entry found!")


def get_abaccus_taxo(tax_dict, out_file, set_cols=True):
    whole_tax_dict = get_taxonomic_df(tax_dict, set_cols=set_cols, fill=True)
    pdf = get_dataframe(whole_tax_dict)
    ncbi = NCBITaxa()
    tax_resolved = ncbi.get_taxid_translator([sp for sp in tax_dict.values()])
    tax_dict = {k: v for k, v in tax_dict.items() if int(v) in tax_resolved.keys()}
    mnemo_sp = {}
    for mnemo in tax_dict:
        id = int(tax_dict[mnemo])
        try:
            species = tax_resolved[id]
        except KeyError:
            species = mnemo
        mnemo_sp[species] = mnemo
    new_mnemo = {}
    for el in pdf["species"]:
        for key in mnemo_sp:
            if key.startswith(el):
                new_mnemo[el] = mnemo_sp[key]

    mnemo_col = pdf["species"].map(new_mnemo)
    pdf.insert(0, "mnemo", mnemo_col)
    tax_col = pdf["mnemo"].map(tax_dict)
    pdf.insert(0, "taxaid", tax_col)
    pdf.to_csv(out_file, index=False, sep=";")
    print(
        "Taxonomy file created. Remember to check if all is okay or eventually add euka entries manually."
    )


def run_abaccus(
    abaccus_path,
    best_trees,
    out_dir,
    taxonomy,
    jumps=2,
    losses=3,
    collapse_dup=True,
    proka=None,
    greasyfile=None,
):
    cmds = []
    with open(best_trees) as bt:
        for line in bt:
            line = line.split("\t")
            main = line[0]
            tree = line[3].strip()
            cmd = (
                "python "
                + abaccus_path
                + ' -i "'
                + tree
                + '"'
                + " -t "
                + taxonomy
                + " -m "
                + main
                + " -o "
                + out_dir
                + "/abaccus_output.txt"
                + " -J "
                + str(jumps)
                + " -L "
                + str(losses)
            )
            if collapse_dup:
                cmd = cmd + " -c "  # collapse_dup
            if proka is not None:
                cmd = cmd + " -p " + proka

            cmds.append(cmd)

            if greasyfile is None:
                run_command(cmd)

    if greasyfile is not None:
        out_greasy = out_dir + "/" + greasyfile
        with open(out_greasy, "w") as gf:
            for item in cmds:
                gf.write("%s\n" % item)
