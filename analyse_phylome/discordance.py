# -*- coding: utf-8 -*-
import pandas as pd
import os
from ete3 import Tree, NCBITaxa
import plotly.express as px
import plotly.figure_factory as ff

# import plotly.figure_factory as ff
# from plotly.subplots import make_subplots
import glob
from .utils import run_command, get_taxonomic_df, get_dataframe

from .ReconciledTreeIO import recPhyloXML_parser


def get_generax_df(res_dir):
    """Get a Pandas dataframe from Generax output directory.

    Parameters
    ----------
    res_dir : str
        Generax output directory.

    Returns
    -------
    df
        A pandas dataframe with Generax results.

    """
    pattern = res_dir + "/reconciliations/*eventCounts*"
    list_files = glob.glob(pattern)

    main_df = pd.DataFrame()

    for file in list_files:
        id = os.path.basename(file).replace("_eventCounts.txt", "")
        df = pd.read_csv(file, delimiter=":", header=None).transpose()
        df = df.rename(columns=df.iloc[0]).drop(df.index[0])
        cols = df.columns
        df[cols] = df[cols].apply(pd.to_numeric, errors="coerce")
        df["Index"] = id
        main_df = main_df.append(df, ignore_index=True)
    # main_df["Leaf"] = pd.to_numeric(main_df["Leaf"])
    for i in ["S", "SL", "D", "T", "TL", "L"]:
        new_nm = i + "_norm"
        main_df[new_nm] = main_df[i] / main_df["Leaf"]
    # main_df.to_csv(args.outfile, index=False)
    return main_df


def get_ecce_df(treeFile, ecce_out):  # , out_dir):
    """Get a Pandas dataframe from ecceTERA output directory.

    Parameters
    ----------
    treeFile : str
        Path to best trees file.
    ecce_out : str
        Path to eccetera results dir.

    Returns
    -------
    df
        A pandas dataframe with ecceTERA results.

    """

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


def ternary_ecce_plot(
    df_ecce,
    renderer=None,
    template="plotly_white",
    show=True,
    out_img=None,
    interactive=True,
):
    fig = px.scatter_ternary(
        df_ecce,
        a="D_norm",
        b="T_norm",
        c="L_norm",
        color="Leaf",
        hover_name="Index",
        template=template,
    )  # , title="Ternary plots of each gene family ecceTera DTL normalized rates colored by number of leaf in the family")
    if show:
        if renderer is not None:
            fig.show(renderer=renderer)
        else:
            fig.show()
    if not show and out_img is not None:
        if interactive:
            fig.write_html(out_img)
        else:
            fig.write_image(out_img)


def dist_ecce_plot(
    df,
    cols=["D_norm", "T_norm", "L_norm"],
    show_hist=False,
    show_rug=False,
    show_curve=True,
    bin_size=0.03,
    template="plotly_white",
    renderer=None,
    show=True,
    out_img=None,
    interactive=True,
):
    # option for clums and for rug and hist
    fig = ff.create_distplot(
        [df[c] for c in cols],
        cols,
        show_hist=show_hist,
        show_rug=show_rug,
        show_curve=show_curve,
        bin_size=bin_size,
    )

    fig.update_layout(template=template)
    if show:
        if renderer is not None:
            fig.show(renderer=renderer)
        else:
            fig.show()
    if not show and out_img is not None:
        if interactive:
            fig.write_html(out_img)
        else:
            fig.write_image(out_img)


def ternary_grax_plot(
    df_grax,
    axes=["D_norm", "T_norm", "SL_norm"],
    renderer=None,
    template="plotly_white",
    show=True,
    out_img=None,
    interactive=True,
):
    fig = px.scatter_ternary(
        df_grax,
        a=axes[0],
        b=axes[1],
        c=axes[2],
        color="Leaf",
        hover_name="Index",
        template=template,
    )  # , title="Ternary plots of each gene family ecceTera DTL normalized rates colored by number of leaf in the family")
    if show:
        if renderer is not None:
            fig.show(renderer=renderer)
        else:
            fig.show()
    if not show and out_img is not None:
        if interactive:
            fig.write_html(out_img)
        else:
            fig.write_image(out_img)


def dist_grax_plot(
    df,
    cols=["D_norm", "S_norm", "SL_norm", "T_norm"],
    show_hist=False,
    show_rug=False,
    show_curve=True,
    bin_size=0.03,
    template="plotly_white",
    renderer=None,
    show=True,
    out_img=None,
    interactive=True,
):
    # option for clums and for rug and hist
    fig = ff.create_distplot(
        [df[c] for c in cols],
        cols,
        show_hist=show_hist,
        show_rug=show_rug,
        show_curve=show_curve,
        bin_size=bin_size,
    )

    fig.update_layout(template=template)
    if show:
        if renderer is not None:
            fig.show(renderer=renderer)
        else:
            fig.show()
    if not show and out_img is not None:
        if interactive:
            fig.write_html(out_img)
        else:
            fig.write_image(out_img)


# SCATTERPLOT MISSING VS COVERAGE


def read_missingcov(fracmiss, spcov):
    miss_df = pd.read_csv(fracmiss, header=None, sep=" ", names=["SPECIES", "MISSING"])
    cov_df = pd.read_csv(spcov, sep=":")
    merge_df = pd.merge(miss_df, cov_df)
    return merge_df


def plot_missingcov(
    df,
    text=True,
    template="simple_white",
    renderer=None,
    show=True,
    out_img=None,
    interactive=True,
):
    if not text:
        fig = px.scatter(
            df,
            x=" FAMILY_COVERAGE",
            y="MISSING",
            hover_data=["SPECIES"],
            template=template,
        )
        fig.update_traces(textposition="top center")
    else:
        fig = px.scatter(
            df, x=" FAMILY_COVERAGE", y="MISSING", text="SPECIES", template=template
        )
        fig.update_traces(textposition="top center")
    if show:
        if renderer is not None:
            fig.show(renderer=renderer)
        else:
            fig.show()
    if not show and out_img is not None:
        if interactive:
            fig.write_html(out_img)
        else:
            fig.write_image(out_img)


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


def get_ternary_phylome_data(out_df, tree, method):

    phy_id = int("".join(filter(str.isdigit, tree)))
    # for each phylome in files in directory
    # file must cointain either grax or ecce and phylome_id!
    if method not in ["grax", "ecce"]:
        raise ValueError("Method must be either grax or ecce")

    df = pd.read_csv(out_df)
    tree = Tree(tree, format=1)
    num_sp = len(tree.get_leaves())
    num_genes = len(df)
    # Sum of each ctegory will be normalized by num genes * num of species
    norm_factor = num_genes * num_sp

    D_norm = sum(df["D"]) / norm_factor
    T_norm = sum(df["T"]) / norm_factor
    if method == "grax":
        L_norm = sum(df["SL"]) / norm_factor
    elif method == "ecce":
        L_norm = sum(df["L"]) / norm_factor

    norm_to_one = 1 / sum([D_norm, T_norm, L_norm])
    D_norm_one = D_norm * norm_to_one
    T_norm_one = T_norm * norm_to_one
    L_norm_one = L_norm * norm_to_one

    row = [
        phy_id,
        D_norm,
        T_norm,
        L_norm,
        D_norm_one,
        T_norm_one,
        L_norm_one,
        method,
        num_sp,
        num_genes,
    ]
    #     data.append(row)
    # colnames = [
    #     "phy_id",
    #     "D_norm",
    #     "T_norm",
    #     "L_norm",
    #     "D_norm_one",
    #     "T_norm_one",
    #     "L_norm_one",
    #     "method",
    #     "num_sp",
    #     "num_genes",
    # ]
    #
    # df = pd.DataFrame(data, columns=colnames)
    # df["phylome"] = df["phy_id"].astype(str)
    return row


def plot_ternary_phylome(
    df,
    renderer=None,
    template="plotly_white",
    show=True,
    out_img=None,
    interactive=True,
):
    fig = px.scatter_ternary(
        df,
        a="D_norm_one",
        b="T_norm_one",
        c="L_norm_one",
        hover_name="phy_id",
        symbol="method",
        # think of a way to color it if many phylomes
        # color="phylome",
        template=template,
    )
    if show:
        if renderer is not None:
            fig.show(renderer=renderer)
        else:
            fig.show()
    if not show and out_img is not None:
        if interactive:
            fig.write_html(out_img)
        else:
            fig.write_image(out_img)


def get_translation_dict(ecce_sp, grax_sp):
    ecce_names = []
    for node in ecce_sp.traverse():
        ecce_names.append(node.name)
    # prune unsampled leaf as generax does not have it
    ecce_sp.prune([el for el in ecce_names if el != "UNSAMPLED"])

    ecce_dict = {}
    for node in ecce_sp.traverse():
        desc = [el.name for el in node.get_descendants() if el.is_leaf()]
        ecce_dict[node.name] = desc

    grax_dict = {}
    for node in grax_sp.traverse():
        desc = [el.name for el in node.get_descendants() if el.is_leaf()]
        grax_dict[node.name] = desc
    # find keys with same descendants so that internal nodes may be translated
    matches = {}
    for k in ecce_dict:
        for i in grax_dict:
            if len(ecce_dict[k]) == 0:
                matches[k] = k
            elif sorted(ecce_dict[k]) == sorted(grax_dict[i]):
                matches[k] = i
    return matches


def compare_rectree(ecce_rec, grax_rec, matches):
    ecce_sum = ecce_rec.getEventsSummary()
    grax_sum = grax_rec.getEventsSummary()

    new_ecce_sum = {}
    for key, val in ecce_sum.items():
        new_tuple = [matches[el] for el in val if el != "UNSAMPLED"]
        new_ecce_sum[key] = new_tuple

    only_ecce = {}
    only_grax = {}
    both = {}

    for key in new_ecce_sum:
        both[key] = set(new_ecce_sum[key]).intersection(set(grax_sum[key]))
        only_ecce[key] = set(new_ecce_sum[key]).difference(set(grax_sum[key]))
        only_grax[key] = set(grax_sum[key]).difference(set(new_ecce_sum[key]))

    return only_ecce, only_grax, both


def get_ecce_counts(out_dir, pattern="canonical_symmetric"):
    rec_ecce_files = glob.glob(out_dir + "/*" + pattern + "*")
    parser = recPhyloXML_parser()

    # dict_long = {}
    dup_dict = {}
    loss_dict = {}
    trans_dict = {}

    for file in rec_ecce_files:
        ecce_rec = parser.parse(file)
        events_dict = ecce_rec.getEventsSummary()
        for key in events_dict:
            if key == "loss":
                loss = {
                    value: events_dict[key].count(value) for value in events_dict[key]
                }
                loss_dict = {
                    k: loss.get(k, 0) + loss_dict.get(k, 0)
                    for k in set(loss) | set(loss_dict)
                }
            if key == "duplication":
                dup = {
                    value: events_dict[key].count(value) for value in events_dict[key]
                }
                dup_dict = {
                    k: dup.get(k, 0) + dup_dict.get(k, 0)
                    for k in set(dup) | set(dup_dict)
                }
            if key == "transferReception":
                trans = {
                    value: events_dict[key].count(value) for value in events_dict[key]
                }
                trans_dict = {
                    k: trans.get(k, 0) + trans_dict.get(k, 0)
                    for k in set(trans) | set(trans_dict)
                }

    ecce_names = [node.name for node in ecce_rec.spTree.traverse()]

    sp_ev_counts_ecce = []

    for key in ecce_names:
        D = dup_dict.get(key)
        T = trans_dict.get(key)
        L = loss_dict.get(key)
        if D is None:
            D = 0
        if T is None:
            T = 0
        if L is None:
            L = 0
        row = [key, D, T, L]
        sp_ev_counts_ecce.append(row)

    counts_df = pd.DataFrame(sp_ev_counts_ecce, columns=["node", "D", "T", "L"])
    counts_df = counts_df.set_index("node")
    return counts_df
