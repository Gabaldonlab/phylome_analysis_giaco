import pandas as pd
import os
from ete3 import Tree
import plotly.express as px
import plotly.figure_factory as ff
from plotly.subplots import make_subplots
import glob


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
    for i in ["S","SL","D","T","TL","L"]:
        new_nm = i + "_norm"
        main_df[new_nm] = main_df[i] / main_df['Leaf']
    # main_df.to_csv(args.outfile, index=False)
    return main_df

def get_ecce_df(treeFile, ecce_out):#, out_dir):

    outfile_nm = os.path.basename(ecce_out)
    outfile_nm_noex = os.path.splitext(outfile_nm)[0]

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

    vals = [m+n for m,n in zip(list_lens,vals)]

    df = pd.DataFrame(vals, columns=["Index","Leaf","Leaves_int","D", "T", "L"])
    for i in ["D","T","L"]:
        new_nm = i + "_norm"
        df[new_nm] = df[i] / df['Leaves_int']
    # df.reset_index(inplace=True)
    # df.to_csv(output, index=False)
    # if mult_factor > 1:
    #     print("Warning: genes were repeated " + str(mult_factor) + " times!")
    return df


def ternary_ecce_plot(df_ecce):
    fig = px.scatter_ternary(df_ecce, a="D_norm", b="T_norm", c="L_norm", color="Leaf", hover_name="Index", template="simple_white")#, title="Ternary plots of each gene family ecceTera DTL normalized rates colored by number of leaf in the family")
    fig.show()

def ternary_grax_plot(df_grax, axes=["D_norm","T_norm","SL_norm"]):
    fig = px.scatter_ternary(df_grax, a=axes[0], b=axes[1], c=axes[2], color="Leaf", hover_name="Index", template="simple_white")#, title="Ternary plots of each gene family ecceTera DTL normalized rates colored by number of leaf in the family")
    fig.show()
