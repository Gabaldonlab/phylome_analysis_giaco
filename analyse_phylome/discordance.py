import pandas as pd
import os


def analyze_ecce(gene_trees, ecce_out, out_dir):

    outfile_nm = os.path.basename(ecce_out)
    outfile_nm_noex = os.path.splitext(outfile_nm)[0]

    output = out_dir + "/" + outfile_nm_noex + "_results.csv"

    gene = []

    with open(gene_trees) as t:
        for line in t:
            line = line.split()
            gene.append(line[0])

    vals = []

    with open(ecce_out) as o:
        for line in o:
            if "nAm" in line:
                line = line.split(",")
                dups = int(line[1].replace("duplications=", "").strip())
                transf = int(line[2].replace("transfers=", "").strip())
                losses = int(line[3].replace("losses=", "").strip())
                vals.append([dups, transf, losses])

    # why????????? THIS IS AN ISSUE
    mult_factor = int(len(vals) / len(gene))

    gene = gene * mult_factor

    df = pd.DataFrame(vals, index=gene, columns=["D", "T", "L"])
    df.reset_index(inplace=True)
    df.to_csv(output, index=False)
    if mult_factor > 1:
        print("Warning: genes were repeated " + str(mult_factor) + " times!")
