# %% codecell
import analyse_phylome as ap
# %% markdown
# Concatenate the alignments of phylome 13.
#
# First of all, we set all the useful variables. In our project we
# %% codecell


readalPath = (
    "/home/giacomo/master-thesis/Second_proj/obtain_phylome_data_scripts/readal"
)


outDir = "out_dir15"
path = "test_data/all_algs"
pathsFile = outDir + "/paths15.txt"
treeFile = "test_data/best_trees15.txt"
trees121File = outDir + "/phy_15_trees121.txt"
tag = "alg_aa"
# %% markdown
#
#
#
# %% codecell
ap.create_folder(outDir)
# help(ap.create_pathfile)

# remember to tar the algs file
ap.create_pathfile(pathsFile, path, tag)
ap.obtain_121_trees(treeFile, trees121File)

sptree = ap.load_species_tree("test_data/rooted_phylome_15_sp_duptree.txt")
spe2age = ap.build_sp2age(sptree, "341454")

# spe2age = ap.get_all_species(treeFile)

help(ap.build_concatenated_alg)
ap.build_concatenated_alg(trees121File, spe2age, pathsFile, outDir, readalPath, prop = 0.99, max_length=10000, partition=True, treeFile=treeFile)

help(ap.build_extra_concatenated_alg2)
ap.build_extra_concatenated_alg2(treeFile, pathsFile, outDir, readalPath, spe2age, min=17, max_length=10000, partition=True)

help(ap.build_extra_concatenated_alg3)
ap.build_extra_concatenated_alg3(treeFile, pathsFile, outDir, readalPath, spe2age, min=15 , at_least=200, max_length=10000, partition=True)
