import analyse_phylome as ap


outDir = "out_dir15"
treeFile = "test_data/best_trees15.txt"
aln_dir = "test_data/all_algs"
sptree = "test_data/rooted_phylome_15_sp_duptree.txt"



# generax
# help(ap.get_generax_data)
out_sprax = outDir + "/sprax15"

# scanforus
ap.get_generax_data(treeFile, out_sprax, aln_dir)

a = get_generax_df(res_dir)

ternary_grax_plot(a)

# then mpiexec -np $np generax -s /gpfs/projects/bsc40/gmutti/Second_proj/all_phylomes/rooted_phylome_$i\_sp_duptree.txt --seed 42 --families /gpfs/scratch/bsc40/bsc40122/all_generax/families/$i\_generax_family.txt --strategy EVAL --rec-model UndatedDTL --per-family-rates --prefix /gpfs/scratch/bsc40/bsc40122/all_generax/generax_phylome_$i

# then analyse



# ecceTERA
ecce_out = "test_data/ecce_15.out"
ap.get_ecce_data(treeFile, sptree)

# ecceTERA species.file=/gpfs/projects/bsc40/gmutti/Second_proj/all_phylomes/rooted_phylome_$i\_sp_duptree.txt gene.file=/gpfs/scratch/bsc40/bsc40122/all_ecce/ecce_best_trees$i.txt dated=0 verbose=1 print.newick=1 output.dir=/gpfs/scratch/bsc40/bsc40122/all_ecce/ecce_phylome_$i/ print.reconciliations=1 print.info=1 orthology.output=1 print.newick=1 print.graph=1 > /gpfs/scratch/bsc40/bsc40122/all_ecce/ecce_$i.out

ecce_df = ap.get_ecce_df(treeFile, ecce_out)

ap.ternary_ecce_plot(ecce_df)


# tomorrow 26 ago finish doc 
