Analyse_phylome
--------

This python package will allow to:

* Concatenate alignments
* Build and visualize species tree
* Prepare data for discordance softwares (GeneRax and ecceTERA) and analyze the results

Starting from a phylome stored in PhylomeDB.

It is mostly copied from Marina Marcet script in [obtain_phylome_data.BSC.py](https://github.com/Gabaldonlab/projects/blob/master/phylome_scripts/obtain_phylome_data.BSC.py) but it's structured as a python package (python3 compatible) and it has different species tree building methods and a gene-tree/species-tree discordance module.

### Installation

`git clone https://github.com/Gabaldonlab/phylome_analysis_giaco`

Then go in phylome_analysis_giaco dir and run:

`pip install .`

If you are in the cluster. Clone the directory then cd into it and run:

`python setup.py install --user`

## Usage

See jupyter notebooks in "examples" directory to browse different cases. Alternatively most function have documentation accessible through help(function)

#### DOUBTS:

* duptree decide if rooting gene trees (before it was midpoint rooting, now no rooting).

#### TODOs:

* maybe add species and annotate monophyletyc taxonomy in species tree generax viz. Single gene trees have per species count data!!!! Maybe a plotly dashboard will be perfect although overkill!
* eccetera plots. why in stdout number differ?? recphyloxml once updated then use https://github.com/WandrilleD/recPhyloXML ----> add warning that unsampled is not counted + add example with recPhyloXML.
* add instructions to get files to reproduce example and add files that cant be obtained.
* Documentation + errors mgmt + eventually automatic testing
* add this workflow https://docs.github.com/en/actions/guides/building-and-testing-python
* add consensus tree + compare trees function.
* improve box layout consistency + better code for viz

global variables for outdir, readal, duptree etc how do they work??

Eventually

* Merge with other stuff from Marina scripts you did not use
* add angst and other discordance softwares data creation scripts
