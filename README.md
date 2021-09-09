Analyse_phylome
--------

This python package will allow to:

* concatenate the alignments
* build and visualize species tree
* Prepare data for discordance softwares (GeneRax and ecceTERA) and analyze the results

Starting from a phylome stored in PhylomeDB:


It is mostly copied from Marina Marcet script in [obtain_phylome_data.BSC.py](https://github.com/Gabaldonlab/projects/blob/master/phylome_scripts/obtain_phylome_data.BSC.py) but it's structured as a python package (python3 compatible) and it has different species tree building methods and a gene-tree/species-tree discordance module.

### Installation



### Concatenation module

##### Method 1

diff w/ Marina:

* User now can set threshold of species with prop

##### Method 2

diff w/ Marina:

* now load_sequence does not root every time.
* User can set both min number of sequence and minimum number or length of concat alns

##### Method 3

diff w/ Marina:

* User can set both min number of sequence and minimum number or length of concat alns


### Species Tree module



### Discordance module



#### DOUBTS:

* duptree decide if rooting gene trees (before it was midpoint rooting, now no rooting).
* wonder if dependencies will be installed if they are in requirements
* how to normalize DTL values?
* add compaison even though already easy in ete?
* add functions to run software or leave it to the user?


#### TODOs:

* Documentation + errors mgmt + eventually automatic testing
* add download phylome script.
* add installation instructions
* misc example
* add this workflow https://docs.github.com/en/actions/guides/building-and-testing-python
* add consensus tree
* viz species tree can be imporved with dynami rectangles width

global variables for outdir, readal, duptree etc how do they work??

Eventually

* Merge with other stuff from Marina scripts you did not use
* add angst and other discordance softwares data creation scripts
* greasifier???
* Merge into Marina pipeline
