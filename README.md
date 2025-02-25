# autoddp



## Getting Started


## Usage

- Create a directory for your project and place autoddp in the root:
```bash
$ wget 
```
- See the available flags:
```bash
$ python3 autoddp.py --help
```
```bash
usage: autoddp.py [-h] --receptor RECEPTOR --config CONFIG [--pH PH] [--ligands LIGANDS] [--complexes COMPLEXES]

autoddp: automate the boring stuff

options:
  -h, --help            show this help message and exit
  --receptor RECEPTOR, -r RECEPTOR
                        Name of receptor file in PDBQT format. Target protein for docking.
  --config CONFIG, -c CONFIG
                        Name of configuration file. Contains parameters for Vina.
  --pH PH, -p PH        pH for ligand protonation. Defaults to 7.4. Can be disabled by specifying pH to 0
  --ligands LIGANDS, -l LIGANDS
                        Name of ligand library in SDF format previously prepared.
  --complexes COMPLEXES, -mc COMPLEXES
                        Make complexes of the top X hits. Defaults to None.
```


- For most virtual screening setups, you can call as follows:
```bash
$ python3 main.py --receptor your_protein.pdb --
```

## Roadmap


## License
This project is licensed under the MIT License.
