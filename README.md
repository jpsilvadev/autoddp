# autoddp

autoddp automates protein-ligand virtual screening workflows

## Getting Started
- The recommended way to use `autoddp` is in a conda environment. If you do not have conda (or alternative)
get [miniconda](https://docs.anaconda.com/miniconda/install/#quick-command-line-install)

- Create a environment and install the dependencies:
```bash
$ conda create -n autoddp -c conda-forge pymol-open-source vina openbabel
$ conda activate autoddp
```
- autoddp uses:
  - `vina` as the molecular docking engine
  - `openbabel` for file conversions and to assign ligand protonation
  - `pymol` to assemble the protein-ligand complexes

## Usage

- Create a directory for your project and place autoddp in the root:
```bash
$ wget https://raw.githubusercontent.com/jpsilvadev/autoddp/refs/heads/main/autoddp.py
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
- `autoddp` expects to find in the root directory:
  - prepared receptor in `pdbqt` format, 
  - prepared ligand library in `sdf` format
  - config file to be passed to `vina`

- If you do not have a project ready or just wish to test `autoddp`, you can use the example provided in `examples/` by running the example.sh script:
```bash
$ git clone https://github.com/jpsilvadev/autoddp.git
$ ./example.sh
```
> Note: This will pin your cpu to 100% for the duration of the molecular docking process.
> If you wish to manually specify, add the `num_cpus` arg to `conf.txt` and run the command below

- For most virtual screening setups, you can follow this structure:

```bash
$ python3 autoddp.py --receptor receptor.pdb --config conf.txt --ligands ligands.sdf --pH 7.4 --complexes 10
```

- autoddp will run the virtual screening, aggregate the results and make the protein-ligand complexes for the top hits, in an organized manner.
```bash
├── autoddp.py
├── backup # individual ligand sdf files
├── complexes # writes the protein-ligand complexes for the number specified in --complexes
├── conf.txt
├── docking_log.txt # logs errors
├── inputs # stores the pdbqt of the ligands passed to vina
├── logs # logs for the vina run of each ligand
├── outputs # written poses
├── receptor.pdbqt
└── results
```
- `results`, and more specifically, `results_sorted.txt` will contain the ranking of the ligands based on the top pose
> With this info and the protein-ligand complexes, you can proceed with 3D visualization/inspection and with 2D interaction plot analysis



## Roadmap
autoddp is WIP. Currently, its a single file script for portability reasons.
The most immediate additions will be:
- Create 2D protein-ligand interactions plots with [prolif](https://prolif.readthedocs.io/en/latest/index.html)
- Refactoring: Replace syscalls  with python bindings for various dependencies and modularize. 
- Create conda, deb and rpm packages
- Use [meeko](https://meeko.readthedocs.io/en/release-doc/) to allow protein and ligand preparation
- Seek out a solution for MD system preparation. Current option is using [pyCHARMM](https://academiccharmm.org/documentation/version/c47b1/pycharmm#Examples)
- Add molecular dynamics functionality with [GROMACS](https://manual.gromacs.org/current/gmxapi/userguide/usage.html)


## License
This project is licensed under the MIT License.
