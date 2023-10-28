"""
Assumptions made:
 - ligands are in library format (e.g., a single file containing all ligands)
"""
import os
import time
import shutil
import subprocess
import logging
import argparse
from contextlib import suppress
import pymol

PATH = os.getcwd()
LOGS = os.path.join(PATH, "logs")
OUTS = os.path.join(PATH, "outputs")
RES = os.path.join(PATH, "results")
BACKUP = os.path.join(PATH, "backup")
INPUTS = os.path.join(PATH, "inputs")
COMPS = os.path.join(PATH, "complexes")


def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.FileHandler("docking_log.txt"), logging.StreamHandler()],
    )


def move_results(extension, path, position):
    try:
        if not os.path.exists(path):
            os.mkdir(path)
        for file in os.listdir(os.getcwd()):
            if getattr(file, position)(extension) and file != RECEPTOR_FILE:
                with suppress(shutil.Error):
                    shutil.move(file, path)
    except Exception as err:
        logging.error(f"Error moving files: {err}")


def run_subprocess(cmd):
    try:
        completed_process = subprocess.run(
            cmd, shell=True, text=True, capture_output=True, check=True
        )
        return completed_process.stdout
    except subprocess.CalledProcessError as err:
        logging.error(f"Error executing command: {err.stderr}")
        return ""


def convert_sdf_to_pdbqt(lig_library, pH=7.4):
    if pH == 0:
        os.system(f"obabel {lig_library} -O prep_subs.sdf --unique cansmiNS")
    else:
        os.system(f"obabel {lig_library} -O prep_subs.sdf -p {pH} --unique cansmiNS")

    os.system("obabel -isdf prep_subs.sdf -osdf -O *.sdf --split --unique")
    os.system("obabel -isdf *.sdf -opdbqt -O*.pdbqt")

    for file in ["prep_subs.pdbqt", "prep_subs.sdf", "conformers.pdbqt"]:
        with suppress(FileNotFoundError):
            os.remove(file)


def run_vina_on_ligands():
    for file in os.listdir(PATH):
        if file == RECEPTOR_FILE:
            continue
        elif file.endswith(".pdbqt"):
            cmd = f"vina --config {CONFIG_FILE} --ligand {file}"
            output = run_subprocess(cmd)
            print(output)
            with open(f"{file}_log.log", "w", encoding="utf-8") as log:
                log.write(output)


def extract_and_sort_results():
    os.system("tail -n14 *.log > results.txt")

    print("\nStarting analysis...")
    dct = {}

    for file in os.listdir(PATH):
        if file.endswith("out.pdbqt"):
            with open(file, "r", encoding="utf-8") as f:
                try:
                    for i, line in enumerate(f):
                        if i == 1:
                            dct[file] = float(line.split(":")[1].split()[0])
                            break
                except Exception as e:
                    logging.error(f"Error processing {file}. \nError: {e}")

    ordered_dict = {k: v for k, v in sorted(dct.items(), key=lambda item: item[1])}

    with open("results_sorted.txt", "w", encoding="utf-8") as out:
        out.write("Sorted Docking Results\n\n")
        for k, v in ordered_dict.items():
            out.write(f"{k}: {v}\n")

    print("\n")
    print("Analysis complete. See your results in the results_sorted.txt file.")


def assemble_complexes_list(comp_num):
    with open(os.path.join(RES, "results_sorted.txt"), "r", encoding="utf-8") as res:
        lines = res.readlines()[2:]
        complexes = []
        for i in range(comp_num):
            complexes.append(lines[i].split(":")[0].strip())

    os.mkdir(COMPS)
    for comp in complexes:
        src_file = os.path.join(OUTS, comp)
        dest_file = os.path.join(COMPS, comp)
        shutil.copy(src_file, dest_file)
    shutil.copy(RECEPTOR_FILE, COMPS)


def make_complexes():

    os.chdir(COMPS)
    pymol.pymol_argv = ["pymol", "-qc"]
    pymol.finish_launching()

    for lig in os.listdir():
        pymol.cmd.load(RECEPTOR_FILE, "receptor")
        pymol.cmd.load(lig, "lig")
        complex_name = f"complex_{lig}.pdb"
        pymol.cmd.create(complex_name, "receptor or lig")
        pymol.cmd.save(complex_name, complex_name)
        pymol.cmd.reinitialize()

    pymol.cmd.quit()


def main():
    setup_logging()
    convert_sdf_to_pdbqt(LIBRARY_FILE, pH=PH)

    time1 = time.time()

    move_results(".sdf", BACKUP, "endswith")
    run_vina_on_ligands()
    extract_and_sort_results()

    time2 = time.time()
    runtime = time2 - time1
    print("\n")
    print(str(runtime / 60) + " mins runtime.")

    move_results(".log", LOGS, "endswith")
    move_results("_out.pdbqt", OUTS, "endswith")
    move_results(".pdbqt", INPUTS, "endswith")
    move_results("results", RES, "startswith")

    if COMPLEXES is not None:
        assemble_complexes_list(COMPLEXES)
        make_complexes()


if __name__ == "__main__":
    # Set up argparse to handle command-line arguments
    parser = argparse.ArgumentParser(description="Script for Vina Docking")
    parser.add_argument(
        "--receptor",
        "-r",
        required=True,
        type=str,
        help="Name of receptor file in PDBQT format. Target protein for docking.",
    )
    parser.add_argument(
        "--config",
        "-c",
        required=True,
        type=str,
        help="Name of configuration file. Contains parameters for Vina.",
    )
    parser.add_argument(
        "--pH",
        "-p",
        type=float,
        default=7.4,
        help="pH for ligand protonation. Defaults to 7.4. Can be disabled by specifying pH to 0",
    )
    parser.add_argument(
        "--ligands",
        "-l",
        type=str,
        help="Name of ligand library in SDF format previously prepared.",
    )
    parser.add_argument(
        "--complexes",
        "-mc",
        type=int,
        default=None,
        help="Make complexes of the top X hits. Defaults to None.",
    )
    args = parser.parse_args()

    # Bind args values to global variables
    RECEPTOR_FILE = args.receptor
    CONFIG_FILE = args.config
    LIBRARY_FILE = args.ligands
    PH = args.pH
    COMPLEXES = args.complexes

    main()
