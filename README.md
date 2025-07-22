# README

This package provides the function `run_lrip(**kwargs)`, which enables you to directly run the LRIP scoring function. All the modules in the package are presented as APIs. Users can customize scoring functions themselves or link the LRIP to their programs.

## Prerequisite

To run the `lrip_py` package, please install **AMBER molecular dynamics** with version **>= 24**. Then make sure you can execute the following from the command line:

- `antechamber`
- `pmemd`
- `sander`
- `parmchk2`
- `pmemd.MPI`
- `pmemd.cuda`

**Schrodinger Maestro** should also be installed before using this package for molecular docking. `SCHRODINGER` and `AMBERHOME` system variables **should be** set properly.

## How to install

1. Download the code
```bash
git clone git@github.com:nnnty/ipsf_py.git
cd ipsf_py
conda create -n <replace with your environment name> python=3.12
pip install .
```

2. After finished successfully, try to import the module lrip_py in another directory which is different from the ipsf_py/ directory:
    ```python
    import lrip_py
    ```
  This import behavior should not return any thing if lrip_py was successfully installed.

## How to use

Python script for running a job:
```python
import lrip_py
lrip_py.run_lrip('<path of a ligands list file，str>',
                 '<path of a configuration JSON file，str>',
                 '<path of the directory contains the ligands structure file，str>',
                 '<Glide grid file (.zip)，str>',
                 '<job type key word，str>',
                 '<path of the experimental values file，str>',
                 log_file="<log file，str>",
                 over_write=<If the log file will be overwirtten, bool>,
                 )
```

### Format of ligands list file

This file **must** contain ligand names as the first column. We highly recommend providing the net charge for each ligand. For example:
```
<name1>,<net charge>
```

### Suggested configuration JSON file
```json
{
  "job_root": "<job name>",
  "lig_format": "<ligand structure files format, .mol2 format is highly recommended>",
  "receptor_path": "<receptor structure file, must be PDB format>",
  "num_mpi": "<MPI option for the jobs. Could be int or str. If an int is given, the MD will be run through pmemd.MPI. If a cuda str is given, the MD will be performed by pmemd.cuda.>"
}
```

### Job type keywords

1. `train`: Train 7 different ML models using experimental data provided.
2. `pred`: Load the pre-trained ML model to make predictions. For this type, the configuration JSON **must** include a `model_root_path` keyword pointing to the directory storing the pre-trained model.
3. `bs`: Perform bootstrap train-validation 10 times (default). This will write several `evaluation_ML_<i>.csv` files to `<execution directory>/<job name>/BS/`. Each file contains statistical metrics for each ML model. The model order is described in `MODEL_LIST.csv` in the same directory.

### Format of experimental values file

This file (provided as the first argument to `lrip_py.run_lrip()`) **must** contain two key columns: `ID` and `exp`, separated by commas. Example:
```
ID,<other columns>,exp,<more columns>
<name1>,<...>,<value1>,<...>
<name2>,<...>,<value2>,<...>
```

## Important info before running the job

1. Each ligand's structural file **must** have the same name as listed in the ligands list file. For example, if the list includes `lig1`, then a structural file named `lig1.mol2` must be present.
2. Receptor PDB files **must** end with a TER line.。Although the package provides simple handling of PDB files, it is **highly recommended** to manually preprocess them.
3. The ligand name in the mol2 file (i.e., the first line under @<TRIPOS>MOLECULE) **must** exactly match the name in the ligands list.
4. **Ligands containing iodine atoms are not supported. Please remove such ligands in advance.**
