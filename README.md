# NUISANCE Neutrino Data Repository

A repository for NUISANCE-maintained HepData records.

## Conventions

This is meant as a list of optional extras and code snippets that might be useful for rapidly converting existing NUISANCE data releases to the standard HepData format.

The [NUISANCE HepData Specification](https://github.com/NUISANCEMC/HEPData), must always be adhered to when converting or creating HepData releases for automatic consumption.

### Directory Structure

Data releases should be stored within the `data` directory in the repository root. We strongly suggest storing new release directories under subdirectories named according to the experiment that the measurement was made with. Further subdirectories are optional, but can help partition experiment-specific directories that contain a lot of measurements. A release directory does not have to be named anything specifically, we recommend naming them something relatively short that corresponds to the type of measurement(s) contained and includes the year of publication. For example, `data/T2K/CH-target/CC0Pi_TKI_2018`, might refer to the release for [this publication](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.032003).

The directory structure can be relatively freeform because release directories are discovered by searching for `ToHepData.py` scripts. The means that you must always name your submission creation script `ToHepData.py`.

### `ToHepData.py`

Every record directory must contain a `ToHepData.py` script. It should take no arguments and be executable in the venv described [below](#hepdata-lib-venv). It should output a `submission.tar.gz` file that can be uploaded to HepData, and produce a directory, `submission-<inspireid>`, containing the unpacked HepData submission files. A template follows:

```python
#!/usr/bin/env python3

from hepdata_lib import Submission, Table, Variable, Uncertainty, RootFileReader

paper_short_ref = "EPJC.83.782"
INSPIRE_id=2638628

# Submission creation goes here

submission.add_additional_resource(description="Python conversion script used to build this submisson. Part of NUISANCE.",
    location="ToHepData.py",
    copy_file=True)

submission.create_files(f"submission-{INSPIRE_id}", remove_old=True)
```

Input files for conversion, such as official data release files, can be downloaded as part of the submission creation process, or can be manually downloaded and stored in relevant submission directory.

More examples of submission creation are given in the main [NUISANCE HEPData specification respository](https://github.com/NUISANCEMC/HEPData#building-a-submission).

### What to commit

Everything that is needed to build the data releases should be committed, this includes manually downloaded official data releases and the ToHepData.py scripts. The submission directory and the submission.tar.gz files themselves should not be comitted and are ignored by default in the `.gitignore` for this repository.

## Utilities

### Updating the local record store

Run `update_local_store` to build a local directory structure of `inspireid`-type reference records that can be used with NUISANCE.

This local record store should be committed to this repository.

### `hepdata_lib` venv

A `requirements.txt` file is included in env/ to help you set up a venv for running the ToHepData.py scripts in this repository. You can initialize the venv from the repository root directory like:

```bash
export NUISANCE_DATA_ROOT=$(readlink -f .)
cd env
python3 -m venv .
. bin/activate 
pip install -r requirements.txt
```

On returning, the `setup.sh` script in the repository root can be sourced to activate the venv.