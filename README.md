# Sigmond Scripts #

This repository contains a set of scripts used for the analysis of 2-point temporal correlation functions.
They depend on the sigmond analysis package ([https://bitbucket.org/ahanlon/sigmond](https://bitbucket.org/ahanlon/sigmond)),
and you should first have sigmond installed before using these scripts.

The analysis scripts use YAML configuration files to specify what you want to do.

### Directory structure ###

- analysis - contains the main scripts for running the analysis
    - run_sigmond.py - the driving script (pass -h to see options)
- data_conversion - contains various scripts for converting data to various formats (e.g. hdf5, LapH-binary, sigmond-bins).
- example_yamls - contains some example yaml files that can be passed to run_sigmond.py
- doc - contains documenation for how to use the analysis scripts (not up to date at all)

### What software requirements are there? ###

- Python 3.8 is required
- Many Python modules
    - wheel
    - cython
    - pybind11
    - pyyaml
    - progressbar
    - sortedcontainers
    - pylatex
    - numpy
    - uncertainties
    - aenum
    - tqdm
    - h5py
    - matplotlib

### How do I get started? ###

- You first need data in a format that sigmond can read (either LapH-binary or sigmond-bins).
  The data_conversion scripts can be used for this
- Next you need some yaml configuration files to pass to run_sigmond.py using the -c option.

### YAML files ###

For more information about YAML files, visit [The Official YAML Web Site](https://yaml.org/)

Note that you can pass various yaml files using the -c option in run_sigmond.py.
These will be combined to form one yaml file that actually gets run.
Be careful about identical sections, because the last section encountered is the one that actually gets run.
Using multiple yaml files allows for aliases and anchors between these files.

### C103 Example ###

- You first need the data in a format that sigmond can read.
It can read the raw LapH data files produced from `last_laph`.
However, the keys are identical in the files for each source time, so sigmond cannot distinguish these.
We must first average all sources that have the same keys, then put the result in a separate file.
Either these files will be provided for you, or you can generate them using the script `data_conversion/C103/to_sigmond.py`.
Make sure to update the `base_data_dir` in `data_conversion/C103/defs.py`. `base_data_dir` should contain the following
    - r005/src0
    - r005/src1
    - r005/src2
    - r005/src3

  where the raw correlator files are located in each of the `src*` directories.
Note that these scripts take 4.5 hours on my local computer to run.

- `cd example_yaml/C103_NN/`
- edit `ensemble_info/C103_all.yml`:

        Initialize:
          ensembles_file: <location of ensembles.xml file>  # can be left blank to use default
          project_directory: <location on local machine where you want all project files to be created>
          raw_data_directories:
            - <location of converted data from first step above>
          precompute: true   # doesn't really matter

    You can also change the resampling mode to Bootstrap.
          
          MCSamplingInfo:
            mode: bootstrap
            number_resampling: 1000
            seed: 3103
            boot_skip: 301

- edit `tasks/averaged_data.yml`. Change `sigmond_batch` to the location of the `sigmond` binary produced when you installed sigmond.
- Next run

        run_sigmond.py -d -c ensemble_info/C103_all.yml tasks/averaged_data.yml

  This run will average over all equivalent momentum frames and irrep rows, then produce a PDF showing the averaged data.
Note that reading the large data files can take some time (not more than 5 minutes on my computer)
One this task has been completed, subsequent tasks will be much easier.

- edit `ensemble_info/C103_none.yml` to point to the correct file locations.
edit `tasks/rotations/isosinglet_nonstrange_nucleonnucleon.yml` to change `sigmond_batch` location.

- Then run

        run_sigmond.py -d -c ensemble_info/C103_none.yml other_info/pivot_info.yml tasks/rotations/isosinglet_nonstrange_nucleonnucleon.yml

  This will run the GEVP using the parameters in `other_info/pivot_info.yml` and the operators in `tasks/rotations/isosinglet_nonstrange_nucleonnucleon.yml`

- edit `tasks/spectrum/isosinglet_nonstrange_nucleonnucleon.yml` to change `sigmond_batch` location.

- Finally, run
        
        run_sigmond.py -d -c ensemble_info/C103_none.yml other_info/pivot_info.yml other_info/minimizer_info.yml non_interacting_levels/isosinglet_nonstrange_nucleonnucleon.yml tasks/spectrum/isosinglet_nonstrange_nucleonnucleon.yml

  This will produce the final spectrum. The non-interacting levels used to form the ratio are defined in `non_interacting_levels/isosinglet_nonstrange_nucleonnucleon.yml`,
and the fit choices are defined in `tasks/spectrum/isosinglet_nonstrange_nucleonnucleon.yml`.
