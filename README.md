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
