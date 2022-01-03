#!/usr/bin/env python

import os

import numpy as np
import h5py

data_type_dirs = ['isodoublet_strange_nucleonlambda', 'isosinglet_nonstrange_nucleonnucleon', 'isoquartet_strange_nucleonsigma', 'isodoublet_nonstrange', 'isoquartet_nonstrange_fermionic', 'isosinglet_strange_fermionic', 'single_hadrons', 'isotriplet_nonstrange_nucleonnucleon', 'isosinglet_doublystrange_lambdalambda', 'isoquintet_doublystrange_sigmasigma']

ensemble_name = "cls21_d200_r000"
omissions = {2000}

INPUT_DIR = "/disk1/research/projects/D200_project/data/hdf5/cls21_d200_r000"


def main():

  missing_channels = list()
  for data_type_dir in data_type_dirs:
    has_missing = False
    input_data_dir = os.path.join(INPUT_DIR, data_type_dir)
    for h5_file in os.listdir(input_data_dir):
      h5_full_file = os.path.join(input_data_dir, h5_file)
      h5_h = h5py.File(h5_full_file, 'r')

      for channel, channel_group in h5_h.items():
        channel_data = channel_group['data']

        op_list = channel_group.attrs['op_list']
        occurences = channel_group.attrs['op_occurences']

        for op_str, occurence in zip(op_list, occurences):
          if occurence <= 1:
            print(f'    - "{op_str}"')
            has_missing = True

    print()
    if has_missing:
      missing_channels.append(data_type_dir)

  print()
  print()
  print("missing")
  for missing_channel in missing_channels:
    print(f"  {missing_channel}")
    


if __name__ == "__main__":
  main()
