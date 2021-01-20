#!/usr/bin/env python

import os
import tarfile

import argparse

replicas = ['r005', 'r006', 'r007', 'r008']
sources = [0, 12, 24, 36, 48, 60, 72, 84]
parities = ['fwd', 'bwd']

def main():

  parser = argparse.ArgumentParser(description="Untar C103 data")
  parser.add_argument("-i", "--input", type=str, required=True, metavar='input directory',
                      help="Specify directory containing raw data")
  parser.add_argument("-o", "--output", type=str, required=True, metavar='output directory',
                      help="Specify output directory to write averaged data")

  args = parser.parse_args()

  for replica in replicas:
    for src_i, source in enumerate(sources):
      for parity in parities:
        tarfile_name = f"merged_{replica}_{parity}_src{src_i}.tar.gz"
        tarfile_path = os.path.join(args.input, tarfile_name)
        print(f"Extracting {tarfile_path}...", end='', flush=True)
        tar_h = tarfile.open(tarfile_path)
        rel_output_dir = os.path.join(replica, f"src{source}", parity)
        full_output_dir = os.path.join(args.output, rel_output_dir)
        os.makedirs(full_output_dir)
        tar_h.extractall(full_output_dir)
        print("done")

if __name__ == "__main__":
  main()
