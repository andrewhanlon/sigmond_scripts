#!/usr/bin/env python

import numpy as np
import h5py
import argparse

def main():

  parser = argparse.ArgumentParser(description="Determine different momentum")
  parser.add_argument("--to-merge", type=str, nargs="+", required=True,
                      help="Specify the files to merge")
  parser.add_argument("--output", type=str, required=True,
                      help="Specify the output file")

  args = parser.parse_args()

  to_merge_handlers = list()
  for filename in args.to_merge:
    to_merge_handlers.append(h5py.File(filename, 'r'))

  output_handler = h5py.File(args.output, 'w')

  output_handler.attrs.create("spatExt", to_merge_handlers[0].attrs['spatExt'])
  configs = list()
  for to_merge_handler in to_merge_handlers:
    configs.extend(to_merge_handler.attrs['cfgList'])

  output_handler.attrs.create("cfgList", configs)

  for channel in to_merge_handlers[0].keys():
    output_group = output_handler.create_group(channel)
    output_group.attrs.create('opList', to_merge_handlers[0][channel].attrs['opList'])

    data_to_merge = list()
    for to_merge_handler in to_merge_handlers:
      data_to_merge.append(to_merge_handler[channel]['t0avg'])

    merged_data = np.concatenate(data_to_merge, axis=0)
    output_group.create_dataset('t0avg', data=merged_data)

  output_handler.close()

  for to_merge_handler in to_merge_handlers:
    to_merge_handler.close()

if __name__ == "__main__":
  main()
