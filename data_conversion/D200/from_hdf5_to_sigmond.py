#!/usr/bin/env python

import os

import numpy as np
import h5py

import sigmond

data_type_dirs = ['isodoublet_strange_nucleonlambda', 'isosinglet_nonstrange_nucleonnucleon', 'isoquartet_strange_nucleonsigma', 'isodoublet_nonstrange', 'isoquartet_nonstrange_fermionic', 'isosinglet_strange_fermionic', 'single_hadrons', 'isotriplet_nonstrange_nucleonnucleon', 'isosinglet_doublystrange_lambdalambda', 'isoquintet_doublystrange_sigmasigma']

ensemble_name = "cls21_d200_r000"
omissions = {2000}

INPUT_DIR = "/disk1/research/projects/D200_project/data/hdf5/cls21_d200_r000"
OUTPUT_DIR = "/disk1/research/projects/D200_project/data/sigmond/cls21_d200_r000"


def main():

  ensemble_info = sigmond.MCEnsembleInfo(ensemble_name, 'ensembles.xml')
  bins_info = sigmond.MCBinsInfo(ensemble_info)
  bins_info.addOmissions(omissions)
      
  for data_type_dir in data_type_dirs:
    print(data_type_dir)
    input_data_dir = os.path.join(INPUT_DIR, data_type_dir)
    output_data_dir = os.path.join(OUTPUT_DIR, data_type_dir)
    os.makedirs(output_data_dir, exist_ok=True)
    for h5_file in os.listdir(input_data_dir):
      h5_full_file = os.path.join(input_data_dir, h5_file)
      h5_h = h5py.File(h5_full_file, 'r')

      bin_file = os.path.splitext(h5_file)[0] + '.bin'
      bin_full_file = os.path.join(output_data_dir, bin_file)
      bins_handler = sigmond.BinsPutHandler(bins_info, bin_full_file, sigmond.WriteMode.Protect, False)

      for channel, channel_group in h5_h.items():
        channel_data = channel_group['data']

        op_list = channel_group.attrs['op_list']
        for src_i, src_opstr in enumerate(op_list):
          src_op = sigmond.OperatorInfo(src_opstr, sigmond.OpKind.BasicLapH)
          for snk_i, snk_opstr in enumerate(op_list[src_i:], start=src_i):
            snk_op = sigmond.OperatorInfo(snk_opstr, sigmond.OpKind.BasicLapH)
            corr_t = sigmond.CorrelatorAtTimeInfo(snk_op, src_op, 0, True, False)
            for t in range(channel_data.shape[-1]):
              corr_t.resetTimeSeparation(t)

              if len(channel_data.shape) == 2:
                real_data = channel_data[:,t].real
                imag_data = channel_data[:,t].imag

              else:
                real_data = channel_data[:,snk_i,src_i,t].real
                imag_data = channel_data[:,snk_i,src_i,t].imag

              if np.any(real_data):
                obs_info = sigmond.MCObsInfo(corr_t, sigmond.ComplexArg.RealPart)
                bins_handler.putData(obs_info, sigmond.RVector(list(real_data)))

              if np.any(imag_data):
                obs_info = sigmond.MCObsInfo(corr_t, sigmond.ComplexArg.ImaginaryPart)
                bins_handler.putData(obs_info, sigmond.RVector(list(imag_data)))

      bins_handler.close()


if __name__ == "__main__":
  main()
