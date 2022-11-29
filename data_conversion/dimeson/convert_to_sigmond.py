#!/usr/bin/env python

import os
import sys

import h5py
import numpy as np
import regex
import tqdm

import sigmond as sig
import defs

FORCE_HERM = True
FORCE_REAL = True

base_data_dir = "/disk3/research/data/dimeson_correlators/"
output_dir = "data"

ensembles_to_do = ["B450", "N300", "N202", "H101"]

def main():

  for ensemble in defs.ensembles:
    if ensemble.name not in ensembles_to_do:
      continue

    print(f"Processing ensemble {ensemble.name}")
    try:
      convert_dimesons(ensemble)
      convert_mesons(ensemble)
    except Exception as e:
      print("Exception:")
      print("\tEnsemble: {}".format(ensemble.name))
      print("\tMessage: {}".format(e))

def convert_dimesons(ensemble):
  ensemble_info = sig.MCEnsembleInfo(f"{ensemble.type}_{ensemble.name}", 'ensembles.xml')
  bins_info = sig.MCBinsInfo(ensemble_info)
  sampling_info = sig.MCSamplingInfo()
  xml_obs = sig.XMLHandler("MCObservables", "")
  obs_get_handler = sig.MCObsGetHandler(xml_obs, bins_info, sampling_info)
  obs_handler = sig.MCObsHandler(obs_get_handler, False)

  raw_data_dir = os.path.join(base_data_dir, f"{ensemble.dir_name}")

  data_handlers = list()
  for replica in ensemble.replica:
    data_filename = f"{ensemble.dir_name}{replica}_{ensemble.modes}modes_DD.hdf5"
    data_file = os.path.join(raw_data_dir, data_filename)
    data_handler = h5py.File(data_file, 'r')
    data_handlers.append(data_handler)

  primary_data_handler = data_handlers[0]

  for isospin in primary_data_handler['/c2_DD'].keys():
    for pref in primary_data_handler[f'/c2_DD/{isospin}'].keys():
      for irrep in primary_data_handler[f'/c2_DD/{isospin}/{pref}'].keys():
        dataset_path = f'/c2_DD/{isospin}/{pref}/{irrep}'
        opsets = [list(data_handler[dataset_path].attrs['operators']) for data_handler in data_handlers]
        if not all([opset == opsets[0] for opset in opsets]):
          print(f"not all opsets identical for {dataset_path}")
          sys.exit()

        opset = opsets[0]
        datasets = [data_handler[dataset_path] for data_handler in data_handlers]

        bin_dir = os.path.join(output_dir, ensemble.name, "dimesons", isospin)
        os.makedirs(bin_dir, exist_ok=True)
        bin_filename = f"{ensemble.name}{ensemble.replica_str}_{ensemble.modes}modes_{pref}_{irrep}_{isospin}.bin"
        bin_file = os.path.join(bin_dir, bin_filename)
        obs_keys = set()

        for opsrc_i, opsrc_str in enumerate(opset):
          opsrc_str = defs.get_opstr(pref, irrep, isospin, opsrc_str)
          opsnks_enum = enumerate(opset[opsrc_i:], start=opsrc_i) if FORCE_HERM or FORCE_REAL else enumerate(opset)
          for opsnk_i, opsnk_str in opsnks_enum:
            opsnk_str = defs.get_opstr(pref, irrep, isospin, opsnk_str)
            for t in range(ensemble.t0, ensemble.t0+ensemble.ts):
              bin_datas = list()
              for dataset in datasets:
                bin_data = np.zeros(dataset.shape[0], dtype=np.complex128)
                for tsrc in range(ensemble.srcs):
                  for par in range(2):
                    if skip_source(ensemble.name, tsrc, par):
                      continue

                    if FORCE_REAL:
                      bin_data += 0.5*(dataset[:,tsrc,par,opsnk_i,opsrc_i,t-ensemble.t0].real + dataset[:,tsrc,par,opsrc_i,opsnk_i,t-ensemble.t0].real)
                    elif FORCE_HERM:
                      bin_data += 0.5*(dataset[:,tsrc,par,opsnk_i,opsrc_i,t-ensemble.t0] + np.conj(dataset[:,tsrc,par,opsrc_i,opsnk_i,t-ensemble.t0]))
                    else:
                      bin_data += dataset[:,tsrc,par,opsnk_i,opsrc_i,t-ensemble.t0]

                bin_datas.append(bin_data)

              opsrc = sig.OperatorInfo(opsrc_str, sig.OpKind.GenIrrep)
              opsnk = sig.OperatorInfo(opsnk_str, sig.OpKind.GenIrrep)
              corr = sig.CorrelatorAtTimeInfo(opsnk, opsrc, t, FORCE_HERM or FORCE_REAL, False)

              corr_re = sig.MCObsInfo(corr, sig.ComplexArg.RealPart)
              re_bin_data = sig.RVector(sum([[x.real for x in bins] for bins in bin_datas], []))
              obs_handler.putBins(corr_re, re_bin_data)
              obs_keys.add(corr_re)

              if (not FORCE_REAL and not FORCE_HERM) or opsrc != opsnk:
                corr_im = sig.MCObsInfo(corr, sig.ComplexArg.ImaginaryPart)
                im_bin_data = sig.RVector(sum([[x.imag for x in bins] for bins in bin_datas], []))
                obs_handler.putBins(corr_im, im_bin_data)
                obs_keys.add(corr_im)

        xml_out = sig.XMLHandler("output", "")
        obs_handler.writeBinsToFile(obs_keys, bin_file, xml_out, sig.WriteMode.Protect, 'D')
        obs_handler.clearData()

  for data_handler in data_handlers:
    data_handler.close()


meson_flavors = ['ll', 'lc', 'cc']
meson_isospin = {
    'll': 'I1',
    'lc': 'I1h',
    'cc': 'I1',
}
def convert_mesons(ensemble):
  ensemble_info = sig.MCEnsembleInfo(f"{ensemble.type}_{ensemble.name}", 'ensembles.xml')
  bins_info = sig.MCBinsInfo(ensemble_info)
  sampling_info = sig.MCSamplingInfo()
  xml_obs = sig.XMLHandler("MCObservables", "")
  obs_get_handler = sig.MCObsGetHandler(xml_obs, bins_info, sampling_info)
  obs_handler = sig.MCObsHandler(obs_get_handler, False)

  raw_data_dir = os.path.join(base_data_dir, f"{ensemble.dir_name}")

  data_handlers = list()
  for replica in ensemble.replica:
    data_filename = f"{ensemble.dir_name}{replica}_{ensemble.modes}modes_DD.hdf5"
    data_file = os.path.join(raw_data_dir, data_filename)
    data_handler = h5py.File(data_file, 'r')
    data_handlers.append(data_handler)

  primary_data_handler = data_handlers[0]

  for pref in primary_data_handler[f'/c2_meson'].keys():
    for irrep in primary_data_handler[f'/c2_meson/{pref}'].keys():
      dataset_path = f'/c2_meson/{pref}/{irrep}'
      opsets = [list(data_handler[dataset_path].attrs['operators']) for data_handler in data_handlers]
      if not all([opset == opsets[0] for opset in opsets]):
        print(f"not all opsets identical for {dataset_path}")
        sys.exit()

      opset = opsets[0]
      datasets = [data_handler[dataset_path] for data_handler in data_handlers]

      bin_dir = os.path.join(output_dir, ensemble.name, "single_mesons")
      os.makedirs(bin_dir, exist_ok=True)
      bin_filename = f"{ensemble.name}{ensemble.replica_str}_{ensemble.modes}modes_{pref}_{irrep}_single-mesons.bin"
      bin_file = os.path.join(bin_dir, bin_filename)
      obs_keys = set()

      for meson_flavor_i, meson_flavor in enumerate(meson_flavors):
        for opsrc_i, opsrc_str in enumerate(opset):
          opsrc_str = defs.get_opstr(pref, irrep, meson_isospin[meson_flavor], f"{opsrc_str}_{meson_flavor}")
          opsnks_enum = enumerate(opset[opsrc_i:], start=opsrc_i) if FORCE_HERM or FORCE_REAL else enumerate(opset)
          for opsnk_i, opsnk_str in opsnks_enum:
            opsnk_str = defs.get_opstr(pref, irrep, meson_isospin[meson_flavor], f"{opsnk_str}_{meson_flavor}")
            for t in range(ensemble.t0, ensemble.t0+ensemble.ts):
              bin_datas = list()
              for dataset in datasets:
                bin_data = np.zeros(dataset.shape[0], dtype=np.complex128)
                for tsrc in range(ensemble.srcs):
                  for par in range(2):
                    if skip_source(ensemble.name, tsrc, par):
                      continue

                    if FORCE_REAL:
                      bin_data += 0.5*(dataset[:,tsrc,par,meson_flavor_i,opsnk_i,opsrc_i,t-ensemble.t0].real + dataset[:,tsrc,par,meson_flavor_i,opsrc_i,opsnk_i,t-ensemble.t0].real)
                    elif FORCE_HERM:
                      bin_data += 0.5*(dataset[:,tsrc,par,meson_flavor_i,opsnk_i,opsrc_i,t-ensemble.t0] + np.conj(dataset[:,tsrc,par,meson_flavor_i,opsrc_i,opsnk_i,t-ensemble.t0]))
                    else:
                      bin_data += dataset[:,tsrc,par,meson_flavor_i,opsnk_i,opsrc_i,t-ensemble.t0]

                bin_datas.append(bin_data)

              opsrc = sig.OperatorInfo(opsrc_str, sig.OpKind.GenIrrep)
              opsnk = sig.OperatorInfo(opsnk_str, sig.OpKind.GenIrrep)
              corr = sig.CorrelatorAtTimeInfo(opsnk, opsrc, t, FORCE_HERM or FORCE_REAL, False)

              corr_re = sig.MCObsInfo(corr, sig.ComplexArg.RealPart)
              re_bin_data = sig.RVector(sum([[x.real for x in bins] for bins in bin_datas], []))
              obs_handler.putBins(corr_re, re_bin_data)
              obs_keys.add(corr_re)

              if (not FORCE_REAL and not FORCE_HERM) or opsrc != opsnk:
                corr_im = sig.MCObsInfo(corr, sig.ComplexArg.ImaginaryPart)
                im_bin_data = sig.RVector(sum([[x.imag for x in bins] for bins in bin_datas], []))
                obs_handler.putBins(corr_im, im_bin_data)
                obs_keys.add(corr_im)

      xml_out = sig.XMLHandler("output", "")
      obs_handler.writeBinsToFile(obs_keys, bin_file, xml_out, sig.WriteMode.Protect, 'D')
      obs_handler.clearData()

  for data_handler in data_handlers:
    data_handler.close()


def skip_source(ensemble_name, tsrc, par):
  if par == 0 and tsrc in defs.forward_prop_skip[ensemble_name]:
    return True
  elif par == 1 and tsrc in defs.backward_prop_skip[ensemble_name]:
    return True
  return False


if __name__ == "__main__":
  main()
