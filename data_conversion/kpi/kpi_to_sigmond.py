#!/usr/bin/env python

import os
import h5py
import numpy as np

import defs

import sigmondbind as sig

FORCE_HERM = True

def main():

  for ensemble in defs.ensembles:
    print(f"{ensemble}")
    ensemble_info = defs.data_info[ensemble]
    for isospin in defs.isospsins:
      print(f"Converting {isospin}")
      convert_data(ensemble_info["ensemble_name"], ensemble_info["combined_name"], isospin)

    for single_hadron in defs.single_hadrons:
      print(f"Converting {single_hadron}")
      convert_single_hadron(ensemble_info["ensemble_name"], ensemble_info["combined_name"], single_hadron)

def convert_single_hadron(replicum, ensemble_name, single_hadron):
  file_handlers = list()
  for replica in replicum:
    filename = os.path.join(defs.output_dir, f"{replica}_single_hadrons.hdf5")
    file_handlers.append(h5py.File(filename, 'r'))

  ensemble_info = sig.MCEnsembleInfo(ensemble_name, 'ensembles.xml')
  bins_info = sig.MCBinsInfo(ensemble_info)
  bins_info.addOmissions(defs.omissions[ensemble_name])
  sampling_info = sig.MCSamplingInfo()
  xml_obs = sig.XMLHandler("MCObservables", "")
  obs_get_handler = sig.MCObsGetHandler(xml_obs, bins_info, sampling_info)
  obs_handler = sig.MCObsHandler(obs_get_handler, False)

  bin_dir = os.path.join(defs.output_dir, ensemble_name, "single_hadrons")

  op_list = list(file_handlers[0][single_hadron].attrs['opList'])
  tmax = file_handlers[0][single_hadron]['t0avg'].shape[1]
  for file_handler in file_handlers:
    if list(file_handler[single_hadron].attrs['opList']) != op_list:
      print("Op lists differ between replicum")
      exit()
    if file_handler[single_hadron]['t0avg'].shape[1] != tmax:
      print("tmaxes differ between replicum")
      exit()

  os.makedirs(bin_dir, exist_ok=True)
  data_filename = f"{ensemble_name}_{single_hadron}.bin"
  bin_file = os.path.join(bin_dir, data_filename)
  obs_keys = set()

  for op_i, psq_str in enumerate(op_list):
    for tsep in range(tmax):
      bin_datas = list()
      for file_handler in file_handlers:
        bin_data = file_handler[single_hadron]['t0avg'][:,tsep,op_i]
        bin_datas.append(bin_data)

      if not np.any(bin_datas[0]):
        continue

      opstr = f"{defs.single_hadron_isospin[single_hadron]} S={defs.single_hadron_strangeness[single_hadron]}"
      psq_str = psq_str.replace('pSq', 'PSQ=')
      psq = int(psq_str[4:])
      opstr += f" {psq_str} {defs.single_hadron_irreps[single_hadron][psq]} {single_hadron}"

      op = sig.OperatorInfo(opstr, sig.OpKind.GenIrrep)
      corr = sig.CorrelatorAtTimeInfo(op, op, tsep, FORCE_HERM, False)

      corr_re = sig.MCObsInfo(corr, sig.ComplexArg.RealPart)
      re_bin_data = sig.RVector(sum([[x.real for x in bins] for bins in bin_datas], []))
      obs_handler.putBins(corr_re, re_bin_data)
      obs_keys.add(corr_re)

      if not FORCE_HERM:
        corr_im = sig.MCObsInfo(corr, sig.ComplexArg.ImaginaryPart)
        im_bin_data = sig.RVector(sum([[x.imag for x in bins] for bins in bin_datas], []))
        obs_handler.putBins(corr_im, im_bin_data)
        obs_keys.add(corr_im)

    xml_out = sig.XMLHandler("output", "")
    obs_handler.writeBinsToFile(obs_keys, bin_file, xml_out, sig.WriteMode.Protect)
    obs_handler.clearData()

  for file_handler in file_handlers:
    file_handler.close()

def convert_data(replicum, ensemble_name, isospin):
  file_handlers = list()
  for replica in replicum:
    filename = os.path.join(defs.output_dir, f"{replica}_{isospin}.hdf5")
    file_handlers.append(h5py.File(filename, 'r'))

  ensemble_info = sig.MCEnsembleInfo(ensemble_name, 'ensembles.xml')
  bins_info = sig.MCBinsInfo(ensemble_info)
  bins_info.addOmissions(defs.omissions[ensemble_name])
  sampling_info = sig.MCSamplingInfo()
  xml_obs = sig.XMLHandler("MCObservables", "")
  obs_get_handler = sig.MCObsGetHandler(xml_obs, bins_info, sampling_info)
  obs_handler = sig.MCObsHandler(obs_get_handler, False)

  bin_dir = os.path.join(defs.output_dir, ensemble_name, "kpi")
  for channel in file_handlers[0].keys():
    print(f"Working on {channel}")
    psq, irrep = channel.split('-')
    psq = psq.replace('pSq', 'PSQ=')
    op_list = list(file_handlers[0][channel].attrs['opList'])
    tmax = file_handlers[0][channel]['t0avg'].shape[1]
    for file_handler in file_handlers:
      if list(file_handler[channel].attrs['opList']) != op_list:
        print("Op lists differ between replicum")
        exit()
      if file_handler[channel]['t0avg'].shape[1] != tmax:
        print("tmaxes differ between replicum")
        exit()

    data_dir = os.path.join(bin_dir, f"{isospin}", f"{psq.replace('=','')}", irrep)
    os.makedirs(data_dir, exist_ok=True)
    data_filename = f"{ensemble_name}_kpi_{isospin}_{psq.replace('=','')}_{irrep}.bin"
    bin_file = os.path.join(data_dir, data_filename)
    obs_keys = set()

    for opsrc_i, opsrc_str in enumerate(op_list):
      opsnks_enum = enumerate(op_list[opsrc_i:], start=opsrc_i) if FORCE_HERM else enumerate(op_list)
      for opsnk_i, opsnk_str in opsnks_enum:
        for tsep in range(tmax):
          bin_datas = list()
          for file_handler in file_handlers:
            if FORCE_HERM:
              bin_data = 0.5*(file_handler[channel]['t0avg'][:,tsep,opsnk_i,opsrc_i] + np.conj(file_handler[channel]['t0avg'][:,tsep,opsrc_i,opsnk_i]))
            else:
              bin_data = file_handler[channel]['t0avg'][:,tsep,opsnk_i,opsrc_i]

            bin_datas.append(bin_data)

          if not np.any(bin_datas[0]):
            continue

          opsrc = sig.OperatorInfo(f"{isospin} S=1 {psq} {irrep} {opsrc_str.replace(' ', '-')}", sig.OpKind.GenIrrep)
          opsnk = sig.OperatorInfo(f"{isospin} S=1 {psq} {irrep} {opsnk_str.replace(' ', '-')}", sig.OpKind.GenIrrep)
          corr = sig.CorrelatorAtTimeInfo(opsnk, opsrc, tsep, FORCE_HERM, False)

          corr_re = sig.MCObsInfo(corr, sig.ComplexArg.RealPart)
          re_bin_data = sig.RVector(sum([[x.real for x in bins] for bins in bin_datas], []))
          obs_handler.putBins(corr_re, re_bin_data)
          obs_keys.add(corr_re)

          if not FORCE_HERM or opsrc_i != opsnk_i:
            corr_im = sig.MCObsInfo(corr, sig.ComplexArg.ImaginaryPart)
            im_bin_data = sig.RVector(sum([[x.imag for x in bins] for bins in bin_datas], []))
            obs_handler.putBins(corr_im, im_bin_data)
            obs_keys.add(corr_im)

    xml_out = sig.XMLHandler("output", "")
    obs_handler.writeBinsToFile(obs_keys, bin_file, xml_out, sig.WriteMode.Protect)
    obs_handler.clearData()
    print("done")

  for file_handler in file_handlers:
    file_handler.close()




if __name__ == "__main__":
  main()
