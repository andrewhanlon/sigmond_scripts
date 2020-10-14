#!/usr/bin/env python

import os

import h5py
import numpy as np
import tqdm

import sigmond as sig
import defs

FORCE_HERM = True

base_data_dir = "/media/ext2/research/data/raw_dibaryon/"

#ensembles_to_do = ['U102', 'U103', 'B450', 'H101', 'H200', 'A653', 'N300', 'N202', 'E1', 'E5']
ensembles_to_do = ['B450', 'N300']

def main():
  for ensemble in defs.ensembles:
    if ensemble.name not in ensembles_to_do:
      continue

    print(f"Processing ensemble {ensemble.name}")
    for channel, ops in tqdm.tqdm(defs.dibaryon_ops.items()):
      try:
        convert_dibaryons(ensemble, channel, ops)
      except Exception as e:
        print("Exception:")
        print("\tEnsemble: {}".format(ensemble.name))
        print("\tChannel: {}".format(channel))
        print("\tMessage: {}".format(e))
        exit()

    convert_baryons(ensemble)
    convert_pseudoscalar(ensemble)
    if ensemble.name in defs.decuplet_ensembles:
      convert_decuplets(ensemble)


def convert_dibaryons(ensemble, channel, ops):
  if channel.flavor is None and ensemble.su3:
    return

  if channel.flavor is not None and not ensemble.su3:
    return

  if not ops:
    return

  data_dir = os.path.join(base_data_dir, f"analysis_{ensemble.dir_name}")
  data_str = 'c2_dibaryon'
  if ensemble.su3:
    data_str += f"/{channel.P}/{channel.irrep}/{channel.flavor}"
  else:
    data_str += f"/{channel.strangeness}/{channel.isospin}/{channel.P}/{channel.irrep}"


  data_handlers = list()
  raw_datas = list()
  for replica in ensemble.replica:
    if ensemble.su3:
      data_filename = f"{ensemble.dir_name}{replica}_{ensemble.modes}modes.hdf5"
    else:
      data_filename = f"{ensemble.dir_name}{replica}_{ensemble.modes}modes_{channel.strangeness}.hdf5"
    data_file = os.path.join(data_dir, data_filename)
    f_handler = h5py.File(data_file, 'r')
    data_handlers.append(f_handler)
    raw_datas.append(f_handler[data_str])

  ensemble_info = sig.MCEnsembleInfo(f"cls_{ensemble.name}", 'ensembles.xml')
  bins_info = sig.MCBinsInfo(ensemble_info)
  sampling_info = sig.MCSamplingInfo()
  xml_obs = sig.XMLHandler("MCObservables", "")
  obs_get_handler = sig.MCObsGetHandler(xml_obs, bins_info, sampling_info)
  obs_handler = sig.MCObsHandler(obs_get_handler, False)

  if ensemble.su3:
    bin_dir = os.path.join("data", ensemble.name, "dibaryons", channel.P, channel.irrep,
                           channel.flavor)
    os.makedirs(bin_dir, exist_ok=True)
    bin_filename = "{}{}_{}modes_{}_{}_F{}.bin".format(
        ensemble.name, ensemble.replica_str, ensemble.modes, channel.P, channel.irrep,
        channel.flavor)
  else:
    bin_dir = os.path.join("data", ensemble.name, "dibaryons", f"{channel.isospin}_{channel.strangeness}",
                           channel.P, channel.irrep)
    os.makedirs(bin_dir, exist_ok=True)
    bin_filename = "{}{}_{}modes_{}_{}_{}_{}.bin".format(
        ensemble.name, ensemble.replica_str, ensemble.modes, channel.isospin,
        channel.strangeness, channel.P, channel.irrep)

  bin_file = os.path.join(bin_dir, bin_filename)
  obs_keys = set()

  for opsrc_i, opsrc_str in enumerate(ops):
    opsnks_enum = enumerate(ops[opsrc_i:], start=opsrc_i) if FORCE_HERM else enumerate(ops)
    for opsnk_i, opsnk_str in opsnks_enum:
      for t in range(ensemble.t0, ensemble.t0+ensemble.ts):
        bin_datas = list()
        for raw_data in raw_datas:
          bin_data = np.zeros(raw_data.shape[0], dtype=np.complex128)
          for tsrc in range(ensemble.srcs):
            for par in range(2):
              if skip_source(ensemble.name, tsrc, par):
                continue

              if FORCE_HERM:
                bin_data += 0.5*(raw_data[:,tsrc,par,opsnk_i,opsrc_i,t-ensemble.t0] + np.conj(raw_data[:,tsrc,par,opsrc_i,opsnk_i,t-ensemble.t0]))
              else:
                bin_data += raw_data[:,tsrc,par,opsnk_i,opsrc_i,t-ensemble.t0]

          bin_datas.append(bin_data)

        opsrc = sig.OperatorInfo(opsrc_str, sig.OpKind.GenIrrep)
        opsnk = sig.OperatorInfo(opsnk_str, sig.OpKind.GenIrrep)
        corr = sig.CorrelatorAtTimeInfo(opsnk, opsrc, t, FORCE_HERM, False)

        corr_re = sig.MCObsInfo(corr, sig.ComplexArg.RealPart)
        re_bin_data = sig.RVector(sum([[x.real for x in bins] for bins in bin_datas], []))
        obs_handler.putBins(corr_re, re_bin_data)
        obs_keys.add(corr_re)

        if not FORCE_HERM or opsrc != opsnk:
          corr_im = sig.MCObsInfo(corr, sig.ComplexArg.ImaginaryPart)
          im_bin_data = sig.RVector(sum([[x.imag for x in bins] for bins in bin_datas], []))
          obs_handler.putBins(corr_im, im_bin_data)
          obs_keys.add(corr_im)

  xml_out = sig.XMLHandler("output", "")
  obs_handler.writeBinsToFile(obs_keys, bin_file, xml_out, sig.WriteMode.Protect)
  obs_handler.clearData()

  for data_handler in data_handlers:
    data_handler.close()


def convert_decuplets(ensemble):
  data_dir = os.path.join(base_data_dir, f"analysis_{ensemble.dir_name}")
  data_handlers = list()
  raw_datas = list()
  for replica in ensemble.replica:
    data_filename = "{}{}_{}modes_decuplet.hdf5".format(ensemble.dir_name, replica, ensemble.modes)
    data_file = os.path.join(data_dir, data_filename)
    f_handler = h5py.File(data_file, 'r')
    data_handlers.append(f_handler)
    raw_datas.append(f_handler['c2_decuplet'])

  flavors = ['decuplet']
  num_flavors = 1
  if not ensemble.su3:
    flavors = list(map(lambda x: x.decode('utf-8'), f_handler['decuplets'][:]))
    num_flavors = len(flavors)

  ensemble_info = sig.MCEnsembleInfo(f"cls_{ensemble.name}", 'ensembles.xml')
  bins_info = sig.MCBinsInfo(ensemble_info)
  sampling_info = sig.MCSamplingInfo()
  xml_obs = sig.XMLHandler("MCObservables", "")
  obs_get_handler = sig.MCObsGetHandler(xml_obs, bins_info, sampling_info)
  obs_handler = sig.MCObsHandler(obs_get_handler, False)

  bin_dir = os.path.join("data", ensemble.name, "single_particles")
  os.makedirs(bin_dir, exist_ok=True)
  if ensemble.su3:
    bin_filename = "{}{}_{}modes_decuplet.bin".format(ensemble.name, ensemble.replica_str, ensemble.modes)
  else:
    bin_filename = "{}{}_{}modes_decuplets.bin".format(ensemble.name, ensemble.replica_str, ensemble.modes)
  bin_file = os.path.join(bin_dir, bin_filename)
  obs_keys = set()

  for psq in range(1):
    irrep = defs.spin_three_half_irreps[psq]
    for flavor in range(num_flavors):
      flav_name = flavors[flavor]
      isospin = defs.baryon_isospin[flav_name]
      strangeness = defs.baryon_strangeness[flav_name]
      op_str = f"{isospin} S={strangeness} PSQ={psq} {irrep} {flav_name} 0"
      op = sig.OperatorInfo(op_str, sig.OpKind.GenIrrep)

      for t in range(ensemble.t0, ensemble.t0+ensemble.ts):
        bin_datas = list()
        for raw_data in raw_datas:
          bin_data = np.zeros(raw_data.shape[0], dtype=np.complex128)
          for tsrc in range(ensemble.srcs):
            for par in range(2):
              if skip_source(ensemble.name, tsrc, par):
                continue

              bin_data += raw_data[:,tsrc,par,psq,flavor,t-ensemble.t0]
          bin_datas.append(bin_data)

        corr = sig.CorrelatorAtTimeInfo(op, op, t, FORCE_HERM, False)

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

  for data_handler in data_handlers:
    data_handler.close()


def convert_baryons(ensemble):
  data_dir = os.path.join(base_data_dir, f"analysis_{ensemble.dir_name}")
  data_handlers = list()
  raw_datas = list()
  for replica in ensemble.replica:
    data_filename = "{}{}_{}modes_baryon.hdf5".format(ensemble.dir_name, replica, ensemble.modes)
    data_file = os.path.join(data_dir, data_filename)
    f_handler = h5py.File(data_file, 'r')
    data_handlers.append(f_handler)
    raw_datas.append(f_handler['c2_baryon'])

  flavors = ['lambda']
  num_flavors = 1
  if not ensemble.su3:
    flavors = list(map(lambda x: x.decode('utf-8'), f_handler['baryons'][:]))
    num_flavors = len(flavors)

  ensemble_info = sig.MCEnsembleInfo(f"cls_{ensemble.name}", 'ensembles.xml')
  bins_info = sig.MCBinsInfo(ensemble_info)
  sampling_info = sig.MCSamplingInfo()
  xml_obs = sig.XMLHandler("MCObservables", "")
  obs_get_handler = sig.MCObsGetHandler(xml_obs, bins_info, sampling_info)
  obs_handler = sig.MCObsHandler(obs_get_handler, False)

  bin_dir = os.path.join("data", ensemble.name, "single_particles")
  os.makedirs(bin_dir, exist_ok=True)
  if ensemble.su3:
    bin_filename = "{}{}_{}modes_lambda.bin".format(ensemble.name, ensemble.replica_str, ensemble.modes)
  else:
    bin_filename = "{}{}_{}modes_baryons.bin".format(ensemble.name, ensemble.replica_str, ensemble.modes)
  bin_file = os.path.join(bin_dir, bin_filename)
  obs_keys = set()

  for psq in range(4):
    irrep = defs.spin_half_irreps[psq]
    for flavor in range(num_flavors):
      flav_name = flavors[flavor]
      isospin = defs.baryon_isospin[flav_name]
      strangeness = defs.baryon_strangeness[flav_name]
      op_str = f"{isospin} S={strangeness} PSQ={psq} {irrep} {flav_name} 0"
      op = sig.OperatorInfo(op_str, sig.OpKind.GenIrrep)

      for t in range(ensemble.t0, ensemble.t0+ensemble.ts):
        bin_datas = list()
        for raw_data in raw_datas:
          bin_data = np.zeros(raw_data.shape[0], dtype=np.complex128)
          for tsrc in range(ensemble.srcs):
            for par in range(2):
              if skip_source(ensemble.name, tsrc, par):
                continue

              bin_data += raw_data[:,tsrc,par,psq,flavor,t-ensemble.t0]
          bin_datas.append(bin_data)

        corr = sig.CorrelatorAtTimeInfo(op, op, t, FORCE_HERM, False)

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

  for data_handler in data_handlers:
    data_handler.close()


def convert_pseudoscalar(ensemble):
  data_dir = os.path.join(base_data_dir, f"analysis_{ensemble.dir_name}")
  data_handlers = list()
  tsrc_list = dict()
  raw_datas = dict()
  for ps_name in defs.pseudoscalar_names[ensemble.name]:
    raw_datas[ps_name] = list()

  for rep_num, replica in enumerate(ensemble.replica):
    data_filename = "{}{}_{}modes_pseudoscalar.hdf5".format(ensemble.dir_name, replica, defs.pseudoscalar_modes[ensemble.name])
    data_file = os.path.join(data_dir, data_filename)
    f_handler = h5py.File(data_file, 'r')
    data_handlers.append(f_handler)
    for ps_name in defs.pseudoscalar_names[ensemble.name]:
      raw_datas[ps_name].append(f_handler[ps_name])

    if ensemble.name in defs.tsrc_files:
      tsrc_list[rep_num] = dict()
      tsrc_file = os.path.join(data_dir, defs.tsrc_files[ensemble.name] + replica)
      with open(tsrc_file, 'r') as tsrc_file_h:
        for line in tsrc_file_h:
          conf, tsrc = line.split()
          tsrc_list[rep_num][int(conf)-1] = int(tsrc)

  ensemble_info = sig.MCEnsembleInfo(f"cls_{ensemble.name}", 'ensembles.xml')
  bins_info = sig.MCBinsInfo(ensemble_info)
  sampling_info = sig.MCSamplingInfo()
  xml_obs = sig.XMLHandler("MCObservables", "")
  obs_get_handler = sig.MCObsGetHandler(xml_obs, bins_info, sampling_info)
  obs_handler = sig.MCObsHandler(obs_get_handler, False)

  bin_dir = os.path.join("data", ensemble.name, "single_particles")
  os.makedirs(bin_dir, exist_ok=True)
  bin_filename = "{}{}_{}modes_pseudoscalar.bin".format(ensemble.name, ensemble.replica_str, defs.pseudoscalar_modes[ensemble.name])
  bin_file = os.path.join(bin_dir, bin_filename)
  obs_keys = set()
  for ps_name in defs.pseudoscalar_names[ensemble.name]:
    for t in range(ensemble.Nt):
      bin_datas = list()
      for rep_num, raw_data in enumerate(raw_datas[ps_name]):
        bin_data = np.zeros(raw_data.shape[0], dtype=np.complex128)
        for tsrc in range(ensemble.srcs):

          if tsrc_list:
            for config in range(raw_data.shape[0]):
              t0 = defs.pseudoscalar_sources[ensemble.name][tsrc] + tsrc_list[rep_num][config]
              ts_forward = (t0+t) % ensemble.Nt
              ts_reverse = (t0-t) % ensemble.Nt

              if not ensemble.open:
                bin_data[config] += raw_data[config,tsrc,ts_forward]
              elif use_forward_pseudoscalar(ensemble.name, tsrc) and use_backward_pseudoscalar(ensemble.name, tsrc):
                bin_data[config] += raw_data[config,tsrc,ts_forward]
                bin_data[config] += raw_data[config,tsrc,ts_reverse]
              elif use_forward_pseudoscalar(ensemble.name, tsrc):
                bin_data[config] += raw_data[config,tsrc,ts_forward]
              elif use_backward_pseudoscalar(ensemble.name, tsrc):
                bin_data[config] += raw_data[config,tsrc,ts_reverse]

          else:
            t0 = defs.pseudoscalar_sources[ensemble.name][tsrc]
            ts_forward = (t0+t) % ensemble.Nt
            ts_reverse = (t0-t) % ensemble.Nt

            if not ensemble.open:
              bin_data += raw_data[:,tsrc,ts_forward]
            elif use_forward_pseudoscalar(ensemble.name, tsrc) and use_backward_pseudoscalar(ensemble.name, tsrc):
              bin_data += raw_data[:,tsrc,ts_forward]
              bin_data += raw_data[:,tsrc,ts_reverse]
            elif use_forward_pseudoscalar(ensemble.name, tsrc):
              bin_data += raw_data[:,tsrc,ts_forward]
            elif use_backward_pseudoscalar(ensemble.name, tsrc):
              bin_data += raw_data[:,tsrc,ts_reverse]

        bin_datas.append(bin_data)

      op_str = defs.pseudoscalar_op_strs[ensemble.name][ps_name]
      op = sig.OperatorInfo(op_str, sig.OpKind.GenIrrep)
      corr = sig.CorrelatorAtTimeInfo(op, op, t, FORCE_HERM, False)

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

  for data_handler in data_handlers:
    data_handler.close()

def skip_source(ensemble_name, tsrc, par):
  if par == 0 and tsrc in defs.forward_prop_skip[ensemble_name]:
    return True
  elif par == 1 and tsrc in defs.backward_prop_skip[ensemble_name]:
    return True
  return False

def use_forward_pseudoscalar(ensemble_name, tsrc):
  if tsrc not in defs.pseudoscalar_forward_prop_skip[ensemble_name]:
    return True
  return False

def use_backward_pseudoscalar(ensemble_name, tsrc):
  if tsrc not in defs.pseudoscalar_backward_prop_skip[ensemble_name]:
    return True
  return False

if __name__ == "__main__":
  main()
