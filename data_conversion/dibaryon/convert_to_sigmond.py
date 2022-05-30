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

base_data_dir = "/disk3/research/data/dibaryon/"
output_dir = "data"

#ensembles_to_do = ["A653", "B450", "B451", "B452", "H101", "H200", "J500", "N200", "N202", "N300", "U102", "U103", "E1", "E5"]
ensembles_to_do = ["a064_m400_mL6.4_trMc"]

particle_map = {
    'Σ': 'S',
    'Ξ': 'X',
    'Λ': 'L',
}

def main():

  dibaryon_ops = read_op_files()

  '''
  import json
  print(json.dumps(dibaryon_ops, indent=2))
  sys.exit()
  '''

  for ensemble in defs.ensembles:
    if ensemble.name not in ensembles_to_do:
      continue

    print(f"Processing ensemble {ensemble.name}")
    try:
      convert_dibaryons(ensemble, dibaryon_ops)
      convert_baryons(ensemble)
      if ensemble.name in defs.decuplet_ensembles:
        convert_decuplet(ensemble)

      convert_pseudoscalar(ensemble)
    except Exception as e:
      print("Exception:")
      print("\tEnsemble: {}".format(ensemble.name))
      print("\tMessage: {}".format(e))

def read_op_files():
  op_data = dict()
  for strangeness, filename in defs.op_files.items():
    with open(filename, 'r') as fh:
      isospin = None
      flavor = None
      pref_str = None
      pref_key = None
      irrep_key = None
      irrep_str = None
      for line in fh:
        line = line.rstrip()
        whitespace = len(line) - len(line.lstrip())
        if whitespace == 0:
          isospin = line
          flavor = f"{isospin}_{strangeness}"
          op_data[flavor] = dict()

        elif whitespace == 2:
          pref, irrep = line.split(')')
          pref_str = pref.replace(' ', '') + ")"
          pref_key = "P" + "".join(pref_str.replace('(', '').replace(')', '').split(','))
          psq = int(pref_key[1])**2 + int(pref_key[2])**2 + int(pref_key[3])**2
          if pref_key not in op_data[flavor]:
            op_data[flavor][pref_key] = dict()

          irrep_key = irrep.strip()
          irrep_str = defs.convert_irrep(irrep_key, psq)
          op_data[flavor][pref_key][irrep_key] = list()

        elif whitespace == 3:
          op_index, op_str = line.split(': ')
          op_index = int(op_index)
          if len(op_data[flavor][pref_key][irrep_key]) != op_index:
            print("op_index problem")
            sys.exit()

          op_str = op_str.replace(' ', '_')
          for un_str, safe_str in particle_map.items():
            op_str = op_str.replace(un_str, safe_str)

          op_data[flavor][pref_key][irrep_key].append(op_str)

  # get irrep-spin ordering
  op_data_ordering = dict()
  for pref, irreps in op_data[defs.irrep_spin_ordering].items():
    op_data_ordering[pref] = dict()
    for irrep, op_strs in irreps.items():
      op_data_ordering[pref][irrep] = list()
      for op_str in op_strs:
        irrep_spin_str = "_".join(op_str.split(' ')[-1].split('_')[1:])
        if irrep_spin_str not in op_data_ordering[pref][irrep]:
          op_data_ordering[pref][irrep].append(irrep_spin_str)

  # get SU(3)-flavor irreps
  for flavor in defs.SU3_flavors:
    if flavor in defs.symmetric_SU3_flavors and flavor in defs.anti_symmetric_SU3_flavors:
      op_data[flavor] = dict()
      for pref, irreps in op_data_ordering.items():
        psq = int(pref[1])**2 + int(pref[2])**2 + int(pref[3])**2
        pref_str = f"Pref=({pref[1]},{pref[2]},{pref[3]})"
        op_data[flavor][pref] = dict()
        for irrep, irrep_spin_strs in irreps.items():
          irrep_str = defs.convert_irrep(irrep, psq)
          op_data[flavor][pref][irrep] = list()
          for irrep_spin_str in irrep_spin_strs:
            symmetric_op_str = f"{defs.symmetric_SU3_flavor_dibaryon}_{irrep_spin_str}"
            if symmetric_op_str in op_data[defs.symmetric_SU2_flavor][pref][irrep]:
              op_data[flavor][pref][irrep].append(symmetric_op_str)

            anti_symmetric_op_str = f"{defs.anti_symmetric_SU3_flavor_dibaryon}_{irrep_spin_str}"
            if anti_symmetric_op_str in op_data[defs.anti_symmetric_SU2_flavor][pref][irrep]:
              op_data[flavor][pref][irrep].append(anti_symmetric_op_str)

    elif flavor in defs.symmetric_SU3_flavors:
      op_data[flavor] = op_data[defs.symmetric_SU2_flavor]
    elif flavor in defs.anti_symmetric_SU3_flavors:
      op_data[flavor] = op_data[defs.anti_symmetric_SU2_flavor]

  return op_data

def convert_dibaryons(ensemble, ops):
  ensemble_info = sig.MCEnsembleInfo(f"{ensemble.type}_{ensemble.name}", 'ensembles.xml')
  bins_info = sig.MCBinsInfo(ensemble_info)
  sampling_info = sig.MCSamplingInfo()
  xml_obs = sig.XMLHandler("MCObservables", "")
  obs_get_handler = sig.MCObsGetHandler(xml_obs, bins_info, sampling_info)
  obs_handler = sig.MCObsHandler(obs_get_handler, False)

  raw_data_dir = os.path.join(base_data_dir, f"{ensemble.dir_name}")
  
  flavors = defs.SU3_flavors if ensemble.su3 else defs.SU2_flavors

  for flavor in flavors:
    data_handlers = list()
    for replica in ensemble.replica:
      if ensemble.type == "exp":
        data_filename = f"{ensemble.dir_name}{replica}_{ensemble.modes}modes_dibaryon.hdf5"
      elif ensemble.su3:
        data_filename = f"{ensemble.dir_name}{replica}_{ensemble.modes}modes.hdf5"
      else:
        isospin, strangeness = flavor.split('_')
        data_filename = f"{ensemble.dir_name}{replica}_{ensemble.modes}modes_{strangeness}.hdf5"

      data_file = os.path.join(raw_data_dir, data_filename)
      if not os.path.isfile(data_file):
        data_filename = f"{ensemble.name}{replica}_SU3_{ensemble.modes}modes.hdf5"
        data_file = os.path.join(raw_data_dir, data_filename)

      data_handler = h5py.File(data_file, 'r')
      data_handlers.append(data_handler)

    for pref, irreps in ops[flavor].items():
      for irrep, op_strs in irreps.items():
        if not op_strs:
          continue

        datasets = list()
        for data_handler in data_handlers:
          if ensemble.su3:
            dataset_name = f"/c2_dibaryon/{pref}/{irrep}/{flavor}"
          else:
            isospin, strangeness = flavor.split('_')
            dataset_name = f"/c2_dibaryon/{strangeness}/{isospin}/{pref}/{irrep}"

          dataset = data_handler[dataset_name]
          op_size = dataset.shape[3]
          if op_size != len(op_strs):
            print(f"Mismatch in operator size: {dataset_name}")
            sys.exit()

          datasets.append(dataset)

        # make bin files
        bin_dir = os.path.join(output_dir, ensemble.name, "dibaryons", flavor)
        os.makedirs(bin_dir, exist_ok=True)
        bin_filename = f"{ensemble.name}{ensemble.replica_str}_{ensemble.modes}modes_{pref}_{irrep}_{flavor}.bin"
        bin_file = os.path.join(bin_dir, bin_filename)
        obs_keys = set()

        for opsrc_i, opsrc_str in enumerate(op_strs):
          opsrc_str = defs.get_opstr(ensemble.su3, pref, irrep, flavor, opsrc_str)
          opsnks_enum = enumerate(op_strs[opsrc_i:], start=opsrc_i) if FORCE_HERM else enumerate(op_strs)
          for opsnk_i, opsnk_str in opsnks_enum:
            opsnk_str = defs.get_opstr(ensemble.su3, pref, irrep, flavor, opsnk_str)
            for t in range(ensemble.t0, ensemble.t0+ensemble.ts):
              bin_datas = list()
              for dataset in datasets:
                bin_data = np.zeros(dataset.shape[0], dtype=np.complex128)
                for tsrc in range(ensemble.srcs):
                  for par in range(2):
                    if skip_source(ensemble.name, tsrc, par):
                      continue

                    if FORCE_HERM:
                      bin_data += 0.5*(dataset[:,tsrc,par,opsnk_i,opsrc_i,t-ensemble.t0] + np.conj(dataset[:,tsrc,par,opsrc_i,opsnk_i,t-ensemble.t0]))
                    else:
                      bin_data += dataset[:,tsrc,par,opsnk_i,opsrc_i,t-ensemble.t0]

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


def convert_baryons(ensemble):
  data_dir = os.path.join(base_data_dir, f"{ensemble.dir_name}")
  data_handlers = list()
  datasets = list()
  for replica in ensemble.replica:
    data_filename = f"{ensemble.dir_name}{replica}_{ensemble.modes}modes_baryon.hdf5"
    data_file = os.path.join(data_dir, data_filename)
    if not os.path.isfile(data_file):
      if ensemble.type == "exp":
        data_filename = f"{ensemble.name}{replica}_{ensemble.modes}modes_dibaryon.hdf5"
        data_file = os.path.join(data_dir, data_filename)
      elif ensemble.su3:
        data_filename = f"{ensemble.name}{replica}_{ensemble.modes}modes_baryon.hdf5"
        data_file = os.path.join(data_dir, data_filename)
      else:
        data_filename = f"{ensemble.dir_name}{replica}_{ensemble.modes}modes_S-2.hdf5"
        data_file = os.path.join(data_dir, data_filename)

    f_handler = h5py.File(data_file, 'r')
    data_handlers.append(f_handler)
    datasets.append(f_handler['c2_baryon'])

  flavors = ['octet']
  num_flavors = 1
  if not ensemble.su3:
    flavors = list(map(lambda x: x.decode('utf-8'), f_handler['baryons'][:]))
    num_flavors = len(flavors)

  ensemble_info = sig.MCEnsembleInfo(f"{ensemble.type}_{ensemble.name}", 'ensembles.xml')
  bins_info = sig.MCBinsInfo(ensemble_info)
  sampling_info = sig.MCSamplingInfo()
  xml_obs = sig.XMLHandler("MCObservables", "")
  obs_get_handler = sig.MCObsGetHandler(xml_obs, bins_info, sampling_info)
  obs_handler = sig.MCObsHandler(obs_get_handler, False)

  bin_dir = os.path.join(output_dir, ensemble.name, "single_particles")
  os.makedirs(bin_dir, exist_ok=True)
  bin_filename = f"{ensemble.name}{ensemble.replica_str}_{ensemble.modes}modes_baryon.bin"
  bin_file = os.path.join(bin_dir, bin_filename)
  obs_keys = set()

  for psq, pref in defs.total_prefs.items():
    irrep = defs.spin_half_irreps[psq]
    for flavor_i, flavor_name in enumerate(flavors):
      op_str = f"Pref={pref} {irrep} flavor={defs.baryon_flavor[flavor_name]} {flavor_name} 0"
      op = sig.OperatorInfo(op_str, sig.OpKind.GenIrrep)

      for t in range(ensemble.t0, ensemble.t0+ensemble.ts):
        bin_datas = list()
        for dataset in datasets:
          bin_data = np.zeros(dataset.shape[0], dtype=np.complex128)
          for tsrc in range(ensemble.srcs):
            for par in range(2):
              if skip_source(ensemble.name, tsrc, par):
                continue

              bin_data += dataset[:,tsrc,par,psq,flavor_i,t-ensemble.t0]
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


def convert_decuplet(ensemble):
  data_dir = os.path.join(base_data_dir, f"{ensemble.dir_name}")
  data_handlers = list()
  datasets = list()
  for replica in ensemble.replica:
    data_filename = f"{ensemble.dir_name}{replica}_{ensemble.modes}modes_decuplet.hdf5"
    data_file = os.path.join(data_dir, data_filename)
    f_handler = h5py.File(data_file, 'r')
    data_handlers.append(f_handler)
    datasets.append(f_handler['c2_decuplet'])

  flavor_name = "decuplet"
  num_flavors = 1 if ensemble.su3 else 4

  ensemble_info = sig.MCEnsembleInfo(f"{ensemble.type}_{ensemble.name}", 'ensembles.xml')
  bins_info = sig.MCBinsInfo(ensemble_info)
  sampling_info = sig.MCSamplingInfo()
  xml_obs = sig.XMLHandler("MCObservables", "")
  obs_get_handler = sig.MCObsGetHandler(xml_obs, bins_info, sampling_info)
  obs_handler = sig.MCObsHandler(obs_get_handler, False)

  bin_dir = os.path.join(output_dir, ensemble.name, "single_particles")
  os.makedirs(bin_dir, exist_ok=True)
  bin_filename = f"{ensemble.name}{ensemble.replica_str}_{ensemble.modes}modes_decuplet.bin"
  bin_file = os.path.join(bin_dir, bin_filename)
  obs_keys = set()

  for psq in range(1):
    pref = defs.total_prefs[psq]
    irrep = defs.spin_three_half_irreps[psq]
    for flavor_i in range(num_flavors):
      op_str = f"Pref={pref} {irrep} flavor={defs.baryon_flavor[flavor_name]} {flavor_name} {flavor_i}"
      op = sig.OperatorInfo(op_str, sig.OpKind.GenIrrep)

      for t in range(ensemble.t0, ensemble.t0+ensemble.ts):
        bin_datas = list()
        for dataset in datasets:
          bin_data = np.zeros(dataset.shape[0], dtype=np.complex128)
          for tsrc in range(ensemble.srcs):
            for par in range(2):
              if skip_source(ensemble.name, tsrc, par):
                continue

              bin_data += dataset[:,tsrc,par,psq,flavor_i,t-ensemble.t0]
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
  data_dir = os.path.join(base_data_dir, f"{ensemble.dir_name}")
  data_handlers = list()
  tsrc_list = dict()
  datasets = dict()
  for ps_name in defs.pseudoscalar_names[ensemble.name]:
    datasets[ps_name] = list()

  for rep_num, replica in enumerate(ensemble.replica):
    data_filename = f"{ensemble.dir_name}{replica}_{defs.pseudoscalar_modes[ensemble.name]}modes_pseudoscalar.hdf5"
    data_file = os.path.join(data_dir, data_filename)
    if not os.path.isfile(data_file) and ensemble.su3:
      data_filename = f"{ensemble.name}{replica}_{defs.pseudoscalar_modes[ensemble.name]}modes_pseudoscalar.hdf5"
      data_file = os.path.join(data_dir, data_filename)

    f_handler = h5py.File(data_file, 'r')
    data_handlers.append(f_handler)
    for ps_name in defs.pseudoscalar_names[ensemble.name]:
      datasets[ps_name].append(f_handler[ps_name])

    if ensemble.name in defs.tsrc_files:
      tsrc_list[rep_num] = dict()
      tsrc_file = os.path.join(data_dir, defs.tsrc_files[ensemble.name] + replica)
      with open(tsrc_file, 'r') as tsrc_file_h:
        for line in tsrc_file_h:
          conf, tsrc = line.split()
          tsrc_list[rep_num][int(conf)-1] = int(tsrc)

  ensemble_info = sig.MCEnsembleInfo(f"{ensemble.type}_{ensemble.name}", 'ensembles.xml')
  bins_info = sig.MCBinsInfo(ensemble_info)
  sampling_info = sig.MCSamplingInfo()
  xml_obs = sig.XMLHandler("MCObservables", "")
  obs_get_handler = sig.MCObsGetHandler(xml_obs, bins_info, sampling_info)
  obs_handler = sig.MCObsHandler(obs_get_handler, False)

  bin_dir = os.path.join(output_dir, ensemble.name, "single_particles")
  os.makedirs(bin_dir, exist_ok=True)
  bin_filename = f"{ensemble.name}{ensemble.replica_str}_{defs.pseudoscalar_modes[ensemble.name]}modes_pseudoscalar.bin"
  bin_file = os.path.join(bin_dir, bin_filename)
  obs_keys = set()
  for ps_name in defs.pseudoscalar_names[ensemble.name]:
    for t in range(ensemble.Nt):
      bin_datas = list()
      for rep_num, dataset in enumerate(datasets[ps_name]):
        bin_data = np.zeros(dataset.shape[0], dtype=np.complex128)
        for tsrc in range(ensemble.srcs):

          if tsrc_list:
            for config in range(dataset.shape[0]):
              t0 = defs.pseudoscalar_sources[ensemble.name][tsrc] + tsrc_list[rep_num][config]
              ts_forward = (t0+t) % ensemble.Nt
              ts_reverse = (t0-t) % ensemble.Nt

              if not ensemble.open:
                bin_data[config] += dataset[config,tsrc,ts_forward]
              elif use_forward_pseudoscalar(ensemble.name, tsrc) and use_backward_pseudoscalar(ensemble.name, tsrc):
                bin_data[config] += dataset[config,tsrc,ts_forward]
                bin_data[config] += dataset[config,tsrc,ts_reverse]
              elif use_forward_pseudoscalar(ensemble.name, tsrc):
                bin_data[config] += dataset[config,tsrc,ts_forward]
              elif use_backward_pseudoscalar(ensemble.name, tsrc):
                bin_data[config] += dataset[config,tsrc,ts_reverse]

          else:
            t0 = defs.pseudoscalar_sources[ensemble.name][tsrc]
            ts_forward = (t0+t) % ensemble.Nt
            ts_reverse = (t0-t) % ensemble.Nt

            if not ensemble.open:
              bin_data += dataset[:,tsrc,ts_forward]
            elif use_forward_pseudoscalar(ensemble.name, tsrc) and use_backward_pseudoscalar(ensemble.name, tsrc):
              bin_data += dataset[:,tsrc,ts_forward]
              bin_data += dataset[:,tsrc,ts_reverse]
            elif use_forward_pseudoscalar(ensemble.name, tsrc):
              bin_data += dataset[:,tsrc,ts_forward]
            elif use_backward_pseudoscalar(ensemble.name, tsrc):
              bin_data += dataset[:,tsrc,ts_reverse]

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
