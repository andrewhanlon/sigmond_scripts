#!/usr/bin/env python

import os

import h5py
import numpy as np
import xml.etree.ElementTree as ET

import defs

import sigmondbind as sig

base_data_dir = "/disk2/research/data/kpi/raw/"

FORCE_HERM = True

def main():
  for ensemble in defs.ensembles:
    ensemble_info = defs.data_info[ensemble]
    for ensemble_name in ensemble_info["ensemble_name"]:
      replica = ensemble_name.split('_')[-1]
      for isospin in defs.isospsins:
        search_dir = os.path.join(base_data_dir, ensemble, f"{isospin}_{replica}")

        sources = list()
        for source in ensemble_info["sources"]:
          sources.append(find_data(ensemble_name, search_dir, source))

        print("Averaging over sources")
        averaged_data = average_sources(sources)
        print("Writing to file")
        write_to_file(averaged_data, ensemble_name, isospin)
        print("done")
      for single_hadron in defs.single_hadrons:
        search_dir = os.path.join(base_data_dir, ensemble, f"{single_hadron}_{replica}")

        sources = list()
        for source in ensemble_info["sources"]:
          sources.append(find_data(ensemble_name, search_dir, source))

        print("Averaging over sources")
        averaged_data = average_sources(sources)
        # Change names
        new_op_list, new_data = update_single_hadrons(ensemble_name, averaged_data)
        print("Writing to file")
        write_sh_to_file(new_data, new_op_list, ensemble_name, single_hadron)
        print("done")

def update_single_hadrons(ensemble_name, averaged_data):
  tmax = 0
  num_ops = 0
  for channel, (op_list, data) in averaged_data.items():
    if len(op_list) != 1:
      print("Single hadrons cannot have multiple operators")
      exit()
    if data.shape[1] > tmax:
      tmax = data.shape[1]
    num_ops += 1

  new_data = np.zeros((len(defs.configs[ensemble_name]), tmax, num_ops), dtype=np.complex128)
  new_op_list = list()

  op_count = 0
  for channel, (op_list, data) in sorted(averaged_data.items()):
    new_op_list.append(f"pSq{channel[0]}")

    new_data[:,:,op_count] = data[:,:,0,0]
    op_count += 1

  return (new_op_list, new_data)

def write_sh_to_file(data, op_list, ensemble_name, single_hadron):
  os.makedirs(defs.output_dir, exist_ok=True)
  filename = os.path.join(defs.output_dir, f"{ensemble_name}_single_hadrons.hdf5")
  h5_file = h5py.File(filename, 'a')
  replica = ensemble_name.split('_')[-1]
  configs = [f"{replica}n{c_num}" for c_num in defs.configs[ensemble_name]]
  h5_file.attrs.create("cfgList", configs)
  h5_file.attrs.create("spatExt", defs.spatial_extent[ensemble_name])

  sh_group = h5_file.create_group(single_hadron)
  sh_group.attrs.create('opList', op_list)
  sh_group.create_dataset("t0avg", data=data)

  h5_file.close()

def write_to_file(averaged_data, ensemble_name, isospin):
  os.makedirs(defs.output_dir, exist_ok=True)
  filename = os.path.join(defs.output_dir, f"{ensemble_name}_{isospin}.hdf5")
  h5_file = h5py.File(filename, 'w')
  replica = ensemble_name.split('_')[-1]
  configs = [f"{replica}n{c_num}" for c_num in defs.configs[ensemble_name]]
  h5_file.attrs.create("cfgList", configs)
  h5_file.attrs.create("spatExt", defs.spatial_extent[ensemble_name])

  for channel, (op_list, data) in averaged_data.items():
    channel_group = h5_file.create_group(f"pSq{channel[0]}-{channel[1]}")
    channel_group.attrs.create('opList', op_list)
    channel_group.create_dataset("t0avg", data=data)

  h5_file.close()


def average_sources(sources):
  averaged_sources = dict()
  for channel in sources[0].keys():
    datas = list()
    for source in sources:
      datas.append(source[channel][1])

    averaged_data = sum(datas) / len(datas)
    averaged_sources[channel] = (sources[0][channel][0], averaged_data)

  return averaged_sources


def find_data(ensemble_name, search_dir, source):
  print(f"starting ensemble {ensemble_name} for source time {source}")
  ensemble_info = sig.MCEnsembleInfo(ensemble_name, 'ensembles.xml')

  print(f"Searching for correlators and time separations in {search_dir}")

  file_list_infos = dict()
  for root, dirs, files in os.walk(search_dir, topdown=True):
    for filename in files:
      base, ext = os.path.splitext(filename)
      if base != f"correlator_t0{source}":
        continue

      full_base = os.path.join(root, base)

      suffix = int(ext[1:])

      if full_base in file_list_infos:
        if suffix < file_list_infos[full_base][0]:
          file_list_infos[full_base][0] = suffix
        elif suffix > file_list_infos[full_base][1]:
          file_list_infos[full_base][1] = suffix
      else:
        file_list_infos[full_base] = [suffix, suffix]

  file_list_infos = [sig.FileListInfo(stub, suffix[0], suffix[1], False)
                     for stub, suffix in file_list_infos.items()]


  mcobs_xml = ET.Element("MCObservables")
  corr_data_xml = ET.SubElement(mcobs_xml, "BLCorrelatorData")
  for file_list_info in file_list_infos:
    corr_data_xml.append(file_list_info.xml())

  mcobs_xml_handler = sig.XMLHandler()
  mcobs_xml_handler.set_from_string(ET.tostring(mcobs_xml))
  sampling_info = sig.MCSamplingInfo()
  bins_info = sig.MCBinsInfo(ensemble_info)
  bins_info.addOmissions(defs.omissions[ensemble_name])
  obs_get_handler = sig.MCObsGetHandler(mcobs_xml_handler, bins_info, sampling_info)
  obs_handler = sig.MCObsHandler(obs_get_handler, False)

  corr_handler = sig.BLCorrelatorDataHandler(file_list_infos, set(), set(), ensemble_info)
  corrs = corr_handler.getCorrelatorSet()

  operators = dict()
  max_tsep = 0
  time_seps = dict()
  for corr in corrs:
    keys = corr_handler.getOrderedKeys(corr)
    tmin = keys[0].getTimeIndex()
    tmax = keys[-1].getTimeIndex()
    if tmax > max_tsep:
      max_tsep = tmax
    time_seps[corr] = (tmin, tmax)

    op_snk = corr.getSink()
    op_src = corr.getSource()
    op_snk_bl = op_snk.getBasicLapH()
    op_src_bl = op_src.getBasicLapH()
    P = (op_src_bl.getXMomentum(), op_src_bl.getYMomentum(), op_src_bl.getZMomentum())
    Psq = P[0]**2 + P[1]**2 + P[2]**2
    irrep = op_src_bl.getLGIrrep()
    irrep_row = op_src_bl.getLGIrrepRow()
    if (Psq, irrep) not in operators:
      operators[(Psq, irrep)] = dict()
    if (P, irrep_row) not in operators[(Psq, irrep)]:
      operators[(Psq, irrep)][(P, irrep_row)] = dict()
    op_snk_short_name = get_short_name(op_snk_bl)
    op_src_short_name = get_short_name(op_src_bl)
    if op_snk_short_name in operators[(Psq, irrep)][(P, irrep_row)] and operators[(Psq, irrep)][(P, irrep_row)][op_snk_short_name] != op_snk:
      print("short name conflict")
      print(operators[(Psq, irrep)][(P, irrep_row)][op_snk_short_name])
      print(op_snk)
      exit()
    if op_src_short_name in operators[(Psq, irrep)][(P, irrep_row)] and operators[(Psq, irrep)][(P, irrep_row)][op_src_short_name] != op_src:
      print("short name conflict")
      print(operators[(Psq, irrep)][(P, irrep_row)][op_src_short_name])
      print(op_src)
      exit()

    operators[(Psq, irrep)][(P, irrep_row)][op_snk_short_name] = op_snk
    operators[(Psq, irrep)][(P, irrep_row)][op_src_short_name] = op_src

  print("Done searching for all correlators and time separations")
  print("Now finding data")

  data = dict()
  for channel, channel_ops in operators.items():
    print(f"Working on {channel}")
    bin_datas = list()
    op_list = list()
    for equivalanet_channel, operator_dict in channel_ops.items():
      print(f"\tWorking on {equivalanet_channel}")

      if op_list:
        op_list_check = sorted(operator_dict.keys())
        if op_list_check != op_list:
          print("bad op list")
          exit()
      else:
        op_list = sorted(operator_dict.keys())

      bin_data = np.zeros((len(defs.configs[ensemble_name]), max_tsep+1, len(operator_dict), len(operator_dict)), dtype=np.complex128)
      for op_snk_i, (op_snk_short_name, op_snk) in enumerate(sorted(operator_dict.items())):
        snk_op_phase = defs.phases.get(op_snk.op_str(), 1)
        for op_src_i, (op_src_short_name, op_src) in enumerate(sorted(operator_dict.items())[op_snk_i:], start=op_snk_i):
          src_op_phase = defs.phases.get(op_src.op_str(), 1)
          phase = snk_op_phase*np.conj(src_op_phase)

          corr_info = sig.CorrelatorInfo(op_snk, op_src)
          corr_info_opposite = sig.CorrelatorInfo(op_src, op_snk)
          if corr_info_opposite not in time_seps:
            trange = time_seps[corr_info]
          elif corr_info not in time_seps:
            trange = time_seps[corr_info_opposite]
          else:
            trange1 = time_seps[corr_info]
            trange2 = time_seps[corr_info_opposite]
            trange = (max(trange1[0], trange2[0]), min(trange1[1], trange2[1]))

          corr_time_info = sig.CorrelatorAtTimeInfo(corr_info, 0, False, False)
          corr_time_info_opposite = sig.CorrelatorAtTimeInfo(corr_info_opposite, 0, False, False)
          for tsep in range(trange[0], trange[1]+1):
            corr_time_info.resetTimeSeparation(tsep)
            corr_time_info_opposite.resetTimeSeparation(tsep)

            obs_info_re = sig.MCObsInfo(corr_time_info, sig.ComplexArg.RealPart)
            obs_info_im = sig.MCObsInfo(corr_time_info, sig.ComplexArg.ImaginaryPart)

            obs_info_opposite_re = sig.MCObsInfo(corr_time_info_opposite, sig.ComplexArg.RealPart)
            obs_info_opposite_im = sig.MCObsInfo(corr_time_info_opposite, sig.ComplexArg.ImaginaryPart)

            if corr_info_opposite not in time_seps:
              the_data = np.array([phase*(re + im*1j) for re, im in zip(obs_handler.getBins(obs_info_re).array(), obs_handler.getBins(obs_info_im).array())])

              bin_data[:,tsep,op_snk_i,op_src_i] = the_data
              bin_data[:,tsep,op_src_i,op_snk_i] = np.conj(the_data)

            elif corr_info not in time_seps:
              the_data = np.array([np.conj(phase)*(re + im*1j) for re, im in zip(obs_handler.getBins(obs_info_opposite_re).array(), obs_handler.getBins(obs_info_opposite_im).array())])

              bin_data[:,tsep,op_snk_i,op_src_i] = np.conj(the_data)
              bin_data[:,tsep,op_src_i,op_snk_i] = the_data

            elif corr_info == corr_info_opposite:
              if FORCE_HERM:
                the_data = np.array(obs_handler.getBins(obs_info_re).array())
              else:
                the_data = np.array([(re + im*1j) for re, im in zip(obs_handler.getBins(obs_info_re).array(), obs_handler.getBins(obs_info_im).array())])

              bin_data[:,tsep,op_snk_i,op_src_i] = the_data

            else:
              the_data = np.array([phase*(re + im*1j) for re, im in zip(obs_handler.getBins(obs_info_re).array(), obs_handler.getBins(obs_info_im).array())])
              the_data_opposite = np.array([np.conj(phase)*(re + im*1j) for re, im in zip(obs_handler.getBins(obs_info_opposite_re).array(), obs_handler.getBins(obs_info_opposite_im).array())])

              if FORCE_HERM:
                bin_data[:,tsep,op_snk_i,op_src_i] = 0.5*(the_data + np.conj(the_data_opposite))
                bin_data[:,tsep,op_src_i,op_snk_i] = 0.5*(np.conj(the_data) + the_data_opposite)
              else:
                bin_data[:,tsep,op_snk_i,op_src_i] = the_data
                bin_data[:,tsep,op_src_i,op_snk_i] = the_data_opposite

      bin_datas.append(bin_data)
    averaged_data = sum(bin_datas) / len(bin_datas)

    data[channel] = (op_list, averaged_data)
    obs_handler.clearData()

  return data

def get_short_name(operator):
  if operator.getNumberOfHadrons() == 1:
    return f"{defs.name_map[operator.getFlavor()]} {operator.getHadronSpatialType(1)}_{operator.getHadronSpatialIdNumber(1)}"
  else:
    short_name = ""
    for had_num in range(1, operator.getNumberOfHadrons()+1):
      had_name = defs.name_map[operator.getHadronFlavor(had_num)]
      had_psq = operator.getHadronXMomentum(had_num)**2 + operator.getHadronYMomentum(had_num)**2 + operator.getHadronZMomentum(had_num)**2
      short_name += f"{had_name}({had_psq}) "

    short_name = short_name[:-1]
    if (cgid := operator.getLGClebschGordonIdNum()) > 0:
      short_name += f" {cgid}"

  return short_name


if __name__ == "__main__":
  main()
