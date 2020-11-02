#!/usr/bin/env python

import os

import h5py
import numpy as np
import xml.etree.ElementTree as ET
from sortedcontainers import SortedSet

import sigmond

import defs
import excluded_ops
import missing_sh_ops

base_data_dir = "/latticeQCD/raid0/ahanlon/data/D200/"

FORCE_HERM = False
AVERAGE_DATA = False
EXCLUDE_OPERATORS = True

def main():
  for ensemble in defs.ensembles:
    print(f"starting for ensemble {ensemble}")
    for channel in os.listdir(base_data_dir):
      search_dir = os.path.join(base_data_dir, channel)
      if channel == "single_hadrons":
        channel_data = find_sh_data(ensemble, search_dir)
        if AVERAGE_DATA:
          channel_data = average_sh_data(channel_data)
          write_averaged_sh_to_file(channel_data, ensemble)
        else:
          write_sh_to_file(channel_data, ensemble)
      
      else:
        channel_data = find_data(ensemble, search_dir)
        if AVERAGE_DATA:
          channel_data = average_data(channel_data)
          write_averaged_to_file(channel_data, ensemble, channel)
        else:
          write_to_file(channel_data, ensemble, channel)

      print("done")


def write_averaged_sh_to_file(data, ensemble_name):
  os.makedirs(defs.output_dir, exist_ok=True)
  filename = os.path.join(defs.output_dir, f"{ensemble_name}_single_hadrons.hdf5")
  h5_file = h5py.File(filename, 'w')
  replica = ensemble_name.split('_')[-1]
  configs = [f"{replica}n{c_num}" for c_num in defs.configs[ensemble_name]]
  h5_file.attrs.create("cfgList", configs)
  h5_file.attrs.create("spatExt", defs.spatial_extent[ensemble_name])
  op_list = list()
  for flavor, averaged_channels in data.items():
    flavor_group = h5_file.create_group(flavor)
    bin_datas = list()
    op_list = list()
    for averaged_channel, bin_data in averaged_channels.items():
      bin_datas.append(bin_data)
      op_list.append(f"pSq{averaged_channel.psq}")
    flavor_group.attrs.create('opList', op_list)
    flavor_group.create_dataset("data", data=np.stack(bin_datas,axis=-1))

  h5_file.close()

def write_sh_to_file(data, ensemble_name):
  os.makedirs(defs.output_dir, exist_ok=True)
  filename = os.path.join(defs.output_dir, f"{ensemble_name}_single_hadrons.hdf5")
  h5_file = h5py.File(filename, 'w')
  replica = ensemble_name.split('_')[-1]
  configs = [f"{replica}n{c_num}" for c_num in defs.configs[ensemble_name]]
  h5_file.attrs.create("cfgList", configs)
  h5_file.attrs.create("spatExt", defs.spatial_extent[ensemble_name])
  for flavor, averaged_channels in data.items():
    flavor_group = h5_file.create_group(flavor)
    bin_datas = list()
    for averaged_channel, equivalent_frames in averaged_channels.items():
      average_channel_group = flavor_group.create_group(repr(averaged_channel))
      for equivalent_frame, bin_data in equivalent_frames.items():
        equivalent_group = average_channel_group.create_group(repr(equivalent_frame))
        equivalent_group.create_dataset("data", data=bin_data)

  h5_file.close()

def write_averaged_to_file(the_data, ensemble_name, the_channel):
  os.makedirs(defs.output_dir, exist_ok=True)
  filename = os.path.join(defs.output_dir, f"{ensemble_name}_{the_channel}.hdf5")
  h5_file = h5py.File(filename, 'w')
  replica = ensemble_name.split('_')[-1]
  configs = [f"{replica}n{c_num}" for c_num in defs.configs[ensemble_name]]
  h5_file.attrs.create("cfgList", configs)
  h5_file.attrs.create("spatExt", defs.spatial_extent[ensemble_name])

  for channel, (op_list, data) in the_data.items():
    channel_group = h5_file.create_group(repr(channel))
    channel_group.attrs.create('opList', op_list)
    channel_group.create_dataset("data", data=data)

  h5_file.close()

def write_to_file(the_data, ensemble_name, the_channel):
  os.makedirs(defs.output_dir, exist_ok=True)
  filename = os.path.join(defs.output_dir, f"{ensemble_name}_{the_channel}.hdf5")
  h5_file = h5py.File(filename, 'w')
  replica = ensemble_name.split('_')[-1]
  configs = [f"{replica}n{c_num}" for c_num in defs.configs[ensemble_name]]
  h5_file.attrs.create("cfgList", configs)
  h5_file.attrs.create("spatExt", defs.spatial_extent[ensemble_name])

  for channel, equivalent_frames in the_data.items():
    channel_group = h5_file.create_group(repr(channel))
    for equivalent_frame, (op_list, data) in equivalent_frames.items():
      op_list = [op.op_str() for op in op_list]
      equivalent_frame_group = channel_group.create_group(repr(equivalent_frame))
      equivalent_frame_group.attrs.create('opList', op_list)
      equivalent_frame_group.create_dataset("data", data=data)

  h5_file.close()

def find_data(ensemble_name, search_dir):
  ensemble_info = sigmond.MCEnsembleInfo(ensemble_name, 'ensembles.xml')

  # search for data files
  print(f"Searching for correlators and time separations in {search_dir}")
  file_list_infos = dict()
  for root, dirs, files in os.walk(search_dir, topdown=True):
    for filename in files:
      base, ext = os.path.splitext(filename)
      full_base = os.path.join(root, base)

      try:
        suffix = int(ext[1:])
      except ValueError:
        continue

      if full_base in file_list_infos:
        if suffix < file_list_infos[full_base][0]:
          file_list_infos[full_base][0] = suffix
        elif suffix > file_list_infos[full_base][1]:
          file_list_infos[full_base][1] = suffix
      else:
        file_list_infos[full_base] = [suffix, suffix]

  file_list_infos = [sigmond.FileListInfo(stub, suffix[0], suffix[1], False)
                     for stub, suffix in file_list_infos.items()]


  # create sigmond handlers
  mcobs_xml = ET.Element("MCObservables")
  corr_data_xml = ET.SubElement(mcobs_xml, "BLCorrelatorData")
  for file_list_info in file_list_infos:
    corr_data_xml.append(file_list_info.xml())

  mcobs_xml_handler = sigmond.XMLHandler()
  mcobs_xml_handler.set_from_string(ET.tostring(mcobs_xml))
  sampling_info = sigmond.MCSamplingInfo()
  bins_info = sigmond.MCBinsInfo(ensemble_info)
  bins_info.addOmissions(defs.omissions[ensemble_name])
  obs_get_handler = sigmond.MCObsGetHandler(mcobs_xml_handler, bins_info, sampling_info)
  obs_handler = sigmond.MCObsHandler(obs_get_handler, False)

  corr_handler = sigmond.BLCorrelatorDataHandler(file_list_infos, set(), set(), ensemble_info)
  corrs = corr_handler.getCorrelatorSet()

  # search through found data
  operators = dict()
  max_tsep = 0
  time_seps = dict()
  for corr in corrs:
    op_snk = corr.getSink()
    op_src = corr.getSource()
    if EXCLUDE_OPERATORS and (op_snk.op_str() in excluded_ops.excluded_operators or op_src.op_str() in excluded_ops.excluded_operators):
      continue

    keys = corr_handler.getOrderedKeys(corr)
    tmin = keys[0].getTimeIndex()
    tmax = keys[-1].getTimeIndex()
    if tmax > max_tsep:
      max_tsep = tmax
    time_seps[corr] = (tmin, tmax)

    op_snk_bl = op_snk.getBasicLapH()
    op_src_bl = op_src.getBasicLapH()
    P = (op_src_bl.getXMomentum(), op_src_bl.getYMomentum(), op_src_bl.getZMomentum())
    Psq = P[0]**2 + P[1]**2 + P[2]**2
    irrep = op_src_bl.getLGIrrep()
    irrep_row = op_src_bl.getLGIrrepRow()
    averaged_channel = defs.AveragedChannel(Psq, irrep)
    equivalent_frame = defs.EquivalentFrame(P, irrep_row)
    if averaged_channel not in operators:
      operators[averaged_channel] = dict()
    if equivalent_frame not in operators[averaged_channel]:
      operators[averaged_channel][equivalent_frame] = SortedSet()

    operators[averaged_channel][equivalent_frame].add(op_snk)
    operators[averaged_channel][equivalent_frame].add(op_src)

  print("Done searching for all correlators and time separations")
  print("Now finding data")

  data = dict()
  for averaged_channel, channel_ops in operators.items():
    print(f"Working on {averaged_channel}")
    data[averaged_channel] = dict()
    for equivalent_frame, operators in channel_ops.items():
      print(f"\tWorking on {equivalent_frame}")
      operator_list = list(operators)

      bin_data = np.zeros((len(defs.configs[ensemble_name]), max_tsep+1, len(operator_list), len(operator_list)), dtype=np.complex128)
      for op_snk_i, op_snk in enumerate(operator_list):
        snk_op_phase = defs.phases.get(op_snk.op_str(), 1)
        for op_src_i, op_src in enumerate(operator_list[op_snk_i:], start=op_snk_i):
          src_op_phase = defs.phases.get(op_src.op_str(), 1)
          phase = snk_op_phase*np.conj(src_op_phase)

          corr_info = sigmond.CorrelatorInfo(op_snk, op_src)
          corr_info_opposite = sigmond.CorrelatorInfo(op_src, op_snk)
          if corr_info_opposite not in time_seps:
            trange = time_seps[corr_info]
          elif corr_info not in time_seps:
            trange = time_seps[corr_info_opposite]
          else:
            trange1 = time_seps[corr_info]
            trange2 = time_seps[corr_info_opposite]
            trange = (max(trange1[0], trange2[0]), min(trange1[1], trange2[1]))

          corr_time_info = sigmond.CorrelatorAtTimeInfo(corr_info, 0, False, False)
          corr_time_info_opposite = sigmond.CorrelatorAtTimeInfo(corr_info_opposite, 0, False, False)
          for tsep in range(trange[0], trange[1]+1):
            corr_time_info.resetTimeSeparation(tsep)
            corr_time_info_opposite.resetTimeSeparation(tsep)

            obs_info_re = sigmond.MCObsInfo(corr_time_info, sigmond.ComplexArg.RealPart)
            obs_info_im = sigmond.MCObsInfo(corr_time_info, sigmond.ComplexArg.ImaginaryPart)

            obs_info_opposite_re = sigmond.MCObsInfo(corr_time_info_opposite, sigmond.ComplexArg.RealPart)
            obs_info_opposite_im = sigmond.MCObsInfo(corr_time_info_opposite, sigmond.ComplexArg.ImaginaryPart)

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

      data[averaged_channel][equivalent_frame] = (operator_list, bin_data)
      obs_handler.clearData()

  return data

def find_sh_data(ensemble_name, search_dir):
  ensemble_info = sigmond.MCEnsembleInfo(ensemble_name, 'ensembles.xml')

  # search for data files
  print(f"Searching for correlators and time separations in {search_dir}")
  file_list_infos = dict()
  for root, dirs, files in os.walk(search_dir, topdown=True):
    for filename in files:
      base, ext = os.path.splitext(filename)
      full_base = os.path.join(root, base)

      try:
        suffix = int(ext[1:])
      except ValueError:
        continue

      if full_base in file_list_infos:
        if suffix < file_list_infos[full_base][0]:
          file_list_infos[full_base][0] = suffix
        elif suffix > file_list_infos[full_base][1]:
          file_list_infos[full_base][1] = suffix
      else:
        file_list_infos[full_base] = [suffix, suffix]

  file_list_infos = [sigmond.FileListInfo(stub, suffix[0], suffix[1], False)
                     for stub, suffix in file_list_infos.items()]


  # create sigmond handlers
  mcobs_xml = ET.Element("MCObservables")
  corr_data_xml = ET.SubElement(mcobs_xml, "BLCorrelatorData")
  for file_list_info in file_list_infos:
    corr_data_xml.append(file_list_info.xml())

  mcobs_xml_handler = sigmond.XMLHandler()
  mcobs_xml_handler.set_from_string(ET.tostring(mcobs_xml))
  sampling_info = sigmond.MCSamplingInfo()
  bins_info = sigmond.MCBinsInfo(ensemble_info)
  bins_info.addOmissions(defs.omissions[ensemble_name])
  obs_get_handler = sigmond.MCObsGetHandler(mcobs_xml_handler, bins_info, sampling_info)
  obs_handler = sigmond.MCObsHandler(obs_get_handler, False)

  corr_handler = sigmond.BLCorrelatorDataHandler(file_list_infos, set(), set(), ensemble_info)
  corrs = corr_handler.getCorrelatorSet()

  # search through found data
  operators = dict()
  max_tsep = 0
  time_seps = dict()
  for corr in corrs:
    if not corr.isSinkSourceSame():
      print("aren't sh ops diagonal?")
      exit()

    op_src = corr.getSource()
    if op_src.op_str() in missing_sh_ops.missing_operators:
      continue

    keys = corr_handler.getOrderedKeys(corr)
    tmin = keys[0].getTimeIndex()
    tmax = keys[-1].getTimeIndex()
    if tmax > max_tsep:
      max_tsep = tmax
    time_seps[corr] = (tmin, tmax)

    op_src = corr.getSource()
    op_src_bl = op_src.getBasicLapH()
    P = (op_src_bl.getXMomentum(), op_src_bl.getYMomentum(), op_src_bl.getZMomentum())
    Psq = P[0]**2 + P[1]**2 + P[2]**2
    irrep = op_src_bl.getLGIrrep()
    irrep_row = op_src_bl.getLGIrrepRow()
    flavor = op_src_bl.getFlavor()
    averaged_channel = defs.AveragedChannel(Psq, irrep)
    equivalent_frame = defs.EquivalentFrame(P, irrep_row)
    if flavor not in operators:
      operators[flavor] = dict()
    if averaged_channel not in operators[flavor]:
      operators[flavor][averaged_channel] = dict()
    if equivalent_frame in operators[flavor][averaged_channel]:
      print("aren't sh ops 1x1?")
      exit()

    operators[flavor][averaged_channel][equivalent_frame] = op_src

  print("Done searching for all correlators and time separations")
  print("Now finding data")

  data = dict()
  for flavor, averaged_channels in operators.items():
    print(f"Working on {flavor}")
    data[flavor] = dict()
    for averaged_channel, channel_ops in averaged_channels.items():
      print(f"\tWorking on {averaged_channel}")
      data[flavor][averaged_channel] = dict()
      for equivalent_frame, operator in channel_ops.items():
        print(f"\t\tWorking on {equivalent_frame}")
        bin_data = np.zeros((len(defs.configs[ensemble_name]), max_tsep+1), dtype=np.complex128)
        corr_info = sigmond.CorrelatorInfo(operator, operator)
        trange = time_seps[corr_info]
        corr_time_info = sigmond.CorrelatorAtTimeInfo(corr_info, 0, False, False)
        for tsep in range(trange[0], trange[1]+1):
          corr_time_info.resetTimeSeparation(tsep)

          if FORCE_HERM:
            obs_info_re = sigmond.MCObsInfo(corr_time_info, sigmond.ComplexArg.RealPart)

            bin_data[:,tsep] = np.array([obs_handler.getBins(obs_info_re).array()])
          else:
            obs_info_re = sigmond.MCObsInfo(corr_time_info, sigmond.ComplexArg.RealPart)
            obs_info_im = sigmond.MCObsInfo(corr_time_info, sigmond.ComplexArg.ImaginaryPart)

            bin_data[:,tsep] = np.array([(re + im*1j) for re, im in zip(obs_handler.getBins(obs_info_re).array(), obs_handler.getBins(obs_info_im).array())])


        data[flavor][averaged_channel][equivalent_frame] = bin_data
        obs_handler.clearData()

  return data

def average_data(data):
  averaged_data = dict()
  for averaged_channel, equivalent_frames in data.items():
    averaged_data[averaged_channel] = dict()
    bin_datas = list()
    averaged_op_list = list()
    for equivalent_frame, (operators, bin_data) in equivalent_frames.items():
      averaged_op_list, reorder_list = get_averaged_op_list(operators, averaged_op_list)
      bin_data = bin_data[:,:,reorder_list,:]
      bin_data = bin_data[:,:,:,reorder_list]
      bin_datas.append(bin_data)

    temp_averaged_data = sum(bin_datas) / len(bin_datas)
    averaged_data[averaged_channel] = (averaged_op_list, temp_averaged_data)

  return averaged_data

def average_sh_data(data):
  averaged_data = dict()
  for flavor, averaged_channels in data.items():
    averaged_data[flavor] = dict()
    for averaged_channel, equivalent_frames in averaged_channels.items():
      averaged_data[flavor][averaged_channel] = dict()
      bin_datas = list()
      for equivalent_frame, bin_data in equivalent_frames.items():
        bin_datas.append(bin_data)

      temp_averaged_data = sum(bin_datas) / len(bin_datas)
      averaged_data[flavor][averaged_channel] = temp_averaged_data

  return averaged_data

def get_averaged_op_list(operators, current_averaged_list):
  averaged_op_list = list()
  for operator in operators:
    short_op_name = get_short_name(operator)
    if short_op_name in averaged_op_list:
      averaged_op_list = get_alt_short_names(operators)
      break

    averaged_op_list.append(short_op_name)

  reorder_list = sorted(range(len(averaged_op_list)), key=lambda k: averaged_op_list[k])
  averaged_op_list = sorted(averaged_op_list)
  if current_averaged_list and averaged_op_list != current_averaged_list:
    print("mismatch")
    exit()

  return averaged_op_list, reorder_list

def get_alt_short_names(operators):
  short_name_list = list()
  for operator in operators:
    short_op_name = get_alt_short_name(operator)
    if short_op_name in short_name_list:
      print("name conflict")
      exit()

    short_name_list.append(short_op_name)

  return short_name_list
      

def get_short_name(operator):
  operator = operator.getBasicLapH()
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

def get_alt_short_name(operator):
  operator = operator.getBasicLapH()
  if operator.getNumberOfHadrons() == 1:
    return f"{defs.name_map[operator.getFlavor()]} {operator.getHadronSpatialType(1)}_{operator.getHadronSpatialIdNumber(1)}"
  else:
    short_name = ""
    for had_num in range(1, operator.getNumberOfHadrons()+1):
      had_name = defs.name_map[operator.getHadronFlavor(had_num)]
      had_psq = operator.getHadronXMomentum(had_num)**2 + operator.getHadronYMomentum(had_num)**2 + operator.getHadronZMomentum(had_num)**2
      spat_type = operator.getHadronSpatialType(had_num)
      spat_id = operator.getHadronSpatialIdNumber(had_num)

      short_name += f"{had_name}({had_psq},{spat_type}{spat_id}) "

    short_name = short_name[:-1]
    if (cgid := operator.getLGClebschGordonIdNum()) > 0:
      short_name += f" {cgid}"

  return short_name



if __name__ == "__main__":
  main()
