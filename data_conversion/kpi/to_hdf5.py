#!/usr/bin/env python

import os

import numpy as np
import xml.etree.ElementTree as ET
import h5py

import tqdm

import defs

import sigmond

import sys
sys.path.insert(1, "../../analysis/")
import operator_info.operator

COMPLEX_ARGS = [sigmond.ComplexArg.RealPart, sigmond.ComplexArg.ImaginaryPart]

MAX_CORRS = 50

AVE_CORRS = True

def main():
  for ensemble in defs.ensembles:
    ensemble_name = ensemble.name
    print(f"Working on {ensemble_name}")

    corr_files = [[0 for replica in ensemble.replica] for tsrc in ensemble.sources]
    channels = list()
    for replica_i, replica in enumerate(ensemble.replica):
      replica_ensemble_name = f"{ensemble_name}_{replica}"
      for tsrc_i, (tsrc, parity) in enumerate(ensemble.sources):
        data_dir = os.path.join(ensemble_name, replica, f"src{tsrc}", parity)
        search_dir = os.path.join(defs.base_data_dir, data_dir)
        print(f"Searching in {search_dir}")
        corr_files[tsrc_i][replica_i] = get_corr_files(replica_ensemble_name, search_dir)
        channels.append(corr_files[tsrc_i][replica_i].keys())

    if not all(chans == channels[0] for chans in channels):
      print("Not all channels equal")
      sys.exit()

    channels = channels[0]

    for channel in tqdm.tqdm(channels):
      corrs_to_extend = list()
      op_lists = list()
      for replica_i, replica in enumerate(ensemble.replica):
        replica_ensemble_name = f"{ensemble_name}_{replica}"

        if AVE_CORRS:
          corrs_to_average = list()
          for tsrc_i, (tsrc, parity) in enumerate(ensemble.sources):
            correlators = corr_files[tsrc_i][replica_i][channel]
            correlator_data, op_list = get_data(correlators, replica_ensemble_name, tsrc, parity)
            op_lists.append(op_list)
            corrs_to_average.append(correlator_data)

          averaged_corr_data = np.mean(corrs_to_average, axis=0)

          all_equal = all(op_list==op_lists[0] for op_list in op_lists)
          if not all_equal:
            print("not all op lists equal")
            exit()
          
          hdf5_file = get_ave_hdf5_file(replica_ensemble_name, channel)
          write_data(averaged_corr_data, channel, op_lists[0], hdf5_file, replica, replica_ensemble_name)

        else:
          for tsrc_i, (tsrc, parity) in enumerate(ensemble.sources):
            correlators = corr_files[tsrc_i][replica_i][channel]
            correlator_data, op_list = get_data(correlators, replica_ensemble_name, tsrc, parity)

            hdf5_file = get_hdf5_file(replica_ensemble_name, channel, tsrc, parity)
            write_data(correlator_data, channel, op_list, hdf5_file, replica, replica_ensemble_name)


def write_data(data, channel, op_list, hdf5_file, replica, replica_ensemble_name):
  exists = os.path.isfile(hdf5_file)

  f_hand = h5py.File(hdf5_file, 'a')
  if not exists:
    configs = [f"{replica}n{c_num}" for c_num in defs.configs[replica_ensemble_name]]
    f_hand.attrs.create("cfgList", configs)

  channel_group = f_hand.create_group(channel.data_channel_str())
  channel_group.attrs.create('opList', op_list)
  channel_group.create_dataset('data', data=data)

  f_hand.close()


def get_data(correlators, ensemble_name, tsrc, parity):
  mcobs_xml = ET.Element("MCObservables")
  corr_data_xml = ET.SubElement(mcobs_xml, "BLCorrelatorData")

  file_list_infos = list()

  for data_file in correlators['data_files']:
    stub, ext = os.path.splitext(data_file)
    suffix = int(ext[1:])
    file_list_info = sigmond.FileListInfo(stub, suffix, suffix, False)
    corr_data_xml.append(file_list_info.xml())
    file_list_infos.append(file_list_info)

  mcobs_xml_handler = sigmond.XMLHandler()
  mcobs_xml_handler.set_from_string(ET.tostring(mcobs_xml))
  sampling_info = sigmond.MCSamplingInfo()
  ensemble_info = sigmond.MCEnsembleInfo(ensemble_name, 'ensembles.xml')
  bins_info = sigmond.MCBinsInfo(ensemble_info)
  bins_info.addOmissions(defs.omissions[ensemble_name])
  obs_get_handler = sigmond.MCObsGetHandler(mcobs_xml_handler, bins_info, sampling_info)
  obs_handler = sigmond.MCObsHandler(obs_get_handler, False)

  corr_handler = sigmond.BLCorrelatorDataHandler(file_list_infos, set(), set(), ensemble_info)

  operators = set()
  tmin = None
  tmax = None
  for corr in correlators['correlators']:
    src_op_str = corr.getSource()
    snk_op_str = corr.getSink()
    operators.add(src_op_str)
    operators.add(snk_op_str)

    if tmin is None:
      keys = corr_handler.getOrderedKeys(corr)
      try:
        tmin = keys[0].getTimeIndex()
        tmax = keys[-1].getTimeIndex()
      except IndexError:
        print(f"no keys for {corr.corr_str()}")
        print("in files:")
        print(file_list_infos)
        #sys.exit()
        continue

  operators = sorted(list([operator_info.operator.Operator(op) for op in operators]))
  operators = [op.operator_info for op in operators]

  single_correlator = True if len(correlators['correlators']) == 1 else False

  if single_correlator:
    corr_data = np.zeros((bins_info.getNumberOfBins(), tmax+1), dtype=np.complex128)
  else:
    corr_data = np.zeros((bins_info.getNumberOfBins(), len(operators), len(operators), tmax+1), dtype=np.complex128)


  for snk_i, snk_op in enumerate(operators):
    for src_i, src_op in enumerate(operators[snk_i:], start=snk_i):

      correlator = sigmond.CorrelatorInfo(snk_op, src_op)
      correlator_opp = sigmond.CorrelatorInfo(src_op, snk_op)
      for tsep in range(tmin, tmax+1):
        correlator_time = sigmond.CorrelatorAtTimeInfo(correlator, tsep, False, False)
        correlator_time_re_obsinfo = sigmond.MCObsInfo(correlator_time, sigmond.ComplexArg.RealPart)
        correlator_time_im_obsinfo = sigmond.MCObsInfo(correlator_time, sigmond.ComplexArg.ImaginaryPart)

        correlator_opp_time = sigmond.CorrelatorAtTimeInfo(correlator_opp, tsep, False, False)
        correlator_opp_time_re_obsinfo = sigmond.MCObsInfo(correlator_opp_time, sigmond.ComplexArg.RealPart)
        correlator_opp_time_im_obsinfo = sigmond.MCObsInfo(correlator_opp_time, sigmond.ComplexArg.ImaginaryPart)

        has_re = obs_handler.queryBins(correlator_time_re_obsinfo)
        has_im = obs_handler.queryBins(correlator_time_im_obsinfo)
        if (has_re and not has_im) or (has_im and not has_re):
          print("only real or imaginary")
          sys.exit()

        has_corr = has_re and has_im

        has_opp_re = obs_handler.queryBins(correlator_opp_time_re_obsinfo)
        has_opp_im = obs_handler.queryBins(correlator_opp_time_im_obsinfo)
        if (has_opp_re and not has_opp_im) or (has_opp_im and not has_opp_re):
          print("only real or imaginary")
          sys.exit()

        has_opp_corr = has_opp_re and has_opp_im

        if has_corr and has_opp_corr:
          data_not_opp = np.array(obs_handler.getBins(correlator_time_re_obsinfo).array()) + 1j*np.array(obs_handler.getBins(correlator_time_im_obsinfo).array())
          data_opp = np.array(obs_handler.getBins(correlator_opp_time_re_obsinfo).array()) - 1j*np.array(obs_handler.getBins(correlator_opp_time_im_obsinfo).array())

          data = 0.5*(data_not_opp + data_opp)

        elif has_corr:
          data = np.array(obs_handler.getBins(correlator_time_re_obsinfo).array()) + 1j*np.array(obs_handler.getBins(correlator_time_im_obsinfo).array())

        elif has_opp_corr:
          data = np.array(obs_handler.getBins(correlator_opp_time_re_obsinfo).array()) - 1j*np.array(obs_handler.getBins(correlator_opp_time_im_obsinfo).array())
        else:
          continue
          

        if single_correlator:
          corr_data[:,tsep] = data
        else:
          corr_data[:,snk_i,src_i,tsep] = data
          corr_data[:,src_i,snk_i,tsep] = np.conj(data)

  op_list = list()
  for operator in operators:
    operator.setForwards()
    op_list.append(operator.op_str())

  return corr_data, op_list


def get_corr_files(ensemble_name, search_dir):
  file_list_infos = dict()
  for root, dirs, files in os.walk(search_dir, topdown=True):
    for filename in files:
      full_filename = os.path.join(root, filename)
      try:
        file_type = sigmond.getFileID(full_filename)
      except ValueError:
        continue

      if file_type != sigmond.FileType.Correlator:
        continue


      base, ext = os.path.splitext(full_filename)

      suffix = int(ext[1:])

      if base in file_list_infos:
        if suffix < file_list_infos[base][0]:
          file_list_infos[base][0] = suffix
        elif suffix > file_list_infos[base][1]:
          file_list_infos[base][1] = suffix
      else:
        file_list_infos[base] = [suffix, suffix]

  file_list_infos = [sigmond.FileListInfo(stub, suffix[0], suffix[1], False)
                     for stub, suffix in file_list_infos.items()]

  mcobs_xml = ET.Element("MCObservables")
  corr_data_xml = ET.SubElement(mcobs_xml, "BLCorrelatorData")
  for file_list_info in file_list_infos:
    corr_data_xml.append(file_list_info.xml())

  mcobs_xml_handler = sigmond.XMLHandler()
  mcobs_xml_handler.set_from_string(ET.tostring(mcobs_xml))
  ensemble_info = sigmond.MCEnsembleInfo(ensemble_name, 'ensembles.xml')
  sampling_info = sigmond.MCSamplingInfo()
  bins_info = sigmond.MCBinsInfo(ensemble_info)
  bins_info.addOmissions(defs.omissions[ensemble_name])
  obs_get_handler = sigmond.MCObsGetHandler(mcobs_xml_handler, bins_info, sampling_info)
  obs_handler = sigmond.MCObsHandler(obs_get_handler, False)

  corr_handler = sigmond.BLCorrelatorDataHandler(file_list_infos, set(), set(), ensemble_info)

  corr_files = dict()
  for corr in corr_handler.getCorrelatorSet():
    '''
    if corr.getSource().isBackwards() or corr.getSink().isBackwards():
      print(corr.corr_str())
      print(corr_handler.getFileName(corr))
      #sys.exit()
    '''

    channel = get_channel(corr)
    if channel not in corr_files:
      corr_files[channel] = dict()
      corr_files[channel]['data_files'] = list()
      corr_files[channel]['correlators'] = list()

    corr_files[channel]['data_files'].append(corr_handler.getFileName(corr))
    corr_files[channel]['correlators'].append(corr)

  '''
  # DEBUG
  for channel, corrs in corr_files.items():
    print(f"Channel {channel}")
    for x in range(len(corrs['data_files'])):
      print(f"  Data file: {corrs['data_files'][x]}")
      print(f"  corr: {corrs['correlators'][x].corr_str()}")
      print()
  '''

  return corr_files

def get_ave_hdf5_file(ensemble_name, channel):
  output_dir = os.path.join(defs.output_dir, 'hdf5', ensemble_name)
  os.makedirs(output_dir, exist_ok=True)

  hdf5_file = os.path.join(output_dir, f"{ensemble_name}_{channel.iso_strange_str()}_ave.hdf5")
  return hdf5_file

def get_hdf5_file(ensemble_name, channel, tsrc, parity):
  output_dir = os.path.join(defs.output_dir, 'hdf5', ensemble_name)
  os.makedirs(output_dir, exist_ok=True)

  hdf5_file = os.path.join(output_dir, f"{ensemble_name}_{channel.iso_strange_str()}_t0{tsrc}_{parity}.hdf5")
  return hdf5_file

def get_channel(correlator):
  op = correlator.getSink().getBasicLapH()
  P = (op.getXMomentum(), op.getYMomentum(), op.getZMomentum())
  irrep = op.getLGIrrep()
  irrep_row = op.getLGIrrepRow()
  isospin = op.getIsospin()
  isospin = defs.name_map.get(isospin, isospin)
  strangeness = op.getStrangeness()
  channel = defs.Channel(P, irrep, irrep_row, isospin, strangeness)
  return channel


if __name__ == "__main__":
  main()
