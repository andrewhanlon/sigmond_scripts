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

def main():
  for ensemble in defs.ensembles:
    ensemble_name = ensemble.name

    corr_files = [[0 for replica in ensemble.replica] for tsrc in ensemble.sources]
    for replica_i, replica in enumerate(ensemble.replica):
      replica_ensemble_name = f"{ensemble_name}_{replica}"
      for tsrc_i, tsrc in enumerate(ensemble.sources):
        data_dir = os.path.join(replica, f"t0{tsrc}")
        search_dir = os.path.join(defs.base_data_dir, data_dir)
        corr_files[tsrc_i][replica_i] = get_corr_files(replica_ensemble_name, search_dir)

    # TODO: check all elements of correlators have same keys
    channels = corr_files[0][0].keys()

    for channel in tqdm.tqdm(channels):
      corrs_to_extend = list()
      op_lists = list()
      for replica_i, replica in enumerate(ensemble.replica):
        replica_ensemble_name = f"{ensemble_name}_{replica}"
        '''
        corrs_to_average = list()
        for tsrc_i, tsrc in enumerate(ensemble.sources):
          correlators = corr_files[tsrc_i][replica_i][channel]
          correlator_data, op_list = get_data(correlators, replica_ensemble_name, ensemble.Nt, tsrc)
          op_lists.append(op_list)
          corrs_to_average.append(correlator_data)

        averaged_corr_data = np.mean(corrs_to_average, axis=0)

        all_equal = all(op_list==op_lists[0] for op_list in op_lists)
        if not all_equal:
          print("not all op lists equal")
          exit()
        
        hdf5_file = get_hdf5_file(replica_ensemble_name)
        write_data(averaged_corr_data, channel, op_lists[0], hdf5_file)
        '''
        for tsrc_i, tsrc in enumerate(ensemble.sources):
          correlators = corr_files[tsrc_i][replica_i][channel]
          correlator_data, op_list = get_data(correlators, replica_ensemble_name, ensemble.Nt, tsrc)

          hdf5_file = get_hdf5_file(replica_ensemble_name, tsrc)
          write_data(correlator_data, channel, op_list, hdf5_file)


def write_data(data, channel, op_list, hdf5_file):
  f_hand = h5py.File(hdf5_file, 'a')
  channel_group = f_hand.create_group(channel.data_channel_str())
  channel_group.attrs.create('opList', op_list)
  channel_group.create_dataset('data', data=data)

  f_hand.close()



def get_data(correlators, ensemble_name, ensemble_Nt, tsrc):
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
      tmin = keys[0].getTimeIndex()
      tmax = keys[-1].getTimeIndex()

  operators = sorted(list([operator_info.operator.Operator(op) for op in operators]))
  operators = [op.operator_info for op in operators]

  single_correlator = True if len(correlators['correlators']) == 1 else False

  if single_correlator:
    corr_data = np.zeros((bins_info.getNumberOfBins(), tmax+1), dtype=np.complex128)
  else:
    corr_data = np.zeros((bins_info.getNumberOfBins(), len(operators), len(operators), tmax+1), dtype=np.complex128)


  for snk_i, snk_op in enumerate(operators):
    for src_i, src_op in enumerate(operators):

      correlator = sigmond.CorrelatorInfo(snk_op, src_op)
      correlator_opposite = sigmond.CorrelatorInfo(src_op, snk_op)
      for tsep in range(tmin, tmax+1):
        correlator_time = sigmond.CorrelatorAtTimeInfo(correlator, tsep, False, False)
        correlator_time_re_obsinfo = sigmond.MCObsInfo(correlator_time, sigmond.ComplexArg.RealPart)
        correlator_time_im_obsinfo = sigmond.MCObsInfo(correlator_time, sigmond.ComplexArg.ImaginaryPart)

        has_re = obs_handler.queryBins(correlator_time_re_obsinfo)
        has_im = obs_handler.queryBins(correlator_time_im_obsinfo)

        if has_re and has_im:
          data = np.array(obs_handler.getBins(correlator_time_re_obsinfo).array()) + 1j*np.array(obs_handler.getBins(correlator_time_im_obsinfo).array())
        elif has_re:
          data = np.array(obs_handler.getBins(correlator_time_re_obsinfo).array())
        elif has_im:
          data = 1j*np.array(obs_handler.getBins(correlator_time_im_obsinfo).array())
        else:
          continue

        if single_correlator:
          corr_data[:,tsep] = data
        else:
          corr_data[:,snk_i,src_i,tsep] = data

  operators = [op.op_str() for op in operators]
  return corr_data, operators


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
    channel = get_channel(corr)
    if channel not in corr_files:
      corr_files[channel] = dict()
      corr_files[channel]['data_files'] = list()
      corr_files[channel]['correlators'] = list()

    corr_files[channel]['data_files'].append(corr_handler.getFileName(corr))
    corr_files[channel]['correlators'].append(corr)

  return corr_files


def get_hdf5_file(ensemble_name, tsrc):
  output_dir = os.path.join(defs.output_dir, ensemble_name)
  os.makedirs(output_dir, exist_ok=True)

  hdf5_file = os.path.join(output_dir, f"{ensemble_name}_t0{tsrc}.hdf5")
  return hdf5_file

def get_channel(correlator):
  op = correlator.getSink().getBasicLapH()
  P = (op.getXMomentum(), op.getYMomentum(), op.getZMomentum())
  irrep = op.getLGIrrep()
  irrep_row = op.getLGIrrepRow()
  isospin = op.getIsospin()
  if isospin=="pion":
    isospin="triplet"

  if isospin=="kaon":
    isospin="doublet"

  strangeness = op.getStrangeness()
  channel = defs.Channel(P, irrep, irrep_row, isospin, strangeness)
  return channel


if __name__ == "__main__":
  main()
