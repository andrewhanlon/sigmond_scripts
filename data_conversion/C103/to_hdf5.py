#!/usr/bin/env python

import os

import numpy as np
import xml.etree.ElementTree as ET
import h5py

import tqdm

import defs
from sourceTimeList import srcTList

import sigmond

COMPLEX_ARGS = [sigmond.ComplexArg.RealPart, sigmond.ComplexArg.ImaginaryPart]

MAX_CORRS = 50

def main():
  for ensemble in defs.ensembles:
    ensemble_name = ensemble.name

    corr_files = [[0 for replica in ensemble.replica] for tsrc in ensemble.sources]
    for replica_i, replica in enumerate(ensemble.replica):
      replica_ensemble_name = f"{ensemble_name}_{replica}"
      for tsrc_i, tsrc in enumerate(ensemble.sources):
        #data_dir = os.path.join(ensemble_name, replica, f"src{tsrc_i}")
        data_dir = os.path.join(replica, f"src{tsrc_i}")
        search_dir = os.path.join(defs.base_data_dir, data_dir)
        corr_files[tsrc_i][replica_i] = get_corr_files(replica_ensemble_name, search_dir)

    # TODO: check all elements of correlators have same keys
    channels = corr_files[0][0].keys()

    for channel in tqdm.tqdm(channels):
      corrs_to_extend = list()
      op_lists = list()
      for replica_i, replica in enumerate(ensemble.replica):
        replica_ensemble_name = f"{ensemble_name}_{replica}"
        corrs_to_average = list()
        for tsrc_i, tsrc in enumerate(ensemble.sources):
          correlators = corr_files[tsrc_i][replica_i][channel]
          correlator_data, op_list = get_data(correlators, replica_ensemble_name, ensemble.Nt, tsrc)
          op_lists.append(op_list)
          corrs_to_average.append(correlator_data)

        averaged_corr_data = np.mean(corrs_to_average, axis=0)
        corrs_to_extend.append(averaged_corr_data)

      if len(corrs_to_extend) > 1:
        extended_corr_data = np.concatenate(corrs_to_extend)
      else:
        extended_corr_data = corrs_to_extend[0]

      all_equal = all(op_list==op_lists[0] for op_list in op_lists)
      if not all_equal:
        print("not all op lists equal")
        exit()
      
      hdf5_file = get_hdf5_file(channel, replica_ensemble_name)
      write_data(extended_corr_data, channel, op_lists[0], hdf5_file)


def write_data(data, channel, op_list, hdf5_file):
  f_hand = h5py.File(hdf5_file, 'w')
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

  operators = sorted(list(operators))


  single_correlator = True if len(correlators['correlators']) == 1 else False

  if single_correlator:
    corr_data = np.zeros((bins_info.getNumberOfBins(), tmax+1), dtype=np.complex128)
  else:
    corr_data = np.zeros((bins_info.getNumberOfBins(), len(operators), len(operators), tmax+1), dtype=np.complex128)


  for snk_i, snk_op in enumerate(operators):
    bl_op = snk_op.getBasicLapH()
    # TODO: replace with isFermionic()
    change_sign = bl_op.isBaryon() or bl_op.isMesonBaryon()

    for src_i, src_op in enumerate(operators[snk_i:], start=snk_i):

      correlator = sigmond.CorrelatorInfo(snk_op, src_op)
      correlator_opposite = sigmond.CorrelatorInfo(src_op, snk_op)
      for tsep in range(tmin, tmax+1):
        correlator_time = sigmond.CorrelatorAtTimeInfo(correlator, tsep, True, False)
        correlator_time_re_obsinfo = sigmond.MCObsInfo(correlator_time, sigmond.ComplexArg.RealPart)
        correlator_time_im_obsinfo = sigmond.MCObsInfo(correlator_time, sigmond.ComplexArg.ImaginaryPart)
        data_re = np.array(obs_handler.getBins(correlator_time_re_obsinfo).array())
        data_im = np.array(obs_handler.getBins(correlator_time_im_obsinfo).array())

        if snk_i != src_i:
          correlator_opposite_time = sigmond.CorrelatorAtTimeInfo(correlator_opposite, tsep, True, False)
          correlator_opposite_time_re_obsinfo = sigmond.MCObsInfo(correlator_opposite_time, sigmond.ComplexArg.RealPart)
          correlator_opposite_time_im_obsinfo = sigmond.MCObsInfo(correlator_opposite_time, sigmond.ComplexArg.ImaginaryPart)
          data_opposite_re = np.array(obs_handler.getBins(correlator_opposite_time_re_obsinfo).array())
          data_opposite_im = np.array(obs_handler.getBins(correlator_opposite_time_im_obsinfo).array())

          data_re = 0.5*(data_re + data_opposite_re)
          data_im = 0.5*(data_im - data_opposite_im)

        data = data_re + 1j*data_im

        if change_sign:
          for config in range(bins_info.getNumberOfBins()):
            T = tsep + tsrc + srcTList[config]
            if T >= ensemble_Nt:
              data[config] = -data[config]

        if single_correlator:
          corr_data[:,tsep] = data
        else:
          corr_data[:,snk_i,src_i,tsep] = data
          if snk_i != src_i:
            corr_data[:,src_i,snk_i,tsep] = np.conj(data)

  operators = [operator.op_str() for operator in operators]
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


def get_hdf5_file(channel, ensemble_name):
  output_dir = os.path.join(defs.output_dir, ensemble_name)
  os.makedirs(output_dir, exist_ok=True)

  hdf5_file = os.path.join(output_dir, f"{channel.iso_strange_str()}.hdf5")
  return hdf5_file

def get_channel(correlator):
  operator = correlator.getSink().getBasicLapH()
  P = (operator.getXMomentum(), operator.getYMomentum(), operator.getZMomentum())
  irrep = operator.getLGIrrep()
  irrep_row = operator.getLGIrrepRow()
  isospin = operator.getIsospin()
  strangeness = operator.getStrangeness()
  channel = defs.Channel(P, irrep, irrep_row, isospin, strangeness)
  return channel


if __name__ == "__main__":
  main()