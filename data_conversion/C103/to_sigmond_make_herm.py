#!/usr/bin/env python

import os
import psutil

import numpy as np
import xml.etree.ElementTree as ET

import tqdm

import defs
from sourceTimeList import srcTList

import sigmond

COMPLEX_ARGS = [sigmond.ComplexArg.RealPart, sigmond.ComplexArg.ImaginaryPart]

MAX_CORRS = 50

def main():
  process = psutil.Process(os.getpid())
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
    correlators = corr_files[0][0].keys()

    data_to_write = dict()

    for corr_num, correlator in enumerate(tqdm.tqdm(correlators), 1):
      corrs_to_extend = list()
      for replica_i, replica in enumerate(ensemble.replica):
        replica_ensemble_name = f"{ensemble_name}_{replica}"
        corrs_to_average = list()
        for tsrc_i, tsrc in enumerate(ensemble.sources):
          data_files = corr_files[tsrc_i][replica_i][correlator]
          correlator_data = get_data(correlator, data_files, replica_ensemble_name, ensemble.Nt, tsrc)
          corrs_to_average.append(correlator_data)

        averaged_corr_data = average_data(corrs_to_average)
        corrs_to_extend.append(averaged_corr_data)

      extended_corr_data = extend_data(corrs_to_extend)
      binfile = get_binfile(correlator, ensemble_name)
      if binfile not in data_to_write:
        data_to_write[binfile] = list()

      data_to_write[binfile].append(extended_corr_data)

      if corr_num % MAX_CORRS == 0:
        write_multiple_data(data_to_write, ensemble_name)

    write_multiple_data(data_to_write, ensemble_name)

def write_multiple_data(data_to_write, ensemble_name):
  ensemble_info = sigmond.MCEnsembleInfo(ensemble_name, 'ensembles.xml')
  bins_info = sigmond.MCBinsInfo(ensemble_info)
  bins_info.addOmissions(defs.omissions[ensemble_name])

  for bin_file, datas in data_to_write.items():
    bins_handler = sigmond.BinsPutHandler(bins_info, bin_file, sigmond.WriteMode.Protect, False)
    for data in datas:
      for obs_info, data_bins in data.items():
        bins_handler.putData(obs_info, sigmond.RVector(data_bins))

    bins_handler.close()
      
  data_to_write.clear()


def write_single_data(data, ensemble_name, bin_file):
  ensemble_info = sigmond.MCEnsembleInfo(ensemble_name, 'ensembles.xml')
  mcobs_xml = ET.Element("MCObservables")
  mcobs_xml_handler = sigmond.XMLHandler()
  mcobs_xml_handler.set_from_string(ET.tostring(mcobs_xml))
  sampling_info = sigmond.MCSamplingInfo()
  bins_info = sigmond.MCBinsInfo(ensemble_info)
  bins_info.addOmissions(defs.omissions[ensemble_name])
  obs_get_handler = sigmond.MCObsGetHandler(mcobs_xml_handler, bins_info, sampling_info)
  obs_handler = sigmond.MCObsHandler(obs_get_handler, False)

  for obs_info, data_bins in data.items():
    obs_handler.putBins(obs_info, sigmond.RVector(data_bins))
    xml_out = sigmond.XMLHandler("output", "")
    obs_handler.writeBinsToFile({obs_info}, bin_file, xml_out, sigmond.WriteMode.Protect)


def extend_data(data):
  extended_data = dict()

  #TODO: make sure all keys match between data

  obs_infos = data[0].keys()

  for obs_info in obs_infos:
    to_extend = [data_el[obs_info] for data_el in data]
    extend = np.concatenate(to_extend)
    extended_data[obs_info] = extend

  return extended_data


def average_data(data):
  averaged_data = dict()

  #TODO: make sure all keys match between data

  obs_infos = data[0].keys()

  for obs_info in obs_infos:
    to_average = [data_el[obs_info] for data_el in data]
    average = sum(to_average) / len(to_average)
    averaged_data[obs_info] = average

  return averaged_data


def get_data(correlator, data_files, ensemble_name, ensemble_Nt, tsrc):
  mcobs_xml = ET.Element("MCObservables")
  corr_data_xml = ET.SubElement(mcobs_xml, "BLCorrelatorData")

  bl_op = correlator.getSource().getBasicLapH()
  # TODO: replace with isFermionic()
  change_sign = bl_op.isBaryon() or bl_op.isMesonBaryon()

  file_list_infos = list()
  has_opposite = data_files[1] is not None
  for data_file in data_files:
    if data_file is None:
      continue
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

  keys = corr_handler.getOrderedKeys(correlator)
  tmin = keys[0].getTimeIndex()
  tmax = keys[-1].getTimeIndex()

  if has_opposite:
    correlator_opposite = sigmond.CorrelatorInfo(correlator.getSource(), correlator.getSink())
    keys = corr_handler.getOrderedKeys(correlator_opposite)
    tmin_opp = keys[0].getTimeIndex()
    tmax_opp = keys[-1].getTimeIndex()

    if tmin != tmin_opp:
      print("mismatch between tmin on corresponding correlators")
      exit()

    if tmax != tmax_opp:
      print("mismatch between tmax on corresponding correlators")
      exit()

  correlator_data = dict()
  corr_time_info = sigmond.CorrelatorAtTimeInfo(correlator, 0, False, False)
  if has_opposite:
    corr_time_info_opp = sigmond.CorrelatorAtTimeInfo(correlator_opposite, 0, False, False)

  corr_time_info_herm = sigmond.CorrelatorAtTimeInfo(correlator, 0, True, False)
  for tsep in range(tmin, tmax+1):
    for cmp_i, complex_arg in enumerate(COMPLEX_ARGS):
      if correlator.isSinkSourceSame() and complex_arg is sigmond.ComplexArg.ImaginaryPart:
        continue

      corr_time_info.resetTimeSeparation(tsep)
      corr_time_info_obs_info = sigmond.MCObsInfo(corr_time_info, complex_arg)
      data = np.array(obs_handler.getBins(corr_time_info_obs_info).array())

      if has_opposite:
        corr_time_info_opp.resetTimeSeparation(tsep)
        corr_time_info_opp_obs_info = sigmond.MCObsInfo(corr_time_info_opp, complex_arg)
        data_opp = np.array(obs_handler.getBins(corr_time_info_opp_obs_info).array())
        data = 0.5*(data + np.conj(data_opp))

      if change_sign:
        for cfg_ind, config in enumerate(defs.configs[ensemble_name]):
          T = tsep + tsrc + srcTList[config]
          if T >= ensemble_Nt:
            data[cfg_ind] = -data[cfg_ind]

      corr_time_info_herm.resetTimeSeparation(tsep)
      corr_time_info_herm_obs_info = sigmond.MCObsInfo(corr_time_info_herm, complex_arg)
      correlator_data[corr_time_info_herm_obs_info] = data

  return correlator_data


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
    corr_opposite = sigmond.CorrelatorInfo(corr.getSource(), corr.getSink())
    if corr_opposite in corr_files:
      corr_files[corr_opposite][1] = corr.getFileName(corr)
    else:
      corr_files[corr] = (corr_handler.getFileName(corr), None)

  return corr_files


def get_binfile(correlator, ensemble_name):
  operator = correlator.getSink().getBasicLapH()
  P = (operator.getXMomentum(), operator.getYMomentum(), operator.getZMomentum())
  irrep = operator.getLGIrrep()
  irrep_row = operator.getLGIrrepRow()
  isospin = operator.getIsospin()
  strangeness = operator.getStrangeness()
  channel = defs.Channel(P, irrep, irrep_row, isospin, strangeness)

  output_dir = os.path.join(defs.output_dir, ensemble_name)
  os.makedirs(output_dir, exist_ok=True)

  bin_file = os.path.join(output_dir, f"{channel.iso_strange_str()}.bin")
  return bin_file


if __name__ == "__main__":
  main()
