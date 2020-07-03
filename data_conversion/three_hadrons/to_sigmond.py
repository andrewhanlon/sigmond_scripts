#!/usr/bin/env python

import os

import numpy as np
import xml.etree.ElementTree as ET

import tqdm

import defs

import sigmond

COMPLEX_ARGS = [sigmond.ComplexArg.RealPart, sigmond.ComplexArg.ImaginaryPart]

def main():
  for ensemble in defs.ensembles:
    ensemble_name = ensemble.name
    for flavor_channel in ensemble.flavor_channels:

      corr_files = [[0 for replica in ensemble.replica] for tsrc in ensemble.sources]
      for replica_i, replica in enumerate(ensemble.replica):
        replica_ensemble_name = f"{ensemble_name}_{replica}"
        for tsrc_i, tsrc in enumerate(ensemble.sources):
          data_dir = os.path.join(ensemble_name, replica, f"{flavor_channel}_t0{tsrc}")
          search_dir = os.path.join(defs.base_data_dir, data_dir)
          corr_files[tsrc_i][replica_i] = get_corr_files(replica_ensemble_name, search_dir)

      # TODO: check all elements of correlators have same keys
      correlators = corr_files[0][0].keys()

      for correlator in tqdm.tqdm(correlators):
        corrs_to_extend = list()
        for replica_i, replica in enumerate(ensemble.replica):
          replica_ensemble_name = f"{ensemble_name}_{replica}"
          corrs_to_average = list()
          for tsrc_i, tsrc in enumerate(ensemble.sources):
            data_file = corr_files[tsrc_i][replica_i][correlator]
            correlator_data = get_data(correlator, data_file, replica_ensemble_name)
            corrs_to_average.append(correlator_data)

          averaged_corr_data = average_data(corrs_to_average)
          corrs_to_extend.append(averaged_corr_data)

        extended_corr_data = extend_data(corrs_to_extend)
        channel = get_channel(correlator)
        write_data(extended_corr_data, ensemble_name, channel, flavor_channel)


def write_data(data, ensemble_name, channel, flavor_channel):
  output_dir = os.path.join(defs.output_dir, ensemble_name, flavor_channel)
  os.makedirs(output_dir, exist_ok=True)

  ensemble_info = sigmond.MCEnsembleInfo(ensemble_name, 'ensembles.xml')
  mcobs_xml = ET.Element("MCObservables")
  mcobs_xml_handler = sigmond.XMLHandler()
  mcobs_xml_handler.set_from_string(ET.tostring(mcobs_xml))
  sampling_info = sigmond.MCSamplingInfo()
  bins_info = sigmond.MCBinsInfo(ensemble_info)
  bins_info.addOmissions(defs.omissions[ensemble_name])
  obs_get_handler = sigmond.MCObsGetHandler(mcobs_xml_handler, bins_info, sampling_info)
  obs_handler = sigmond.MCObsHandler(obs_get_handler, False)

  bin_file = os.path.join(output_dir, f"{channel!r}.bin")
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


def get_data(correlator, data_file, ensemble_name):
  stub, ext = os.path.splitext(data_file)
  suffix = int(ext[1:])
  file_list_info = sigmond.FileListInfo(stub, suffix, suffix, False)

  mcobs_xml = ET.Element("MCObservables")
  corr_data_xml = ET.SubElement(mcobs_xml, "BLCorrelatorData")
  corr_data_xml.append(file_list_info.xml())

  mcobs_xml_handler = sigmond.XMLHandler()
  mcobs_xml_handler.set_from_string(ET.tostring(mcobs_xml))
  sampling_info = sigmond.MCSamplingInfo()
  ensemble_info = sigmond.MCEnsembleInfo(ensemble_name, 'ensembles.xml')
  bins_info = sigmond.MCBinsInfo(ensemble_info)
  bins_info.addOmissions(defs.omissions[ensemble_name])
  obs_get_handler = sigmond.MCObsGetHandler(mcobs_xml_handler, bins_info, sampling_info)
  obs_handler = sigmond.MCObsHandler(obs_get_handler, False)

  corr_handler = sigmond.BLCorrelatorDataHandler([file_list_info], set(), set(), ensemble_info)

  keys = corr_handler.getOrderedKeys(correlator)
  tmin = keys[0].getTimeIndex()
  tmax = keys[-1].getTimeIndex()

  correlator_data = dict()
  corr_time_info = sigmond.CorrelatorAtTimeInfo(correlator, 0, True, False)
  for tsep in range(tmin, tmax+1):
    for cmp_i, complex_arg in enumerate(COMPLEX_ARGS):
      corr_time_info.resetTimeSeparation(tsep)
      obs_info = sigmond.MCObsInfo(corr_time_info, complex_arg)
      data = obs_handler.getBins(obs_info)
      correlator_data[obs_info] = np.array(data.array())

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
    corr_files[corr] = corr_handler.getFileName(corr)

  return corr_files


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
