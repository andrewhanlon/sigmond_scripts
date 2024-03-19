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
    print(f"Ensemble: {ensemble_name}")

    for flavor_channel in ensemble.flavor_channels:
      print(f"Flavor: {flavor_channel}")

      corr_files = dict()
      correlators = list()
      for replica in ensemble.replica:
        replica_ensemble_name = f"{ensemble_name}_{replica}"
        corr_files[replica] = dict()
        for tsrc in ensemble.sources:
          data_dir = os.path.join(ensemble_name, replica, f"{flavor_channel}_t0{tsrc[0]}", tsrc[1])
          print(f"Searching {data_dir}...", end='', flush=True)
          search_dir = os.path.join(defs.base_data_dir, data_dir)
          corr_files[replica][tsrc] = get_corr_files(replica_ensemble_name, search_dir)
          correlators.append(corr_files[replica][tsrc].keys())
          print("done")

      # check all elements of correlators have same keys
      all_equal = all(corrs==correlators[0] for corrs in correlators)
      if not all_equal:
        print("not all equal\n\n")
        for corrs_num, corrs in enumerate(correlators):
          print(f"corrs num {corrs_num}")
          print("--------------------")
          for corr in sorted(corrs):
            print(f"\t{corr}")

          print()

      correlators = correlators[0]

      for correlator in tqdm.tqdm(correlators):
        corrs_to_extend = list()
        for replica in ensemble.replica:
          replica_ensemble_name = f"{ensemble_name}_{replica}"
          corrs_to_average = list()
          for tsrc in ensemble.sources:
            data_file, data_file_opposite, is_backwards = corr_files[replica][tsrc][correlator]
            correlator_data = get_data(correlator, data_file, data_file_opposite, is_backwards, replica_ensemble_name)
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
    obs_handler.writeBinsToFile({obs_info}, bin_file, xml_out, sigmond.WriteMode.Protect, 'D')


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


def get_data(correlator, data_file, data_file_opposite, is_backwards, ensemble_name):
  if is_backwards:
    correlator.setBackwards()

  mcobs_xml = ET.Element("MCObservables")
  corr_data_xml = ET.SubElement(mcobs_xml, "BLCorrelatorData")

  file_list_infos = list()
  has_opposite = data_file_opposite is not None
  if has_opposite:
    data_files = [data_file, data_file_opposite]
  else:
    data_files = [data_file]

  for d_file in data_files:
    stub, ext = os.path.splitext(d_file)
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
    for complex_arg in COMPLEX_ARGS:
      if correlator.isSinkSourceSame() and complex_arg is sigmond.ComplexArg.ImaginaryPart:
        continue

      sign = 1. if complex_arg is sigmond.ComplexArg.RealPart else -1.

      corr_time_info.resetTimeSeparation(tsep)
      corr_time_info_obs_info = sigmond.MCObsInfo(corr_time_info, complex_arg)
      data = np.array(obs_handler.getBins(corr_time_info_obs_info).array())

      if has_opposite:
        print('has_opposite')
        sys.exit("has opposite")
        corr_time_info_opp.resetTimeSeparation(tsep)
        corr_time_info_opp_obs_info = sigmond.MCObsInfo(corr_time_info_opp, complex_arg)
        data_opp = sign*np.array(obs_handler.getBins(corr_time_info_opp_obs_info).array())
        data = 0.5*(data + data_opp)

      corr_time_info_herm.resetTimeSeparation(tsep)
      if is_backwards:
        corr_time_info_herm.setForwards()

      corr_time_info_herm_obs_info = sigmond.MCObsInfo(corr_time_info_herm, complex_arg)
      correlator_data[corr_time_info_herm_obs_info] = data

  if is_backwards:
    correlator.setForwards()

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
    corr_file_name = corr_handler.getFileName(corr)
    is_backwards = corr.isBackwards()
    if is_backwards:
      corr.setForwards()

    corr_opposite = sigmond.CorrelatorInfo(corr.getSource(), corr.getSink())
    if corr_opposite in corr_files:
      #corr_files[corr_opposite][1] = corr_file_name
      corr_files[corr_opposite] = (corr_files[corr_opposite][0], corr_file_name, corr_files[corr_opposite][2])
    else:
      corr_files[corr] = (corr_file_name, None, is_backwards)

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
