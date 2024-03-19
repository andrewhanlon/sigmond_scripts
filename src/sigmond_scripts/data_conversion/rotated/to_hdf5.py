#!/usr/bin/env python

import os

import h5py
import numpy as np
import xml.etree.ElementTree as ET
from sortedcontainers import SortedSet

import sigmond

import defs

rel_rot_data_dir = ".sigmond/data/rotated_correlators"

ensembles = ["cls21_c103"]

BINS = True

def main():
  for ensemble in ensembles:
    print(f"starting for ensemble {ensemble}")
    rot_data_dir = os.path.join(defs.input_dirs[ensemble], rel_rot_data_dir)
    if not os.path.isdir(rot_data_dir):
      continue

    rot_data = find_data(ensemble, rot_data_dir)
    write_to_file(rot_data, ensemble)
    print("done")


def write_to_file(the_data, ensemble_name):
  os.makedirs(defs.output_dir, exist_ok=True)
  filename = os.path.join(defs.output_dir, f"rotated_{ensemble_name}.hdf5")
  h5_file = h5py.File(filename, 'w')

  for channel, data in the_data.items():
    h5_file.create_dataset(repr(channel), data=data)

  h5_file.close()

def find_data(ensemble_name, search_dir):
  ensemble_info = sigmond.MCEnsembleInfo(ensemble_name, 'ensembles.xml')

  # search for data files
  bin_files = set()
  smp_files = set()
  print(f"Searching for correlators and time separations in {search_dir}")
  for root, dirs, files in os.walk(search_dir, topdown=True):
    for filename in files:
      full_filename = os.path.join(root, filename)
      try:
        file_type = sigmond.getFileID(full_filename)
      except ValueError:
        continue
      
      if file_type == sigmond.FileType.Bins:
        bin_files.add(full_filename)
      elif file_type == sigmond.FileType.Samplings:
        smp_files.add(full_filename)
      else:
        continue

  if not BINS and bin_files:
    print("ignoring bin files found")
  elif BINS and smp_files:
    print("ignoring smp files found")
      
  # create sigmond handlers
  mcobs_xml = ET.Element("MCObservables")
  if BINS:
    bins_data_xml = ET.SubElement(mcobs_xml, "BinData")
    for bin_file in bin_files:
      ET.SubElement(bins_data_xml, "FileName").text = bin_file
  else:
    smp_data_xml = ET.SubElement(mcobs_xml, "SamplingData")
    for smp_file in smp_files:
      ET.SubElement(smp_data_xml, "FileName").text = smp_file

  mcobs_xml_handler = sigmond.XMLHandler()
  mcobs_xml_handler.set_from_string(ET.tostring(mcobs_xml))
  sampling_info = sigmond.MCSamplingInfo()
  bins_info = sigmond.MCBinsInfo(ensemble_info)
  bins_info.addOmissions(defs.omissions[ensemble_name])
  obs_get_handler = sigmond.MCObsGetHandler(mcobs_xml_handler, bins_info, sampling_info)
  obs_handler = sigmond.MCObsHandler(obs_get_handler, False)

  if BINS:
    data_handler = sigmond.BinsGetHandler(bins_info, bin_files)
    data_size = bins_info.getNumberOfBins()
  else:
    data_handler = sigmond.SamplingsGetHandler(bins_info, sampling_info, smp_files)
    data_size = sampling_info.getNumberOfReSamplings(bins_info) + 1

  obs_infos = data_handler.getKeys()

  data_dict = dict()
  max_level = dict()
  max_time = dict()
  data_size

  for obs_info in obs_infos:
    if obs_info.isImaginaryPart():
      continue
    operator = obs_info.getCorrelatorSourceInfo().getGenIrrep()
    channel = get_channel(operator)
    if channel not in data_dict:
      data_dict[channel] = dict()
    if channel not in max_level:
      max_level[channel] = 0
    if channel not in max_time:
      max_time[channel] = 0

    level = operator.getIDIndex()
    if level not in data_dict[channel]:
      data_dict[channel][level] = dict()

    max_level[channel] = max(max_level[channel], level)

    time = obs_info.getCorrelatorTimeIndex()
    #print(f"{channel!r} - {level} - {time}")
    if time in data_dict[channel][level]:
      print("another problem")
      exit()
    max_time[channel] = max(max_time[channel], time)

    if BINS:
      data = obs_handler.getBins(obs_info).array()
    else:
      data = obs_handler.getFullAndSamplingValues(obs_info, sampling_info.getSamplingMode()).array()

    data_dict[channel][level][time] = data

  final_data_dict = dict()
  for channel, levels in data_dict.items():
    data_arr = np.zeros((data_size, max_level[channel]+1, max_time[channel]+1), dtype=np.float64)
    for level, times in levels.items():
      for time, data in times.items():
        data_arr[:,level,time] = data

    final_data_dict[channel] = data_arr

  return final_data_dict


def get_channel(operator):
  PSQ = operator.getMomentumSquared()
  irrep = operator.getLGIrrep()
  isospin = operator.getIsospin()
  strangeness = operator.getStrangeness()
  channel = defs.Channel(PSQ, irrep, isospin, strangeness)
  return channel


if __name__ == "__main__":
  main()
