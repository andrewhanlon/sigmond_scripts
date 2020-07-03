#!/usr/bin/env python

import os

import numpy as np
import xml.etree.ElementTree as ET

import tqdm

import defs

import sigmond

REL_TOL=1e-17
ABS_TOL=1e-17

ensemble_name = "cls21_c103_r005"
new_data_dir = "/media/ext2/research/data/raw_C103/r005/src0/"
old_data_dir = "/media/ext2/research/data/raw_C103/r005.old/src0/"

def main():

  new_corr_files = get_corr_files(ensemble_name, new_data_dir)
  old_corr_files = get_corr_files(ensemble_name, old_data_dir)

  success = 0
  fails = 0

  for corr, corr_file in new_corr_files.items():
    if corr not in old_corr_files:
      continue

    new_data = get_data(corr, corr_file, ensemble_name)
    old_data = get_data(corr, old_corr_files[corr], ensemble_name)

    if check_data(new_data, old_data):
      success += 1
    else:
      print(f"FAILURE: {corr.corr_str()}")
      fails += 1

  print(f"{fails} FAILS out of {fails+success}")

def check_data(new_data, old_data):
  if not np.allclose(new_data, old_data, rtol=REL_TOL, atol=ABS_TOL):
    for new_val, old_val in zip(new_data, old_data):
      if not np.isclose(new_val, old_val, rtol=REL_TOL, atol=ABS_TOL):
        print(f"{new_val} != {old_val}")

    return False

  else:
    return True



def get_data(correlator, data_file, ensemble_name):
  mcobs_xml = ET.Element("MCObservables")
  corr_data_xml = ET.SubElement(mcobs_xml, "BLCorrelatorData")

  stub, ext = os.path.splitext(data_file)
  suffix = int(ext[1:])
  file_list_info = sigmond.FileListInfo(stub, suffix, suffix, False)
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

  data = np.zeros(tmax-tmin+1, dtype=np.complex128)

  corr_time_info = sigmond.CorrelatorAtTimeInfo(correlator, 0, False, False)
  for tsep in range(tmin, tmax+1):
    corr_time_info.resetTimeSeparation(tsep)
    obs_info_re = sigmond.MCObsInfo(corr_time_info, sigmond.ComplexArg.RealPart)
    obs_info_im = sigmond.MCObsInfo(corr_time_info, sigmond.ComplexArg.ImaginaryPart)

    re_val = obs_handler.getBin(obs_info_re, 0)
    im_val = obs_handler.getBin(obs_info_im, 0)

    data[tsep-tmin] = re_val + im_val*1j

  return data


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



if __name__ == "__main__":
  main()
