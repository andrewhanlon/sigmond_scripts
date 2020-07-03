#!/usr/bin/env python

import os

import h5py
import numpy as np
import xml.etree.ElementTree as ET
from sortedcontainers import SortedSet

import sigmond

import defs

ref_elements = {
    0: defs.EquivalentFrame((0,0,0), 1),
    1: defs.EquivalentFrame((0,0,1), 1),
    2: defs.EquivalentFrame((0,1,1), 1),
    3: defs.EquivalentFrame((1,1,1), 1),
    4: defs.EquivalentFrame((0,0,2), 1),
    5: defs.EquivalentFrame((0,1,2), 1),
    6: defs.EquivalentFrame((1,1,2), 1),
    8: defs.EquivalentFrame((0,2,2), 1),
}

base_data_dir = "/media/ext2/research/data/raw_C103/r005/src0/"

def main():
  for ensemble in defs.ensembles:
    for replica in ensemble.replica:
      #search_dir = os.path.join(base_data_dir, channel)
      search_dir = base_data_dir
      ensemble_name = f"{ensemble.name}_{replica}"
      channel_data = find_data(ensemble_name, search_dir)
      channel_data = average_data(channel_data)


def find_data(ensemble_name, search_dir):
  ensemble_info = sigmond.MCEnsembleInfo(ensemble_name, 'ensembles.xml')

  # search for data files
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
    isospin = op_src_bl.getIsospin()
    strangeness = op_src_bl.getStrangeness()
    averaged_channel = defs.AveragedChannel(Psq, irrep, isospin, strangeness)
    equivalent_frame = defs.EquivalentFrame(P, irrep_row)
    if averaged_channel not in operators:
      operators[averaged_channel] = dict()
    if equivalent_frame not in operators[averaged_channel]:
      operators[averaged_channel][equivalent_frame] = SortedSet()

    operators[averaged_channel][equivalent_frame].add(op_snk)
    operators[averaged_channel][equivalent_frame].add(op_src)

  return operators


def average_data(data):
  for averaged_channel, equivalent_frames in data.items():
    ref_channel = ref_elements[averaged_channel.psq]
    full_channel = defs.Channel(ref_channel.mom, averaged_channel.irrep, ref_channel.row, averaged_channel.isospin, averaged_channel.strangeness)
    ref_operators = equivalent_frames[ref_channel]
    ref_ops_map = _getOperatorsMap(ref_operators, averaged_channel)
    result_operators = sorted(list(ref_ops_map.keys()))
    mismatch = False
    for equivalent_frame, operators in equivalent_frames.items():
      ops_map = _getOperatorsMap(operators, averaged_channel)
      if sorted(list(ops_map.keys())) != result_operators:

        if not mismatch:
          print(f"Mismatch found in {full_channel!s}")
          print("  All reference frame operators:")
          for averaged_op, raw_op in sorted(ref_ops_map.items()):
            print(f"    {raw_op}")
          print()

          mismatch = True

    if mismatch:
      for equivalent_frame, operators in equivalent_frames.items():
        ops_map = _getOperatorsMap(operators, averaged_channel)
        print("  Mismatched operators in equivalent frame:")
        for averaged_op, raw_op in sorted(ops_map.items()):
          if averaged_op in ref_ops_map:
            print(f"    {raw_op}")
          else:
            print(f"  * {raw_op}")
        print()
      
def _getOperatorsMap(operators, averaged_channel):
  op_map = dict()
  for operator in operators:
    averaged_op = _getAveragedOperator(operator, averaged_channel)
    if averaged_op in op_map:
      print(f"Conflicting operators {operator} and {op_map[averaged_op]}")
      print()
      print(averaged_op)
      print()
      for op in operators:
        print(op)
      exit()

    op_map[averaged_op] = operator

  return op_map

def _getAveragedOperator(operator, averaged_channel):
  op_info = operator.getBasicLapH()
  if op_info.getNumberOfHadrons() == 1:
    obs_name = f"{NAME_MAP[op_info.getFlavor()]}-{op_info.getHadronSpatialType(1)}_{op_info.getHadronSpatialIdNumber(1)}"
    obs_id = 0
  else:
    obs_name = ""
    for had_num in range(1, op_info.getNumberOfHadrons()+1):
      had_name = NAME_MAP[op_info.getHadronFlavor(had_num)]
      had_psq = op_info.getHadronXMomentum(had_num)**2 + op_info.getHadronYMomentum(had_num)**2 + op_info.getHadronZMomentum(had_num)**2
      spat_type = op_info.getHadronSpatialType(had_num)
      spat_id = op_info.getHadronSpatialIdNumber(had_num)
      irrep = op_info.getHadronLGIrrep(had_num)
      obs_name += f"{had_name}({had_psq}_{spat_type}_{spat_id}_{irrep})"

    obs_id = op_info.getLGClebschGordonIdNum()

  return f"{obs_name} {obs_id}"
  

# TODO: use flavor_map in operator_info/operator.py
NAME_MAP = {
    'pion': 'pi',
    'eta': 'e',
    'phi': 'p',
    'kaon': 'k',
    'kbar': 'kb',
    'nucleon': 'N',
    'delta': 'D',
    'sigma': 'S',
    'lambda': 'L',
    'xi': 'X',
    'omega': 'O',
}


if __name__ == "__main__":
  main()
