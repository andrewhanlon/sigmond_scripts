from typing import NamedTuple
import copy
import itertools

import logging

from data_handling.data_files import DataFiles, FileInfo
from operator_info.operator_set import OperatorSet
from operator_info.operator import Operator

import sigmond

class CorrelatorDataInfo(NamedTuple):
  tmin: int
  tmax: int
  fileinfo: FileInfo

  def getUpdatedTsepRange(self, tmin, tmax):
    if self.tmin < tmin:
      tmin = self.tmin

    if self.tmax > tmax:
      tmax = self.tmax

    return CorrelatorDataInfo(tmin, tmax, self.fileinfo)


class CorrelatorData:

  def __init__(self):
    self._data_files = dict()
    self._correlators = dict()
    self._vevs = dict()
    self._operator_set = OperatorSet()

  def addCorrelator(self, correlator_info, tmin, tmax, fileinfo):
    op_snk = Operator(correlator_info.getSink())
    op_src = Operator(correlator_info.getSource())
    self._operator_set.addOperator(op_snk)
    self._operator_set.addOperator(op_src)

    if correlator_info not in self._correlators:
      self._correlators[correlator_info] = CorrelatorDataInfo(tmin, tmax, fileinfo)
      self._addFileInfo(fileinfo, op_snk.channel)
    elif self._correlators[correlator_info].fileinfo != fileinfo:
      logging.warning(f"Correlator {correlator_info} already added with different file, skipping...")
    else:
      self._correlators[correlator_info] = self._correlators[correlator_info].getUpdatedTsepRange(tmin, tmax)

  def addVEV(self, vev, fileinfo):
    op = Operator(vev)
    self._operator_set.addOperator(op)

    if vev not in self._vevs:
      self._vevs[vev] = fileinfo
      self._addFileInfo(fileinfo, op.channel)
    elif fileinfo != self._vevs[vev]:
      logging.warning(f"VEV {vev} already added with different file, skipping...")

  @property
  def channels(self):
    return self._operator_set.channels

  @property
  def channel(self):
    if len(self.channels) != 1:
      logging.error("Requires one and only one channel")
    return self.channels.pop()

  @property
  def operators(self):
    return self._operator_set.operators

  @property
  def operator_set(self):
    return self._operator_set

  def getChannelOperators(self, channel):
    return self._operator_set.getOperators(channel)

  def getChannelDataFiles(self, channel):
    return self._data_files[channel]

  def hasCorrelator(self, correlator):
    return correlator in self._correlators

  def getCorrInfo(self, corr):
    if self.hasCorrelator(corr):
      return self._correlators[corr]
    else:
      return None

  def getChannelsLargestTRange(self, max_time, *channels):
    tmax = 0
    tmin = max_time
    found = False
    for corr, info in self._correlators.items():
      channel = Operator(corr.getSource()).channel
      if channel not in channels:
        continue

      found = True

      if info.tmin < tmin:
        tmin = info.tmin
      if info.tmax > tmax:
        tmax = info.tmax

    if found:
      return (tmin, tmax)
    else:
      return (-1, -1)

  def getOperatorSetSmallestTRange(self, max_time, operator_set):
    tmax = max_time
    tmin = 0
    for two_ops in itertools.combinations_with_replacement(operator_set.operators, 2):
      corr = sigmond.CorrelatorInfo(two_ops[0].operator_info, two_ops[1].operator_info)
      corr_rev = sigmond.CorrelatorInfo(two_ops[1].operator_info, two_ops[0].operator_info)

      corr_info = self.getCorrInfo(corr)
      corr_info_rev = self.getCorrInfo(corr_rev)

      if corr_info is not None and corr_info_rev is not None:
        tmin_try = min(corr_info.tmin, corr_info_rev.tmin)
        tmax_try = max(corr_info.tmax, corr_info_rev.tmax)

      elif corr_info is not None:
        tmin_try = corr_info.tmin
        tmax_try = corr_info.tmax

      elif corr_info_rev is not None:
        tmin_try = corr_info_rev.tmin
        tmax_try = corr_info_rev.tmax

      else:
        return (-1, -1)

      if tmin_try > tmin:
        tmin = tmin_try
      if tmax_try < tmax:
        tmax = tmax_try

    return (tmin, tmax)

  def _addFileInfo(self, fileinfo, channel):
    if channel not in self._data_files:
      self._data_files[channel] = DataFiles()

    if fileinfo.filetype is sigmond.FileType.Correlator:
      self._data_files[channel].addCorrelatorFiles(fileinfo.filename)
    elif fileinfo.filetype is sigmond.FileType.VEV:
      self._data_files[channel].addVEVFiles(fileinfo.filename)
    elif fileinfo.filetype is sigmond.FileType.Bins:
      self._data_files[channel].addBinFiles(fileinfo.filename)
    elif fileinfo.filetype is sigmond.FileType.Samplings:
      self._data_files[channel].addSamplingFiles(fileinfo.filename)
    else:
      logging.error(f"Unrecognized filetype {fileinfo.filetype}")

  def __add__(self, other):
    if isinstance(other, self.__class__):
      new_instance = self.__class__()

      for channel in list(self._data_files.keys()) + list(other._data_files.keys()):
        new_instance._data_files[channel] = self._data_files.get(channel, DataFiles()) + other._data_files.get(channel, DataFiles())

      new_instance._correlators = self._correlators.copy()
      new_instance._vevs = self._vevs.copy()
      for op in self._operator_set.operators:
        new_instance._operator_set.addOperator(op)

      for corr_info, corr_data_info in other._correlators.items():
        new_instance.addCorrelator(corr_info, **corr_data_info._asdict())

      for vev, fileinfo in other._vevs.items():
        new_instance.addVEV(vev, fileinfo)

      return new_instance

    return NotImplemented

