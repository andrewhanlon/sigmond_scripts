import os

import xml.etree.ElementTree as ET
from typing import NamedTuple

import sigmond
import sigmond_scripts.analysis.utils.util as util


class FileInfo(NamedTuple):
  filename: str
  filetype: sigmond.FileType


class DataFiles:

  def __init__(self, checksums=False):
    self._bl_corr_files = dict()
    self._bl_vev_files = dict()
    self._bin_files = set()
    self._sampling_files = set()

    self._checksums = checksums

  def addCorrelatorFiles(self, *filenames):
    file_list_infos = list()
    for filename in filenames:
      base, ext = os.path.splitext(filename)
      suffix = int(ext[1:])
      file_list_info = sigmond.FileListInfo(base, suffix, suffix, False)
      file_list_infos.append(file_list_info)

    self.addCorrelators(*file_list_infos)
    
  def addVEVFiles(self, *filenames):
    file_list_infos = list()
    for filename in filenames:
      base, ext = os.path.splitext(filename)
      suffix = int(ext[1:])
      file_list_info = sigmond.FileListInfo(base, suffix, suffix, False)
      file_list_infos.append(file_list_info)

    self.addVEVs(*file_list_infos)

  def addCorrelators(self, *file_lists):
    for file_list in file_lists:
      file_stub = file_list.getFileStub()
      if file_stub in self._bl_corr_files:
        min_file = self._bl_corr_files[file_stub].getMinFileNumber()
        max_file = self._bl_corr_files[file_stub].getMaxFileNumber()
        old_min_file = min_file
        old_max_file = max_file
        if file_list.getMinFileNumber() < min_file:
          min_file = file_list.getMinFileNumber()
        if file_list.getMaxFileNumber() > max_file:
          max_file = file_list.getMaxFileNumber()
        
        min_file = 0

        overwrite = self._bl_corr_files[file_stub].isModeOverwrite() and file_list.isModeOverwrite()
        new_file_list = sigmond.FileListInfo(file_stub, min_file, max_file, overwrite)
        if new_file_list != self._bl_corr_files[file_stub] and (min_file < old_min_file or max_file > old_max_file):
          self._bl_corr_files[file_stub] = new_file_list

      else:
        self._bl_corr_files[file_stub] = file_list

  def addVEVs(self, *file_lists):
    for file_list in file_lists:
      file_stub = file_list.getFileStub()
      if file_stub in self._bl_vev_files:

        min_file = self._bl_vev_files[file_stub].getMinFileNumber()
        max_file = self._bl_vev_files[file_stub].getMaxFileNumber()
        old_min_file = min_file
        old_max_file = max_file
        if file_list.getMinFileNumber() < min_file:
          min_file = file_list.getMinFileNumber()
        if file_list.getMaxFileNumber() > max_file:
          max_file = file_list.getMaxFileNumber()

        overwrite = self._bl_vev_files[file_stub].isModeOverwrite() and file_list.isModeOverwrite()
        new_file_list = sigmond.FileListInfo(file_stub, min_file, max_file, overwrite)

        if new_file_list != self._bl_vev_files[file_stub] and (min_file < old_min_file or max_file > old_max_file):
          self._bl_vev_files[file_stub] = new_file_list

      else:
        self._bl_vev_files[file_stub] = file_list

  def addBinFiles(self, *filenames):
    for filename in filenames:
      self._bin_files.add(filename)

  def addSamplingFiles(self, *filenames):
    for filename in filenames:
      self._sampling_files.add(filename)

  @property
  def bl_corr_files(self):
    return [file_list_info for file_list_info in self._bl_corr_files.values()]

  @property
  def bl_vev_files(self):
    return [file_list_info for file_list_info in self._bl_vev_files.values()]
  
  @property
  def bin_files(self):
    return self._bin_files

  @property
  def sampling_files(self):
    return self._sampling_files

  @property
  def checksums(self):
    return self._checksums

  @property
  def empty(self):
    return not (self.bl_corr_files or self.bl_vev_files or self.bin_files or self.sampling_files)

  def xml(self):
    _xml = ET.Element("MCObservables")
    if self.checksums:
      ET.SubElement(_xml, "UseCheckSums")

    if self.bl_corr_files:
      corr_tag = ET.SubElement(_xml, "BLCorrelatorData")
      for corr_file in self.bl_corr_files:
        corr_tag.append(corr_file.xml())

    if self.bl_vev_files:
      vev_tag = ET.SubElement(_xml, "BLVEVData")
      for vev_file in self.bl_vev_files:
        vev_tag.append(vev_file.xml())

    if self.bin_files:
      bin_tag = ET.SubElement(_xml, "BinData")
      for bin_file in self.bin_files:
        ET.SubElement(bin_tag, "FileName").text = bin_file

    if self.sampling_files:
      sampling_tag = ET.SubElement(_xml, "SamplingData")
      for sampling_file in self.sampling_files:
        ET.SubElement(sampling_tag, "FileName").text = sampling_file

    return _xml

  def __str__(self):
    return util.xmltostr(self.xml())

  def __add__(self, other):
    """
    Note: order matters due to the insertion order being preserved in dicts used for
          the basic laph files. 'self' will have its files ordered before 'other'.
    """
    if isinstance(other, self.__class__):
      checksums = other.checksums or self.checksums
      new_data_files = self.__class__(checksums)
      new_data_files._bl_corr_files = self._bl_corr_files.copy()
      new_data_files.addCorrelators(*other.bl_corr_files)
      new_data_files._bl_vev_files = self._bl_vev_files.copy()
      new_data_files.addVEVs(*other.bl_vev_files)

      new_data_files._bin_files = self.bin_files | other.bin_files
      new_data_files._sampling_files = self.sampling_files | other.sampling_files

      return new_data_files

    return NotImplemented
