import os
import h5py

from typing import NamedTuple
import logging
from sortedcontainers import SortedSet
import tqdm

import sigmond_scripts.analysis.utils.util as util
from sigmond_scripts.analysis.data_handling.data_files import DataFiles, FileInfo
from sigmond_scripts.analysis.data_handling.correlator_data import CorrelatorData
from sigmond_scripts.analysis.operator_info.operator import Operator

import sigmond

class DataHandler(metaclass=util.Singleton):
  """ The all important Data Handler!  """

  rel_averaged_datadir = [""]
  rel_rotated_datadir = [""]
  # rel_rotated_pivotdir = ""

  def __init__(self, project_info):
    self.project_info = project_info
    
    self._raw_data = CorrelatorData()
    self._averaged_data = CorrelatorData()
    self._rotated_data = CorrelatorData()

    self.raw_data_files = DataFiles()
    self.averaged_data_files = DataFiles()
    self.rotated_data_files = DataFiles()

    self._findRawData()
    if self.rel_averaged_datadir[0]:
      self.findAveragedData()
    if self.rel_rotated_datadir[0]:
      self.findRotatedData()

  @property
  def project_dir(self):
    return self.project_info.project_dir

  @property
  def raw_data_dirs(self):
    return self.project_info.raw_data_dirs

  @property
  def bins_info(self):
    return self.project_info.bins_info

  @property
  def sampling_info(self):
    return self.project_info.sampling_info

  @property
  def ensemble_info(self):
    return self.bins_info.getMCEnsembleInfo()

  @property
  def data_files(self):
    return self.project_info.data_files
  
  # def raw_data_files(self):
  #   return self.raw_data_files
  
  @property
  def precompute(self):
    return self.project_info.precompute
  
  @property
  def raw_channels(self):
    return self._raw_data.channels

  @property
  def averaged_channels(self):
    return self._averaged_data.channels

  @property
  def rotated_channels(self):
    return self._rotated_data.channels
  
  @property
  def averaged_datadir(self):
    datadir = os.path.join(self.project_dir, self.rel_averaged_datadir[0])
    # os.makedirs(datadir, exist_ok=True)
    return datadir
  
  @property
  def averaged_datadirs(self):
    datadirs = [os.path.join(self.project_dir, datadir) for datadir in self.rel_averaged_datadir]
    # os.makedirs(datadir, exist_ok=True)
    return datadirs

  @property
  def rotated_datadir(self):
    datadir = os.path.join(self.project_dir, self.rel_rotated_datadir[0])
    # os.makedirs(datadir, exist_ok=True)
    return datadir

  @property
  def rotated_datadirs(self):
    datadirs = [os.path.join(self.project_dir, datadir) for datadir in self.rel_rotated_datadir]
    return datadirs

  def rotated_datafile(self, rotated_operator_set):
    datafile = f"{rotated_operator_set!r}.dat"
    return os.path.join(self.rotated_datadir, datafile)

  @property
  def pivot_datadir(self):
    datadir = os.path.join(self.project_dir, self.rel_rotated_pivotdir)
    os.makedirs(datadir, exist_ok=True)
    return datadir

  def pivotfile(self, rotated_operator_set):
    datafile = f"{rotated_operator_set!r}.piv"
    return os.path.join(self.pivot_datadir, datafile)

  # def getRotatedDataFiles(self, rotated_basis=None):
  #   corr_data = self._readRotatedData(rotated_basis)
  #   return corr_data.getChannelDataFiles(corr_data.channel)

  # def getRotatedOperators(self, rotated_basis=None):
  #   corr_data = self._readRotatedData(rotated_basis)
  #   return corr_data.getChannelOperators(corr_data.channel)
  
  # def getRotatedTRange(self, rotated_basis=None):
  #   corr_data = self._readRotatedData(rotated_basis)
  #   time_extent = self.ensemble_info.getLatticeTimeExtent()
  #   return corr_data.getOperatorSetSmallestTRange(time_extent, corr_data.operator_set)
  
  # def getRotatedChannels(self):
  #   corr_data = self._readRotatedData()
  #   return corr_data.keys()


  def getChannelOperators(self, channel):
    if channel in self.averaged_channels:
      return self._averaged_data.getChannelOperators(channel)
    
    elif channel in self.rotated_channels:
      return self._rotated_data.getChannelOperators(channel)

    elif channel in self.raw_channels:
      return self._raw_data.getChannelOperators(channel)

    return SortedSet()
  
  def getAveragedOperators(self, channel):
    if channel in self.averaged_channels:
      return self._averaged_data.getChannelOperators(channel)

    return SortedSet()
  
  def getRotatedOperators(self, channel):
    if channel in self.rotated_channels:
      return self._rotated_data.getChannelOperators(channel)

    return SortedSet()

  def getChannelDataFiles(self, channel):
    if channel in self.averaged_channels:
      return self._averaged_data.getChannelDataFiles(channel)

    if channel in self.raw_channels:
      return self._raw_data.getChannelDataFiles(channel)

    return DataFiles()
  
  def getAveragedDataFiles(self, channel):
    if channel in self.averaged_channels:
      return self._averaged_data.getChannelDataFiles(channel)

    return DataFiles()
  
  def getRotatedDataFiles(self, channel):
    if channel in self.rotated_channels:
      return self._rotated_data.getChannelDataFiles(channel)

    return DataFiles()

  def hasCorrelator(self, correlator):
    return self._averaged_data.hasCorrelator(correlator) or self._raw_data.hasCorrelator(correlator)

  def getChannelsLargestTRange(self, *channels):
    time_extent = self.ensemble_info.getLatticeTimeExtent()
    ave_trange = self._averaged_data.getChannelsLargestTRange(time_extent, *channels)
    raw_trange = self._raw_data.getChannelsLargestTRange(time_extent, *channels)

    if ave_trange != (-1, -1):
      return ave_trange
    else:
      return raw_trange

  def getOperatorSetSmallestTRange(self, operator_set):
    time_extent = self.ensemble_info.getLatticeTimeExtent()
    ave_trange = self._averaged_data.getOperatorSetSmallestTRange(time_extent, operator_set)
    raw_trange = self._raw_data.getOperatorSetSmallestTRange(time_extent, operator_set)
    if ave_trange != (-1, -1):
      return ave_trange
    elif raw_trange != (-1, -1):
      return raw_trange
    else:
      logging.warning( f"Operator set {operator_set.name} cannot find data files for all operators")
      return raw_trange

  # def _readRotatedData(self, rotated_basis):
  #   if rotated_basis in self._rotated_data:
  #     return self._rotated_data[rotated_basis]

  #   data_file = self.rotated_datafile(rotated_basis)
  #   if not os.path.isfile(data_file):
  #     logging.error(f"No basis {rotated_basis.name} with {rotated_basis.pivot_info}")

  #   file_type = sigmond.getFileID(data_file)
  #   if file_type == sigmond.FileType.Bins:
  #     sigmond_handler = sigmond.BinsGetHandler(self.bins_info, {data_file})
  #     fileinfo = FileInfo(data_file, sigmond.FileType.Bins)
  #   elif file_type == sigmond.FileType.Samplings:
  #     sigmond_handler = sigmond.SamplingsGetHandler(self.bins_info, self.sampling_info, {data_file})
  #     fileinfo = FileInfo(data_file, sigmond.FileType.Samplings)
  #   else:
  #     logging.critical(f"Invalid rotated data file {data_file}")

  #   rotated_data = CorrelatorData()
  #   for mc_obs in sigmond_handler.getKeys():
  #     if mc_obs.isCorrelatorAtTime():
  #       corr = mc_obs.getCorrelatorInfo()
  #       t = mc_obs.getCorrelatorTimeIndex()
  #       rotated_data.addCorrelator(corr, t, t, fileinfo)

  #     elif mc_obs.isVEV():
  #       rotated_data.addVEV(vev, fileinfo)

  #     else:
  #       logging.warning(f"Ignoring observable '{mc_obs}'")

  #   self._rotated_data[rotated_basis] = rotated_data
  #   return rotated_data

  def _findRawData(self):
    data_files = DataFiles()
    for data_dir in self.raw_data_dirs:
      data_files += _find_data_files(data_dir)

    logging.info("Searching through found raw data files...")
    self._raw_data += self._findLaphData(data_files)

    if data_files.bin_files:
      logging.info("Reading bin data")
      for bin_file in tqdm.tqdm(data_files.bin_files):
        self._raw_data += self._findSigmondData(bin_file, sigmond.FileType.Bins)

    if data_files.sampling_files:
      logging.info("Reading sampling data")
      for smp_file in tqdm.tqdm(data_files.sampling_files):
        self._raw_data += self._findSigmondData(smp_file, sigmond.FileType.Samplings)

    self.raw_data_files = data_files
    logging.info("done")

  def findAveragedData(self, rotated=False):
    if rotated:
      data_dirs = self.rotated_datadirs
      dtype = "rotated"
    else:
      data_dirs = self.averaged_datadirs
      dtype = "averaged"

    data_files = _find_data_files(data_dirs[0])
    for data_dir in data_dirs[1:]:
      data_files += _find_data_files(data_dir)

    logging.info(f"Searching through found {dtype} data files...")

    self._averaged_data += self._findLaphData(data_files)

    if data_files.bin_files:
      logging.info("Reading bin data")
      for bin_file in tqdm.tqdm(data_files.bin_files):
        if rotated:
          self._rotated_data += self._findSigmondData(bin_file, sigmond.FileType.Bins)
        else:
          self._averaged_data += self._findSigmondData(bin_file, sigmond.FileType.Bins)

    if data_files.sampling_files:
      logging.info("Reading sampling data")
      for sampling_file in tqdm.tqdm(data_files.sampling_files):
        if rotated:
          self._rotated_data += self._findSigmondData(sampling_file, sigmond.FileType.Samplings)
        else:
          self._averaged_data += self._findSigmondData(sampling_file, sigmond.FileType.Samplings)

    if rotated:
      self.rotated_data_files = data_files
    else:
      self.averaged_data_files = data_files
    logging.info("done")
    return data_files
  
  def findRotatedData(self):
    return self.findAveragedData(True)

  def _findLaphData(self, data_files):
    laph_corr_handler = sigmond.BLCorrelatorDataHandler(data_files.bl_corr_files, set(), set(), self.ensemble_info)
    laph_vev_handler = sigmond.BLVEVDataHandler(data_files.bl_vev_files, set(), self.ensemble_info)

    data = CorrelatorData()
    if (corr_set := laph_corr_handler.getCorrelatorSet()):
      logging.info("Reading LapH correlators")
      for correlator_info in tqdm.tqdm(corr_set):
        tmin, tmax = laph_corr_handler.getTimeSepRange(correlator_info)
        filename = laph_corr_handler.getFileName(correlator_info)
        fileinfo = FileInfo(filename, sigmond.FileType.Correlator)
        data.addCorrelator(correlator_info, tmin, tmax, fileinfo)

    if (vev_set := laph_vev_handler.getOperatorSet()):
      logging.info("Reading LapH vevs")
      for vev in tqdm.tqdm(vev_set):
        filename = laph_vev_handler.getFileName(vev)
        fileinfo = FileInfo(filename, sigmond.FileType.VEV)
        data.addVEV(vev, fileinfo)

    return data

  def _findSigmondData(self, filename, filetype):
    if filetype is sigmond.FileType.Bins:
      sigmond_handler = sigmond.BinsGetHandler(self.bins_info, {filename})
    elif filetype is sigmond.FileType.Samplings:
      sigmond_handler = sigmond.SamplingsGetHandler(self.bins_info, self.sampling_info, {filename})
    else:
      logging.error(f"Invalid filetype {filetype} passed to DataHandler::_getSigmondData")

    fileinfo = FileInfo(filename, filetype)

    data = CorrelatorData()
    for mc_obs in sigmond_handler.getKeys():
      if mc_obs.isCorrelatorAtTime():
        corr = mc_obs.getCorrelatorInfo()
        t = mc_obs.getCorrelatorTimeIndex()
        data.addCorrelator(corr, t, t, fileinfo)

      elif mc_obs.isVEV():
        data.addVEV(vev, fileinfo)

      else:
        logging.warning(f"Ignoring observable '{mc_obs}'")
        
    return data


def _find_data_files(data_dir):

  logging.info(f"Searching '{data_dir}' for data...")
  
  data_files = DataFiles()

  if os.path.isdir(data_dir):
    for root, dirs, files in os.walk(data_dir, topdown=True):
      for filename in files:
        full_filename = os.path.join(root, filename)
        _find_data_file(full_filename, data_files)

  elif os.path.isfile(data_dir):
      _find_data_file(data_dir, data_files)

  if data_files.empty:
    logging.info(f"Found no data in '{data_dir}'")

  return data_files

def _find_data_file(full_filename, data_files):
  file_list_infos = dict()
  bl_corr_files = set()
  bl_vev_files = set()
  bin_files = set()
  smp_files = set()
  roots = []
  try:
    with h5py.File(full_filename, 'r') as h5py_file:
      if h5py_file['Info']['FIdentifier'][()].decode()=="Sigmond--SamplingsFile":
        file_type = sigmond.FileType.Samplings
      elif h5py_file['Info']['FIdentifier'][()].decode()=="Sigmond--BinsFile":
        file_type = sigmond.FileType.Bins
      roots = list(h5py_file.keys())
      roots.remove("Info")
  except OSError as err:
    pass

  if not roots: 
    try:
      file_type = sigmond.getFileID(full_filename)
    except ValueError:
      logging.warning(f"Invalid file '{full_filename}'")
      return

  if file_type == sigmond.FileType.Correlator:
    logging.info(f"Found BasicLaph Correlator file '{full_filename}'")
    bl_corr_files.add(full_filename)
  elif file_type == sigmond.FileType.VEV:
    logging.info(f"Found BasicLaph VEV file '{full_filename}'")
    bl_vev_files.add(full_filename)
  elif file_type == sigmond.FileType.Bins:
    logging.info(f"Found Bins file '{full_filename}'")
    if roots:
        for hroot in roots:
            bin_files.add(f"{full_filename}[{hroot}]")
    else:
        bin_files.add(full_filename)
  elif file_type == sigmond.FileType.Samplings:
    logging.info(f"Found Samplings file '{full_filename}'")
    if roots:
        for hroot in roots:
            smp_files.add(f"{full_filename}[{hroot}]")
    else:
        smp_files.add(full_filename)
  else:
    logging.error(f"Unrecognized filetype {file_type}")

  data_files.addCorrelatorFiles(*bl_corr_files)
  data_files.addVEVFiles(*bl_vev_files)
  data_files.addBinFiles(*bin_files)
  data_files.addSamplingFiles(*smp_files)
