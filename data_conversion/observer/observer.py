import os
from enum import Enum, auto
import regex
import logging

import h5py
import numpy as np

class ObserverData:

  _data_files = []

  def __init__(self, *dirs, reg):
    for direc in dirs:
      for filename in glob.glob(os.path.join(direc, "*{}*".format(reg))):
        self._data_files.append(ObserverDataFile(filename))


  def write2PtData(self, channel, momentum, polarization, prop_src, prop_snk):





class DataType(Enum):
  meson_2pt = auto()
  baryon_2pt = auto()
  meson_3pt = auto()
  baryon_3pt = auto()

class Baryon(Enum):
  Nucleon = auto()
  Lambda = auto()
  Sigma = auto()
  Cascade = auto()

class Propagator(Enum):
  smeared = auto()
  local = auto()

class Polarization(Enum):
  G0 = auto()
  G1 = auto()
  G2 = auto()
  G3 = auto()


class ObserverDataFile:
  """
  TODO: close file at close
  """

  def __init__(self, filename):
    self._filename = filename
    self._h5_file_obj = h5py.File(filename, 'r')

  @property
  def filename(self):
    return self._filename

  @property
  def file_type(self):
    data_types = list()
    for data_key in self._h5_file_obj.keys():
      try:
        data_types.append(DataType[data_key])
      except KeyError as e:
        pass

    if len(data_types) != 1:
      raise TypeError("Did not find exactly one data type in : {}".format(self._filename))

    return data_types[0]

  @property
  def propagators(self):
    ...

  def getSection(self, section_title):
    """
    Used for getting sections in the infile. The argument 'section_title' is meant as
    a substring for the section name

    Args:
      section_title (str): a substring of the section name to be searched for.

    Returns: 
      a dictionary of sections that have titles in which 'section_title' is a substring.
      The keys of the dictionary correspond to the section title, and the value is another
      dictionary with the key-value pairs of the parameters in the section.
    """
    
    infile = self.infile
    matches = regex.finditer(r"\[(.*{}.*)\]\n".format(section_title), infile)

    sections = dict()
    for match in matches:
      start = match.start()
      end = infile.find('[', start+1)
      section_params = infile[start:end].split('\n')[1:]

      section = dict()

      for section_param in section_params:
        param_match = regex.search(r"^(\S+)\s+(\S.+)$", section_param.strip())

        if param_match is None:
          continue

        section[param_match.group(1)] = param_match.group(2)

      sections[match.group(1)] = section

    if not sections:
      raise ValueError("Section with substring '{}' not found".format(section))

    return sections

  @property
  def flavors(self):
    kappa = self.getSection("Lattice parameters")["Lattice parameters"]["kappa"].split()
    return len(kappa)

  @property
  def config(self):
    """
    TODO: Might want to just use the self.infile (which decodes the numpy bytes string
    """

    confs = self.getSection("Configurations")['Configurations']
    if confs['first'] != confs['last']:
      raise TypeError("More than one configuration contained in: {}".format(self._filename))

    return int(confs['first'])


  @property
  def axes(self):
    return list(self._h5_file_obj["axes"].keys())

  def printAxes(self):
    for k, v in self._h5_file_obj["axes"].items():
      print("{}:".format(k))
      ind = 0
      for el in v:
        print("\t{} - {}".format(ind, el))
        ind += 1
      print()

  @property
  def infile(self):
    return self._h5_file_obj["infile"][()].decode('UTF-8')

  def getPropagatorIds(self, prop_src, prop_snk):

    propagators = dict()
    for prop_title, prop_section in self.getSection("Propagator").items():
      prop_ind_match = regex.search(r"Propagator\s+(\d+)$", prop_title.strip())
      prop_ind = int(prop_ind_match.group(1))
      prop = list()
      if "smearing_parameter_source" in prop_section:
        prop.append(Propagator.smeared)
      else:
        prop.append(Propagator.local)

      if "smearing_parameter_sink" in prop_section:
        prop.append(Propagator.smeared)
      else:
        prop.append(Propagator.local)

      propagators[tuple(prop)] = prop_ind

    prop_ind = propagators[(prop_src, prop_snk)]

    prop_ids = list()
    for i in range(self.flavors):
      prop_ids.append(prop_ind)

    return prop_ids



  def get2PtData(self, channel, momentum, polarization, prop_src, prop_snk):
    """
    Args:
      channel (Baryon): specifies the flavor of Baryon to use
      momentum (int or list): specifying either the psq (averaging over all
          equivalent frames), or a particular momentum.
      polarization (Polarization): specifies the polarization
      prop_src (Propagator): specifies the propagator to use for the source.
      prop_snk (Propagator): specifies the propagator to use for the sink.
    """

    if self.file_type != DataType.baryon_2pt:
      raise TypeError("Only Baryon 2pt functions are supported currently")

    if momentum != [0,0,0]:
      raise TypeError("Only P=[0,0,0] currently supported")

    propagator_ids = self.getPropagatorIds(prop_src, prop_snk)

    if polarization == Polarization.G0:
      data_I_p = self._data(reverse=False, channel=channel.name, momentum=momentum, polarizations="1", propagator_list=propagator_ids)
      data_I_n = self._data(reverse=True, channel=channel.name, momentum=momentum, polarizations="1", propagator_list=propagator_ids)
   
      data_g0_p = self._data(reverse=False, channel=channel.name, momentum=momentum, polarizations="g4", propagator_list=propagator_ids)
      data_g0_n = self._data(reverse=True, channel=channel.name, momentum=momentum, polarizations="g4", propagator_list=propagator_ids)

      data_p = 0.5*(data_I_p + data_g0_p)
      data_n = 0.5*(data_I_n - data_g0_n)

      data = 0.5*(data_p + data_n)

    else:
      raise TypeError("Only G0 polarization supported currently")

    return data


  def _data(self, reverse=False, tmin=0, tmax=0, **kw_args):
    """
    Args:
      reverse (bool): determines wheter we want the backward (True) or forward (False) propagating
          correlators.
      tmin (int): specifies the smallest time source to include in the average
      tmax (int): specifies the largest time source to include in the average (a value of
          zero means to not consider a largest time source)
      **source_position (list): specifies the source position to use for retreiving data
          (shifted to zero). If absent, all source positions found are used, shifted to 0,
          and averaged over.
      kw_args (dict): containes values for in the 'axes' to return. If any given 'axis' is
          misssing, then all values for that 'axis' are returned.
    TODO:
      I think this code is pretty hacky.
      Does reversal also need complex conjugation?
    """

    _axes = self._h5_file_obj["axes"]
    _data = self._h5_file_obj[self.file_type.name]
    inds = [slice(None)] * len(_data.shape)

    dims = {}
    for dim in _data.dims:
      dims[dim.label.replace(' ', '_')] = dim._dimension

    for k, v in kw_args.items():

      if k not in _axes:
        continue

      if isinstance(v, str):
        v = v.encode()

      if len(_axes[k][()].shape) == 1:
        ind = np.where(_axes[k][()] == v)[0]
      elif len(_axes[k][()].shape) == 2:
        ind = np.where((_axes[k][()] == v).all(axis=1))[0]
      else:
        raise TypeError("Unexpected axes in {}".format(self._filename))

      if len(ind):
        ind = ind[0]
      else:
        logging.warning("Could not find {}={}".format(k, v))
        continue

      inds[dims[k]] = ind

    # shift source to zero and average over all source_positions

    if 'source_position' in kw_args:
      t_src = kw_args['source_position'][3]

      if t_src < tmin or (tmax and t_src > tmax):
        logging.warning("only source time is not inside allowed range, return None")
        return None

      if reverse:
        t_src += 1

      _data = np.roll(_data, -t_src, axis=dims['T'])
      if reverse:
        _data = np.flip(_data, axis=dims['T'])

    else:
      del inds[dims['source_position']]
      _data_srcs = []
      for ind, src_pos in enumerate(_axes['source_position']):
        t_src = src_pos[3]

        if t_src < tmin or (tmax and t_src > tmax):
          logging.warning("t_src {} outside of allowed range, skipping source {}".format(t_src, src_pos))

        if reverse:
          t_src += 1

        _data_src = np.squeeze(np.take(_data, [ind], axis=dims['source_position']), axis=dims['source_position'])
        t_ind = dims['T']
        if dims['source_position'] < t_ind:
          t_ind -= 1
        _data_src = np.roll(_data_src, -t_src, axis=t_ind)
        if reverse:
          _data_src = np.flip(_data_src, axis=t_ind)

        _data_srcs.append(_data_src)

      if _data_srcs:
        _data = np.stack(tuple(_data_srcs), axis=dims['source_position'])
        _data = np.mean(_data, axis=dims['source_position'])
      else:
        logging.warning("No sources, returning None")
        return None

    _data = _data[tuple(inds)]

    return _data

