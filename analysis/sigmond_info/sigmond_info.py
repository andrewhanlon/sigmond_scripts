import logging
from enum import Enum
from aenum import MultiValueEnum
from typing import NamedTuple

import utils.util as util
import sigmond

class Goodness(Enum):
  chisq = 'chisq'
  quality = 'qual'

class SymbolColor(Enum):
  black = "black"
  red = "red"
  blue = "blue"
  green = "green"
  yellow = "yellow"
  magenta = "magenta"
  cyan = "cyan"
  orange = "orange"
  violet = "violet"
  maroon = "maroon"

class SymbolType(Enum):
  none = "none"
  circle = "circle"
  square = "square"
  diamond = "diamond"
  triangleup = "triangleup"
  triangleleft = "triangleleft"
  triangledown = "triangledown"
  plus = "plus"
  X = "X"
  star = "star"

class EffEnergyType(MultiValueEnum):
  TimeForward = 0, "time_forward"
  TimeSymmetric = 1, "time_symmetric"
  TimeForwardPlusConst = 2, "time_forward_const"
  TimeSymmetricPlusConst = 3, "time_symmetric_const"

class PlotInfo(NamedTuple):
  corrname: str = "standard"
  timestep: int = 3
  rescale: float = 1.0
  symbol_color: SymbolColor = SymbolColor.blue
  symbol_type: SymbolType = SymbolType.circle
  max_relative_error: float = 0.0
  eff_energy_type: EffEnergyType = EffEnergyType.TimeForward

  @classmethod
  def createFromConfig(cls, options):
    plot_info = options.pop('plot_info', dict())
    util.updateOption(plot_info, 'symbol_color', SymbolColor)
    util.updateOption(plot_info, 'symbol_type', SymbolType)
    util.updateOption(plot_info, 'eff_energy_type', EffEnergyType)

    try:
      plot_info = cls(**plot_info)
    except TypeError as err:
      logging.error(f"Invalid 'plot_info' config: {err}")

    return plot_info

class FitPlotInfo(NamedTuple):
  corrname: str = "standard"
  timestep: int = 3
  symbol_color: SymbolColor = SymbolColor.blue
  symbol_type: SymbolType = SymbolType.circle
  show_approach: bool = True
  goodness: Goodness = Goodness.chisq
  max_relative_error: float = 0.0

  @classmethod
  def createFromConfig(cls, options):
    fit_plot_info = options.pop('fit_plot_info', dict())
    util.updateOption(fit_plot_info, 'symbol_color', SymbolColor)
    util.updateOption(fit_plot_info, 'symbol_type', SymbolType)
    util.updateOption(fit_plot_info, 'goodness', Goodness)

    try:
      fit_plot_info = cls(**fit_plot_info)
    except TypeError as err:
      logging.error(f"Invalid 'fit_plot_info' config: {err}")

    return fit_plot_info

class AnisotropyPlotInfo(NamedTuple):
  goodness: Goodness = Goodness.chisq
  symbol_color: SymbolColor = SymbolColor.blue
  symbol_type: SymbolType = SymbolType.circle

  @classmethod
  def createFromConfig(cls, options):
    anisotropy_plot_info = options.pop('anitostropy_plot_info', dict())
    util.updateOption(anisotropy_plot_info, 'goodness', Goodness)
    util.updateOption(anisotropy_plot_info, 'symbol_color', SymbolColor)
    util.updateOption(anisotropy_plot_info, 'symbol_type', SymbolType)

    try:
      anisotropy_plot_info = cls(**anisotropy_plot_info)
    except TypeError as err:
      logging.error(f"Invalid 'anisotropy_plot_info' config: {err}")

    return anisotropy_plot_info

class DispersionPlotInfo(NamedTuple):
  goodness: Goodness = Goodness.chisq
  symbol_color: SymbolColor = SymbolColor.blue
  symbol_type: SymbolType = SymbolType.circle

  @classmethod
  def createFromConfig(cls, options):
    dispersion_plot_info = options.pop('dispersion_plot_info', dict())
    util.updateOption(dispersion_plot_info, 'goodness', Goodness)
    util.updateOption(dispersion_plot_info, 'symbol_color', SymbolColor)
    util.updateOption(dispersion_plot_info, 'symbol_type', SymbolType)

    try:
      dispersion_plot_info = cls(**dispersion_plot_info)
    except TypeError as err:
      logging.error(f"Invalid 'dispersion_plot_info' config: {err}")

    return dispersion_plot_info


class TMinPlotInfo(NamedTuple):
  obsname: str = "standard"
  symbol_type: SymbolType = SymbolType.circle
  goodfit_color: SymbolColor = SymbolColor.blue
  badfit_color: SymbolColor = SymbolColor.red
  correlatedfit_hollow: bool = True 
  uncorrelatedfit_hollow: bool = False
  quality_threshold: float = 0.1
  correlated_threshold: float = 1.0

  @classmethod
  def createFromConfig(cls, options):
    tmin_plot_info = options.pop('tmin_plot_info', dict())
    util.updateOption(tmin_plot_info, 'symbol_type', SymbolType)
    util.updateOption(tmin_plot_info, 'goodfit_color', SymbolColor)
    util.updateOption(tmin_plot_info, 'badfit_color', SymbolColor)

    try:
      tmin_plot_info = cls(**tmin_plot_info)
    except TypeError as err:
      logging.error(f"Invalid 'tmin_plot_info' config: {err}")

    return tmin_plot_info

  def setChosenFit(self, fit_info):
    new_plot_info = TMinPlotInfo(
        self.obsname, self.symbol_type, self.goodfit_color, self.badfit_color,
        self.correlatedfit_hollow, self.uncorrelatedfit_hollow, self.quality_threshold,
        self.correlated_threshold)

    return new_plot_info

class RotateMode(Enum):
  bins = 'bins'
  samplings = 'samplings'
  samplings_unsubt = 'samplings_unsubt'
  samplings_all = 'samplings_all'

class DataFormat(MultiValueEnum):
  bins = 'bin', 'bins'
  samplings = 'smp', 'samplings'

class PivotType(Enum):
  SinglePivot = 'single_pivot'
  SingleTimePivot = 'single_time_pivot'
  RollingPivot = 'rolling_pivot'
  PrincipalAxes = 'principal_axes'

class PivotInfo:
  """ specifies all the needed info about the pivot used """

  def __init__(self, pivot_type, **extra_info):
    self.pivot_type = PivotType(pivot_type)
    
    try:
      self.norm_time = extra_info['norm_time']
      self.metric_time = extra_info['metric_time']
      self.diagonalize_time = extra_info['diagonalize_time']
      self.max_condition_number = extra_info['max_condition_number']

    except KeyError as err:
      logging.error(f"Missing required key {err} for {self.pivot_type}")

  def __repr__(self):
    return f"{self.pivot_type.value}_n{self.norm_time}_m{self.metric_time}_d{self.diagonalize_time}_c{self.max_condition_number}"

  def __str__(self):
    return f"{self.pivot_type.name} - ({self.norm_time},{self.metric_time},{self.diagonalize_time},{self.max_condition_number})"

  def __hash__(self):
    return hash(repr(self))

  def __eq__(self, other):
    return repr(self) == repr(other)

  def __ne__(self, other):
    return not self.__eq__(other)

class ScatteringParticle:
  """ Scattering Particle

  This class is used for representing a single particle with some 
  P^2. This is convenient for strings like 'name(3)' for representing
  a single particle called 'name' with P^2=3.
  """

  def __init__(self, name, psq, irrep=None):
    self.name = name
    self.psq = psq
    self.irrep = irrep

  @classmethod
  def create(cls, particle):
    try:
      name, arg = particle[:-1].split('(')
      if '_' in arg:
        psq, irrep = arg.split('_')
        psq = int(psq)
      else:
        psq = int(arg)
        irrep = None

    except ValueError:
      logging.error(f"Invalid ScatteringParticle '{particle}'")

    return cls(name, psq, irrep)

  def __str__(self):
    if self.irrep is None:
      return f"{self.name}({self.psq})"
    else:
      return f"{self.name}({self.psq}_{self.irrep})"

  def __repr__(self):
    if self.irrep is None:
      return f"{self.name}({self.psq})"
    else:
      return f"{self.name}({self.psq}_{self.irrep})"

  def __hash__(self):
    return hash(self.__repr__())

  def __eq__(self, other):
    if isinstance(other, self.__class__):
      return self.__repr__() == other.__repr__()
    return NotImplemented

  def __ne__(self, other):
    return not self.__eq__(other)


class NonInteractingLevel:
  """ Non-interacting Level

  Basically just a list of ScatteringParticle objects
  that corresponds to a single non-interacting level.
  """
  
  def __init__(self, *particles):
    self._particles = list(particles)

  def addParticle(self, particle):
    self._particles.append(particle)

  def names_str(self):
    return ''.join(self._particles)

  def __iter__(self):
    yield from self._particles


class NonInteractingOperators(NamedTuple):
  non_interacting_level: NonInteractingLevel
  operators: list

  @classmethod
  def create(cls, scattering_particles, non_interacting_level):
    the_non_interacting_level = NonInteractingLevel()
    operators = list()
    for particle in non_interacting_level:
      particle = ScatteringParticle.create(particle)
      the_non_interacting_level.addParticle(particle)
      try:
        operator = scattering_particles[particle]
        operators.append(operator)
      except KeyError:
        logging.error("Scattering particle {particle} not specified")

    return cls(the_non_interacting_level, operators)


def getMinimizerInfo(options):
  if (minimizer_info_conf := options.pop('minimizer_info', False)):
    minimizer = minimizer_info_conf.pop("minimizer", "lmder")[0].upper()
    parameter_rel_tol = float(minimizer_info_conf.pop("parameter_rel_tol", 1e-6))
    chisquare_rel_tol = float(minimizer_info_conf.pop("chisquare_rel_tol", 1e-4))
    max_iterations = minimizer_info_conf.pop("max_iterations", 1024)
    verbosity = minimizer_info_conf.pop("verbosity", "low")[0].upper()

    util.check_extra_keys(minimizer_info_conf, "minimizer_info")

    minimizer_info = sigmond.MinimizerInfo(minimizer, parameter_rel_tol, chisquare_rel_tol, max_iterations, verbosity)
  else:
    minimizer_info = sigmond.MinimizerInfo()

  return minimizer_info

