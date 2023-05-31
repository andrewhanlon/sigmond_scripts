import xml.etree.ElementTree as ET
import logging
from typing import NamedTuple
from sortedcontainers import SortedSet
from aenum import MultiValueEnum
import regex

import sigmond
import utils.util as util
import operator_info.operator



class FitModel(MultiValueEnum):
  TimeForwardSingleExponential = 0, "1-exp"
  TimeSymSingleExponential = 1, "1-exp-sym"
  TimeForwardSingleExponentialPlusConstant = 2, "1-exp-const"
  TimeSymSingleExponentialPlusConstant = 3, "1-exp-sym-const"
  TimeForwardTwoExponential = 4, "2-exp"
  TimeSymTwoExponential = 5, "2-exp-sym"
  TimeForwardTwoExponentialPlusConstant = 6, "2-exp-const"
  TimeSymTwoExponentialPlusConstant = 7, "2-exp-sym-const"
  TimeForwardGeomSeriesExponential = 8, "geom"
  TimeSymGeomSeriesExponential = 9, "geom-sym"
  LogTimeForwardSingleExponential = 10, "log-1-exp"
  LogTimeForwardTwoExponential = 11, "log-2-exp"
  TimeForwardMultiExponential = 12, "multi-exp"
  TimeForwardSTIGeomSeriesExponential = 13, "sti-geom"
  TimeForwardTruncGeomSeriesExponential = 14, "trunc-geom"
  TimeForwardDoubleExpRatio = 15, "two-exp-ratio"
  TimeForwardTwoIndExp = 16, "two-ind-exp"
  TimeForwardGeomSeriesSTI = 17, "geom-sti" #colin wrote this one
  TimeForwardThreeIndExp = 18, "3-ind-exp"

  @property
  def short_name(self):
    return FIT_MODEL_SHORT_NAMES[self]

  @property
  def has_gap(self):
    return 'SqrtGapToSecondEnergy' in FitInfo.PARAMETERS[self]

  @property
  def has_const(self):
    return 'AddedConstant' in FitInfo.PARAMETERS[self]


FIT_MODEL_SHORT_NAMES = {
    FitModel.TimeForwardSingleExponential: "1-exp",
    FitModel.TimeSymSingleExponential: "1-exp-sym",
    FitModel.TimeForwardSingleExponentialPlusConstant: "1-exp-const",
    FitModel.TimeSymSingleExponentialPlusConstant: "1-exp-sym-const",
    FitModel.TimeForwardTwoExponential: "2-exp",
    FitModel.TimeSymTwoExponential: "2-exp-sym",
    FitModel.TimeForwardTwoExponentialPlusConstant: "2-exp-const",
    FitModel.TimeSymTwoExponentialPlusConstant: "2-exp-sym-const",
    FitModel.TimeForwardGeomSeriesExponential: "geom",
    FitModel.TimeSymGeomSeriesExponential: "geom-sym",
    FitModel.LogTimeForwardSingleExponential: "log-1-exp",
    FitModel.LogTimeForwardTwoExponential: "log-2-exp",
    FitModel.TimeForwardMultiExponential: "multi-exp",
    FitModel.TimeForwardSTIGeomSeriesExponential: "sti-geom",
    FitModel.TimeForwardTruncGeomSeriesExponential: "trunc-geom",
    FitModel.TimeForwardDoubleExpRatio: "two-exp-ratio",
    FitModel.TimeForwardTwoIndExp: "two-ind-exp",
    FitModel.TimeForwardGeomSeriesSTI: "geom-sti",
    FitModel.TimeForwardThreeIndExp: "3-ind-exp",
}


class FitInfo:

  PARAMETERS = {
      FitModel.TimeForwardSingleExponential: [
          "Energy",
          "Amplitude",
      ],
      FitModel.TimeSymSingleExponential: [
          "Energy",
          "Amplitude",
      ],
      FitModel.TimeForwardSingleExponentialPlusConstant: [
          "Energy",
          "Amplitude",
          "AddedConstant",
      ],
      FitModel.TimeSymSingleExponentialPlusConstant: [
          "Energy",
          "Amplitude",
          "AddedConstant",
      ],
      FitModel.TimeForwardTwoExponential: [
          "FirstEnergy",
          "FirstAmplitude",
          "SqrtGapToSecondEnergy",
          "SecondAmplitudeRatio",
      ],
      FitModel.TimeSymTwoExponential: [
          "FirstEnergy",
          "FirstAmplitude",
          "SqrtGapToSecondEnergy",
          "SecondAmplitudeRatio",
      ],
      FitModel.TimeForwardTwoExponentialPlusConstant: [
          "FirstEnergy",
          "FirstAmplitude",
          "SqrtGapToSecondEnergy",
          "SecondAmplitudeRatio",
          "AddedConstant",
      ],
      FitModel.TimeSymTwoExponentialPlusConstant: [
          "FirstEnergy",
          "FirstAmplitude",
          "SqrtGapToSecondEnergy",
          "SecondAmplitudeRatio",
          "AddedConstant",
      ],
      FitModel.TimeForwardGeomSeriesExponential: [
          "FirstEnergy",
          "FirstAmplitude",
          "SqrtGapToSecondEnergy",
          "SecondAmplitudeRatio",
      ],
      FitModel.TimeSymGeomSeriesExponential: [
          "FirstEnergy",
          "FirstAmplitude",
          "SqrtGapToSecondEnergy",
          "SecondAmplitudeRatio",
      ],
      FitModel.LogTimeForwardSingleExponential: [
          "Energy",
          "LogAmplitude",
      ],
      FitModel.LogTimeForwardTwoExponential: [
          "FirstEnergy",
          "LogFirstAmplitude",
          "SqrtGapToSecondEnergy",
          "SecondAmplitudeRatio",
      ],
      FitModel.TimeForwardMultiExponential: [
          "E0","E1","E2","E3","E4",
          "A0","A1","A2","A3","A4",
      ],
      FitModel.TimeForwardSTIGeomSeriesExponential: [
          "FirstEnergy",
          "FirstAmplitude",
          "SqrtGapToSecondEnergy",
          "SecondAmplitudeRatio",
          "STIAmplitudeRatio1",
          "STIAmplitudeRatio2"
      ],
      FitModel.TimeForwardTruncGeomSeriesExponential: [
          "FirstEnergy",
          "FirstAmplitude",
          "SqrtGapToSecondEnergy",
          "SecondAmplitudeRatio",
      ],
      FitModel.TimeForwardDoubleExpRatio: [
          "Energy",
          "Amplitude",
          "NumGap",
          "NumGapAmp",
          "SH1Gap",
          "SH1GapAmp",
          "SH2Gap",
          "SH2GapAmp",
      ],
      FitModel.TimeForwardTwoIndExp: [
          "Energy",
          "Amplitude",
          "Energy1",
          "Amplitude1",
      ],
      FitModel.TimeForwardThreeIndExp: [
          "Energy",
          "Amplitude",
          "Energy1",
          "Amplitude1",
          "Energy2",
          "Amplitude2",
      ],
      FitModel.TimeForwardGeomSeriesSTI: [
          "FirstEnergy",
          "FirstAmplitude",
          "SqrtGapToSecondEnergy",
          "SecondAmplitudeRatio",
          "STIAmplitude1",
          "STIAmplitude2",
          "STIAmplitude3",
          "STIAmplitude4",
      ],
  }

  def __init__(self, operator, model, tmin, tmax, subtractvev=False, ratio=False,
               exclude_times=[], noise_cutoff=0.0, non_interacting_operators=None, 
               tmin_max=-1, tmax_min=-1, max_level = 6, initial_gap=1.0, repeating_gap=1.0,
              sim_fit = False, initial_params = {}):
    """
    Args:
      operator (sigmondbind.OperatorInfo):
      model (FitModel): the model to use
      tmin (int): minimum time slice
      tmax (int): maximum time slice
      subtractvev (bool): whether the subtracted vev should be used
      exclude_times (list of ints): A list of times to exclude from
          the fit
      noise_cutoff (float): A error cutoff in the fit
      tmin_max (int): If doing a TminVary Fit, this is needed
      max_level (int): for multiseries fit, max fit form level
    """

    self.operator = operator
    self.model = model
    self.subtractvev = subtractvev
    self.tmin = tmin
    self.tmax = tmax
    self.ratio = ratio
    self.exclude_times = exclude_times
    self.noise_cutoff = noise_cutoff
    self.max_level = max_level #for multiseries fit
    self.initial_gap = initial_gap #for multiseries fit
    self.repeating_gap = repeating_gap #for multiseries fit
    self.sim_fit = sim_fit #for simultaneous fits #not actually set up yet
    self.initial_params = initial_params #for TimeForwardDoubleExpRatio

    self.non_interacting_operators = non_interacting_operators

    self.tmin_max = tmin_max
    self.tmax_min = tmax_min #for tmax vary plot

  @classmethod
  def createFromObservable(cls, obs_info):
    pattern = r"^(?P<comp_op_str>\S+)T(?P<tmin>\d+)-(?P<tmax>\d+)(-(?P<tmin_max>\d+))?" \
              r"(?P<subtractvev>S)?(?P<ratio>R)?(C(?P<cutoff>\d*\.?\d+))?$"
    matches = regex.match(pattern, obs_info.getObsName())
    if matches is None:
      print(obs_info)
      logging.warning("Failed at constructing FitInfo")
      return None

    op = operator_info.operator.Operator.createFromCompact(matches.group('comp_op_str'))
    min_time = int(matches.group('tmin'))
    max_time = int(matches.group('tmax'))
    tmin_max = -1 if matches.group('tmin_max') is None else int(matches.group('tmin_max'))
    ratio = False if matches.group('ratio') is None else True
    subtractvev = False if matches.group('subtractvev') is None else True

    noise_cutoff = 0.0
    if matches.group('cutoff') is not None:
      noise_cutoff = float(matches.group('cutoff')) / 10.

    id_index = obs_info.getObsIndex()
    model_num = id_index % 100
    model_type = FitModel(model_num)

    exclude_times_int = id_index // 1000

    exclude_times = SortedSet()
    if exclude_times_int > 0:
      excludes_time_str = str(exclude_times_int)
      exclude_time = ""
      for pos, digit in enumerate(excludes_time_str):
        if digit == "0" and pos + 1 < len(excludes_time_str) and excludes_time_str[pos+1] != "0":
          exclude_times.add(int(current_time))
          current_time = ""
        else:
          exclude_time += digit

      if exclude_time != "":
        exclude_times.add(int(exclude_time))

    exclude_times = list(exclude_times)

    return cls(op, model_type, min_time, max_time, subtractvev, ratio, exclude_times,
               noise_cutoff, tmin_max=tmin_max)

  @classmethod
  def createFromConfig(cls, config):
    try:
      operator = operator_info.operator.Operator(config.pop('operator'))
      model = FitModel(config.pop('model'))

      return cls(operator, model, **config)

    except KeyError as err:
      logging.error("Missing {err} in fit_info")

  @property
  def has_gap(self):
    return 'SqrtGapToSecondEnergy' in self.PARAMETERS[self.model]

  @property
  def has_const(self):
    return 'AddedConstant' in self.PARAMETERS[self.model]

  @property
  def num_params(self):
    return len(self.PARAMETERS[self.model])

  @property
  def model_value(self):
    return self.model.value

  @property
  def model_name(self):
    return self.model.name

  @property
  def fit_type(self):
    _fit_type = "TemporalCorrelator"
    if self.is_log_fit:
      _fit_type = f"Log{_fit_type}"

    if self.ratio:
      _fit_type = f"{_fit_type}InteractionRatio"

    
    if self.is_tmin_vary:
      _fit_type = f"{_fit_type}TminVary"

    if self.is_tmax_vary:
      _fit_type = f"{_fit_type}TmaxVary"

    if self.sim_fit:
      _fit_type = f"NSimTemporalCorrelator"
    
    return _fit_type

  @property
  def is_tmin_vary(self):
    return self.tmin_max > 0

  @property
  def is_tmax_vary(self):
    return self.tmax_min > 0

  @property
  def is_log_fit(self):
    return self.model.name.startswith("Log")

  @property
  def plotfile(self):
    return f"fit_{self.operator!r}_{self.obs_name}_{self.obs_id(0)}"

  def xml(self):
    xml = ET.Element(f"{self.fit_type}Fit")

    if self.ratio:
      ratio_xml = ET.SubElement(xml, "Ratio")
      ratio_xml.append(self.operator.ratio_op.xml())
      interacting_xml = ET.SubElement(xml, "InteractingOperator")
      interacting_xml.append(self.operator.xml())
      if self.subtractvev and self.operator.vev:
        ET.SubElement(interacting_xml, "SubtractVEV")

      if self.non_interacting_operators:
        for non_interacting_operator in self.non_interacting_operators.operators:
          non_interacting_xml = ET.SubElement(xml, "NonInteractingOperator")
          non_interacting_xml.append(non_interacting_operator.xml())
          if self.subtractvev and non_interacting_operator.vev:
            ET.SubElement(non_interacting_xml, "SubtractVEV")

      else:
        logging.error("No non-interacting operators passed to xml for ratio fit")

    else:
      xml.append(self.operator.xml())
      if self.subtractvev and self.operator.vev:
        ET.SubElement(xml, "SubtractVEV")

    if self.is_tmin_vary > 0:
      ET.SubElement(xml, "TminFirst").text = str(self.tmin)
      ET.SubElement(xml, "TminLast").text = str(self.tmin_max)
      ET.SubElement(xml, "Tmax").text = str(self.tmax)
    elif self.is_tmax_vary > 0:
      ET.SubElement(xml, "TmaxFirst").text = str(self.tmax_min)
      ET.SubElement(xml, "TmaxLast").text = str(self.tmax)
      ET.SubElement(xml, "Tmin").text = str(self.tmin)
    else:
      ET.SubElement(xml, "MinimumTimeSeparation").text = str(self.tmin)
      ET.SubElement(xml, "MaximumTimeSeparation").text = str(self.tmax)

    if self.exclude_times:
      ET.SubElement(xml, "ExcludeTimes").text = " ".join(str(t) for t in self.exclude_times)

    if self.noise_cutoff and not self.is_tmin_vary and not self.is_tmax_vary:
      ET.SubElement(xml, "LargeTimeNoiseCutoff").text = str(self.noise_cutoff)

    if self.is_log_fit:
      model_xml = ET.SubElement(xml, "LogModel")
    else:
      model_xml = ET.SubElement(xml, "Model")

    ET.SubElement(model_xml, "Type").text = self.model_name
    if self.model_name == "TimeForwardMultiExponential":
        ET.SubElement(model_xml, "MaxLevel").text = str(self.max_level)
        ET.SubElement(model_xml, "InitialGap").text = str(self.initial_gap)
        ET.SubElement(model_xml, "RepeatingGap").text = str(self.repeating_gap)
        
    if self.model_name == "TimeForwardGeomSeriesSTI":
        ET.SubElement(model_xml, "NumberOfSTIEnergies").text = str(self.max_level)
        ET.SubElement(model_xml, "STIEnergyGap").text = str(self.initial_gap)
        ET.SubElement(model_xml, "STIEnergyStep").text = str(self.repeating_gap)
        
    param_count = 0
    for param in self.PARAMETERS[self.model]:
      param_xml = ET.SubElement(model_xml, param)
      ET.SubElement(param_xml, "Name").text = self.obs_name
      ET.SubElement(param_xml, "IDIndex").text = str(self.obs_id(param_count))
      if self.model_name == "TimeForwardDoubleExpRatio" and self.initial_params:
        if param in self.initial_params.keys():
            ET.SubElement(param_xml, "InitialValue").text = str(self.initial_params[param])
#       if self.model_name == "TimeForwardMultiExponential":
#         ET.SubElement(param_xml, "PriorValue").text = str(self.initial_gap)
#         ET.SubElement(param_xml, "PriorWidth").text = str(self.repeating_gap)
            
      param_count += 1

    return xml

  @property
  def obs_name(self):
    name = f"{self.operator.compact_str}T{self.tmin}-{self.tmax}"
    if self.is_tmin_vary > 0:
      name += f"-{self.tmin_max}"
    
    if self.is_tmax_vary > 0:
      name += f"-{self.tmax_min}" #is not set up properly so that when reading the logfile
                                    #sigmond scripts thinks this is a tmin plot, but I have no good ideas rn
                                    #to fix it

    if self.subtractvev:
      name += "S"

    if self.ratio:
      name += "R"

    if self.noise_cutoff > 0.0:
      noise_cutoff = self.noise_cutoff * 10
      name += "C{val:.0f}".format(val=noise_cutoff)

    return name

  def obs_id(self, param_num):
    id_index = ""
    if self.exclude_times:
      id_index += "0".join(map(str, self.exclude_times))
    id_index += f"{param_num}{self.model_value:02d}"
    return int(id_index)

  @property
  def energy_index(self):
    return 0

  @property
  def amplitude_index(self):
    if (self.model==FitModel.TimeForwardMultiExponential): # or (self.model==FitModel.TimeForwardDoubleExpRatio)
        return 5
    return 1

  @property
  def energy_observable(self):
    return sigmond.MCObsInfo(self.obs_name, self.obs_id(self.energy_index))

  def fit_param_obs(self, index):
    return sigmond.MCObsInfo(self.obs_name, self.obs_id(index))

  @property
  def amplitude_observable(self):
    return sigmond.MCObsInfo(self.obs_name, self.obs_id(self.amplitude_index))

  def __str__(self):
    return util.xmltostr(self.xml())

  def __repr__(self):
    return str(self.__cmp())

  def __cmp(self):
    tup_rep = (self.operator, self.model_value, self.tmin, self.tmax,
               ":".join(map(str, self.exclude_times)), self.ratio, self.noise_cutoff)

    return tup_rep

  def __hash__(self):
    return hash(self.__repr__())

  def __eq__(self, other):
    return self.__repr__() == other.__repr__()

  def __ne__(self, other):
    return not self.__eq__(other)

  def __lt__(self, other):
    if isinstance(other, self.__class__):
      return self.__cmp() < other.__cmp()
    return NotImplemented

  def __gt__(self, other):
    if isinstance(other, self.__class__):
      return self.__cmp() > other.__cmp()
    return NotImplemented

  def __ge__(self, other):
    if isinstance(other, self.__class__):
      return self.__cmp() >= other.__cmp()
    return NotImplemented

  def __le__(self, other):
    if isinstance(other, self.__class__):
      return self.__cmp() <= other.__cmp()
    return NotImplemented
