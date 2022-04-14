import os
import logging
import xml.etree.ElementTree as ET
from abc import ABCMeta, abstractmethod
import regex

from typing import NamedTuple
from sortedcontainers import SortedSet, SortedDict

import sigmond
import utils.util as util
from sigmond_info.fit_info import FitInfo
from operator_info.operator import Operator


class SigmondLog(metaclass=ABCMeta):

  def __init__(self, logfile):
    self.logfile = logfile
    try:
      log_xml_root = ET.parse(logfile).getroot()
    except ET.ParseError:
      logging.critical("Bad logfile: {}".format(logfile))
      return

    self.parse(log_xml_root)

  @abstractmethod
  def parse(self, log_xml):
    pass


class DiagFracErrors(NamedTuple):
  metric: float
  matrix: float

class DeviationFromZero(NamedTuple):
  status: bool
  max: float
  one: float
  two: float
  three: float
  four: float


class RotationLog(SigmondLog):

  def parse(self, log_xml_root):
    self.failed = False
    if log_xml_root.find("Task/Error") is not None:
      self.failed = True
      return

    rotation_tasks_xml = log_xml_root.findall("Task/DoCorrMatrixRotation")
    if len(rotation_tasks_xml) != 1:
      logging.warning("Could not find single <DoCorrMatrixRotation> tag")
      return

    rotation_task_xml = rotation_tasks_xml[0]

    if rotation_task_xml.find("SinglePivot") is None:
      logging.warning("Reading of non 'SinglePivot' log files not currently supported")
      return NotImplemented

    pivot_xml = rotation_task_xml.find("SinglePivot/InitiateNew/CreatePivot")
    self.diag_corr_errors_xml = pivot_xml.find("DiagonalCorrelatorFractionalErrors")
    self.analyze_metric_xml = pivot_xml.find("AnalyzeMetric")
    self.analyze_matrix_xml = pivot_xml.find("AnalyzeMatrix")
    self.do_rotation_xml = rotation_task_xml.find("DoRotation")
    self.transformation_matrix_xml = rotation_task_xml.find("SinglePivot/InitiateNew/TransformationMatrix")

  @property
  def metric_null_space_message(self):
    _message = self.analyze_matrix_xml.findtext(
        "CheckNullSpaceCommonality/MetricNullSpace")
    if _message:
      return _message
    else:
      return "Passed"

  @property
  def diagonal_correlator_errors(self):
    corr_errors = SortedDict()
    for diag_corr_error_xml in self.diag_corr_errors_xml.iter("DiagonalCorrelator"):
      op_str = [xml.text for xml in diag_corr_error_xml.iter()
                if 'Operator' in xml.tag][0]

      metric_error = float(diag_corr_error_xml.findtext("MetricTimeFractionalError"))
      matrix_error = float(diag_corr_error_xml.findtext("MatrixTimeFractionalError"))

      corr_errors[op_str] = DiagFracErrors(metric_error, matrix_error)

    return corr_errors

  def metric_condition(self, retained=True):
    if retained:
      eigenvalues = [float(eig.text) for eig in self.analyze_metric_xml.find(
                                                   "MetricRetainedEigenvalues").iter("Value")]
      largest_eigenvalue = max(eigenvalues)
      smallest_eigenvalue = min(eigenvalues)
    else:
      eigenvalues = [float(eig.text) for eig in self.analyze_metric_xml.find(
                                                   "MetricAllEigenvalues").iter("Value")]
      largest_eigenvalue = max(eigenvalues)
      smallest_eigenvalue = min(eigenvalues)

    return round(largest_eigenvalue / smallest_eigenvalue, 2)

  def matrix_condition(self, retained=True):
    if retained:
      eigenvalues = [float(eig.text) for eig in self.analyze_matrix_xml.find(
                                                   "GMatrixRetainedEigenvalues").iter("Value")]
      largest_eigenvalue = max(eigenvalues)
      smallest_eigenvalue = min(eigenvalues)
    else:
      eigenvalues = [float(eig.text) for eig in self.analyze_matrix_xml.find(
                                                   "GMatrixAllEigenvalues").iter("Value")]
      largest_eigenvalue = max(eigenvalues)
      smallest_eigenvalue = min(eigenvalues)

    return round(largest_eigenvalue / smallest_eigenvalue, 2)

  @property
  def deviations_from_zero(self):
    deviations = SortedDict()
    for rotation_xml in self.do_rotation_xml.iter("CorrelatorRotation"):
      time = int(rotation_xml.findtext("TimeValue"))
      status = rotation_xml.findtext("Status")

      off_diagonal_xml = rotation_xml.find("OffDiagonalChecks")
      if off_diagonal_xml is not None:
        max_err = float(off_diagonal_xml.findtext("MaximumDeviationFromZero/RelativeToError"))
        percent_xml = off_diagonal_xml.find("PercentDeviationsFromZero")
        one_sigma = round(float(percent_xml.findtext("GreaterThanOneSigma")), 2)
        two_sigma = round(float(percent_xml.findtext("GreaterThanTwoSigma")), 2)
        three_sigma = round(float(percent_xml.findtext("GreaterThanThreeSigma")), 2)
        four_sigma = round(float(percent_xml.findtext("GreaterThanFourSigma")), 2)
        deviation = DeviationFromZero(status, max_err, one_sigma, two_sigma, three_sigma, four_sigma)
      else:
        deviation = DeviationFromZero(status, '', '', '', '', '')

      deviations[time] = deviation

    return deviations

  @property
  def number_levels(self):
    return len(list(self.analyze_matrix_xml.find("GMatrixRetainedEigenvalues")))

  @property
  def improved_operators(self):
    improved_ops = list()
    for op in self.transformation_matrix_xml.find("ImprovedOperators").iter("ImprovedOperator"):
      op_info = list()
      op_info.append( op.find("OpName").findtext("GIOperatorString") )
      for term in op.iter("OpTerm"):
        op_info.append( term.findtext("GIOperatorString") )
        op_info.append( term.findtext("Coefficient") ) #.replace('(', '').replace(')', '').split(',') )
      improved_ops.append( op_info )
    # [ [name, term1, coeff1, term2, coeff2...], [name, term1, coeff1...  
    return improved_ops


class FitResult(NamedTuple):
  chisq: float
  quality: float
  energy: str
  reconstructed_energy: str #if ratio, otherwise it's just the energy
  amplitude: str
  gap: str
  const: str
  covariance_condition: float

class FitLog(SigmondLog):

  def parse(self, log_xml_root):
    self.fits = SortedDict()
    for task_xml in log_xml_root.findall("Task"):
      count = task_xml.findtext("Count")
      fit_xml = task_xml.find("DoFit")
      if fit_xml is None or fit_xml.find("Error") is not None:
        continue

      try:
        eigenvalues = list()
        for eigenvalue in fit_xml.find("CovarianceMatrixEigenvalues"):
          eigenvalues.append(float(eigenvalue.text))

        eigenvalues.sort()
        cov_cond = float(fit_xml.findtext("CovarianceMatrixConditionNumber"))
        fit_results = fit_xml.find("BestFitResult")
        chisq_dof = float(fit_results.findtext("ChiSquarePerDof"))
        quality = float(fit_results.findtext("FitQuality"))
        energy_fit = fit_results.find("FitParameter0")
        energy_obs_str = energy_fit.findtext("MCObservable/Info")
        pattern = r"^(?P<obsname>\S+) (?P<obsid>\d+) (?P<simple>s|n) (?P<complex_arg>re|im)$"
        match = regex.match(pattern, energy_obs_str.strip())
        if match.group('simple') != 'n' or match.group('complex_arg') != 're':
          logging.error("Energies are supposed to be simple and real")

        energy_obs = sigmond.MCObsInfo(match.group('obsname'), int(match.group('obsid')))
        fit_info = FitInfo.createFromObservable(energy_obs)
        energy_value = float(energy_fit.findtext("MCEstimate/FullEstimate"))
        energy_error = float(energy_fit.findtext("MCEstimate/SymmetricError"))
        energy = util.nice_value(energy_value, energy_error)
        amplitude_fit = fit_results.find("FitParameter1")
        amplitude_value = float(amplitude_fit.findtext("MCEstimate/FullEstimate"))
        amplitude_error = float(amplitude_fit.findtext("MCEstimate/SymmetricError"))
        amplitude = util.nice_value(amplitude_value, amplitude_error)
        gap = "---"
        if fit_info.has_gap:
          sqrt_gap_fit = fit_results.find("FitParameter2")
          sqrt_gap_value = float(sqrt_gap_fit.findtext("MCEstimate/FullEstimate"))
          gap_value = sqrt_gap_value**2
          gap_error = 2.*abs(sqrt_gap_value)*float(sqrt_gap_fit.findtext("MCEstimate/SymmetricError"))
          gap = util.nice_value(gap_value, gap_error)

        const = '---'
        if fit_info.has_const:
          const_fit_num = fit_info.num_params - 1
          const_fit = fit_results.find(f"FitParameter{const_fit_num}")
          const_value = float(const_fit.findtext("MCEstimate/FullEstimate"))
          const_err = float(const_fit.findtext("MCEstimate/SymmetricError"))
          const = util.nice_value(const_value, const_err)
        
        if fit_info.ratio: #changing all output to be in the same lab frame units
            reconstructed_energy = energy
            for task_xml2 in log_xml_root.findall("Task"):
                reconstructed_energy_xml = task_xml2.find("DoObsFunction")
                if reconstructed_energy_xml is None or reconstructed_energy_xml.find("Error") is not None:
                    continue
                if reconstructed_energy_xml.findtext("Type") != "ReconstructEnergy":
                    continue
                if energy_obs_str != reconstructed_energy_xml.findtext("EnergyDifference/MCObservable/Info"):
                    continue
                reconstructed_value = float(reconstructed_energy_xml.findtext("MCEstimate/FullEstimate"))
                reconstructed_err = float(reconstructed_energy_xml.findtext("MCEstimate/SymmetricError"))
                reconstructed_energy = util.nice_value(reconstructed_value,reconstructed_err)
                break

            fit_result = FitResult(chisq_dof, quality, energy, reconstructed_energy, amplitude, gap, const, cov_cond)
        else:
            fit_result = FitResult(chisq_dof, quality, energy, energy, amplitude, gap, const, cov_cond)

        if fit_info in self.fits:
          logging.warning(f"Found two identical fits in {self.logfile}, ignoring...")
          continue

        self.fits[fit_info] = fit_result

      except AttributeError as err:
        logging.warning(f"{err} in DoFit task {count} in {self.logfile}")


class Level(NamedTuple):
  new: int
  original: int

class SpectrumLog(SigmondLog):

  def parse(self, log_xml_root):
    self.energies = SortedDict()
    self.reorder = True
    energy_level_xmls = log_xml_root.findall("Task/DoRotCorrMatInsertFitInfos/SinglePivot/ReorderEnergies/EnergyLevel")
    if not energy_level_xmls:
      self.reorder = False
      energy_level_xmls = log_xml_root.findall("Task/GetFromPivot/Energies/EnergyLevel")

    for energy_level_xml in energy_level_xmls:
      if self.reorder:
        new_level = int(energy_level_xml.findtext("LevelIndex"))
        orig_level = int(energy_level_xml.findtext("OriginalIndex"))
      else:
        new_level = int(energy_level_xml.findtext("Level"))
        orig_level = int(energy_level_xml.findtext("Level"))

      level = Level(new_level, orig_level)
      energy_obs_str = energy_level_xml.findtext("MCObservable/Info")
      pattern = r"^(?P<obsname>\S+) (?P<obsid>\d+) (?P<simple>s|n) (?P<complex_arg>re|im)$"
      match = regex.match(pattern, energy_obs_str.strip())
      if match.group('simple') != 'n' or match.group('complex_arg') != 're':
        logging.error("Energies are supposed to be simple and real")

      energy_obs = sigmond.MCObsInfo(match.group('obsname'), int(match.group('obsid')))
      fit_info = FitInfo.createFromObservable(energy_obs)
      self.energies[level] = fit_info

  
    self.zfactors = SortedDict()
    for operator_zfactor_xml in log_xml_root.findall(
        "Task/DoCorrMatrixZMagSquares/OperatorZMagnitudeSquares"):
      op_zfactors = SortedDict()
      op_str = [xml.text for xml in operator_zfactor_xml.iter()
                if 'OperatorString' in xml.tag][0]

      op = Operator(op_str)

      for zfactor_xml in operator_zfactor_xml.findall("ZMagSquare"):
        level = int(zfactor_xml.findtext("Level"))
        zfactor_value = float(zfactor_xml.findtext("Value/MCEstimate/FullEstimate"))
        op_zfactors[level] = zfactor_value

      self.zfactors[op] = op_zfactors

