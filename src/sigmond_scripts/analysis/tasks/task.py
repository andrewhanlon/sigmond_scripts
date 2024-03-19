import os
from abc import ABCMeta, abstractmethod
import pylatex
from typing import NamedTuple
from enum import Enum

import operator_info.operator
import data_handling.data_handler
import sigmond_info.sigmond_input
import utils.util as util
import sigmond

class Task(metaclass=ABCMeta):
  """The Task abstract base class

  All tasks in run-sigmond derive from this base class
  """

  relative_logdir = ".sigmond/logs"
  relative_inputdir = ".sigmond/inputs"
  relative_plotdir = ".sigmond/plots"

  def __init__(self, project_info, task_name):
    """The __init__ method for Task

    Args:
      project_info (launcher.ProjectInfo): This is a named tuple
          which holds all the needed information about the current
          project.
    """
    self.project_info = project_info
    self.task_name = task_name
    self.data_handler = data_handling.data_handler.DataHandler(project_info)

  @abstractmethod
  def readConfig(self, **task_options):
    pass

  @abstractmethod
  def getSigmondInputs(self):
    pass

  @abstractmethod
  def finalize(self):
    pass

  @property
  def results_dir(self):
    results_dir = os.path.join(self.project_dir, self.task_type)
    os.makedirs(results_dir, exist_ok=True)
    return results_dir

  def project_name(self, name):
    return f"{self.task_type}_{self.task_name}_{name}"

  @property
  def logdir(self):
    logdir = os.path.join(self.project_dir, self.relative_logdir, self.task_type)
    os.makedirs(logdir, exist_ok=True)
    return logdir

  def logfile(self, name):
    logdir = self.logdir
    logfile = f"{self.task_name}.{name}.log"
    return os.path.join(logdir, logfile)

  @property
  def inputdir(self):
    inputdir = os.path.join(self.project_dir, self.relative_inputdir, self.task_type)
    os.makedirs(inputdir, exist_ok=True)
    return inputdir

  def inputfile(self, name):
    inputdir = self.inputdir
    inputfile = f"{self.task_name}.{name}.input"
    return os.path.join(inputdir, inputfile)

  @property
  def plotdir(self):
    plotdir = os.path.join(self.project_dir, self.relative_plotdir, self.task_type)
    os.makedirs(plotdir, exist_ok=True)
    return plotdir
  
  @property
  def ensemble_name(self):
    return self.bins_info.getMCEnsembleInfo().getId().split('|')[0]

  @property
  def ensemble_time_extent(self):
    return self.bins_info.getLatticeTimeExtent()

  @property
  def ensemble_spatial_extent(self):
    x_extent = self.bins_info.getLatticeXExtent()
    y_extent = self.bins_info.getLatticeYExtent()
    z_extent = self.bins_info.getLatticeZExtent()

    if x_extent != y_extent or x_extent != z_extent:
      logging.error("Spatial extents of lattice not all equal")
      return (x_extent, y_extent, z_extent)

    return x_extent

  @property
  def project_dir(self):
    return self.project_info.project_dir

  @property
  def raw_data_dirs(self):
    return self.project_info.raw_data_dirs

  @property
  def ensembles_file(self):
    return self.project_info.ensembles_file

  @property
  def echo_xml(self):
    return self.project_info.echo_xml

  @property
  def bins_info(self):
    return self.project_info.bins_info

  @property
  def rebin(self):
    return self.bins_info.getRebinFactor()

  @property
  def sampling_info(self):
    return self.project_info.sampling_info

  @property
  def sampling_mode(self):
    return self.sampling_info.getSamplingMode()

  @property
  def precompute(self):
    return self.project_info.precompute

  @property
  def data_files(self):
    return self.project_info.data_files

  @property
  def latex_compiler(self):
    return self.project_info.latex_compiler

  def new_sigmond_input(self, project_name, inputfile, logfile, data_files):
    return sigmond_info.sigmond_input.SigmondInput(
        project_name, self.bins_info, self.sampling_info, self.ensembles_file, data_files,
        inputfile, logfile, self.precompute, self.echo_xml)


  def insertSigmondPlotTasks(self, sigmond_input, name, operators): 
    for op_src in operators:
      src_channel = op_src.channel
      subtractvev = self.subtractvev and src_channel.vev

      for op_snk in operators:
        snk_channel = op_snk.channel

        if src_channel != snk_channel:
          continue

        corr = sigmond.CorrelatorInfo(op_snk.operator_info, op_src.operator_info)
        '''
        The plot task will simply fail if the data isn't found
        if not self.data_handler.hasCorrelator(corr):
          continue
        '''

        if corr.isSinkSourceSame():
          corr_plotfile = self.correlator_plotfile(corr, name)
          energy_plotfile = self.energy_plotfile(op_src, name)

          sigmond_input.doCorrelatorPlot(
              corr, sigmond.ComplexArg.RealPart, corr_plotfile, hermitian=self.hermitian,
              subtractvev=subtractvev, corrname=self.plot_info.corrname,
              symbol_color=self.plot_info.symbol_color, symbol_type=self.plot_info.symbol_type,
              rescale=self.plot_info.rescale)
          sigmond_input.doEnergyPlot(
              corr, sigmond.ComplexArg.RealPart, energy_plotfile,
              eff_energy_type=self.plot_info.eff_energy_type, timestep=self.plot_info.timestep,
              hermitian=self.hermitian, subtractvev=subtractvev, corrname=self.plot_info.corrname,
              symbol_color=self.plot_info.symbol_color, symbol_type=self.plot_info.symbol_type)

        elif self.off_diagonal:
          re_corr_plotfile = self.correlator_plotfile(
              corr, name, complex_arg=sigmond.ComplexArg.RealPart)
          sigmond_input.doCorrelatorPlot(
              corr, sigmond.ComplexArg.RealPart, re_corr_plotfile, hermitian=self.hermitian,
              subtractvev=subtractvev, corrname=self.plot_info.corrname,
              symbol_color=self.plot_info.symbol_color, symbol_type=self.plot_info.symbol_type,
              rescale=self.plot_info.rescale)

          im_corr_plotfile = self.correlator_plotfile(
              corr, name, complex_arg=sigmond.ComplexArg.ImaginaryPart)
          sigmond_input.doCorrelatorPlot(
              corr, sigmond.ComplexArg.ImaginaryPart, im_corr_plotfile, hermitian=self.hermitian,
              subtractvev=subtractvev, corrname=self.plot_info.corrname,
              symbol_color=self.plot_info.symbol_color, symbol_type=self.plot_info.symbol_type,
              rescale=self.plot_info.rescale)


  def addPlotsToPDF(self, doc, data_files, operators, name):
    obs_handler, _ = util.get_obs_handlers(data_files, self.bins_info, self.sampling_info)

    corr_plotsdir = self.correlator_plotdir(name)
    energy_plotsdir = self.energy_plotdir(name)
    util.dirGrace2pdf(corr_plotsdir)
    util.dirGrace2pdf(energy_plotsdir)

    off_diag_corrs = list()

    for op_src in operators:
      for op_snk in operators:
        if op_src == op_snk:
          continue

        corr = sigmond.CorrelatorInfo(op_snk.operator_info, op_src.operator_info)
        if not self.data_handler.hasCorrelator(corr):
          continue

        off_diag_corrs.append(corr)

    with doc.create(pylatex.Subsection("Diagonal Correlators")):
      for operator in operators:
        corr = sigmond.CorrelatorInfo(operator.operator_info, operator.operator_info)
        if self.data_handler.hasCorrelator(corr):
          with doc.create(pylatex.Subsubsection(str(operator))):
            util.add_correlator(doc, self, corr, name, obs_handler)

    if self.off_diagonal and off_diag_corrs:
      with doc.create(pylatex.Subsection("Off-Diagonal Correlators")):
        for corr in off_diag_corrs:
          with doc.create(pylatex.Subsubsection(corr.corr_str())):
            util.add_correlator(doc, self, corr, name, obs_handler)

  def correlator_plotdir(self, name):
    plotdir = os.path.join(self.project_dir, self.relative_plotdir, f"{self.task_type}_correlators", name)
    os.makedirs(plotdir, exist_ok=True)
    return plotdir

  def correlator_plotfile(self, correlator, name, **options):
    """get a correlator plotfile
      
      correlator (sigmond.CorrelatorInfo): the correlator to plot
      name (str): ...
      **complex_arg (ComplexArg): Required if the correlator
          is not diagonal. If it is diagonal, this is
          ignored
      **extension (PlotExtension): 
    """
    extension = options.pop('extension', util.PlotExtension.grace)
    complex_arg = options.pop('complex_arg', sigmond.ComplexArg.RealPart)

    plotdir = self.correlator_plotdir(name)
    if correlator.isSinkSourceSame():
      opsrc = operator_info.operator.Operator(correlator.getSource())
      plotfile = f"corr_{opsrc!r}.{extension.value}"
    else:
      complex_arg = "re" if complex_arg == sigmond.ComplexArg.RealPart else "im"
      corr_name = util.str_to_file(correlator.corr_str())
      plotfile = f"corr_{corr_name}_{complex_arg}.{extension.value}"

    return os.path.join(plotdir, plotfile)
  
  def energy_plotdir(self, name):
    plotdir = os.path.join(self.project_dir, self.relative_plotdir, f"{self.task_type}_eff_energies", name)
    os.makedirs(plotdir, exist_ok=True)
    return plotdir

  def energy_plotfile(self, operator, name, extension=util.PlotExtension.grace):
    plotdir = self.energy_plotdir(name)
    plotfile = f"eff-energy_{operator!r}.{extension.value}"
    return os.path.join(plotdir, plotfile)


