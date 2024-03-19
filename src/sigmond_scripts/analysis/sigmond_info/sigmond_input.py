"""sigmond_input module

This module defines the SigmondInput class which is used for
creating SigMonD input XML files.
"""

import logging
import xml.etree.ElementTree as ET
import copy
import xml.dom.minidom as minidom
import numpy as np

import sigmond_scripts.analysis.sigmond_info.sigmond_info as sigmond_info
import sigmond_scripts.analysis.sigmond_info.fit_info as fit_info_lib
import sigmond_scripts.analysis.utils.util as util

import sigmond


class SigmondInput:
  """The SigmondInput class 
  
  The __init__ method creates the <Initialize> section of the XML
  file, and tasks are added by calling the various methods starting
  with 'do'. When you are finished creating your sigmond input XML
  file, you can call the 'writeInput' method to write the input file
  to disk.

  For more details on these XML files, see the SigMonD documentation.
  """

  def __init__(self, project_name, bins_info, sampling_info, ensembles_file, data_files,
               filename, logfile, precompute, echo_xml):
    """The __init__ method for SigmondInput.

    This method creates the <Initialize> section of the SigMonD input
    XML file. It also creates an empty <TaskSequence> section.

    Args:
      project_name (str): A string that denotes the ProjectName for
          the job.  bins_info (sigmond_xml.BinsInfo): An object of
          type BinsInfo that describes the details of the gauge
          configurations being used.
      sampling_info (sigmondbind.MCSamplingInfo): An object of type
          SamplingInfo that describes the type of resampling to use.
      ensembles_file (str): A filename containing all known ensembles
      data_files (data_handler.DataFiles): An object of type DataFiles
          that specifies the data files that are required in the
          sigmond run.
      filename (str): The filename where the input will be written.
      logfile (str): The filename where the log will be written.
      precompute (bool): Specifies whether the bootstrap samples
          should be precomputed.
    """
    self.filename = filename

    self.sigmondXML = ET.Element("SigMonD")
    init_tag = ET.SubElement(self.sigmondXML, "Initialize")
    ET.SubElement(init_tag, "ProjectName").text = project_name
    ET.SubElement(init_tag, "LogFile").text = logfile
    if echo_xml:
      ET.SubElement(init_tag, "EchoXML")
    if ensembles_file:
      ET.SubElement(init_tag, "KnownEnsemblesFile").text = ensembles_file
    init_tag.append(bins_info.xml())
    sampling_info_xml = sampling_info.xml()
    if precompute and sampling_info.isBootstrapMode():
      bootstrap_xml = sampling_info_xml.find("Bootstrapper")
      ET.SubElement(sampling_info_xml, "Precompute")
    init_tag.append(sampling_info_xml)
    init_tag.append(data_files.xml())

    self.taskSequence_tag = ET.SubElement(self.sigmondXML, "TaskSequence")

  def _addTask(self, task):
    self.taskSequence_tag.append(task)

  def write(self):
    """
    Writes the input XML file to disk
    """
    xmlstr = minidom.parseString(
        ET.tostring(self.sigmondXML, 'utf-8')).toprettyxml(indent="  ")
    xmllines = xmlstr.split('\n')
    with open(self.filename, 'w') as f:
      for xmlline in xmllines:
        xmlline = xmlline.rstrip()
        if xmlline != '':
          f.write('{}\n'.format(xmlline))

    logging.info("Wrote XML to: {}".format(self.filename))

  def to_str(self):
    return minidom.parseString(
      ET.tostring(self.sigmondXML, 'utf-8')).toprettyxml(indent="  ")

  def doCorrelatorCheck(self, operators, mintime, maxtime, **extra_options):
    """Adds a 'DoChecks' task of type 'TemporalCorrelatorMatrix'
       to the SimgondInput object.

    Args:
      operators (list): A list of objects of type 'Operator' that are
          to form a correlator matrix to check
      mintime (int): the mininum time slice considered for the
          correlators
      maxtime (int): the maximum time slice considered for the
          correlators
      **outlier (float): the outlier parameter as described in
          the SigMonD documenation
      **hermitian (bool): specifies whether the correlator matrix
          should be assumed to be Hermitian.
      **subtractvev (bool): specifies whether the correlators should
          subtract their VEVs (if non-zero VEVs exist).
      **verbose (bool): specifies whether the output to the logfile
          from sigmond should add more information.
    """

    xml = ET.Element("Task")
    ET.SubElement(xml, "Action").text = "DoChecks"
    ET.SubElement(xml, "Type").text = "TemporalCorrelatorMatrix"
    correlator_matrix_tag = ET.SubElement(xml, "CorrelatorMatrixInfo")
    for operator in operators:
      correlator_matrix_tag.append(operator.xml())

    if extra_options.get('hermitian'):
      ET.SubElement(correlator_matrix_tag, "HermitianMatrix")
    if extra_options.get('subtractvev'):
      ET.SubElement(correlator_matrix_tag, "SubtractVEV")

    ET.SubElement(xml, "MinTimeSep").text = str(mintime)
    ET.SubElement(xml, "MaxTimeSep").text = str(maxtime)

    if extra_options.get('verbose'):
      ET.SubElement(xml, "Verbose")

    if 'outlier' in extra_options:
      ET.SubElement(xml, "OutlierScale").text = str(extra_options['outlier'])

    self._addTask(xml)

  def doHermitianCheck(self, operators, mintime, maxtime, **extra_options):
    """Adds a 'DoChecks' task of type
       'TemporalCorrelatorMatrixIsHermitian' to the SigmondInput
       object.

    Args:
      operators (list): A list of objects of type 'Operator' that
          are to form a correlator matrix to be checked for
          Hermiticity.
      mintime (int): the mininum time slice considered for the
          correlators
      maxtime (int): the maximum time slice considered for the
          correlators
      **verbose (bool): specifies whether the output to the logfile
          from sigmond should add more information.
    """

    xml = ET.Element("Task")
    ET.SubElement(xml, "Action").text = "DoChecks"
    ET.SubElement(xml, "Type").text = "TemporalCorrelatorMatrixIsHermitian"
    correlator_matrix_tag = ET.SubElement(xml, "CorrelatorMatrixInfo")
    for operator in operators:
      correlator_matrix_tag.append(operator.xml())

    ET.SubElement(xml, "MinTimeSep").text = str(mintime)
    ET.SubElement(xml, "MaxTimeSep").text = str(maxtime)

    if extra_options.get('verbose'):
      ET.SubElement(xml, "Verbose")

    self._addTask(xml)

  def doCorrelatorPrint(self, correlator, arg, **extra_options):
    """Adds a 'PrintXML' task of type 'TemporalCorrelator' to the
       SigmondInput object

    Args:
      correlator (sigmondbind.Correlator): The correlator to be printed.
      arg (sigmondbind.ComplexArg): Specifies whether the real or
          imaginary part should be printed.
      **hermitian (bool): specifies whether the correlator matrices
          should be assumed to be Hermitian.
      **subtractvev (bool): specifies whether the correlators should
          subtract their VEVs (if non-zero VEVs exist).
      **sampling_mode (sigmondbind.SamplingMode): specifies a sampling
          mode to use
    """

    xml = ET.Element("Task")
    ET.SubElement(xml, "Action").text = "PrintXML"
    ET.SubElement(xml, "Type").text = "TemporalCorrelator"
    xml.append(correlator.xml())
    ET.SubElement(xml, "Arg").text = str(arg)
    if extra_options.get('hermitian'):
      ET.SubElement(xml, "HermitianMatrix")
    if extra_options.get('subtractvev'):
      ET.SubElement(xml, "SubtractVEV")
    if 'sampling_mode' in extra_options:
      ET.SubElement(xml, "SamplingMode").text = str(extra_options['sampling_mode'])

    self._addTask(xml)

  def doEnergyPrint(self, correlator, arg, **extra_options):
    """Adds a 'PrintXML' task of type 'EffectiveEnergy' to the
        SigmondInput object

    Args:
      correlator (sigmondbind.Correlator): The correlator to be printed.
      arg (sigmondbind.ComplexArg): Specifies whether the real or
          imaginary part should be printed.
      **eff_energy_type (EffEnergyType): Specifies the type of effective
          energy to compute.
      **timestep (int): Specifies the time step to be used to compute
          the effective energy.
      **hermitian (bool): specifies whether the correlator matrices
          should be assumed to be Hermitian.
      **subtractvev (bool): specifies whether the correlators should
          subtract their VEVs (if non-zero VEVs exist).
      **sampling_mode (sigmond_xml.SamplingMode): specifies a sampling
          mode to use
    """

    xml = ET.Element("Task")
    ET.SubElement(xml, "Action").text = "PrintXML"
    ET.SubElement(xml, "Type").text = "EffectiveEnergy"
    if 'eff_energy_type' in extra_options:
      ET.SubElement(xml, "EffEnergyType").text = extra_options['eff_energy_type'].name
    if 'timestep' in extra_options:
      ET.SubElement(xml, "TimeStep").text = str(extra_options['timestep'])
    xml.append(correlator.xml())
    ET.SubElement(xml, "Arg").text = str(arg)
    if extra_options.get('hermitian'):
      ET.SubElement(xml, "HermitianMatrix")
    if extra_options.get('subtractvev'):
      ET.SubElement(xml, "SubtractVEV")
    if 'sampling_mode' in extra_options:
      ET.SubElement(xml, "SamplingMode").text = str(extra_options['sampling_mode'])

    self._addTask(xml)

  def doCorrelatorPlot(self, correlator, arg, plotfile, **extra_options):
    """Adds a 'DoPlot' task of type 'TemporalCorrelator' to the
        SigmondInput object

    Args:
      correlator (sigmondbind.Correlator): The correlator to be printed.
      arg (sigmondbind.ComplexArg): Specifies whether the real or
      imaginary part should be printed.
      plotfile (str): Specifies the file to store the plot in.
      **hermitian (bool): specifies whether the correlator matrices
          should be assumed to be Hermitian.
      **subtractvev (bool): specifies whether the correlators should
          subtract their VEVs (if non-zero VEVs exist).
      **sampling_mode (sigmond_xml.SamplingMode): specifies a sampling
          mode to use
      **corrname (str): specifies the name to print on the plot.
          'standard' produces the standard print.
      **symbol_color (SymbolColor): specifies the color to use for the data
          points.
      **symbol_type (SymbolType): specifies the shape of the data points to use.
      **rescale (float): specifies a rescaling factors for the correlator.
    """

    xml = ET.Element("Task")
    ET.SubElement(xml, "Action").text = "DoPlot"
    ET.SubElement(xml, "Type").text = "TemporalCorrelator"
    xml.append(correlator.xml())
    ET.SubElement(xml, "Arg").text = str(arg)
    if extra_options.get('hermitian'):
      ET.SubElement(xml, "HermitianMatrix")
    if extra_options.get('subtractvev'):
      ET.SubElement(xml, "SubtractVEV")
    if 'sampling_mode' in extra_options:
      ET.SubElement(xml, "SamplingMode").text = str(extra_options['sampling_mode'])

    ET.SubElement(xml, "PlotFile").text = plotfile
    if 'corrname' in extra_options:
      ET.SubElement(xml, "CorrName").text = extra_options['corrname']
    if 'symbol_color' in extra_options:
      ET.SubElement(xml, "SymbolColor").text = extra_options['symbol_color'].value
    if 'symbol_type' in extra_options:
      ET.SubElement(xml, "SymbolType").text = extra_options['symbol_type'].value
    if 'rescale' in extra_options:
      ET.SubElement(xml, "Rescale").text = str(extra_options['rescale'])

    self._addTask(xml)

  def doEnergyPlot(self, correlator, arg, plotfile, **extra_options):
    """Adds a 'DoPlot' task of type 'EffectiveEnergy' to the
        SigmondInput object

    Args:
      correlator (sigmondbind.Correlator): The correlator to be printed.
      arg (sigmondbind.ComplexArg): Specifies whether the real or
          imaginary part should be printed.
      plotfile (str): Specifies the file to store the plot in.
      **eff_energy_type (EffEnergyType): Specifies the type of
          effective energy to compute.
      **timestep (int): Specifies the time step to be used to compute
          the effective energy.
      **hermitian (bool): specifies whether the correlator matrices
          should be assumed to be Hermitian.
      **subtractvev (bool): specifies whether the correlators should
          subtract their VEVs (if non-zero VEVs exist).
      **sampling_mode (sigmond_xml.SamplingMode): specifies a sampling
          mode to use
      **corrname (str): specifies the name to print on the plot.
          'standard' produces the standard print.
      **symbol_color (SymbolColor): specifies the color to use for the data points.
      **symbol_type (SymbolType): specifies the shape of the data points to use.
      **max_error (float): specifies the max error to plot
    """

    xml = ET.Element("Task")
    ET.SubElement(xml, "Action").text = "DoPlot"
    ET.SubElement(xml, "Type").text = "EffectiveEnergy"
    if 'eff_energy_type' in extra_options:
      ET.SubElement(xml, "EffEnergyType").text = extra_options['eff_energy_type'].name
    if 'timestep' in extra_options:
      ET.SubElement(xml, "TimeStep").text = str(extra_options['timestep'])
    xml.append(correlator.xml())
    ET.SubElement(xml, "Arg").text = str(arg)
    if extra_options.get('hermitian'):
      ET.SubElement(xml, "HermitianMatrix")
    if extra_options.get('subtractvev'):
      ET.SubElement(xml, "SubtractVEV")
    if 'sampling_mode' in extra_options:
      ET.SubElement(xml, "SamplingMode").text = str(extra_options['sampling_mode'])

    ET.SubElement(xml, "PlotFile").text = plotfile
    if 'corrname' in extra_options:
      ET.SubElement(xml, "CorrName").text = extra_options['corrname']
    if 'symbol_color' in extra_options:
      ET.SubElement(xml, "SymbolColor").text = extra_options['symbol_color'].value
    if 'symbol_type' in extra_options:
      ET.SubElement(xml, "SymbolType").text = extra_options['symbol_type'].value
    if 'max_error' in extra_options:
      ET.SubElement(xml, "MaxErrorToPlot").text = str(extra_options['max_error'])

    self._addTask(xml)

  def doAnisotropyFromDispersion(self, scattering_particle_energies, spatial_extent, **extra_options):
    """ Adds a 'DoAnisotropyFromDispersion' task to the
        SigmondInput object

    Args:
      scattering_particle_energies (dict): Specifies the single
          hadron fits
      spatial_extent (int): Need this if you want to get the
          absolute/difference energy from a tmin fit
      **name (str): used for naming the fit parameters
      **minimizer_info (sigmond_xml.MinimizerInfo): Specifies the
          info for the Minimizer.
      **sampling_mode (sigmond_xml.SamplingMode): specifies a
          sampling mode to use for the fit.
      **cov_sampling_mode (sigmond_xml.SamplingMode): specifies the
          sampling mode to use in the calculation of the covariance
          matrix.
      **plotfile (str): specifies the plot filename
      **plot_info (AnisotropyPlotInfo): all the information
          for the plotting.
    """

    xml = ET.Element("Task")
    ET.SubElement(xml, "Action").text = "DoFit"
    ET.SubElement(xml, "Type").text = "AnisotropyFromDispersion"
    if 'minimizer_info' in extra_options:
      xml.append(extra_options['minimizer_info'].xml())
    if 'sampling_mode' in extra_options:
      ET.SubElement(xml, "SamplingMode").text = str(extra_options['sampling_mode'])
    if 'cov_sampling_mode' in extra_options:
      ET.SubElement(xml, "CovMatCalcSamplingMode").text = str(extra_options['cov_sampling_mode'])

    fit_tag = ET.SubElement(xml, "AnisotropyFromDispersionFit")
    ET.SubElement(fit_tag, "SpatialExtentNumSites").text = str(spatial_extent)
    for scat_energy_obs, psq in scattering_particle_energies:
      scat_energy_xml = ET.SubElement(fit_tag, "Energy")
      ET.SubElement(scat_energy_xml, "Name").text = scat_energy_obs.getObsName()
      ET.SubElement(scat_energy_xml, "IDIndex").text = str(scat_energy_obs.getObsIndex())
      ET.SubElement(scat_energy_xml, "IntMomSquared").text = str(psq)

    name = extra_options.get('name', '')
    aniso_xml = ET.SubElement(fit_tag, "Anisotropy")
    ET.SubElement(aniso_xml, "Name").text = f"{name}_aniso"
    ET.SubElement(aniso_xml, "IDIndex").text = "0"
    mass_xml = ET.SubElement(fit_tag, "RestMassSquared")
    ET.SubElement(mass_xml, "Name").text = f"{name}_mass_squared"
    ET.SubElement(mass_xml, "IDIndex").text = "0"

    if 'plotfile' in extra_options:
      plot_tag = ET.SubElement(fit_tag, "DoPlot")
      ET.SubElement(plot_tag, "PlotFile").text = extra_options['plotfile']
      if 'name' in extra_options:
        ET.SubElement(plot_tag, "ParticleName").text = extra_options['name']
      if 'plot_info' in extra_options:
        plot_info = extra_options['plot_info']
        ET.SubElement(plot_tag, "SymbolColor").text = plot_info.symbol_color.value
        ET.SubElement(plot_tag, "SymbolType").text = plot_info.symbol_type.value
        ET.SubElement(plot_tag, "Goodness").text = plot_info.goodness.value

    self._addTask(xml)

  def doDispersion(self, scattering_particle_energies, spatial_extent, **extra_options):
    """ Adds a 'DoDispersion' task to the
        SigmondInput object

    Args:
      scattering_particle_energies (dict): Specifies the single
          hadron fits
      spatial_extent (int): Need this if you want to get the
          absolute/difference energy from a tmin fit
      **name (str): used for naming the fit parameters
      **minimizer_info (sigmond_xml.MinimizerInfo): Specifies the
          info for the Minimizer.
      **sampling_mode (sigmond_xml.SamplingMode): specifies a
          sampling mode to use for the fit.
      **cov_sampling_mode (sigmond_xml.SamplingMode): specifies the
          sampling mode to use in the calculation of the covariance
          matrix.
      **plotfile (str): specifies the plot filename
      **plot_info (DispersionPlotInfo): all the information
          for the plotting.
    """

    xml = ET.Element("Task")
    ET.SubElement(xml, "Action").text = "DoFit"
    ET.SubElement(xml, "Type").text = "Dispersion"
    if 'minimizer_info' in extra_options:
      xml.append(extra_options['minimizer_info'].xml())
    if 'sampling_mode' in extra_options:
      ET.SubElement(xml, "SamplingMode").text = str(extra_options['sampling_mode'])
    if 'cov_sampling_mode' in extra_options:
      ET.SubElement(xml, "CovMatCalcSamplingMode").text = str(extra_options['cov_sampling_mode'])

    fit_tag = ET.SubElement(xml, "DispersionFit")
    ET.SubElement(fit_tag, "SpatialExtentNumSites").text = str(spatial_extent)
    for scat_energy_obs, psq in scattering_particle_energies:
      scat_energy_xml = ET.SubElement(fit_tag, "Energy")
      ET.SubElement(scat_energy_xml, "Name").text = scat_energy_obs.getObsName()
      ET.SubElement(scat_energy_xml, "IDIndex").text = str(scat_energy_obs.getObsIndex())
      ET.SubElement(scat_energy_xml, "IntMomSquared").text = str(psq)

    name = extra_options.get('name', '')
    coeff_xml = ET.SubElement(fit_tag, "Coefficient")
    ET.SubElement(coeff_xml, "Name").text = f"{name}_coeff"
    ET.SubElement(coeff_xml, "IDIndex").text = "0"
    mass_xml = ET.SubElement(fit_tag, "RestMassSquared")
    ET.SubElement(mass_xml, "Name").text = f"{name}_mass_squared"
    ET.SubElement(mass_xml, "IDIndex").text = "0"

    if 'plotfile' in extra_options:
      plot_tag = ET.SubElement(fit_tag, "DoPlot")
      ET.SubElement(plot_tag, "PlotFile").text = extra_options['plotfile']
      if 'name' in extra_options:
        ET.SubElement(plot_tag, "ParticleName").text = extra_options['name']
      if 'plot_info' in extra_options:
        plot_info = extra_options['plot_info']
        ET.SubElement(plot_tag, "SymbolColor").text = plot_info.symbol_color.value
        ET.SubElement(plot_tag, "SymbolType").text = plot_info.symbol_type.value
        ET.SubElement(plot_tag, "Goodness").text = plot_info.goodness.value

    self._addTask(xml)

  def doTemporalCorrelatorFit(self, fit_info, **extra_options):
    """ Adds a 'DoFit' task of type 'TemporalCorrelator' to the
        SigmondInput object

    Args:
      fit_info (sigmond_info.FitInfo): Specifies what is
          being fit
      **minimizer_info (sigmond_xml.MinimizerInfo): Specifies the
          info for the Minimizer.
      **sampling_mode (sigmond_xml.SamplingMode): specifies a
          sampling mode to use for the fit.
      **cov_sampling_mode (sigmond_xml.SamplingMode): specifies the
          sampling mode to use in the calculation of the covariance
          matrix.
      **plotfile (str): specifies the plot filename
      **plot_info (FitPlotInfo or TMinPlotInfo): all the information
          for the plotting.
      **chosen_fit_info (MCObsInfo): The chosen fit info
      **reference_energy (sigmondbind.MCOsInfo): an observable
          that contains the reference energy.
      **spatial_extent (int): Need this if you want to get the
          absolute/difference energy from a tmin fit
      **non_interacting_level ([(MCObsInfo, psq)]): also need this
          if you want the absolute/difference energy from tmin fit
      **anisotropy (MCObsInfo): If wanting to do absolute/difference
          energy and have an anisotropy
    """

    conspiracy_fits = [
      fit_info_lib.FitModel.TwoExpConspiracy,
      # fit_info_lib.FitModel.TimeForwardTwoExponential,
      # fit_info_lib.FitModel.TimeForwardThreeExponential,
      # fit_info_lib.FitModel.TimeForwardFourExponential,
      # fit_info_lib.FitModel.TimeForwardThreeIndExponential,
    ]
    deg_conspiracy_fits = [
      fit_info_lib.FitModel.DegTwoExpConspiracy,
      fit_info_lib.FitModel.DegThreeExpConspiracy,
    ]

    xml = ET.Element("Task")
    ET.SubElement(xml, "Action").text = "DoFit"
    ET.SubElement(xml, "Type").text = fit_info.fit_type
    if 'minimizer_info' in extra_options:
      xml.append(extra_options['minimizer_info'].xml())
    if 'sampling_mode' in extra_options:
      ET.SubElement(xml, "SamplingMode").text = str(extra_options['sampling_mode'])
    if 'cov_sampling_mode' in extra_options:
      ET.SubElement(xml, "CovMatCalcSamplingMode").text = str(extra_options['cov_sampling_mode'])
    
    
    if 'plotfile' in extra_options:
      if fit_info.is_tmin_vary or fit_info.is_tmax_vary:
        plot_tag = ET.Element("PlotInfo")
        ET.SubElement(plot_tag, "PlotFile").text = extra_options['plotfile']
        if 'plot_info' in extra_options:
          plot_info = extra_options['plot_info']
          ET.SubElement(plot_tag, "CorrName").text = plot_info.obsname
          ET.SubElement(plot_tag, "SymbolType").text = plot_info.symbol_type.value
          ET.SubElement(plot_tag, "GoodFitSymbolColor").text = plot_info.goodfit_color.value
          ET.SubElement(plot_tag, "BadFitSymbolColor").text = plot_info.badfit_color.value
          if plot_info.correlatedfit_hollow:
            ET.SubElement(plot_tag, "CorrelatedFitSymbolHollow")
          if plot_info.uncorrelatedfit_hollow:
            ET.SubElement(plot_tag, "UncorrelatedFitSymbolHollow")
          ET.SubElement(plot_tag, "QualityThreshold").text = str(plot_info.quality_threshold)
          ET.SubElement(plot_tag, "CorrelatedThreshold").text = str(plot_info.correlated_threshold)

      else:
        plot_tag = ET.Element("DoEffectiveEnergyPlot")
        ET.SubElement(plot_tag, "PlotFile").text = extra_options['plotfile']
        if 'plot_info' in extra_options:
          plot_info = extra_options['plot_info']
          ET.SubElement(plot_tag, "CorrName").text = plot_info.corrname
          ET.SubElement(plot_tag, "TimeStep").text = str(plot_info.timestep)
          ET.SubElement(plot_tag, "SymbolColor").text = plot_info.symbol_color.value
          ET.SubElement(plot_tag, "SymbolType").text = plot_info.symbol_type.value
          ET.SubElement(plot_tag, "Goodness").text = plot_info.goodness.value
          ET.SubElement(plot_tag, "ShowApproach")
          if plot_info.max_relative_error > 0.:
            ET.SubElement(plot_tag, "MaxRelativeErrorToPlot").text = str(plot_info.max_relative_error)

        if 'reference_energy' in extra_options and extra_options['reference_energy'] is not None:
          reference_energy = extra_options['reference_energy']
          ref_energy_tag = ET.SubElement(plot_tag, "ReferenceEnergy")
          ET.SubElement(ref_energy_tag, "Name").text = reference_energy.getObsName()
          ET.SubElement(ref_energy_tag, "IDIndex").text = str(reference_energy.getObsIndex())

    
    if 'chosen_fit_info' in extra_options:
      chosen_fit_xml = ET.Element("ChosenFitInfo")
      ET.SubElement(chosen_fit_xml, "Name").text = extra_options['chosen_fit_info'].getObsName()
      ET.SubElement(chosen_fit_xml, "IDIndex").text = str(extra_options['chosen_fit_info'].getObsIndex())
            
    if fit_info.sim_fit: 
        tmin_min = fit_info.tmin
        tmin_max = fit_info.tmin_max
        tmax_min = fit_info.tmax_min
        tmax_max = fit_info.tmax
        
        fit_tag2 = ET.Element("Fits")
        
        if fit_info.is_tmin_vary:
            xml.find("Type").text = f"{fit_info.fit_type}TminVary"
            fit_tag = ET.Element(f"{fit_info.fit_type}TminVaryFit")
            ET.SubElement(fit_tag, "TminFirst").text = str(tmin_min)
            ET.SubElement(fit_tag, "TminLast").text = str(tmin_max)
            ET.SubElement(fit_tag, "Tmax").text = str(tmax_max)
            tmin_fit_tag = ET.SubElement(fit_tag, f"{fit_info.fit_type}Fit")
            tmin_fit_tag.append(fit_tag2)
        elif fit_info.is_tmax_vary:
            xml.find("Type").text = f"{fit_info.fit_type}TmaxVary"
            fit_tag = ET.Element(f"{fit_info.fit_type}TmaxVaryFit")
            ET.SubElement(fit_tag, "TmaxFirst").text = str(tmax_min)
            ET.SubElement(fit_tag, "TmaxLast").text = str(tmax_max)
            ET.SubElement(fit_tag, "Tmin").text = str(tmin_min)
            tmin_fit_tag = ET.SubElement(fit_tag, f"{fit_info.fit_type}Fit")
            tmin_fit_tag.append(fit_tag2)
        else:
            fit_tag = ET.Element(f"{fit_info.fit_type}Fit")
            fit_tag.append(fit_tag2)
        
        FinalEnergy = ET.SubElement(xml, "FinalEnergy")
        ET.SubElement(FinalEnergy, "Name").text = fit_info.obs_name
        ET.SubElement(FinalEnergy, "IDIndex").text = str(fit_info.obs_id(fit_info.energy_index))
        FinalAmplitude = ET.SubElement(xml, "FinalAmplitude")
        ET.SubElement(FinalAmplitude, "Name").text = fit_info.obs_name
        ET.SubElement(FinalAmplitude, "IDIndex").text = str(fit_info.obs_id(fit_info.amplitude_index))
        
        fit_info.sim_fit = False
        fit_info.tmin_max = -1
        fit_info.tmax_min = -1
        fit_xml = fit_info.xml()
        if 'chosen_fit_info' in extra_options:
            fit_xml.append(chosen_fit_xml)
        fit_info.tmin_max = tmin_max
        fit_info.tmax_min = tmax_min
        if 'plotfile' in extra_options:
            fit_xml.append(plot_tag)
        if (fit_info.is_tmin_vary or fit_info.is_tmax_vary):
            fit_xml.remove( fit_xml.find("MinimumTimeSeparation") )
            fit_xml.remove( fit_xml.find("MaximumTimeSeparation") )
        for name in fit_xml.findall("Model/*/Name"):
            name.text = fit_info.obs_name
                                           
        fit_tag2.append(fit_xml)
        if 'scattering_fit_info' in extra_options:
            sh_priors = extra_options['sh_priors']
            twothree = ["Second","Third"]
            scat_xmls= [copy.deepcopy(scat_fit_info.xml()) for scat_fit_info in extra_options['scattering_fit_info']]
            if len(scat_xmls)==2:
              for i, this_xml in enumerate(scat_xmls):
                  param = this_xml.find("Model/SqrtGapToSecondEnergy")
                  if i==0:
                    shift1 = sh_priors[param.find("Name").text][param.tag]["Mean"]
                  else:
                    shift1 = np.sqrt(sh_priors[param.find("Name").text][param.tag]["Mean"]**2-shift1*shift1)
              if np.isnan(shift1):
                scat_xmls.reverse()
                extra_options['scattering_fit_info'].reverse()

            for i, this_xml in enumerate(scat_xmls):
                # this_xml = copy.deepcopy(scat_fit_info.xml())
                if (fit_info.is_tmin_vary or fit_info.is_tmax_vary): #add shift energy if shift plot
                    chosen_name = this_xml.find('Model/FirstEnergy/Name').text
                    chosen_id = this_xml.find('Model/FirstEnergy/IDIndex').text
                    chosen_fit_xml = ET.SubElement( this_xml, "ChosenFitInfo")
                    ET.SubElement( chosen_fit_xml, "Name").text = chosen_name
                    ET.SubElement( chosen_fit_xml, "IDIndex").text = chosen_id

                new_priors_xml = ET.SubElement( this_xml, "Priors")
                for param in this_xml.findall('Model/*'):
                    if param.tag!="Type":
                      new_prior_xml = ET.SubElement( new_priors_xml, param.tag)
                      if param.tag=="SqrtGapToSecondEnergy":
                        if i==0:
                          shift1 = sh_priors[param.find("Name").text][param.tag]["Mean"]
                          err1 = sh_priors[param.find("Name").text][param.tag]["Error"]
                        else:
                          shift1 = np.sqrt(sh_priors[param.find("Name").text][param.tag]["Mean"]**2-shift1*shift1)
                          err1 = 0.8*shift1 #prior_width_multiplier*sh_priors[param.find("Name").text][param.tag]["Error"]+err1
                        ET.SubElement( new_prior_xml, "Mean").text = str(shift1)
                        ET.SubElement( new_prior_xml, "Error").text = str(err1)
                      else:
                        ET.SubElement( new_prior_xml, "Mean").text = str(sh_priors[param.find("Name").text][param.tag]["Mean"])
                        ET.SubElement( new_prior_xml, "Error").text = str(sh_priors[param.find("Name").text][param.tag]["Error"])
                      param.find("Name").text = fit_info.obs_name+"-"+param.find("Name").text
                     
                if (fit_info.is_tmin_vary or fit_info.is_tmax_vary):
                    ET.SubElement(this_xml, "Fixed")    

                if fit_info.model in conspiracy_fits: #insert other energy
                    ratio_name = f'TemporalCorrelatorFit/Model/SqrtGapTo{twothree[i]}Energy/'
                    fit_tag2.find(ratio_name+"Name").text=this_xml.find('Model/SqrtGapToSecondEnergy/Name').text
                    fit_tag2.find(ratio_name+"IDIndex").text=this_xml.find('Model/SqrtGapToSecondEnergy/IDIndex').text
                elif fit_info.model in deg_conspiracy_fits:
                  for i in range(deg_conspiracy_fits.index(fit_info.model)+1):
                    ratio_name = f'TemporalCorrelatorFit/Model/SqrtGapTo{twothree[i]}Energy/'
                    fit_tag2.find(ratio_name+"Name").text=this_xml.find(f'Model/SqrtGapTo{twothree[i]}Energy/Name').text
                    fit_tag2.find(ratio_name+"IDIndex").text=this_xml.find(f'Model/SqrtGapTo{twothree[i]}Energy/IDIndex').text
                    
                if fit_info.model==fit_info_lib.FitModel.TimeForwardThreeIndExponential:
                    ratio_name = f'TemporalCorrelatorFit/Model/SingleHadronEnergy/'
                    fit_tag2.find(ratio_name+"Name").text=this_xml.find('Model/FirstEnergy/Name').text
                    fit_tag2.find(ratio_name+"IDIndex").text=this_xml.find('Model/FirstEnergy/IDIndex').text

                if fit_info.model==fit_info_lib.FitModel.TimeForwardDoubleExpRatio1:
                    ratio_name = f'TemporalCorrelatorInteractionRatioFit/Model/SH{i+1}Gap/'
                    fit_tag2.find(ratio_name+"Name").text=this_xml.find('Model/SqrtGapToSecondEnergy/Name').text
                    fit_tag2.find(ratio_name+"IDIndex").text=this_xml.find('Model/SqrtGapToSecondEnergy/IDIndex').text
                    
                if fit_info.model==fit_info_lib.FitModel.TimeForwardDoubleExpRatio1:
                    ratio_name = f'TemporalCorrelatorInteractionRatioFit/Model/SHGap/'
                    fit_tag2.find(ratio_name+"Name").text=this_xml.find('Model/SqrtGapToSecondEnergy/Name').text
                    fit_tag2.find(ratio_name+"IDIndex").text=this_xml.find('Model/SqrtGapToSecondEnergy/IDIndex').text
                    # ratio_name = f'TemporalCorrelatorInteractionRatioFit/Model/SHGapAmp/'
                    # fit_tag2.find(ratio_name+"Name").text=this_xml.find('Model/SecondAmplitudeRatio/Name').text
                    # fit_tag2.find(ratio_name+"IDIndex").text=this_xml.find('Model/SecondAmplitudeRatio/IDIndex').text

                if fit_info.model==fit_info_lib.FitModel.TwoExpConspiracy and i==0:
                  shift_parameter = ET.Element("SqrtGapToSecondEnergyShift")
                  for item in this_xml.findall('Model/SqrtGapToSecondEnergy/*'):
                    ET.SubElement(shift_parameter, item.tag).text = item.text

                #2exp to shifted 2exp
                if fit_info.model==fit_info_lib.FitModel.TwoExpConspiracy and i==1:
                  this_xml.find('Model/Type').text = "TimeForwardTwoExponentialForCons"
                  this_xml.find('Model').append(shift_parameter)

                # if 'plotfile' in extra_options:
                #     ni_plot_tag = copy.deepcopy(plot_tag)
                #     ni_plot_tag.find('CorrName').text = ni_plot_tag.find('CorrName').text+f"-{scat_fit_info.operator.compact_str}"
                #     ni_plot_tag.find('PlotFile').text = ni_plot_tag.find('PlotFile').text.replace(".agr",f"-{scat_fit_info.operator.compact_str}.agr")
                #     this_xml.append(ni_plot_tag)

                fit_tag2.append(this_xml)
                if fit_info.model==fit_info_lib.FitModel.TimeForwardTwoExponential:
                  break
                if fit_info.model in deg_conspiracy_fits:
                  break
                if fit_info.model==fit_info_lib.FitModel.TimeForwardThreeIndExponential:
                  break
                        
        
        fit_info.sim_fit = True
        
    else:
        fit_tag = fit_info.xml()
        
    if (fit_info.is_tmin_vary or fit_info.is_tmax_vary) and 'spatial_extent' in extra_options and 'non_interacting_level' in extra_options:
      if fit_info.ratio:
        energy_xml = ET.SubElement(xml, "DoReconstructEnergy")
      else:
        energy_xml = ET.SubElement(xml, "DoEnergyDifference")

      ET.SubElement(energy_xml, "SpatialExtentNumSites").text = str(extra_options['spatial_extent'])
      if 'anisotropy' in extra_options:
        aniso_xml = ET.SubElement(energy_xml, "Anisotropy")
        ET.SubElement(aniso_xml, "Name").text = extra_options['anisotropy'].getObsName()
        ET.SubElement(aniso_xml, "IDIndex").text = str(extra_options['anisotropy'].getObsIndex())

      for scat_energy_obs, psq in extra_options['non_interacting_level']:
        scat_energy_xml = ET.SubElement(energy_xml, "ScatteringParticleEnergyFit")
        ET.SubElement(scat_energy_xml, "IntMomSquared").text = str(psq)
        ET.SubElement(scat_energy_xml, "Name").text = scat_energy_obs.getObsName()
        ET.SubElement(scat_energy_xml, "IDIndex").text = str(scat_energy_obs.getObsIndex())

    xml.append(fit_tag)
    
    
    if 'chosen_fit_info' in extra_options and not fit_info.sim_fit:
      xml.append(chosen_fit_xml)

    if 'plotfile' in extra_options:
      if not fit_info.sim_fit and (fit_info.is_tmin_vary or fit_info.is_tmax_vary):
          xml.append(plot_tag)
      elif not fit_info.sim_fit:
          fit_tag.append(plot_tag)
            
    self._addTask(xml)

  def doCorrMatrixRotation(self, pivot_info, rotate_mode, correlator_matrix, resulting_operator, 
                           mintime, maxtime, **extra_options):
    """Adds a 'DoCorrMatrixRotations' task to the SigmondInput object

    Args:
      pivot_info (PivotInfo): The pivot info to be used in the rotation
      rotate_mode (RotateMode): Specifies whether to rotate by bins or samplings
      correlator_matrix (sigmonbind.CorrelatorMatrix): The correlator
          matrix to be rotated.
      resulting_operator (sigmondbind.GIOperator): The GIOperator stub
          to construct the resulting rotated correlator matrix.
      mintime (int): the minimum timeslice to rotate.
      maxtime (int): the maximum timeslice to rotate.
      **neg_eig_alarm (float): the negative eigenvalue alarm value.
          For details see the run-sigmond and SigMonD documentation.
      **check_metric_errors (bool): Specifies whether estimates for
          the correlator matrix at the metric time should be output.
      **check_common_nullspace (bool): See run-sigmond and SigMonD
          documentation for details
      **pivot_name (str): specifies a name for the pivot, so that it
          can be referenced later.
      **improved_ops (list): a list of 'ImprovedOperator' objects.
          If some of the operators in correlator_matrix are improved
          operators, then these need to specified here.
      **show_transformation (bool): Specifies whether the transformation
          matrix computed to perfrom the rotation is printed to the log.
      **pivot_filename (str): specifies a file to write the pivot to.
      **pivot_overwrite (bool): if a pivot_filename is given, this
          specifies whether overwrite mode should be used.
      **set_imaginary_parts_zero (bool): if present, it will set the
          imaginary parts of the correlators to zero
      **set_to_zero (list): A list of correlators to set to zero
      **rotated_corrs_filename (str): specifies a file to write the
          rotated correlators to.
      **file_mode (sigmondbind.WriteMode): Specifies the write mode
          to use
      **corr_plotstub (str): a plotstub to use for the correlator plots.
          This is required if plotting of the rotated correlators is
          desired.
      **energy_plotstub (str): a plotstub to use for the effective
          energy plots. This is required if plotting of the rotated
          effective energies is desired.
      **sampling_mode (sigmond_xml.SamplingMode): specifies a sampling
          mode to use for the plots of the rotated correlators and
          effective energies.
      **eff_energy_type (EffEnergyType): Specifies the type of effective
          energy to compute for the plots.
      **timestep (int): Specifies the time step to be used to compute the
          effective energy for the plots.
      **symbol_color (SymbolColor): specifies the color to use for the data points.
      **symbol_type (SymbolType): specifies the shape of the data points to use.
      **rescale (float): specifies a rescaling factors for the correlator
          plots.
      **remove_off_diag (bool): specifies whether off diagonal elements
          should be removed
      **diagonly_time (int): specifies time when only diagonal rotated
          correlators are calculated
    """

    xml = ET.Element("Task")
    ET.SubElement(xml, "Action").text = "DoCorrMatrixRotation"
    ET.SubElement(xml, "MinTimeSep").text = str(mintime)
    ET.SubElement(xml, "MaxTimeSep").text = str(maxtime)
    ET.SubElement(xml, "RotateMode").text = rotate_mode.name
    ET.SubElement(xml, "Type").text = pivot_info.pivot_type.name

    if 'remove_off_diag' in extra_options and extra_options['remove_off_diag']:
      ET.SubElement(xml, "RemoveOffDiag")
    if 'diagonly_time' in extra_options:
      ET.SubElement(xml, "DiagonalOnlyTime").text = str(extra_options['diagonly_time'])

    pivot_init_tag = ET.SubElement(xml, "{}Initiate".format(pivot_info.pivot_type.name))
    ET.SubElement(pivot_init_tag, "RotatedCorrelator").append(resulting_operator.xml())
    if 'pivot_name' in extra_options:
      ET.SubElement(pivot_init_tag, "AssignName").text = extra_options['pivot_name']
    #pivot_init_tag.append(correlator_matrix.xml())
    
    if 'improved_ops' in extra_options:
      if extra_options['improved_ops']:
        #make list of replacements
        ops_to_replace = list()
        replacements = list()
        for improved_op_set in extra_options['improved_ops']:
          for improved_op in improved_op_set:
            replacements.append(improved_op[0])
          for i in range(1,len(improved_op_set[0]),2):
            ops_to_replace.append(improved_op_set[0][i])
          if len(ops_to_replace)!=len(replacements):
            logging.warning("Mismatch in improved operators")
            ops_to_replace = list()
            replacements = list()
      
        #make replacements in correlator matrix
        improved_correlator_matrix = correlator_matrix.xml()
        for child in improved_correlator_matrix.findall('GIOperatorString'):
          for i, (op) in enumerate(ops_to_replace):
            if child.text==op:
              child.text = replacements[i]
        pivot_init_tag.append(improved_correlator_matrix)
      else:
        pivot_init_tag.append(correlator_matrix.xml())
    else:
      pivot_init_tag.append(correlator_matrix.xml())

    if 'improved_ops' in extra_options:
      if extra_options['improved_ops']:

        #add improved operators to pivot
        improved_ops_tag = ET.SubElement(pivot_init_tag, "ImprovedOperators")
        for improved_op_set in extra_options['improved_ops']:
          for improved_op in improved_op_set:
            improved_op_tag = ET.SubElement(improved_ops_tag,"ImprovedOperator")
            opname_tag = ET.SubElement(improved_op_tag, "OpName")
            ET.SubElement(opname_tag, "GIOperatorString").text = improved_op[0]
            for i in range(1,len(improved_op),2):
              opterm_tag = ET.SubElement(improved_op_tag, "OpTerm")
              ET.SubElement(opterm_tag, "GIOperatorString").text = improved_op[i]
              ET.SubElement(opterm_tag, "Coefficient").text = improved_op[i+1]


    if pivot_info.pivot_type is sigmond_info.PivotType.SinglePivot:
      ET.SubElement(pivot_init_tag, "NormTime").text = str(pivot_info.norm_time)
      ET.SubElement(pivot_init_tag, "MetricTime").text = str(pivot_info.metric_time)
      ET.SubElement(pivot_init_tag, "DiagonalizeTime").text = str(pivot_info.diagonalize_time)
      ET.SubElement(pivot_init_tag, "MinimumInverseConditionNumber").text = str(1./pivot_info.max_condition_number)
    elif pivot_info.pivot_type is sigmond_info.PivotType.RollingPivot:
      ET.SubElement(pivot_init_tag, "NormTime").text = str(pivot_info.norm_time)
      ET.SubElement(pivot_init_tag, "MetricTime").text = str(pivot_info.metric_time)
      ET.SubElement(pivot_init_tag, "ZMagSqTime").text = str(pivot_info.diagonalize_time)
      ET.SubElement(pivot_init_tag, "MinimumInverseConditionNumber").text = str(1./pivot_info.max_condition_number)
    else:
      logging.error(f"pivot of type {pivot_info.pivot_type} currently not supported")

    if 'neg_eig_alarm' in extra_options:
      ET.SubElement(pivot_init_tag, "NegativeEigenvalueAlarm").text = str(extra_options['neg_eig_alarm'])
    if extra_options.get('check_metric_errors'):
      ET.SubElement(pivot_init_tag, "CheckMetricErrors")
    if extra_options.get('check_common_nullspace'):
      ET.SubElement(pivot_init_tag, "CheckCommonMetricMatrixNullSpace")
    if 'pivot_filename' in extra_options:
      write_pivot_tag = ET.SubElement(pivot_init_tag, "WritePivotToFile")
      ET.SubElement(write_pivot_tag, "PivotFileName").text = extra_options['pivot_filename']
      if extra_options.get('pivot_overwrite'):
        ET.SubElement(write_pivot_tag, "Overwrite")
      else:
        ET.SubElement(write_pivot_tag, "Update")
    if extra_options.get('show_transformation'):
      ET.SubElement(pivot_init_tag, "PrintTransformationMatrix")
    if extra_options.get('set_imaginary_parts_zero'):
      ET.SubElement(pivot_init_tag, "SetImaginaryPartsZero")
    if 'set_to_zero' in extra_options:
      for corr in extra_options['set_to_zero']:
        set_zero_xml = ET.SubElement(pivot_init_tag, "SetToZero")
        set_zero_xml.append(corr.xml())

    if 'rotated_corrs_filename' in extra_options:
      write_rotated_corr_tag = ET.SubElement(xml, "WriteRotatedCorrToFile")
      ET.SubElement(write_rotated_corr_tag, "RotatedCorrFileName").text = \
          extra_options['rotated_corrs_filename']
      if 'file_mode' in extra_options:
        ET.SubElement(write_rotated_corr_tag, "WriteMode").text = str(extra_options['file_mode'])

    if 'corr_plotstub' in extra_options:
      plot_tag = ET.SubElement(xml, "PlotRotatedCorrelators")
      ET.SubElement(plot_tag, "PlotFileStub").text = extra_options['corr_plotstub']
      if 'sampling_mode' in extra_options:
        ET.SubElement(plot_tag, "SamplingMode").text = str(extra_options['sampling_mode'])
      if 'symbol_color' in extra_options:
        ET.SubElement(plot_tag, "SymbolColor").text = extra_options['symbol_color'].value
      if 'symbol_type' in extra_options:
        ET.SubElement(plot_tag, "SymbolType").text = extra_options['symbol_type'].value
      if 'rescale' in extra_options:
        ET.SubElement(plot_tag, "Rescale").text = str(extra_options['rescale'])

    if 'energy_plotstub' in extra_options:
      plot_tag = ET.SubElement(xml, "PlotRotatedEffectiveEnergies")
      ET.SubElement(plot_tag, "PlotFileStub").text = extra_options['energy_plotstub']
      if 'eff_energy_type' in extra_options:
        ET.SubElement(plot_tag, "EffEnergyType").text = extra_options['eff_energy_type'].name
      if 'timestep' in extra_options:
        ET.SubElement(plot_tag, "TimeStep").text = str(extra_options['timestep'])
      if 'sampling_mode' in extra_options:
        ET.SubElement(plot_tag, "SamplingMode").text = str(extra_options['sampling_mode'])
      if 'symbol_color' in extra_options:
        ET.SubElement(plot_tag, "SymbolColor").text = extra_options['symbol_color'].value
      if 'symbol_type' in extra_options:
        ET.SubElement(plot_tag, "SymbolType").text = extra_options['symbol_type'].value

    self._addTask(xml)

  def doRotCorrMatInsertFitInfos(self, pivot_type, energies, amplitudes,
                                 **extra_options):
    """Adds a 'DoRotCorrMatInsertFitInfos' task to the SigmondInput object

    Args:
      pivot_type (PivotType): The pivot type to be used in the rotation
      energies (list of sigmondbind.MCObservable): energies to insert
      amplitudes (list of sigmondbind.MCObservable): amplitudes to insert
      **pivot_name (str): A name for the pivot so that it can be used in
          other tasks.
      **pivot_filename (str): the filename with the stored pivot
      **reorder (bool): specifies whether the energies should be reordered
    """

    xml = ET.Element("Task")
    ET.SubElement(xml, "Action").text = "DoRotCorrMatInsertFitInfos"
    ET.SubElement(xml, "Type").text = pivot_type.name
    pivot_init_tag = ET.SubElement(xml, f"{pivot_type.name}Initiate")
    if 'pivot_filename' in extra_options:
      read_tag = ET.SubElement(pivot_init_tag, "ReadPivotFromFile")
      ET.SubElement(read_tag, "PivotFileName").text = extra_options['pivot_filename']
      if 'pivot_name' in extra_options:
        ET.SubElement(pivot_init_tag, "AssignName").text = extra_options['pivot_name']

    elif 'pivot_name' in extra_options:
      memory_tag = ET.SubElement(pivot_init_tag, "GetFromMemory")
      ET.SubElement(memory_tag, "IDName").text = extra_options['pivot_name']
    else:
      logging.warning("No pivot given, task can not be added")
      return

    if 'reorder' in extra_options and extra_options['reorder']:
      ET.SubElement(xml, "ReorderByFitEnergy")

    level = 0
    for energy in energies:
      energy_fit_tag = ET.SubElement(xml, "EnergyFit")
      ET.SubElement(energy_fit_tag, "Level").text = str(level)
      ET.SubElement(energy_fit_tag, "Name").text = energy.getObsName()
      ET.SubElement(energy_fit_tag, "IDIndex").text = str(energy.getObsIndex())
      level += 1

    level = 0
    for amplitude in amplitudes:
      amplitude_fit_tag = ET.SubElement(xml, "RotatedAmplitude")
      ET.SubElement(amplitude_fit_tag, "Level").text = str(level)
      ET.SubElement(amplitude_fit_tag, "Name").text = amplitude.getObsName()
      ET.SubElement(amplitude_fit_tag, "IDIndex").text = str(amplitude.getObsIndex())
      level += 1

    self._addTask(xml)

  def doCorrMatrixZMagSquares(self, pivot_type, **extra_options):
    """Adds a 'DoCorrMatrixZMagSquares' task to the SigmondInput object

    Args:
      pivot_type (PivotType): The pivot type to be used in the rotation
      **pivot_name (str): A name for the pivot so that it can be used in
          other tasks.
      **pivot_filename (str): the filename with the stored pivot
      **plot_stub (str): a file stub to use for the zfactor plots
      **operators (list): a list of sigmondbind.OperatorInfo specifying
          the operators to calculate z-factors for
      **bar_color (str): a color to use for the bar plots
      **obsname (str): a name for the plots
    """
    xml = ET.Element("Task")
    ET.SubElement(xml, "Action").text = "DoCorrMatrixZMagSquares"
    ET.SubElement(xml, "Type").text = pivot_type.name
    pivot_init_tag = ET.SubElement(xml, f"{pivot_type.name}Initiate")
    if 'pivot_filename' in extra_options:
      read_tag = ET.SubElement(pivot_init_tag, "ReadPivotFromFile")
      ET.SubElement(read_tag, "PivotFileName").text = extra_options['pivot_filename']
      if 'pivot_name' in extra_options:
        ET.SubElement(pivot_init_tag, "AssignName").text = extra_options['pivot_name']

    elif 'pivot_name' in extra_options:
      memory_tag = ET.SubElement(pivot_init_tag, "GetFromMemory")
      ET.SubElement(memory_tag, "IDName").text = extra_options['pivot_name']
    else:
      logging.warning("No pivot given, task can not be added")
      return

    if 'plot_stub' in extra_options:
      plot_tag = ET.SubElement(xml, "DoPlots")
      ET.SubElement(plot_tag, "PlotFileStub").text = extra_options['plot_stub']
      if 'bar_color' in extra_options:
        ET.SubElement(plot_tag, "BarColor").text = extra_options['bar_color'].value

      for operator in extra_options.get('operators', []):
        zmag_tag = ET.SubElement(plot_tag, "ZMagSqPlot")
        zmag_tag.append(operator.xml())
        if 'obsname' in extra_options:
          ET.SubElement(zmag_tag, "ObsName").text = extra_options['obsname']
        ET.SubElement(zmag_tag, "FileSuffix").text = \
            util.str_to_file(operator.op_str())

    self._addTask(xml)


  def doGEVPCheck(self, pivot_type, **extra_options):
    """Adds a 'doGEVPCheck' task to the SigmondInput object

    Args:
      pivot_type (PivotType): The pivot type to be used in the rotation
      **pivot_name (str): A name for the pivot so that it can be used in
          other tasks.
      **pivot_filename (str): the filename with the stored pivot
      **plot_stub (str): a file stub to use for the zfactor plots
      **operators (list): a list of sigmondbind.OperatorInfo specifying
          the operators to calculate z-factors for
      **bar_color (str): a color to use for the bar plots
      **obsname (str): a name for the plots
    """
    xml = ET.Element("Task")
    ET.SubElement(xml, "Action").text = "DoGEVPCheck"
    ET.SubElement(xml, "Type").text = pivot_type.name
    ET.SubElement(xml, "MinTimeSep").text = str(extra_options.pop("tmin",2))
    ET.SubElement(xml, "MaxTimeSep").text = str(extra_options.pop("tmin",25))
    # ET.SubElement(xml, "NormTime").text = "tN"

    if 'minimizer_info' in extra_options:
      xml.append(extra_options['minimizer_info'].xml())
    pivot_init_tag = ET.SubElement(xml, f"{pivot_type.name}Initiate")
    if 'pivot_filename' in extra_options:
      read_tag = ET.SubElement(pivot_init_tag, "ReadPivotFromFile")
      ET.SubElement(read_tag, "PivotFileName").text = extra_options['pivot_filename']
      if 'pivot_name' in extra_options:
        ET.SubElement(pivot_init_tag, "AssignName").text = extra_options['pivot_name']

    elif 'pivot_name' in extra_options:
      memory_tag = ET.SubElement(pivot_init_tag, "GetFromMemory")
      ET.SubElement(memory_tag, "IDName").text = extra_options['pivot_name']
    else:
      logging.warning("No pivot given, task can not be added")
      return

    if 'plot_stub' in extra_options:
      # plot_tag = ET.SubElement(xml, "DoPlots")
      ET.SubElement(xml, "PlotFileStub").text = extra_options['plot_stub']
      # if 'bar_color' in extra_options:
      #   ET.SubElement(plot_tag, "BarColor").text = extra_options['bar_color'].value

      # for operator in extra_options.get('operators', []):
      #   zmag_tag = ET.SubElement(plot_tag, "ZMagSqPlot")
      #   zmag_tag.append(operator.xml())
      #   if 'obsname' in extra_options:
      #     ET.SubElement(zmag_tag, "ObsName").text = extra_options['obsname']
      #   ET.SubElement(zmag_tag, "FileSuffix").text = \
      #       util.str_to_file(operator.op_str())

    self._addTask(xml)

  def getFromPivot(self, pivot_type, energy_name, amplitude_name, **extra_options):
    """Adds a 'DoCorrMatrixZMagSquares' task to the SigmondInput object

    Args:
      pivot_type (PivotType): The pivot type to be used in the rotation
      energy_name (str): the obs_name for the reordered energies
      amplitude_name (str): the obs_name for the reordered amplitudes
      **pivot_name (str): A name for the pivot so that it can be used in
          other tasks.
      **pivot_filename (str): the filename with the stored pivot
    """
    xml = ET.Element("Task")
    ET.SubElement(xml, "Action").text = "GetFromPivot"
    ET.SubElement(xml, "Type").text = pivot_type.name
    pivot_init_tag = ET.SubElement(xml, f"{pivot_type.name}Initiate")
    if 'pivot_filename' in extra_options:
      read_tag = ET.SubElement(pivot_init_tag, "ReadPivotFromFile")
      ET.SubElement(read_tag, "PivotFileName").text = extra_options['pivot_filename']
      if 'pivot_name' in extra_options:
        ET.SubElement(pivot_init_tag, "AssignName").text = extra_options['pivot_name']

    elif 'pivot_name' in extra_options:
      memory_tag = ET.SubElement(pivot_init_tag, "GetFromMemory")
      ET.SubElement(memory_tag, "IDName").text = extra_options['pivot_name']
    else:
      logging.warning("No pivot given, task can not be added")
      return
    
    ET.SubElement(xml, "EnergyName").text = energy_name
    ET.SubElement(xml, "AmplitudeName").text = amplitude_name

    self._addTask(xml)

  def doReconstructEnergy(self, result_obs, energy_diff_obs, spatial_extent,
                          scattering_particle_energies, **extra_options):

    xml = ET.Element("Task")
    ET.SubElement(xml, "Action").text = "DoObsFunction"
    ET.SubElement(xml, "Type").text = "ReconstructEnergy"
    ET.SubElement(xml, "SpatialExtentNumSites").text = str(spatial_extent)
    if 'anisotropy' in extra_options:
      aniso_xml = ET.SubElement(xml, "Anisotropy")
      ET.SubElement(aniso_xml, "Name").text = extra_options['anisotropy'].getObsName()
      ET.SubElement(aniso_xml, "IDIndex").text = str(extra_options['anisotropy'].getObsIndex())

    result_xml = ET.SubElement(xml, "Result")
    ET.SubElement(result_xml, "Name").text = result_obs.getObsName()
    ET.SubElement(result_xml, "IDIndex").text = str(result_obs.getObsIndex())

    energy_diff_fit_xml = ET.SubElement(xml, "EnergyDifferenceFit")
    ET.SubElement(energy_diff_fit_xml, "Name").text = energy_diff_obs.getObsName()
    ET.SubElement(energy_diff_fit_xml, "IDIndex").text = str(energy_diff_obs.getObsIndex())

    for scat_energy_obs, psq in scattering_particle_energies:
      scat_energy_xml = ET.SubElement(xml, "ScatteringParticleEnergyFit")
      ET.SubElement(scat_energy_xml, "IntMomSquared").text = str(psq)
      ET.SubElement(scat_energy_xml, "Name").text = scat_energy_obs.getObsName()
      ET.SubElement(scat_energy_xml, "IDIndex").text = str(scat_energy_obs.getObsIndex())

    if 'sampling_mode' in extra_options:
      ET.SubElement(xml, "SamplingMode").text = str(extra_options['sampling_mode'])

    self._addTask(xml)

  def doReconstructAmplitude(self, result_obs, energy_diff_amp_obs,
                             scattering_particle_amps, **extra_options):

    xml = ET.Element("Task")
    ET.SubElement(xml, "Action").text = "DoObsFunction"
    ET.SubElement(xml, "Type").text = "ReconstructAmplitude"

    result_xml = ET.SubElement(xml, "Result")
    ET.SubElement(result_xml, "Name").text = result_obs.getObsName()
    ET.SubElement(result_xml, "IDIndex").text = str(result_obs.getObsIndex())

    energy_diff_amp_fit_xml = ET.SubElement(xml, "EnergyDifferenceAmplitudeFit")
    ET.SubElement(energy_diff_amp_fit_xml, "Name").text = energy_diff_amp_obs.getObsName()
    ET.SubElement(energy_diff_amp_fit_xml, "IDIndex").text = str(energy_diff_amp_obs.getObsIndex())

    for scat_amp_obs in scattering_particle_amps:
      scat_amp_xml = ET.SubElement(xml, "ScatteringParticleAmplitudeFit")
      ET.SubElement(scat_amp_xml, "Name").text = scat_amp_obs.getObsName()
      ET.SubElement(scat_amp_xml, "IDIndex").text = str(scat_amp_obs.getObsIndex())

    if 'sampling_mode' in extra_options:
      ET.SubElement(xml, "SamplingMode").text = str(extra_options['sampling_mode'])

    self._addTask(xml)

  def doEnergyDifference(self, result_obs, energy_obs, spatial_extent,
                         scattering_particle_energies, **extra_options):

    xml = ET.Element("Task")
    ET.SubElement(xml, "Action").text = "DoObsFunction"
    ET.SubElement(xml, "Type").text = "EnergyDifference"
    ET.SubElement(xml, "SpatialExtentNumSites").text = str(spatial_extent)
    if 'anisotropy' in extra_options:
      aniso_xml = ET.SubElement(xml, "Anisotropy")
      ET.SubElement(aniso_xml, "Name").text = extra_options['anisotropy'].getObsName()
      ET.SubElement(aniso_xml, "IDIndex").text = str(extra_options['anisotropy'].getObsIndex())

    result_xml = ET.SubElement(xml, "Result")
    ET.SubElement(result_xml, "Name").text = result_obs.getObsName()
    ET.SubElement(result_xml, "IDIndex").text = str(result_obs.getObsIndex())

    energy_diff_fit_xml = ET.SubElement(xml, "EnergyFit")
    ET.SubElement(energy_diff_fit_xml, "Name").text = energy_obs.getObsName()
    ET.SubElement(energy_diff_fit_xml, "IDIndex").text = str(energy_obs.getObsIndex())

    for scat_energy_obs, psq in scattering_particle_energies:
      scat_energy_xml = ET.SubElement(xml, "ScatteringParticleEnergyFit")
      ET.SubElement(scat_energy_xml, "IntMomSquared").text = str(psq)
      ET.SubElement(scat_energy_xml, "Name").text = scat_energy_obs.getObsName()
      ET.SubElement(scat_energy_xml, "IDIndex").text = str(scat_energy_obs.getObsIndex())

    if 'sampling_mode' in extra_options:
      ET.SubElement(xml, "SamplingMode").text = str(extra_options['sampling_mode'])

    self._addTask(xml)

  def doBoost(self, result, mom_squared, spatial_extent, frame_energy, **extra_options):
    '''
      anisotropy=None,
              boost_to_cm=True, sampling_mode=None, reference_energy=None):
    '''

    xml = ET.Element("Task")
    ET.SubElement(xml, "Action").text = "DoObsFunction"
    ET.SubElement(xml, "Type").text = "BoostEnergy"
    if 'boost_to_cm' not in extra_options or extra_options['boost_to_cm']:
      ET.SubElement(xml, "BoostToCM")

    result_xml = ET.SubElement(xml, "Result")
    ET.SubElement(result_xml, "Name").text = result.getObsName()
    ET.SubElement(result_xml, "IDIndex").text = str(result.getObsIndex())
    ET.SubElement(xml, "IntMomSquared").text = str(mom_squared)
    ET.SubElement(xml, "SpatialExtentNumSites").text = str(spatial_extent)
    frame_energy_xml = ET.SubElement(xml, "FrameEnergy")
    frame_energy_xml.append(frame_energy.xml())
    if 'anisotropy' in extra_options:
      aniso_xml = ET.SubElement(xml, "Anisotropy")
      aniso_xml.append(extra_options['anisotropy'].xml())

    if 'sampling_mode' in extra_options:
      ET.SubElement(xml, "SamplingMode").text = str(extra_options['sampling_mode'])

    self._addTask(xml)

  def doCopy(self, result_obs, in_obs, **extra_options):

    xml = ET.Element("Task")
    ET.SubElement(xml, "Action").text = "DoObsFunction"
    ET.SubElement(xml, "Type").text = "Copy"

    from_xml = ET.SubElement(xml, "InObservable")
    from_xml.append(in_obs.xml())

    to_xml = ET.SubElement(xml, "Result")
    ET.SubElement(to_xml, "Name").text = result_obs.getObsName()
    ET.SubElement(to_xml, "IDIndex").text = str(result_obs.getObsIndex())

    if 'sampling_mode' in extra_options:
      ET.SubElement(xml, "SamplingMode").text = str(extra_options['sampling_mode'])

    self._addTask(xml)

  def doRatio(self, result, numerator, denominator, **extra_options):

    xml = ET.Element("Task")
    ET.SubElement(xml, "Action").text = "DoObsFunction"
    ET.SubElement(xml, "Type").text = "Ratio"
    result_xml = ET.SubElement(xml, "Result")
    ET.SubElement(result_xml, "Name").text = result.getObsName()
    ET.SubElement(result_xml, "IDIndex").text = str(result.getObsIndex())
    numerator_xml = ET.SubElement(xml, "Numerator")
    numerator_xml.append(numerator.xml())
    denominator_xml = ET.SubElement(xml, "Denominator")
    denominator_xml.append(denominator.xml())

    if 'sampling_mode' in extra_options:
      ET.SubElement(xml, "SamplingMode").text = str(extra_options['sampling_mode'])

    self._addTask(xml)

  def doReconstructRatioFromMultiExpFit(self, result1, result2, numerator, denominator1, denominator2, **extra_options):

    xml = ET.Element("Task")
    ET.SubElement(xml, "Action").text = "DoFit"
    ET.SubElement(xml, "Type").text = "ReconstructRatioFromMultiExpFit"
    if 'minimizer_info' in extra_options:
      xml.append(extra_options['minimizer_info'].xml())
    
    numerator_xml = ET.SubElement(xml, "Numerator")
    numerator_xml.append(numerator.xml())
    denominator1_xml = ET.SubElement(xml, "Denominator1")
    denominator1_xml.append(denominator1.xml())
    denominator2_xml = ET.SubElement(xml, "Denominator2")
    denominator2_xml.append(denominator2.xml())
    ratio_xml = ET.SubElement(xml, "Ratio")
    ratio_xml.append(result1.xml())
    fit_ratio_xml = ET.SubElement(xml, "FitRatio")
    fit_ratio_xml.append(result2.xml())

    if 'plotfile' in extra_options:
        plot_tag = ET.SubElement(xml, "DoEffectiveEnergyPlot")
        ET.SubElement(plot_tag, "PlotFile").text = extra_options['plotfile']
        if 'plot_info' in extra_options:
          plot_info = extra_options['plot_info']
          ET.SubElement(plot_tag, "CorrName").text = plot_info.corrname
          ET.SubElement(plot_tag, "TimeStep").text = str(plot_info.timestep)
          ET.SubElement(plot_tag, "SymbolColor").text = plot_info.symbol_color.value
          ET.SubElement(plot_tag, "SymbolType").text = plot_info.symbol_type.value
          ET.SubElement(plot_tag, "Goodness").text = plot_info.goodness.value
          ET.SubElement(plot_tag, "ShowApproach")
          if plot_info.max_relative_error > 0.:
            ET.SubElement(plot_tag, "MaxRelativeErrorToPlot").text = str(plot_info.max_relative_error)

    self._addTask(xml)
#     print(ET.tostring(xml))


  def doLinearSuperposition(self, result, summands, **extra_options):

    xml = ET.Element("Task")
    ET.SubElement(xml, "Action").text = "DoObsFunction"
    ET.SubElement(xml, "Type").text = "LinearSuperposition"
    result_xml = ET.SubElement(xml, "Result")
    ET.SubElement(result_xml, "Name").text = result.getObsName()
    ET.SubElement(result_xml, "IDIndex").text = str(result.getObsIndex())

    for summand in summands:
      summand_xml = ET.SubElement(xml, "Summand")
      summand_xml.append(summand['observable'].xml())
      ET.SubElement(summand_xml, "Coefficient").text = str(summand['coefficient'])

    if 'sampling_mode' in extra_options:
      ET.SubElement(xml, "SamplingMode").text = str(extra_options['sampling_mode'])

    self._addTask(xml)

  def doExp(self, result_obs, in_obs, **extra_options):

    xml = ET.Element("Task")
    ET.SubElement(xml, "Action").text = "DoObsFunction"
    ET.SubElement(xml, "Type").text = "Exp"

    from_xml = ET.SubElement(xml, "InObservable")
    from_xml.append(in_obs.xml())

    to_xml = ET.SubElement(xml, "Result")
    ET.SubElement(to_xml, "Name").text = result_obs.getObsName()
    ET.SubElement(to_xml, "IDIndex").text = str(result_obs.getObsIndex())

    if 'sampling_mode' in extra_options:
      ET.SubElement(xml, "SamplingMode").text = str(extra_options['sampling_mode'])

    self._addTask(xml)

  def doLog(self, result_obs, in_obs, **extra_options):

    xml = ET.Element("Task")
    ET.SubElement(xml, "Action").text = "DoObsFunction"
    ET.SubElement(xml, "Type").text = "Log"

    from_xml = ET.SubElement(xml, "InObservable")
    from_xml.append(in_obs.xml())

    to_xml = ET.SubElement(xml, "Result")
    ET.SubElement(to_xml, "Name").text = result_obs.getObsName()
    ET.SubElement(to_xml, "IDIndex").text = str(result_obs.getObsIndex())

    if 'sampling_mode' in extra_options:
      ET.SubElement(xml, "SamplingMode").text = str(extra_options['sampling_mode'])

    self._addTask(xml)

  def doCorrelatorMatrixSuperposition(self, result_operators, original_operators,
                                      mintime, maxtime, **extra_options):
    """Adds a 'CorrelatorMatrixSuperposition' task to the SigmondInput
       object

    Args:
      result_operators (list): A list of objects of type 'GIOperator'
          that are used to form the resulting correlator matrix.
      original_operators (list(dict())): Each list element
          corresponds to a correlator matrix to be averaged. Each
          key corresponds to a corresponding operator in result_operators,
          and the value is the actual operator to average.
      mintime (int): the mininum time slice considered for the
          correlators
      maxtime (int): the maximum time slice considered for the
          correlators
      **coefficients (dict): This is a dictionary which uses operators
          as keys. The value of each key corresponds to the coefficient
          to be used for this operator in the average. Any operators
          not present in the dictionary will assume a value of '1.0'
          for the coefficient.
      **hermitian (bool): specifies whether the correlator matrices
          should be assumed to be Hermitian.
      **filename (str): specifies the filename to write the resulting
          correlator matrix bins to.
      **file_type (DataFormat): Specifies whether to save by bins or
          samplings
      **file_mode (sigmondbind.WriteMode): Specifies the write mode
          to use
    """

    xml = ET.Element("Task")
    ET.SubElement(xml, "Action").text = "DoObsFunction"
    ET.SubElement(xml, "Type").text = "CorrelatorMatrixSuperposition"
    result_op_tag = ET.SubElement(xml, "ResultOperatorOrderedList")
    for result_operator in result_operators:
      result_op_tag.append(result_operator.xml())

    coefficients = extra_options.get('coefficients', dict())
    for operators in original_operators:
      operators_tag = ET.SubElement(xml, "OperatorOrderedList")
      for result_operator in result_operators:
        raw_operator = operators[result_operator]
        coeff = coefficients.get(raw_operator, 1.)
        item_tag = ET.SubElement(operators_tag, "Item")
        item_tag.append(raw_operator.xml())
        ET.SubElement(item_tag, "Coefficient").text = str(coeff)

    ET.SubElement(xml, "MinimumTimeSeparation").text = str(mintime)
    ET.SubElement(xml, "MaximumTimeSeparation").text = str(maxtime)
    if extra_options.get('hermitian'):
      ET.SubElement(xml, "HermitianMatrix")

    if 'filename' in extra_options:
      write_xml = ET.SubElement(xml, "WriteToFile")
      ET.SubElement(write_xml, "FileName").text = extra_options['filename']
      ET.SubElement(write_xml, "FileType").text = extra_options['file_type'].name
      if 'file_mode' in extra_options:
        ET.SubElement(write_xml, "WriteMode").text = str(extra_options['file_mode'])

    self._addTask(xml)

  def writeToFile(self, file_name, observables,
                  file_type=sigmond_info.DataFormat.samplings, 
                  file_mode=sigmond.WriteMode.Overwrite, hdf5 = False):
    """Adds a 'WriteToFile' task to the SigmondInput object

    Args:
      file_name (str): Specifies the file to write to
      observables (list): A list of observables of type MCObsInfo to write to the file
      **file_type (DataFormat): whether the data should be written as bins or samplings
      **file_mode (sigmondbind.WriteMode): Specifies the write mode
          to use
    """

    if not observables:
      logging.warning("No observables passed to writeToFile")

    xml = ET.Element("Task")
    ET.SubElement(xml, "Action").text = "WriteToFile"
    ET.SubElement(xml, "FileName").text = file_name
    ET.SubElement(xml, "FileType").text = file_type.name
    ET.SubElement(xml, "WriteMode").text = str(file_mode)
    if hdf5:
        ET.SubElement(xml, "FileFormat").text = "hdf5" 
        
    for observable in observables:
      xml.append(observable.xml())

    self._addTask(xml)

  def writeCorrMatToFile(self, file_name, corrs,
                  file_type=sigmond_info.DataFormat.samplings, 
                  file_mode=sigmond.WriteMode.Overwrite, tmin=None, tmax=None):
    """Adds a 'WriteToFile' task to the SigmondInput object

    Args:
      file_name (str): Specifies the file to write to
      corrs (sigmond.CorrelatorMatrixInfo):The correlation Matrix info
      **file_type (DataFormat): whether the data should be written as samplings or samplings
      **file_mode (sigmondbind.WriteMode): Specifies the write mode
          to use
      **tmin (int): min t to write
      **tmax (int): max t to write
    """

    if not corrs:
      logging.warning("No correlators passed to writeCorrMatToFile")

    xml = ET.Element("Task")
    ET.SubElement(xml, "Action").text = "WriteCorrMatToFile"
    ET.SubElement(xml, "FileName").text = file_name
    ET.SubElement(xml, "FileType").text = "samplings" #file_type.name
    ET.SubElement(xml, "FileFormat").text = "hdf5"
    ET.SubElement(xml, "WriteMode").text = str(file_mode)
    xml.append(corrs.xml())
    
    if tmin:
        ET.SubElement(xml, "MinTimeSep").text = str(tmin)
    else:
        ET.SubElement(xml, "MinTimeSep").text = str(0)
    if tmax:
        ET.SubElement(xml, "MaxTimeSep").text = str(tmax)
    else:
        ET.SubElement(xml, "MaxTimeSep").text = str(64)
    self._addTask(xml)

