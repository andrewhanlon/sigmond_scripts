import os
import shutil
import logging
import xml.etree.ElementTree as ET
import h5py
import pylatex
import numpy as np
import regex

import tasks.task
import utils.plotting
import data_handling.data_files
import sigmond_info.sigmond_input
import sigmond_info.sigmond_info
import operator_info.operator_set
import utils.util as util

import sigmond

class Spectrum(tasks.task.Task):
  """ The Spectrum Task!

  TODO: 
    - put all the PSQ{psq}/{irrep} stuff into variables at the top
    - The non-interacting levels are calculated in _getNonInteractingLevels,
      and these rely on the FitInfos in self.spectrum having the
      non_interacting levels in them. Should this be done this way?
    - tmin plots for single hadrons
  """

  task_type = "spectrum"

  sh_name = "single_hadrons"
  ref_name = "ref"
  ordered_energy = "ordered_energy"
  ordered_amplitude = "ordered_amplitude"

  def initialize(self, spectrum, **extra_options):
    """ initializes the needed objects to perform fits to correlators
        with sigmond

    Args:
      spectrum ({OperatorSet: [FitInfos]): a dictionary
          mapping each operator basis to a list of fit infos
          for each level
      **flagged_levels ({OperatorSet: [bool]): a dictionary
          mapping each operator to its 'flag'
      **reorder (bool): whether energies should be reordered
      **tmin_infos ({OperatorSet: [[FitInfo]]): a dictionary
          containing the differnt tmin fits to do.
      **minimizer_info (sigmondbind.MinimizerInfo): Speicifies the 
          minimizer to use and the info to pass to the minimizer
      **fit_plot_info (FitPlot): Contains information for how to make
          the fit plots.
      **default_tmin_plot_info (TMinPlotInfo): Contains information on
          defaults for tmin plots
      **reference_fit_info (FitInfo): The fit to be used as a
          reference energy.
      **scattering_particles ({ScatteringParticle: FitInfo}): A dictionary
          of the scattering particles to their fit infos.
      **thresholds (list): a list of thresholds to draw on spectrum
      **latex_map (dict): a mapping of the latex to use
      **bar_color (SymbolColor): color for the zfactor plots
      **rotate_labels (bool): whether the irrep labels on the spectrum plot
          should be rotated
      **non_interacting_energy_labels (bool): whether it prints the
          non interacting levels on the spectrum tikz plot or not.
      **single_hadrons (list): list of the two single hadron names that 
          make up the correlators for this channel to calculate the q^2_cm.
          Not set up for anything other than two operator interactions.
    """

    self.spectrum = spectrum
    self.flagged_levels = extra_options.pop('flagged_levels', dict())
    self.reorder = extra_options.pop('reorder', True)
    self.tmin_infos = extra_options.pop('tmin_infos', dict())
    self.minimizer_info = extra_options.pop('minimizer_info', sigmond.MinimizerInfo())
    self.fit_plot_info = extra_options.pop('fit_plot_info', sigmond_info.sigmond_info.FitPlotInfo())
    self.default_tmin_plot_info = extra_options.pop('default_tmin_plot_info', sigmond_info.sigmond_info.TMinPlotInfo())
    self.scattering_particles = extra_options.pop('scattering_particles', dict())
    self.reference_fit_info = extra_options.pop('reference_fit_info', None)
    self.thresholds = extra_options.pop('thresholds', list())
    self.latex_map = extra_options.pop('latex_map', dict())
    self.bar_color = extra_options.pop('bar_color', sigmond_info.sigmond_info.SymbolColor.cyan)
    self.rotate_labels = extra_options.pop('rotate_labels', False)
    self.plot_width_factor = extra_options.pop('plot_width_factor', 1.)
    self.non_interacting_energy_labels = extra_options.pop('non_interacting_energy_labels', False)
    self.single_hadrons = extra_options.pop('single_hadrons',list())
    self.write_params_to_file = extra_options.pop('write_params_to_file',False)

    util.check_extra_keys(extra_options, 'Spectrum.initialize')
    
    self.ni_fit_info = {}

  def readConfig(self, **task_options):
    """readConfig function for Spectrum

    YAML:
      subtractvev: true   # optional
      noise_cutoff: 3.    # optional
      non_interacting_energy_labels: true #optional
      reorder: false   # optional
      write_params_to_file: false #optional
      
      # optional
      minimizer_info:
        minimizer: lmder
        parameter_rel_tol: 1e-6
        chisquare_rel_tol: 1e-4
        max_iterations: 1024
        verbosity: high

      # optional
      fit_plot_info:
        timestep: 3
        show_approach: true
        goodness: chisq
        corrname: standard
        symbol_color: blue
        symbol_type: circle
        max_relative_error: 0.0

      bar_color: cyan

      rotate_labels: False
      plot_width_factor: 1.

      # optional (no tmin plots will be made if missing)
      tmin_plot_info:
        obsname: standard
        symbol_type: circle
        goodfit_color: blue
        badfit_color: red
        goodfit_hollow: false
        badfit_hollow: false
        quality_threshold: 0.1

      # optional
      reference_fit_info:
        # option 1: if defined in "scattering_particles", uses that info instead
        scattering_particle: pi
        # option 2: uses this info for reference fit
        operator: "isotriplet S=0 P=(0,0,0) A1um P 0"
        model: 1-exp
        tmin: 14
        tmax: 40

      # optional (used for ratio fits)
      scattering_particles:  
        - name: pi
          fits:
            - operator: isotriplet S=0 P=(0,0,0) A1um P 0
              model: 1-exp
              tmin: 12
              tmax: 30
              tmin_info:
                - model: 1-exp
                  tmin_min: 3
                  tmin_max: 10
                  extra_tmaxes: [18, 19]
            - operator: isotriplet S=0 PSQ=1 A2m P 0
              model: 2-exp
              tmin: 5
              tmax: 15
            - operator: isotriplet S=0 PSQ=2 A2m P 0
              model: log-1-exp
              tmin: 8
              tmax: 12
            ...
        - name: lambda
          fits:
            - ...
            ...
        ...

      # optional (used for drawing thresholds on spectrum plot).
      # requires scattering particles
      thresholds:
        - [lambda, lambda]
        - [pi, pi]
        - [pi, pi, pi, pi]
        - ...
        ...

      # optional (used to translates scattering particle names and
      # irreps to latex
      latex_map:
        ref: \pi
        pi: \pi
        lambda: \Lambda
        A1: 'A_1'

      spectrum:
        - name: op_list1
          operators:
            - isodoublet S=0 P=(0,0,0) A1um P 0
            - isodoublet S=0 P=(0,0,0) A1um P 1
          non_interacting_levels:    # required for ratio fits
            - [pi(0), pi(2), pi(3)]
            - [pi(1), pi(2), pi(3)]
          levels:
            - model: 2-exp
              tmin: 10
              tmax: 20
            - model: 1-exp
              tmin: 14
              tmax: 20
          tmin_info: (optional)
            - fit_infos:
              - model: 1-exp
                tmin_min: 3
                tmin_max: 10
                extra_tmaxes: [18, 19]

        - name: rotated_op_basis
          pivot_info:
            pivot_type: single_pivot
            norm_time: 5
            metric_time: 5
            diagonalize_time: 10
            max_condition_number: 100
          non_interacting_levels:    # required for ratio fits
            - [pi(0), pi(2), pi(3)]
            - [pi(1), pi(2), pi(3)]
          levels:
            - model: 2-exp
              tmin: 12
              tmax: 25
              flag: True        # default is False
            - model: 1-exp
              tmin: 16
              tmax: 23
        - ...
        ...
    """

    subtractvev = task_options.pop('subtractvev', True)
    global_noise_cutoff = task_options.pop('noise_cutoff', 0.)
    task_options['minimizer_info'] = sigmond_info.sigmond_info.getMinimizerInfo(task_options)
    task_options['fit_plot_info'] = sigmond_info.sigmond_info.FitPlotInfo.createFromConfig(task_options)
    task_options['default_tmin_plot_info'] = sigmond_info.sigmond_info.TMinPlotInfo.createFromConfig(task_options)

    ref_fit_info = task_options.pop('reference_fit_info', None)
    if ref_fit_info is not None:
      if 'scattering_particle' not in ref_fit_info.keys():
        ref_fit_info['noise_cutoff'] = ref_fit_info.get('noise_cutoff', global_noise_cutoff)
        ref_fit_info = sigmond_info.fit_info.FitInfo.createFromConfig(ref_fit_info)
    task_options['reference_fit_info'] = ref_fit_info

    util.updateOption(task_options, 'bar_color', sigmond_info.sigmond_info.SymbolColor)

    # check for scattering_particles
    tmin_infos = dict()
    tmin_infos[self.sh_name] = dict()
    scattering_particles = dict()
    for scattering_particle in task_options.pop('scattering_particles', []):
      try:
        name = scattering_particle.pop('name')
        fits = scattering_particle.pop('fits')
        for fit in fits:
          tmin_info_conf = fit.pop('tmin_info', [])
          operator = operator_info.operator.Operator(fit.pop('operator'))
          use_irrep = fit.pop('use_irrep', False)
          fit_model = sigmond_info.fit_info.FitModel(fit.pop('model'))
          noise_cutoff = fit.pop('noise_cutoff', global_noise_cutoff)
          fit_info = sigmond_info.fit_info.FitInfo(operator, fit_model, **fit, noise_cutoff=noise_cutoff)
          psq = operator.psq
          if use_irrep:
            irrep = operator.getLGIrrep()
            scattering_particle = sigmond_info.sigmond_info.ScatteringParticle(name, psq, False, irrep)
          else:
            scattering_particle = sigmond_info.sigmond_info.ScatteringParticle(name, psq, False)

          if scattering_particle in scattering_particles:
            logging.error(f"Scattering particle '{scattering_particle}' encountered twice")

          scattering_particles[scattering_particle] = fit_info

          tmin_infos[self.sh_name][scattering_particle] = list()
          for tmin_fit_info in tmin_info_conf:
            fit_model = sigmond_info.fit_info.FitModel(tmin_fit_info['model'])
            tmin_min = fit_info.tmin
            tmin_max = -1
            tmax_min = -1
            tmax_max = fit_info.tmax
            if ('tmin_min' in tmin_fit_info.keys()) and ('tmin_max' in tmin_fit_info.keys()):
                tmin_min = tmin_fit_info['tmin_min']
                tmin_max = tmin_fit_info['tmin_max']
            elif ('tmax_min' in tmin_fit_info.keys()) and ('tmax_max' in tmin_fit_info.keys()):
                tmax_min = tmin_fit_info['tmax_min']
                tmax_max = tmin_fit_info['tmax_max']
            else:
                logging.error(f"Invalid TminVary or TmaxVary config")

            if 'extra_tmaxes' in tmin_fit_info:
              tmaxes = tmin_fit_info.get('extra_tmaxes')
              if isinstance(tmaxes, int):
                tmaxes = [int(tmaxes)]
              elif '-' in tmaxes:
                tmaxes = list(map(int, tmaxes.split('-')))
                tmaxes = list(range(tmaxes[0], tmaxes[-1]+1))
              elif isinstance(tmaxes, list):
                tmaxes = list(map(int, tmaxes))
              else:
                logging.error("Invalid 'extra_tmaxes' parameter")

              tmaxes.append(fit_info.tmax)

            else:
              tmaxes = [fit_info.tmax]

#             for tmin_tmax in sorted(tmaxes):
#               tmin_fit_info = sigmond_info.fit_info.FitInfo(
#                   operator, fit_model, tmin, tmin_tmax, subtractvev, False,
#                   fit_info.exclude_times, fit_info.noise_cutoff, None, tmin_max)
#               tmin_infos[self.sh_name][scattering_particle].append(tmin_fit_info)
            
            if tmax_min<0:
                for tmin_tmax in sorted(tmaxes):
                    tmin_fit_info = sigmond_info.fit_info.FitInfo(
                        operator, fit_model, tmin_min, tmin_tmax, subtractvev, False, fit_info.exclude_times, 
                        fit_info.noise_cutoff, None, tmin_max)

                    tmin_infos[self.sh_name][scattering_particle].append(tmin_fit_info)
            else:
                tmin_fit_info = sigmond_info.fit_info.FitInfo(
                        operator, fit_model, tmin_min, tmax_max, subtractvev, False, fit_info.exclude_times, 
                        fit_info.noise_cutoff, None, tmin_max, tmax_min)
                tmin_infos[self.sh_name][scattering_particle].append(tmin_fit_info)

      except KeyError as err:
        logging.error(f"Missing required key in 'scattering_particles': {err}")

    scattering_particle_operators = dict()
    for scattering_particle, fit_info in scattering_particles.items():
      scattering_particle_operators[scattering_particle] = fit_info.operator
            
    task_options['scattering_particles'] = scattering_particles
    if ref_fit_info:
      if type(task_options['reference_fit_info'])==dict:
        scattering_particle = sigmond_info.sigmond_info.ScatteringParticle(ref_fit_info['scattering_particle'], 0, False)
        task_options['reference_fit_info'] = task_options['scattering_particles'][scattering_particle]

    # get the spectrum now
    spectrum = dict()
    flagged_levels = dict()
    try:
      for energy_set in task_options.pop('spectrum', []):
        operator_set = operator_info.operator_set.getOperatorSet(energy_set)
        operators = operator_set.getRotatedOperators() if operator_set.is_rotated else operator_set.operators
        non_interacting_levels = energy_set.pop('non_interacting_levels', [None]*len(operators))
        levels = energy_set.pop('levels')
        tmin_info_confs = energy_set.pop('tmin_info', [None]*len(operators))
        ratio_fit_log = energy_set.pop('ratio_log', None )

        spectrum[operator_set] = list()
        flagged_levels[operator_set] = list()
        tmin_infos[operator_set] = list()
        
        if len(non_interacting_levels) != len(operators):
            non_interacting_levels += [None]*(len(operators)-len(non_interacting_levels))
            
        if ratio_fit_log:
          ratio_xml = ET.parse(ratio_fit_log)
          ratio_fit_xmls = {item.find('InteractingOperator/GIOperatorString').text:item for item in ratio_xml.findall(f"./Task/DoFit/[Type='TemporalCorrelatorInteractionRatio']")}
          sh_fit_xmls = {item.find('TemporalCorrelatorFit/GIOperatorString').text:item for item in ratio_xml.findall(f"./Task/DoFit/[Type='TemporalCorrelator']")}
        
        for operator, level, non_interacting_level, tmin_info in zip(operators, levels, non_interacting_levels, tmin_info_confs):
          fit_model = sigmond_info.fit_info.FitModel(level.pop('model'))
          tmin = int(level.pop('tmin'))
          tmax = int(level.pop('tmax'))
          exclude_times = level.pop('exclude_times', [])
          noise_cutoff = level.pop('noise_cutoff', global_noise_cutoff)
          ratio = level.pop('ratio', False)
          flag = level.pop('flag', False)
          max_level = level.pop('max_level',6)
          initial_gap = level.pop('initial_gap',1.0)
          repeating_gap = level.pop('repeating_gap',1.0)
          sim_fit = level.pop('sim_fit',False)
          initial_params = level.pop('initial_params',{})
          if ratio_fit_log:
            for key, val in [('Energy','FitParameter0'),('Amplitude','FitParameter1')]:
              initial_params[key] = float(ratio_fit_xmls[str(operator)].find(f'BestFitResult/{val}/MCEstimate/FullEstimate').text)
            for i, sh in enumerate(non_interacting_level):
                sh0 = sigmond_info.sigmond_info.ScatteringParticle.create(sh)
                this_xml = sh_fit_xmls[str(scattering_particles[sh0].operator)]
                for key, val in [(f'SH{i+1}Gap','FitParameter2'),(f'SH{i+1}GapAmp','FitParameter3')]:
                    initial_params[key] = float(this_xml.find(f'BestFitResult/{val}/MCEstimate/FullEstimate').text)
            initial_params["NumGap"] = 1.0
            initial_params["NumGapAmp"] = 0.2 #make input
            
          util.check_extra_keys(level, "spectrum.level")

          if non_interacting_level is None:
            non_interacting_operators = None
          else:
            non_interacting_operators = sigmond_info.sigmond_info.NonInteractingOperators.create(
                scattering_particle_operators, non_interacting_level)

          
          fit_info = sigmond_info.fit_info.FitInfo(
              operator, fit_model, tmin, tmax, subtractvev, ratio, exclude_times, noise_cutoff,
              non_interacting_operators, -1, -1, max_level,initial_gap,repeating_gap,sim_fit,initial_params)

          spectrum[operator_set].append(fit_info)
          flagged_levels[operator_set].append(flag)

          if tmin_info is not None:
            tmin_info_list = list()
            for tmin_fit_info in tmin_info['fit_infos']:
              fit_model = sigmond_info.fit_info.FitModel(tmin_fit_info['model'])
              ratio = tmin_fit_info.get('ratio', False)
              ratio = tmin_fit_info['ratio']
              sim_fit = tmin_fit_info.get('sim_fit', False)
#               if 'sim_fit' in tmin_fit_info:
#                 sim_fit = tmin_fit_info['sim_fit']
#               initial_params = level.pop('initial_params',{})
              tmin_min = tmin
              tmin_max = -1
              tmax_min = -1
              tmax_max = tmax
              if ('tmin_min' in tmin_fit_info.keys()) and ('tmin_max' in tmin_fit_info.keys()):
                  tmin_min = tmin_fit_info['tmin_min']
                  tmin_max = tmin_fit_info['tmin_max']
              elif ('tmax_min' in tmin_fit_info.keys()) and ('tmax_max' in tmin_fit_info.keys()):
                  tmax_min = tmin_fit_info['tmax_min']
                  tmax_max = tmin_fit_info['tmax_max']
              else:
                  logging.error(f"Invalid TminVary or TmaxVary config")
              if 'extra_tmaxes' in tmin_fit_info:
                tmaxes = tmin_fit_info.get('extra_tmaxes')
                if isinstance(tmaxes, int):
                  tmaxes = [int(tmaxes)]
                elif '-' in tmaxes:
                  tmaxes = list(map(int, tmaxes.split('-')))
                  tmaxes = list(range(tmaxes[0], tmaxes[-1]+1))
                elif isinstance(tmaxes, list):
                  tmaxes = list(map(int, tmaxes))
                else:
                  logging.error("Invalid 'extra_tmaxes' parameter")

                tmaxes.append(tmax)

              else:
                tmaxes = [tmax]

              if tmax_min<0:
                for tmin_tmax in sorted(tmaxes):
                    fit_info = sigmond_info.fit_info.FitInfo(
                        operator, fit_model, tmin_min, tmin_tmax, subtractvev, ratio, exclude_times, noise_cutoff,
                        non_interacting_operators, tmin_max, sim_fit=sim_fit,initial_params =initial_params)

                    if fit_info not in tmin_info_list:
                      tmin_info_list.append(fit_info)
              else:
                fit_info = sigmond_info.fit_info.FitInfo(
                        operator, fit_model, tmin_min, tmax_max, subtractvev, ratio, exclude_times, noise_cutoff,
                        non_interacting_operators, tmin_max, tmax_min, sim_fit=sim_fit,initial_params =initial_params)
                if fit_info not in tmin_info_list:
                  tmin_info_list.append(fit_info)

            tmin_infos[operator_set].append(tmin_info_list)

    except KeyError as err:
      logging.error(f"Missing required key in 'spectrum': {err}")

    task_options['tmin_infos'] = tmin_infos

    self.initialize(spectrum, flagged_levels=flagged_levels, **task_options)

  @property
  def results_dir(self):
    results_dir = os.path.join(super().results_dir, self.task_name)
    os.makedirs(results_dir, exist_ok=True)
    return results_dir

  @property
  def samplings_dir(self):
    samp_dir = os.path.join(self.results_dir, "samplings")
    os.makedirs(samp_dir, exist_ok=True)
    return samp_dir

  @property
  def hdf5_filename(self):
    filename = f"energy_samplings_{self.task_name}_rebin{self.rebin}.hdf5"
    return os.path.join(self.samplings_dir, filename)

  @property
  def qsqr_filename(self):
    filename = f"qsqr_samplings_{self.task_name}_rebin{self.rebin}.hdf5"
    return os.path.join(self.samplings_dir, filename)

  @property
  def estimates_filename(self):
    filename = f"energy_estimates_{self.task_name}_rebin{self.rebin}_{self.sampling_mode}.csv"
    return os.path.join(self.results_dir, filename)
  
  @property
  def samplings_filename(self):
    filename = f"energy_samplings_{self.task_name}_rebin{self.rebin}.smp"
    return os.path.join(self.samplings_dir, filename)

  def params_filename(self, tag="", with_root=False):
    if tag:
        filename = f"param_samplings_{tag}_{self.task_name}_rebin{self.rebin}.hdf5"
    else:
        filename = f"param_samplings_{self.task_name}_rebin{self.rebin}.hdf5"
    if with_root:
        filename+="[params]"
    return os.path.join(self.samplings_dir, filename)

  def zfactor_plotdir(self, operator_set_name):
    plotdir = os.path.join(self.plotdir, "zfactors", self.task_name, f"rebin{self.rebin}", operator_set_name)
    os.makedirs(plotdir, exist_ok=True)
    return plotdir

  def zfactor_plotstub(self, operator_set_name):
    plotdir = self.zfactor_plotdir(operator_set_name)
    plotstub = "zfactor"
    return os.path.join(plotdir, plotstub)

  def zfactor_plotfile(self, operator_set_name, operator, extension):
    plotstub = self.zfactor_plotstub(operator_set_name)
    plotfile = f"{plotstub}_{operator!r}.{extension.value}"
    return plotfile

  def fit_plotdir(self, operator_set_name):
    plotdir = os.path.join(self.plotdir, "fit_plots", self.task_name, f"rebin{self.rebin}", operator_set_name)
    os.makedirs(plotdir, exist_ok=True)
    return plotdir

  def fit_plotfile(self, operator_set_name, level, extension):
    plotdir = self.fit_plotdir(operator_set_name)
    plotfile = f"fit_{level}.{extension.value}"
    return os.path.join(plotdir, plotfile)

  def fit_sh_plotdir(self):
    plotdir = os.path.join(self.plotdir, "fit_plots", self.task_name, f"rebin{self.rebin}", self.sh_name)
    os.makedirs(plotdir, exist_ok=True)
    return plotdir

  def fit_sh_plotfile(self, name, extension):
    plotdir = self.fit_sh_plotdir()
    plotfile = f"fit_{name}.{extension.value}"
    return os.path.join(plotdir, plotfile)

  def tmin_sh_plotdir(self):
    plotdir = os.path.join(self.plotdir, "tmin_plots", self.task_name, f"rebin{self.rebin}", self.sh_name)
    os.makedirs(plotdir, exist_ok=True)
    return plotdir

  def tmin_sh_plotfile(self, fit_info, extension):
    plotdir = self.tmin_sh_plotdir()
    plotfile = f"tmin_{fit_info.plotfile}.{extension.value}"
    return os.path.join(plotdir, plotfile)

  def tmin_plotdir(self, operator_set_name):
    plotdir = os.path.join(self.plotdir, "tmin_plots", self.task_name, f"rebin{self.rebin}", operator_set_name)
    os.makedirs(plotdir, exist_ok=True)
    return plotdir

  def tmin_relabelled_plotdir(self, operator_set_name):
    plotdir = os.path.join(self.plotdir, "tmin_plots", self.task_name, f"rebin{self.rebin}_relabelled", operator_set_name)
    os.makedirs(plotdir, exist_ok=True)
    return plotdir

  def tmin_plotfile(self, operator_set_name, level, fit_info, show_shift, extension):
    plotdir = self.tmin_plotdir(operator_set_name)
    shift_str = "D_" if show_shift else ""
    plotfile = f"tmin_{fit_info.plotfile}_{shift_str}{level}.{extension.value}"
    return os.path.join(plotdir, plotfile)

  def tmin_relabelled_plotfile(self, operator_set_name, level, fit_info, show_shift, extension):
    plotdir = self.tmin_relabelled_plotdir(operator_set_name)
    shift_str = "D_" if show_shift else ""
    plotfile = f"tmin_{fit_info.plotfile}_{shift_str}{level}.{extension.value}"
    return os.path.join(plotdir, plotfile)

  def getSigmondInputs(self):
    sigmond_inputs = list()

    if os.path.exists(self.samplings_filename):
      os.remove(self.samplings_filename)
    
    default_data_files = self.data_handler.data_files
    # get reference observable data files
    if self.reference_fit_info is not None:
      default_data_files += self.data_handler.getChannelDataFiles(self.reference_fit_info.operator.channel)

    # get scattering observables
    for fit_info in self.scattering_particles.values():
      default_data_files += self.data_handler.getChannelDataFiles(fit_info.operator.channel)
      
    shutil.rmtree(self.fit_sh_plotdir())
    # Add the single particle energies (i.e. the ref. and the scattering particles)
    if self.reference_fit_info is not None or self.scattering_particles:
      project_name = self.project_name(self.sh_name)
      inputfile = self.inputfile(self.sh_name)
      logfile = self.logfile(self.sh_name)

      sigmond_input = self.new_sigmond_input(project_name, inputfile, logfile, default_data_files)

      energy_samplings_obs = self._addReferenceAndScatteringFits(sigmond_input, True)
      sigmond_input.writeToFile(
          self.samplings_filename, energy_samplings_obs,
          file_type=sigmond_info.sigmond_info.DataFormat.samplings, file_mode=sigmond.WriteMode.Protect)

      sigmond_input.write()

      sigmond_inputs.append(sigmond_input)

    # Now the spectra
    for operator_set, fit_infos in self.spectrum.items():
      tmin_plotdir = self.tmin_plotdir(repr(operator_set))
      fit_plotdir = self.fit_plotdir(repr(operator_set))
      shutil.rmtree(tmin_plotdir)
      shutil.rmtree(fit_plotdir)

      if operator_set.is_rotated:
        data_files = default_data_files + self.data_handler.getRotatedDataFiles(operator_set)
      else:
        data_files = default_data_files
        for operator in operator_set.operators:
          data_files += self.data_handler.getChannelDataFiles(operator.channel)

      project_name = self.project_name(repr(operator_set))
      inputfile = self.inputfile(repr(operator_set))
      logfile = self.logfile(repr(operator_set))

      sigmond_input = self.new_sigmond_input(project_name, inputfile, logfile, data_files)
      self._addReferenceAndScatteringFits(sigmond_input, False)
        
      energies = list()
      amplitudes = list()
      denergies = list()
      tmin_fit_infos = self.tmin_infos.get(operator_set, [None]*len(fit_infos))
      if not tmin_fit_infos:
        tmin_fit_infos = [None]*len(fit_infos)
    
      for level, (fit_info, tmin_fit_info) in enumerate(zip(fit_infos, tmin_fit_infos)):
        energy_obs = fit_info.energy_observable
        amplitude_obs = fit_info.amplitude_observable


        non_interacting_level = list()
        non_interacting_amp = list()
        non_interacting_scattering_fit_info = []
        if fit_info.non_interacting_operators is not None: 
          for scattering_particle in fit_info.non_interacting_operators.non_interacting_level:
            scattering_particle_fit_info = self.scattering_particles[scattering_particle]
            non_interacting_scattering_fit_info.append(scattering_particle_fit_info)
            at_rest_scattering_particle = sigmond_info.sigmond_info.ScatteringParticle(scattering_particle.name, 0, False)
            at_rest_scattering_particle_fit_info = self.scattering_particles[at_rest_scattering_particle]

            non_interacting_level.append((at_rest_scattering_particle_fit_info.energy_observable, scattering_particle.psq))
            non_interacting_amp.append(scattering_particle_fit_info.amplitude_observable)
            
        
        plotfile = self.fit_plotfile(repr(operator_set), level, util.PlotExtension.grace)
        sigmond_input.doTemporalCorrelatorFit(
            fit_info, minimizer_info=self.minimizer_info, plotfile=plotfile,
            plot_info=self.fit_plot_info, scattering_fit_info=non_interacting_scattering_fit_info)

        if fit_info.is_log_fit:
          sigmond_input.doExp(amplitude_obs, amplitude_obs, sampling_mode=self.sampling_mode)
        
            
        if tmin_fit_info is not None:
          tmin_plot_info = self.default_tmin_plot_info.setChosenFit(fit_info)

          if fit_info.ratio:
            shift_energy_obs = fit_info.energy_observable
            noshift_energy_obs = sigmond.MCObsInfo("noshift_energy", 20000)
            sigmond_input.doReconstructEnergy(
                noshift_energy_obs, shift_energy_obs, self.ensemble_spatial_extent, non_interacting_level,
                sampling_mode=self.sampling_mode)

          else:
            noshift_energy_obs = fit_info.energy_observable
            shift_energy_obs = sigmond.MCObsInfo("shift_energy", 20000)
            sigmond_input.doEnergyDifference(
                shift_energy_obs, fit_info.energy_observable, self.ensemble_spatial_extent, non_interacting_level,
                sampling_mode=self.sampling_mode)

          for shift in [True, False]:
            for a_tmin_fit_info in tmin_fit_info:
              plotfile = self.tmin_plotfile(repr(operator_set), level, a_tmin_fit_info, shift, util.PlotExtension.grace)

              if a_tmin_fit_info.ratio and shift:
                sigmond_input.doTemporalCorrelatorFit(
                    a_tmin_fit_info, minimizer_info=self.minimizer_info, plotfile=plotfile,
                    plot_info=tmin_plot_info, chosen_fit_info=shift_energy_obs, 
                    scattering_fit_info=non_interacting_scattering_fit_info)
              elif a_tmin_fit_info.ratio and not shift:
                sigmond_input.doTemporalCorrelatorFit(
                    a_tmin_fit_info, minimizer_info=self.minimizer_info, plotfile=plotfile,
                    plot_info=tmin_plot_info, chosen_fit_info=noshift_energy_obs,
                    non_interacting_level=non_interacting_level, spatial_extent=self.ensemble_spatial_extent, 
                    scattering_fit_info=non_interacting_scattering_fit_info)
              elif not a_tmin_fit_info.ratio and shift:
                sigmond_input.doTemporalCorrelatorFit(
                    a_tmin_fit_info, minimizer_info=self.minimizer_info, plotfile=plotfile,
                    plot_info=tmin_plot_info, chosen_fit_info=shift_energy_obs,
                    non_interacting_level=non_interacting_level, spatial_extent=self.ensemble_spatial_extent)
              else:
                sigmond_input.doTemporalCorrelatorFit(
                    a_tmin_fit_info, minimizer_info=self.minimizer_info, plotfile=plotfile,
                    plot_info=tmin_plot_info, chosen_fit_info=noshift_energy_obs)
        
        if fit_info.ratio:
          channel = fit_info.operator.channel
          samp_dir = f"PSQ{channel.psq}/{channel.irrep}"
          delab = sigmond.MCObsInfo(f"{samp_dir}/dElab_{level}", 0)
          sigmond_input.doCopy(delab, energy_obs, sampling_mode=self.sampling_mode)
          denergies.append(delab)
          sigmond_input.doReconstructEnergy(
              energy_obs, energy_obs, self.ensemble_spatial_extent, non_interacting_level,
              sampling_mode=self.sampling_mode)

          sigmond_input.doReconstructAmplitude(
              amplitude_obs, amplitude_obs, non_interacting_amp, sampling_mode=self.sampling_mode)
        
        energies.append(energy_obs)
        amplitudes.append(amplitude_obs)
        
      if operator_set.is_rotated:
        zfactor_plotdir = self.zfactor_plotdir(repr(operator_set))
        shutil.rmtree(zfactor_plotdir)
        pivot_name = f"pivot_{operator_set!r}"
        pivot_file = self.data_handler.pivotfile(operator_set)
        zfactor_stub = self.zfactor_plotstub(repr(operator_set))

        sigmond_input.doRotCorrMatInsertFitInfos(
            operator_set.pivot_info.pivot_type, energies, amplitudes, pivot_name=pivot_name,
            pivot_filename=pivot_file, reorder=self.reorder)

        sigmond_input.doCorrMatrixZMagSquares(
            operator_set.pivot_info.pivot_type, pivot_name=pivot_name, plot_stub=zfactor_stub,
            operators=operator_set.operators, bar_color=self.bar_color,
            obsname='standard')

        sigmond_input.getFromPivot(
            operator_set.pivot_info.pivot_type, self.ordered_energy, self.ordered_amplitude,
            pivot_name=pivot_name)

      else:
        for level, energy in enumerate(energies):
          ordered_energy = sigmond.MCObsInfo(self.ordered_energy, level)
          sigmond_input.doCopy(ordered_energy, energy, sampling_mode=self.sampling_mode)

      # Get all the energies
      energy_samplings_obs = list()
      all_fit_param_obs = list()
      for level, fit_info in enumerate(fit_infos):
        channel = fit_info.operator.channel
        samp_dir = f"PSQ{channel.psq}/{channel.irrep}"
        
        ordered_energy = sigmond.MCObsInfo(self.ordered_energy, level)
        elab = sigmond.MCObsInfo(f"{samp_dir}/elab_{level}", 0)
        ecm = sigmond.MCObsInfo(f"{samp_dir}/ecm_{level}", 0)

        sigmond_input.doCopy(elab, ordered_energy, sampling_mode=self.sampling_mode)
        sigmond_input.doBoost(ecm, channel.psq, self.ensemble_spatial_extent,
                              elab, sampling_mode=self.sampling_mode)

        energy_samplings_obs.append(elab)
        energy_samplings_obs.append(ecm)
        for i in range(fit_info.num_params):
            all_fit_param_obs.append( fit_info.fit_param_obs(i) )

      
        if self.reference_fit_info is not None:
          elab_ref = sigmond.MCObsInfo(f"{samp_dir}/elab_{level}_{self.ref_name}", 0)
          ecm_ref = sigmond.MCObsInfo(f"{samp_dir}/ecm_{level}_{self.ref_name}", 0)

          sigmond_input.doRatio(elab_ref, elab, self.reference_fit_info.energy_observable,
                                sampling_mode=self.sampling_mode)
          sigmond_input.doRatio(ecm_ref, ecm, self.reference_fit_info.energy_observable,
                                sampling_mode=self.sampling_mode)

          energy_samplings_obs.append(elab_ref)
          energy_samplings_obs.append(ecm_ref)
            
#       energy_samplings_obs.extend(denergies)

      sigmond_input.writeToFile(
          self.samplings_filename, energy_samplings_obs,
          file_type=sigmond_info.sigmond_info.DataFormat.samplings, file_mode=sigmond.WriteMode.Protect)
      
      param_tag = f"{channel.irrep}({channel.psq})"
      if os.path.exists(self.params_filename(param_tag)):
        os.remove(self.params_filename(param_tag))
    
      if self.write_params_to_file:
        sigmond_input.writeToFile(
          self.params_filename(param_tag,True), all_fit_param_obs,
          file_type=sigmond_info.sigmond_info.DataFormat.samplings, file_mode=sigmond.WriteMode.Overwrite,
          hdf5 = True
        )

      sigmond_input.write()
      sigmond_inputs.append(sigmond_input)

    return sigmond_inputs

  def _addReferenceAndScatteringFits(self, sigmond_input, plot=True):
    energy_samplings_obs = list()
    all_fit_param_obs = list()
    if self.reference_fit_info is not None:
      if plot:
        plotfile = self.fit_sh_plotfile(self.ref_name, util.PlotExtension.grace)
        sigmond_input.doTemporalCorrelatorFit(
            self.reference_fit_info, minimizer_info=self.minimizer_info, plotfile=plotfile,
            plot_info=self.fit_plot_info)
        for i in range(self.reference_fit_info.num_params):
            all_fit_param_obs.append( self.reference_fit_info.fit_param_obs(i) )
      else:
        sigmond_input.doTemporalCorrelatorFit(self.reference_fit_info, minimizer_info=self.minimizer_info)

      ref_info = sigmond.MCObsInfo(f"{self.sh_name}/{self.ref_name}", 0)
      sigmond_input.doCopy(ref_info, self.reference_fit_info.energy_observable,
                             sampling_mode=self.sampling_mode)
      energy_samplings_obs.append(ref_info)

    sh_tmin_fit_infos = self.tmin_infos.get(self.sh_name)
    for scattering_particle, scat_fit_info in self.scattering_particles.items():
      if scat_fit_info != self.reference_fit_info and plot:
        plotfile = self.fit_sh_plotfile(repr(scattering_particle), util.PlotExtension.grace)
        sigmond_input.doTemporalCorrelatorFit(
            scat_fit_info, minimizer_info=self.minimizer_info, plotfile=plotfile,
            plot_info=self.fit_plot_info)
        
        for i in range(scat_fit_info.num_params):
            all_fit_param_obs.append( scat_fit_info.fit_param_obs(i) )

        for tmin_fit_info in sh_tmin_fit_infos[scattering_particle]:
          tmin_plot_info = self.default_tmin_plot_info.setChosenFit(scat_fit_info)
          tmin_plotfile = self.tmin_sh_plotfile(tmin_fit_info, util.PlotExtension.grace)
          sigmond_input.doTemporalCorrelatorFit(
              tmin_fit_info, minimizer_info=self.minimizer_info, plotfile=tmin_plotfile,
              plot_info=tmin_plot_info, chosen_fit_info=scat_fit_info.energy_observable)

      elif scat_fit_info != self.reference_fit_info:
        sigmond_input.doTemporalCorrelatorFit(scat_fit_info, minimizer_info=self.minimizer_info)
      elif plot:
        ref_file = self.fit_sh_plotfile(self.ref_name, util.PlotExtension.grace)
        scat_file = self.fit_sh_plotfile(repr(scattering_particle), util.PlotExtension.grace)
        self.copy_file = (ref_file, scat_file)

        for tmin_fit_info in sh_tmin_fit_infos[scattering_particle]:
          tmin_plot_info = self.default_tmin_plot_info.setChosenFit(scat_fit_info)
          tmin_plotfile = self.tmin_sh_plotfile(tmin_fit_info, util.PlotExtension.grace)
          sigmond_input.doTemporalCorrelatorFit(
              tmin_fit_info, minimizer_info=self.minimizer_info, plotfile=tmin_plotfile,
              plot_info=tmin_plot_info, chosen_fit_info=scat_fit_info.energy_observable)


      scat_info = sigmond.MCObsInfo(f"{self.sh_name}/{scattering_particle!r}", 0)
      sigmond_input.doCopy(scat_info, scat_fit_info.energy_observable, sampling_mode=self.sampling_mode)
      energy_samplings_obs.append(scat_info)

      if self.reference_fit_info is not None:
        scat_info_ref = sigmond.MCObsInfo(f"{self.sh_name}/{scattering_particle!r}_{self.ref_name}", 0)
        sigmond_input.doRatio(scat_info_ref, scat_info, self.reference_fit_info.energy_observable,
                              sampling_mode=self.sampling_mode)
        energy_samplings_obs.append(scat_info_ref)

    if plot:
      #delete param file
      if os.path.exists(self.params_filename("sh")):
        os.remove(self.params_filename("sh"))
      sigmond_input.writeToFile(
          self.params_filename("sh",True), all_fit_param_obs,
          file_type=sigmond_info.sigmond_info.DataFormat.samplings, file_mode=sigmond.WriteMode.Overwrite,
          hdf5 = True
      )
        
    return energy_samplings_obs


  def finalize(self):
    self.spectrum_logs = dict()
    
    logfile = self.logfile('single_hadrons')
    spectrum_log = sigmond_info.sigmond_log.SpectrumLog(logfile)
    #self.spectrum_logs['single_hadrons'] = spectrum_log
    
    for operator_set in self.spectrum.keys():
      if operator_set.is_rotated:
        logfile = self.logfile(repr(operator_set))
        spectrum_log = sigmond_info.sigmond_log.SpectrumLog(logfile)
        self.spectrum_logs[operator_set] = spectrum_log

    self._setEnergies()
    if self.write_params_to_file:
      self._combine_param_files()
      self._include_fit_model()
    doc = util.create_doc(f"Spectrum: {self.ensemble_name} - {self.task_name}", True)
    if self.spectrum:
      with doc.create(pylatex.Section("Spectrum Results")):
        spectrum_tikz_file = self._makePlot()

        doc.append(pylatex.NoEscape(r"\resizebox{\textwidth}{!}{"))
        doc.append(pylatex.NoEscape(rf"\input{{{spectrum_tikz_file}.tikz}}"))
        doc.append(pylatex.NoEscape(r"}"))
        self._addFitTable(doc)
        #if non interacting levels add Noninteracting levels fit table
        doc.append(pylatex.NoEscape(r"\newpage"))

    self._addFittedEnergies(doc)
    doc.append(pylatex.NoEscape(r"\newpage"))
    self._addTminPlots(doc)
    doc.append(pylatex.NoEscape(r"\newpage"))

    if self.spectrum:
      self._addZfactors(doc)

    filename = os.path.join(self.results_dir, f"spectrum_{self.ensemble_name}_{self.task_name}_rebin{self.rebin}")
    util.compile_pdf(doc, filename, self.latex_compiler)
    self._reorderTminFitPlots()
          
  def _addZfactors(self, doc):
    with doc.create(pylatex.Section("Overlap factors")):
      for operator_set, spectrum_log in self.spectrum_logs.items():
        zfactor_dir = self.zfactor_plotdir(repr(operator_set))
        util.dirGrace2pdf(zfactor_dir)
        grouped_ops = [operator_set.operators[n:n+3] for n in range(0, len(operator_set.operators), 3)]
        with doc.create(pylatex.Subsection(str(operator_set))):
          for ops in grouped_ops:
            with doc.create(pylatex.Figure(position='H')) as fig:
              for op in ops:
                zfactor_plot = self.zfactor_plotfile(repr(operator_set), op, util.PlotExtension.pdf)
                util.add_image(fig, self.results_dir, zfactor_plot, width="0.3", view=False)

  def _addTminPlots(self, doc):
    with doc.create(pylatex.Section(pylatex.NoEscape(r"$t_{\rm min}$ plots"))):
      for operator_set, tmin_fit_infos in self.tmin_infos.items():
        # This is a stupid hack to get the sh particles...
        if isinstance(operator_set, str):
          if operator_set != self.sh_name:
            logging.critical("What?")

          util.dirGrace2pdf(self.tmin_sh_plotdir())
          tmin_fit_infos = self.tmin_infos[self.sh_name]
          with doc.create(pylatex.Subsection("Scattering Particles")):
            for scattering_particle, tmin_fit_info in tmin_fit_infos.items():
              list_of_tmin_plots = list()
              for fit_info in tmin_fit_info:
                plotfile = self.tmin_sh_plotfile(fit_info, util.PlotExtension.pdf)
                caption = fit_info.operator.op_str().replace('_', '\_')
                caption += f"\\newline {fit_info.model.short_name}"
                caption += f", $t_{{\\rm max}} = {fit_info.tmax}$"
                if fit_info.noise_cutoff > 0.:
                  caption += f" ($\\sigma_{{\\rm cut}} = {round(fit_info.noise_cutoff, 2)}$)"
                list_of_tmin_plots.append((plotfile, caption))

              grouped_plots = [list_of_tmin_plots[n:n+3] for n in range(0, len(list_of_tmin_plots), 3)]
              if grouped_plots:
                with doc.create(pylatex.Subsubsection(str(scattering_particle))):
                  for plot_group in grouped_plots:
                    with doc.create(pylatex.Figure(position='H')) as subfig:
                      for plot in plot_group:
                        plot_file = plot[0]
                        cap = plot[1]
                        if len(plot_group) == 1:
                          util.add_image(subfig, self.results_dir, plot_file, width='0.33', caption=cap)
                        else:
                          with doc.create(pylatex.SubFigure(position='b', width=pylatex.NoEscape(r'0.33\linewidth'))) as fig:
                            util.add_image(fig, self.results_dir, plot_file, width='1.0', caption=cap)
        else:
          util.dirGrace2pdf(self.tmin_plotdir(repr(operator_set)))
          with doc.create(pylatex.Subsection(str(operator_set))):
            for level, tmin_fit_info in enumerate(tmin_fit_infos):
              list_of_tmin_plots = list()
              list_of_tmin_shift_plots = list()
              for shift in [True, False]:
                for fit_info in tmin_fit_info:
                  plotfile = self.tmin_plotfile(repr(operator_set), level, fit_info, shift, util.PlotExtension.pdf)
                  caption = fit_info.operator.op_str().replace('_', '\_')
                  if shift:
                    caption += r"\newline $a_t \Delta E_{\rm lab}$, "
                  else:
                    caption += r"\newline $a_t E_{\rm lab}$, "
                  caption += fit_info.model.short_name
                  if fit_info.ratio:
                    caption += " - ratio"
                  caption += f", $t_{{\\rm max}} = {fit_info.tmax}$"
                  if fit_info.noise_cutoff > 0.:
                    caption += f" ($\\sigma_{{\\rm cut}} = {round(fit_info.noise_cutoff, 2)}$)"
                  if fit_info.model==sigmond_info.fit_info.FitModel.TimeForwardDoubleExpRatio and not shift:
                    continue
                  if shift:
                    list_of_tmin_shift_plots.append((plotfile, caption))
                  else:
                    list_of_tmin_plots.append((plotfile, caption))

              #add other plots to this list
              sim_plots = []
              pattern = rf'^tmin_fit_\S+-ROT-{level}_\S+_D_\S+-(?<particle_name>\S+)\[(?P<spatial>[a-zA-Z]+)(?P<spatial_id>\d+)\]-0.pdf$'
              for test_str in os.listdir(self.tmin_plotdir(repr(operator_set))):
                  match = regex.match(pattern, test_str)
                  if match:
                      matches = match.groupdict()
                      matches['file'] = os.path.join(self.tmin_plotdir(repr(operator_set)),test_str)
                      matches['level'] = level
                      sim_plots.append(matches)
                    
              if sim_plots:
                sim_plots.sort( key=util.sort_plots )
                for item in sim_plots:
                  list_of_tmin_shift_plots.append( (item['file'], fr"\newline $a_t \Delta E_{{\rm lab}}$, {sigmond_info.fit_info.FitModel.TimeForwardDoubleExpRatio.short_name} {item['particle_name']} fit" ) )
                                             
              grouped_shift_plots = [list_of_tmin_shift_plots[n:n+3] for n in range(0, len(list_of_tmin_shift_plots), 3)]
              grouped_plots = [list_of_tmin_plots[n:n+3] for n in range(0, len(list_of_tmin_plots), 3)]
              if grouped_plots:
                with doc.create(pylatex.Subsubsection(f"ROT {level}")):
                  for plot_group in grouped_shift_plots:
                    with doc.create(pylatex.Figure(position='H')) as subfig:
                      for plot in plot_group:
                        plot_file = plot[0]
                        cap = plot[1]
                        if len(plot_group) == 1:
                          util.add_image(subfig, self.results_dir, plot_file, width='0.33', caption=cap)
                        else:
                          with doc.create(pylatex.SubFigure(position='b', width=pylatex.NoEscape(r'0.33\linewidth'))) as fig:
                            util.add_image(fig, self.results_dir, plot_file, width='1.0', caption=cap)

                  doc.append(pylatex.NoEscape(r"\hrule"))
                  for plot_group in grouped_plots:
                    with doc.create(pylatex.Figure(position='H')) as subfig:
                      for plot in plot_group:
                        plot_file = plot[0]
                        cap = plot[1]
                        if len(plot_group) == 1:
                          util.add_image(subfig, self.results_dir, plot_file, width='0.33', caption=cap)
                        else:
                          with doc.create(pylatex.SubFigure(position='b', width=pylatex.NoEscape(r'0.33\linewidth'))) as fig:
                            util.add_image(fig, self.results_dir, plot_file, width='1.0', caption=cap)

  def _reorderTminFitPlots(self):
    for operator_set, tmin_fit_infos in self.tmin_infos.items():
        if isinstance(operator_set, str):
          if operator_set != self.sh_name:
            logging.critical("What?")
        else:
            new_level = 0
            for level, energy in self.energies[operator_set].items():
                tmin_fit_info = tmin_fit_infos[level]
                for fit_info in tmin_fit_info:
                    for shift in [True, False]:
                        plotfile = self.tmin_plotfile(repr(operator_set), level, fit_info, shift, util.PlotExtension.grace)
                        new_plotfile = self.tmin_relabelled_plotfile(repr(operator_set), new_level, fit_info, shift, util.PlotExtension.grace)
                        shutil.copyfile(plotfile,new_plotfile)
                new_level+=1

  def _addFittedEnergies(self, doc):
    with doc.create(pylatex.Section("Fitted Energies")):
      if hasattr(self, "copy_file"):
        shutil.copy(self.copy_file[0], self.copy_file[1])
      util.dirGrace2pdf(self.fit_sh_plotdir())

      # first the reference
      if self.reference_fit_info is not None:
        plotfile = self.fit_sh_plotfile(self.ref_name, util.PlotExtension.pdf)
        with doc.create(pylatex.Subsection("Reference Energy")):
          with doc.create(pylatex.Figure(position='H')) as fig:
            cap = self.reference_fit_info.operator.op_str().replace('_', '\_')
            util.add_image(fig, self.results_dir, plotfile, width='0.5', caption=cap)

      # now scattering particles
      if self.scattering_particles:
        with doc.create(pylatex.Subsection("Scattering Energies")):
          scattering_particles = list(self.scattering_particles.keys())

          grouped_scattering_particles = [scattering_particles[n:n+2] for n in range(0, len(scattering_particles), 2)]
          for scattering_particles_group in grouped_scattering_particles:
            with doc.create(pylatex.Figure(position='H')) as subfig:
              for scattering_particle in scattering_particles_group:
                plot_file = self.fit_sh_plotfile(repr(scattering_particle), util.PlotExtension.pdf)
                cap = self.scattering_particles[scattering_particle].operator.op_str().replace('_', '\_')
                if len(scattering_particles_group) == 1:
                  util.add_image(subfig, self.results_dir, plot_file, width='0.5', caption=cap)
                else:
                  with doc.create(pylatex.SubFigure(position='b', width=pylatex.NoEscape(r'0.5\linewidth'))) as fig:
                    util.add_image(fig, self.results_dir, plot_file, width='1.0', caption=cap)



      # And spectrum
      for operator_set, fit_infos in self.spectrum.items():
        util.dirGrace2pdf(self.fit_plotdir(repr(operator_set)))
        with doc.create(pylatex.Subsection(str(operator_set))):
          levels = sorted(list(self.energies[operator_set].keys()))
          grouped_levels = [levels[n:n+3] for n in range(0, len(levels), 3)]
          for grouped_level in grouped_levels:
            with doc.create(pylatex.Figure(position='H')) as subfig:
              for level in grouped_level:
                plot_file = self.fit_plotfile(repr(operator_set), level, util.PlotExtension.pdf)
                cap = fit_infos[level].operator.op_str().replace('_', '\_')
                if len(grouped_level) == 1:
                  util.add_image(subfig, self.results_dir, plot_file, width='0.33', caption=cap)
                else:
                  with doc.create(pylatex.SubFigure(position='b', width=pylatex.NoEscape(r'0.33\linewidth'))) as fig:
                    util.add_image(fig, self.results_dir, plot_file, width='1.0', caption=cap)
        #collect other graphs
#         pattern = r'fit_[0-9]-[A-Za-z0-9]*\.agr' #fit_3-3-1p1G1-S[SS0]-0.agr
        
        sim_plots = []
        pattern = r'^fit_(?P<level>\d+)-\S+-(?<particle_name>\S+)\[(?P<spatial>[a-zA-Z]+)(?P<spatial_id>\d+)\]-0.pdf$' #p(?P<psq>\d+)(?P<irrep>\S+)
        for test_str in os.listdir(self.fit_plotdir(repr(operator_set))):
            match = regex.match(pattern, test_str)
            if match:
                matches = match.groupdict()
                matches['file'] = os.path.join(self.fit_plotdir(repr(operator_set)),test_str)
                sim_plots.append(matches)
#                 print( test_str, matches )
          
        if sim_plots:
            sim_plots.sort( key=util.sort_plots )
            grouped_sim_plots = [sim_plots[n:n+3] for n in range(0, len(sim_plots), 3)]
            for group in grouped_sim_plots:
              with doc.create(pylatex.Figure(position='H')) as subfig:
                  for plot in group:
                    cap = f"ROT {plot['level']}, {plot['particle_name']} fit"
                    if len(group) == 1:
                      util.add_image(subfig, self.results_dir, plot['file'], width='0.33', caption=cap)
                    else:
                      with doc.create(pylatex.SubFigure(position='b', width=pylatex.NoEscape(r'0.33\linewidth'))) as fig:
                        util.add_image(fig, self.results_dir, plot['file'], width='1.0', caption=cap)


  def _addFitTable(self, doc):
    #make table for spectrum fits
    with doc.create(pylatex.Center()) as centered:
      long_tabu = "X[0.1,c] X[0.5,c] X[0.5,c] X[2,c] X[3,c]  X[4,c] X[3,c] X[4,c] X[3,c] X[1.5,c] X[2,c]"
      with centered.create(pylatex.LongTabu(long_tabu, to=r"\linewidth")) as data_table:

        if self.reference_fit_info is None:
          ref_energy_header = r"$a_t E_{\rm cm}$"
        else:
          ref_latex = self.latex_map.get(self.ref_name, self.ref_name)
          ref_energy_header = fr"$E_{{\rm cm}} / E_{{{ref_latex}}}$"

        header_row = [
            pylatex.NoEscape(""),
            pylatex.NoEscape(r"$d^2$"),
            pylatex.NoEscape(r"$\Lambda$"),
            pylatex.NoEscape(r"Level"),
            pylatex.NoEscape(ref_energy_header),
            pylatex.NoEscape(r"$a_t (\Delta) E_{\rm lab}$"),
            pylatex.NoEscape(r"$a_t E_{\rm lab}$"),
            pylatex.NoEscape("Fit model"),
            pylatex.NoEscape(r"$(t_{\rm min}, t_{\rm max})$"),
#             pylatex.NoEscape(r"$\chi^2/\text{dof}$"),
            pylatex.NoEscape(r"$p$-val."),
            pylatex.NoEscape(r"NI"),
        ]

        data_table.add_row(header_row, mapper=[pylatex.utils.bold])
        data_table.end_table_header()

        for operator_set, fit_infos in self.spectrum.items():
          data_table.add_hline()
          new_level = 0
          for level, energy in self.energies[operator_set].items():
            flag = ""
            if operator_set in self.flagged_levels:
              flag = "X" if self.flagged_levels[operator_set][level] else ""

            fit_info = fit_infos[level]
            noninteracting_level = ""
            if fit_info.ratio:
                if fit_info.non_interacting_operators is not None:
                  for scattering_particle in fit_info.non_interacting_operators.non_interacting_level:
                    noninteracting_level+=str(scattering_particle)
            channel = operator_set.channel
            logfile = self.logfile(repr(operator_set))
            fit_log = sigmond_info.sigmond_log.FitLog(logfile)

            fit_result = fit_log.fits[fit_info]

            energy = util.nice_value(energy.getFullEstimate(), energy.getSymmetricError())
            at_energy = fit_result.energy
            recon_energy = fit_result.reconstructed_energy
            fit_model = fit_info.model.short_name
            if fit_info.ratio:
              fit_model += " - ratio"

            irrep = channel.irrep
            if irrep in self.latex_map:
              irrep = pylatex.NoEscape(rf"${self.latex_map[irrep]}$")

            data_row = [
                flag,
                channel.psq,
                irrep,
                pylatex.NoEscape(rf"${level} \to {new_level}$"),
                energy,
                at_energy,
                recon_energy,
                fit_model,
                pylatex.NoEscape(rf"$({fit_info.tmin}, {fit_info.tmax})$"),
#                 round(fit_result.chisq,2),
                round(fit_result.quality,2),
                noninteracting_level,
            ]

            data_table.add_row(data_row)
            new_level += 1
    
    #make table for scattering particles
    if self.scattering_particles:
      with doc.create(pylatex.Center()) as centered:
        long_tabu = "X[8,c] X[4,c] X[4,c] X[4,c] X[3,c] X[2,c] X[2,c]"
        with centered.create(pylatex.LongTabu(long_tabu, to=r"\linewidth")) as data_table2:
          header_row = [
                pylatex.NoEscape(r"Scattering Particles($d^2$)"),
                pylatex.NoEscape(ref_energy_header),
                pylatex.NoEscape(r"$a_t (\Delta) E_{\rm lab}$"),
                pylatex.NoEscape("Fit model"),
                pylatex.NoEscape(r"$(t_{\rm min}, t_{\rm max})$"),
                pylatex.NoEscape(r"$\chi^2/\text{dof}$"),
                pylatex.NoEscape(r"$p$-val."),
          ]

          data_table2.add_row(header_row, mapper=[pylatex.utils.bold])
          data_table2.end_table_header()
          
          data_table2.add_hline()
          new_level = 0
          scattering_particles = list(self.scattering_particles.keys())
          grouped_scattering_particles = [scattering_particles[n:n+2] for n in range(0, len(scattering_particles), 2)]
          logfile = self.logfile('single_hadrons')
          fit_log = sigmond_info.sigmond_log.FitLog(logfile)
          for scattering_particles_group in grouped_scattering_particles:
            for scattering_particle in scattering_particles_group:
              fit_info = self.scattering_particles[scattering_particle]
              energy = self.ni_energies[scattering_particle]
              energy = util.nice_value(energy.getFullEstimate(), energy.getSymmetricError())
              fit_result = fit_log.fits[fit_info]
              at_energy = fit_result.energy
              fit_model = fit_info.model.short_name
              if fit_info.ratio:
                fit_model += " - ratio"
              data_row = [
                      scattering_particle,
                      energy,
                      at_energy,
                      fit_model,
                      pylatex.NoEscape(rf"$({fit_info.tmin}, {fit_info.tmax})$"),
                      round(fit_result.chisq,2),
                      round(fit_result.quality,2),
              ]
    
              data_table2.add_row(data_row)
              new_level += 1

  def _combine_param_files(self):
    
    if os.path.exists(self.params_filename()):
        os.remove(self.params_filename())
        
    outfile = h5py.File(self.params_filename(), "a")
    #single hadrons
    if os.path.exists(self.params_filename("sh")):
        hdf5_h = h5py.File(self.params_filename("sh"),"r+")
        outfile.create_group('Info')
        for info in hdf5_h['Info'].keys():
            data=hdf5_h['Info'][info][()]
            outfile['Info'].create_dataset(info,data.shape,data=data )
            
        if 'params' not in outfile.keys():
            outfile.create_group('params')
        for irrep in hdf5_h['params'].keys():
            if type(hdf5_h['params'][irrep])==h5py._hl.dataset.Dataset:
                if irrep not in outfile['params'].keys():
                    data=hdf5_h['params'][irrep][()]
                    outfile['params'].create_dataset(irrep,data.shape,data=data)
            else:
                if irrep not in outfile['params'].keys():
                    outfile['params'].create_group(irrep)
                for value in hdf5_h['params'][irrep].keys():
                    data=hdf5_h['params'][irrep][value][()]
                    outfile['params'][irrep].create_dataset(value,data.shape,data=data)
            
        hdf5_h.close()
        os.remove(self.params_filename("sh"))
        
    for operator_set, fit_infos in self.spectrum.items():
        channel = fit_infos[0].operator.channel
        param_tag = f"{channel.irrep}({channel.psq})"
        hdf5_h = h5py.File(self.params_filename(param_tag),"r+")
        if 'params' not in outfile.keys():
            outfile.create_group('params')
        for irrep in hdf5_h['params'].keys():
            if type(hdf5_h['params'][irrep])==h5py._hl.dataset.Dataset:
                if irrep not in outfile['params'].keys():
                    data=hdf5_h['params'][irrep][()]
                    outfile['params'].create_dataset(irrep,data.shape,data=data)
            else:
                if irrep not in outfile['params'].keys():
                    outfile['params'].create_group(irrep)
                for value in hdf5_h['params'][irrep].keys():
                    data=hdf5_h['params'][irrep][value][()]
                    outfile['params'][irrep].create_dataset(value,data.shape,data=data)
        hdf5_h.close()
        os.remove(self.params_filename(param_tag))
    outfile.close()
            
  def _include_fit_model(self):
    #include single hadrons
    hdf5_h = h5py.File(self.params_filename(),"r+")
    for operator_set, fit_infos in self.spectrum.items():
      
      param_info = hdf5_h.create_group("param_info")
      logfile = self.spectrum_logs[operator_set].logfile
      log_xml_root = ET.parse(logfile).getroot()
      for item in log_xml_root.findall('./Task/DoFit'):
        if item.find("./Type").text =="TemporalCorrelator" and item.find("./TemporalCorrelatorFit/Model").text =="TimeForwardMultiExponential":
          corr = item.find('./TemporalCorrelatorFit/GIOperatorString').text
          if corr not in param_info.keys():
              param_info_this_corr = param_info.create_group(corr)
              if item.findall('./BestFitResult'):
                  fit_level = int(item.find('./BestFitResult/FitLevel').text)
                  final_tmin = int(item.find('./BestFitResult/FinalTmin').text)
                  chisqr = float(item.find('./BestFitResult/ChiSquarePerDof').text)
                  final_tmax = int(item.find('./TemporalCorrelatorFit/TimeSeparations').text.split(" ")[-1])
                  n = []
                  try:
                    for i in range(2,6):
                        n.append(float(item.find(f'./BestFitResult/N{i}').text))
                  except:
                    pass

                  mc_observables = [mcobs.text for mcobs in item.findall('./BestFitResult/*/MCObservable/Info')]
                  param_info_this_corr.attrs.create("FitLevel",fit_level)
                  param_info_this_corr.attrs.create("FinalTmin",final_tmin)
                  param_info_this_corr.attrs.create("FinalTmax",final_tmax)
                  param_info_this_corr.attrs.create("ChiSquarePerDof",chisqr)
                  for i,this_n in enumerate(n):
                      param_info_this_corr.attrs.create(f"N{i+2}",this_n)
                  param_info_this_corr.attrs.create("FitParams",mc_observables)
    
    hdf5_h.close()
    

  def _setEnergies(self):
    # output energies to file
    data_files = data_handling.data_files.DataFiles()
    data_files.addSamplingFiles(self.samplings_filename)
    obs_handler, _ = util.get_obs_handlers(data_files, self.bins_info, self.sampling_info)

    samplings_handler = sigmond.SamplingsGetHandler(
        self.bins_info, self.sampling_info, set([self.samplings_filename]))

#     hdf5_filename = self.hdf5_filename
#     if os.path.exists(hdf5_filename):
#       os.remove(hdf5_filename)
#     hdf5_h = h5py.File(hdf5_filename, 'w')
    
    est_filename = self.estimates_filename
    fests = open(est_filename, 'w+')
    fests.write(f"obs,val,err\n")
    
    #print non interacting levels for each level to their correponding irrp into the csv
    for operator_set, fit_infos in self.spectrum.items():
      non_interacting_energies,non_interacting_energy_names = self._getNonInteractingLevels(operator_set)
      for ni_name, ni_energy in zip(non_interacting_energy_names,non_interacting_energies):
        obsname = f"PSQ{operator_set.channel.psq}/{operator_set.channel.irrep}/"
        for particle, mom in zip(ni_name[0],ni_name[1]):
          obsname+=f"{particle}({mom})"
        obsname+="_ref"
        fests.write(f"{obsname},{ni_energy.getFullEstimate()},{ni_energy.getSymmetricError()}\n")

    for obs_info in samplings_handler.getKeys():
      np_data = util.get_samplings(obs_handler, obs_info)
      val = obs_handler.getEstimate(obs_info).getFullEstimate()
      err = obs_handler.getEstimate(obs_info).getSymmetricError()
      fests.write(f"{obs_info.getObsName()},{val},{err}\n")

#     hdf5_h.close()
    fests.close()
    
    #write qsqr to file
    if len(self.single_hadrons)==2:
        sh1ref_obs_info = sigmond.MCObsInfo(f'single_hadrons/{self.single_hadrons[0]}(0)_ref',0)
        sh2ref_obs_info = sigmond.MCObsInfo(f'single_hadrons/{self.single_hadrons[1]}(0)_ref',0)
        if sh1ref_obs_info in samplings_handler.getKeys() and sh2ref_obs_info in samplings_handler.getKeys():
            qsqr_filename = self.qsqr_filename
            if os.path.exists(qsqr_filename):
              os.remove(qsqr_filename)
            hdf5_qsqr = h5py.File(qsqr_filename, 'w')
            sh1ref_data = util.get_samplings(obs_handler, sh1ref_obs_info)
            sh2ref_data = util.get_samplings(obs_handler, sh2ref_obs_info)
            for obs_info in samplings_handler.getKeys():
                if '_ref' in obs_info.getObsName() and 'cm' in obs_info.getObsName():
                    np_data = util.get_samplings(obs_handler, obs_info)
                    qsqr_data = np_data*np_data/4.0-(sh1ref_data*sh1ref_data+sh2ref_data*sh2ref_data)/2.0 + (sh1ref_data*sh1ref_data-sh2ref_data*sh2ref_data)*(sh1ref_data*sh1ref_data-sh2ref_data*sh2ref_data)/(4.0*np_data*np_data)
                    hdf5_qsqr.create_dataset(obs_info.getObsName().replace('ecm','q2cm'), data=qsqr_data)

            hdf5_qsqr.close()
    elif self.single_hadrons:
        logging.warning("q^2 calculations are only set up for two hadron correlators") 

    # save energies to self
    self.energies = dict()
    for operator_set, fit_infos in self.spectrum.items():
      self.energies[operator_set] = dict()
      if operator_set.is_rotated:
        for level, fit_info in self.spectrum_logs[operator_set].energies.items():
          operator = fit_info.operator
          if operator.channel.irrep not in operator_set.channel.irrep:
            continue
          obs_name = f"PSQ{operator.psq}/{operator.channel.irrep}/ecm_{level.new}"
          if self.reference_fit_info is not None:
            obs_name += f"_{self.ref_name}"

          try:
            energy = obs_handler.getEstimate(sigmond.MCObsInfo(obs_name, 0))
          except RuntimeError as error:
            logging.warning(f"{error}")
#             logging.critical(f"Failed to get samplings for {obs_name}")
            logging.warning(f"Failed to get samplings for {obs_name}")
            continue

          self.energies[operator_set][level.original] = energy

      else:
        for level, fit_info in enumerate(fit_infos):
          operator = fit_info.operator
          obs_name = f"PSQ{operator.psq}/{operator.channel.irrep}/ecm_{level}"
          if self.reference_fit_info is not None:
            obs_name += f"_{self.ref_name}"

          try:
            energy = obs_handler.getEstimate(sigmond.MCObsInfo(obs_name, 0))
          except RuntimeError as error:
            logging.warning(f"{error}")
            logging.critical(f"Failed to get samplings for {obs_name}")

          self.energies[operator_set][level] = energy

    samplings_handler = sigmond.SamplingsGetHandler(
        self.bins_info, self.sampling_info, set([self.samplings_filename]))

    hdf5_filename = self.hdf5_filename
    if os.path.exists(hdf5_filename):
      os.remove(hdf5_filename)
    hdf5_h = h5py.File(hdf5_filename, 'w')

    for obs_info in samplings_handler.getKeys():
      np_data = util.get_samplings(obs_handler, obs_info)
      hdf5_h.create_dataset(obs_info.getObsName(), data=np_data)

    # Add non-interacting level
    for operator_set, fit_infos in self.spectrum.items():
      channel = operator_set.channel
      group_name = f"/PSQ{channel.psq}/{channel.irrep}"
      group = hdf5_h[group_name]
      
      # get reorder map
      reorder = dict()
      for new_level, level in enumerate(self.energies[operator_set].keys()):
        reorder[level] = new_level

      # old
      '''
      particles = set()
      moms = [None]*len(fit_infos)
      for level, fit_info in enumerate(fit_infos):
        particle_set = list()
        mom_set = list()

        for scattering_particle in fit_info.non_interacting_operators.non_interacting_level.particles:
          particle_set.append(scattering_particle.name)
          mom_set.append(scattering_particle.psq)
        
        particles.add(tuple(particle_set))
        moms[reorder[level]] = mom_set

      if len(particles) != 1:
        logging.critical("Ordering of particles is not constant")

      group.attrs.create('particles', list(particles.pop()))
      group.attrs.create('free_levels', moms)
      '''
      # new
      try:
        free_levels = [None]*len(fit_infos)
        for level, fit_info in enumerate(fit_infos):
          free_levels_set = list()

          for scattering_particle in fit_info.non_interacting_operators.non_interacting_level.particles:
            free_levels_set.append(str(scattering_particle))
          
          free_levels[reorder[level]] = free_levels_set

        group.attrs.create('free_levels', free_levels)
      except Exception as e:
        logging.warning(f"Could not add free levels: {e}")


    hdf5_h.close()

    #save single hadrons to self
    if self.scattering_particles:
      self.ni_energies = dict()
      for operator, fit_info in self.scattering_particles.items():
        #operator = fit_info.operator
        obs_name = f"single_hadrons/{operator}"
        if self.reference_fit_info is not None:
          obs_name += f"_{self.ref_name}"
        try:
          energy = obs_handler.getEstimate(sigmond.MCObsInfo(obs_name, 0))
        except RuntimeError:
          logging.critical(f"Failed to get samplings for {obs_name}")
        self.ni_energies[operator] = energy

  def _getThresholds(self):
    data_files = data_handling.data_files.DataFiles()
    data_files.addSamplingFiles(self.samplings_filename)
    obs_handler, _ = util.get_obs_handlers(data_files, self.bins_info, self.sampling_info)

    returned_thresholds = dict()
    for threshold in self.thresholds:
      obs_infos = list()
      coeffs = list()
      latexs = list()
      for scattering_particle_name in threshold:
        scattering_particle = sigmond_info.sigmond_info.ScatteringParticle(scattering_particle_name, 0, False)
        energy_obs = sigmond.MCObsInfo(f"{self.sh_name}/{scattering_particle!r}", 0)
        obs_infos.append(energy_obs)
        coeffs.append(1.0)
        latexs.append(self.latex_map.get(scattering_particle_name, scattering_particle.name))

      threshold_obs = util.linear_superposition_obs(obs_handler, obs_infos, coeffs)

      if self.reference_fit_info is not None:
        ref_obs = sigmond.MCObsInfo(f"{self.sh_name}/{self.ref_name}", 0)
        threshold_obs = util.ratio_obs(obs_handler, threshold_obs, ref_obs)

      threshold_latex = " ".join(latexs)
      threshold_latex = rf"${threshold_latex}$"
      returned_thresholds[threshold_latex] = obs_handler.getEstimate(threshold_obs)

    return returned_thresholds

      
  def _getNonInteractingLevels(self, operator_basis):
    data_files = data_handling.data_files.DataFiles()
    data_files.addSamplingFiles(self.samplings_filename)
    obs_handler, _ = util.get_obs_handlers(data_files, self.bins_info, self.sampling_info)

    non_interacting_levels = list()      #    level 0     ,   level 1...
    non_interacting_level_names = list() #[ [[N,pi].[0,0]], [[N,pi].[1,0]], ...
    for fit_info in self.spectrum[operator_basis]:
      single_particle_names = list()
      single_particle_moms = list()
      single_particles = list()

      if fit_info.non_interacting_operators is None:
        continue

      for scattering_particle in fit_info.non_interacting_operators.non_interacting_level:
        at_rest_scattering_particle = sigmond_info.sigmond_info.ScatteringParticle(scattering_particle.name, 0, False)
        energy_obs = sigmond.MCObsInfo(f"{self.sh_name}/{at_rest_scattering_particle!r}", 0)
        single_particle = util.boost_obs(obs_handler, energy_obs, scattering_particle.psq, self.ensemble_spatial_extent)
        single_particles.append(single_particle)
        single_particle_names.append(scattering_particle.name)
        single_particle_moms.append(scattering_particle.psq)

      total = util.linear_superposition_obs(obs_handler, single_particles, [1.0]*len(single_particles))
      total_cm = util.boost_obs_to_cm(obs_handler, total, fit_info.operator.psq, self.ensemble_spatial_extent)

      if self.reference_fit_info is not None:
        ref_obs = sigmond.MCObsInfo(f"{self.sh_name}/{self.ref_name}", 0)
        total_cm = util.ratio_obs(obs_handler, total_cm, ref_obs)

      non_interacting_levels.append(obs_handler.getEstimate(total_cm))
      non_interacting_level_names.append([single_particle_names,single_particle_moms])

    return non_interacting_levels,non_interacting_level_names

  def _makePlot(self):

    thresholds = self._getThresholds()

    # get energies and non_interacting energies
    energies = dict()
    non_interacting_energies = dict()
    non_interacting_energy_names = dict()
    
    for operator_set, fit_infos in self.spectrum.items():
      if operator_set.is_rotated:
        channel = operator_set.channel
        irrep = self.latex_map.get(channel.irrep, channel.irrep)
        psq = channel.psq
        label = rf"${irrep} ({psq})$"
        energies[label] = list()
        for energy in self.energies[operator_set].values():
          energies[label].append(energy)

      else:
        for level, fit_info in enumerate(fit_infos):
          operator = fit_info.operator
          channel = operator.channel
          irrep = self.latex_map.get(channel.irrep, channel.irrep)
          psq = channel.psq
          label = rf"${irrep} ({psq})$"
          energy = self.energies[operator_set][level]
          if label not in energies:
            energies[label] = list()

          energies[label].append(energy)

      non_interacting_energies[label],non_interacting_energy_names[label] = self._getNonInteractingLevels(operator_set)

    plot_directory = os.path.join(self.results_dir, "spectrum_plot")
    os.makedirs(plot_directory, exist_ok=True)
    plot_filename = os.path.join(plot_directory, "spectrum")
    relative_plot_directory = os.path.join("spectrum_plot","spectrum")

    if self.non_interacting_energy_labels:
      utils.plotting.spectrum(thresholds, energies, non_interacting_energies, self.latex_map,
                            self.rotate_labels, self.plot_width_factor, self.ref_name, plot_filename,
                            self.latex_compiler, non_interacting_energy_names)
    else:
      utils.plotting.spectrum(thresholds, energies, non_interacting_energies, self.latex_map,
                            self.rotate_labels, self.plot_width_factor, self.ref_name, plot_filename,
                            self.latex_compiler)
    return relative_plot_directory


