import logging
from typing import NamedTuple
from sortedcontainers import SortedSet
import os
import pylatex
import itertools

import sigmond

import tasks.task
import utils.util as util
import sigmond_info.fit_info
import sigmond_info.sigmond_input
import sigmond_info.sigmond_log
import operator_info.operator
import operator_info.operator_set
import operator_info.channel


class TRange(NamedTuple):
  tmin: int
  tmax: int

class Fits(NamedTuple):
  name: str
  model: sigmond_info.fit_info.FitModel
  tranges: tuple
  ratio: bool
  exclude_times: tuple
  noise_cutoff: float


class FitCorrelators(tasks.task.Task):
  """ Fit Correlators Task

  TODO: Decide on better way to do tmin plots?
  """

  task_type = "fit_corrs"

  def initialize(self, operator_fits, **extra_options):
    """ initializes the needed objects to perform fits to correlators
        with sigmond

    Args:
      operator_fits ({OperatorSet: {Operator: {Fits: [FitInfo]}}}): 
          Specifies all the fits to be done.
      **minimizer_info (MinimizerInfo): Speicifies the 
          minimizer to use and the info to pass to the minimizer
      **fit_plots (bool): whether fit plots should be made
      **fit_plot_info (FitPlotInfo): Contains information for how to make
          the fit plots.
      **tmin_plots (bool): whether tmin plots should be made
      **tmin_plot_info (TMinPlotInfo): Contains information for how to make
          the tmin plots.
      **scattering_particles ({ScatteringParticle: Operator}): A dict
          of all the scattering particles used for ratio fits.
      **reference_fit_info (FitInfo): The fit to be used as a
          reference energy.
    """
    self.operator_fits = operator_fits
    self.minimizer_info = extra_options.pop('minimizer_info', sigmond.MinimizerInfo())
    self.fit_plots = extra_options.pop('fit_plots', True)
    self.fit_plot_info = extra_options.pop('fit_plot_info', sigmond_info.sigmond_info.FitPlotInfo())
    self.tmin_plots = extra_options.pop('tmin_plots', True)
    self.tmin_plot_info = extra_options.pop('tmin_plot_info', sigmond_info.sigmond_info.TMinPlotInfo())
    self.scattering_particle = extra_options.pop('scattering_particles', dict())
    self.reference_fit_info = extra_options.pop('reference_fit_info', None)

    util.check_extra_keys(extra_options, 'FitCorrelators.initialize')

  def readConfig(self, **task_options):
    """readConfig function for FitCorrelators

    YAML:
      subtractvev: true   # optional

      # optional
      minimizer_info:
        minimizer: lmder
        parameter_rel_tol: 1e-6
        chisquare_rel_tol: 1e-4
        max_iterations: 1024
        verbosity: high

      # optional
      fit_plots: true
      fit_plot_info:
        timestep: 3
        show_approach: true
        goodness: chisq
        corrname: standard
        symbol_color: blue
        symbol_type: circle
        max_relative_error: 0.0

      # optional
      tmin_plots: true
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
        operator: isotriplet S=0 P=(0,0,0) A1um P 0
        model: 1-exp
        tmin: 14
        tmax: 40

      # optional (used for ratio fits)
      scattering_particles:  
        - name: pi
          operators:
            - isotriplet S=0 P=(0,0,0) A1um P 0
            - isotriplet S=0 PSQ=1 A2m P 0
            - isotriplet S=0 PSQ=2 A2m P 0
            - isotriplet S=0 PSQ=3 A2m P 0
            ...
        - ...
        ...

      fits:
        - name: fit_type1
          model: 2-exp
          tmin: 7-12
          tmax: 15-30
          ratio: false             # optional
          exclude_times: [10,12]   # optional
          noise_cutoff: 1.2        # optional
        - ...
        ...

      ### operator_sets
      operator_sets:
        - name: pion_0
          operators: 
            - isotriplet S=0 P=(0,0,0) A1um P 0
          fits:                 # optional (if missing, does all fits)
            - fit_type1
        - name: op_basis1
          operators:
            - isodoublet S=0 P=(0,0,0) A1um P 0
            - isodoublet S=0 P=(0,0,0) A1um P 1
          non_interacting_levels:    # required for ratio fits
            - [pi(0), pi(2), pi(3)]
            - [pi(1), pi(2), pi(3)]
            - ...
          fits:                 # optional (if missing, does all fits)
            - fit_type1
            - ...
            ...
        - name: rotated_op_basis
          pivot_info:
            pivot_type: single_pivot
            norm_time: 5
            metric_time: 5
            diagonalize_time: 10
            max_condition_number: 100
        - ...
        ...
    """

    task_options['minimizer_info'] = sigmond_info.sigmond_info.getMinimizerInfo(task_options)
    task_options['fit_plot_info'] = sigmond_info.sigmond_info.FitPlotInfo.createFromConfig(task_options)
    task_options['tmin_plot_info'] = sigmond_info.sigmond_info.TMinPlotInfo.createFromConfig(task_options)

    ref_fit_info = task_options.pop('reference_fit_info', None)
    if ref_fit_info is not None:
      ref_fit_info = sigmond_info.fit_info.FitInfo.createFromConfig(ref_fit_info)
    task_options['reference_fit_info'] = ref_fit_info

    # check for scattering_particles
    scattering_particles = dict()
    for scattering_particle in task_options.pop('scattering_particles', []):
      try:
        name = scattering_particle.pop('name')
        op_strs = scattering_particle.pop('operators')
      except KeyError as err:
        logging.error(f"Missing required key in 'scattering_particles': {err}")

      operators = dict()
      for op_str in op_strs:
        operator = operator_info.operator.Operator(op_str)
        psq = operator.psq
        scattering_particle = sigmond_info.sigmond_info.ScatteringParticle(name, psq)
        if scattering_particle in scattering_particles:
          logging.error(f"Scattering particle '{scattering_particle}' encountered twice")

        scattering_particles[scattering_particle] = operator

    task_options['scattering_particles'] = scattering_particles

    # Read fits
    fits = dict()
    try:
      for fit in task_options.pop('fits'):
        name = fit.pop('name')
        if name in fits:
          logging.warning(f"Fit with name '{name}' appeared more than once...overwriting")

        model = sigmond_info.fit_info.FitModel(fit.pop('model'))
        exclude_times = fit.pop('exclude_times', tuple())
        noise_cutoff = fit.pop('noise_cutoff', 0.0)
        ratio = fit.pop('ratio', False)

        tmins = fit.pop('tmin')
        if isinstance(tmins, int):
          tmins = [int(tmins)]
        else:
          tmins = list(map(int, tmins.split('-')))
          tmins = list(range(tmins[0], tmins[-1]+1))

        tmaxs = fit.pop('tmax')
        if isinstance(tmaxs, int):
          tmaxs = [int(tmaxs)]
        else:
          tmaxs = list(map(int, tmaxs.split('-')))
          tmaxs = list(range(tmaxs[0], tmaxs[-1]+1))

        tranges = list()

        n_params = len(sigmond_info.fit_info.FitInfo.PARAMETERS[model])
        for tmin in tmins:
          for tmax in tmaxs:
            n_dof = tmax - tmin + 1 - n_params
            if n_dof < 1:
              continue

            tranges.append(TRange(tmin, tmax))

        util.check_extra_keys(fit, "fits")

        fits[name] = Fits(name, model, tuple(tranges), ratio, tuple(exclude_times), noise_cutoff)

    except KeyError as err:
      err = str(err)
      if err == "'fits'":
        logging.error("Fit tasks need a 'fits' section")
      elif err == "'name'":
        logging.error(f"No 'name' for fit in 'fits' section of task '{self.task_name}'")
      else:
        logging.error(f"Missing required key in '{name}' fit of task '{self.task_name}': {err}")

    # read operator sets
    subtractvev = task_options.pop('subtractvev', True)
    operator_fits = dict()
    non_interacting_operators_lists = dict()
    try:
      for operator_set_conf in task_options.pop('operator_sets'):
        operator_set = operator_info.operator_set.getOperatorSet(operator_set_conf)
        operator_fits[operator_set] = dict()
        operators = operator_set.getRotatedOperators() if operator_set.is_rotated else operator_set.operators
        non_interacting_levels = operator_set_conf.pop('non_interacting_levels', [None]*len(operators))

        fit_keys = operator_set_conf.pop('fits', list(fits.keys()))
        for operator, non_interacting_level in zip(operators, non_interacting_levels):
          if non_interacting_level is None:
            non_interacting_operators = None
          else:
            non_interacting_operators = sigmond_info.sigmond_info.NonInteractingOperators.create(
                scattering_particles, non_interacting_level)

          operator_fits[operator_set][operator] = dict()
          for fit_key in fit_keys:
            if fit_key not in fits:
              logging.error(f"Invalid fit name '{fit_key}")

            fits_to_do = fits[fit_key]
            operator_fits[operator_set][operator][fits_to_do] = {'normal': [], 'tmin': []}
            for trange in fits_to_do.tranges:
              fit_info = sigmond_info.fit_info.FitInfo(
                  operator, fits_to_do.model, trange.tmin, trange.tmax, subtractvev, fits_to_do.ratio,
                  list(fits_to_do.exclude_times), fits_to_do.noise_cutoff, non_interacting_operators)

              operator_fits[operator_set][operator][fits_to_do]['normal'].append(fit_info)

            for tmax, tranges in itertools.groupby(sorted(fits_to_do.tranges, key=lambda trange: trange.tmax), lambda trange: trange.tmax):
              tmins = list(trange.tmin for trange in tranges)
              tmin_min = min(tmins)
              tmin_max = max(tmins)
              fit_info = sigmond_info.fit_info.FitInfo(
                  operator, fits_to_do.model, tmin_min, tmax, subtractvev, fits_to_do.ratio,
                  list(fits_to_do.exclude_times), fits_to_do.noise_cutoff, non_interacting_operators,
                  tmin_max)

              operator_fits[operator_set][operator][fits_to_do]['tmin'].append(fit_info)

        util.check_extra_keys(operator_set_conf, "operator_sets")
    except KeyError as err:
      logging.error(f"Missing required key in 'operator_bases': {err}")

    self.initialize(operator_fits, **task_options)

  def project_name(self, operator_set, operator, fit_name):
    return super().project_name(f"{operator_set!r}_{operator!r}_{util.str_to_file(fit_name)}")

  def tmin_project_name(self, operator_set, operator, fit_name):
    return super().project_name(f"tmin_{operator_set!r}_{operator!r}_{util.str_to_file(fit_name)}")

  def logfile(self, operator_set, operator, fit_name):
    logdir = os.path.join(self.logdir, repr(operator_set), 'normal')
    os.makedirs(logdir, exist_ok=True)
    logfile = f"{self.task_name}_{operator!r}_{util.str_to_file(fit_name)}.log"
    return os.path.join(logdir, logfile)

  def tmin_logfile(self, operator_set, operator, fit_name):
    logdir = os.path.join(self.logdir, repr(operator_set), 'tmin')
    os.makedirs(logdir, exist_ok=True)
    logfile = f"{self.task_name}_{operator!r}_{util.str_to_file(fit_name)}.log"
    return os.path.join(logdir, logfile)

  def inputfile(self, operator_set, operator, fit_name):
    inputdir = os.path.join(self.inputdir, repr(operator_set), 'normal')
    os.makedirs(inputdir, exist_ok=True)
    inputfile = f"{self.task_name}_{operator!r}_{util.str_to_file(fit_name)}.input"
    return os.path.join(inputdir, inputfile)

  def tmin_inputfile(self, operator_set, operator, fit_name):
    inputdir = os.path.join(self.inputdir, repr(operator_set), 'tmin')
    os.makedirs(inputdir, exist_ok=True)
    inputfile = f"{self.task_name}_{operator!r}_{util.str_to_file(fit_name)}.input"
    return os.path.join(inputdir, inputfile)

  def fit_plotdir(self, operator_set, operator, fit_name):
    plotdir = os.path.join(self.plotdir, repr(operator_set), 'normal', repr(operator), util.str_to_file(fit_name))
    os.makedirs(plotdir, exist_ok=True)
    return plotdir

  def tmin_fit_plotdir(self, operator_set, operator, fit_name):
    plotdir = os.path.join(self.plotdir, repr(operator_set), 'tmin', repr(operator), util.str_to_file(fit_name))
    os.makedirs(plotdir, exist_ok=True)
    return plotdir

  def fit_plotfile(self, operator_set, fit_name, fit_info, extension):
    operator = fit_info.operator
    plotdir = self.fit_plotdir(operator_set, operator, fit_name)
    plotfile = f"{fit_info.plotfile}.{extension.value}"
    return os.path.join(plotdir, plotfile)

  def tmin_fit_plotfile(self, operator_set, fit_name, fit_info, extension):
    operator = fit_info.operator
    plotdir = self.tmin_fit_plotdir(operator_set, operator, fit_name)
    plotfile = f"{fit_info.plotfile}.{extension.value}"
    return os.path.join(plotdir, plotfile)

  def getSigmondInputs(self):
    sigmond_inputs = list()

    default_data_files = self.data_handler.data_files
    # get reference observable data files
    reference_observable = None
    if self.reference_fit_info is not None:
      reference_observable = self.reference_fit_info.energy_observable
      default_data_files += self.data_handler.getChannelDataFiles(self.reference_fit_info.operator.channel)

    for scattering_operator in self.scattering_particle.values():
      default_data_files += self.data_handler.getChannelDataFiles(scattering_operator.channel)

    for operator_set, operator_fits in self.operator_fits.items():
      if operator_set.is_rotated:
        data_files = default_data_files + self.data_handler.getRotatedDataFiles(operator_set)

      for operator, fits in operator_fits.items():
        if not operator_set.is_rotated:
          data_files = default_data_files + self.data_handler.getChannelDataFiles(operator.channel)

        for fit, fit_infos in fits.items():
          # normal fits
          project_name = self.project_name(operator_set, operator, fit.name)
          inputfile = self.inputfile(operator_set, operator, fit.name)
          logfile = self.logfile(operator_set, operator, fit.name)

          sigmond_input = self.new_sigmond_input(project_name, inputfile, logfile, data_files)
          if self.reference_fit_info is not None:
            sigmond_input.doTemporalCorrelatorFit(self.reference_fit_info, minimizer_info=self.minimizer_info)

          for fit_info in fit_infos['normal']:
            plotfile = self.fit_plotfile(operator_set, fit.name, fit_info, util.PlotExtension.grace)
            if self.fit_plots:
              sigmond_input.doTemporalCorrelatorFit(
                  fit_info, minimizer_info=self.minimizer_info, plotfile=plotfile,
                  plot_info=self.fit_plot_info, reference_energy=reference_observable)
            else:
              sigmond_input.doTemporalCorrelatorFit(
                  fit_info, minimizer_info=self.minimizer_info)

          sigmond_input.write()
          sigmond_inputs.append(sigmond_input)

          # tmin fits
          if self.tmin_plots:
            project_name = self.tmin_project_name(operator_set, operator, fit.name)
            inputfile = self.tmin_inputfile(operator_set, operator, fit.name)
            logfile = self.tmin_logfile(operator_set, operator, fit.name)

            sigmond_input = self.new_sigmond_input(project_name, inputfile, logfile, data_files)
            for fit_info in fit_infos['tmin']:
              plotfile = self.tmin_fit_plotfile(operator_set, fit.name, fit_info, util.PlotExtension.grace)
              sigmond_input.doTemporalCorrelatorFit(
                  fit_info, minimizer_info=self.minimizer_info, plotfile=plotfile,
                  plot_info=self.tmin_plot_info)

            sigmond_input.write()
            sigmond_inputs.append(sigmond_input)

    return sigmond_inputs


  def finalize(self):
    doc = util.create_doc(f"Correlator Fits: {self.ensemble_name} - {self.task_name}")


    for operator_set, operator_fits in self.operator_fits.items():
      with doc.create(pylatex.Section(str(operator_set))):
        for operator, fits in operator_fits.items():
          with doc.create(pylatex.Subsection(str(operator))):
            for fit, fit_infos in fits.items():

              # normal fits
              logfile = self.logfile(operator_set, operator, fit.name)
              fit_log = sigmond_info.sigmond_log.FitLog(logfile)
              if self.fit_plots:
                plotdir = self.fit_plotdir(operator_set, operator, fit.name)
                util.dirGrace2pdf(plotdir)

              section_title = f"{fit.name} - Model: {fit.model.short_name}"
              with doc.create(pylatex.Subsubsection(section_title)):
                self._add_fits(doc, fit_log, fit.name, operator_set, fit.ratio, fit.model.has_gap, fit.model.has_const)

              # tmin fits
              if self.tmin_plots:
                plotdir = self.tmin_fit_plotdir(operator_set, operator, fit.name)
                util.dirGrace2pdf(plotdir)
                tmin_fit_infos = list()
                for fit_info in fit_infos['tmin']:
                  plotfile = self.tmin_fit_plotfile(operator_set, fit.name, fit_info, extension=util.PlotExtension.pdf)
                  if os.path.isfile(plotfile):
                    tmin_fit_infos.append(fit_info)

                if len(tmin_fit_infos) == 0:
                  continue

                tmin_fit_infos.sort(key=lambda fit_info: fit_info.tmax)
                section_title = f"$t_{{\\rm min}}$ plots - {fit.name} - Model: {fit.model.short_name}"
                with doc.create(pylatex.Subsubsection(pylatex.NoEscape(section_title))):
                  self._add_tmins(doc, tmin_fit_infos, fit.name, operator_set, fit.ratio)

    results_dir = self.results_dir
    os.makedirs(results_dir, exist_ok=True)
    filename = os.path.join(results_dir, self.task_name)
    util.compile_pdf(doc, filename)

  def _add_tmins(self, doc, fit_infos, fit_name, operator_set, ratio):
    for fit_info1, fit_info2 in zip(*[iter(fit_infos)]*2):
      fit_plot1 = self.tmin_fit_plotfile(operator_set, fit_name, fit_info1, extension=util.PlotExtension.pdf)
      fit_plot2 = self.tmin_fit_plotfile(operator_set, fit_name, fit_info2, extension=util.PlotExtension.pdf)

      with doc.create(pylatex.Figure(position='H')) as fig:
        with doc.create(pylatex.SubFigure(position='b', width=pylatex.NoEscape(r'0.5\linewidth'))) as fig1:
          caption = f"{fit_info1.model.short_name}, $t_{{\\rm max}} = {fit_info1.tmax}$"
          if fit_info1.ratio:
            caption += " Ratio"
          if fit_info1.exclude_times:
            caption += f" $t_{{\\rm exc}}$ = {fit_info1.exclude_times}"
          util.add_image(fig1, self.results_dir, fit_plot1, width="1.0", caption=caption)
        with doc.create(pylatex.SubFigure(position='b', width=pylatex.NoEscape(r'0.5\linewidth'))) as fig2:
          caption = f"{fit_info2.model.short_name}, $t_{{\\rm max}} = {fit_info2.tmax}$"
          if fit_info1.ratio:
            caption += " Ratio"
          if fit_info2.exclude_times:
            caption += f" $t_{{\\rm exc}}$ = {fit_info2.exclude_times}"
          util.add_image(fig2, self.results_dir, fit_plot2, width="1.0", caption=caption)

    if len(fit_infos) % 2 != 0:
      fit_info = list(fit_infos)[-1]
      fit_plot = self.tmin_fit_plotfile(operator_set, fit_name, fit_info, extension=util.PlotExtension.pdf)
      with doc.create(pylatex.Figure(position='H')) as fig:
        caption = f"{fit_info.model.short_name}, $t_{{\\rm max}} = {fit_info.tmax}$"
        if fit_info.ratio:
          caption += " Ratio"
        if fit_info.exclude_times:
          caption += f" $t_{{\\rm exc}}$ = {fit_info.exclude_times}"
        util.add_image(fig, self.results_dir, fit_plot, width="0.5", caption=caption)

    doc.append(pylatex.NoEscape(r"\newpage"))


  def _add_fits(self, doc, fit_log, fit_name, operator_set, ratio, gap=False, const=False):
    with doc.create(pylatex.Center()) as centered:
      if gap and const and util.ERR_PREC <= 4:
        cols = "X[c] X[c] X[4,c] X[4,c] X[4,c] X[4,c] X[2,c] X[2,c] X[2,c] X[c]"
      elif (gap or const) and util.ERR_PREC <= 4:
        cols = "X[c] X[c] X[4,c] X[4,c] X[4,c] X[2,c] X[2,c] X[2,c] X[c]"
      else:
        cols = "X[c] X[c] X[4,c] X[4,c] X[2,c] X[2,c] X[2,c] X[c]"
      with centered.create(pylatex.LongTabu(cols, to=r"\linewidth")) as fit_table:
        if ratio:
          energy_header = r"$a_t \Delta E_{\rm fit}$"
        else:
          energy_header = r"$a_t E_{\rm fit}$"
        header_row = [
            pylatex.NoEscape(r"$t_{\rm min}$"),
            pylatex.NoEscape(r"$t_{\rm max}$"),
            pylatex.NoEscape(energy_header),
            "A",
            pylatex.NoEscape(r"$\chi^2 / \text{dof}$"),
            pylatex.NoEscape(r"$p$-value"),
            pylatex.NoEscape(r"$t_{\rm exc}$"),
            pylatex.NoEscape(r"$\sigma_{\rm cut}$"),
        ]
        if const and util.ERR_PREC <= 4:
          header_row.insert(4, pylatex.NoEscape(r"const"))
        if gap and util.ERR_PREC <= 4:
          header_row.insert(4, pylatex.NoEscape(r"$a \Delta$"))

        fit_table.add_row(header_row, mapper=[pylatex.utils.bold])
        fit_table.add_hline()
        fit_table.end_table_header()
        for fit_info, fit_result in fit_log.fits.items():
          value_row = [
              fit_info.tmin,
              fit_info.tmax,
              fit_result.energy,
              fit_result.amplitude,
              round(fit_result.chisq, 2),
              round(fit_result.quality, 2),
              fit_info.exclude_times,
              round(fit_info.noise_cutoff, 1),
          ]
          if const and util.ERR_PREC <= 4:
            value_row.insert(4, fit_result.const)
          if gap and util.ERR_PREC <= 4:
            value_row.insert(4, fit_result.gap)
          fit_table.add_row(value_row)

    if self.fit_plots:
      for fit_info1, fit_info2 in zip(*[iter(fit_log.fits.keys())]*2):
        fit_plot1 = self.fit_plotfile(operator_set, fit_name, fit_info1, extension=util.PlotExtension.pdf)
        fit_plot2 = self.fit_plotfile(operator_set, fit_name, fit_info2, extension=util.PlotExtension.pdf)

        with doc.create(pylatex.Figure(position='H')) as fig:
          with doc.create(pylatex.SubFigure(position='b', width=pylatex.NoEscape(r'0.5\linewidth'))) as fig1:
            caption = f"{fit_info1.model.short_name}, $t_{{\\rm fit}} = {fit_info1.tmin}, {fit_info1.tmax}$"
            if fit_info1.ratio:
              caption += " Ratio"
            if fit_info1.exclude_times:
              caption += f" $t_{{\\rm exc}}$ = {fit_info1.exclude_times}"
            if fit_info1.noise_cutoff:
              caption += f" $\sigma_{{\\rm cut}} = {fit_info1.noise_cutoff}$"
            util.add_image(fig1, self.results_dir, fit_plot1, width="1.0", caption=caption)
          with doc.create(pylatex.SubFigure(position='b', width=pylatex.NoEscape(r'0.5\linewidth'))) as fig2:
            caption = f"{fit_info2.model.short_name}, $t_{{\\rm fit}} = {fit_info2.tmin}, {fit_info2.tmax}$"
            if fit_info1.ratio:
              caption += " Ratio"
            if fit_info2.exclude_times:
              caption += f" $t_{{\\rm exc}}$ = {fit_info2.exclude_times}"
            if fit_info2.noise_cutoff:
              caption += f" $\sigma_{{\\rm cut}} = {fit_info2.noise_cutoff}$"
            util.add_image(fig2, self.results_dir, fit_plot2, width="1.0", caption=caption)

      if len(fit_log.fits) % 2 != 0:
        fit_info = list(fit_log.fits.keys())[-1]
        fit_plot = self.fit_plotfile(operator_set, fit_name, fit_info, extension=util.PlotExtension.pdf)
        with doc.create(pylatex.Figure(position='H')) as fig:
          caption = f"{fit_info.model.short_name}, $t_{{\\rm fit}} = {fit_info.tmin}, {fit_info.tmax}$"
          if fit_info.ratio:
            caption += " Ratio"
          if fit_info.exclude_times:
            caption += f" $t_{{\\rm exc}}$ = {fit_info.exclude_times}"
          if fit_info.noise_cutoff:
            caption += f" $\sigma_{{\\rm cut}} = {fit_info.noise_cutoff}$"
          util.add_image(fig, self.results_dir, fit_plot, width="0.5", caption=caption)

    doc.append(pylatex.NoEscape(r"\newpage"))

