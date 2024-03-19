import os
import time

import tasks.task
import sigmond_info.sigmond_info
import sigmond_info.fit_info
import operator_info.operator
import utils.util as util

import sigmond


class Anisotropy(tasks.task.Task):

  task_type = "anisotropy"

  def initialize(self, scattering_particles, **options):
    """sets the needed paramters for the view data task

    Args:
      scattering_particles ({ScatteringParticle: FitInfo}): A dictionary
          of the scattering particles to their fit infos.
      **minimizer_info (sigmondbind.MinimizerInfo): Speicifies the 
          minimizer to use and the info to pass to the minimizer
      **anisotropy_plot_info (AnisotropyPlotInfo): Contains information for
        how to make the plots.
    """

    self.scattering_particles = scattering_particles
    self.minimizer_info = options.pop('minimizer_info', sigmond.MinimizerInfo())
    self.anisotropy_plot_info = options.pop('anisotropy_plot_info', sigmond_info.sigmond_info.AnisotropyPlotInfo())

    util.check_extra_keys(options, self.task_name)

  def readConfig(self, **task_options):
    """readConfig function for ViewData

    YAML:
      # optional
      minimizer_info:
        minimizer: lmder
        parameter_rel_tol: 1e-6
        chisquare_rel_tol: 1e-4
        max_iterations: 1024
        verbosity: high

      particles: [pi, kaon]   # optional (default to include all)

      # optional
      anisotropy_plot_info:
        goodness: chisq
        symbol_color: blue
        symbol_type: circle

      scattering_particles:  
        - name: pi
          fits:
            - operator: isotriplet S=0 P=(0,0,0) A1um P 0
              model: 1-exp
              tmin: 12
              tmax: 30
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
    """

    task_options['minimizer_info'] = sigmond_info.sigmond_info.getMinimizerInfo(task_options)
    task_options['anisotropy_plot_info'] = sigmond_info.sigmond_info.AnisotropyPlotInfo.createFromConfig(task_options)

    particles = task_options.pop('particles', list())
    
    # check for scattering_particles
    scattering_particles = dict()
    try:
      for scattering_particle in task_options.pop('scattering_particles'):
        name = scattering_particle.pop('name')
        if particles and name not in particles:
          continue

        scattering_particles[name] = list()
        fits = scattering_particle.pop('fits')
        for fit in fits:
          fit.pop('tmin_info', None) # removes tmin_info if someone includes it
          operator = operator_info.operator.Operator(fit.pop('operator'))
          fit_model = sigmond_info.fit_info.FitModel(fit.pop('model'))
          fit_info = sigmond_info.fit_info.FitInfo(operator, fit_model, **fit)
          scattering_particles[name].append(fit_info)

    except KeyError as err:
      logging.error(f"Missing required key in 'scattering_particles': {err}")

    self.initialize(scattering_particles, **task_options)

  def getSigmondInputs(self):
    sigmond_inputs = list()
    for particle_name, fit_infos in self.scattering_particles.items():
      data_files = self.data_handler.data_files
      for fit_info in fit_infos:
        data_files += self.data_handler.getChannelDataFiles(fit_info.operator.channel)

      project_name = self.project_name(particle_name)
      inputfile = self.inputfile(particle_name)
      logfile = self.logfile(particle_name)

      sigmond_input = self.new_sigmond_input(project_name, inputfile, logfile, data_files)

      scattering_particle_energies = list()
      for fit_info in fit_infos:
        sigmond_input.doTemporalCorrelatorFit(fit_info, minimizer_info=self.minimizer_info)
        obs_info = fit_info.energy_observable
        psq = fit_info.operator.psq
        scattering_particle_energies.append((obs_info,psq))

      pdf_plotfile = f"{self.task_name}_{particle_name}.pdf"
      pdf_plotfile = os.path.join(self.results_dir, pdf_plotfile)
      if os.path.isfile(pdf_plotfile):
        os.remove(pdf_plotfile)
      plotfile = f"{self.task_name}_{particle_name}.agr"
      plotfile = os.path.join(self.results_dir, plotfile)
      if os.path.isfile(plotfile):
        os.remove(plotfile)
      sigmond_input.doAnisotropyFromDispersion(
          scattering_particle_energies, self.ensemble_spatial_extent, name=particle_name,
          minimizer_info=self.minimizer_info, plotfile=plotfile,
          plot_info=self.anisotropy_plot_info)
      
      sigmond_input.write()
      sigmond_inputs.append(sigmond_input)

    return sigmond_inputs


  def finalize(self):
    time.sleep(1) # needed or else grace file might not be create yet...I guess.
    util.dirGrace2pdf(self.results_dir)

