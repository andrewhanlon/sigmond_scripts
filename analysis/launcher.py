import subprocess
import os
import logging
from typing import NamedTuple
import concurrent.futures
import multiprocessing
import xml.etree.ElementTree as ET
import datetime
import math, time

import yaml
import progressbar
from simple_slurm import Slurm

import tasks.view_data
import tasks.average_corrs
import tasks.rotate_correlators
import tasks.fit_correlators
import tasks.spectrum
import tasks.reconstruct_fit_ratio
import tasks.anisotropy
import tasks.dispersion
import data_handling.data_files
import utils.util as util
import sigmond

task_map = {
    tasks.view_data.ViewData.task_type                   : tasks.view_data.ViewData,
    tasks.average_corrs.AverageCorrelators.task_type     : tasks.average_corrs.AverageCorrelators,
    tasks.rotate_correlators.RotateCorrelators.task_type : tasks.rotate_correlators.RotateCorrelators,
    tasks.fit_correlators.FitCorrelators.task_type       : tasks.fit_correlators.FitCorrelators,
    tasks.spectrum.Spectrum.task_type                    : tasks.spectrum.Spectrum,
    tasks.reconstruct_fit_ratio.ReconstructFitRatio.task_type : tasks.reconstruct_fit_ratio.ReconstructFitRatio,
    tasks.anisotropy.Anisotropy.task_type                : tasks.anisotropy.Anisotropy,
    tasks.dispersion.Dispersion.task_type                : tasks.dispersion.Dispersion,
}

            
class ProjectInfo(NamedTuple):
  project_dir: str
  raw_data_dirs: list
  ensembles_file: str
  echo_xml: bool
  bins_info: sigmond.MCBinsInfo
  sampling_info: sigmond.MCSamplingInfo
  data_files: data_handling.data_files.DataFiles
  precompute: bool
  latex_compiler: str


def launch(config_filenames):
  executor, tasks_found = read_config(config_filenames)

  for task_name, task in tasks_found.items():
    logging.info(f"Starting task '{task_name}'")
    sigmond_inputs = task.getSigmondInputs()
    executor.execute(sigmond_inputs)
    task.finalize()


def read_config(config_filenames):

  # Get the config file
  if config_filenames is None:
    config_filenames = [create_config()]

  config = _read_config_files(config_filenames)
  yaml_file_contents = yaml.dump(config)

  # Construct the Executor
  executor = None
  if 'Execute' in config:
    try:
      executor = Executor(**config.pop('Execute'))
    except TypeError as err:
      logging.critical("Invalid arguments passed to 'Execute' block")

  # Read 'Initialize' section which includes
  # 'precompute', project directory and possibly some raw data directories
  try:
    init_conf = config.pop('Initialize')
    project_dir = init_conf['project_directory']
    logging.info(f"Project directory: {project_dir}")
  except KeyError as err:
    logging.critical(f"Missing required key {err}")

  ensembles_file = init_conf.get('ensembles_file', '')
  if ensembles_file:
    logging.info(f"Ensembles File: {ensembles_file}")
  else:
    logging.info(f"Ensembles File: default")

  precompute = init_conf.get('precompute_bootstraps', True)
  logging.info(f"Precompute Bootstraps: {precompute}")
  echo_xml = init_conf.get('echo_xml', False)
  logging.info(f"Echo XML: {echo_xml}")

  raw_data_dirs = init_conf.get('raw_data_directories', list())
  if raw_data_dirs:
    _raw_dir_str = "Raw data directories:"
    _check_is_not_parent_dir(project_dir, raw_data_dirs)
    for raw_data_dir in raw_data_dirs:
      _raw_dir_str += f"\n        {raw_data_dir}"

    logging.info(_raw_dir_str)
  else:
    logging.info("No raw data directories given")

  # check precision
  util.ERR_PREC = config.pop('precision', 2)
  logging.info(f"Error precision set to {util.ERR_PREC}")

  # get latex compiler
  latex_compiler = config.pop('latex_compiler', None)

  # Read in the MCBinsInfo (which includes ensemble info)
  try:
    bins_info_config = config.pop('MCBinsInfo')
    bins_info = _get_bins_info(bins_info_config, ensembles_file)

  except KeyError as err:
    logging.critical(f"Missing required key {err}")

  # Read in the MCSamplingInfo
  sampling_info = config.pop('MCSamplingInfo', dict())
  sampling_info = _get_sampling_info(sampling_info, bins_info)

  # Read in any custom MCObservables specified by the user
  if 'MCObservables' in config:
    obs_conf = config.pop('MCObservables')
    checksums = obs_conf.get('checksums', False)
    data_files = data_handling.data_files.DataFiles(checksums)
    if 'BLCorrelatorData' in obs_conf:
      file_lists = _get_file_list_infos(obs_conf.pop('BLCorrelatorData', []))
      data_files.addCorrelators(*file_lists)

    if 'BLVEVData' in obs_conf:
      file_lists = obs_conf['BLVEVData']
      file_lists = _get_file_list_infos(obs_conf.pop('BLVEVData', []))
      data_files.addVEVs(*file_lists)

    if 'BinData' in obs_conf:
      data_files.addBins(*obs_conf['BinData'])

    if 'SamplingData' in obs_conf:
      data_files.addSamplings(*obs_conf['SamplingData'])

    logging.info(f"Custom data files:\n{data_files}")

  else:
    data_files = data_handling.data_files.DataFiles()

  # Collect info for the Tasks
  project_info = ProjectInfo(
      project_dir=project_dir, raw_data_dirs=raw_data_dirs, ensembles_file=ensembles_file,
      echo_xml=echo_xml, bins_info=bins_info, sampling_info=sampling_info, data_files=data_files,
      precompute=precompute, latex_compiler=latex_compiler)

  if not config:
    logging.error("No tasks blocks found")

  # output used config to file
  yaml_config_dir = os.path.join(project_dir, 'yaml_confs')
  os.makedirs(yaml_config_dir, exist_ok=True)
  yaml_filename = f"input_{util.datestr()}.yml"
  yaml_file = os.path.join(yaml_config_dir, yaml_filename)
  with open(yaml_file, 'w') as yaml_file_h:
    yaml_file_h.write(yaml_file_contents)
    logging.info(f"Wrote YAML input file {yaml_filename}")

  tasks_found = dict()
  for task_name, task_options in config.items():
    try:
      task_type = task_options.pop("task_type")
    except KeyError:
      logging.error(f"Task '{task_name}' missing task_type")
    except TypeError:
      logging.error(f"Task '{task_name}' has invalid task_type: '{task_options[-1]}'")

    if task_type not in task_map:
      logging.error(f"Task type '{task_type}' not recognized")

    task_obj = task_map[task_type](project_info, task_name)
    task_obj.readConfig(**task_options)
    tasks_found[task_name] = task_obj

  return executor, tasks_found


class Executor:

  sigmond_return_string = "Error: batch mode requires a file name as the only argument\n"

  def __init__(self, mode, sigmond_batch, **exec_options):
    self.mode = mode
    self.sigmond_batch = sigmond_batch
    self.exec_options = exec_options
    
    #for slurm:
    self.email = None
    self.max_hours = 10
    self.partition = 'blue'
    self.max_mem_in_mb = 3000

    # check sigmond_batch
    try:
      result = subprocess.run(self.sigmond_batch, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      if result.stdout.decode('utf-8') != self.sigmond_return_string:
        logging.error(f"sigmond_batch is invalid '{self.sigmond_batch}'")

    except FileNotFoundError:
      logging.error(f"sigmond_batch not found at '{self.sigmond_batch}'")

    if self.mode == "local" or self.mode == "slurm":
      try: #if task=rotate and write_to_file=True, max_sim should be 1
        self.simultaneous_jobs = self.exec_options.get('max_simultaneous', multiprocessing.cpu_count())
      except ValueError as err:
        logging.error("Invalid value passed to 'max_simultaneous'")
      if self.mode == "slurm":
        self.email = self.exec_options.get('email', None)
        self.max_hours = self.exec_options.get('max_hours', 10)
        self.partition = self.exec_options.get('partition', 'blue') #figure out how to detect this?
        self.max_mem_in_mb = self.exec_options.get('max_mem_in_mb', 3000) #figure out how to detect this?

      logging.info("Local execution with {} max simultaneous jobs".format(self.simultaneous_jobs))
      logging.info("Sigmond batch binary: {}".format(self.sigmond_batch))
        

    elif self.mode == "PBS":
      logging.critical("PBS execution not yet supported...")
    else:
      logging.critical("Invalid execution mode '{}'".format(self.mode))

  def execute(self, sigmond_inputs):
    if len(sigmond_inputs) == 0:
      logging.info("No sigmond inputs passed to executor...")
      return

    logging.info(f"Executing {len(sigmond_inputs)} job(s)...")
    if self.mode == "local":
      executor = concurrent.futures.ThreadPoolExecutor(max_workers=self.simultaneous_jobs)
      futures = []
      for sig_input in sigmond_inputs:
        futures.append(executor.submit(subprocess.run, [self.sigmond_batch, sig_input.filename]))

      executor.shutdown(False)
      status = concurrent.futures.wait(futures, return_when=concurrent.futures.FIRST_COMPLETED,
                                       timeout=0.1)
      pbar_widgets = [progressbar.SimpleProgress(), ' job(s) finished ', progressbar.Percentage(),
                      ' ', progressbar.Bar()]
      pbar = progressbar.ProgressBar(widgets=pbar_widgets, maxval=len(futures)).start()
      while len(status.not_done):
        pbar.update(len(status.done))
        status = concurrent.futures.wait(futures, return_when=concurrent.futures.FIRST_COMPLETED)
      pbar.finish()
    
    if self.mode == "slurm":
      if len(sigmond_inputs) < self.simultaneous_jobs:
        self.simultaneous_jobs = len(sigmond_inputs)
        
      n_tasks_per_node = math.ceil( len(sigmond_inputs)/self.simultaneous_jobs )
        
      jobs = []
    
      groups = [sigmond_inputs[n:n+n_tasks_per_node] for n in range(0, len(sigmond_inputs), n_tasks_per_node)]
    
      for group in groups:
          slurm = Slurm(job_name='run_sigmond')
          slurm.add_arguments( partition=self.partition )
          slurm.add_arguments( time=datetime.timedelta(hours=self.max_hours) )
          slurm.add_arguments( nodes=1 ) 
          slurm.add_arguments( ntasks_per_node=1 ) 
          slurm.add_arguments( mem=f'{self.max_mem_in_mb}mb' ) 
          if group==groups[-1] and self.email:
              slurm.add_arguments( mail_type='END' ) 
              slurm.add_arguments( mail_user=self.email )
        
          command = ""
          for sig_input in group:
              command+=f"{self.sigmond_batch} {sig_input.filename};"
          jobs.append(slurm.sbatch(command))
            
      #wait for job to finish
      for job in jobs:
          wait_for_job_completion(job)
        

def wait_for_job_completion(job_id):
    slurm = Slurm()
    while True:
        output = subprocess.check_output(["squeue", "-j", str(job_id), "-h"], universal_newlines=True)
        lines = output.strip().split("\n")
        if '' in lines:
            lines.remove('')
        if not lines:
            print(f"Job {job_id} completed."+" "*50)
            break
        else:
            output = subprocess.check_output(["squeue", "-j", str(job_id), "-h"], universal_newlines=True)
            print(output.strip(), end="\r" ) 
            time.sleep(2)  # Adjust the sleep time as needed
            
def create_config():
  logging.critical("Interactive mode not currently supported")

def _read_config_files(config_filenames):
  """ read the config file """

  config_file_contents = ""
  for config_filename in config_filenames:
    logging.info(f"Reading config file '{config_filename}'")
    with open(config_filename) as config_file:
      config_file_contents += f"\n{config_file.read()}"

  config = yaml.safe_load(config_file_contents.strip())

  # remove all 'x-' keys
  config = {k: v for k, v in config.items() if not k.startswith('x-')}

  return config

def _get_bins_info(bins_info_config, ensembles_file):
  if 'num_measurements' in bins_info_config:
    ensemble_info = sigmond.MCEnsembleInfo(
        bins_info_config['ensemble_id'], bins_info_config['num_measurements'],
        bins_info_config['num_streams'], bins_info_config['n_x'],
        bins_info_config['n_y'], bins_info_config['n_z'], bins_info_config['n_t'])
  elif ensembles_file:
    ensemble_info = sigmond.MCEnsembleInfo(bins_info_config['ensemble_id'], ensembles_file)
  else:
    ensemble_info = sigmond.MCEnsembleInfo(bins_info_config['ensemble_id'])
    
  logging.info("Ensemble: {}".format(ensemble_info))

  if 'keep_first' in bins_info_config:
    new_bins_info = ET.Element('MCBinsInfo')
    new_bins_info.append(ensemble_info.xml())
    tweaks = ET.SubElement(new_bins_info,'TweakEnsemble')
    ET.SubElement(tweaks,'KeepFirst').text = str(bins_info_config['keep_first'])
    ET.SubElement(tweaks,'KeepLast').text = str(bins_info_config['keep_last'])
    bins_info = sigmond.MCBinsInfo(sigmond.XMLHandler().set_from_string(ET.tostring(new_bins_info)))
  else:
    bins_info = sigmond.MCBinsInfo(ensemble_info)
  bins_info.setRebin(bins_info_config.get('rebin', 1))
  bins_info.addOmissions(set(bins_info_config.get("omissions", [])))
  logging.info("Rebin: {}".format(bins_info.getRebinFactor()))
  logging.info("Config omissions: {}".format(bins_info.getOmissions()))

  return bins_info

def _get_sampling_info(sampling_info_config, bins_info):
  if 'mode' in sampling_info_config:
    try:
      sampling_mode = sigmond.SamplingMode.create(sampling_info_config['mode'].lower())
    except KeyError as err:
      logging.critical("Unknown sampling mode {}".format(err))

  else:
    sampling_mode = sigmond.SamplingMode.Jackknife

  logging.info(str(sampling_mode).replace(".", ": "))

  if sampling_mode == sigmond.SamplingMode.Bootstrap:
    try:
      sampling_info = sigmond.MCSamplingInfo(
          sampling_info_config['number_resampling'], sampling_info_config['seed'],
          sampling_info_config['boot_skip'])

      logging.info("Number Of Resamplings: {}".format(sampling_info.getNumberOfReSamplings(bins_info)))
      logging.info("Seed: {}".format(sampling_info.getRNGSeed()))
      logging.info("Boot skip: {}".format(sampling_info.getSkipValue()))

    except KeyError as err:
      logging.critical("Missing required key {}".format(err))
  else:
    sampling_info = sigmond.MCSamplingInfo()

  return sampling_info

def _get_file_list_infos(conf_infos):
  file_list_infos = list()
  for file_list_conf in conf_infos:
    try:
      file_stub = file_list_conf.pop('file_stub')
      min_suffix = file_list_conf.pop('min_suffix')
      max_suffix = file_list_conf.pop('max_suffix')
      overwrite = file_list_conf.pop('overwrite', False)
    except KeyError as err:
      logging.error(f"Missing required key in file list info {err}")

    util.check_extra_keys(file_list_conf, 'file_list_info')
    file_list_info = sigmond.FileListInfo(file_stub, min_suffix, max_suffix, overwrite)
    file_list_infos.append(file_list_info)

  return file_list_infos

def _check_is_not_parent_dir(project_dir, data_dirs):
  parent_path = os.path.abspath(project_dir)
  for data_dir in data_dirs:
    child_path = os.path.abspath(data_dir)

    if parent_path == os.path.commonpath([parent_path, child_path]):
      logging.critical(f"Data directory '{child_path}' cannot be a subdirectory of project directory '{parent_path}'")
