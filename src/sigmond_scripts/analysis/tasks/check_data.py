from sortedcontainers import SortedSet

from tasks.task import Task
from operator_info.channel import Channel


class CheckData(Task):
  """The check data task

  TODO:
    Have the finalize put results in a PDF?
  """

  task_type = "check_data"

  def initiliaze(self, **options):
    """sets the needed parameters for the task

    Args:
      **channels (SortedSet(Channel)): The set of channels
          to be checked. If not provided, all raw channels
          found will be used.
      **outlierscale (int): the outlier scale (see sigmond)
      **hermitian (bool): whether the hermitian option
          should be used for the correlator check
      **subtractvev (bool): whether the subtractvev flag
          should be used for operators with non-zero VEV
    """
    self.channels = options.pop('channels', self.data_handler.getRawChannels())
    self.outlierscale = options.pop('outlierscale', 10)
    self.hermitian = options.pop('hermitian', True)
    self.subtractvev = options.pop('subtractvev', True)

    util.check_extra_keys(options, self.task_name)

  def readConfig(self, **task_options):
    """readConfig function for CheckData

    YAML:
      outlierscale: 10  # optional
      hermitian: true   # optional
      subtractvev: true # optional

      channels:         # optional
        - isospin: quintet
          strangeness: 0
          momentum: [0,0,0]
          irrep: A1gp
          irreprow: 1
        - ...
        ...
    """
    if 'channels' in task_options:
      channels = SortedSet()
      for channel in task_options.pop('channels'):
        channels.add(Channel(**channel))
      task_options['channels'] = channels

    self.initiliaze(**task_options)

  def getSigmondInputs(self):

    sigmond_inputs = list()
    for channel in self.channels:
      logging.info(f"{self.task_name} sigmond input for channel {channel.nicename}")
      operator_set = self.data_handler.getRawOperatorSet(channel)
      min_time, max_time = self.data_handler.getMaximumTRange(operator_set)
      data_files = self.data_handler.getCorrelatorDataFiles(operator_set)
      project_name = f"{self.task_type}_{self.task_name}_{channel.filename}"
      logfile = self.logfile(channel.filename)
      inputfile = self.inputfile(channel.filename)

      sigmond_input = self.new_sigmond_input(project_name, inputfile, logfile, data_files)
      sigmond_input.doHermitianCheck(operator_set.operators, min_time, max_time, verbose=True)

      subtractvev = self.subtractvev and channel.vev
      sigmond_input.doCorrelatorCheck(
          operator_set.operators, min_time, max_time, hermitian=self.hermitian,
          subtractvev=subtractvev, outlier=self.outlierscale, verbose=True)

      sigmond_input.write()
      sigmond_inputs.append(sigmond_input)

    return sigmond_inputs

  def finalize(self):
    for channel in self.channels:
      logfile = self.logfile(channel.filename)
      result_filename = f"checkdata.{channel.filename}.log"
      result_filename = os.path.join(self.results_dir(), result_filename)
      shutil.copyfile(logfile, result_filename)
      logging.info(f"Log file copied to: {result_filename}")
