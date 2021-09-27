import os
import logging
from sortedcontainers import SortedSet
import pylatex

import sigmond

import tasks.task
import sigmond_info.sigmond_input
import operator_info.operator
import operator_info.operator_set
import utils.util as util
import utils.menu

class AverageCorrelators(tasks.task.Task):
  """ Average Correlators task

  TODO:
    Do a finalize in which averaged correlators are put into a PDF?
  """

  task_type = "average_corrs"

  def initiliaze(self, **options):
    """sets the needed parameters for the task

    Args:
      **averaged_channels (SortedSet(Channel)): The set of
          channels to be created after averaging over
          the available raw channels.
      **coefficients ({Operator: float}): the coefficients
          other than 1.0 to use for operators. All operators
          not present are assumed to have a coefficient
          of 1.0
      **excluded_operators ({Operator}): A set of operators to not consider
      **file_mode (sigmond.WriteMode): the file mode to use (see
          sigmond docs)
      **off_diagonal (bool): specifies whether the off diagonal
          correlators should be printed
      **hermitian (bool): whether we should be looking for hermitian
          correlators
      **use_spatial_info (bool): specifies if spatial info should be
          given in averaged op name
      **use_irrep_info (bool): specifies if irrep info should be
          given in averaged op name
      **plot_info (PlotInfo): info about the plot
      **write_operators (bool): determines whether the operators
          should be written to file
    """
    self.subtractvev = False

    self.hermitian = options.pop('hermitian', True)
    self.off_diagonal = options.pop('off_diagonal', False)
    self.plot_info = options.pop('plot_info', sigmond_info.sigmond_info.PlotInfo())
    self.write_operators = options.pop('write_operators', False)
    self.excluded_operators = options.pop('excluded_operators', set())

    self.use_spatial_info = options.pop('use_spatial_info', False)
    self.use_irrep_info = options.pop('use_irrep_info', False)

    raw_channels = self.data_handler.raw_channels
    self.averaged_channels = dict()
    for raw_channel in raw_channels:
      if raw_channel.is_averaged:
        logging.warning(f"Channel '{raw_channel}' is averaged already")
        continue

      averaged_channel = raw_channel.averaged
      if averaged_channel not in self.averaged_channels:
        self.averaged_channels[averaged_channel] = SortedSet()

      self.averaged_channels[averaged_channel].add(raw_channel)

    if (user_averaged_channels := options.pop('averaged_channels', SortedSet())):
      try:
        self.averaged_channels = {averaged_channel: self.averaged_channels[averaged_channel] for averaged_channel in user_averaged_channels}
      except KeyError as err:
        logging.error("Channel {err} not an averaged channel of the found raw channels")

    self.coefficients = options.pop('coefficients', dict())
    self.file_mode = options.pop('file_mode', sigmond.WriteMode.Overwrite)

    util.check_extra_keys(options, self.task_name)

  def readConfig(self, **task_options):
    """readConfig function for AverageCorrelators

    YAML:
      file_mode: overwrite  # optional

      off_diagonal: false
      hermitian: false

      # used for giving more details in averaged op name
      use_spatial_info: true
      use_irrep_info: true

      write_operators: false

      plot_info:
        corrname: standard                   # optional
        timestep: 3                          # optional
        rescale: 1.0                         # optional
        symbol_color: blue                   # optional
        symbol_type: circle                  # optional
        max_relative_error: 2.0              # optional
        effective_energy_type: time_forward  # optional

      coefficients:         # optional
        - operator: op_string
          coefficient: -92
        - ...
        ...

      # list of operators to exclude from averaging
      excluded_operators:
        - kbar S=1 P=(0,0,1) A2 k 0
        - ...

      averaged_channels:    # optional
        - isospin: quintet
          strangeness: 0
          momentum: [0,0,0]
          irrep: A1gp
        - ...
        ...
    """
    if 'averaged_channels' in task_options:
      averaged_channels = SortedSet()
      for averaged_channel in task_options.pop('averaged_channels'):
        averaged_channels.add(operator_info.channel.Channel(**averaged_channel))
      task_options['averaged_channels'] = averaged_channels

    if 'coefficients' in task_options:
      coefficients = dict()
      for coefficient in task_options.pop('coefficients'):
        try:
          operator = operator_info.operator.Operator(coefficient.get('operator'))
          coefficient = float(coefficient.get('coefficient'))
        except KeyError as err:
          logging.error(f"coefficients item missing required key '{err}'")

        coefficients[operator] = coefficient

      task_options['coefficients'] = coefficients

    if 'file_mode' in task_options:
      try:
        task_options['file_mode'] = sigmond.WriteMode.create(task_options['file_mode'])
      except ValueError:
        logging.error(f"Invalid file_mode in '{self.task_name}'")

    excluded_operators = set()
    for excluded_operator in task_options.pop('excluded_operators', []):
      excluded_operators.add(operator_info.operator.Operator(excluded_operator))

    task_options['excluded_operators'] = excluded_operators

    task_options['plot_info'] = sigmond_info.sigmond_info.PlotInfo.createFromConfig(task_options)

    self.initiliaze(**task_options)

  def datafile(self, project_dir, averaged_channel):
    datafile = f"corrs_{averaged_channel!r}.bin"
    return os.path.join(self.data_handler.averaged_datadir, datafile)

  def getSigmondInputs(self):
    self.final_info = dict()
    sigmond_inputs = list()
    for averaged_channel, raw_channels in self.averaged_channels.items():
      logging.info(f"Averaged channel {averaged_channel!s}")

      coefficients = dict()
      original_operators = list()
      result_operators = None
      last_ops_map = None

      data_files = self.data_files
      for raw_channel in raw_channels:
        logging.info(f"    Raw Channel {raw_channel!s}")
        data_files += self.data_handler.getChannelDataFiles(raw_channel)
        operators = [op for op in self.data_handler.getChannelOperators(raw_channel) if op not in self.excluded_operators]

        if len(operators) == 1:
          operator = list(operators)[0]
          if operator.psq == 9 and abs(operator.getXMomentum()) != 0 and abs(operator.getXMomentum()) != 3:
            continue

        for operator in operators:
          coeff = 1.
          if operator in self.coefficients:
            coefficients[operator] = self.coefficients[operator]
            coeff = coefficients[operator]

          logging.info(f"        {coeff} {operator!s}")

        ops_map = _getOperatorsMap(operators, averaged_channel, self.use_spatial_info, self.use_irrep_info)

        if result_operators is None:
          result_operators = sorted(list(ops_map.keys()))
        elif (check_result_operators := sorted(list(ops_map.keys()))) != result_operators:
          print("Mismatch between operator sets to average")
          print("Previous mappings:")
          choice_ops = list()
          for averaged_op, raw_op in sorted(last_ops_map.items()):
            if averaged_op in ops_map:
              print(f"  {averaged_op} <= {raw_op}")
            else:
              choice_ops.append(averaged_op)
              print(f"* {averaged_op} <= {raw_op}")

          print()
          print("Current mappings:")
          for averaged_op, raw_op in sorted(ops_map.items()):
            if averaged_op in last_ops_map:
              print(f"  {averaged_op} <= {raw_op}")
            else:
              print(f"* {averaged_op} <= {raw_op}")

          new_ops_map = dict()
          for averaged_op, raw_op in ops_map.items():
            if averaged_op in result_operators:
              new_ops_map[averaged_op] = raw_op
              continue

            if not choice_ops:
              break

            op_menu = utils.menu.Menu("Choice?", title=f"{averaged_op} <= {raw_op}")
            op_menu.addItems(*choice_ops)
            chosen_op = op_menu.getInput(multi=False, no_items=True)
            if chosen_op is None:
              continue

            choice_ops.remove(chosen_op)
            new_ops_map[chosen_op] = raw_op

          for not_chosen_op in choice_ops:
            result_operators.remove(not_chosen_op)

          ops_map = new_ops_map

        last_ops_map = ops_map
        original_operators.append(ops_map)

      for original_operator_sets in original_operators:
        for k in list(original_operator_sets.keys()):
          if k not in result_operators:
            del original_operator_sets[k]

      project_name = self.project_name(repr(averaged_channel))
      logfile = self.logfile(repr(averaged_channel))
      inputfile = self.inputfile(repr(averaged_channel))
      sigmond_input = self.new_sigmond_input(project_name, inputfile, logfile, data_files)

      datafile = self.datafile(self.project_dir, averaged_channel)
      min_time, max_time = self.data_handler.getChannelsLargestTRange(*raw_channels)
      logging.info(f"Time separations [{min_time},{max_time}]")

      sigmond_input.doCorrelatorMatrixSuperposition(
          result_operators, original_operators, min_time, max_time, 
          coefficients=coefficients, hermitian=self.hermitian, filename=datafile,
          file_type=sigmond_info.sigmond_info.DataFormat.bins, file_mode=self.file_mode)
      
      self.final_info[averaged_channel] = (result_operators, original_operators, coefficients)

      self.insertSigmondPlotTasks(sigmond_input, repr(averaged_channel), result_operators)

      sigmond_input.write()
      sigmond_inputs.append(sigmond_input)

    return sigmond_inputs

  def op_file(self, channel):
    filename = f"{channel!r}.ops"
    dirname = os.path.join(self.results_dir, "operators")
    os.makedirs(dirname, exist_ok=True)
    return os.path.join(dirname, filename)

  @property
  def op_yaml_file(self):
    filename = os.path.join(self.results_dir, f"{util.str_to_file(self.task_name)}.yml")
    return filename

  def _suggest_rotation_yml_file(self):
    yaml_file = os.path.join(self.results_dir, "rotate_suggestion.yml")
    file_mode = 'w+'
    f_handler = open(yaml_file, file_mode)
    
    f_handler.write(f"#fill/alter the necessary header information\n")
    f_handler.write("Execute:\n\n")
    proj_name = self.project_dir.split('/')[-2]
    f_handler.write(f"rotate_{proj_name}:\n")
    f_handler.write(f"  task_type: rotate_corrs\n  show_transformation: true\n  negative_eigenvalue_alarm: -0.10\n  subtractvev: false\n\n")
    f_handler.write(f"  plot_info:\n    corrname: standard\n    symbol_color: blue\n    symbol_type: circle\n    eff_energy_type: time_forward\n    timestep: 1\n\n")
    f_handler.write(f"  operator_bases:\n")
    
    for channel in self.averaged_channels:
      data_files = self.data_files + self.data_handler.getChannelDataFiles(channel)
      operators = self.data_handler.getChannelOperators(channel)
      if len(operators) > 1:
        f_handler.write(f"    - name: {repr(channel)}\n")
        f_handler.write(f"      pivot_info:\n")
        f_handler.write(f"        <<: *PIVOT_INFO\n")
        f_handler.write(f"      operators:\n")
        for operator in operators:
          f_handler.write(f"        - {operator.op_str()}\n")
        f_handler.write("\n")

    f_handler.close()
    logging.info(f"Suggested rotation yaml: {yaml_file}")
    return None

  def finalize(self):
    if self.write_operators and os.path.isfile(self.op_yaml_file):
      os.remove(self.op_yaml_file)  

    self.data_handler.findAveragedData()
    doc = util.create_doc(
        f"Averaged Correlators and Effective Energies: {self.ensemble_name} - {self.task_name}")

    for channel in self.averaged_channels:
      data_files = self.data_files + self.data_handler.getChannelDataFiles(channel)
      operators = self.data_handler.getChannelOperators(channel)
      if self.write_operators:
        operator_info.operator_set.write_operators(self.op_file(channel), operators, True, False)
        operator_info.operator_set.write_operators_to_yaml(self.op_yaml_file, repr(channel), operators, True)

      result_operators, original_operators, coefficients = self.final_info[channel]
      with doc.create(pylatex.Section(str(channel))):
        with doc.create(pylatex.Subsection("Operators averaged")):
          for result_operator in result_operators:
            with doc.create(pylatex.Itemize()) as itemize:
              itemize.add_item(result_operator.op_str())
              with doc.create(pylatex.Itemize()) as sub_itemize:
                for original_operator in original_operators:
                  the_op = original_operator[result_operator]
                  the_coeff = coefficients.get(the_op, 1.)
                  sub_itemize.add_item(f"{the_coeff} {the_op.op_str()}")

          doc.append(pylatex.NoEscape(r"\newpage"))

        self.addPlotsToPDF(doc, data_files, operators, repr(channel))

    self._suggest_rotation_yml_file()
    filename = os.path.join(self.results_dir, util.str_to_file(self.task_name))
    util.compile_pdf(doc, filename, self.latex_compiler)


def _getOperatorsMap(operators, averaged_channel, get_had_spat=False, get_had_irrep=False):
  op_map = dict()
  for operator in operators:
    averaged_op = _getAveragedOperator(operator, averaged_channel, get_had_spat, get_had_irrep)
    if averaged_op in op_map:
      logging.critical(f"Conflicting operators {operator} and {op_map[averaged_op]}")

    op_map[averaged_op] = operator

  return op_map

def _getAveragedOperator(operator, averaged_channel, get_had_spat=False, get_had_irrep=False):
  if operator.operator_type is sigmond.OpKind.GenIrrep:
    logging.critical("Averaging of GIOperators not currently supported")

  op_info = operator.operator_info.getBasicLapH()
  if op_info.getNumberOfHadrons() == 1:
    obs_name = f"{NAME_MAP[op_info.getFlavor()]}-{op_info.getHadronSpatialType(1)}_{op_info.getHadronSpatialIdNumber(1)}"
    obs_id = 0
  else:
    obs_name = ""
    for had_num in range(1, op_info.getNumberOfHadrons()+1):
      had_name = NAME_MAP[op_info.getHadronFlavor(had_num)]
      had_psq = op_info.getHadronXMomentum(had_num)**2 + op_info.getHadronYMomentum(had_num)**2 + op_info.getHadronZMomentum(had_num)**2
      had_str = str(had_psq)
      if get_had_spat:
        had_spat_type = op_info.getHadronSpatialType(had_num)
        had_spat_id = op_info.getHadronSpatialIdNumber(had_num)
        had_str += f"_{had_spat_type}{had_spat_id}"

      if get_had_irrep:
        had_irrep_str = op_info.getHadronLGIrrep(had_num)
        had_str += f"_{had_irrep_str}" 

      obs_name += f"{had_name}({had_str})"

    obs_id = op_info.getLGClebschGordonIdNum()

  return operator_info.operator.Operator(averaged_channel.getGIOperator(obs_name, obs_id))

# TODO: use flavor_map in operator_info/operator.py
NAME_MAP = {
    'pion': 'pi',
    'eta': 'e',
    'phi': 'p',
    'kaon': 'k',
    'kbar': 'kb',
    'nucleon': 'N',
    'delta': 'D',
    'sigma': 'S',
    'lambda': 'L',
    'xi': 'X',
    'omega': 'O',
}
