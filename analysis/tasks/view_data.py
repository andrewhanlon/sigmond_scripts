import os
from enum import Enum

from sortedcontainers import SortedSet
import logging
import pylatex
import xml.etree.ElementTree as ET

import utils.util as util
import tasks.task
import operator_info.operator_set
import sigmond_info.sigmond_input
import sigmond_info.sigmond_info

import sigmond


class ViewData(tasks.task.Task):

  task_type = "view_data"

  def initialize(self, **options):
    """sets the needed paramters for the view data task

    Args:
      **channels (SortedSet(Channel)): The set of channels to print
      **operator_sets ([NamedOperatorSet]): A list of named operator sets
          to print
      **excluded_operators ({Operator}): A set of operators to not consider
      **auto_add (bool): Specifies whether channels should automatically
          be added.
      **off_diagonal (bool): specifies whether the off diagonal
          correlators should be printed
      **hermitian (bool): specifies whether the Hermitian tag
          should be used
      **subtractvev (bool): specifies whether operators with non-zero VEV
          should have the VEV subtracted
      **plot_info (PlotInfo): info about the plot
      **write_operators (bool): determines whether the operators
          should be written to file
      **split_pdfs (bool): specifies if PDFs should be split by channel,
          in which case it will group by irrep and P^2
    """

    self.channels = options.pop('channels', SortedSet())
    self.operator_sets = options.pop('operator_sets', dict())
    self.excluded_operators = options.pop('excluded_operators', set())

    self.off_diagonal = options.pop('off_diagonal', False)
    self.hermitian = options.pop('hermitian', True)
    self.subtractvev = options.pop('subtractvev', True)
    self.write_operators = options.pop('write_operators', False)
    self.split_pdfs = options.pop('split_pdfs', False)

    self.plot_info = options.pop('plot_info', sigmond_info.sigmond_info.PlotInfo())

    if options.pop('auto_add', False):
      self.channels |= self.data_handler.raw_channels

    util.check_extra_keys(options, self.task_name)

  def readConfig(self, **task_options):
    """readConfig function for ViewData

    YAML:
      off_diagonal: false                  # optional
      hermitian: true                      # optional
      subtractvev: false                   # optional

      write_operators: false
      split_pdfs: true                 # optional

      auto_add: false                      # optional

      plot_info:
        corrname: standard                   # optional
        timestep: 3                          # optional
        rescale: 1.0                         # optional
        symbol_color: blue                   # optional
        symbol_type: circle                  # optional
        max_relative_error: 2.0              # optional
        effective_energy_type: time_forward  # optional


      channels:
        - isospin: quintet
          strangeness: 0
          momentum: [0,0,0]
          irrep: A1gp
          irreprow: 1
        - ...
        ...

      excluded_operators:
        - kbar S=1 PSQ=1 A2 k 0
        - pion S=0 P=(0,0,0) A1um P 0

      operator_sets:
        - name: pions
          operators:
            - isotriplet S=0 P=(0,0,0) A1um P 0
            - isotriplet S=0 PSQ=1 A2m P 0
            - isotriplet S=0 PSQ=2 A2m P 0
            - isotriplet S=0 PSQ=3 A2m P 0
            - isotriplet S=0 PSQ=4 A2m P 0
        - ...
        ...
    """

    channels = SortedSet()
    for channel in task_options.pop('channels', []):
      channels.add(operator_info.channel.Channel(**channel))

    task_options['channels'] = channels

    operator_sets = list()
    for operator_set in task_options.pop('operator_sets', []):
      operator_sets.append(operator_info.operator_set.getOperatorSet(operator_set))

    excluded_operators = set()
    for excluded_operator in task_options.pop('excluded_operators', []):
      excluded_operators.add(operator_info.operator.Operator(excluded_operator))

    task_options['excluded_operators'] = excluded_operators

    task_options['operator_sets'] = operator_sets

    task_options['plot_info'] = sigmond_info.sigmond_info.PlotInfo.createFromConfig(task_options)

    self.initialize(**task_options)


  def getSigmondInputs(self):

    sigmond_inputs = list()

    all_operator_set_ops = SortedSet()
    if self.operator_sets:
      all_operator_set_ops = SortedSet.union(*[op_set.operators for op_set in self.operator_sets])

    for channel in self.channels:
      operators = [op for op in self.data_handler.getChannelOperators(channel) if op not in all_operator_set_ops and op not in self.excluded_operators]
      if not operators:
        continue

      data_files = self.data_files + self.data_handler.getChannelDataFiles(channel)
      project_name = self.project_name(repr(channel))
      logfile = self.logfile(repr(channel))
      inputfile = self.inputfile(repr(channel))
      sigmond_input = self.new_sigmond_input(project_name, inputfile, logfile, data_files)

      self.insertSigmondPlotTasks(sigmond_input, repr(channel), operators)

      sigmond_input.write()
      sigmond_inputs.append(sigmond_input)

    for operator_set in self.operator_sets:
      operators = operator_set.operators
      data_files = self.data_files
      for channel in operator_set.channels:
        channel_data_files = self.data_handler.getChannelDataFiles(channel)
        data_files += self.data_handler.getChannelDataFiles(channel)

      project_name = self.project_name(operator_set.name)
      logfile = self.logfile(operator_set.name)
      inputfile = self.inputfile(operator_set.name)
      sigmond_input = self.new_sigmond_input(project_name, inputfile, logfile, data_files)

      self.insertSigmondPlotTasks(sigmond_input, operator_set.name, operators)

      sigmond_input.write()
      sigmond_inputs.append(sigmond_input)

    return sigmond_inputs

  def op_file(self, opset_name):
    filename = f"{opset_name}.ops"
    dirname = os.path.join(self.results_dir, "operators")
    os.makedirs(dirname, exist_ok=True)
    return os.path.join(dirname, filename)

  @property
  def op_yaml_file(self):
    filename = os.path.join(self.results_dir, f"{util.str_to_file(self.task_name)}.yml")
    return filename

  def finalize(self):
    if self.write_operators and os.path.isfile(self.op_yaml_file):
      os.remove(self.op_yaml_file)

    doc = util.create_doc(
        f"Correlators and Effective Energies: {self.ensemble_name} - {self.task_name}")

    all_operator_set_ops = SortedSet()
    if self.operator_sets:
      all_operator_set_ops = SortedSet.union(*[op_set.operators for op_set in self.operator_sets])

    for channel in self.channels:
      data_files = self.data_files + self.data_handler.getChannelDataFiles(channel)

      operators = [op for op in self.data_handler.getChannelOperators(channel) if op not in all_operator_set_ops and op not in self.excluded_operators]
      if not operators:
        continue

      if self.write_operators:
        operator_info.operator_set.write_operators(self.op_file(repr(channel)), operators, True, False)
        operator_info.operator_set.write_operators_to_yaml(self.op_yaml_file, repr(channel), operators, True)


      with doc.create(pylatex.Section(str(channel))):
        self.addPlotsToPDF(doc, data_files, operators, repr(channel))

    for operator_set in self.operator_sets:
      if self.write_operators:
        operator_info.operator_set.write_operators(self.op_file(operator_set.name), operator_set.operators, True, False)
        operator_info.operator_set.write_operators_to_yaml(self.op_yaml_file, operator_set.name, operator_set.operators, True)

      data_files = self.data_files
      for channel in operator_set.channels:
        channel_data_files = self.data_handler.getChannelDataFiles(channel)
        data_files += self.data_handler.getChannelDataFiles(channel)

      with doc.create(pylatex.Section(operator_set.name)):
        self.addPlotsToPDF(doc, data_files, operator_set.operators, operator_set.name)

    filename = os.path.join(self.results_dir, f"{util.str_to_file(self.task_name)}_rebin{self.rebin}")
    util.compile_pdf(doc, filename)

  def finalize(self):
    if self.write_operators and os.path.isfile(self.op_yaml_file):
      os.remove(self.op_yaml_file)

    all_operator_set_ops = SortedSet()
    if self.operator_sets:
      all_operator_set_ops = SortedSet.union(*[op_set.operators for op_set in self.operator_sets])

    # create docs
    if self.split_pdfs:
      docs = dict()
      for channel in self.channels:
        operators = [op for op in self.data_handler.getChannelOperators(channel) if op not in all_operator_set_ops and op not in self.excluded_operators]
        if not operators:
          continue
        
        if channel.irrep_psq_key not in docs:
          docs[channel.irrep_psq_key] = util.create_doc(f"Correlators and Effective Energies: {self.ensemble_name} - {self.task_name} - {channel.irrep_psq_key}")

      for operator_set in self.operator_sets:
        docs[operator_set.name] = util.create_doc(f"Correlators and Effective Energies: {self.ensemble_name} - {self.task_name} - {operator_set.name}")

    else:
      doc = util.create_doc(f"Correlators and Effective Energies: {self.ensemble_name} - {self.task_name}")

    # create content
    for channel in self.channels:
      if self.split_pdfs:
        doc = docs[channel.irrep_psq_key]

      data_files = self.data_files + self.data_handler.getChannelDataFiles(channel)

      operators = [op for op in self.data_handler.getChannelOperators(channel) if op not in all_operator_set_ops and op not in self.excluded_operators]
      if not operators:
        continue

      if self.write_operators:
        operator_info.operator_set.write_operators(self.op_file(repr(channel)), operators, True, False)
        operator_info.operator_set.write_operators_to_yaml(self.op_yaml_file, repr(channel), operators, True)

      with doc.create(pylatex.Section(str(channel))):
        self.addPlotsToPDF(doc, data_files, operators, repr(channel))

    for operator_set in self.operator_sets:
      if self.write_operators:
        operator_info.operator_set.write_operators(self.op_file(operator_set.name), operator_set.operators, True, False)
        operator_info.operator_set.write_operators_to_yaml(self.op_yaml_file, operator_set.name, operator_set.operators, True)

      data_files = self.data_files
      for channel in operator_set.channels:
        channel_data_files = self.data_handler.getChannelDataFiles(channel)
        data_files += self.data_handler.getChannelDataFiles(channel)

      if self.split_pdfs:
        self.addPlotsToPDF(doc, data_files, operator_set.operators, operator_set.name)

      else:
        with doc.create(pylatex.Section(operator_set.name)):
          self.addPlotsToPDF(doc, data_files, operator_set.operators, operator_set.name)

    # compile
    if self.split_pdfs:
      for split_key, doc in docs.items():
        filename = os.path.join(self.results_dir, f"{util.str_to_file(self.task_name)}_{split_key}_rebin{self.rebin}")
        util.compile_pdf(doc, filename)
    else:
      filename = os.path.join(self.results_dir, f"{util.str_to_file(self.task_name)}_rebin{self.rebin}")
      util.compile_pdf(doc, filename)
