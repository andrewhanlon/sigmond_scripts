import os
import stat
import logging
from sortedcontainers import SortedSet
from typing import NamedTuple

import utils.util as util
import operator_info.operator
import operator_info.channel
import sigmond_info.sigmond_input
import data_handling.data_handler

import sigmond

def getOperatorSet(options):
  operators = [operator_info.operator.Operator(op) for op in options.pop('operators', list())]
  if (pivot_info := options.pop('pivot_info', False)):
    try:
      name = options.pop('name')
      pivot_info = sigmond_info.sigmond_info.PivotInfo(**pivot_info)
      operator_set = operator_info.operator_set.RotatedOperatorSet(name, pivot_info, *operators)
    except KeyError as err:
      logging.error("RotatedOperatorSet missing key {err}")

  elif (name := options.pop('name', False)):
    operator_set = operator_info.operator_set.NamedOperatorSet(name, *operators)
  else:
    operator_set = operator_info.operator_set.OperatorSet(*operators)

  return operator_set


class OperatorSet:

  def __init__(self, *operators):
    """OperatorBasis __init__ method

    Args:
      *operators (list of Operators): a list of operators to use.
    """

    self._operators = dict()
    for operator in operators:
      self._add_operator(operator)


  def addOperator(self, operator):
    self._add_operator(operator)

  def _add_operator(self, operator):
    if operator.channel not in self._operators:
      self._operators[operator.channel] = SortedSet()

    self._operators[operator.channel].add(operator)

  '''
  def items(self):
    return self._operators.items()
  '''

  @property
  def is_rotated(self):
    return hasattr(self, "_pivot_info")

  @property
  def channel(self):
    if self.num_channels == 1:
      return self.channels.pop()

    logging.critical("OperatorSet::channel requires only one channel")

  @property
  def channels(self):
    return SortedSet(self._operators.keys())

  @property
  def operators(self):
    return SortedSet().union(*self._operators.values())

  def getOperators(self, *channels):
    operator_set = SortedSet()
    if channels:
      for channel in channels:
        operator_set.update(self._operators[channel])
    else:
      operator_set.update(self.operators)
    
    return operator_set

  @property
  def num_operators(self):
    return len(self.operators)

  @property
  def num_channels(self):
    return len(self.channels)

  def __str__(self):
    return '\n'.join([str(op) for op in self.operators])

  def __repr__(self):
    return "[" + ','.join([repr(op) for op in self.operators]) + "]"

  def __eq__(self, other):
    return self.operators == other.operators

  def __ne__(self, other):
    return not self.__eq__(other)

  def __lt__(self, other):
    return tuple(self.operators) < tuple(other.operators)

  def __le__(self, other):
    return tuple(self.operators) <= tuple(other.operators)

  def __gt__(self, other):
    return tuple(self.operators) > tuple(other.operators)

  def __ge__(self, other):
    return tuple(self.operators) >= tuple(other.operators)


class NamedOperatorSet(OperatorSet):

  def __init__(self, name, *operators):
    self._name = name
    super().__init__(*operators)

  def addOperator(self, operator):
    logging.error("NamedOperatorSet is not mutable")

  @property
  def name(self):
    return self._name

  def __str__(self):
    return self.name

  def __repr__(self):
    return util.str_to_file(self.name)

  def __cmp(self):
    return (self.name, tuple(self.operators))

  def __hash__(self):
    return hash(self.__cmp())

  def __eq__(self, other):
    return self.__cmp() == other.__cmp()

  def __ne__(self, other):
    return not self.__eq__(other)

  def __lt__(self, other):
    return self.__cmp() < other.__cmp()

  def __le__(self, other):
    return self.__cmp() <= other.__cmp()

  def __gt__(self, other):
    return self.__cmp() > other.__cmp()

  def __ge__(self, other):
    return self.__cmp() >= other.__cmp()


class RotatedOperatorSet(NamedOperatorSet):

  def __init__(self, name, pivot_info, *operators):
    self._name = name
    self._pivot_info = pivot_info

    operators = self._verifyOperators(operators)

    super().__init__(name, *operators)

  def _verifyOperators(self, operators):
    operators = SortedSet(operators)
    pivot_file = data_handling.data_handler.DataHandler().pivotfile(self)
    if (pivot_exists := os.path.isfile(pivot_file)):
      original_operators = SortedSet([operator_info.operator.Operator(op) for op in sigmond.getOperatorBasis(pivot_file)])

    if pivot_exists and operators and original_operators != operators:
      logging.error("RotatedOperatorSet operators do not match operators on disk")
    elif pivot_exists and operators:
      return operators
    elif pivot_exists and not operators:
      return original_operators
    elif operators:
      return operators
    else:
      logging.error("RotatedOperatorSet cannot find operators")

  def getRotatedOperators(self):
    return data_handling.data_handler.DataHandler().getRotatedOperators(self)

  @property
  def name(self):
    return self._name
  
  @property
  def pivot_info(self):
    return self._pivot_info

  def __str__(self):
    return f"{self.name} - {self.pivot_info!s}"

  def __repr__(self):
    return f"{util.str_to_file(self.name)}_{self.pivot_info!r}"

  def __cmp(self):
    return (self.name, self.pivot_info, tuple(self.operators))

  def __hash__(self):
    return hash(self.__cmp())

  def __eq__(self, other):
    return self.__cmp() == other.__cmp()

  def __ne__(self, other):
    return not self.__eq__(other)

  def __lt__(self, other):
    return self.__cmp() < other.__cmp()

  def __le__(self, other):
    return self.__cmp() <= other.__cmp()

  def __gt__(self, other):
    return self.__cmp() > other.__cmp()

  def __ge__(self, other):
    return self.__cmp() >= other.__cmp()


def write_operators(op_file, operators, read_only=True, check_file=True):
  operators = SortedSet(operators)
  if check_file and os.path.isfile(op_file):
    ops_from_file = SortedSet(read_operators(opset_name))
    if ops_from_file != operators:
      logging.error(f"Different operator set already exists in '{op_file}'")
    else:
      logging.info(f"operator set found in '{op_file}'")

  else:
    f_handler = open(op_file, 'w')
    for op in operators:
      f_handler.write(f"{op.op_str()}\n")
    f_handler.close()

    if read_only:
      os.chmod(op_file, stat.S_IREAD|stat.S_IRGRP|stat.S_IROTH)

    logging.info(f"Wrote operator file to '{op_file}'")

def write_operators_to_yaml(yaml_file, name, operators, append=True):
  """ Add logging info """
  operators = SortedSet(operators)

  file_mode = 'a' if append else 'w'
  f_handler = open(yaml_file, file_mode)
  f_handler.write(f"x-{name}: &{name}\n")
  f_handler.write(f"  name: {name}\n")
  if len(operators) > 1:
    f_handler.write(f"  pivot_info:\n")
    f_handler.write(f"    <<: *PIVOT_INFO\n")
  f_handler.write(f"  operators:\n")
  for operator in operators:
    f_handler.write(f"    - {operator.op_str()}\n")
  f_handler.write("\n")


def read_operators(opset_name):
  op_file = operators_file(opset_name)
  if not os.path.isfile(op_file):
    logging.error(f"Cannot read operators, no such operator set: {opset_name}")

  f_handler = open(op_file, 'r')
  ops = list()
  for op_str in f_handler:
    op = operator_info.operator.Operator(op_str)
    ops.append(op)

  return ops
