import copy
from sortedcontainers import SortedSet

import sigmond
from operator_info.isospin import Isospin

irreprows = {
    'A1': 1,
    'A2': 1,
    'E' : 2,
    'T1': 3,
    'T2': 3,
    'G1': 2,
    'G2': 2,
    'H' : 4,
    'B1': 1,
    'B2': 1,
    'G' : 2,
    'F1': 1,
    'F2': 1,
}


class Channel:

  def __init__(self, momentum, flavor, irrep="NONE", irreprow=0, ref_momentum=False):
    """Channel __init__ method

    Args:
      momentum (tuple of 3 ints): the definite momentum for
          a channel.
      flavor (tuple): the flavor of the channel
      irrep (str): the irrep for the channel
      irreprow (int): the irrep row for the channel. If missing, it
          is assumed that there is only one irrep row or the irrep row
          has been averaged over.
      ref_momentum (bool): assumed False if absent

    TODO:
      - Make use of sigmondbind.Momentum ?
    """

    self.momentum = tuple(momentum)
    self.flavor = flavor
    self.irrep = irrep
    self.irreprow = irreprow
    self.ref_momentum = ref_momentum

    if self.ref_momentum:
      self.momentum = self.refP
    
  @classmethod
  def createFromOperator(cls, operator):
    if operator.isBasicLapH():
      bl_op = operator.getBasicLapH()

      momentum=(bl_op.getXMomentum(), bl_op.getYMomentum(), bl_op.getZMomentum())

      isospin = bl_op.getIsospin().split('_')[0]  # Takes care of tetraquarks - Something better?
      if isospin.startswith("iso"):
        isospin = isospin[len("iso"):]
      isospin = Isospin(isospin).to_str
      strangeness = str(bl_op.getStrangeness())
      flavor = (isospin, strangeness)

      return cls(momentum=momentum, flavor=flavor, irrep=bl_op.getLGIrrep(), irreprow=bl_op.getLGIrrepRow())

    else:
      gi_op = operator.getGenIrrep()
      momentum = (gi_op.getXMomentum(), gi_op.getYMomentum(), gi_op.getZMomentum())
      return cls(momentum=momentum, flavor=gi_op.getFlavor(), irrep=gi_op.getLGIrrep(),
                 irreprow=gi_op.getLGIrrepRow(), ref_momentum=gi_op.isReferenceMomentum())

  @property
  def averaged(self):
    aver_chan = copy.copy(self)
    aver_chan.irreprow = 0
    aver_chan.momentum = self.refP
    aver_chan.ref_momentum = True
    return aver_chan

  @property
  def is_averaged(self):
    return self == self.averaged

  @property
  def at_rest(self):
    return self.psq == 0

  @property
  def refP(self):
    return tuple(sorted([abs(pi) for pi in self.momentum]))

  @property
  def refP_str(self):
     ref_p = self.refP
     return f"Pref{ref_p[0]}{ref_p[1]}{ref_p[2]}"

  @property
  def psq(self):
    return self.momentum[0]**2 + self.momentum[1]**2 + self.momentum[2]**2

  @property
  def vev(self):
    return (self.irrep == "A1g" or self.irrep == "A1gp") and self.flavor == ("0","0")

  def getRotatedOp(self, level=0):
    return self.getGIOperator("ROT", level)

  def getGIOperator(self, obs_name, obs_id=0):
    op_str = f"{self!s} {obs_name} {obs_id}"
    return sigmond.GenIrrepOperatorInfo(op_str)

  @property
  def mom_str(self):
    _mom_str = "Pref" if self.ref_momentum else "P"
    _mom_str += "{self.momentum[0]}{self.momentum[1]}{self.momentum[2]}"
    return _mom_str

  @property
  def irrep_refP_key(self):
    return f"{self.irrep}_Pref{self.refP}"

  def __str__(self):
    _str = "Pref=" if self.ref_momentum else "P="
    _str += f"({self.momentum[0]},{self.momentum[1]},{self.momentum[2]})"

    if self.irrep != "NONE":
      _str += f" {self.irrep}"
      if self.irreprow != 0:
        _str += f"_{self.irreprow}"

    _str += " flavor="
    for f_i in self.flavor:
      _str += f"{f_i},"
    _str = _str[:-1]

    return _str
  

  def __repr__(self):
    _str = "Pr" if self.ref_momentum else "P"
    _str += f"{self.momentum[0]}_{self.momentum[1]}_{self.momentum[2]}"

    if self.irrep != "NONE":
      _str += f"_{self.irrep}"
      if self.irreprow != 0:
        _str += f"_{self.irreprow}"

    _str += "_F"
    for f_i in self.flavor:
      _str += f_i

    return _str.replace('-', 'm')
  
  def __cmp(self):
    return (self.flavor, self.psq, self.irrep, self.irreprow, self.momentum)

  def __hash__(self):
    return hash(repr(self))

  def __eq__(self, other):
    return repr(self) == repr(other)

  def __ne__(self, other):
    return not self.__eq__(other)

  def __lt__(self, other):
    if isinstance(other, self.__class__):
      return self.__cmp() < other.__cmp()
    return NotImplemented

  def __le__(self, other):
    if isinstance(other, self.__class__):
      return self.__cmp() <= other.__cmp()
    return NotImplemented

  def __gt__(self, other):
    if isinstance(other, self.__class__):
      return self.__cmp() > other.__cmp()
    return NotImplemented

  def __ge__(self, other):
    if isinstance(other, self.__class__):
      return self.__cmp() >= other.__cmp()
    return NotImplemented
