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

  EXTRA_INFO_KEYS = ['irreprow', 'momentum', 'momentum_squared']

  def __init__(self, isospin, strangeness, irrep, **extra_info):
    """Channel __init__ method

    Args:
      isospin (str): the isospin for the channel
      strangeness (int): the strangeness of the channel
      irrep (str): the irrep for the channel
      **irreprow (int): the irrep row for the channel. If missing, it
          is assumed that there is only one irrep row or the irrep row
          has been averaged over.
      **momentum (tuple of 3 ints): the definite momentum for
          a channel.
      **momentum_squared (int): the momentum squared for the channel,
          which assumes the equivalent momentum frames have been
          averaged over. Note that a momentum squared of 0 is treated
          as a momentum channel with P = (0, 0, 0)

    TODO:
      - Make use of sigmondbind.Momentum ?
    """

    self.isospin = Isospin(isospin).value
    self.strangeness = strangeness
    self.irrep = irrep

    for attr, value in extra_info.items():
      if attr not in self.EXTRA_INFO_KEYS:
        logging.error(f"Unrecognized key {attr} in Channel")

      setattr(self, attr, value)

    if hasattr(self, 'momentum'):
      self.momentum = tuple(self.momentum)
    
    if hasattr(self, 'momentum') and hasattr(self, 'momentum_squared'):
      raise ValueError("Channel can not have both 'momentum' and 'momentum_squared'")

    if hasattr(self, 'momentum') and self.momentum == (0,0,0):
      del self.momentum
      self.momentum_squared = 0

  @classmethod
  def createFromOperator(cls, operator):
    if operator.isBasicLapH():
      bl_op = operator.getBasicLapH()

      # Takes care of tetraquarks - Something better?
      isospin = bl_op.getIsospin().split('_')[0]
      if isospin.startswith("iso"):
        isospin = isospin[len("iso"):]

      return cls(isospin=isospin, strangeness=bl_op.getStrangeness(),
                 irrep=bl_op.getLGIrrep(), irreprow=bl_op.getLGIrrepRow(),
                 momentum=(bl_op.getXMomentum(), bl_op.getYMomentum(),
                           bl_op.getZMomentum()))

    else:
      gi_op = operator.getGenIrrep()
      kw_args = {
          "isospin"     : gi_op.getIsospin(),
          "strangeness" : gi_op.getStrangeness(),
          "irrep"       : gi_op.getLGIrrep(),
      }

      if gi_op.getLGIrrepRow() > 0:
        kw_args["irreprow"] = gi_op.getLGIrrepRow()

      if gi_op.hasDefiniteMomentum():
        kw_args["momentum"] = (gi_op.getXMomentum(), gi_op.getYMomentum(),
                               gi_op.getZMomentum())
      else:
        kw_args["momentum_squared"] = gi_op.getMomentumSquared()

      return cls(**kw_args)

  @property
  def averaged(self):
    aver_chan = copy.copy(self)
    if hasattr(aver_chan, 'irreprow'):
      del aver_chan.irreprow
    if hasattr(aver_chan, 'momentum'):
      del aver_chan.momentum
      aver_chan.momentum_squared = self.psq

    return aver_chan

  @property
  def is_averaged(self):
    return self == self.averaged

  @property
  def at_rest(self):
    return self.psq == 0

  @property
  def psq(self):
    if hasattr(self, 'momentum_squared'):
      return self.momentum_squared
    else:
      return self.momentum[0]**2 + self.momentum[1]**2 + self.momentum[2]**2

  @property
  def vev(self):
    return (self.irrep == "A1g" or self.irrep == "A1gp") and self.isospin == "singlet" and self.strangeness == 0

  def getRotatedOp(self, level=0):
    return self.getGIOperator("ROT", level)

  def getGIOperator(self, obs_name, obs_id=0):
    op_str = "iso{} S={}".format(self.isospin, self.strangeness)
    if hasattr(self, "momentum"):
      op_str += " P=({},{},{})".format(self.momentum[0], self.momentum[1],
                                       self.momentum[2])

    if hasattr(self, "momentum_squared"):
      op_str += " PSQ={}".format(self.momentum_squared)

    op_str += " {}".format(self.irrep)
    if hasattr(self, "irreprow"):
      op_str += "_{}".format(self.irreprow)

    op_str += " {} {}".format(obs_name, obs_id)

    return sigmond.GenIrrepOperatorInfo(op_str)

  @property
  def mom_str(self):
    if hasattr(self, "momentum"):
      return f"P{self.momentum[0]}{self.momentum[1]}{self.momentum[2]}"
    else:
      return f"PSQ{self.momentum_squared}"

  @property
  def irrep_psq_key(self):
    return f"{self.irrep}_PSQ{self.psq}"

  def __str__(self):
    _str = "iso{} S={} {}".format(self.isospin, self.strangeness, self.irrep)
    if hasattr(self, "irreprow"):
      _str += "_{}".format(self.irreprow)

    if hasattr(self, "momentum"):
      _str += " P=({},{},{})".format(self.momentum[0], self.momentum[1], self.momentum[2])
    elif hasattr(self, "momentum_squared"):
      _str += " PSQ={}".format(self.momentum_squared)

    return _str
  

  def __repr__(self):
    _str = "iso{}_S{}_{}".format(self.isospin, self.strangeness, self.irrep)
    if hasattr(self, "irreprow"):
      _str += "_{}".format(self.irreprow)

    if hasattr(self, "momentum"):
      _str += "_P{}{}{}".format(self.momentum[0], self.momentum[1], self.momentum[2])
    elif hasattr(self, "momentum_squared"):
      _str += "_P{}".format(self.momentum_squared)

    return _str.replace('-', 'm')
  
  def __cmp(self):
    irreprow = getattr(self, 'irreprow', 0)
    momentum = getattr(self, 'momentum', None)

    return (self.isospin, self.strangeness, self.psq, self.irrep, irreprow, momentum)

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
