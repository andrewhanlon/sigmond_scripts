import logging
import regex

from operator_info.isospin import Isospin
from operator_info.channel import Channel
import sigmond
import utils.util as util

flavor_map = {
    'G': "glueball",
    'P': "pion",
    'E': "eta",
    'F': "phi",
    'K': "kaon",
    'k': "kbar",
    'N': "nucleon",
    'D': "delta",
    'S': "sigma",
    'L': "lambda",
    'X': "xi",
    'W': "omega",
    'EE00': "tquuuu1m",
    'EE01': "tquuuu1p",
    'EP20': "tquudu3m",
    'EP21': "tquudu3p",
    'PP00': "tqdudu1m",
    'PP01': "tqdudu1p",
    'PP20': "tqdudu3m",
    'PP21': "tqdudu3p",
    'PP40': "tqdudu5m",
    'PP41': "tqdudu5p",
    'KE10': "tqsuuu2m",
    'KE11': "tqsuuu2p",
    'KP10': "tqsudu2m",
    'KP11': "tqsudu2p",
    'KP30': "tqsudu4m",
    'KP31': "tqsudu4p",
    'FP20': "tqssdu3m",
    'FP21': "tqssdu3p",
    'EF00': "tquuss1m",
    'EF01': "tquuss1p",
    'KF10': "tqsuss2m",
    'KF11': "tqsuss2p",
    'FF00': "tqssss1m",
    'FF01': "tqssss1p",
}


class Operator:

  def __init__(self, operator):
    if isinstance(operator, str):
      self.operator_type = Operator.operator_type(operator)
      self.operator_info = sigmond.OperatorInfo(operator, self.operator_type)
    elif isinstance(operator, sigmond.OperatorInfo):
      self.operator_info = operator
      self.operator_type = sigmond.OpKind.GenIrrep if operator.isGenIrrep() else sigmond.OpKind.BasicLapH
    elif isinstance(operator, sigmond.BasicLapHOperatorInfo):
      self.operator_type = sigmond.OpKind.BasicLapH
      self.operator_info = sigmond.OperatorInfo(operator)
    elif isinstance(operator, sigmond.GenIrrepOperatorInfo):
      self.operator_type = sigmond.OpKind.GenIrrep
      self.operator_info = sigmond.OperatorInfo(operator)
    else:
      raise TypeError(f"Invalid operator type: {type(operator)}")

    self._channel = Channel.createFromOperator(self)

  @classmethod
  def createFromCompact(cls, compact_str):
    compact_str = compact_str.strip()
    if compact_str[0].isdigit(): # is gi_operator
      pattern = r"^(?P<iso_int>\d)(?P<strange>-?\d)(p(?P<psq>\d)|(?P<px>-?\d)" \
                r"(?P<py>-?\d)(?P<pz>-?\d))(?P<irrep>[A-Z]\d?[gu]?[mp]?)" \
                r"(_(?P<irreprow>\d))?-(?P<id_name>\S+)-(?P<id_index>\d+)$"

      match = regex.match(pattern, compact_str)
      if match is None:
        logging.critical(f"Invalid compact operator string passed: {compact_str}")

      matches = match.groupdict()
      isospin = Isospin(int(matches['iso_int'])).name
      op_str = f"iso{isospin} S={matches['strange']}"
      if matches['psq'] is not None:
        op_str += f" PSQ={matches['psq']}"
      else:
        op_str += f" P=({matches['px']},{matches['py']},{matches['pz']})"

      op_str += f" {matches['irrep']}"
      if matches['irreprow'] is not None:
        op_str += f"_{matches['irreprow']}"

      op_str += f" {matches['id_name']} {matches['id_index']}"

    else: # is bl_operator
      if compact_str[3].isdigit() or compact_str[3] == '-': # single hadron
        if compact_str[1].isdigit() or compact_str[1] == '-': # not a tetraquark
          pattern = r"^(?P<flavor>\S)"
        else: # a tetraquark
          pattern = r"^(?P<flavor>\S\S\d\d)"

        pattern += r"(?P<px>-?\d)(?P<py>-?\d)(?P<pz>-?\d)(?<irrep>[A-Z]\d?[gu]?[mp]?)" \
                   r"(?P<irreprow>\d)(?P<spat_type>[^\s\d]+)(?P<spat_id>\d+)?$"

        match = regex.match(pattern, compact_str)
        if match is None:
          logging.critical(f"Invalid compact operator string passed: {compact_str}")
        
        matches = match.groupdict()
        flavor = flavor_map[matches['flavor']]
        op_str = f"{flavor} P=({matches['px']},{matches['py']},{matches['pz']})" \
                 f" {matches['irrep']}_{matches['irreprow']} {matches['spat_type']}"

        if matches['spat_id'] is not None:
          op_str += f"_{matches['spat_id']}"

      else: # multi-hadron
        orig_compact_str = compact_str
        iso_index = int(compact_str[regex.search("\d", compact_str).start()])
        op_str = "iso" + Isospin(iso_index+1).name
        for flavor_ind in range(iso_index):
          flavor = flavor_map[compact_str[flavor_ind]]
          op_str += "_{flavor}"

        compact_str = compact_str[iso_index+1:]
        pattern = r"^(?P<irrep>[A-Z]\d?[gu]?[mp]?)(?P<irreprow>\d)(C(?P<cg>\d))?"
        for ind in range(iso_index):
          pattern += rf"(?P<px{ind}>-?\d)(?P<py{ind}>-?\d)(?P<pz{ind}>-?\d)" \
                     rf"(?P<irrep{ind}>[A-Z]\d?[gu]?[mp]?)(?P<spat_type{ind}>[^\s\d]+)" \
                     rf"(?P<spat_id{ind}>\d+)?"

        pattern += r"$"

        match = regex.match(pattern, compact_str)
        if match is None:
          logging.critical(f"Invalid compact operator string passed: {org_compact_str}")

        matches = match.groupdict()
        op_str += f" {matches['irrep']}_{matches['irreprow']}"
        if matches['cg'] is not None:
          op_str += f" CG_{matches['cg']}"

        for ind in range(iso_index):
          op_str += " [P=({px},{py},{pz}) {irrep} {spat_type}".format(
              px=matches[f'px{ind}'], py=matches[f'py{ind}'], pz=matches[f'pz{ind}'],
              irrep=matches[f'irrep{ind}'], spat_type=matches[f'spat_type{ind}'])

          if matches[f'spat_id{ind}'] is not None:
            op_str += "_{}".format(matches[f'spat_id{ind}'])

          op_str += "]"

    return cls(op_str)

  @property
  def compact_str(self):
    if self.operator_type == sigmond.OpKind.GenIrrep:
      operator = self.operator_info.getGenIrrep()
      if self.psq >= 10:
        logging.critical("PSQ >= 10 not supported")

      if abs(operator.getStrangeness()) >= 10:
        logging.critical("Strangeness >= 10 not supported")

      isospin_int = Isospin(operator.getIsospin()).to_int()
      compact_str = f"{isospin_int}{operator.getStrangeness()}"
      if operator.hasDefiniteMomentum():
        compact_str += f"{operator.getXMomentum()}" \
                       f"{operator.getYMomentum()}" \
                       f"{operator.getZMomentum()}"
      else:
        compact_str += f"p{operator.getMomentumSquared()}"

      compact_str += operator.getLGIrrep()
      if operator.getLGIrrepRow():
        compact_str += f"_{operator.getLGIrrepRow()}"

      compact_str += f"-{operator.getIDName()}-{operator.getIDIndex()}"

    else:
      operator = self.operator_info.getBasicLapH()
      if self.psq >= 10:
        logging.critical("PSQ >= 10 not supported")

      if abs(operator.getStrangeness()) >= 10:
        logging.critical("Strangeness >= 10 not supported")

      compact_str = operator.getFlavorCode()
      if operator.isTetraquark():
        compact_str += str(operator.getTetraquarkColorType())
      if operator.getNumberOfHadrons() == 1:
        compact_str += f"{operator.getXMomentum()}" \
                       f"{operator.getYMomentum()}" \
                       f"{operator.getZMomentum()}"

        compact_str += f"{operator.getLGIrrep()}" \
                       f"{operator.getLGIrrepRow()}" \
                       f"{operator.getHadronSpatialType(1)}"

        if not operator.isGlueball():
          compact_str += str(operator.getHadronSpatialIdNumber(1))

      else:
        compact_str += f"{operator.getLGIrrep()}{operator.getLGIrrepRow()}"

        if operator.getLGClebschGordonIdNum():
          compact_str += f"C{operator.getLGClebschGordonIdNum()}"

        for hadron in range(operator.getNumberOfHadrons()):
          compact_str += f"{operator.getHadronXMomentum(hadron+1)}" \
                         f"{operator.getHadronYMomentum(hadron+1)}" \
                         f"{operator.getHadronZMomentum(hadron+1)}"

          compact_str += f"{operator.getHadronLGIrrep(hadron+1)}" \
                         f"{operator.getHadronSpatialType(hadron+1)}"

          if not operator.isHadronGlueball(hadron+1):
            compact_str += str(operator.getHadronSpatialIdNumber(hadron+1))

    return compact_str

  @staticmethod
  def operator_type(op_str):
    isospin = op_str.split()[0]
    if isospin.startswith("iso") and "_" not in isospin:
      return sigmond.OpKind.GenIrrep
    else:
      return sigmond.OpKind.BasicLapH

  def __getattr__(self, name):
    try:
      return getattr(self.operator_info, name)
    except AttributeError:
      try:
        if self.operator_type is sigmond.OpKind.BasicLapH:
          return getattr(self.operator_info.getBasicLapH(), name)
        else:
          return getattr(self.operator_info.getGenIrrep(), name)
      except AttributeError:
        logging.critical(f"Operator has no attribute '{name}'")

  @property
  def channel(self):
    return self._channel

  @property
  def id_index(self):
    if self.operator_type is sigmond.OpKind.GenIrrep:
      return self.getIDIndex()
    else:
      logging.warning("BasicLaphOperatorInfo has no IDIndex")
      return 0

  @id_index.setter
  def id_index(self, in_index):
    if self.operator_type is sigmond.OpKind.GenIrrep:
      gen_op = self.getGenIrrep()
      gen_op.resetIDIndex(in_index)
      self.operator_info = sigmond.OperatorInfo(gen_op)
    else:
      logging.warning("Cannot update IDIndex of BasicLaphOperatorInfo")

  @property
  def ratio_op(self):
    if self.operator_type is sigmond.OpKind.GenIrrep:
      op_str_tokens = self.op_str().split()
      op_str_tokens[4] += "r"
      op_str = ' '.join(op_str_tokens)
      ratio_op = Operator(op_str)
      return ratio_op
    else:
      logging.warning("Cannot update IDIndex of BasicLaphOperatorInfo")
      return self

  @property
  def psq(self):
    if self.operator_type == sigmond.OpKind.BasicLapH:
      return self.getXMomentum()**2 + self.getYMomentum()**2 + self.getZMomentum()**2
    else:
      return self.getMomentumSquared()

  @property
  def vev(self):
    return self.channel.vev

  @property
  def is_rotated(self):
    return self.operator_type is sigmond.OpKind.GenIrrep and self.getIDName() == "ROT"

  @property
  def level(self):
    if self.is_rotated:
      return self.getIDIndex()
    
    return 0

  def __repr__(self):
    return util.str_to_file(self.__str__())

  def __str__(self):
    return self.operator_info.op_str()

  def __hash__(self):
    return hash(self.__repr__())

  def __eq__(self, other):
    return self.__repr__() == other.__repr__()

  def __ne__(self, other):
    return not self.__eq__(other)

  # @ADH - This is poorly written...
  def __comp_op(self):
    if self.operator_type is sigmond.OpKind.BasicLapH:
      _comp_list = [
          self.getIsospin(),
          self.getStrangeness(),
          self.psq,
          self.getLGIrrep(),
          self.getLGIrrepRow(),
          self.getNumberOfHadrons(),
      ]

      if self.getNumberOfHadrons() == 1:
        _comp_list.extend([self.getFlavorCode(), self.getHadronSpatialType(1)])

        if not self.isGlueball():
          _comp_list.append(self.getHadronSpatialIdNumber(1))

      else:
        flavors = []
        total_psq = 0
        psqs = []
        lg_irreps = []
        spatial_types = []
        spatial_ids = []
        moms = []

        for hadron_i in range(1, self.getNumberOfHadrons()+1):
          flavors.append(self.getHadronFlavor(hadron_i))
          psq = self.getHadronXMomentum(hadron_i)**2 + self.getHadronYMomentum(hadron_i)**2 + self.getHadronZMomentum(hadron_i)**2
          total_psq += psq
          psqs.append(psq)
          lg_irreps.append(self.getHadronLGIrrep(hadron_i))
          spatial_types.append(self.getHadronSpatialType(hadron_i))
          if not self.isHadronGlueball(hadron_i):
            spatial_ids.append(self.getHadronSpatialIdNumber(hadron_i))

          mom = (self.getHadronXMomentum(hadron_i), self.getHadronYMomentum(hadron_i), self.getHadronZMomentum(hadron_i))
          moms.append(mom)

        _comp_list.extend(flavors)
        _comp_list.extend(lg_irreps)
        _comp_list.extend(spatial_types)
        _comp_list.extend(spatial_ids)
        _comp_list.append(total_psq)
        _comp_list.extend(psqs)

        _comp_list.append(self.getIsospinClebschGordonIdNum())
        _comp_list.append(self.getLGClebschGordonIdNum())

        _comp_list.extend(moms) # @ADH - should this even be ranked at all?

      mom = (self.getXMomentum(), self.getYMomentum(), self.getZMomentum())
      _comp_list.append(mom)

    else:
      _comp_list = [
          self.getIsospin(),
          self.getStrangeness(),
          self.psq,
          self.getLGIrrep(),
          self.getLGIrrepRow(),
          self.getIDName(),
          self.getIDIndex(),
      ]

      if self.hasDefiniteMomentum():
        mom = (self.getXMomentum(), self.getYMomentum(), self.getZMomentum())
        _comp_list.append(mom)

    return tuple(_comp_list)

  def __cmp(self):
    return (str(self.operator_type), self.__comp_op(), self.operator_info)

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

