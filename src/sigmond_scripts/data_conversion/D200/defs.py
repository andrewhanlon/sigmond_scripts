from typing import NamedTuple


class Ensemble(NamedTuple):
  name: str
  Nt: int
  replica: list
  sources: list


parity_name = {
    True: 'fwd',
    False: 'bwd',
}

ensemble = Ensemble("cls21_d200", 128, ['r000'], [(35, True), (64, True), (64, False), (92, False)])
ensemble_name = "cls21_s64_t128_D200"
search_source = (92, False)

config_indices = {
    'r000': list(range(1, 2001, 1)),
}

omissions = {
    'r000': set([2000]),
}


isospin_map = {
    'nucleon': 'doublet',
    'pion': 'triplet',
    'kaon': 'doublet',
    'kbar': 'doublet',
    'delta': 'quartet',
    'lambda': 'singlet',
    'sigma': 'triplet',
    'xi': 'doublet',
}

class Channel(NamedTuple):
  momentum: tuple
  irrep: str
  irrep_row: int
  isospin: str
  strangeness: int

  def __repr__(self):
    mom_str = f"P{self.momentum[0]}{self.momentum[1]}{self.momentum[2]}".replace('-', 'm')
    strangeness_str = f"S{self.strangeness}".replace('-', 'm')

    return f"{mom_str}_{self.irrep}_{self.irrep_row}_iso{self.isospin}_{strangeness_str}"

  def data_channel_str(self):
    mom_str = f"P{self.momentum[0]}{self.momentum[1]}{self.momentum[2]}".replace('-', 'm')

    return f"{mom_str}_{self.irrep}_{self.irrep_row}"

  def averaged_channel_str(self):
    psq = self.momentum[0]**2 + self.momentum[1]**2 + self.momentum[2]**2
    strangeness_str = f"S{self.strangeness}".replace('-', 'm')
    return f"PSQ{psq}_{self.irrep}_iso{self.isospin}_{strangeness_str}"

  def iso_strange_str(self):
    strangeness_str = f"S{self.strangeness}".replace('-', 'm')
    return f"iso{self.isospin}_{strangeness_str}"




class AveragedChannel(NamedTuple):
  psq: int
  irrep: str
  isospin: str
  strangeness: int

  def iso_strange_str(self):
    strangeness_str = f"S{self.strangeness}".replace('-', 'm')
    return f"iso{self.isospin}_{strangeness_str}"

  def data_channel_str(self):
    return f"PSQ{self.psq}_{self.irrep}"

  def __str__(self):
    return f"(psq={self.psq}, irrep={self.irrep}, isospin={self.isospin}, S={self.strangeness})"

  def __repr__(self):
    return f"PSQ{self.psq}-{self.irrep}-{self.isospin}-{self.strangeness}"

  def __hash__(self):
    return hash(repr(self))

  def __eq__(self, other):
    return repr(self) == repr(other)
