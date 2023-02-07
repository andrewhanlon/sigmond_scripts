from typing import NamedTuple

base_data_dir = "/disk2/research/data/rho/laph/"
output_dir = "data"

class Ensemble(NamedTuple):
  name: str
  Nt: int
  replica: list
  sources: list

parity_name = {
    True: 'fwd',
    False: 'bwd',
}

ensembles = [
    #Ensemble(
    #  "cls21_d200",
    #  128, 
    #  ['r000'],
    #  [
    #    (35, True),
    #    (92, False),
    #  ]
    #),
    #Ensemble(
    #  "cls21_j303",
    #  192, 
    #  ['r003'],
    #  [
    #    (48, True),
    #    (74, True),
    #    (118, False),
    #    (144, False),
    #  ]
    #),
    Ensemble(
      "cls21_e250",
      192, 
      ['r001'],
      [
        (0, True), (0, False),
        (1, True), (1, False),
        (2, True), (2, False),
        (3, True), (3, False),
      ]
    ),
]

source_lists = {}

configs = {
    "cls21_j303_r003": list(range(1, 1073, 2)),
    "cls21_d200_r000": list(range(1, 2001, 1)),
    "cls21_e250_r001": list(range(1, 502, 4)) + list(range(503, 902, 2)) + list(range(905, 1010, 4)),
}

omissions = {
    "cls21_d200": set(),
    "cls21_d200_r000": set(),

    "cls21_j303": set(list(range(1, 1073, 2))),
    "cls21_j303_r003": set(list(range(1, 1073, 2))),

    "cls21_e250": set(list(range(0, 1009))) -  set(list(range(0, 501, 4)) + list(range(502, 901, 2)) + list(range(904, 1009, 4))),
    "cls21_e250_r001": set(list(range(0, 1010))) -  set(list(range(0, 501, 4)) + list(range(502, 901, 2)) + list(range(904, 1009, 4))),
}


config_indices = {
}

isospin_map = {
    'nucleon': 'doublet',
    'pion': 'triplet',
    'kaon': 'doublet',
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

    return f"{mom_str}_{self.irrep}_{self.irrep_row}_{self.isospin}_{strangeness_str}"

  def data_channel_str(self):
    mom_str = f"P{self.momentum[0]}{self.momentum[1]}{self.momentum[2]}".replace('-', 'm')

    return f"{mom_str}_{self.irrep}_{self.irrep_row}"

  def averaged_channel_str(self):
    psq = self.momentum[0]**2 + self.momentum[1]**2 + self.momentum[2]**2
    strangeness_str = f"S{self.strangeness}".replace('-', 'm')
    return f"PSQ{psq}_{self.irrep}_{self.isospin}_{strangeness_str}"

  def iso_strange_str(self):
    strangeness_str = f"S{self.strangeness}".replace('-', 'm')
    return f"{self.isospin}_{strangeness_str}"




class AveragedChannel(NamedTuple):
  psq: int
  irrep: str
  isospin: str
  strangeness: int

  def iso_strange_str(self):
    strangeness_str = f"S{self.strangeness}".replace('-', 'm')
    return f"{self.isospin}_{strangeness_str}"

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
