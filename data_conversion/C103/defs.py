from typing import NamedTuple

import sourceTimeList

base_data_dir = "/media/ext2/research/data/raw_C103/v2/"
output_dir = "data"

class Ensemble(NamedTuple):
  name: str
  Nt: int
  replica: list
  sources: list
  parities: list

ensembles = [
    #Ensemble("cls21_c103", 96, ['r005', 'r006', 'r007'], [0, 24, 48, 72]),
    Ensemble(
      "cls21_c103",
      96,
      ['r007'],
      [
        (0, True), (0, False),
        (12, True), (12, False),
        (24, True), (24, False),
        (36, True), (36, False),
        (48, True), (48, False),
        (60, True), (60, False),
        (72, True), (72, False),
        (84, True), (84, False),
      ],
      ['fwd', 'bwd']
    ),
]

source_lists = {
    "cls21_c103_r005": sourceTimeList.r005,
    "cls21_c103_r006": sourceTimeList.r006,
    "cls21_c103_r007": sourceTimeList.r007,
}

config_indices = {
    "cls21_c103_r005": list(range(401)),
    "cls21_c103_r006": list(range(401)),
    "cls21_c103_r007": list(range(10)),
}

omissions = {
    "cls21_c103_r005": set(),
    "cls21_c103_r006": set(),
    "cls21_c103_r007": set(),
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

  def __str__(self):
    return f"(psq={self.psq}, irrep={self.irrep}, isospin={self.isospin}, S={self.strangeness})"

  def __repr__(self):
    return f"pSq{self.psq}-{self.irrep}-{self.isospin}-{self.strangeness}"

class EquivalentFrame(NamedTuple):
  mom: tuple
  row: int

  def __str__(self):
    return f"(P={self.mom}, row={self.row})"

  def __repr__(self):
    return f"P{''.join(list(map(str, self.mom)))}_{self.row}".replace('-', 'm')

