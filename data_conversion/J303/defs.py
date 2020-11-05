from typing import NamedTuple

base_data_dir = "/media/ext2/research/data/J303/"
output_dir = "data"

class Ensemble(NamedTuple):
  name: str
  Nt: int
  replica: list
  sources: list

ensembles = [
    Ensemble("cls21_j303", 192, ['r003'], [48, 74, 118, 144]),
]

config_indices = {
    "cls21_j303_r003": list(range(1, 9, 2)) + list(range(11, 82, 2)),
}

omissions = {
    "cls21_j303_r003": set(list(range(1, 8, 2)) + [8] + list(range(9, 48, 2))),
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
    strangeness_str = f"S{self.strangeness}".replace('-', 'm')
    mom_str = f"P{self.momentum[0]}{self.momentum[1]}{self.momentum[2]}".replace('-', 'm')

    return f"{self.isospin}_{strangeness_str}/{mom_str}_{self.irrep}_{self.irrep_row}"
