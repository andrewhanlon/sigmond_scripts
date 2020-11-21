from typing import NamedTuple

base_data_dir = "/media/ext2/research/data/ben/kpi_ins100/"
output_dir = "data/kpi_ins100/"

class Ensemble(NamedTuple):
  name: str
  Nt: int
  replica: list
  sources: list

ensembles = [
    Ensemble("cls21_n200", 128, ['r001'], [32, 52]),
]

omissions = {
    "cls21_n200_r001": set(),
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
