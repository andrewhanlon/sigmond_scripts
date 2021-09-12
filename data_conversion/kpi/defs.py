from typing import NamedTuple

base_data_dir = "/disk2/research/data/kpi/raw"
output_dir = "/disk1/research/projects/kpi/data"

class Ensemble(NamedTuple):
  name: str
  replica: list
  sources: list


ensembles = [
    Ensemble(
        "cls21_n203",
        ['r000', 'r001'],
        [(32, 'fwd'), (52, 'fwd')]
    ),
    Ensemble(
        "cls21_n200",
        ['r000'],
        [(32, 'fwd'), (52, 'fwd')]
    ),
    Ensemble(
        "cls21_d200",
        ['r000'],
        [(35, 'fwd'), (92, 'bwd')]
    ),
]

omissions = {
    'cls21_n203_r000': set(range(1, 756, 2)),
    'cls21_n203_r001': set(range(0, 787, 2)),
    'cls21_n203': set(range(1, 756, 2)).union(set(range(756, 1544, 2))),
    "cls21_n200_r000": set(range(1, 856, 2)),
    "cls21_n200": set(range(1, 856, 2)),
    'cls21_d200_r000': set([2000]),
    'cls21_d200': set([2000]),
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

    isospin = name_map.get(self.isospin, self.isospin)

    return f"{mom_str}_{self.irrep}_{self.irrep_row}_{isospin}_{strangeness_str}"

name_map = {
    'kaon': 'doublet',
    'pion': 'triplet',
    'eta': 'singlet',
    'phi': 'singlet',
}
