from typing import NamedTuple

base_data_dir = "/latticeQCD/raid0/ahanlon/data/raw_three_hadrons/"
output_dir = "data"

class Ensemble(NamedTuple):
  name: str
  replica: list
  sources: list
  flavor_channels: list

'''
ensembles = [
    Ensemble("cls21_n203", ['r000', 'r001'], [32, 52], ['kaons','pions']),
    Ensemble("cls21_n200", ['r000', 'r001'], [32, 52], ['kaons','pions']),
    Ensemble("cls21_nd00", ['r000'], [35], ['kaons','pions']),
]
'''
ensembles = [
    Ensemble("cls21_d200", ['r000'], [35], ['kaons','pions']),
]

omissions = {
    'cls21_n203_r000': set(range(1, 756, 2)),
    'cls21_n203_r001': set(range(0, 787, 2)),
    'cls21_n203': set(range(1, 756, 2)).union(set(range(756, 1544, 2))),
    'cls21_n200_r000': set(),
    'cls21_n200_r001': set(),
    'cls21_n200': set(),
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

    return f"{mom_str}_{self.irrep}_{self.irrep_row}_{self.isospin}_{strangeness_str}"
