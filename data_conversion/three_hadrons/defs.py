from typing import NamedTuple

base_data_dir = "/disk2/research/data/three_hadrons/raw"
output_dir = "/disk2/research/data/three_hadrons/sigmond"

class Ensemble(NamedTuple):
  name: str
  replica: list
  sources: list
  flavor_channels: list

'''
ensembles = [
    Ensemble("cls21_n203", ['r000', 'r001'], [32, 52], ['kaons', 'pions', 'KKp_Kpp']),
    Ensemble("cls21_n200", ['r000', 'r001'], [32, 52], ['kaons', 'pions', 'KKp_Kpp']),
    Ensemble("cls21_d200", ['r000'], [35, 92], ['kaons', 'pions', 'KKp_Kpp']),
]
'''
ensembles = [
'''
    Ensemble(
        "cls21_n203",
        ['r000', 'r001'],
        [(32, 'fwd'), (52, 'fwd')],
        ['kaons', 'pions', 'KKp_Kpp']
    ),
    Ensemble(
        "cls21_n200",
        ['r000', 'r001'],
        [(32, 'fwd'), (52, 'fwd')],
        ['pions', 'kaons', 'KKp_Kpp']
    ),
    Ensemble(
        "cls21_d200",
        ['r000'],
        [(35, 'fwd'), (92, 'bwd')],
        ['pions', 'kaons', 'KKp_Kpp']
    ),
'''
    Ensemble(
        "cls21_e250",
        ['r001'],
        [(0, 'fwd'), (0, 'bwd'), (1, 'fwd'), (1, 'bwd'), (2, 'fwd'), (2, 'bwd'), (3, 'fwd'), (3, 'bwd')],
        ['pions', 'kaons', 'KKp', 'Kpp']
    ),
]
'''
ensembles = [
    Ensemble(
        "cls21_e250",
        ['r001'],
        [(0, 'fwd'), (0, 'bwd'), (1, 'fwd'), (1, 'bwd'), (2, 'fwd'), (2, 'bwd'), (3, 'fwd'), (3, 'bwd')],
        ['pions', 'kaons', 'KKp', 'Kpp']
    ),
]
'''

omissions = {
    'cls21_n203_r000': set(range(1, 756, 2)),
    'cls21_n203_r001': set(range(0, 787, 2)),
    'cls21_n203': set(range(1, 756, 2)).union(set(range(756, 1544, 2))),
    'cls21_n200_r000': set(),
    'cls21_n200_r001': set(),
    'cls21_n200': set(),
    'cls21_d200_r000': set([2000]),
    'cls21_d200': set([2000]),
    "cls21_e250_r001": set(list(range(0, 1009))) -  set(list(range(0, 501, 4)) + list(range(502, 901, 2)) + list(range(904, 1009, 4))),
    "cls21_e250": set(list(range(0, 1009))) -  set(list(range(0, 501, 4)) + list(range(502, 901, 2)) + list(range(904, 1009, 4))),
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
