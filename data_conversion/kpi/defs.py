from typing import NamedTuple

base_data_dir = "/disk2/research/data/kpi/raw/"
output_dir = "/disk2/research/data/kpi/hdf5/"

class Ensemble(NamedTuple):
  name: str
  replica: list
  sources: list

ensembles = [
    #Ensemble("cls21_d200", ['r000'], [35]),
    Ensemble("cls21_d200", ['r000'], [35, 92]),
    #Ensemble("cls21_n200", ['r000'], [32, 52]),
    #Ensemble("cls21_n203", ['r000', 'r001'], [32, 52]),
]

omissions = {
    "cls21_d200_r000": set([2000]),
    "cls21_n200_r000": set(range(1, 856, 2)),
    "cls21_n203_r000": set(range(1, 757, 2)),
    "cls21_n203_r001": set(range(0, 788, 2)),
}

configs = {
    "cls21_d200_r000": list(range(1, 2001, 1)),
    "cls21_n200_r000": list(range(1, 856, 2)),
    "cls21_n203_r000": list(range(1, 757, 2)),
    "cls21_n203_r001": list(range(2, 787, 2)),
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

  def iso_strange_str(self):
    strangeness_str = f"S{self.strangeness}".replace('-', 'm')
    return f"{self.isospin}_{strangeness_str}"


name_map = {
    'kaon': 'doublet',
    'pion': 'triplet',
    'eta': 'singlet',
    'phi': 'singlet',
}


allMom = {
    0: [(0,0,0),],
    1: [(0,0,1),(0,0,-1),(0,1,0),(0,-1,0),(1,0,0),(-1,0,0),],
    2: [(0, 1, 1), (0, 1, -1), (0, -1, 1), (0, -1, -1), (1, 0, 1), (1, 0, -1), (-1, 0, 1), (-1, 0, -1), (1, 1, 0), (1, -1, 0), (-1, 1, 0), (-1, -1, 0)],
    3: [(1, 1, 1), (1, -1, 1), (-1, 1, 1), (-1, -1, 1), (1, 1, -1), (1, -1, -1), (-1, 1, -1), (-1, -1, -1)],
    4: [(1, 1, 1), (1, -1, 1), (-1, 1, 1), (-1, -1, 1), (1, 1, -1), (1, -1, -1), (-1, 1, -1), (-1, -1, -1)],
}

phases = {}
phases.update({'kaon P=('+','.join(map(str, aM))+') B1_1 SS_2' : -1.j for aM in allMom[2]})
phases.update({'kaon P=('+','.join(map(str, aM))+') B1_1 SS_3' : 1.j for aM in allMom[2]})
phases.update({'kaon P=('+','.join(map(str, aM))+') B2_1 SS_0' : 1.j for aM in allMom[2]})
phases.update({'kaon P=('+','.join(map(str, aM))+') B2_1 SS_1' : 1.j for aM in allMom[2]})
phases.update({'kaon P=('+','.join(map(str, aM))+') E_'+str(iR)+' SS_3' : 1.j for aM in allMom[3] for iR in [1,2]})
phases.update({'kaon P=('+','.join(map(str, aM))+') E_'+str(iR)+' SS_2' : 1.j for aM in allMom[3] for iR in [1,2]})
phases.update({'kaon P=('+','.join(map(str, aM))+') E_'+str(iR)+' SS_1' : 0.5**0.5*(1.+1.j) for aM in allMom[3] for iR in [1,2]})
phases.update({'kaon P=('+','.join(map(str, aM))+') E_'+str(iR)+' SS_0' : 0.5**0.5*(1.+1.j) for aM in allMom[3] for iR in [1,2]})
