from collections import namedtuple

Ensemble = namedtuple('Ensemble', ['name', 'dir_name', 'type', 'su3', 'open', 'Nt', 'replica', 'replica_str', 'modes', 'srcs', 't0', 'ts'])

# ensemble info
ensembles = [
    Ensemble("B450", "B450", 'cls', True, False, 64, ('r000',), 'r000', 32, 8, 0, 33),
    Ensemble("N300", "N300", 'cls', True, True, 128, ('r001','r002'), 'r001-002', 32, 12, 0, 48),
    Ensemble("N202", "N202", 'cls', True, True, 128, ('r001',), 'r001', 68, 8, 0, 40),
]

backward_prop_skip = {
    'a064_m400_mL6.4_trMc': [],
    'a094_m400_mL6.0_trMc': [],
    'a12_m400_mL6.0_trMc': [],
    'A653': [],
    'B450': [],
    'B451': [],
    'B452': [],
    'H101': [0],
    'H102': [0],
    'H107': [0],
    'H200': [0,1,2,3],
    'J500': [0,1,2,3],
    'N200': [0,1],
    'N202': [0,1],
    'N300': [0,1,2,3,4,5],
    'N451': [],
    'U102': [0],
    'U103': [0],
    'U103_60modes': [0],
    'U103_4': [0],
    'E1':   [],
    'E5':   [],
}

forward_prop_skip = {
    'a064_m400_mL6.4_trMc': [],
    'a094_m400_mL6.0_trMc': [],
    'a12_m400_mL6.0_trMc': [],
    'A653': [],
    'B450': [],
    'B451': [],
    'B452': [],
    'H101': [3],
    'H102': [3],
    'H107': [3],
    'H200': [4,5,6,7],
    'J500': [8,9,10,11],
    'N200': [6,7],
    'N202': [6,7],
    'N300': [6,7,8,9,10,11],
    'N451': [],
    'U102': [4],
    'U103': [4],
    'U103_60modes': [4],
    'U103_4': [4],
    'E1':   [],
    'E5':   [],
}


def convert_irrep(irrep, psq):
  if psq == 2 and irrep == "B1":
    return "B2"
  elif psq == 2 and irrep == "B2":
    return "B1"
  elif irrep == "E2":
    return "E"

  return irrep.replace("+", "g").replace("-", "u")


def get_opstr(pref, irrep, isospin, op_str):
  psq = int(pref[1])**2 + int(pref[2])**2 + int(pref[3])**2
  pref = f"Pref=({pref[1]},{pref[2]},{pref[3]})"
  irrep = convert_irrep(irrep, psq)
  flavor = f"{isospin[1:]},0"
  op_str = op_str.replace(' ', '_').replace('*', 'S')

  return f"{pref} {irrep} flavor={flavor} {op_str}"
