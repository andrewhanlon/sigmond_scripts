output_dir = "data"

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

data_info = {
    "N203": {
      "sources": [32, 52],
      "ensemble_name": ["cls21_n203_r000", "cls21_n203_r001"],
      "combined_name": "cls21_n203",
    },
    "D200": {
      "sources": [35],
      "ensemble_name": ["cls21_d200_r000"],
      "combined_name": "cls21_d200_r000",
    },
}

ensembles = [
    #"N203",
    "D200",
]

isospsins = [
    "isodoublet",
    "isoquartet",
]

single_hadrons = [
    "pion",
    "kaon",
]

omissions = {
    "cls21_n203_r000": set(range(1, 757, 2)),
    "cls21_n203_r001": set(range(0, 788, 2)),
    "cls21_n203": set(range(1, 757, 2)).union(set(range(756, 1544, 2))),
    "cls21_d200_r000": set([2000])
}

configs = {
    "cls21_n203_r000": list(range(1, 757, 2)),
    "cls21_n203_r001": list(range(2, 787, 2)),
    "cls21_d200_r000": list(range(1, 2001, 1)),
}

spatial_extent = {
    "cls21_n203_r000": 48,
    "cls21_n203_r001": 48,
    "cls21_d200_r000": 64,
}

single_hadron_irreps = {
    'kaon': {
        0: 'A1u',
        1: 'A2',
        2: 'A2',
        3: 'A2',
        4: 'A2',
    },
    'pion': {
        0: 'A1um',
        1: 'A2m',
        2: 'A2m',
        3: 'A2m',
        4: 'A2m',
    },
}

single_hadron_isospin = {
    'kaon': 'isodoublet',
    'pion': 'isotriplet',
}

single_hadron_strangeness = {
    'kaon': 1,
    'pion': 0,
}


name_map = {
    'kaon': 'K',
    'pion': 'pi',
    'eta': 'e',
    'phi': 'f',
}
