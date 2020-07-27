from collections import namedtuple

Ensemble = namedtuple('Ensemble', ['name', 'su3', 'open', 'Nt', 'replica', 'replica_str', 'modes', 'srcs', 't0', 'ts'])
Channel = namedtuple('Channel', ['P', 'irrep', 'isospin', 'strangeness', 'flavor'])

ensembles = [
    Ensemble("A653", True, False, 48, ('r000',), 'r000', 32, 4, 0, 24),
    Ensemble("B450", True, False, 64, ('r000',), 'r000', 32, 4, 0, 32),
    Ensemble("H101", True, True, 96, ('r000','r001',), 'r000-001', 48, 4, 0, 24),
    Ensemble("U103", True, True, 128, ('r001','r002','r003',), 'r001-003', 20, 5, 0, 24),
    Ensemble("H200", True, True, 96, ('r000','r001'), 'r000-001', 20, 8, 0, 32),
    Ensemble("N300", True, True, 128, ('r001','r002'), 'r001-002', 32, 6, 0, 40),
    Ensemble("N202", True, True, 128, ('r001',), 'r001', 68, 8, 0, 32),
    Ensemble("U102", False, True, 128, ('r001', 'r002',), 'r001-002', 20, 5, 0, 24),
    Ensemble("E1", True, False, 64, ('',), '', 30, 8, 0, 32),
    Ensemble("E5", True, False, 64, ('f', 'g',), 'fg', 30, 4, 0, 32),
]

backward_prop_skip = {
    'A653': [],
    'B450': [],
    'H101': [0],
    'H200': [0,1,2,3],
    'N300': [0,1,2],
    'N202': [0],
    'U103': [0],
    'U102': [0],
    'E1':   [],
    'E5':   [],
}

forward_prop_skip = {
    'A653': [],
    'B450': [],
    'H101': [3],
    'H200': [4,5,6,7],
    'N300': [3,4,5],
    'N202': [7],
    'U103': [4],
    'U102': [4],
    'E1':   [],
    'E5':   [],
}

pseudoscalar_backward_prop_skip = {
    'A653': [],
    'B450': [],
    'H101': [0,1],
    'H200': [0,1,2,3,4,5],
    'N300': [0,1,2,3],
    'N202': [0,1,2],
    'U103': [0,1],
    'U102': [0,1],
    'E1':   [],
    'E5':   [],
}

pseudoscalar_forward_prop_skip = {
    'A653': [],
    'B450': [],
    'H101': [2,3],
    'H200': [2,3,4,5,6,7],
    'N300': [2,3,4,5],
    'N202': [5,6,7],
    'U103': [3,4],
    'U102': [3,4],
    'E1':   [],
    'E5':   [],
}


pseudoscalar_modes = {
    'A653': 32,
    'B450': 32,
    'H101': 144,
    'H200': 20,
    'U103': 20,
    'N300': 32,
    'N202': 68,
    'U102': 20,
    'E1':   56,
    'E5':   60,
}

pseudoscalar_sources = {
    'A653': [0, 0, 0, 0],
    'B450': [0, 0, 0, 0],
    'H101': [24, 40, 55, 71],
    'H200': [0, 4, 8, 12, 16, 20, 24, 28],
    'N300': [0, 8, 16, 24, 32, 40],
    'N202': [0, 6, 12, 18, 24, 30, 36, 42],
    'U103': [32, 48, 64, 79, 95],
    'U102': [32, 48, 64, 79, 95],
    'E1':   [0, 0, 0, 0, 0, 0, 0, 0],
    'E5':   [0, 0, 0, 0],
}

pseudoscalar_names = {
    'A653': ['ps'],
    'B450': ['ps'],
    'H101': ['ps'],
    'H200': ['ps'],
    'N300': ['ps'],
    'N202': ['ps'],
    'U103': ['ps'],
    'U102': ['pion', 'kaon', 'eta_s'],
    'E1':   ['ps'],
    'E5':   ['pion'],
}

pseudoscalar_op_strs = {
    'A653': {
      'ps': "isotriplet S=0 PSQ=0 A1um pion 0",
    },
    'B450': {
      'ps': "isotriplet S=0 PSQ=0 A1um pion 0",
    },
    'H101': {
      'ps': "isotriplet S=0 PSQ=0 A1um pion 0",
    },
    'H200': {
      'ps': "isotriplet S=0 PSQ=0 A1um pion 0",
    },
    'N300': {
      'ps': "isotriplet S=0 PSQ=0 A1um pion 0",
    },
    'N202': {
      'ps': "isotriplet S=0 PSQ=0 A1um pion 0",
    },
    'U103': {
      'ps': "isotriplet S=0 PSQ=0 A1um pion 0",
    },
    'U102': {
      'pion': "isotriplet S=0 PSQ=0 A1um pion 0",
      'kaon': "isodoublet S=1 PSQ=0 A1u kaon 0",
      'eta_s': "isosinglet S=0 PSQ=0 A1up eta 0",
    },
    'E1': {
      'ps': "isotriplet S=0 PSQ=0 A1um pion 0",
    },
    'E5': {
      'pion': "isotriplet S=0 PSQ=0 A1um pion 0",
    },
}

tsrc_files = {
    'H200': 'tsrc_list_H200',
    'N300': 'tsrc_list_N300',
    'N202': 'tsrc_list_N202',
}

decuplet_ensembles = [
    'A653',
    'B450',
    'H101',
    'H200',
    'U103',
    'N300',
    'N202',
]

baryon_strangeness = {
    'nucleon':  '0',
    'lambda':  '-1',
    'sigma':   '-1',
    'xi':      '-2',
    'decuplet': '-2',
}

baryon_isospin = {
    'nucleon': 'isodoublet',
    'lambda':  'isosinglet',
    'sigma':   'isotriplet',
    'xi':      'isodoublet',
    'decuplet': 'isodoublet',
}

spin_half_irreps = {
    0: 'G1g',
    1: 'G1',
    2: 'G',
    3: 'G',
    4: 'G1',
}

spin_three_half_irreps = {
    0: 'Hg',
}

dibaryon_ops = {
    Channel('P000', 'A1+', 'I0', 'S-2', None) : [
        "isosinglet S=-2 PSQ=0 A1g L(0)L(0)s_S0 0",
        "isosinglet S=-2 PSQ=0 A1g S(0)S(0)s_S0 0",
        "isosinglet S=-2 PSQ=0 A1g N(0)X(0)s_S0 0",
        "isosinglet S=-2 PSQ=0 A1g L(1)L(1)s_S0 0",
        "isosinglet S=-2 PSQ=0 A1g S(1)S(1)s_S0 0",
        "isosinglet S=-2 PSQ=0 A1g N(1)X(1)s_S0 0",
        "isosinglet S=-2 PSQ=0 A1g L(2)L(2)s_S0 0",
        "isosinglet S=-2 PSQ=0 A1g S(2)S(2)s_S0 0",
        "isosinglet S=-2 PSQ=0 A1g N(2)X(2)s_S0 0",
    ],
    Channel('P001', 'A1', 'I0', 'S-2',  None) : [
        "isosinglet S=-2 PSQ=1 A1 L(1)L(0)s_S0 0",
        "isosinglet S=-2 PSQ=1 A1 S(1)S(0)s_S0 0",
        "isosinglet S=-2 PSQ=1 A1 N(1)X(0)s_S0 0",
        "isosinglet S=-2 PSQ=1 A1 N(1)X(0)a_S0 0",
        "isosinglet S=-2 PSQ=1 A1 L(2)L(1)s_S0 0",
        "isosinglet S=-2 PSQ=1 A1 S(2)S(1)s_S0 0",
        "isosinglet S=-2 PSQ=1 A1 N(2)X(1)s_S0 0",
        "isosinglet S=-2 PSQ=1 A1 N(2)X(1)a_S0 0",
        "isosinglet S=-2 PSQ=1 A1 L(2)L(1)s_S1 0",
        "isosinglet S=-2 PSQ=1 A1 S(2)S(1)s_S1 0",
        "isosinglet S=-2 PSQ=1 A1 N(2)X(1)s_S1 0",
        "isosinglet S=-2 PSQ=1 A1 N(2)X(1)a_S1 0",
    ],
    Channel('P011', 'A1', 'I0', 'S-2',  None) : [
        "isosinglet S=-2 PSQ=2 A1 L(2)L(0)s_S0 0",
        "isosinglet S=-2 PSQ=2 A1 S(2)S(0)s_S0 0",
        "isosinglet S=-2 PSQ=2 A1 N(2)X(0)s_S0 0",
        "isosinglet S=-2 PSQ=2 A1 N(2)X(0)a_S0 0",
        "isosinglet S=-2 PSQ=2 A1 L(1)L(1)s_S0 0",
        "isosinglet S=-2 PSQ=2 A1 S(1)S(1)s_S0 0",
        "isosinglet S=-2 PSQ=2 A1 N(1)X(1)s_S0 0",
        "isosinglet S=-2 PSQ=2 A1 L(1)L(1)s_S1 0",
        "isosinglet S=-2 PSQ=2 A1 S(1)S(1)s_S1 0",
        "isosinglet S=-2 PSQ=2 A1 N(1)X(1)s_S1 0",
    ],
    Channel('P111', 'A1', 'I0', 'S-2',  None) : [
        "isosinglet S=-2 PSQ=3 A1 L(3)L(0)s_S0 0",
        "isosinglet S=-2 PSQ=3 A1 S(3)S(0)s_S0 0",
        "isosinglet S=-2 PSQ=3 A1 N(3)X(0)s_S0 0",
        "isosinglet S=-2 PSQ=3 A1 N(3)X(0)a_S0 0",
        "isosinglet S=-2 PSQ=3 A1 L(2)L(1)s_S0 0",
        "isosinglet S=-2 PSQ=3 A1 S(2)S(1)s_S0 0",
        "isosinglet S=-2 PSQ=3 A1 N(2)X(1)s_S0 0",
        "isosinglet S=-2 PSQ=3 A1 N(2)X(1)a_S0 0",
        "isosinglet S=-2 PSQ=3 A1 L(2)L(1)s_S1 0",
        "isosinglet S=-2 PSQ=3 A1 S(2)S(1)s_S1 0",
        "isosinglet S=-2 PSQ=3 A1 N(2)X(1)s_S1 0",
        "isosinglet S=-2 PSQ=3 A1 N(2)X(1)a_S1 0",
    ],
    Channel('P002', 'A1', 'I0', 'S-2',  None) : [
        "isosinglet S=-2 PSQ=4 A1 L(1)L(1)s_S0 0",
        "isosinglet S=-2 PSQ=4 A1 S(1)S(1)s_S0 0",
        "isosinglet S=-2 PSQ=4 A1 N(1)X(1)s_S0 0",
    ],



    Channel('P000', 'A1+', 'I0', 'S-2', '1') : [
        "isosinglet S=-2 PSQ=0 A1g L(0)L(0)s_S0 1",
        "isosinglet S=-2 PSQ=0 A1g L(1)L(1)s_S0 1",
        "isosinglet S=-2 PSQ=0 A1g L(2)L(2)s_S0 1",
    ],
    Channel('P000', 'A1+', 'I0', 'S-2', '27') : [
        "isosinglet S=-2 PSQ=0 A1g L(0)L(0)s_S0 27",
        "isosinglet S=-2 PSQ=0 A1g L(1)L(1)s_S0 27",
        "isosinglet S=-2 PSQ=0 A1g L(2)L(2)s_S0 27",
    ],
    Channel('P000', 'A1+', 'I0', 'S-2', '8') : [
        "isosinglet S=-2 PSQ=0 A1g L(0)L(0)s_S0 8",
        "isosinglet S=-2 PSQ=0 A1g L(1)L(1)s_S0 8",
        "isosinglet S=-2 PSQ=0 A1g L(2)L(2)s_S0 8",
    ],
    Channel('P000', 'A2+', 'I0', 'S-2', '10')  : [],
    Channel('P000', 'A2+', 'I0', 'S-2', '10b') : [],
    Channel('P000', 'A2+', 'I0', 'S-2', '8')   : [],
    Channel('P000', 'E+' , 'I0', 'S-2', '1')   : [],
    Channel('P000', 'E+' , 'I0', 'S-2', '10')  : [],
    Channel('P000', 'E+' , 'I0', 'S-2', '10b') : [],
    Channel('P000', 'E+' , 'I0', 'S-2', '27')  : [],
    Channel('P000', 'E+' , 'I0', 'S-2', '8')   : [],
    Channel('P000', 'T1+', 'I0', 'S-2', '10')  : [],
    Channel('P000', 'T1+', 'I0', 'S-2', '10b') : [],
    Channel('P000', 'T1+', 'I0', 'S-2', '8')   : [],
    Channel('P000', 'T2+', 'I0', 'S-2', '1')   : [],
    Channel('P000', 'T2+', 'I0', 'S-2', '10')  : [],
    Channel('P000', 'T2+', 'I0', 'S-2', '10b') : [],
    Channel('P000', 'T2+', 'I0', 'S-2', '27')  : [],
    Channel('P000', 'T2+', 'I0', 'S-2', '8')   : [],
    Channel('P000', 'A1-', 'I0', 'S-2', '1')   : [],
    Channel('P000', 'A1-', 'I0', 'S-2', '27')  : [],
    Channel('P000', 'A1-', 'I0', 'S-2', '8')   : [],
    Channel('P000', 'A2-', 'I0', 'S-2', '1')   : [],
    Channel('P000', 'A2-', 'I0', 'S-2', '27')  : [],
    Channel('P000', 'A2-', 'I0', 'S-2', '8')   : [],
    Channel('P000', 'E-' , 'I0', 'S-2', '1')   : [],
    Channel('P000', 'E-' , 'I0', 'S-2', '27')  : [],
    Channel('P000', 'E-' , 'I0', 'S-2', '8')   : [],
    Channel('P000', 'T1-', 'I0', 'S-2', '1')   : [],
    Channel('P000', 'T1-', 'I0', 'S-2', '10')  : [],
    Channel('P000', 'T1-', 'I0', 'S-2', '10b') : [],
    Channel('P000', 'T1-', 'I0', 'S-2', '27')  : [],
    Channel('P000', 'T1-', 'I0', 'S-2', '8')   : [],
    Channel('P000', 'T2-', 'I0', 'S-2', '1')   : [],
    Channel('P000', 'T2-', 'I0', 'S-2', '10')  : [],
    Channel('P000', 'T2-', 'I0', 'S-2', '10b') : [],
    Channel('P000', 'T2-', 'I0', 'S-2', '27')  : [],
    Channel('P000', 'T2-', 'I0', 'S-2', '8')   : [],
    Channel('P001', 'A1', 'I0', 'S-2',  '1')   : [
        "isosinglet S=-2 PSQ=1 A1 L(1)L(0)s_S0 1",
        "isosinglet S=-2 PSQ=1 A1 L(2)L(1)s_S0 1",
        "isosinglet S=-2 PSQ=1 A1 L(2)L(1)s_S1 1",
    ],
    Channel('P001', 'A1', 'I0', 'S-2',  '10') : [],
    Channel('P001', 'A1', 'I0', 'S-2',  '10b') : [],
    Channel('P001', 'A1', 'I0', 'S-2',  '27')  : [
        "isosinglet S=-2 PSQ=1 A1 L(1)L(0)s_S0 27",
        "isosinglet S=-2 PSQ=1 A1 L(2)L(1)s_S0 27",
        "isosinglet S=-2 PSQ=1 A1 L(2)L(1)s_S1 27",
    ],
    Channel('P001', 'A1', 'I0', 'S-2',  '8')   : [
        "isosinglet S=-2 PSQ=1 A1 L(1)L(0)s_S0 8",
        "isosinglet S=-2 PSQ=1 A1 L(1)L(0)a_S0 8",
        "isosinglet S=-2 PSQ=1 A1 L(2)L(1)s_S0 8",
        "isosinglet S=-2 PSQ=1 A1 L(2)L(1)a_S0 8",
        "isosinglet S=-2 PSQ=1 A1 L(2)L(1)s_S1 8",
        "isosinglet S=-2 PSQ=1 A1 L(2)L(1)a_S1 8",
    ],
    Channel('P001', 'A2', 'I0', 'S-2',  '1')   : [],
    Channel('P001', 'A2', 'I0', 'S-2',  '10',) : [],
    Channel('P001', 'A2', 'I0', 'S-2',  '10b') : [],
    Channel('P001', 'A2', 'I0', 'S-2',  '27')  : [],
    Channel('P001', 'A2', 'I0', 'S-2',  '8')   : [],
    Channel('P001', 'B1', 'I0', 'S-2',  '1')   : [],
    Channel('P001', 'B1', 'I0', 'S-2',  '10',) : [],
    Channel('P001', 'B1', 'I0', 'S-2',  '10b') : [],
    Channel('P001', 'B1', 'I0', 'S-2',  '27')  : [],
    Channel('P001', 'B1', 'I0', 'S-2',  '8')   : [],
    Channel('P001', 'B2', 'I0', 'S-2',  '1')   : [],
    Channel('P001', 'B2', 'I0', 'S-2',  '10',) : [],
    Channel('P001', 'B2', 'I0', 'S-2',  '10b') : [],
    Channel('P001', 'B2', 'I0', 'S-2',  '27')  : [],
    Channel('P001', 'B2', 'I0', 'S-2',  '8')   : [],
    Channel('P001', 'E', 'I0', 'S-2',   '1')   : [],
    Channel('P001', 'E', 'I0', 'S-2',   '10',) : [],
    Channel('P001', 'E', 'I0', 'S-2',   '10b') : [],
    Channel('P001', 'E', 'I0', 'S-2',   '27')  : [],
    Channel('P001', 'E', 'I0', 'S-2',   '8')   : [],
    Channel('P011', 'A1', 'I0', 'S-2',  '1')   : [
        "isosinglet S=-2 PSQ=2 A1 L(2)L(0)s_S0 1",
        "isosinglet S=-2 PSQ=2 A1 L(1)L(1)s_S0 1",
        "isosinglet S=-2 PSQ=2 A1 L(1)L(1)s_S1 1",
    ],
    Channel('P011', 'A1', 'I0', 'S-2',  '10') : [],
    Channel('P011', 'A1', 'I0', 'S-2',  '10b') : [],
    Channel('P011', 'A1', 'I0', 'S-2',  '27')  : [
        "isosinglet S=-2 PSQ=2 A1 L(2)L(0)s_S0 27",
        "isosinglet S=-2 PSQ=2 A1 L(1)L(1)s_S0 27",
        "isosinglet S=-2 PSQ=2 A1 L(1)L(1)s_S1 27",
    ],
    Channel('P011', 'A1', 'I0', 'S-2',  '8')   : [
        "isosinglet S=-2 PSQ=2 A1 L(2)L(0)s_S0 8",
        "isosinglet S=-2 PSQ=2 A1 L(2)L(0)a_S0 8",
        "isosinglet S=-2 PSQ=2 A1 L(1)L(1)s_S0 8",
        "isosinglet S=-2 PSQ=2 A1 L(1)L(1)s_S1 8",
    ],
    Channel('P011', 'A2', 'I0', 'S-2',  '1')   : [],
    Channel('P011', 'A2', 'I0', 'S-2',  '10',) : [],
    Channel('P011', 'A2', 'I0', 'S-2',  '10b') : [],
    Channel('P011', 'A2', 'I0', 'S-2',  '27')  : [],
    Channel('P011', 'A2', 'I0', 'S-2',  '8')   : [],
    Channel('P011', 'B2', 'I0', 'S-2',  '1')   : [],
    Channel('P011', 'B2', 'I0', 'S-2',  '10',) : [],
    Channel('P011', 'B2', 'I0', 'S-2',  '10b') : [],
    Channel('P011', 'B2', 'I0', 'S-2',  '27')  : [],
    Channel('P011', 'B2', 'I0', 'S-2',  '8')   : [],
    Channel('P011', 'B1', 'I0', 'S-2',  '1')   : [],
    Channel('P011', 'B1', 'I0', 'S-2',  '10',) : [],
    Channel('P011', 'B1', 'I0', 'S-2',  '10b') : [],
    Channel('P011', 'B1', 'I0', 'S-2',  '27')  : [],
    Channel('P011', 'B1', 'I0', 'S-2',  '8')   : [],
    Channel('P111', 'A1', 'I0', 'S-2',  '1')   : [
        "isosinglet S=-2 PSQ=3 A1 L(3)L(0)s_S0 1",
        "isosinglet S=-2 PSQ=3 A1 L(2)L(1)s_S0 1",
        "isosinglet S=-2 PSQ=3 A1 L(2)L(1)s_S1 1",
    ],
    Channel('P111', 'A1', 'I0', 'S-2',  '10') : [],
    Channel('P111', 'A1', 'I0', 'S-2',  '10b') : [],
    Channel('P111', 'A1', 'I0', 'S-2',  '27')  : [
        "isosinglet S=-2 PSQ=3 A1 L(3)L(0)s_S0 27",
        "isosinglet S=-2 PSQ=3 A1 L(2)L(1)s_S0 27",
        "isosinglet S=-2 PSQ=3 A1 L(2)L(1)s_S1 27",
    ],
    Channel('P111', 'A1', 'I0', 'S-2',  '8')   : [
        "isosinglet S=-2 PSQ=3 A1 L(3)L(0)s_S0 8",
        "isosinglet S=-2 PSQ=3 A1 L(3)L(0)a_S0 8",
        "isosinglet S=-2 PSQ=3 A1 L(2)L(1)s_S0 8",
        "isosinglet S=-2 PSQ=3 A1 L(2)L(1)a_S0 8",
        "isosinglet S=-2 PSQ=3 A1 L(2)L(1)s_S1 8",
        "isosinglet S=-2 PSQ=3 A1 L(2)L(1)a_S1 8",
    ],
    Channel('P111', 'A2', 'I0', 'S-2',  '1')   : [],
    Channel('P111', 'A2', 'I0', 'S-2',  '10',) : [],
    Channel('P111', 'A2', 'I0', 'S-2',  '10b') : [],
    Channel('P111', 'A2', 'I0', 'S-2',  '27')  : [],
    Channel('P111', 'A2', 'I0', 'S-2',  '8')   : [],
    Channel('P111', 'E', 'I0', 'S-2',   '1')   : [],
    Channel('P111', 'E', 'I0', 'S-2',   '10',) : [],
    Channel('P111', 'E', 'I0', 'S-2',   '10b') : [],
    Channel('P111', 'E', 'I0', 'S-2',   '27')  : [],
    Channel('P111', 'E', 'I0', 'S-2',   '8')   : [],
    Channel('P002', 'A1', 'I0', 'S-2',  '1')   : [
        "isosinglet S=-2 PSQ=4 A1 L(1)L(1)s_S0 1",
    ],
    Channel('P002', 'A1', 'I0', 'S-2',  '27')  : [
        "isosinglet S=-2 PSQ=4 A1 L(1)L(1)s_S0 27",
    ],
    Channel('P002', 'A1', 'I0', 'S-2',  '8')   : [
        "isosinglet S=-2 PSQ=4 A1 L(1)L(1)s_S0 8",
    ],
    Channel('P002', 'A2', 'I0', 'S-2',  '10',) : [],
    Channel('P002', 'A2', 'I0', 'S-2',  '10b') : [],
    Channel('P002', 'A2', 'I0', 'S-2',  '8')   : [],
    Channel('P002', 'E', 'I0', 'S-2',   '10',) : [],
    Channel('P002', 'E', 'I0', 'S-2',   '10b') : [],
    Channel('P002', 'E', 'I0', 'S-2',   '8')   : [],
}

