from collections import namedtuple

Ensemble = namedtuple('Ensemble', ['name', 'dir_name', 'type', 'su3', 'open', 'flow', 'Nt', 'replica', 'replica_fwd', 'replica_str', 'modes', 'srcs', 't0', 'ts'])
Channel = namedtuple('Channel', ['P', 'irrep', 'flavor'])

# ensemble info
ensembles = [
    Ensemble("clover_T1.0_cls_A653", "A653", 'clover_T1.0_cls', True, False, 'T1.0_', 48, ('r000',), (True,), 'r000', 32, 8, 0, 24),
    Ensemble("clover_T1.0_cls_B450", "B450", 'clover_T1.0_cls', True, False, 'T1.0_', 64, ('r000',), (True,), 'r000', 32, 8, 0, 32),
    Ensemble("clover_T0.125_cls_A653", "A653", 'clover_T0.125_cls', True, False, 'T0.125_', 48, ('r000',), (True,), 'r000', 32, 8, 0, 24),
    Ensemble("mdwf_T1.0_cls_A653", "A653", 'mdwf_T1.0_cls', True, False, '', 48, ('r000',), (True,), 'r000', 32, 8, 0, 24),
    Ensemble("mdwf_T1.0_cls_B450", "B450", 'mdwf_T1.0_cls', True, False, '', 64, ('r000',), (True,), 'r000', 32, 8, 0, 32),
    Ensemble("mdwf_T1.0_cls_H101", "H101", 'mdwf_T1.0_cls', True, True, '', 96, ('r000','r001',), (True, True,), 'r000-001', 48, 4, 0, 24),
    Ensemble("mdwf_T1.0_cls_H200", "H200", 'mdwf_T1.0_cls', True, True, '', 96, ('r000','r001'), (True, True,), 'r000-001', 20, 8, 0, 32),
    Ensemble("mdwf_T0.4_cls_A653", "A653", 'mdwf_T0.4_cls', True, False, '', 48, ('r000',), (True,), 'r000', 32, 8, 0, 24),
    Ensemble("mdwf_T1.0_hisq_a15m400trMc_a", "a15m400trMc_a", 'mdwf_hisq', True, False, '', 48, ('',), (True,), '', 60, 8, 0, 16),
    Ensemble("mdwf_T1.0_hisq_a12m400trMc_a", "a12m400trMc_a", 'mdwf_hisq', True, False, '', 64, ('',), (True,), '', 54, 8, 0, 20),
    Ensemble("mdwf_T1.0_hisq_a09m400trMc_a", "a09m400trMc_a", 'mdwf_hisq', True, False, '', 96, ('',), (True,), '', 54, 8, 0, 28),
    Ensemble("mdwf_T1.0_hisq_a06m400trMc_a", "a06m400trMc_a", 'mdwf_hisq', True, False, '', 144, ('',), (True,), '', 46, 8, 0, 44),
    Ensemble("exp_a12_m400_mL6.0_trMc", "a12_m400_mL6.0_trMc", 'exp', True, False, '', 96, ('',), (True,), '', 54, 8, 0, 20),
    Ensemble("exp_a094_m400_mL6.0_trMc", "a094_m400_mL6.0_trMc", 'exp', True, False, '', 96, ('_ext',), (True,), 'ext', 64, 8, 0, 24),
    Ensemble("exp_a064_m400_mL6.4_trMc", "a064_m400_mL6.4_trMc", 'exp', True, False, '', 96, ('_1','_2',), (True, True), '1-2', 68, 8, 0, 40),
    Ensemble("cls_A653", "A653", 'cls', True, False, '', 48, ('r000',), (True,), 'r000', 32, 4, 0, 24),
    Ensemble("cls_B450", "B450", 'cls', True, False, '', 64, ('r000',), (True,), 'r000', 32, 8, 0, 32),
    Ensemble("cls_H101", "H101", 'cls', True, True, '', 96, ('r000','r001',), (True, True), 'r000-001', 48, 4, 0, 24),
    Ensemble("cls_H200", "H200", 'cls', True, True, '', 96, ('r000','r001'), (True, True), 'r000-001', 20, 8, 0, 32),
    Ensemble("cls_J500", "J500", 'cls', True, True, '', 192, ('r004','r005'), (True, True), 'r004-005', 36, 12, 0, 56),
    Ensemble("cls_J500_sm2", "J500_sm2", 'cls', True, True, '', 192, ('r004','r005','r006'), (True, True, True), 'r004-006', 36, 12, 0, 56),
    Ensemble("cls_J307", "J307", 'cls', True, True, '', 192, ('r001','r000'), (True, False), 'r001-000b', 76, 8, 0, 40),
    Ensemble("cls_N202", "N202", 'cls', True, True, '', 128, ('r001',), (True,), 'r001', 68, 8, 0, 32),
    Ensemble("cls_N300", "N300", 'cls', True, True, '', 128, ('r001','r002'), (True, True), 'r001-002', 32, 12, 0, 40),
    Ensemble("cls_U103", "U103", 'cls', True, True, '', 128, ('r001','r002','r003',), (True, True, True), 'r001-003', 20, 5, 0, 24),
    #Ensemble("cls_E1", "E1", 'cls', True, False, '', 64, ('',), '', 30, 8, 0, 32),
    #Ensemble("cls_E5", "E5_SU3", 'cls', True, False, '', 64, ('f', 'g',), 'fg', 30, 4, 0, 32),
    #Ensemble("cls_B451", "B451", 'cls', False, False, '', 64, ('r000',), 'r000', 32, 4, 0, 32),
    #Ensemble("cls_B452", "B452", 'cls', False, False, '', 64, ('r000',), 'r000', 32, 4, 0, 32),
    #Ensemble("cls_H102", "H102", 'cls', False, True, '', 96, ('r001','r002',), 'r001-002', 48, 4, 0, 24),
    #Ensemble("cls_H107", "H107", 'cls', False, True, '', 96, ('r000','r001','r002','r003','r004','r005'), 'r000-005', 48, 4, 0, 24),
    #Ensemble("cls_N200", "N200", 'cls', False, True, '', 128, ('r000','r001'), 'r000-001', 68, 8, 0, 32),
    #Ensemble("cls_N451", "N451", 'cls', False, False, '', 128, ('r000',), 'r000', 108, 8, 0, 32),
    #Ensemble("cls_U102", "U102", 'cls', False, True, '', 128, ('r001', 'r002',), 'r001-002', 20, 5, 0, 24),
]

backward_prop_skip = {
    'clover_T1.0_cls_A653': [],
    'clover_T1.0_cls_B450': [],
    'clover_T0.125_cls_A653': [],
    'mdwf_T1.0_cls_A653': [],
    'mdwf_T1.0_cls_B450': [],
    'mdwf_T1.0_cls_H101': [0],
    'mdwf_T1.0_cls_H200': [0,1,2,3],
    'mdwf_T0.4_cls_A653': [],
    'mdwf_T1.0_hisq_a15m400trMc_a': [],
    'mdwf_T1.0_hisq_a12m400trMc_a': [],
    'mdwf_T1.0_hisq_a09m400trMc_a': [],
    'mdwf_T1.0_hisq_a06m400trMc_a': [],
    'exp_a12_m400_mL6.0_trMc': [],
    'exp_a094_m400_mL6.0_trMc': [],
    'exp_a064_m400_mL6.4_trMc': [],
    'cls_A653': [],
    'cls_B450': [],
    'cls_H101': [0],
    'cls_H200': [0,1,2,3],
    'cls_J307': [0],
    'cls_J500': [0,1,2],
    'cls_J500_sm2': [0,1,2,3,4,5],
    'cls_N202': [0,1],
    'cls_N300': [0,1,2,3,4,5,6,7],
    'cls_U103': [0],
    'cls_E1':   [],
    'cls_E5':   [],
    'cls_B451': [],
    'cls_B452': [],
    'cls_H102': [0],
    'cls_H107': [0],
    'cls_N200': [0,1],
    'cls_N451': [],
    'cls_U102': [0],
}

forward_prop_skip = {
    'clover_T1.0_cls_A653': [],
    'clover_T1.0_cls_B450': [],
    'clover_T0.125_cls_A653': [],
    'mdwf_T1.0_cls_A653': [],
    'mdwf_T1.0_cls_B450': [],
    'mdwf_T1.0_cls_H101': [3],
    'mdwf_T1.0_cls_H200': [4,5,6,7],
    'mdwf_T0.4_cls_A653': [],
    'mdwf_T1.0_hisq_a15m400trMc_a': [],
    'mdwf_T1.0_hisq_a12m400trMc_a': [],
    'mdwf_T1.0_hisq_a09m400trMc_a': [],
    'mdwf_T1.0_hisq_a06m400trMc_a': [],
    'exp_a12_m400_mL6.0_trMc': [],
    'exp_a094_m400_mL6.0_trMc': [],
    'exp_a064_m400_mL6.4_trMc': [],
    'cls_A653': [],
    'cls_B450': [],
    'cls_H101': [3],
    'cls_H200': [4,5,6,7],
    'cls_J307': [7],
    'cls_J500': [9,10,11],
    'cls_J500_sm2': [6,7,8,9,10,11],
    'cls_N202': [6,7],
    'cls_N300': [4,5,6,7,8,9,10,11],
    'cls_U103': [4],
    'cls_E1':   [],
    'cls_E5':   [],
    'cls_B451': [],
    'cls_B452': [],
    'cls_H102': [3],
    'cls_H107': [3],
    'cls_N200': [6,7],
    'cls_N451': [],
    'cls_U102': [4],
}

pseudoscalar_backward_prop_skip = {
    'clover_T1.0_cls_A653': [],
    'clover_T1.0_cls_B450': [],
    'clover_T0.125_cls_A653': [],
    'mdwf_T1.0_cls_A653': [],
    'mdwf_T1.0_cls_B450': [],
    'mdwf_T1.0_cls_H101': [0,1],
    'mdwf_T1.0_cls_H200': [0,1,2,3,4,5],
    'mdwf_T0.4_cls_A653': [],
    'mdwf_T1.0_hisq_a15m400trMc_a': [],
    'mdwf_T1.0_hisq_a12m400trMc_a': [],
    'mdwf_T1.0_hisq_a09m400trMc_a': [],
    'mdwf_T1.0_hisq_a06m400trMc_a': [],
    'exp_a12_m400_mL6.0_trMc': [],
    'exp_a094_m400_mL6.0_trMc': [],
    'exp_a064_m400_mL6.4_trMc': [],
    'cls_A653': [],
    'cls_B450': [],
    'cls_H101': [0,1],
    'cls_H200': [0,1,2,3,4,5],
    'cls_J307': [0,1,2],
    'cls_J500': [0,1,2,3,4,5],
    'cls_J500_sm2': [0,1,2,3,4,5],
    'cls_N202': [0,1,2,3],
    'cls_N300': [0,1,2,3,4,5,6,7],
    'cls_U103': [0,1],
    'cls_E1':   [],
    'cls_E5':   [],
    'cls_B451': [],
    'cls_B452': [],
    'cls_H102': [0,1],
    'cls_H107': [0,1],
    'cls_N200': [0,1,2],
    'cls_N451': [],
    'cls_U102': [0,1],
}

pseudoscalar_forward_prop_skip = {
    'clover_T1.0_cls_A653': [],
    'clover_T1.0_cls_B450': [],
    'clover_T0.125_cls_A653': [],
    'mdwf_T1.0_cls_A653': [],
    'mdwf_T1.0_cls_B450': [],
    'mdwf_T1.0_cls_H101': [2,3],
    'mdwf_T1.0_cls_H200': [2,3,4,5,6,7],
    'mdwf_T0.4_cls_A653': [],
    'mdwf_T1.0_hisq_a15m400trMc_a': [],
    'mdwf_T1.0_hisq_a12m400trMc_a': [],
    'mdwf_T1.0_hisq_a09m400trMc_a': [],
    'mdwf_T1.0_hisq_a06m400trMc_a': [],
    'exp_a12_m400_mL6.0_trMc': [],
    'exp_a094_m400_mL6.0_trMc': [],
    'exp_a064_m400_mL6.4_trMc': [],
    'cls_A653': [],
    'cls_B450': [],
    'cls_H101': [2,3],
    'cls_H200': [2,3,4,5,6,7],
    'cls_J307': [5,6,7],
    'cls_J500': [6,7,8,9,10,11],
    'cls_J500_sm2': [6,7,8,9,10,11],
    'cls_N202': [4,5,6,7],
    'cls_N300': [4,5,6,7,8,9,10,11],
    'cls_U103': [3,4],
    'cls_E1':   [],
    'cls_E5':   [],
    'cls_B451': [],
    'cls_B452': [],
    'cls_H102': [2,3],
    'cls_H107': [2,3],
    'cls_N200': [5,6,7],
    'cls_N451': [],
    'cls_U102': [3,4],
}


pseudoscalar_modes = {
    'clover_T1.0_cls_A653': 32,
    'clover_T1.0_cls_B450': 32,
    'clover_T0.125_cls_A653': 32,
    'mdwf_T1.0_cls_A653': 32,
    'mdwf_T1.0_cls_B450': 32,
    'mdwf_T1.0_cls_H101': 48,
    'mdwf_T1.0_cls_H200': 20,
    'mdwf_T0.4_cls_A653': 32,
    'mdwf_T1.0_hisq_a15m400trMc_a': 60,
    'mdwf_T1.0_hisq_a12m400trMc_a': 54,
    'mdwf_T1.0_hisq_a09m400trMc_a': 54,
    'mdwf_T1.0_hisq_a06m400trMc_a': 46,
    'exp_a12_m400_mL6.0_trMc': 54,
    'exp_a094_m400_mL6.0_trMc': 64,
    'exp_a064_m400_mL6.4_trMc': 68,
    'cls_A653': 32,
    'cls_B450': 32,
    'cls_H101': 144,
    'cls_H200': 20,
    'cls_J307': 76,
    'cls_J500': 36,
    'cls_J500_sm2': 36,
    'cls_N202': 68,
    'cls_N300': 32,
    'cls_U103': 20,
    'cls_E1':   56,
    'cls_E5':   60,
    'cls_B451': 32,
    'cls_B452': 32,
    'cls_H102': 72,
    'cls_H107': 48,
    'cls_N200': 68,
    'cls_N451': 108,
    'cls_U102': 20,
}

pseudoscalar_sources = {
    'clover_T1.0_cls_A653': [0, 0, 0, 0, 0, 0, 0, 0],
    'clover_T1.0_cls_B450': [0, 0, 0, 0, 0, 0, 0, 0],
    'clover_T0.125_cls_A653': [0, 0, 0, 0, 0, 0, 0, 0],
    'mdwf_T1.0_cls_A653': [0, 0, 0, 0, 0, 0, 0, 0],
    'mdwf_T1.0_cls_B450': [0, 0, 0, 0, 0, 0, 0, 0],
    'mdwf_T1.0_cls_H101': [0, 0, 0, 0],
    'mdwf_T1.0_cls_H200': [0, 0, 0, 0, 0, 0, 0, 0],
    'mdwf_T0.4_cls_A653': [0, 0, 0, 0, 0, 0, 0, 0],
    'mdwf_T1.0_hisq_a15m400trMc_a': [0, 0, 0, 0, 0, 0, 0, 0],
    'mdwf_T1.0_hisq_a12m400trMc_a': [0, 0, 0, 0, 0, 0, 0, 0],
    'mdwf_T1.0_hisq_a09m400trMc_a': [0, 0, 0, 0, 0, 0, 0, 0],
    'mdwf_T1.0_hisq_a06m400trMc_a': [0, 0, 0, 0, 0, 0, 0, 0],
    'exp_a12_m400_mL6.0_trMc': [0, 0, 0, 0, 0, 0, 0, 0],
    'exp_a094_m400_mL6.0_trMc': [0, 0, 0, 0, 0, 0, 0, 0],
    'exp_a064_m400_mL6.4_trMc': [0, 0, 0, 0, 0, 0, 0, 0],
    'cls_A653': [0, 0, 0, 0],
    'cls_B450': [0, 0, 0, 0, 0, 0, 0, 0],
    'cls_H101': [24, 40, 55, 71],
    'cls_H200': [0, 4, 8, 12, 16, 20, 24, 28],
    'cls_J307': [0, 10, 20, 30, 40, 50, 60, 70],
    'cls_J500': [0, 7, 14, 21, 28, 35, 42, 49, 56, 63, 70, 77],
    'cls_J500_sm2': [0, 7, 14, 21, 28, 35, 42, 49, 56, 63, 70, 77],
    'cls_N202': [0, 6, 12, 18, 24, 30, 36, 42],
    'cls_N300': [0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44],
    'cls_U103': [32, 48, 64, 79, 95],
    'cls_E1':   [0, 0, 0, 0, 0, 0, 0, 0],
    'cls_E5':   [0, 0, 0, 0],
    'cls_B451': [0, 0, 0, 0],
    'cls_B452': [0, 0, 0, 0],
    'cls_H102': [24, 40, 55, 71],
    'cls_H107': [24, 40, 55, 71],
    'cls_N200': [0, 6, 12, 18, 24, 30, 36, 42],
    'cls_N451': [0, 0, 0, 0, 0, 0, 0, 0],
    'cls_U102': [32, 48, 64, 79, 95],
}

pseudoscalar_names = {
    'clover_T1.0_cls_A653': ['ps'],
    'clover_T1.0_cls_B450': ['ps'],
    'clover_T0.125_cls_A653': ['ps'],
    'mdwf_T1.0_cls_A653': ['ps'],
    'mdwf_T1.0_cls_B450': ['ps'],
    'mdwf_T1.0_cls_H101': ['ps'],
    'mdwf_T1.0_cls_H200': ['ps'],
    'mdwf_T0.4_cls_A653': ['ps'],
    'mdwf_T1.0_hisq_a15m400trMc_a': ['ps'],
    'mdwf_T1.0_hisq_a12m400trMc_a': ['ps'],
    'mdwf_T1.0_hisq_a09m400trMc_a': ['ps'],
    'mdwf_T1.0_hisq_a06m400trMc_a': ['ps'],
    'exp_a12_m400_mL6.0_trMc': ['ps'],
    'exp_a094_m400_mL6.0_trMc': ['ps'],
    'exp_a064_m400_mL6.4_trMc': ['ps'],
    'cls_A653': ['ps'],
    'cls_B450': ['ps'],
    'cls_H101': ['ps'],
    'cls_H200': ['ps'],
    'cls_J307': ['pion'],
    'cls_J500': ['ps'],
    'cls_J500_sm2': ['ps'],
    'cls_N202': ['ps'],
    'cls_N300': ['ps'],
    'cls_U103': ['ps'],
    'cls_E1':   ['ps'],
    'cls_E5':   ['pion'],
    'cls_B451': ['pion', 'kaon', 'eta_s'],
    'cls_B452': ['pion', 'kaon', 'eta_s'],
    'cls_H102': ['pion', 'kaon', 'eta_s'],
    'cls_H107': ['pion', 'kaon', 'eta_s'],
    'cls_N200': ['pion', 'kaon', 'eta_s'],
    'cls_N451': ['pion', 'kaon', 'eta_s'],
    'cls_U102': ['pion', 'kaon', 'eta_s'],
}

pseudoscalar_op_strs = {
    'clover_T1.0_cls_A653': {
      'ps': "Flavor=1,0 Pref=(0,0,0) A1um pion 0",
    },
    'clover_T1.0_cls_B450': {
      'ps': "Flavor=1,0 Pref=(0,0,0) A1um pion 0",
    },
    'clover_T0.125_cls_A653': {
      'ps': "Flavor=1,0 Pref=(0,0,0) A1um pion 0",
    },
    'mdwf_T1.0_cls_A653': {
      'ps': "Flavor=1,0 Pref=(0,0,0) A1um pion 0",
    },
    'mdwf_T1.0_cls_B450': {
      'ps': "Flavor=1,0 Pref=(0,0,0) A1um pion 0",
    },
    'mdwf_T1.0_cls_H101': {
      'ps': "Flavor=1,0 Pref=(0,0,0) A1um pion 0",
    },
    'mdwf_T1.0_cls_H200': {
      'ps': "Flavor=1,0 Pref=(0,0,0) A1um pion 0",
    },
    'mdwf_T0.4_cls_A653': {
      'ps': "Flavor=1,0 Pref=(0,0,0) A1um pion 0",
    },
    'mdwf_T1.0_hisq_a15m400trMc_a': {
      'ps': "Flavor=1,0 Pref=(0,0,0) A1um pion 0",
    },
    'mdwf_T1.0_hisq_a12m400trMc_a': {
      'ps': "Flavor=1,0 Pref=(0,0,0) A1um pion 0",
    },
    'mdwf_T1.0_hisq_a09m400trMc_a': {
      'ps': "Flavor=1,0 Pref=(0,0,0) A1um pion 0",
    },
    'mdwf_T1.0_hisq_a06m400trMc_a': {
      'ps': "Flavor=1,0 Pref=(0,0,0) A1um pion 0",
    },
    'exp_a12_m400_mL6.0_trMc': {
      'ps': "Flavor=1,0 Pref=(0,0,0) A1um pion 0",
    },
    'exp_a094_m400_mL6.0_trMc': {
      'ps': "Flavor=1,0 Pref=(0,0,0) A1um pion 0",
    },
    'exp_a064_m400_mL6.4_trMc': {
      'ps': "Flavor=1,0 Pref=(0,0,0) A1um pion 0",
    },
    'cls_A653': {
      'ps': "Flavor=1,0 Pref=(0,0,0) A1um pion 0",
    },
    'cls_B450': {
      'ps': "Flavor=1,0 Pref=(0,0,0) A1um pion 0",
    },
    'cls_H101': {
      'ps': "Flavor=1,0 Pref=(0,0,0) A1um pion 0",
    },
    'cls_H200': {
      'ps': "Flavor=1,0 Pref=(0,0,0) A1um pion 0",
    },
    'cls_J307': {
      'pion': "Flavor=1,0 Pref=(0,0,0) A1um pion 0",
    },
    'cls_J500': {
      'ps': "Flavor=1,0 Pref=(0,0,0) A1um pion 0",
    },
    'cls_J500_sm2': {
      'ps': "Flavor=1,0 Pref=(0,0,0) A1um pion 0",
    },
    'cls_N202': {
      'ps': "Flavor=1,0 Pref=(0,0,0) A1um pion 0",
    },
    'cls_N300': {
      'ps': "Flavor=1,0 Pref=(0,0,0) A1um pion 0",
    },
    'cls_U103': {
      'ps': "Flavor=1,0 Pref=(0,0,0) A1um pion 0",
    },
    'cls_E1': {
      'ps': "Flavor=1,0 Pref=(0,0,0) A1um pion 0",
    },
    'cls_E5': {
      'pion': "Flavor=1,0 Pref=(0,0,0) A1um pion 0",
    },
    'cls_B451': {
      'pion': "Flavor=1,0 Pref=(0,0,0) A1um pion 0",
      'kaon': "Flavor=1h,1 Pref=(0,0,0) A1u kaon 0",
      'eta_s': "Flavor=0,0 Pref=(0,0,0) A1up eta_s 0",
    },
    'cls_B452': {
      'pion': "Flavor=1,0 Pref=(0,0,0) A1um pion 0",
      'kaon': "Flavor=1h,1 Pref=(0,0,0) A1u kaon 0",
      'eta_s': "Flavor=0,0 Pref=(0,0,0) A1up eta_s 0",
    },
    'cls_H102': {
      'pion': "Flavor=1,0 Pref=(0,0,0) A1um pion 0",
      'kaon': "Flavor=1h,1 Pref=(0,0,0) A1u kaon 0",
      'eta_s': "Flavor=0,0 Pref=(0,0,0) A1up eta_s 0",
    },
    'cls_H107': {
      'pion': "Flavor=1,0 Pref=(0,0,0) A1um pion 0",
      'kaon': "Flavor=1h,1 Pref=(0,0,0) A1u kaon 0",
      'eta_s': "Flavor=0,0 Pref=(0,0,0) A1up eta_s 0",
    },
    'cls_N200': {
      'pion': "Flavor=1,0 Pref=(0,0,0) A1um pion 0",
      'kaon': "Flavor=1h,1 Pref=(0,0,0) A1u kaon 0",
      'eta_s': "Flavor=0,0 Pref=(0,0,0) A1up eta_s 0",
    },
    'cls_N451': {
      'pion': "Flavor=1,0 Pref=(0,0,0) A1um pion 0",
      'kaon': "Flavor=1h,1 Pref=(0,0,0) A1u kaon 0",
      'eta_s': "Flavor=0,0 Pref=(0,0,0) A1up eta_s 0",
    },
    'cls_U102': {
      'pion': "Flavor=1,0 Pref=(0,0,0) A1um pion 0",
      'kaon': "Flavor=1h,1 Pref=(0,0,0) A1u kaon 0",
      'eta_s': "Flavor=0,0 Pref=(0,0,0) A1up eta_s 0",
    },
}

tsrc_files = {
    'cls_H200': 'tsrc_list_H200',
    'cls_J307': 'tsrc_list_J307',
    'cls_J500': 'tsrc_list_J500',
    'cls_J500_sm2': 'tsrc_list_J500',
    'cls_N202': 'tsrc_list_N202',
    'cls_N300': 'tsrc_list_N300',
    'cls_N200': 'tsrc_list_N200',
}

decuplet_ensembles = [
    'clover_T1.0_cls_A653',
    'clover_T1.0_cls_B450',
    'clover_T0.125_cls_A653',
    'mdwf_T1.0_cls_A653',
    'mdwf_T1.0_cls_B450',
    'mdwf_T1.0_cls_H101',
    'mdwf_T1.0_cls_H200',
    'mdwf_T0.4_cls_A653',
    'mdwf_T1.0_hisq_a15m400trMc_a',
    'mdwf_T1.0_hisq_a12m400trMc_a',
    'mdwf_T1.0_hisq_a09m400trMc_a',
    'mdwf_T1.0_hisq_a06m400trMc_a',
    'exp_a12_m400_mL6.0_trMc',
    'exp_a094_m400_mL6.0_trMc',
    'exp_a064_m400_mL6.4_trMc',
    'cls_A653',
    'cls_B450',
    'cls_H101',
    'cls_H200',
    'cls_J307',
    'cls_J500',
    'cls_J500_sm2',
    'cls_N202',
    'cls_N300',
    'cls_U103',
    'cls_H102',
    'cls_H107',
    'cls_N200',
    'cls_N451',
]

# operator info
SU3_flavors = ["1", "27", "10", "10b", "8"]
SU2_flavors = ["I0_S0", "I1_S0", "I0_S-2", "I1_S-2", "I2_S-2"]
irrep_spin_ordering = "I0_S-2"
symmetric_SU3_flavors = ["1", "27", "8"]
anti_symmetric_SU3_flavors = ["10", "10b", "8"]
symmetric_SU2_flavor = "I1_S0"
anti_symmetric_SU2_flavor = "I0_S0"
symmetric_SU3_flavor_dibaryon = "NNs"
anti_symmetric_SU3_flavor_dibaryon = "NNa"

op_files = {
    'S0': "operators_S0",
    'S-2': "operators_S-2",
}

baryon_flavor = {
    'nucleon':  '1h,0',
    'lambda':   '0,-1',
    'sigma':    '1,-1',
    'xi':       '1h,-2',
    'octet':    '8',
    'decuplet': '10',
}

spin_half_irreps = {
    0: 'G1g',
    1: 'G1',
    2: 'G',
    3: 'G',
    4: 'G1',
    5: 'F1',
}

spin_three_half_irreps = {
    0: 'Hg',
}

total_prefs = {
    0: '(0,0,0)',
    1: '(0,0,1)',
    2: '(0,1,1)',
    3: '(1,1,1)',
    #4: '(0,0,2)',
    #5: '(0,1,2)',
}

def convert_irrep(irrep, psq):
  if psq == 2 and irrep == "B1":
    return "B2"
  elif psq == 2 and irrep == "B2":
    return "B1"
  elif irrep == "E2":
    return "E"

  return irrep.replace("+", "g").replace("-", "u")


def get_opstr(su3, pref, irrep, flavor, op_str):
  psq = int(pref[1])**2 + int(pref[2])**2 + int(pref[3])**2
  pref = f"Pref=({pref[1]},{pref[2]},{pref[3]})"
  irrep = convert_irrep(irrep, psq)
  if not su3:
    isospin, strangeness = flavor.split('_')
    isospin = isospin[1:]
    strangeness = strangeness[1:]
    flavor = f"{isospin},{strangeness}"

  return f"Flavor={flavor} {pref} {irrep} {op_str}"
