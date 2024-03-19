import os

from typing import NamedTuple
import regex
import enum

class Ensemble(NamedTuple):
  name: str
  dir_name: str
  su3: bool
  openBC: bool
  Nt: int
  replica: tuple
  replica_str: str

class Channel(NamedTuple):
  flavor: str
  two_times_I3: int
  irrep: str
  irrep_row: int

  def __repr__(self):
    return f"{self.flavor}_{self.two_times_I3}-{self.irrep}_{self.irrep_row}"

  def __str__(self):
    irrepz_str = str(int(self.two_times_I3/2)) if self.two_times_I3 % 2 == 0 else str(self.two_times_I3) + 'h'
    return f"{short_name[self.flavor]}_Iz{irrepz_str}"

class DataInfo(NamedTuple):
  operator_ids: list
  contributions: dict

def get_source_time(source):
  source = source.decode()
  return int(source[source.rfind(' ')+1:-1])

ensembles = [
    Ensemble("A654", "A654", False, False, 48, ('r000',), 'r000'),
]

omissions = {
    'A654': set(range(5068)) - set(range(0, 5068, 8)),
}

configs = {
    "A654": list(range(1, 5068, 8))
}

class ParticleType(enum.Enum):
  BOSON = enum.auto()
  FERMION = enum.auto()

particle_types = {
    'etaprime': ParticleType.BOSON,
    'eta_pion': ParticleType.BOSON,
    'pion_etaprime': ParticleType.BOSON,
    'etaprime_pion': ParticleType.BOSON,
    'eta_etaprime': ParticleType.BOSON,
    'pion_eta': ParticleType.BOSON,
    'etaprime_eta': ParticleType.BOSON,
    'pion': ParticleType.BOSON,
    'eta': ParticleType.BOSON,
    'kaon_minus': ParticleType.BOSON,
    'kaon_plus': ParticleType.BOSON,

    'omega': ParticleType.FERMION,
    'sigma_lambdaprime': ParticleType.FERMION,
    'lambda_sigma': ParticleType.FERMION,
    'lambdaprime': ParticleType.FERMION,
    'xi': ParticleType.FERMION,
    'lambdaprime_lambda': ParticleType.FERMION,
    'delta': ParticleType.FERMION,
    'sigma_lambda': ParticleType.FERMION,
    'nucleon': ParticleType.FERMION,
    'lambda': ParticleType.FERMION,
    'lambda_lambdaprime': ParticleType.FERMION,
    'sigma': ParticleType.FERMION,
    'lambdaprime_sigma': ParticleType.FERMION,
}

irrep_rows = {
    'G1g': {1:1, -1:2},
    'G1u': {1:1, -1:2},
    'Hg': {1:1, 3:2, -3:3, -1:4},
    'Hu': {1:1, 3:2, -3:3, -1:4},
    'A1u': {1:1},
}

flavors = {
    'etaprime': ('0', '0'),
    'eta_pion': ('0', '0'),
    'pion_etaprime': ('0', '0'),
    'etaprime_pion': ('0', '0'),
    'eta_etaprime': ('0', '0'),
    'pion_eta': ('0', '0'),
    'etaprime_eta': ('0', '0'),
    'pion': ('1', '0'),
    'eta': ('0', '0'),
    'kaon_minus': ('1h', '1'),
    'kaon_plus': ('1h', '1'),

    'omega': ('0', '-3'),
    'sigma_lambdaprime': ('0', '-1'),
    'lambda_sigma': ('0', '-1'),
    'lambdaprime': ('0', '-1'),
    'xi': ('1h', '-2'),
    'lambdaprime_lambda': ('0', '-1'),
    'delta': ('3h', '0'),
    'sigma_lambda': ('0', '-1'),
    'nucleon': ('1h', '0'),
    'lambda': ('0', '-1'),
    'lambda_lambdaprime': ('0', '-1'),
    'sigma': ('1', '-1'),
    'lambdaprime_sigma': ('0', '-1'),
}

short_name = {
    'etaprime': 'ep',
    'eta_pion': 'e_pi',
    'pion_etaprime': 'pi_ep',
    'etaprime_pion': 'ep_pi',
    'eta_etaprime': 'e_ep',
    'pion_eta': 'pi_e',
    'etaprime_eta': 'ep_e',
    'pion': 'pi',
    'eta': 'e',
    'kaon_minus': 'km',
    'kaon_plus': 'kp',

    'omega': 'O',
    'sigma_lambdaprime': 'S_Lp',
    'lambda_sigma': 'L_S',
    'lambdaprime': 'Lp',
    'xi': 'X',
    'lambdaprime_lambda': 'Lp_L',
    'delta': 'D',
    'sigma_lambda': 'S_L',
    'nucleon': 'N',
    'lambda': 'L',
    'lambda_lambdaprime': 'L_Lp',
    'sigma': 'S',
    'lambdaprime_sigma': 'Lp_S',
}

def shorten_diagram(diagram):
  return diagram.replace('mass', 'm').replace('self', 's').replace('vertex', 'v')

def shorten_contribution(contribution):
  return contribution.replace('symmetric', 'symm')

def get_channel_from_str(channel_str):
  pattern = r"^(?P<flavor>\w+)_(?P<two_times_I3>-?\d)-(?P<irrep>\w+)_(?P<irrep_row>-?\d+)$"
  channel_match = regex.match(pattern, channel_str)

  return Channel(channel_match['flavor'], int(channel_match['two_times_I3']),
                 channel_match['irrep'], int(channel_match['irrep_row']))



def raw_data_file(ensemble_name, config, exact):
  if exact:
    return f"{ensemble_name}_ib_two_point_exactn{config}.hdf5"
  else:
    return f"{ensemble_name}_ib_two_point_sloppyn{config}.hdf5"

def raw_data_dir(base_data_dir, ensemble_name):
  return os.path.join(base_data_dir, ensemble_name, 'correlators', 'raw')

def ave_data_file(ensemble_name, channel, config):
  return f"{ensemble_name}_ib_two_point_{channel!r}_n{config}.hdf5"

def ave_data_dir(base_data_dir, ensemble_name):
  return os.path.join(base_data_dir, ensemble_name, 'correlators', 'ave')

def merged_data_file(ensemble_name, channel):
  return f"{ensemble_name}_ib_two_point_{channel!r}.hdf5"

def merged_data_dir(base_data_dir, ensemble_name):
  return os.path.join(base_data_dir, ensemble_name, 'correlators', 'merged')

def sigmond_data_file(ensemble_name, channel, contribution):
  return f"{ensemble_name}_ib_two_point_{channel!r}_{contribution}.bin"

def sigmond_data_dir(base_data_dir, ensemble_name):
  return os.path.join(base_data_dir, ensemble_name)
