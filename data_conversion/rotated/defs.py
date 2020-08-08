from typing import NamedTuple

class Channel(NamedTuple):
  psq: tuple
  irrep: str
  isospin: str
  strangeness: int

  def __repr__(self):
    strangeness_str = f"S{self.strangeness}".replace('-', 'm')

    return f"PSQ{self.psq}_{self.irrep}_{self.isospin}_{strangeness_str}"

output_dir = "data"

input_dirs = {
    "cls21_s64_t128_D200": "/media/ext1/research/projects/D200_project/analysis/",
    "cls21_c103": "/media/ext1/research/projects/C103_project/analysis/",
}

phases = dict()

omissions = {
    "cls21_s64_t128_D200": set(range(529, 2001, 1)),
    "cls21_c103": set(),
}

name_map = {
    'pion': 'pi',
    'eta': 'e',
    'phi': 'p',
    'kaon': 'k',
    'kbar': 'kb',
    'nucleon': 'N',
    'delta': 'D',
    'sigma': 'S',
    'lambda': 'L',
    'xi': 'X',
    'omega': 'O'
}
