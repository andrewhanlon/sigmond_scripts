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

phases = dict()

ensembles = ["cls21_s64_t128_D200"]

omissions = {
    "cls21_s64_t128_D200": set(range(529, 2001, 1)),
}

configs = {
    "cls21_s64_t128_D200": list(range(1, 530, 1)),
}

spatial_extent = {
    "cls21_s64_t128_D200": 64,
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
