from typing import NamedTuple

class AveragedChannel(NamedTuple):
  psq: int
  irrep: str

  def __str__(self):
    return f"(psq={self.psq}, irrep={self.irrep})"

  def __repr__(self):
    return f"pSq{self.psq}-{self.irrep}"

class EquivalentFrame(NamedTuple):
  mom: tuple
  row: int

  def __str__(self):
    return f"(P={self.mom}, row={self.row})"

  def __repr__(self):
    return f"P{''.join(list(map(str, self.mom)))}_{self.row}".replace('-', 'm')

output_dir = "/latticeQCD/raid0/ahanlon/share/D200/data/"

phases = dict()

ensembles = ["cls21_s64_t128_D200"]

omissions = {
    "cls21_s64_t128_D200": set(list(range(0, 1000, 1)) + [1999, 2000]),
}

configs = {
    "cls21_s64_t128_D200": list(range(1001, 2000, 1)),
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
