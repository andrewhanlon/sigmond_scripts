from typing import NamedTuple

import sourceTimeList

class Ensemble(NamedTuple):
  name: str
  Nt: int
  replica: list
  sources: list
  parities: list

ensembles = [
    Ensemble(
      "cls21_c103",
      96,
      ['r005', 'r006', 'r007', 'r008'],
      [
        (0, True), (0, False),
        (12, True), (12, False),
        (24, True), (24, False),
        (36, True), (36, False),
        (48, True), (48, False),
        (60, True), (60, False),
        (72, True), (72, False),
        (84, True), (84, False),
      ],
      ['fwd', 'bwd']
    ),
]

source_lists = {
    "cls21_c103_r005": sourceTimeList.r005,
    "cls21_c103_r006": sourceTimeList.r006,
    "cls21_c103_r007": sourceTimeList.r007,
    "cls21_c103_r008": sourceTimeList.r008,
}

config_indices = {
    "cls21_c103_r005": list(range(401)),
    "cls21_c103_r006": list(range(401)),
    "cls21_c103_r007": list(range(404)),
    "cls21_c103_r008": sorted(list([0,1,2,3,6,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,25,27,30,31,34,35,37,40,41,44,45,47,50,51,54,55,57,59,62,63,66,68,69,72,73,76,78,79,82,83,86,88,89,92,93,96,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,244,245,246,248,249,252,253,256,258,260,261,264,265,267,270,271,274,275,277,280,281,284,285,287,290,291,294,295,297,300,301,304,305,307,310,311,314,315,317,319,322,323,326,328,329,332,333,336,338,339,342,343,346,348,349,352,353,356,358,360,361,364,365,367,370,371,374,375,377,380,381,384,385,387,390,391,394,395,397,399,402,403]))
}

omissions = {
    "cls21_c103_r005": set(list(range(401))) - set(config_indices['cls21_c103_r005']),
    "cls21_c103_r006": set(list(range(401))) - set(config_indices['cls21_c103_r006']),
    "cls21_c103_r007": set(list(range(404))) - set(config_indices['cls21_c103_r007']),
    "cls21_c103_r008": set(list(range(404))) - set(config_indices['cls21_c103_r008']),
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

  def averaged_channel_str(self):
    psq = self.momentum[0]**2 + self.momentum[1]**2 + self.momentum[2]**2
    strangeness_str = f"S{self.strangeness}".replace('-', 'm')
    return f"PSQ{psq}_{self.irrep}_{self.isospin}_{strangeness_str}"

  def iso_strange_str(self):
    strangeness_str = f"S{self.strangeness}".replace('-', 'm')
    return f"{self.isospin}_{strangeness_str}"




class AveragedChannel(NamedTuple):
  psq: int
  irrep: str
  isospin: str
  strangeness: int

  def __str__(self):
    return f"(psq={self.psq}, irrep={self.irrep}, isospin={self.isospin}, S={self.strangeness})"

  def __repr__(self):
    return f"pSq{self.psq}-{self.irrep}-{self.isospin}-{self.strangeness}"

class EquivalentFrame(NamedTuple):
  mom: tuple
  row: int

  def __str__(self):
    return f"(P={self.mom}, row={self.row})"

  def __repr__(self):
    return f"P{''.join(list(map(str, self.mom)))}_{self.row}".replace('-', 'm')

