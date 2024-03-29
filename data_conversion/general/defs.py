from typing import NamedTuple


class Ensemble(NamedTuple):
  name: str
  Nt: int
  replica: list
  sources: list


parity_name = {
    True: 'fwd',
    False: 'bwd',
}

ensembles = {
    "D200":
        Ensemble(
          "cls21_d200",
          128,
          ['r000'],
          [
            (35, True),
            (92, False),
          ]
        ),
    "N200":
        Ensemble(
          "cls21_n200",
          128,
          ['r000'],
          [
            (32, True),
            (52, True),
          ]
        ),
    "N203":
        Ensemble(
          "cls21_n203",
          128,
          ['r000', 'r001'],
          [
            (32, True),
            (52, True),
          ]
        ),
    "E250":
        Ensemble(
          "cls21_e250",
          192,
          ['r001'],
          [
            (0, True), (0, False),
            (1, True), (1, False),
            (2, True), (2, False),
            (3, True), (3, False),
          ]
        ),
    "C103":
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
          ]
        ),
}

source_lists = {
    "cls21_c103_r005": [10, 2, 21, 15, 12, 5, 3, 13, 3, 7, 14, 2, 23, 20, 11, 19, 0, 9, 13, 17, 21, 4, 18, 13, 20, 17, 1, 20, 17, 9, 13, 2, 11, 23, 9, 4, 13, 15, 11, 7, 19, 11, 20, 2, 6, 15, 23, 11, 5, 19, 23, 17, 23, 8, 15, 22, 17, 2, 20, 10, 21, 22, 20, 5, 12, 8, 5, 10, 22, 14, 22, 16, 20, 11, 21, 0, 1, 15, 9, 6, 16, 18, 23, 13, 20, 18, 17, 23, 12, 10, 18, 12, 17, 11, 23, 8, 16, 3, 18, 23, 22, 13, 9, 14, 8, 12, 21, 18, 3, 17, 11, 15, 17, 20, 10, 18, 1, 19, 21, 18, 19, 3, 16, 3, 10, 16, 11, 9, 6, 8, 9, 15, 11, 15, 18, 11, 12, 23, 13, 1, 16, 6, 11, 17, 18, 22, 2, 21, 16, 4, 3, 6, 1, 6, 7, 4, 9, 19, 8, 23, 17, 3, 20, 4, 3, 2, 1, 19, 8, 18, 17, 16, 16, 3, 2, 8, 6, 21, 8, 22, 21, 10, 23, 22, 1, 16, 18, 20, 1, 23, 14, 10, 5, 8, 13, 6, 17, 20, 12, 14, 10, 20, 19, 21, 1, 9, 2, 19, 16, 18, 4, 13, 12, 11, 10, 23, 16, 16, 6, 1, 10, 15, 10, 5, 9, 8, 20, 13, 13, 3, 19, 2, 11, 5, 3, 16, 22, 9, 16, 7, 4, 4, 12, 2, 3, 21, 3, 15, 14, 18, 4, 18, 0, 11, 18, 8, 0, 19, 22, 5, 17, 4, 21, 18, 6, 10, 17, 22, 6, 2, 20, 23, 9, 13, 10, 18, 15, 0, 21, 19, 15, 5, 9, 9, 23, 16, 9, 0, 22, 13, 2, 5, 16, 21, 21, 11, 12, 8, 1, 16, 14, 5, 8, 18, 12, 6, 8, 10, 23, 7, 17, 23, 13, 19, 15, 8, 22, 5, 2, 19, 16, 2, 17, 15, 17, 9, 15, 5, 11, 15, 19, 15, 4, 13, 17, 14, 20, 18, 13, 13, 6, 22, 15, 7, 2, 11, 9, 8, 12, 0, 3, 21, 5, 6, 18, 21, 23, 19, 7, 15, 1, 13, 2, 14, 4, 8, 0, 21, 14, 16, 8, 6, 9, 12, 4, 14, 6, 14, 0, 7, 3, 8, 17, 18, 21, 9, 10, 1, 11, 22, 11, 3, 5, 23, 19, 14, 23, 10, 9, 0, 17],
    "cls21_c103_r006": [11, 1, 12, 16, 15, 0, 17, 23, 5, 15, 2, 13, 3, 23, 12, 17, 6, 3, 8, 1, 12, 18, 12, 18, 14, 19, 17, 12, 15, 5, 0, 23, 0, 21, 0, 8, 18, 3, 13, 2, 4, 17, 9, 7, 14, 13, 6, 13, 8, 0, 17, 2, 11, 7, 4, 17, 9, 2, 3, 19, 14, 21, 5, 9, 7, 23, 0, 18, 9, 21, 17, 22, 6, 11, 13, 6, 6, 21, 13, 2, 13, 21, 5, 22, 7, 13, 7, 15, 18, 18, 4, 18, 14, 22, 14, 20, 10, 19, 4, 12, 2, 8, 18, 6, 10, 6, 3, 11, 12, 19, 17, 7, 19, 2, 5, 3, 7, 12, 2, 19, 17, 17, 18, 19, 5, 21, 11, 3, 11, 23, 23, 11, 9, 23, 18, 5, 12, 3, 2, 0, 13, 4, 7, 22, 16, 6, 0, 20, 21, 0, 21, 14, 15, 5, 7, 19, 1, 22, 6, 1, 8, 22, 11, 14, 5, 4, 9, 6, 14, 9, 15, 18, 9, 4, 3, 23, 5, 15, 16, 17, 13, 16, 20, 19, 5, 10, 15, 10, 18, 2, 2, 21, 15, 0, 2, 12, 16, 17, 8, 22, 17, 17, 21, 23, 4, 23, 23, 0, 7, 20, 20, 10, 3, 2, 23, 22, 6, 13, 23, 17, 21, 0, 16, 20, 14, 17, 20, 15, 12, 4, 3, 6, 3, 9, 20, 19, 4, 0, 0, 22, 6, 7, 17, 4, 16, 16, 13, 1, 20, 10, 22, 17, 18, 5, 3, 11, 13, 0, 19, 22, 7, 8, 10, 18, 5, 21, 9, 11, 20, 1, 12, 22, 0, 9, 0, 15, 12, 6, 19, 14, 16, 14, 22, 23, 2, 19, 2, 15, 12, 3, 21, 19, 12, 20, 2, 0, 2, 5, 5, 5, 1, 20, 16, 2, 6, 23, 23, 11, 5, 11, 6, 3, 23, 3, 2, 10, 14, 2, 13, 6, 9, 15, 2, 9, 22, 3, 9, 22, 9, 16, 0, 9, 21, 3, 12, 10, 2, 0, 19, 15, 1, 19, 1, 13, 7, 15, 20, 23, 18, 10, 1, 11, 5, 17, 7, 11, 18, 18, 3, 9, 23, 2, 4, 14, 3, 19, 11, 21, 23, 19, 7, 5, 9, 22, 11, 22, 7, 15, 4, 18, 4, 23, 23, 21, 3, 13, 10, 1, 19, 9, 16, 11, 3, 22, 10, 5, 1, 1, 10, 21, 11],
    "cls21_c103_r007": [3, 23, 20, 16, 5, 23, 4, 16, 11, 17, 7, 23, 15, 21, 5, 21, 13, 6, 10, 1, 20, 3, 18, 3, 8, 21, 0, 10, 0, 22, 13, 10, 6, 1, 18, 8, 14, 3, 17, 0, 8, 2, 12, 9, 4, 17, 11, 4, 10, 15, 20, 4, 10, 16, 4, 20, 16, 22, 17, 11, 20, 1, 9, 12, 9, 20, 23, 18, 9, 1, 14, 21, 12, 2, 21, 10, 16, 20, 9, 13, 11, 15, 18, 5, 10, 15, 12, 16, 7, 20, 12, 1, 15, 20, 3, 16, 7, 11, 14, 0, 15, 2, 14, 1, 15, 21, 15, 2, 1, 5, 20, 23, 12, 6, 9, 21, 6, 16, 10, 1, 18, 6, 1, 16, 9, 10, 21, 14, 2, 18, 9, 17, 12, 8, 1, 5, 13, 7, 4, 16, 19, 15, 3, 10, 14, 6, 0, 4, 14, 2, 14, 19, 5, 17, 7, 0, 9, 2, 20, 0, 13, 21, 10, 1, 23, 7, 10, 9, 2, 5, 8, 17, 14, 3, 16, 11, 18, 6, 16, 10, 14, 18, 9, 22, 7, 15, 6, 16, 19, 16, 4, 9, 6, 9, 18, 14, 11, 2, 14, 10, 1, 22, 13, 21, 12, 16, 1, 3, 8, 1, 10, 23, 1, 15, 11, 7, 13, 22, 2, 7, 15, 12, 17, 3, 21, 16, 4, 19, 13, 1, 13, 0, 2, 19, 8, 3, 13, 1, 13, 0, 7, 3, 10, 4, 2, 11, 2, 14, 0, 1, 3, 4, 15, 2, 11, 14, 2, 20, 2, 23, 19, 11, 7, 15, 4, 18, 22, 10, 15, 18, 2, 7, 10, 17, 5, 13, 20, 23, 13, 2, 16, 8, 20, 23, 14, 19, 2, 21, 4, 19, 8, 6, 18, 0, 6, 12, 16, 15, 5, 23, 6, 9, 15, 20, 16, 8, 7, 10, 13, 22, 9, 15, 22, 15, 8, 20, 2, 10, 6, 14, 17, 20, 14, 20, 2, 12, 6, 15, 7, 14, 1, 20, 7, 19, 10, 15, 5, 11, 19, 12, 15, 18, 13, 4, 7, 12, 16, 2, 8, 5, 3, 14, 19, 12, 20, 16, 7, 4, 19, 18, 23, 6, 17, 9, 10, 9, 11, 0, 2, 10, 12, 6, 20, 7, 18, 14, 23, 15, 2, 16, 8, 21, 13, 9, 11, 14, 9, 4, 20, 14, 5, 7, 5, 13, 11, 7, 21, 19, 10, 6, 0, 11, 18, 18],
    "cls21_c103_r008": [18, 9, 2, 13, 6, 23, 12, 2, 15, 9, 3, 22, 17, 22, 12, 17, 5, 12, 4, 18, 8, 4, 12, 6, 11, 15, 23, 6, 18, 12, 3, 9, 23, 9, 20, 8, 16, 1, 23, 5, 13, 15, 11, 5, 0, 14, 23, 7, 2, 13, 23, 10, 16, 23, 11, 6, 18, 9, 20, 3, 17, 21, 14, 5, 18, 6, 0, 3, 14, 7, 12, 9, 23, 2, 12, 18, 6, 17, 3, 17, 4, 10, 3, 14, 21, 17, 13, 5, 15, 1, 13, 21, 2, 7, 13, 21, 4, 15, 7, 10, 18, 5, 23, 19, 11, 19, 2, 11, 13, 1, 0, 8, 20, 23, 13, 21, 2, 9, 19, 4, 19, 14, 7, 15, 8, 18, 21, 17, 12, 7, 2, 15, 23, 18, 11, 6, 12, 17, 6, 1, 18, 3, 22, 11, 16, 21, 9, 12, 23, 7, 17, 23, 1, 11, 3, 19, 13, 4, 14, 11, 21, 10, 19, 9, 13, 6, 19, 16, 11, 6, 3, 19, 10, 16, 3, 9, 6, 23, 13, 8, 22, 1, 12, 2, 11, 23, 6, 8, 20, 6, 14, 10, 8, 14, 5, 20, 11, 6, 12, 0, 10, 8, 22, 16, 13, 8, 11, 1, 14, 4, 10, 20, 1, 15, 9, 18, 1, 4, 8, 0, 23, 11, 17, 20, 7, 23, 10, 12, 20, 9, 18, 14, 4, 1, 12, 4, 12, 2, 22, 11, 0, 8, 17, 22, 15, 21, 3, 22, 12, 19, 11, 5, 20, 11, 1, 6, 20, 11, 10, 0, 6, 19, 23, 12, 20, 9, 6, 0, 20, 2, 21, 2, 5, 23, 10, 8, 23, 9, 16, 14, 2, 7, 19, 16, 0, 9, 22, 5, 21, 13, 23, 5, 14, 20, 1, 20, 12, 8, 16, 21, 12, 5, 23, 18, 2, 23, 14, 5, 11, 18, 11, 7, 21, 13, 18, 1, 17, 2, 18, 10, 23, 13, 1, 23, 4, 0, 10, 17, 8, 6, 1, 19, 4, 14, 10, 3, 13, 21, 7, 4, 12, 1, 23, 12, 16, 22, 2, 23, 11, 13, 23, 12, 1, 5, 15, 9, 5, 13, 18, 22, 12, 1, 3, 11, 7, 22, 9, 0, 14, 10, 11, 15, 12, 23, 9, 21, 14, 23, 12, 1, 23, 16, 10, 22, 15, 22, 6, 14, 23, 8, 0, 5, 22, 12, 5, 16, 6, 10, 7, 23, 14, 20, 21, 15],
}

source_mods = {
    "cls21_c103_r005": 12,
    "cls21_c103_r006": 12,
    "cls21_c103_r007": 12,
    "cls21_c103_r008": 12,
}

config_indices = {
    "cls21_c103_r005": list(range(401)),
    "cls21_c103_r006": list(range(401)),
    "cls21_c103_r007": list(range(404)),
    "cls21_c103_r008": sorted(list([0,1,2,3,6,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,25,27,30,31,34,35,37,40,41,44,45,47,50,51,54,55,57,59,62,63,66,68,69,72,73,76,78,79,82,83,86,88,89,92,93,96,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,244,245,246,248,249,252,253,256,258,260,261,264,265,267,270,271,274,275,277,280,281,284,285,287,290,291,294,295,297,300,301,304,305,307,310,311,314,315,317,319,322,323,326,328,329,332,333,336,338,339,342,343,346,348,349,352,353,356,358,360,361,364,365,367,370,371,374,375,377,380,381,384,385,387,390,391,394,395,397,399,402,403])),
    "cls21_d200_r000": list(range(1, 2001, 1)),
    "cls21_n200_r000": list(range(1, 856, 2)),
    "cls21_n203_r000": list(range(1, 757, 2)),
    "cls21_n203_r001": list(range(2, 787, 2)),
    "cls21_e250_r001": list(range(1, 502, 4)) + list(range(503, 902, 2)) + list(range(905, 1010, 4)),
}

omissions = {
    "cls21_c103_r005": set(list(range(401))) - set(config_indices['cls21_c103_r005']),
    "cls21_c103_r006": set(list(range(401))) - set(config_indices['cls21_c103_r006']),
    "cls21_c103_r007": set(list(range(404))) - set(config_indices['cls21_c103_r007']),
    "cls21_c103_r008": set(list(range(404))) - set(config_indices['cls21_c103_r008']),
    "cls21_d200_r000": set([2000]),
    "cls21_n200_r000": set(range(1, 856, 2)),
    "cls21_n203_r000": set(range(1, 757, 2)),
    "cls21_n203_r001": set(range(0, 788, 2)),
    "cls21_e250_r001": set(list(range(0, 1009))) -  set(list(range(0, 501, 4)) + list(range(502, 901, 2)) + list(range(904, 1009, 4))),
}


isospin_map = {
    'nucleon': 'doublet',
    'pion': 'triplet',
    'kaon': 'doublet',
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

    return f"{mom_str}_{self.irrep}_{self.irrep_row}_iso{self.isospin}_{strangeness_str}"

  def data_channel_str(self):
    mom_str = f"P{self.momentum[0]}{self.momentum[1]}{self.momentum[2]}".replace('-', 'm')

    return f"{mom_str}_{self.irrep}_{self.irrep_row}"

  def averaged_channel_str(self):
    psq = self.momentum[0]**2 + self.momentum[1]**2 + self.momentum[2]**2
    strangeness_str = f"S{self.strangeness}".replace('-', 'm')
    return f"PSQ{psq}_{self.irrep}_iso{self.isospin}_{strangeness_str}"

  def iso_strange_str(self):
    strangeness_str = f"S{self.strangeness}".replace('-', 'm')
    return f"iso{self.isospin}_{strangeness_str}"




class AveragedChannel(NamedTuple):
  psq: int
  irrep: str
  isospin: str
  strangeness: int

  def iso_strange_str(self):
    strangeness_str = f"S{self.strangeness}".replace('-', 'm')
    return f"iso{self.isospin}_{strangeness_str}"

  def data_channel_str(self):
    return f"PSQ{self.psq}_{self.irrep}"

  def __str__(self):
    return f"(psq={self.psq}, irrep={self.irrep}, isospin={self.isospin}, S={self.strangeness})"

  def __repr__(self):
    return f"PSQ{self.psq}-{self.irrep}-{self.isospin}-{self.strangeness}"

  def __hash__(self):
    return hash(repr(self))

  def __eq__(self, other):
    return repr(self) == repr(other)
