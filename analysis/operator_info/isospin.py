from aenum import MultiValueEnum

class Isospin(MultiValueEnum):
  singlet = "singlet", 1, "glueball", "eta", "phi", "lambda", "omega"
  doublet = "doublet", 2, "kaon", "kbar", "nucleon", "xi"
  triplet = "triplet", 3, "pion", "sigma"
  quartet = "quartet", 4, "delta"
  quintet = "quintet", 5
  sextet  = "sextet", 6
  septet  = "septet", 7

  def to_int(self):
    return self.values[1]

  def to_str(self):
    val = self.to_int()-1
    if (val % 2) == 0:
      return str(val/2)
    else:
      return f"{str(val)}h"

