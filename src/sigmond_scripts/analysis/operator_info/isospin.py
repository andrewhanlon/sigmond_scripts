from aenum import MultiValueEnum

class Isospin(MultiValueEnum):
  singlet = "singlet", 1, "0", "glueball", "eta", "phi", "lambda", "omega"
  doublet = "doublet", 2, "1h", "kaon", "kbar", "nucleon", "xi"
  triplet = "triplet", 3, "1", "pion", "sigma"
  quartet = "quartet", 4, "3h", "delta"
  quintet = "quintet", 5, "2"
  sextet  = "sextet", 6, "5h"
  septet  = "septet", 7, "3"

  def to_int(self):
    return self.values[1]

