from typing import NamedTuple

base_data_dir = "/disk2/research/data/3K_3pi"
output_dir = "/disk1/research/projects/3K_3pi/data"

class Ensemble(NamedTuple):
  name: str
  replica: list
  sources: list
  flavor_channels: list

'''
ensembles = [
    Ensemble("cls21_n203", ['r000', 'r001'], [32, 52], ['kaons', 'pions', 'KKp_Kpp']),
    Ensemble("cls21_n200", ['r000', 'r001'], [32, 52], ['kaons', 'pions', 'KKp_Kpp']),
    Ensemble("cls21_d200", ['r000'], [35, 92], ['kaons', 'pions', 'KKp_Kpp']),
]
'''
ensembles = [
    Ensemble(
        "cls21_n203",
        ['r000', 'r001'],
        [(32, 'fwd'), (52, 'fwd')],
        ['kaons', 'pions']
    ),
    Ensemble(
        "cls21_n200",
        ['r000', 'r001'],
        [(32, 'fwd'), (52, 'fwd')],
        ['kaons', 'pions']
    ),
    Ensemble(
        "cls21_d200",
        ['r000'],
        [(35, 'fwd'), (92, 'bwd')],
        ['kaons', 'pions']
    ),
    Ensemble(
        "cls21_e250",
        ['r001'],
        [(0, 'fwd'), (0, 'bwd'), (1, 'fwd'), (1, 'bwd'), (2, 'fwd'), (2, 'bwd'), (3, 'fwd'), (3, 'bwd')],
        ['pions']
    ),
]

omissions = {
    'cls21_n203_r000': set(range(1, 756, 2)),
    'cls21_n203_r001': set(range(0, 787, 2)),
    'cls21_n203': set(range(1, 756, 2)).union(set(range(756, 1544, 2))),
    'cls21_n200_r000': set(),
    'cls21_n200_r001': set(),
    'cls21_n200': set(),
    'cls21_d200_r000': set([2000]),
    'cls21_d200': set([2000]),
    "cls21_e250_r001": {1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 13, 14, 15, 17, 18, 19, 21, 22, 23, 25, 26, 27, 29, 30, 31, 33, 34, 35, 37, 38, 39, 41, 42, 43, 45, 46, 47, 49, 50, 51, 53, 54, 55, 57, 58, 59, 61, 62, 63, 65, 66, 67, 69, 70, 71, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 93, 94, 95, 97, 98, 99, 101, 102, 103, 105, 106, 107, 109, 110, 111, 113, 114, 115, 117, 118, 119, 121, 122, 123, 125, 126, 127, 129, 130, 131, 133, 134, 135, 137, 138, 139, 141, 142, 143, 145, 146, 147, 149, 150, 151, 153, 154, 155, 157, 158, 159, 161, 162, 163, 165, 166, 167, 169, 170, 171, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 193, 194, 195, 197, 198, 199, 201, 202, 203, 204, 205, 206, 207, 209, 210, 211, 213, 214, 215, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 229, 230, 231, 233, 234, 235, 236, 237, 238, 239, 241, 242, 243, 245, 246, 247, 249, 250, 251, 253, 254, 255, 257, 258, 259, 261, 262, 263, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 277, 278, 279, 281, 282, 283, 285, 286, 287, 288, 289, 290, 291, 293, 294, 295, 297, 298, 299, 301, 302, 303, 305, 306, 307, 309, 310, 311, 313, 314, 315, 317, 318, 319, 321, 322, 323, 325, 326, 327, 329, 330, 331, 333, 334, 335, 337, 338, 339, 341, 342, 343, 345, 346, 347, 349, 350, 351, 353, 354, 355, 357, 358, 359, 361, 362, 363, 365, 366, 367, 369, 370, 371, 372, 373, 374, 375, 377, 378, 379, 381, 382, 383, 385, 386, 387, 389, 390, 391, 393, 394, 395, 397, 398, 399, 401, 402, 403, 405, 406, 407, 409, 410, 411, 413, 414, 415, 417, 418, 419, 421, 422, 423, 425, 426, 427, 429, 430, 431, 433, 434, 435, 437, 438, 439, 441, 442, 443, 444, 445, 446, 447, 449, 450, 451, 453, 454, 455, 457, 458, 459, 461, 462, 463, 465, 466, 467, 469, 470, 471, 473, 474, 475, 476, 477, 478, 479, 481, 482, 483, 485, 486, 487, 488, 489, 490, 491, 493, 494, 495, 497, 498, 499, 501, 502, 503, 505, 506, 507, 509, 510, 511, 513, 514, 515, 517, 518, 519, 521, 522, 523, 525, 526, 527, 529, 530, 531, 533, 534, 535, 537, 538, 539, 541, 542, 543, 545, 546, 547, 549, 550, 551, 553, 554, 555, 557, 558, 559, 561, 562, 563, 565, 566, 567, 569, 570, 571, 573, 574, 575, 577, 578, 579, 581, 582, 583, 585, 586, 587, 589, 590, 591, 593, 594, 595, 597, 598, 599, 601, 602, 603, 605, 606, 607, 609, 610, 611, 613, 614, 615, 617, 618, 619, 621, 622, 623, 625, 626, 627, 629, 630, 631, 633, 634, 635, 637, 638, 639, 641, 642, 643, 645, 646, 647, 649, 650, 651, 653, 654, 655, 657, 658, 659, 661, 662, 663, 665, 666, 667, 669, 670, 671, 673, 674, 675, 677, 678, 679, 681, 682, 683, 685, 686, 687, 689, 690, 691, 693, 694, 695, 697, 698, 699, 701, 702, 703, 705, 706, 707, 709, 710, 711, 713, 714, 715, 717, 718, 719, 721, 722, 723, 725, 726, 727, 729, 730, 731, 733, 734, 735, 737, 738, 739, 741, 742, 743, 745, 746, 747, 749, 750, 751, 753, 754, 755, 757, 758, 759, 761, 762, 763, 765, 766, 767, 769, 770, 771, 773, 774, 775, 777, 778, 779, 781, 782, 783, 785, 786, 787, 789, 790, 791, 793, 794, 795, 797, 798, 799, 801, 802, 803, 805, 806, 807, 809, 810, 811, 813, 814, 815, 817, 818, 819, 821, 822, 823, 825, 826, 827, 829, 830, 831, 833, 834, 835, 837, 838, 839, 841, 842, 843, 845, 846, 847, 849, 850, 851, 853, 854, 855, 857, 858, 859, 861, 862, 863, 865, 866, 867, 869, 870, 871, 873, 874, 875, 877, 878, 879, 881, 882, 883, 885, 886, 887, 889, 890, 891, 893, 894, 895, 897, 898, 899, 901, 902, 903, 905, 906, 907, 909, 910, 911, 913, 914, 915, 917, 918, 919, 920, 921, 922, 923, 924, 925, 926, 927, 928, 929, 930, 931, 932, 933, 934, 935, 936, 937, 938, 939, 940, 941, 942, 943, 944, 945, 946, 947, 948, 949},
    "cls21_e250": {1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 13, 14, 15, 17, 18, 19, 21, 22, 23, 25, 26, 27, 29, 30, 31, 33, 34, 35, 37, 38, 39, 41, 42, 43, 45, 46, 47, 49, 50, 51, 53, 54, 55, 57, 58, 59, 61, 62, 63, 65, 66, 67, 69, 70, 71, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 93, 94, 95, 97, 98, 99, 101, 102, 103, 105, 106, 107, 109, 110, 111, 113, 114, 115, 117, 118, 119, 121, 122, 123, 125, 126, 127, 129, 130, 131, 133, 134, 135, 137, 138, 139, 141, 142, 143, 145, 146, 147, 149, 150, 151, 153, 154, 155, 157, 158, 159, 161, 162, 163, 165, 166, 167, 169, 170, 171, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 193, 194, 195, 197, 198, 199, 201, 202, 203, 204, 205, 206, 207, 209, 210, 211, 213, 214, 215, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 229, 230, 231, 233, 234, 235, 236, 237, 238, 239, 241, 242, 243, 245, 246, 247, 249, 250, 251, 253, 254, 255, 257, 258, 259, 261, 262, 263, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 277, 278, 279, 281, 282, 283, 285, 286, 287, 288, 289, 290, 291, 293, 294, 295, 297, 298, 299, 301, 302, 303, 305, 306, 307, 309, 310, 311, 313, 314, 315, 317, 318, 319, 321, 322, 323, 325, 326, 327, 329, 330, 331, 333, 334, 335, 337, 338, 339, 341, 342, 343, 345, 346, 347, 349, 350, 351, 353, 354, 355, 357, 358, 359, 361, 362, 363, 365, 366, 367, 369, 370, 371, 372, 373, 374, 375, 377, 378, 379, 381, 382, 383, 385, 386, 387, 389, 390, 391, 393, 394, 395, 397, 398, 399, 401, 402, 403, 405, 406, 407, 409, 410, 411, 413, 414, 415, 417, 418, 419, 421, 422, 423, 425, 426, 427, 429, 430, 431, 433, 434, 435, 437, 438, 439, 441, 442, 443, 444, 445, 446, 447, 449, 450, 451, 453, 454, 455, 457, 458, 459, 461, 462, 463, 465, 466, 467, 469, 470, 471, 473, 474, 475, 476, 477, 478, 479, 481, 482, 483, 485, 486, 487, 488, 489, 490, 491, 493, 494, 495, 497, 498, 499, 501, 502, 503, 505, 506, 507, 509, 510, 511, 513, 514, 515, 517, 518, 519, 521, 522, 523, 525, 526, 527, 529, 530, 531, 533, 534, 535, 537, 538, 539, 541, 542, 543, 545, 546, 547, 549, 550, 551, 553, 554, 555, 557, 558, 559, 561, 562, 563, 565, 566, 567, 569, 570, 571, 573, 574, 575, 577, 578, 579, 581, 582, 583, 585, 586, 587, 589, 590, 591, 593, 594, 595, 597, 598, 599, 601, 602, 603, 605, 606, 607, 609, 610, 611, 613, 614, 615, 617, 618, 619, 621, 622, 623, 625, 626, 627, 629, 630, 631, 633, 634, 635, 637, 638, 639, 641, 642, 643, 645, 646, 647, 649, 650, 651, 653, 654, 655, 657, 658, 659, 661, 662, 663, 665, 666, 667, 669, 670, 671, 673, 674, 675, 677, 678, 679, 681, 682, 683, 685, 686, 687, 689, 690, 691, 693, 694, 695, 697, 698, 699, 701, 702, 703, 705, 706, 707, 709, 710, 711, 713, 714, 715, 717, 718, 719, 721, 722, 723, 725, 726, 727, 729, 730, 731, 733, 734, 735, 737, 738, 739, 741, 742, 743, 745, 746, 747, 749, 750, 751, 753, 754, 755, 757, 758, 759, 761, 762, 763, 765, 766, 767, 769, 770, 771, 773, 774, 775, 777, 778, 779, 781, 782, 783, 785, 786, 787, 789, 790, 791, 793, 794, 795, 797, 798, 799, 801, 802, 803, 805, 806, 807, 809, 810, 811, 813, 814, 815, 817, 818, 819, 821, 822, 823, 825, 826, 827, 829, 830, 831, 833, 834, 835, 837, 838, 839, 841, 842, 843, 845, 846, 847, 849, 850, 851, 853, 854, 855, 857, 858, 859, 861, 862, 863, 865, 866, 867, 869, 870, 871, 873, 874, 875, 877, 878, 879, 881, 882, 883, 885, 886, 887, 889, 890, 891, 893, 894, 895, 897, 898, 899, 901, 902, 903, 905, 906, 907, 909, 910, 911, 913, 914, 915, 917, 918, 919, 920, 921, 922, 923, 924, 925, 926, 927, 928, 929, 930, 931, 932, 933, 934, 935, 936, 937, 938, 939, 940, 941, 942, 943, 944, 945, 946, 947, 948, 949},
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
