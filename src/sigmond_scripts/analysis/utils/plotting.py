import numpy as np

import pylatex
import matplotlib.pyplot as plt

import utils.util as util

PLOT_MARKS = [
    '*',
    'triangle*',
    'halfsquare*',
    'halfcircle*',
    'halfsquare right*',
    'halfsquare left*',
    'square*',
]

COLORS = [
    'red',
    'teal',
    'cyan',
]

INTERPOLATION = 100

NUM_YTICKS = 5
NUM_XTICKS = 8

NUM_MIN_TICKS = 4

TIK_SIZE = 0.1
MIN_TIK_SIZE = 0.05

FRAME_HEIGHT = 4.45
FRAME_WIDTH = 7.2
PERCENT_ROOM = 0.05

SAMPS_TO_PLOT = 100
SAMPS_PERCENT = 0.16


def calculate_position(pos, x_range, y_range):
  x_pos = FRAME_WIDTH*(pos[0] - x_range[0])/(x_range[1] - x_range[0])
  y_pos = FRAME_HEIGHT*(pos[1] - y_range[0])/(y_range[1] - y_range[0])

  return (x_pos, y_pos)

def get_ticks(x_range, y_range):
  ax = plt.gca()
  ax.plot(x_range, y_range)

  x_ticks = list()
  for x_tick in ax.get_xticks():
    if x_range[0] <= x_tick <= x_range[1]:
      x_ticks.append(round(x_tick, 12))

  y_ticks = list()
  for y_tick in ax.get_yticks():
    if y_range[0] <= y_tick <= y_range[1]:
      y_ticks.append(round(y_tick, 12))

  return (x_ticks, y_ticks)

def get_minor_ticks(x_range, y_range):
  x_ticks, y_ticks = get_ticks(x_range, y_range)

  x_ticks.insert(0, x_ticks[0] - (x_ticks[1] - x_ticks[0]))
  x_ticks.append(x_ticks[-1] + (x_ticks[-1] - x_ticks[-2]))

  y_ticks.insert(0, y_ticks[0] - (y_ticks[1] - y_ticks[0]))
  y_ticks.append(y_ticks[-1] + (y_ticks[-1] - y_ticks[-2]))

  x_interp = x_ticks[1] - x_ticks[0]
  x_min_ticks = list()
  for x_tick in x_ticks:
    for x_min_tick in np.linspace(x_tick, x_tick+x_interp, NUM_MIN_TICKS+2):
      if x_range[0] <= x_min_tick <= x_range[1] and x_min_tick not in x_ticks:
        x_min_ticks.append(round(x_min_tick, 12))

  y_interp = y_ticks[1] - y_ticks[0]
  y_min_ticks = list()
  for y_tick in y_ticks:
    for y_min_tick in np.linspace(y_tick, y_tick+y_interp, NUM_MIN_TICKS+2):
      if y_range[0] <= y_min_tick <= y_range[1] and y_min_tick not in y_ticks:
        y_min_ticks.append(round(y_min_tick, 12))
  
  return (x_min_ticks, y_min_ticks)


def spectrum(thresholds, energies, non_interacting_energies, latex_map, rotate_irreps,
             plot_width_factor, ref_name, filename, latex_compiler, non_interacting_energy_names=None):
  global FRAME_WIDTH
  FRAME_WIDTH *= plot_width_factor
  global PERCENT_ROOM
  PERCENT_ROOM /= plot_width_factor
  # START remove
  '''
  global FRAME_HEIGHT
  FRAME_WIDTH = 12.21
  FRAME_HEIGHT = 6.789
  PERCENT_ROOM *= plot_width_factor
  '''
  # END remove

  # create doc and tikz
  doc = util.create_tikz_doc()
  tikz_filename = f"{filename}.tikz"
  doc.append(pylatex.NoEscape(rf"\input{{spectrum.tikz}}"))

  error_bar_options = [
      pylatex.NoEscape(r"lebar/.style={black,postaction={decorate,decoration={markings, mark=at position 0.0 with {\draw (0pt,-2.4pt) -- ++(0,4.8pt);}}}}"),
      pylatex.NoEscape(r"uebar/.style={black,postaction={decorate,decoration={markings, mark=at position 1.0 with {\draw (0pt,-2.4pt) -- ++(0,4.8pt);}}}}"),
  ]

  tikz_pic = pylatex.TikZ(options=error_bar_options)
  # draw boundary
  border = rf"\draw[black] ({FRAME_WIDTH},0.0) -- (0.0,0.0) -- (0.0,{FRAME_HEIGHT});"
  tikz_pic.append(pylatex.NoEscape(border))

  # find energy_range
  energie_vals = list() #change to where it doen't fuck up if there are no non interacting levels decided
  for energy_lists, non_energy_lists in zip(energies.values(), non_interacting_energies.values()):
    for energy in energy_lists:
      energy_mean = energy.getFullEstimate()
      energy_err = energy.getSymmetricError()
      energie_vals.append(energy_mean+energy_err)
      energie_vals.append(energy_mean-energy_err)

    for energy in non_energy_lists:
      energy_mean = energy.getFullEstimate()
      energy_err = energy.getSymmetricError()
      energie_vals.append(energy_mean+energy_err)
      energie_vals.append(energy_mean-energy_err)

  energy_min = min(energie_vals)
  energy_max = max(energie_vals)
  energy_room = (energy_max - energy_min)*PERCENT_ROOM
  energy_range = (energy_min - energy_room, energy_max + energy_room)
  x_range = (0., FRAME_WIDTH)

  # Draw thresholds
  for latex, threshold in thresholds.items():
    threshold_val = threshold.getFullEstimate()
    left_pos = calculate_position((0.0,threshold_val), x_range, energy_range)
    right_pos = calculate_position((FRAME_WIDTH,threshold_val), x_range, energy_range)
    thres_tikz = rf"\draw[densely dotted,color=black,line width=0.2mm] ({left_pos[0]},{left_pos[1]}) -- ({right_pos[0]}, {right_pos[1]});"
    tikz_pic.append(pylatex.NoEscape(thres_tikz))
    thres_tikz = rf"\draw [anchor=west] ({right_pos[0]},{right_pos[1]}) " \
                 rf"node{{\footnotesize {latex}}};"
    tikz_pic.append(pylatex.NoEscape(thres_tikz))

  _, y_ticks = get_ticks(x_range, energy_range)

  # draw energy ticks
  for y_tick in y_ticks:
    x_val = x_range[0]
    pos = calculate_position((x_val, y_tick), x_range, energy_range)

    tikz_tik = rf"\draw[black] ({pos[0]},{pos[1]}) -- +({TIK_SIZE},0.0);"
    tikz_tik_lab = rf'\node[anchor=east, inner sep=0.2em] at ({pos[0]},{pos[1]})' \
                   rf'{{\footnotesize ${{{y_tick}}}$}};'

    tikz_pic.append(pylatex.NoEscape(tikz_tik))
    tikz_pic.append(pylatex.NoEscape(tikz_tik_lab))

  _, y_min_ticks = get_minor_ticks(x_range, energy_range)

  for y_min_tick in y_min_ticks:
    x_val = x_range[0]
    pos = calculate_position((x_val, y_min_tick), x_range, energy_range)

    tikz_tik = rf"\draw[black] ({pos[0]},{pos[1]}) -- +({MIN_TIK_SIZE},0.0);"

    tikz_pic.append(pylatex.NoEscape(tikz_tik))

  # Draw energy label
  energy_label_pos = calculate_position((x_range[0]-.4, energy_range[1]-.1), x_range, energy_range)
  energy_label_tikz = rf"\draw [anchor=north east] ({energy_label_pos[0]}, {energy_label_pos[1]})"
  if ref_name in latex_map:
    ref_mass_latex = latex_map[ref_name]
    energy_label_tikz += rf" node{{\footnotesize $\frac{{E_\mathrm{{cm}}}}{{m_{ref_mass_latex}}}$}};"
  else:
    energy_label_tikz += rf" node{{\footnotesize $E_\mathrm{{cm}}$}};"

  tikz_pic.append(pylatex.NoEscape(energy_label_tikz))


  # draw irrep labels
  num_channels = len(energies)
  label_width = FRAME_WIDTH / (num_channels+1)
  irrep_label_poss = np.linspace(label_width/2.+FRAME_WIDTH*PERCENT_ROOM,
                                 FRAME_WIDTH-label_width/2-FRAME_WIDTH*PERCENT_ROOM,
                                 num_channels)
  '''
  irrep_label_poss = np.linspace(label_width/2., FRAME_WIDTH-label_width/2, num_channels)
  '''

  print("irrep positions")
  for irrep_label_pos, irrep_label in zip(irrep_label_poss, energies.keys()):
    print(f"{irrep_label}: {irrep_label_pos}")
    if rotate_irreps:
      tikz_lab = rf"\node[anchor=north east, rotate=45, inner sep=0.2em] at ({irrep_label_pos},0.)" \
                 rf"{{\footnotesize {{{irrep_label}}}}};"
    else:
      tikz_lab = rf"\node[anchor=north, inner sep=0.2em] at ({irrep_label_pos},0.)" \
                 rf"{{\footnotesize {{{irrep_label}}}}};"

    tikz_pic.append(pylatex.NoEscape(tikz_lab))


  # Draw non-interacting levels
  box_width = label_width*.8
  for counter, (x_pos, non_energy_list) in enumerate(zip(irrep_label_poss, non_interacting_energies.values())):
    for non_energy_est in non_energy_list:
      non_energy_val = non_energy_est.getFullEstimate()
      non_energy_err = non_energy_est.getSymmetricError()
      non_energy_err_p = non_energy_val + non_energy_err
      non_energy_err_m = non_energy_val - non_energy_err

      box_left = calculate_position((x_pos-box_width/2, non_energy_err_m), x_range, energy_range)
      box_right = calculate_position((x_pos+box_width/2, non_energy_err_p), x_range, energy_range)
      line_left = calculate_position((x_pos-box_width/2, non_energy_val), x_range, energy_range)
      line_right = calculate_position((x_pos+box_width/2, non_energy_val), x_range, energy_range)

      non_tikz = rf"\fill[lightgray] ({box_left[0]},{box_left[1]})" \
                 rf" rectangle ({box_right[0]},{box_right[1]});"
      tikz_pic.append(pylatex.NoEscape(non_tikz))

  # I want the non-interacting dashes on top. Should use z-level, but it's confusing :)
  for counter, (x_pos, non_energy_list) in enumerate(zip(irrep_label_poss, non_interacting_energies.values())):
    for non_energy_est in non_energy_list:
      non_energy_val = non_energy_est.getFullEstimate()
      non_energy_err = non_energy_est.getSymmetricError()
      non_energy_err_p = non_energy_val + non_energy_err
      non_energy_err_m = non_energy_val - non_energy_err

      box_left = calculate_position((x_pos-box_width/2, non_energy_err_m), x_range, energy_range)
      box_right = calculate_position((x_pos+box_width/2, non_energy_err_p), x_range, energy_range)
      line_left = calculate_position((x_pos-box_width/2, non_energy_val), x_range, energy_range)
      line_right = calculate_position((x_pos+box_width/2, non_energy_val), x_range, energy_range)

      non_tikz = r"\draw[dash pattern=on 1pt off 1pt, black] " \
                 rf"({line_left[0]},{line_left[1]}) -- " \
                 rf"({line_right[0]},{line_right[1]});"
      tikz_pic.append(pylatex.NoEscape(non_tikz))

  #label non interacting enerrgies.
  if non_interacting_energy_names is not None:
    for counter, (x_pos, non_energy_list, irrep) in enumerate(zip(irrep_label_poss, non_interacting_energies.values(),non_interacting_energy_names)):
      for i, (non_energy_est) in enumerate(non_energy_list):
        non_energy_val = non_energy_est.getFullEstimate()
        line_right = calculate_position((x_pos+box_width/2, non_energy_val), x_range, energy_range)

        level_name = str(non_interacting_energy_names[irrep][i][0][0])+"_"+str(non_interacting_energy_names[irrep][i][1][0])+" " \
                   + str(non_interacting_energy_names[irrep][i][0][1])+"_"+str(non_interacting_energy_names[irrep][i][1][1])
        non_tikz_labels = rf"\node[scale=0.5] at ({line_right[0]}+0.3,{line_right[1]}) {{${level_name}$}};"
        tikz_pic.append(pylatex.NoEscape(non_tikz_labels))

  # Draw energy levels
  box_width = label_width*.9
  for counter, (x_pos, energy_list) in enumerate(zip(irrep_label_poss, energies.values())):
    for energy_est in energy_list:
      energy_val = energy_est.getFullEstimate()
      energy_err = energy_est.getSymmetricError()
      energy_err_p = energy_val + energy_err
      energy_err_m = energy_val - energy_err

      energy_val = calculate_position((x_pos, energy_val), x_range, energy_range)
      energy_err_m = calculate_position((x_pos, energy_err_m), x_range, energy_range)
      energy_err_p = calculate_position((x_pos, energy_err_p), x_range, energy_range)

      en_tikz = rf"\draw[lebar,uebar] ({energy_err_m[0]},{energy_err_m[1]}) -- ({energy_err_p[0]},{energy_err_p[1]});"
      tikz_pic.append(pylatex.NoEscape(en_tikz))
      en_tikz = r"\draw[mark options={color=black, fill=white}, mark size=1.5pt,mark=*] " \
                rf"plot coordinates {{({energy_val[0]},{energy_val[1]})}};"
      tikz_pic.append(pylatex.NoEscape(en_tikz))

  util.write_tikz(tikz_pic, tikz_filename)
  util.compile_pdf(doc, filename, latex_compiler)

  print(f"x_range: {x_range}")
  print(f"energy_range: {energy_range}")
