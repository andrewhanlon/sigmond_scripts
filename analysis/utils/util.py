import os
import string
import logging
import subprocess
import datetime
import math
import numpy as np
import regex
import pylatex
import xml.etree.ElementTree as ET
import xml.dom.minidom as minidom
import uncertainties
from sortedcontainers import SortedSet
from aenum import MultiValueEnum

import operator_info.operator
import sigmond


class PlotExtension(MultiValueEnum):
  grace = "agr", "grace"
  pdf = "pdf"

ERR_PREC=2


class Singleton(type):
  _instances = {}
  def __call__(cls, *args, **kwargs):
    if cls not in cls._instances:
      cls._instances[cls] = super(Singleton, cls).__call__(*args, **kwargs)
    return cls._instances[cls]

#########################################################################################
# For creating nice looking errors on results

class ShorthandFormatter(string.Formatter):

  def format_field(self, value, format_spec):
    if isinstance(value, uncertainties.UFloat):
      return value.format(format_spec+'S')  # Shorthand option added
    # Special formatting for other types can be added here (floats, etc.)
    else:
      # Usual formatting:
      return super(ShorthandFormatter, self).format_field(value, format_spec)


def nice_value(value, error):
  if error > 0.:
    frmtr = ShorthandFormatter()
    return frmtr.format(f"{{0:.{ERR_PREC}u}}", uncertainties.ufloat(value, error))
  else:
    logging.warning("Error <= 0. passed to nice_value...")
    return str(value)

#########################################################################################
# Some utility functions

def str_to_file(st):
  st = regex.sub('[\[\]()=:,]', '', st.strip()).replace(' ', '-')
  return st

def datestr():
  return str(datetime.datetime.now().replace(microsecond=0)).replace(' ', '_')

def xmltostr(xml):
  doc = minidom.Document()
  declaration = doc.toxml()

  xmlstr = minidom.parseString(ET.tostring(xml, encoding='utf-8')).toprettyxml(indent="  ")[len(declaration)+1:]
  return xmlstr

def check_extra_keys(extra_dict, name):
  for extra_key in extra_dict.keys():
    logging.error(f"Unrecognized key '{extra_key}' in {name}")

def updateOption(options, option, option_class):
  if option in options:
    try:
      options[option] = option_class(options[option])
    except ValueError:
      logging.error(f"Invalid '{option}' for {option_class.__name__}")


#########################################################################################
# writing samplings to file

def write_samplings(obs_handler, obs, filename):
  os.makedirs(os.path.dirname(filename), exist_ok=True)
  f_handler = open(filename, 'w')
  obs_handler.setSamplingBegin()
  while not obs_handler.isSamplingEnd():
    sampling = obs_handler.getCurrentSamplingValue(obs)
    f_handler.write(f"{sampling}\n")
    obs_handler.setSamplingNext()

  f_handler.close()
  logging.info(f"Wrote samplings file: {filename}")

#########################################################################################
# Get samplings in numpy format

def get_samplings(obs_handler, obs):
  obs_handler.setSamplingBegin()
  samples = list()
  while not obs_handler.isSamplingEnd():
    samples.append(obs_handler.getCurrentSamplingValue(obs))
    obs_handler.setSamplingNext()

  return np.array(samples)


#########################################################################################
# For manipulating sigmond observables

# @ADH - Why must I also return obs_get_handler?
def get_obs_handlers(data_files, bins_info, sampling_info):
  xml = sigmond.XMLHandler()
  xml.set_from_string(ET.tostring(data_files.xml()))
  obs_get_handler = sigmond.MCObsGetHandler(xml, bins_info, sampling_info)
  obs_handler = sigmond.MCObsHandler(obs_get_handler, False)
  return obs_handler, obs_get_handler

def _get_free_obs(obs_handler, obs_name):
  obs_name = f"{obs_name}"
  obs_ind = 100000
  obs_info = sigmond.MCObsInfo(obs_name, obs_ind)
  while obs_handler.queryBins(obs_info) or obs_handler.queryFullAndSamplings(obs_info):
    obs_ind += 1
    obs_info = sigmond.MCObsInfo(obs_name, obs_ind)

  return obs_info

# @ADH - merge the following two functions since they are nearly identical
def get_absolute_energy(obs_handler, energy_shift_obs_info, non_interacting_level,
                        spatial_sites):

  boosted_obs = [energy_shift_obs_info]
  coeffs = [1.]
  for scat_obs, psq in non_interacting_level:
    coeffs.append(1.)

    boosted_ob = _get_free_obs(obs_handler, scat_obs.getObsName())
    if boosted_ob is None:
      return None

    psqfactor = psq * (2.*math.pi/spatial_sites)**2
    try:
      sigmond.doBoostBySamplings(obs_handler, scat_obs, psqfactor, boosted_ob)
    except RuntimeError:
      logging.critical(f"Could not get samplings for\n{scat_obs!s}")

    boosted_obs.append(boosted_ob)

  absolute_obs = _get_free_obs(obs_handler, energy_shift_obs_info.getObsName())
  if absolute_obs is None:
    return None

  try:
    absolute_obs = linear_superposition_obs(obs_handler, boosted_obs, coeffs)
  except RuntimeError:
    err_str = "Could not get samplings for one of\n"
    for obs_info in boosted_obs:
      err_str += str(obs_info)

    logging.critical(err_str)

  return absolute_obs

def get_energy_difference(obs_handler, abs_energy_obs_info, non_interacting_level,
                          spatial_sites):

  boosted_obs = [abs_energy_obs_info]
  coeffs = [1.]
  for scat_obs, psq in non_interacting_level:
    coeffs.append(-1.)

    boosted_ob = _get_free_obs(obs_handler, scat_obs.getObsName())
    if boosted_ob is None:
      return None

    psqfactor = psq * (2.*math.pi/spatial_sites)**2
    try:
      sigmond.doBoostBySamplings(obs_handler, scat_obs, psqfactor, boosted_ob)
    except RuntimeError:
      logging.critical(f"Could not get samplings for\n{scat_obs!s}")

    boosted_obs.append(boosted_ob)

  energy_diff_obs = _get_free_obs(obs_handler, abs_energy_obs_info.getObsName())
  if energy_diff_obs is None:
    return None

  try:
    diff_obs = linear_superposition_obs(obs_handler, boosted_obs, coeffs)
  except RuntimeError:
    err_str = "Could not get samplings for one of\n"
    for obs_info in boosted_obs:
      err_str += str(obs_info)

    logging.critical(err_str)


  return diff_obs

def ratio_obs(obs_handler, num_obs_info, den_obs_info):
  ratio_obs = _get_free_obs(obs_handler, num_obs_info.getObsName())
  if ratio_obs is None:
    return None

  try:
    sigmond.doRatioBySamplings(obs_handler, num_obs_info, den_obs_info, ratio_obs)
  except RuntimeError:
    logging.critical(f"Could not get samplings for\n{num_obs_info!s}or\n{den_obs_info!s}")

  return ratio_obs


def boost_obs_to_cm(obs_handler, energy_obs, psq, spatial_sites):
  cm_obs = _get_free_obs(obs_handler, energy_obs.getObsName())

  if cm_obs is None:
    return None
 
  psqfactor = -psq * (2.*math.pi/spatial_sites)**2
  try:
    sigmond.doBoostBySamplings(obs_handler, energy_obs, psqfactor, cm_obs)
  except RuntimeError:
    logging.crtical(f"Could not get samplings for\n{energy_obs!s}")

  return cm_obs

def boost_obs(obs_handler, energy_obs, psq, spatial_sites):
  boost_obs = _get_free_obs(obs_handler, energy_obs.getObsName())

  if boost_obs is None:
    return None
 
  psqfactor = psq * (2.*math.pi/spatial_sites)**2
  try:
    sigmond.doBoostBySamplings(obs_handler, energy_obs, psqfactor, boost_obs)
  except RuntimeError:
    logging.crtical(f"Could not get samplings for\n{energy_obs!s}")

  return boost_obs

def linear_superposition_obs(obs_handler, obs_infos, coeffs):
  result_obs = _get_free_obs(obs_handler, obs_infos[0].getObsName())
  if result_obs is None:
    return None

  try:
    sigmond.doLinearSuperpositionBySamplings(obs_handler, obs_infos, coeffs, result_obs)
  except RuntimeError:
    err_str = "Could not get samplings for one of\n"
    for obs_info in obs_infos:
      err_str += str(obs_info)

    logging.critical(err_str)

  return result_obs

#########################################################################################
# PDF documents
def create_doc(title, with_tikz=False):
  logging.info("Creating PDF...")
  doc = pylatex.Document(geometry_options={'margin': '1.5cm'})
  doc.packages.append(pylatex.Package('hyperref'))
  doc.packages.append(pylatex.Package('amssymb'))
  doc.packages.append(pylatex.Package('amsmath'))
  doc.packages.append(pylatex.Package('float'))
  doc.packages.append(pylatex.Package('caption', {'labelformat': 'empty', 'justification': 'centering'}))
  doc.packages.append(pylatex.Package('todonotes'))
  doc.packages.append(pylatex.NoEscape(r"\usepackage{longtable}[=v4.13]"))

  # I don't think these are needed
  #doc.packages.append(pylatex.Package('needspace'))
  #doc.packages.append(pylatex.Package('pgfplots'))

  if with_tikz:
    doc.packages.append(pylatex.Package('tikz'))
    doc.preamble.append(pylatex.Command('usetikzlibrary', 'decorations.markings'))
    doc.preamble.append(pylatex.Command('usetikzlibrary', 'fit'))
    doc.preamble.append(pylatex.Command('usetikzlibrary', 'plotmarks'))

  doc.preamble.append(pylatex.Command('title', title))
  doc.preamble.append(pylatex.Command('date', ''))

  doc.append(pylatex.NoEscape(r"\maketitle"))
  doc.append(pylatex.NoEscape(r"\tableofcontents"))
  doc.append(pylatex.NoEscape(r"\newpage"))
  doc.append(pylatex.NoEscape(r"\captionsetup[subfigure]{labelformat=empty}"))

  return doc

def create_tikz_doc():
  logging.info("Creating PDF...")
  doc = pylatex.Document(documentclass='standalone', document_options=['crop','tikz'])
  doc.packages.append(pylatex.Package('tikz'))
  doc.preamble.append(pylatex.Command('usetikzlibrary', 'decorations.markings'))
  doc.preamble.append(pylatex.Command('usetikzlibrary', 'fit'))
  doc.preamble.append(pylatex.Command('usetikzlibrary', 'plotmarks'))
  return doc

def add_correlator(doc, task_handler, correlator, name, obs_handler):

  operator_src = operator_info.operator.Operator(correlator.getSource())
  subtractvev = task_handler.subtractvev and operator_src.channel.vev
  if correlator.isSinkSourceSame():
    left_pdf_file = task_handler.correlator_plotfile(correlator, name, extension=PlotExtension.pdf)
    right_pdf_file = task_handler.energy_plotfile(operator_src, name, extension=PlotExtension.pdf)

  else:
    left_pdf_file = task_handler.correlator_plotfile(
        correlator, name, complex_arg=sigmond.ComplexArg.RealPart, extension=PlotExtension.pdf)
    right_pdf_file = task_handler.correlator_plotfile(
        correlator, name, complex_arg=sigmond.ComplexArg.ImaginaryPart, extension=PlotExtension.pdf)

  if correlator.isSinkSourceSame():
    left_estimates = sigmond.getCorrelatorEstimates(
        obs_handler, correlator, task_handler.hermitian, subtractvev, sigmond.ComplexArg.RealPart,
        task_handler.sampling_mode)
    right_estimates = sigmond.getEffectiveEnergy(
        obs_handler, correlator, task_handler.hermitian, subtractvev, sigmond.ComplexArg.RealPart,
        task_handler.sampling_mode, task_handler.plot_info.timestep, 
        task_handler.plot_info.eff_energy_type.value, 0.)
  else:
    left_estimates = sigmond.getCorrelatorEstimates(
        obs_handler, correlator, task_handler.hermitian, subtractvev, sigmond.ComplexArg.RealPart,
        task_handler.sampling_mode)
    right_estimates = sigmond.getCorrelatorEstimates(
        obs_handler, correlator, task_handler.hermitian, subtractvev, sigmond.ComplexArg.ImaginaryPart,
        task_handler.sampling_mode)

  if correlator.isSinkSourceSame():
    doc.append(pylatex.NoEscape(rf"Score: {score(left_estimates)}"))

  with doc.create(pylatex.Figure(position='H')):
    with doc.create(pylatex.SubFigure(
        position='b', width=pylatex.NoEscape(r'0.5\linewidth'))) as left_fig:
      add_image(left_fig, task_handler.results_dir, left_pdf_file, width="1.0")
    with doc.create(pylatex.SubFigure(
        position='b', width=pylatex.NoEscape(r'0.5\linewidth'))) as right_fig:
      add_image(right_fig, task_handler.results_dir, right_pdf_file, width="1.0")


  if correlator.isSinkSourceSame():
    header_row = [
        pylatex.NoEscape(r"$t$"),
        pylatex.NoEscape(r"$C(t)$"),
        pylatex.NoEscape(r"$\delta C(t)$"), 
    ]
    if task_handler.plot_info.eff_energy_type.value < 2:
      header_row.extend([
          pylatex.NoEscape(rf"$a_t E_{{\rm eff}} (t + {task_handler.plot_info.timestep}/2)$"),
          pylatex.NoEscape(rf"$\delta a_t E_{{\rm eff}} (t + {task_handler.plot_info.timestep}/2)$"),])
    else:
      header_row.extend([
          pylatex.NoEscape(r"$a_t E_{\rm eff} (t)$"),
          pylatex.NoEscape(r"$\delta a_t E_{\rm eff} (t)$"),])

  else:
      header_row = [
          pylatex.NoEscape(r"$t$"),
          pylatex.NoEscape(r"$Re C(t)$"),
          pylatex.NoEscape(r"$\delta Re C(t)$"), 
          pylatex.NoEscape(r"$Im C(t)$"),
          pylatex.NoEscape(r"$\delta Im C(t)$")
      ]


  with doc.create(pylatex.Center()) as centered:
    with centered.create(pylatex.LongTabu("X[c] X[2,c] X[2,c] X[2,c] X[2,c]", 
                         to=r"\linewidth")) as data_table:
      data_table.add_row(header_row, mapper=[pylatex.utils.bold])
      data_table.add_hline()
      data_table.end_table_header()
      for t in sorted(left_estimates.keys()):
        left_value = left_estimates[t].getFullEstimate()
        left_error = left_estimates[t].getSymmetricError()
        left_est = nice_value(left_value, left_error)
        left_rel_error = round(left_estimates[t].getRelativeError(), 4)
        t_right = t
        if correlator.isSinkSourceSame() and task_handler.plot_info.eff_energy_type.value < 2:
          t_right = t + 0.5*task_handler.plot_info.timestep

        if t_right in right_estimates:
          right_value = right_estimates[t_right].getFullEstimate()
          right_error = right_estimates[t_right].getSymmetricError()
          right_est = nice_value(right_value, right_error)
          right_rel_error = round(right_estimates[t_right].getRelativeError(), 4)
        else:
          right_est = ""
          right_rel_error = ""

        row = [int(t), left_est, left_rel_error, right_est, right_rel_error]
        data_table.add_row(row)

  doc.append(pylatex.NoEscape(r"\newpage"))


def score(corr_estimates):
  _score = 0.
  for t, estimate in corr_estimates.items():
    if estimate.isStatisticallyZero() or estimate.getFullEstimate() < 0.:
      continue

    _score += 1. / estimate.getRelativeError()

  return _score


def add_image(figure, rel_dir, pdf_file, width="1.0", caption="", view=True):
  rel_pdf_file = os.path.relpath(pdf_file, rel_dir)
  rel_grace_file = os.path.splitext(rel_pdf_file)[0] + ".agr"

  width = f"{width}\linewidth"
  placement = "\centering"
  if os.path.isfile(pdf_file):
    if caption and view:
      caption = f"{caption} \\newline \href{{run:{rel_grace_file}}}{{view}}"
    elif view:
      caption = f"\href{{run:{rel_grace_file}}}{{view}}"

    figure.add_image(rel_pdf_file, width=pylatex.NoEscape(width), placement=pylatex.NoEscape(placement))
    if caption:
      figure.add_caption(pylatex.NoEscape(caption))

  else:
    logging.warning("Could not find file: {}".format(pdf_file))
    missing_figure_tex = r"\missingfigure[figwidth=" + width + r"]{missing " + os.path.basename(pdf_file).replace("_", r"\_") + r"}"
    figure.append(pylatex.NoEscape(missing_figure_tex))


def compile_pdf(doc, filename, compiler=None):
  logging.info(f"Compiling LaTeX: {filename}.tex")
  doc.generate_pdf(filename, clean=True, clean_tex=False, compiler=compiler, compiler_args=['-f'])
  logging.info(f"Created PDF: {filename}.pdf")

def write_tikz(tikz, filename):
  with open(filename, 'w', encoding='utf-8') as f_hand:
    tikz.dump(f_hand)


#########################################################################################
# Grace file converters

def dirGrace2pdf(dirname, margins='0'):
  grace2pdf_sh = os.path.join(os.path.dirname(__file__), 'grace2pdf.sh')
  logging.info(f"Converting grace plots to PDF: {dirname}")
  cwd = os.getcwd()
  os.chdir(dirname)
  try:
    ps = list()
    for filename in os.listdir('.'):
      base_filename, ext = os.path.splitext(filename)
      if ext == ".agr":
        new_filename = f"{base_filename}.pdf"

        if not os.path.isfile(new_filename) or os.path.getmtime(filename) > os.path.getmtime(new_filename):
          command_list = [grace2pdf_sh, filename, margins]

          try:
            p = subprocess.Popen(command_list, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

          except BlockingIOError:
            for p in ps:
              p.wait()
            ps = list()
            p = subprocess.Popen(command_list, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
          ps.append(p)

    for p in ps:
      p.wait()

  except FileNotFoundError:
    print("No such directory: {}".format(dirname))

  os.chdir(cwd)
  logging.info("Finished converting grace plots")

def grace2pdf(grace_filename, margins='0'):
  base_filename, ext = os.path.splitext(grace_filename)
  if ext != ".agr":
    logging.error("Grace2pdf 'grace_filename' must have .agr extension")
  cwd = os.getcwd()
  os.chdir(os.path.dirname(grace_filename))

  pdf_filename = f"{base_filename}.pdf"
  if not os.path.isfile(pdf_filename) or os.path.getmtime(grace_filename) > os.path.getmtime(pdf_filename):
    grace2pdf_sh = os.path.join(os.path.dirname(__file__), 'grace2pdf.sh')
    command_list = [grace2pdf_sh, grace_filename, margins]
    process = subprocess.Popen(command_list, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    process.wait()
    if process.returncode != 0:
      logging.error("Failed to convert grace file to pdf '{grace_filename}'")

  os.chdir(cwd)
