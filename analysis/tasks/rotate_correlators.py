import os

import logging
from sortedcontainers import SortedSet
from enum import Enum
import pylatex

import sigmond

import tasks.task
import utils.util as util
import sigmond_info.sigmond_info
import sigmond_info.sigmond_log
import sigmond_info.sigmond_input
import operator_info.operator
import operator_info.operator_set


class RotateCorrelators(tasks.task.Task):
  task_type = "rotate_corrs"

  def initialize(self, operator_bases, **extra_options):
    """sets the needed parameters for the rotate correlators task

    Args:
      operator_bases ([RotatedOperatorSet]): A list of rotated
          operator sets.
      **rotate_mode (sigmond_info.RotateMode): rotate mode to use,
          either bins, samplings, samplings_all, samplings_unsubt
      **show_transformation (bool): whether or not the transformation
          used for the rotation should be shown
      **negative_eigenvalue_alarm (float): see sigmond docs
      **subtractvev (bool): whether operators with non-zero VEVs
          should have the VEV subtracted
      **remove_off_diag (bool): (default: true) specifies whether the off diagonal
          correlators should be printed
      **diagonly_time (int): (default: 0) specifies what time only diagonal correlators
          are calculated
      **plot_info (PlotInfo): info about the plot
      **write_to_file (bool): write correlator data to csv files
          in results directory and the samplings to hdf5.
              Should only be ran with max_simultaneous: 1.
    """
    self.hermitian = True
    self.operator_bases = operator_bases
    self.remove_off_diag = extra_options.pop('remove_off_diag', False)
    self.off_diagonal = not self.remove_off_diag
    self.diagonly_time = extra_options.pop('diagonly_time', 0)

    self.rotate_mode = extra_options.pop('rotate_mode', sigmond_info.sigmond_info.RotateMode.samplings_all)

    self.show_transformation = extra_options.pop('show_transformation', True)

    self.neg_eig_alarm = extra_options.pop('negative_eigenvalue_alarm', -0.01)
    self.subtractvev = extra_options.pop('subtractvev', True)

    self.plot_info = extra_options.pop('plot_info', sigmond_info.sigmond_info.PlotInfo())
    self.write_to_file = extra_options.pop('write_to_file', False)

    util.check_extra_keys(extra_options, self.task_name)

  def readConfig(self, **task_options):
    """readConfig function for RotateCorrelators

    YAML:
      # optional (defaults shown)
      show_transformation: true
      negative_eigenvalue_alarm: -0.10
      subtractvev: true
      remove_off_diag: false
      diagonly_time: 10
      rotate_mode: bins # or samplings_all, samplings_unsubt, samplings
      write_to_file: true # optional, if flag present writes rotated 
              correlators to csv files and the samplings to hdf5.
              Should only be ran with max_simultaneous: 1.

      plot_info:
        corrname: standard
        symbol_color: blue
        symbol_type: circle
        effective_energy_type: time_forward
        timestep: 3
        rescale: 1.0

      operator_bases:
        - name: basis_name
          pivot_info:             # depends on pivot type
            pivot_type: single_pivot
            norm_time: 10
            metric_time: 18
            diagonalize_time: 22
            max_condition_number: 100
          optimized_operators:
            - basis_name_of_opops (uses name to retrieve information from logfile in same project directory)
          operators:       # optional (if basis_name already defined)
            - op_string1
            - op_string2
            - op_string3
            ...

        - name: other_basis_name
          ...
        ...
    """

    try:
      operator_bases = list()
      for operator_basis in task_options.pop('operator_bases'):
        operator_bases.append(operator_info.operator_set.getOperatorSet(operator_basis))

    except KeyError as err:
      logging.error(f"Invalid key '{err}' in task '{self.task_name}'")
    except AttributeError as err:
      logging.error(f"Rotated basis needs a pivot tag")

    task_options['plot_info'] = sigmond_info.sigmond_info.PlotInfo.createFromConfig(task_options)
    task_options['rotate_mode'] = sigmond_info.sigmond_info.RotateMode(task_options.pop('rotate_mode', 'samplings_all'))
    self.initialize(operator_bases, **task_options)

  def energy_plotdir(self, rotated_basis):
    plotdir = os.path.join(self.project_dir, self.relative_plotdir, "rotated_eff_energies",
                           rotated_basis.name, repr(rotated_basis.pivot_info))

    os.makedirs(plotdir, exist_ok=True)
    return plotdir

  def correlator_plotdir(self, rotated_basis):
    plotdir = os.path.join(self.project_dir, self.relative_plotdir, "rotated_correlators",
                           rotated_basis.name, repr(rotated_basis.pivot_info))
    os.makedirs(plotdir, exist_ok=True)
    return plotdir

  def correlator_plotstub(self, rotated_basis):
    plotdir = self.correlator_plotdir(rotated_basis)
    plotstub = f"corr_{rotated_basis.name}_rotated"
    return os.path.join(plotdir, plotstub)

  def energy_plotstub(self, rotated_basis):
    plotdir = self.energy_plotdir(rotated_basis)
    plotstub = f"eff-energy_{rotated_basis.name}_rotated"
    return os.path.join(plotdir, plotstub)

  def correlator_plotfile(self, correlator, rotated_basis, **options):
    """get a correlator plotfile
      
      correlator (sigmond.CorrelatorInfo): the correlator to plot
      name (str): ...
      **complex_arg (ComplexArg): Required if the correlator
          is not diagonal. If it is diagonal, this is
          ignored
      **extension (PlotExtension): 
    """
    extension = options.pop('extension', util.PlotExtension.grace)
    complex_arg = options.pop('complex_arg', sigmond.ComplexArg.RealPart)
    complex_str = "Re" if complex_arg == sigmond.ComplexArg.RealPart else "Im"

    src_level = operator_info.operator.Operator(correlator.getSource()).level
    snk_level = operator_info.operator.Operator(correlator.getSource()).level
    plotdir = self.correlator_plotdir(rotated_basis)
    plotstub = self.correlator_plotstub(rotated_basis)
    plotfile = f"{plotstub}_{src_level}-{snk_level}_{complex_str}.{extension.value}"
    return os.path.join(plotdir, plotfile)

  def energy_plotfile(self, operator, rotated_basis, extension):
    level = operator.level
    plotdir = self.energy_plotdir(rotated_basis)
    plotstub = self.energy_plotstub(rotated_basis)
    plotfile = f"{plotstub}_{level}.{extension.value}"
    return os.path.join(plotdir, plotfile)

  def samplings_file_hdf5(self, irrep, mom):
    return os.path.join(self.results_dir, f"rotated_correlators.hdf5[{irrep}_PSQ{mom}]")

  def getSigmondInputs(self):
    sigmond_inputs = list()
    overwrite = True
    for operator_basis in self.operator_bases:
      logging.info(f"Working on basis {operator_basis.name}")

      if operator_basis.num_channels != 1:
        logging.info("  Basis must only have one channel")
        continue

      subtractvev = self.subtractvev and operator_basis.vev
      logging.info(f"  Subtract VEV = {subtractvev}")

      logging.info("  Operator basis:")
      op_infos = set()
      for operator in operator_basis.operators:
        logging.info(f"    {operator!s}")
        op_infos.add(operator.operator_info)

      corr_mat = sigmond.CorrelatorMatrixInfo(op_infos, True, subtractvev)

      mintime, maxtime = self.data_handler.getOperatorSetSmallestTRange(operator_basis)
      logging.info(f"  Time separations [{mintime},{maxtime}]")
      data_files = self.data_handler.getChannelDataFiles(operator_basis.channel)


      project_name = self.project_name(repr(operator_basis))
      logfile = self.logfile(repr(operator_basis))
      inputfile = self.inputfile(repr(operator_basis))
      sigmond_input = self.new_sigmond_input(project_name, inputfile, logfile, data_files)

      resulting_op = operator_basis.channel.getRotatedOp()
      pivot_file = self.data_handler.pivotfile(operator_basis)
      rotated_file = self.data_handler.rotated_datafile(operator_basis)

      corr_plotstub = self.correlator_plotstub(operator_basis)
      energy_plotstub = self.energy_plotstub(operator_basis)

      sigmond_input.doCorrMatrixRotation(
          operator_basis.pivot_info, self.rotate_mode, corr_mat, resulting_op, mintime, maxtime, 
          neg_eig_alarm=self.neg_eig_alarm, check_metric_errors=True, check_common_nullspace=True,
          show_transformation=self.show_transformation, pivot_filename=pivot_file,
          pivot_overwrite=True, rotated_corrs_filename=rotated_file,
          file_mode=sigmond.WriteMode.Overwrite, corr_plotstub=corr_plotstub, 
          energy_plotstub=energy_plotstub, eff_energy_type=self.plot_info.eff_energy_type, 
          timestep=self.plot_info.timestep, symbol_color=self.plot_info.symbol_color,
          symbol_type=self.plot_info.symbol_type, rescale=self.plot_info.rescale,
          remove_off_diag=self.remove_off_diag, diagonly_time=self.diagonly_time,
          improved_ops=operator_basis.verifyOptimizedOps(self.logdir))
        
      if self.write_to_file:
          if overwrite:
            write_mode = sigmond.WriteMode.Overwrite
            overwrite = False
          else:
            write_mode = sigmond.WriteMode.Update
            
          rot_op_infos = set()
          rot_op0 = operator_info.operator.Operator(resulting_op)
          rot_op_infos.add(rot_op0.operator_info)
          for i in range(1,len(op_infos)):
            this_compact_string = rot_op0.compact_str[:-1]
            this_compact_string+=f"{i}"
            this_rot_op = operator_info.operator.Operator.createFromCompact(this_compact_string);
            rot_op_infos.add(this_rot_op.operator_info)
#           
          rot_corr_mat = sigmond.CorrelatorMatrixInfo(rot_op_infos, True, subtractvev)
          sigmond_input.writeCorrMatToFile(
              self.samplings_file_hdf5(operator_basis.channel.irrep,operator_basis.channel.psq),rot_corr_mat,file_mode=write_mode,
          )

      sigmond_input.write()
      sigmond_inputs.append(sigmond_input)

    return sigmond_inputs

  def finalize(self):
    doc = util.create_doc(f"Rotated Correlators and Effective Energies: {self.task_name} - {self.ensemble_name}")
    
    #Sarah
    results_dir = self.results_dir
    os.makedirs(results_dir, exist_ok=True)
    if self.write_to_file:
        os.makedirs(os.path.join(results_dir,'estimates'), exist_ok=True)
    
    for operator_basis in self.operator_bases:
      logfile = self.logfile(repr(operator_basis))
      rotation_log = sigmond_info.sigmond_log.RotationLog(logfile)
      if rotation_log.failed:
        logging.warning(f"Rotation {operator_basis.name} failed: log located in {logfile}")
        continue
      
      corr_plotsdir = self.correlator_plotdir(operator_basis)
      energy_plotsdir = self.energy_plotdir(operator_basis)
      util.dirGrace2pdf(corr_plotsdir)
      util.dirGrace2pdf(energy_plotsdir)

      data_files = self.data_handler.getRotatedDataFiles(operator_basis)
      obs_handler, _ = util.get_obs_handlers(data_files, self.bins_info, self.sampling_info)

      with doc.create(pylatex.Section(f"{operator_basis.channel!s} - {operator_basis.name}")):
        with doc.create(pylatex.Subsection("Rotation Info")):
          with doc.create(pylatex.Center()) as centered:
            with centered.create(
                  pylatex.LongTabu("X[c]|X[c]|X[c]|X[c]|X[c]|X[3,c]|X[3,c]|X[3,c]|X[3,c]|X[3,c]",
                                 to=r"\linewidth")) as param_table:
              if(operator_basis.pivot_info.pivot_type==sigmond_info.sigmond_info.PivotType.RollingPivot):
                  header_row = [
                      pylatex.NoEscape(r"$N_{op}$"),
                      pylatex.NoEscape(r"$N_{\text{d}}$"),
                      pylatex.NoEscape(r"$\tau_N$"),
                      pylatex.NoEscape(r"$\tau_0$"),
                      pylatex.NoEscape(r"$\tau_Z$"),
                      pylatex.NoEscape(r"$\xi_{cn}$ (max)"),
                      pylatex.NoEscape(r"$\xi_{cn}^C$ (input)"),
                      pylatex.NoEscape(r"$\xi_{cn}^C$ (retain)"),
                      pylatex.NoEscape(r"$\xi_{cn}^G$ (input)"),
                      pylatex.NoEscape(r"$\xi_{cn}^G$ (retain)"),
                  ]
              else:
                  header_row = [
                      pylatex.NoEscape(r"$N_{op}$"),
                      pylatex.NoEscape(r"$N_{\text{d}}$"),
                      pylatex.NoEscape(r"$\tau_N$"),
                      pylatex.NoEscape(r"$\tau_0$"),
                      pylatex.NoEscape(r"$\tau_D$"),
                      pylatex.NoEscape(r"$\xi_{cn}$ (max)"),
                      pylatex.NoEscape(r"$\xi_{cn}^C$ (input)"),
                      pylatex.NoEscape(r"$\xi_{cn}^C$ (retain)"),
                      pylatex.NoEscape(r"$\xi_{cn}^G$ (input)"),
                      pylatex.NoEscape(r"$\xi_{cn}^G$ (retain)"),
                  ]
              param_table.add_row(header_row, mapper=[pylatex.utils.bold])
              param_table.add_hline()
              param_table.end_table_header()
              value_row = [
                  operator_basis.num_operators,
                  operator_basis.num_operators - rotation_log.number_levels,
                  operator_basis.pivot_info.norm_time,
                  operator_basis.pivot_info.metric_time,
                  operator_basis.pivot_info.diagonalize_time,
                  operator_basis.pivot_info.max_condition_number,
                  rotation_log.metric_condition(False),
                  rotation_log.metric_condition(True),
                  rotation_log.matrix_condition(False),
                  rotation_log.matrix_condition(True),
              ]
              param_table.add_row(value_row)

          doc.append(pylatex.NoEscape(r"\textbf{Metric Null Space Check:} " + \
                                      rotation_log.metric_null_space_message))
          
          if (not operator_basis._optimized_ops) and (operator_basis.pivot_info.pivot_type!=sigmond_info.sigmond_info.PivotType.RollingPivot):
            with doc.create(pylatex.Subsubsection("Input Operators")):
              with doc.create(pylatex.Center()) as centered:
                with centered.create(
                    pylatex.LongTabu("X[2,c] X[c] X[c]", row_height=1.5)) as op_table:
                  header_row = [
                      "Operator",
                      pylatex.NoEscape(r"$\delta C(\tau_0)$"),
                      pylatex.NoEscape(r"$\delta C(\tau_D)$")
                  ]
                  op_table.add_row(header_row, mapper=[pylatex.utils.bold])
                  op_table.add_hline()
                  op_table.end_table_header()
                  for op, errors in rotation_log.diagonal_correlator_errors.items():
                    row = [
                        op,
                        errors.metric,
                        errors.matrix,
                    ]
                    op_table.add_row
                    op_table.add_row(row)
                    
          if operator_basis.pivot_info.pivot_type!=sigmond_info.sigmond_info.PivotType.RollingPivot:
            with doc.create(pylatex.Subsubsection("Diagonal Deviations From Zero")):
                with doc.create(pylatex.Center()) as centered:
                  with centered.create(
                      pylatex.LongTabu("X[c] X[4,c] X[3,c] X[3,c] X[3,c] X[3,c] X[2,c]")) as deviation_table:
                    header_row = [
                        "time",
                        pylatex.NoEscape(r"$\delta 0_{max}$"),
                        pylatex.NoEscape(r"$\% > 1 \sigma$"),
                        pylatex.NoEscape(r"$\% > 2 \sigma$"),
                        pylatex.NoEscape(r"$\% > 3 \sigma$"),
                        pylatex.NoEscape(r"$\% > 4 \sigma$"),
                        "Status",
                    ]
                    deviation_table.add_row(header_row, mapper=[pylatex.utils.bold])
                    deviation_table.add_hline()
                    deviation_table.end_table_header()
                    for time, deviation in rotation_log.deviations_from_zero.items():
                      row = [
                          time,
                          deviation.max,
                          deviation.one,
                          deviation.two,
                          deviation.three,
                          deviation.four,
                          deviation.status,
                      ]
                      deviation_table.add_row(row)

            doc.append(pylatex.NoEscape(r"\newpage"))

        operators = self.data_handler.getRotatedOperators(operator_basis)
        #self.addPlotsToPDF(doc, data_files, operators, operator_basis)

        with doc.create(pylatex.Subsection("Rotated Diagonal Correlators/Effective Energies")):
          for operator in operators:
            with doc.create(pylatex.Subsubsection(str(operator))):
              corr = sigmond.CorrelatorInfo(operator.operator_info, operator.operator_info)
              util.add_correlator(doc, self, corr, operator_basis, obs_handler)
              if self.write_to_file:
                  self._write_diag_corr_to_csv(corr, operator_basis, obs_handler)

        doc.append(pylatex.NoEscape(r"\newpage"))

        off_diag_corrs = list()
        op_list = list(operators)
        for opsrc_i, op_src in enumerate(op_list):
          for op_snk in operators[opsrc_i+1:]:
            corr = sigmond.CorrelatorInfo(op_src.operator_info, op_snk.operator_info)
            off_diag_corrs.append(corr)

        with doc.create(pylatex.Subsection("Rotated Off-Diagonal Correlators")):
          for off_diag_corr in off_diag_corrs:
            with doc.create(pylatex.Subsubsection(off_diag_corr.corr_str())):
              util.add_correlator(doc, self, off_diag_corr, operator_basis, obs_handler)

    filename = os.path.join(results_dir, self.task_name)
    util.compile_pdf(doc, filename, self.latex_compiler)
    

  def _write_diag_corr_to_csv( self, correlator, name, obs_handler):
    operator_src = operator_info.operator.Operator(correlator.getSource())
    subtractvev = self.subtractvev and operator_src.channel.vev
    results_dir = self.results_dir
    if correlator.isSinkSourceSame():
        csv_file_stub = str(correlator.getSink()).replace(' ','_').replace('(0,0,0)','0').replace('=','')
        f = open(os.path.join(results_dir,'estimates',csv_file_stub+'.csv'),'w+')
        f.write('t,value,error\n')
        left_estimates = sigmond.getCorrelatorEstimates(
            obs_handler, correlator, self.hermitian, subtractvev, sigmond.ComplexArg.RealPart,
            self.sampling_mode)
        for t in sorted(left_estimates.keys()):
            left_value = left_estimates[t].getFullEstimate()
            left_error = left_estimates[t].getSymmetricError()
            f.write(f"{t},{left_value},{left_error}\n")
        f.close()
