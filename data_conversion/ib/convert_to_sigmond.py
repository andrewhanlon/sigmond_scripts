#!/usr/bin/env python

import os
import sys

import xml.etree.ElementTree as ET
import h5py
import numpy as np
import tqdm

import argparse

import defs

import sigmond

FORCE_HERM = True

base_data_dir = "/disk2/research/data/ib"
proj_data_dir = "/disk1/research/projects/ib/data"
output_dir = "/disk1/research/projects/ib/data"

ensembles_to_do = ['A654']

def main():
  parser = argparse.ArgumentParser(description="Convert data")
  parser.add_argument("-a", "--ama", action="store_true",
                      help="Specify if ama needs to be done")
  parser.add_argument("-m", "--merge", action="store_true",
                      help="Specify if merging should be done")
  parser.add_argument("-s", "--sigmond", action="store_true",
                      help="Specify if conversion to sigmond should be done")

  args = parser.parse_args()

  for ensemble in defs.ensembles:
    if ensemble.name not in ensembles_to_do:
      continue

    print(f"Processing ensemble {ensemble.name}")
    smearings, data_infos = get_data_infos(ensemble.name)

    if args.ama:
      do_ama(ensemble, smearings, data_infos)

    if args.merge:
      do_merge(ensemble, data_infos)

    if args.sigmond:
      convert_to_sigmond(ensemble, smearings, data_infos)

def do_ama(ensemble, smearings, data_infos):
  raw_data_dir = defs.raw_data_dir(base_data_dir, ensemble.name)
  ave_data_dir = defs.ave_data_dir(base_data_dir, ensemble.name)
  os.makedirs(ave_data_dir, exist_ok=True)

  for config in tqdm.tqdm(defs.configs[ensemble.name]):
    sloppy_file = os.path.join(raw_data_dir, defs.raw_data_file(ensemble.name, config, False))
    sloppy_fh = h5py.File(sloppy_file, mode='r')

    sloppy_smearings = sloppy_fh["/axes/quark_smearing"][:]
    if sloppy_smearings != smearings:
      print("Mismatched smearings")
      sys.exit()

    sloppy_sources = sloppy_fh["/axes/source_position"][:]

    exact_file = os.path.join(raw_data_dir, defs.raw_data_file(ensemble.name, config, True))
    exact_fh = h5py.File(exact_file, mode='r')

    exact_smearings = exact_fh["/axes/quark_smearing"][:]
    if exact_smearings != smearings:
      print("Mismatched smearings")
      sys.exit()

    exact_sources = exact_fh["/axes/source_position"][:]
    exact_sources_sloppy_indices = list()
    for source in exact_sources:
      try:
        exact_sources_sloppy_indices.append(list(sloppy_sources).index(source))
      except ValueError:
        print(f"{source} missing from sloppy for config {config}")


    for channel, data_info in data_infos.items():
      ave_file = os.path.join(ave_data_dir, defs.ave_data_file(ensemble.name, channel, config))
      ave_fh = h5py.File(ave_file, mode='w')

      particle_type = defs.particle_types[channel.flavor]

      for contribution, diagrams in data_info.contributions.items():
        for diagram in diagrams:
          data_key = f"/data/{channel!r}/{contribution}/{diagram}"
          sloppy_data = sloppy_fh[data_key][:]
          for src_i, source in enumerate(sloppy_sources):
            tsrc = defs.get_source_time(source)
            if particle_type is defs.ParticleType.FERMION:
              sloppy_data[src_i,:,:,:tsrc] = -sloppy_data[src_i,:,:,:tsrc]

            sloppy_data[src_i,:] = np.roll(sloppy_data[src_i,:], -tsrc, axis=3)


          exact_data = exact_fh[data_key][:]
          for src_i, source in enumerate(exact_sources):
            tsrc = defs.get_source_time(source)
            if particle_type is defs.ParticleType.FERMION:
              exact_data[src_i,:,:,:tsrc] = -exact_data[src_i,:,:,:tsrc]

            exact_data[src_i,:] = np.roll(exact_data[src_i,:], -tsrc, axis=3)

          if sloppy_data.shape[1:] != exact_data.shape[1:]:
            print("Mismatched data shapes")
            sys.exit()

          bias_correction_data = np.zeros(exact_data.shape, dtype=np.complex128)

          for exact_source_i, sloppy_source_i in enumerate(exact_sources_sloppy_indices):
            bias_correction_data[exact_source_i,:] = exact_data[exact_source_i,:] - sloppy_data[sloppy_source_i,:]

          final_data = np.mean(bias_correction_data, axis=0) + np.mean(sloppy_data, axis=0)

          ave_fh.create_dataset(f"/{contribution}/{diagram}", data=final_data)

      ave_fh.close()

    sloppy_fh.close()
    exact_fh.close()


def do_merge(ensemble, data_infos):
  ave_data_dir = defs.ave_data_dir(base_data_dir, ensemble.name)
  merged_data_dir = defs.merged_data_dir(base_data_dir, ensemble.name)
  os.makedirs(merged_data_dir, exist_ok=True)

  for channel, data_info in tqdm.tqdm(data_infos.items()):
    merged_file = os.path.join(merged_data_dir, defs.merged_data_file(ensemble.name, channel))
    merged_fh = h5py.File(merged_file, mode='w')

    ave_fhs = list()
    for config in defs.configs[ensemble.name]:
      ave_file = os.path.join(ave_data_dir, defs.ave_data_file(ensemble.name, channel, config))
      ave_fhs.append(h5py.File(ave_file, mode='r'))
    
    for contribution, diagrams in tqdm.tqdm(data_info.contributions.items(), leave=False):
      for diagram in tqdm.tqdm(diagrams, leave=False):
        data_key = f"/{contribution}/{diagram}"
        data_list = list()
        for ave_fh in tqdm.tqdm(ave_fhs, leave=False):
          data_list.append(ave_fh[data_key][:])

        merged_fh.create_dataset(data_key, data=np.stack(data_list))

    for ave_fh in ave_fhs:
      ave_fh.close()

    merged_fh.close()


def convert_to_sigmond(ensemble, smearings, data_infos):
  ensemble_info = sigmond.MCEnsembleInfo(f"cls21_{ensemble.name}", 'ensembles.xml')
  mcobs_xml = ET.Element("MCObservables")
  mcobs_xml_handler = sigmond.XMLHandler()
  mcobs_xml_handler.set_from_string(ET.tostring(mcobs_xml))
  sampling_info = sigmond.MCSamplingInfo()
  bins_info = sigmond.MCBinsInfo(ensemble_info)
  bins_info.addOmissions(defs.omissions[ensemble.name])
  obs_get_handler = sigmond.MCObsGetHandler(mcobs_xml_handler, bins_info, sampling_info)
  obs_handler = sigmond.MCObsHandler(obs_get_handler, False)

  merged_data_dir = defs.merged_data_dir(base_data_dir, ensemble.name)
  sigmond_data_dir = defs.sigmond_data_dir(proj_data_dir, ensemble.name)
  os.makedirs(sigmond_data_dir, exist_ok=True)
  smearings_file = os.path.join(sigmond_data_dir, "smearings.txt")
  with open(smearings_file, 'w') as fh:
    for smearing_id, smearing in enumerate(smearings):
      os.makedirs(os.path.join(sigmond_data_dir, f"smearing_{smearing_id}"), exist_ok=True)
      fh.write(f"{smearing_id}: {smearing.decode()}\n")

  for channel, data_info in tqdm.tqdm(data_infos.items()):
    merged_file = os.path.join(merged_data_dir, defs.merged_data_file(ensemble.name, channel))
    merged_fh = h5py.File(merged_file, mode='r')

    irrep = channel.irrep
    irrep_row = defs.irrep_rows[irrep][channel.irrep_row]
    flavor = f"{defs.flavors[channel.flavor][0]},{defs.flavors[channel.flavor][1]}"

    for smearing_id in range(len(smearings)):
      for contribution, diagrams in tqdm.tqdm(data_info.contributions.items(), leave=False):
        sigmond_subdata_dir = os.path.join(sigmond_data_dir, f"smearing_{smearing_id}", channel.flavor, contribution)
        os.makedirs(sigmond_subdata_dir, exist_ok=True)

        sigmond_file = os.path.join(sigmond_subdata_dir, defs.sigmond_data_file(ensemble.name, channel, contribution))
        obs_infos = set()
        for diagram in diagrams:
          data_key = f"/{contribution}/{diagram}"
          data = merged_fh[data_key]

          name = f"{channel!s}_{defs.shorten_contribution(contribution)}_{defs.shorten_diagram(diagram)}"
          if len(name) > 26:
            print("Name too long")
            print(name)
            sys.exit()

          for src_i, src_id in enumerate(data_info.operator_ids):
            src_op_str = f"P=(0,0,0) {irrep}_{irrep_row} flavor={flavor} {name} {src_id}"
            src_op = sigmond.OperatorInfo(src_op_str, sigmond.OpKind.GenIrrep)

            for snk_i, snk_id in enumerate(data_info.operator_ids):
              snk_op_str = f"P=(0,0,0) {irrep}_{irrep_row} flavor={flavor} {name} {snk_id}"
              snk_op = sigmond.OperatorInfo(snk_op_str, sigmond.OpKind.GenIrrep)

              for tsep in range(ensemble.Nt):
                re_obs_info = sigmond.MCObsInfo(snk_op, src_op, tsep, True, sigmond.ComplexArg.RealPart, False)
                re_data = data[:,smearing_id,snk_i,src_i,tsep].real
                obs_handler.putBins(re_obs_info, sigmond.RVector(re_data))
                obs_infos.add(re_obs_info)

                im_obs_info = sigmond.MCObsInfo(snk_op, src_op, tsep, True, sigmond.ComplexArg.ImaginaryPart, False)
                im_data = data[:,smearing_id,snk_i,src_i,tsep].imag
                obs_handler.putBins(im_obs_info, sigmond.RVector(im_data))
                obs_infos.add(im_obs_info)

        xml_out = sigmond.XMLHandler("output", "")
        obs_handler.writeBinsToFile(obs_infos, sigmond_file, xml_out, sigmond.WriteMode.Protect)



    merged_fh.close()


def get_data_infos(ensemble_name):
  first_config = defs.configs[ensemble_name][0]
  raw_data_dir = defs.raw_data_dir(base_data_dir, ensemble_name)
  data_file = os.path.join(raw_data_dir, defs.raw_data_file(ensemble_name, first_config, True))

  data_fh = h5py.File(data_file, 'r')

  smearings = data_fh["/axes/quark_smearing"][:]

  data_infos = dict()
  for channel_str, channel_group in data_fh['/data'].items():
    channel = defs.get_channel_from_str(channel_str)
    op_ids = channel_group.attrs['operator_ids']
    contributions = dict()
    for contribution_str, contribution_group in channel_group.items():
      contributions[contribution_str] = list()
      for diagram_str, dataset in contribution_group.items():
        contributions[contribution_str].append(diagram_str)

    data_infos[channel] = defs.DataInfo(op_ids, contributions)

  data_fh.close()

  return smearings, data_infos


if __name__ == "__main__":
  main()
