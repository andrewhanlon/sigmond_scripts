import h5py
import numpy as np
import itertools
import xml.etree.ElementTree as ET

import sigmondbind as sigmond

eta_0_ops = [
    "eta P=(0,0,0) T1um_1 SS_0",
    "eta P=(0,0,0) T1um_2 SS_0",
    "eta P=(0,0,0) T1um_3 SS_0"
]

eta_1_ops = [
    "eta P=(0,0,0) T1um_1 SS_1",
    "eta P=(0,0,0) T1um_2 SS_1",
    "eta P=(0,0,0) T1um_3 SS_1"
]


tmin=1
tmax=40
tseps = list(range(tmin, tmax+1))

cfgmin = 1
cfgmax = 1100
cfginc = 1
cfgs = list(range(cfgmin, cfgmax+1, cfginc))

def main():
  ensemble = sigmond.MCEnsembleInfo("cls21_s64_t128_D200")
  bins_info = sigmond.MCBinsInfo(ensemble)
  sampling_info = sigmond.MCSamplingInfo()

  xml = ET.Element("MCObservables")
  corr_xml = ET.SubElement(xml, "BLCorrelatorData")
  file_info_xml = ET.SubElement(corr_xml, "FileListInfo")
  ET.SubElement(file_info_xml, "FileNameStub").text = "/home/ahanlon/research/data/D200/omega/mergedcorr"
  ET.SubElement(file_info_xml, "MinFileNumber").text = "0"
  ET.SubElement(file_info_xml, "MaxFileNumber").text = "11"
  obs_xml = ET.tostring(xml)
  obs_xml_h = sigmond.XMLHandler()
  obs_xml_h.set_from_string(obs_xml)
  obs_get_handler = sigmond.MCObsGetHandler(obs_xml_h, bins_info, sampling_info)
  obs_handler = sigmond.MCObsHandler(obs_get_handler, False)


  out_file = h5py.File("omega-d200.hdf5", 'w')
  out_file.attrs.create("cfgList", cfgs)
  out_file.attrs.create("spatExt", 64)
  out_file.attrs.create('tseps', tseps)
  omega_grp = out_file.create_group("omega")


  pSq0_T1um_grp = omega_grp.create_group("pSq0-T1um")
  pSq0_T1um_grp.attrs.create('opList', ['eta(0) SS_0', 'eta(0) SS_1'], shape=(2,), dtype="S11")
  dataset = pSq0_T1um_grp.create_dataset("t0avg", (len(cfgs), len(tseps), 2, 2), dtype='complex128')
  average(dataset, eta_0_ops, eta_1_ops, obs_handler)

  out_file.close()


def average(dataset, eta_ops, phi_ops, obs_handler):
  data = np.empty((len(eta_ops), len(cfgs), len(tseps), 2, 2), dtype='complex128')
  for ind, (eta_op, phi_op) in enumerate(zip(eta_ops, phi_ops)):
    ops = [eta_op, phi_op]

    for (i_snk, snk_op_str), (i_src, src_op_str) in itertools.product(enumerate(ops), enumerate(ops)):
      snk_op = sigmond.OperatorInfo(snk_op_str, sigmond.OpKind.BasicLapH)
      src_op = sigmond.OperatorInfo(src_op_str, sigmond.OpKind.BasicLapH)
      corr = sigmond.CorrelatorAtTimeInfo(snk_op, src_op, 0, False, False, False)
      for it, tsep in enumerate(tseps):
        corr.resetTimeSeparation(tsep)
        corr_re = sigmond.MCObsInfo(corr, sigmond.ComplexArg.RealPart)
        corr_im = sigmond.MCObsInfo(corr, sigmond.ComplexArg.ImaginaryPart)
        bins_re = np.array(obs_handler.getBins(corr_re).array())
        bins_im = np.array(obs_handler.getBins(corr_im).array())

        bins = [r + i*1j for r, i in zip(bins_re, bins_im)]
        data[ind,:,it,i_snk,i_src] = bins

  dataset[:,:,:,:] = np.mean(data, axis=0)


if __name__ == "__main__":
  main()
