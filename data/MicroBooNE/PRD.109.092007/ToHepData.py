#!/usr/bin/env python3

import yaml
import csv
import os
import ROOT
import re
import numpy as np
from math import sqrt

from hepdata_lib import Submission, Table, Variable, Uncertainty, RootFileReader

ref = "PRD.109.092007"
INSPIRE_id=2709091


def build_flux_table(hname, tname):
  inFileName = "MicroBooNE_FHC_numu_flux.root"

  reader = RootFileReader(inFileName)
  
  fh = reader.read_hist_1d(hname)

  #### Build Submission
  EnuVar = Variable("e_nu", is_independent=True, is_binned=True, units="GeV")
  EnuVar.values = fh["x_edges"]

  FluxVar = Variable("flux_nu", is_independent=False, is_binned=False, units="$/cm^{2}$")
  FluxVar.values = fh["y"]

  FluxVar.add_qualifier("variable_type", "probe_flux")
  FluxVar.add_qualifier("probe_particle", "numu")
  FluxVar.add_qualifier("bin_content_type", "count")

  FluxTable = Table(tname)

  FluxTable.add_variable(EnuVar)
  FluxTable.add_variable(FluxVar)

  return FluxTable

SerialDeltaPn_DeltaAlpha3Dq_bins = np.array([
  [
    [0, 0.1],
    [0.1, 0.2],
    [0.2, 0.3],
    [0.3, 0.4],
    [0.4, 0.85],
    [0, 0.1],
    [0.1, 0.2],
    [0.2, 0.3],
    [0.3, 0.4],
    [0.4, 0.55],
    [0.55, 0.85],
    [0, 0.08],
    [0.08, 0.15],
    [0.15, 0.23],
    [0.23, 0.3],
    [0.3, 0.38],
    [0.38, 0.45],
    [0.45, 0.85],
    [0, 0.08],
    [0.08, 0.15],
    [0.15, 0.23],
    [0.23, 0.3],
    [0.3, 0.4],
    [0.4, 0.85],
  ],
  [
    [0, 45],
    [0, 45],
    [0, 45],
    [0, 45],
    [0, 45],
    [45, 90],
    [45, 90],
    [45, 90],
    [45, 90],
    [45, 90],
    [45, 90],
    [90,135],
    [90,135],
    [90,135],
    [90,135],
    [90,135],
    [90,135],
    [90,135],
    [135,180],
    [135,180],
    [135,180],
    [135,180],
    [135,180],
    [135,180],
  ]
])

SerialDeltaAlpha3Dq_DeltaPn_bins = np.array([
  [
    [0, 35],
    [35, 70],
    [70, 92],
    [92, 114],
    [114, 136],
    [136, 158],
    [158, 180],
    [0, 25],
    [25, 50],
    [50, 70],
    [70, 90],
    [90, 110],
    [110, 130],
    [130, 150],
    [150, 180],
    [0, 35],
    [35, 70],
    [70, 90],
    [90, 110],
    [110, 130],
    [130, 145],
    [145, 180],
  ],
  [
    [0, 0.2],
    [0, 0.2],
    [0, 0.2],
    [0, 0.2],
    [0, 0.2],
    [0, 0.2],
    [0, 0.2],
    [0.2, 0.4],
    [0.2, 0.4],
    [0.2, 0.4],
    [0.2, 0.4],
    [0.2, 0.4],
    [0.2, 0.4],
    [0.2, 0.4],
    [0.2, 0.4],
    [0.4,100],
    [0.4,100],
    [0.4,100],
    [0.4,100],
    [0.4,100],
    [0.4,100],
    [0.4,100],
  ]
])

def double_diff_measurement(varname, bins, xsunits, covunits, *args, **kwargs):
  binFileName = "release/DataRelease.root"
  reader = RootFileReader(binFileName)

  inHist = reader.read_hist_1d(f"TotalUnc_{varname}")
  inCov = reader.read_hist_2d(f"Cov_{varname}")
  inAc = reader.read_hist_2d(f"Ac_{varname}")

  if varname == "SerialDeltaPn_DeltaAlpha3Dq":
    xvar = "pn"
    yvar = "alpha3d"
  elif varname == "SerialDeltaAlpha3Dq_DeltaPn":
    yvar = "pn"
    xvar = "alpha3d"

  varname = varname.lower()

  XVar = Variable(xvar, is_independent=True, is_binned=True, units="")
  XVar.values = bins[0]

  YVar = Variable(yvar, is_independent=True, is_binned=True, units="")
  YVar.values = bins[1]

  CrossSection = Variable("cross_section", is_independent=False, is_binned=False, units=xsunits)
  CrossSection.values = inHist["y"]

  CrossSection.add_qualifier("variable_type", "cross_section_measurement")
  CrossSection.add_qualifier("selectfunc", "analysis.cxx:MicroBooNE_CC0Pi_GKI_nu_SelectSignal")
  CrossSection.add_qualifier("pn:projectfunc", "analysis.cxx:MicroBooNE_CC0Pi_GKI_nu_pn")
  CrossSection.add_qualifier("alpha3d:projectfunc", "analysis.cxx:MicroBooNE_CC0Pi_GKI_nu_alpha3d")
  CrossSection.add_qualifier("target", "Ar")
  CrossSection.add_qualifier("cross_section_units", "1e-38 cm2|PerTarget|per_bin_width")
  CrossSection.add_qualifier("probe_flux", "microboone_flux_numu")
  CrossSection.add_qualifier("errors", f"covariance-{varname}")
  CrossSection.add_qualifier("smearing", f"smearing-{varname}")

  TotalUncertainty = Uncertainty("total", is_symmetric=True)
  TotalUncertainty.values = inHist["dy"]

  #data release comes without PDFing by bin area for 2D histograms only.
  bin_areas = []
  for i in range(len(XVar.values)):
    bin_areas.append((XVar.values[i][1] - XVar.values[i][0])*(YVar.values[i][1] - YVar.values[i][0]))
  CrossSection.values = np.divide(CrossSection.values,bin_areas)
  TotalUncertainty.values = np.divide(TotalUncertainty.values,bin_areas)

  CrossSection.add_uncertainty(TotalUncertainty)

  xsTable = Table(f"cross_section-{varname}")
  xsTable.description = ""
  xsTable.location = ""

  xsTable.add_variable(XVar)
  xsTable.add_variable(YVar)
  xsTable.add_variable(CrossSection)

  if "observables" in kwargs:
    xsTable.keywords["observables"] = kwargs["observables"]
  if "reactions" in kwargs:
    xsTable.keywords["reactions"] = kwargs["reactions"]
  if "phrases" in kwargs:
    xsTable.keywords["phrases"] = kwargs["phrases"]

##### matrices

  bin_i = Variable("bin_i", is_independent=True, is_binned=False, units="")
  bin_i.values = []

  bin_j = Variable("bin_j", is_independent=True, is_binned=False, units="")
  bin_j.values = []

  #read_hist_2d loops over x then y (column-major) vs, TH2D which is row-major
  for i in range(len(inHist["x"])):
    for j in range(len(inHist["x"])):
      bin_i.values.append(i)
      bin_j.values.append(j)

  Covariance = Variable("covariance", is_independent=False, is_binned=False, units=covunits)
  Covariance.values = inCov["z"]

  Covariance.add_qualifier("variable_type", "error_table")
  Covariance.add_qualifier("error_type", "covariance")

  covTable = Table(f"covariance-{varname}")
  covTable.description = ""
  covTable.location = ""

  covTable.add_variable(bin_i)
  covTable.add_variable(bin_j)
  covTable.add_variable(Covariance)

  bin_true = Variable("bin_true", is_independent=True, is_binned=False, units="")
  bin_true.values = []

  bin_smear = Variable("bin_smear", is_independent=True, is_binned=False, units="")
  bin_smear.values = []

  #read_hist_2d loops over x then y (column-major) vs, TH2D which is row-major
  for i in range(len(inHist["x"])):
    for j in range(len(inHist["x"])):
      bin_true.values.append(i)
      bin_smear.values.append(j)

  WSmear = Variable("wiener_svd-smearing-matrix", is_independent=False, is_binned=False)
  WSmear.values = inAc["z"]

  WSmear.add_qualifier("variable_type", "smearing_table")
  WSmear.add_qualifier("smearing_type", "smeared_by_true_matrix")

  smearTable = Table(f"smearing-{varname}")
  smearTable.description = ""
  smearTable.location = ""

  smearTable.add_variable(bin_true)
  smearTable.add_variable(bin_smear)
  smearTable.add_variable(WSmear)

  return (xsTable,covTable,smearTable)

def single_diff_measurement(varname, projname, xsunits, covunits, *args, **kwargs):
  binFileName = "release/DataRelease.root"
  reader = RootFileReader(binFileName)

  inHist = reader.read_hist_1d(f"TotalUnc_{varname}")
  inCov = reader.read_hist_2d(f"Cov_{varname}")
  inAc = reader.read_hist_2d(f"Ac_{varname}")

  varname = varname.lower()

#### Build Submission
  XVar = Variable(projname, is_independent=True, is_binned=True, units="")
  XVar.values = inHist["x_edges"]

  if "x_pretty_name" in kwargs:
    XVar.add_qualifier("pretty_name", kwargs["x_pretty_name"])

  CrossSection = Variable("cross_section", is_independent=False, is_binned=False, units=xsunits)
  CrossSection.values = inHist["y"]

  CrossSection.add_qualifier("variable_type", "cross_section_measurement")
  CrossSection.add_qualifier("selectfunc", "analysis.cxx:MicroBooNE_CC0Pi_GKI_nu_SelectSignal")
  CrossSection.add_qualifier(f"{projname}:projectfunc", f"analysis.cxx:MicroBooNE_CC0Pi_GKI_nu_{projname}")
  CrossSection.add_qualifier("target", "Ar")
  CrossSection.add_qualifier("cross_section_units", "1e-38 cm2|PerTarget|per_bin_width")
  CrossSection.add_qualifier("probe_flux", "microboone_flux_numu")
  CrossSection.add_qualifier("errors", f"covariance-{varname}")
  CrossSection.add_qualifier("smearing", f"smearing-{varname}")
  if "pretty_name" in kwargs:
    CrossSection.add_qualifier("pretty_name", kwargs["pretty_name"])

  TotalUncertainty = Uncertainty("total", is_symmetric=True)
  TotalUncertainty.values = inHist["dy"]

  CrossSection.add_uncertainty(TotalUncertainty)

  xsTable = Table(f"cross_section-{varname}")
  xsTable.description = ""
  xsTable.location = ""

  xsTable.add_variable(XVar)
  xsTable.add_variable(CrossSection)

  if "observables" in kwargs:
    xsTable.keywords["observables"] = kwargs["observables"]
  if "reactions" in kwargs:
    xsTable.keywords["reactions"] = kwargs["reactions"]
  if "phrases" in kwargs:
    xsTable.keywords["phrases"] = kwargs["phrases"]

##### matrices

  bin_i = Variable("bin_i", is_independent=True, is_binned=False, units="")
  bin_i.values = []

  bin_j = Variable("bin_j", is_independent=True, is_binned=False, units="")
  bin_j.values = []

  #read_hist_2d loops over x then y (column-major) vs, TH2D which is row-major
  for i in range(len(inHist["x"])):
    for j in range(len(inHist["x"])):
      bin_i.values.append(i)
      bin_j.values.append(j)

  Covariance = Variable("covariance", is_independent=False, is_binned=False, units=covunits)
  Covariance.values = inCov["z"]

  Covariance.add_qualifier("variable_type", "error_table")
  Covariance.add_qualifier("error_type", "covariance")

  covTable = Table(f"covariance-{varname}")
  covTable.description = ""
  covTable.location = ""

  covTable.add_variable(bin_i)
  covTable.add_variable(bin_j)
  covTable.add_variable(Covariance)

  bin_true = Variable("bin_true", is_independent=True, is_binned=False, units="")
  bin_true.values = []

  bin_smear = Variable("bin_smear", is_independent=True, is_binned=False, units="")
  bin_smear.values = []

  #read_hist_2d loops over x then y (column-major) vs, TH2D which is row-major
  for i in range(len(inHist["x"])):
    for j in range(len(inHist["x"])):
      bin_true.values.append(i)
      bin_smear.values.append(j)

  WSmear = Variable("wiener_svd-smearing-matrix", is_independent=False, is_binned=False)
  WSmear.values = inAc["z"]

  WSmear.add_qualifier("variable_type", "smearing_table")
  WSmear.add_qualifier("smearing_type", "smeared_by_true_matrix")

  smearTable = Table(f"smearing-{varname}")
  smearTable.description = ""
  smearTable.location = ""

  smearTable.add_variable(bin_true)
  smearTable.add_variable(bin_smear)
  smearTable.add_variable(WSmear)

  return (xsTable,covTable,smearTable)

submission = Submission()
submission.read_abstract("abstract.txt")

pd_tables = single_diff_measurement(varname="DeltaPn", projname="pn", xsunits=r"$\mathrm{cm}^{2}\ c/\mathrm{GeV}\ /\mathrm{Ar}$", covunits=r"$(\mathrm{cm}^{2}\ c/\mathrm{GeV}\ /\mathrm{Ar})^{2}$")
submission.add_table(pd_tables[0])
submission.add_table(pd_tables[1])
submission.add_table(pd_tables[2])

pd_tables = single_diff_measurement(varname="DeltaAlpha3Dq", projname="alpha3d", xsunits=r"$\mathrm{cm}^{2}\ /\mathrm{degrees}\ /\mathrm{Ar}$", covunits=r"$(\mathrm{cm}^{2}\ c/\mathrm{degrees}\ /\mathrm{Ar})^{2}$")
submission.add_table(pd_tables[0])
submission.add_table(pd_tables[1])
submission.add_table(pd_tables[2])

pd_tables = single_diff_measurement(varname="DeltaPhi3D", projname="phi3d", xsunits=r"$\mathrm{cm}^{2}\ /\mathrm{degrees}\ /\mathrm{Ar}$", covunits=r"$(\mathrm{cm}^{2}\ /\mathrm{degrees}\ /\mathrm{Ar})^{2}$")
submission.add_table(pd_tables[0])
submission.add_table(pd_tables[1])
submission.add_table(pd_tables[2])

pd_tables = single_diff_measurement(varname="DeltaPnPar", projname="pn_para", xsunits=r"$\mathrm{cm}^{2}\ c/\mathrm{GeV}\ /\mathrm{Ar}$", covunits=r"$(\mathrm{cm}^{2}\ /\mathrm{GeV}\ /\mathrm{Ar})^{2}$")
submission.add_table(pd_tables[0])
submission.add_table(pd_tables[1])
submission.add_table(pd_tables[2])

pd_tables = single_diff_measurement(varname="DeltaPnPerp", projname="pn_perp", xsunits=r"$\mathrm{cm}^{2}\ c/\mathrm{GeV}\ /\mathrm{Ar}$", covunits=r"$(\mathrm{cm}^{2}\ c/\mathrm{GeV}\ /\mathrm{Ar})^{2}$")
submission.add_table(pd_tables[0])
submission.add_table(pd_tables[1])
submission.add_table(pd_tables[2])

pd_tables = single_diff_measurement(varname="DeltaPnPerpx", projname="pn_perp_x", xsunits=r"$\mathrm{cm}^{2}\ c/\mathrm{GeV}\ /\mathrm{Ar}$", covunits=r"$(\mathrm{cm}^{2}\ c/\mathrm{GeV}\ /\mathrm{Ar})^{2}$")
submission.add_table(pd_tables[0])
submission.add_table(pd_tables[1])
submission.add_table(pd_tables[2])

pd_tables = single_diff_measurement(varname="DeltaPnPerpy", projname="pn_perp_y", xsunits=r"$\mathrm{cm}^{2}\ c/\mathrm{GeV}\ /\mathrm{Ar}$", covunits=r"$(\mathrm{cm}^{2}\ c/\mathrm{GeV}\ /\mathrm{Ar})^{2}$")
submission.add_table(pd_tables[0])
submission.add_table(pd_tables[1])
submission.add_table(pd_tables[2])

pd_tables = double_diff_measurement(varname="SerialDeltaPn_DeltaAlpha3Dq", bins=SerialDeltaPn_DeltaAlpha3Dq_bins, xsunits=r"$\mathrm{cm}^{2}\ c/\mathrm{GeV}\ /\mathrm{degrees} \ /\mathrm{Ar}$", covunits=r"$(\mathrm{cm}^{2}\ c/\mathrm{GeV}\ /\mathrm{degrees} \ /\mathrm{Ar})^{2}$")
submission.add_table(pd_tables[0])
submission.add_table(pd_tables[1])
submission.add_table(pd_tables[2])

pd_tables = double_diff_measurement(varname="SerialDeltaAlpha3Dq_DeltaPn", bins=SerialDeltaAlpha3Dq_DeltaPn_bins, xsunits=r"$\mathrm{cm}^{2}\ c/\mathrm{GeV}\ /\mathrm{degrees} \ /\mathrm{Ar}$", covunits=r"$(\mathrm{cm}^{2}\ c/\mathrm{GeV}\ /\mathrm{degrees} \ /\mathrm{Ar})^{2}$")
submission.add_table(pd_tables[0])
submission.add_table(pd_tables[1])
submission.add_table(pd_tables[2])

submission.add_table(build_flux_table("numu_hist", "microboone_flux_numu"))

submission.add_additional_resource(description="Python conversion script used to build this submisson. Part of NUISANCE.",
    location="ToHepData.py",
    copy_file=True)

submission.add_additional_resource(description="ROOT version of the flux provided by the MicroBooNE collaboration.",
    location="MicroBooNE_FHC_numu_flux.root",
    copy_file=True)

submission.add_additional_resource(description="Selection and projection function examples. Can be executued in the ProSelecta environment v1.0.",
    location="analysis.cxx",
    copy_file=True, file_type="ProSelecta")

submission.add_additional_resource(description="2D Binning scheme",
    location="release/BinScheme.txt",
    copy_file=True)
submission.add_additional_resource(description="Official data release documentation",
    location="release/README.txt",
    copy_file=True)

submission.add_link(description="pre-print", location="https://doi.org/10.48550/arXiv.2310.06082")
submission.add_link(description="publication", location="https://doi.org/10.1103/PhysRevD.108.112009")
submission.add_link(description="Use with NUISANCE3", location="https://github.com/NUISANCEMC/nuisance3")
submission.add_link(description="Adheres to the NUISANCE HEPData Conventions", location="https://github.com/NUISANCEMC/HEPData/tree/main")

submission.create_files(f"submission-{INSPIRE_id}", remove_old=True)
