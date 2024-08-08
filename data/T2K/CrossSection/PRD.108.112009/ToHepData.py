#!/usr/bin/env python3

import yaml
import csv
import os
import ROOT
import re
from math import sqrt
import requests

from hepdata_lib import Submission, Table, Variable, Uncertainty, RootFileReader

if not os.path.exists("onoffaxis_data_release/analysis_flux.root"):
  if not os.path.exists("onoffaxis_data_release.tar.gz"):
    req = requests.get("https://zenodo.org/records/7768255/files/onoffaxis_data_release.tar.gz?download=1", params={"download": 1})
    if req.status_code != requests.codes.ok:
      raise RuntimeError("Failed to download data release from: https://zenodo.org/records/7768255/files/onoffaxis_data_release.tar.gz?download=1")
    with open("onoffaxis_data_release.tar.gz", 'wb') as fd:
      for chunk in req.iter_content(chunk_size=128):
        fd.write(chunk)
  os.system("mkdir -p onoffaxis_data_release && tar -zxvf onoffaxis_data_release.tar.gz -C onoffaxis_data_release")

ref = "PRD.108.112009"
INSPIRE_id=2646102

def build_flux_table(hname, tname):

  reader = RootFileReader("onoffaxis_data_release/analysis_flux.root")
  
  fh = reader.read_hist_1d(hname)

  #### Build Submission
  EnuVar = Variable("e_nu", is_independent=True, is_binned=True, units="GeV")
  EnuVar.values = fh["x_edges"]

  FluxVar = Variable("flux_nu", is_independent=False, is_binned=False, units="$/cm^{2}/50MeV/10^{21}p.o.t$")
  FluxVar.values = fh["y"]

  FluxVar.add_qualifier("variable_type", "probe_flux")
  FluxVar.add_qualifier("probe_particle", "numu")
  FluxVar.add_qualifier("bin_content_type", "count_density")

  FluxTable = Table(tname)

  FluxTable.add_variable(EnuVar)
  FluxTable.add_variable(FluxVar)

  return FluxTable

def nd280_analysis():
  bins = []

  with open("onoffaxis_data_release/nd280_analysis_binning.csv", newline='') as csvfile:
    csvreader = csv.DictReader(csvfile)
    for row in csvreader:
      bins.append({
        "bin_id": row["bin"], "extents": [
        ( float(row["low_momentum"])/1E3, float(row["high_momentum"])/1E3 ),
        ( float(row["low_angle"]), float(row["high_angle"]) ),
      ]})

  data = {
    "data": {},
    "error": {},
    "nominal_mc": {},
  }

  with open("onoffaxis_data_release/xsec_data_mc.csv", newline='') as csvfile:
    csvreader = csv.DictReader(csvfile)
    for row in csvreader:
      #DictReader doesn't strip keys unfortunately
      for k,v in row.items():
        if k.strip() == "data":
          data["data"][row["bin"]] = float(v)
        elif k.strip() == "nominal_mc":
          data["nominal_mc"][row["bin"]] = float(v)

  with open("onoffaxis_data_release/cov_matrix.csv", newline='') as csvfile:
    csvreader = csv.reader(csvfile)
    for i, row in enumerate(csvreader):
      data["error"][str(i+1)] = { str(j+1): float(x)  for j,x in enumerate(row) }

  #### Build Submission
  CosThetaVar = Variable("cos_theta_mu", is_independent=True, is_binned=True, units="")
  CosThetaVar.values = [ x["extents"][1] for x in bins ]

  PVar = Variable("p_mu", is_independent=True, is_binned=True, units="MeV/c")
  PVar.values = [ x["extents"][0] for x in bins ]

  CrossSection = Variable("cross_section", is_independent=False, is_binned=False, units="$cm${}^{2} c/MeV /Nucleon$")
  CrossSection.values = [ data["data"][x["bin_id"]] for x in bins ]

  CrossSection.add_qualifier("selectfunc", "T2K_CC0Pi_onoffaxis_nu_SelectSignal")
  CrossSection.add_qualifier("cos_theta_mu:projectfunc", "T2K_CC0Pi_onoffaxis_nu_Project_CosThetaMu")
  CrossSection.add_qualifier("p_mu:projectfunc", "T2K_CC0Pi_onoffaxis_nu_Project_PMu")

  CrossSection.add_qualifier("target", "CH")
  CrossSection.add_qualifier("probe_spectra", "flux-offaxis-postfit-fine")
  CrossSection.add_qualifier("variable_type", "cross_section_measurement")

  CrossSectionNEUT = Variable("cross_section_neut-prediction", is_independent=False, is_binned=False, units="$cm${}^{2} c/MeV /Nucleon$")
  CrossSectionNEUT.values = [ data["nominal_mc"][x["bin_id"]] for x in bins ]

  CrossSectionNEUT.add_qualifier("variable_type", "cross_section_prediction")

  TotalUncertainty = Uncertainty("total", is_symmetric=True)
  TotalUncertainty.values = [ sqrt(data["error"][x["bin_id"]][x["bin_id"]]) for x in bins ]

  CrossSection.add_uncertainty(TotalUncertainty)

  xsTable = Table("cross_section-offaxis")
  xsTable.description = """Extracted ND280 cross section as a function of muon momentum in angle bins compared to the nominal NEUT MC prediction. Note that the final bin extending to 30 GeV=c has been omitted for clarity."""
  xsTable.location = "FIG. 21. in the publication"

  xsTable.add_variable(PVar)
  xsTable.add_variable(CosThetaVar)
  xsTable.add_variable(CrossSection)
  xsTable.add_variable(CrossSectionNEUT)
  xsTable.add_image("fig21.png")

  xsTable.keywords["observables"] = ["D2SIG/DP/DCOSTHETA"]
  xsTable.keywords["reactions"] = ["NUMU C --> MU- P"]
  xsTable.keywords["phrases"] = ["Neutrino CC0Pi", "Cross Section"]

  return xsTable

def ingrid_analysis():

  bins = []

  with open("onoffaxis_data_release/ingrid_analysis_binning.csv", newline='') as csvfile:
    csvreader = csv.DictReader(csvfile)
    for row in csvreader:
      bins.append({
        "bin_id": row["bin"], "extents": [
        ( float(row["low_momentum"])/1E3, float(row["high_momentum"])/1E3 ),
        ( float(row["low_angle"]), float(row["high_angle"]) ),
      ]})

  data = {
    "data": {},
    "error": {},
    "nominal_mc": {},
  }

  with open("onoffaxis_data_release/xsec_data_mc.csv", newline='') as csvfile:
    csvreader = csv.DictReader(csvfile)
    for row in csvreader:
      #DictReader doesn't strip keys unfortunately
      for k,v in row.items():
        if k.strip() == "data":
          data["data"][row["bin"]] = float(v)
        elif k.strip() == "nominal_mc":
          data["nominal_mc"][row["bin"]] = float(v)

  with open("onoffaxis_data_release/cov_matrix.csv", newline='') as csvfile:
    csvreader = csv.reader(csvfile)
    for i, row in enumerate(csvreader):
      data["error"][str(i+1)] = { str(j+1): float(x)  for j,x in enumerate(row) }

  #### Build Submission
  CosThetaVar = Variable("cos_theta_mu", is_independent=True, is_binned=True, units="")
  CosThetaVar.values = [ x["extents"][1] for x in bins ]

  PVar = Variable("p_mu", is_independent=True, is_binned=True, units="MeV/c")
  PVar.values = [ x["extents"][0] for x in bins ]

  CrossSection = Variable("cross_section", is_independent=False, is_binned=False, units="$cm${}^{2} c/MeV /Nucleon$")
  CrossSection.values = [ data["data"][x["bin_id"]] for x in bins ]

  CrossSection.add_qualifier("selectfunc", "T2K_CC0Pi_onoffaxis_nu_SelectSignal")
  CrossSection.add_qualifier("cos_theta_mu:projectfunc", "T2K_CC0Pi_onoffaxis_nu_Project_CosThetaMu")
  CrossSection.add_qualifier("p_mu:projectfunc", "T2K_CC0Pi_onoffaxis_nu_Project_PMu")

  CrossSection.add_qualifier("target", "CH")
  CrossSection.add_qualifier("probe_spectra", "flux-onaxis-postfit-fine")
  CrossSection.add_qualifier("variable_type", "cross_section_measurement")
  CrossSection.add_qualifier("pretty_name", r"$p_{\mu}$")

  CrossSectionNEUT = Variable("cross_section_neut-prediction", is_independent=False, is_binned=False, units="$cm${}^{2} c/MeV /Nucleon$")
  CrossSectionNEUT.values = [ data["nominal_mc"][x["bin_id"]] for x in bins ]

  CrossSectionNEUT.add_qualifier("variable_type", "cross_section_prediction")

  TotalUncertainty = Uncertainty("total", is_symmetric=True)
  TotalUncertainty.values = [ sqrt(data["error"][x["bin_id"]][x["bin_id"]]) for x in bins ]

  CrossSection.add_uncertainty(TotalUncertainty)

  xsTable = Table("cross_section-onaxis")
  xsTable.description = """Extracted INGRID cross section as a function of muon momentum in angle bins compared to the nominal NEUT MC prediction. Note that the final bin extending to 30 GeV=c has been omitted for clarity."""
  xsTable.location = """FIG. 22. in the publication."""

  xsTable.add_variable(PVar)
  xsTable.add_variable(CosThetaVar)
  xsTable.add_variable(CrossSection)
  xsTable.add_variable(CrossSectionNEUT)
  xsTable.add_image("fig22.png")

  xsTable.keywords["observables"] = ["D2SIG/DP/DCOSTHETA"]
  xsTable.keywords["reactions"] = ["NUMU C --> MU- P"]
  xsTable.keywords["phrases"] = ["Neutrino CC0Pi", "Cross Section"]

  return xsTable

def joint_analysis():

  #### Build Submission
  bin_i = Variable("bin_i", is_independent=True, is_binned=False, units="")
  bin_i.values = []

  bin_j = Variable("bin_j", is_independent=True, is_binned=False, units="")
  bin_j.values = []

  Covariance = Variable("covariance", is_independent=False, is_binned=False, units=r"$(cm${}^{2} c/MeV /Nucleon)^{2}$")
  Covariance.values = []

  with open("onoffaxis_data_release/cov_matrix.csv", newline='') as csvfile:
    csvreader = csv.reader(csvfile)
    for i, row in enumerate(csvreader):
      for j, val in enumerate(row):
        bin_i.values.append(i)
        bin_j.values.append(j)
        Covariance.values.append(val)

  Invcovariance = Variable("inverse_covariance", is_independent=False, is_binned=False, units=r"$(cm${}^{2} c/MeV /Nucleon)^{-2}$")
  Invcovariance.values = []

  with open("onoffaxis_data_release/inv_matrix.csv", newline='') as csvfile:
    csvreader = csv.reader(csvfile)
    for i, row in enumerate(csvreader):
      for j, val in enumerate(row):
        Invcovariance.values.append(val)

  Covariance.add_qualifier("variable_type", "error_table")
  Covariance.add_qualifier("error_type", "covariance")
  
  Invcovariance.add_qualifier("variable_type", "error_table")
  Invcovariance.add_qualifier("error_type", "inverse_covariance")

  covmatTable = Table("covariance-onoffaxis")
  covmatTable.description = """This table contains the covariance and pre-inverted covariance for the joint on/off axis analysis. See the covered measurements for the constituent measurements."""

  covmatTable.add_variable(bin_i)
  covmatTable.add_variable(bin_j)
  covmatTable.add_variable(Covariance)
  covmatTable.add_variable(Invcovariance)

  jointTable = Table("cross_section-onoffaxis")
  CrossSection = Variable("cross_section", is_independent=False)
  CrossSection.add_qualifier("variable_type", "combined_cross_section_measurement")
  CrossSection.add_qualifier("sub_measurements", "cross_section-offaxis,cross_section-onaxis")
  CrossSection.add_qualifier("error", "covariance-onoffaxis")
  jointTable.add_variable(CrossSection)

  return (covmatTable,jointTable)

submission = Submission()
submission.read_abstract("abstract.txt")

submission.add_table(nd280_analysis())
submission.add_table(ingrid_analysis())
cov,ana = joint_analysis()
submission.add_table(cov)
submission.add_table(ana)
submission.add_table(build_flux_table("ingrid_flux_fine_nominal","flux-onaxis-nominal-fine"))
submission.add_table(build_flux_table("ingrid_flux_coarse_nominal","flux-onaxis-nominal-coarse"))
submission.add_table(build_flux_table("ingrid_flux_fine_postfit","flux-onaxis-postfit-fine"))
submission.add_table(build_flux_table("ingrid_flux_coarse_postfit","flux-onaxis-postfit-coarse"))
submission.add_table(build_flux_table("nd280_flux_fine_nominal","flux-offaxis-nominal-fine"))
submission.add_table(build_flux_table("nd280_flux_coarse_nominal","flux-offaxis-nominal-coarse"))
submission.add_table(build_flux_table("nd280_flux_fine_postfit","flux-offaxis-postfit-fine"))
submission.add_table(build_flux_table("nd280_flux_coarse_postfit","flux-offaxis-postfit-coarse"))

submission.add_additional_resource(description="Python conversion script used to build this submisson. Part of NUISANCE.",
    location="ToHepData.py",
    copy_file=True)

submission.add_additional_resource(description="Selection and projection function examples. Can be executued in the ProSelecta environment v1.0.",
    location="analysis.cxx",
    copy_file=True,
    file_type="ProSelecta")

submission.add_link(description="publication", location="https://doi.org/10.1103/PhysRevD.108.112009")
submission.add_link(description="pre-print", location="https://doi.org/10.48550/arXiv.2303.14228")
submission.add_link(description="official data release", location="https://doi.org/10.5281/zenodo.7768255")
submission.add_link(description="use with NUISANCE3", location="https://github.com/NUISANCEMC/nuisance3")

submission.create_files(f"submission-{INSPIRE_id}", remove_old=True)
