#!/usr/bin/env python3

import yaml
import os
import ROOT
import re
from math import sqrt
import requests
import numpy as np

from hepdata_lib import Submission, Table, Variable, Uncertainty, RootFileReader

ref = ""
INSPIRE_id=999999

def flux_CV():
  reader = RootFileReader("release/best_fit_flux.root")
  fluxes =  {
    "flux_numu": reader.read_hist_1d("flux_numu"),
    "flux_numubar": reader.read_hist_1d("flux_numubar"),
    "flux_nue": reader.read_hist_1d("flux_nue"),
    "flux_nuebar": reader.read_hist_1d("flux_nuebar"),
  }

  EnuVar = Variable("e_nu", is_independent=True, is_binned=True, units="GeV")
  EnuVar.values = fluxes["flux_nue"]["x_edges"]

  FluxTable = Table("best_fit_nue_flux")

  FluxTable.add_variable(EnuVar)

  FluxVar = Variable(f"flux_nue", is_independent=False, is_binned=False, units="/cm$^{2}$ /10$^{21}$ p.o.t /50 MeV")
  FluxVar.values = fluxes["flux_nue"]["y"]

  FluxVar.add_qualifier("variable_type", "probe_flux")
  FluxVar.add_qualifier("probe_particle", "nue")
  FluxVar.add_qualifier("bin_content_type", "count_density")

  FluxTable.add_variable(FluxVar)

  return FluxTable

def analysis_pe_costheta_ppi():

  reader = RootFileReader("release/xsec_data_3D_fix.root")
  
  flat_xs =  reader.read_hist_1d("sel_best_fit")

# p_e: [0, 0.350, 1.700, 30.000] GeV
# cos(θ_e): [-1.0, 0.7, 0.94, 1.0]
# p_π: [0, 0.450, 1.500, 30.000] GeV

  pe_bin_edges = [ [0.350,1.700], [1.700,30.000] ]
  costhetae_bin_edges = [ [0.7,0.94], [0.94,1] ]
  ppi_bin_edges = [ [0,0.450], [0.450,1.500] ]

  binedges_3d = { "pe": [],
                  "costhetae": [],
                  "ppi": [] }

  for bz in ppi_bin_edges:
    for by in costhetae_bin_edges:
      for bx in pe_bin_edges:
        binedges_3d["pe"].append(bx)
        binedges_3d["costhetae"].append(by)
        binedges_3d["ppi"].append(bz)

  ElecMomVar = Variable("p_e", is_independent=True, is_binned=True, units="GeV/c")
  ElecMomVar.values = binedges_3d["pe"]

  CosThetaElecVar = Variable("cos_theta_e", is_independent=True, is_binned=True, units="")
  CosThetaElecVar.values = binedges_3d["costhetae"]

  PionMomVar = Variable("p_pi", is_independent=True, is_binned=True, units="GeV/c")
  PionMomVar.values = binedges_3d["ppi"]

  CrossSection = Variable("cross_section", is_independent=False, is_binned=False, 
                          units=r"$\mathrm{cm}^{2}\ c/\mathrm{GeV}\ /\mathrm{Nucleon}$")
  CrossSection.values = flat_xs["y"][3:] #first three bins are OOPS

  CrossSection.add_qualifier("selectfunc", "analysis.cxx:T2K_nueCCPiplus_2024_Select")
  CrossSection.add_qualifier("p_e:projectfunc", "analysis.cxx:T2K_nueCCPiplus_2024_ElectronMomentum")
  CrossSection.add_qualifier("cos_theta_e:projectfunc", "analysis.cxx:T2K_nueCCPiplus_2024_CosThetaElectron")
  CrossSection.add_qualifier("p_pi:projectfunc", "analysis.cxx:T2K_nueCCPiplus_2024_PiPlusMomentum")
  CrossSection.add_qualifier("p_e:prettyname", r"$p_e$")
  CrossSection.add_qualifier("cos_theta_e:prettyname", r"$\cos{\theta_e}$")
  CrossSection.add_qualifier("p_pi:prettyname", r"$p_\pi$")
  CrossSection.add_qualifier("prettyname", r"$\mathrm{d}^{3}\sigma/\mathrm{d}p_e\mathrm{d}cos{\theta_e}\mathrm{d}p_\pi$")
  CrossSection.add_qualifier("cross_section_units", "cm2|PerTargetNucleon|per_bin_width")

  CrossSection.add_qualifier("target", "CH")
  CrossSection.add_qualifier("probe_flux", "best_fit_nue_flux")
  CrossSection.add_qualifier("variable_type", "cross_section_measurement")
  CrossSection.add_qualifier("error", "covariance:covariance")

  rf = ROOT.TFile.Open("xsec_result.root")
  rcov_matrix = rf.Get("TMatrixTSym<double>")

  cov_matrix = np.array([[rcov_matrix[i+3,j+3] for i in range(8)] for j in range(8)])

  TotalUncertainty = Uncertainty("total", is_symmetric=True)
  TotalUncertainty.values = np.sqrt(np.diagonal(cov_matrix))

  CrossSection.add_uncertainty(TotalUncertainty)

  xsTable = Table("cross_section")

  xsTable.add_variable(ElecMomVar)
  xsTable.add_variable(CosThetaElecVar)
  xsTable.add_variable(PionMomVar)
  xsTable.add_variable(CrossSection)

  xsTable.keywords["observables"] = ["D3SIG/DPE/DCOSTHETA/DPPi"]
  xsTable.keywords["reactions"] = ["NUE CH --> E- Pi+"]
  xsTable.keywords["phrases"] = ["Neutrino CCPi+", "Cross Section"]

  bin_i = Variable("bin_i", is_independent=True, is_binned=False, units="")
  bin_i.values = np.array([[i for i in range(8)] for j in range(8)]).ravel()

  bin_j = Variable("bin_j", is_independent=True, is_binned=False, units="")
  bin_j.values = np.array([[j for i in range(8)] for j in range(8)]).ravel()

  Covariance = Variable("covariance", is_independent=False, is_binned=False, units=r"$(\mathrm{cm}^{2}\ (c/\mathrm{GeV})^{2} \/\mathrm{Nucleon})^{2}$")
  Covariance.values = np.ravel(cov_matrix)

  Covariance.add_qualifier("variable_type", "error_table")
  Covariance.add_qualifier("error_type", "covariance")

  cov = Table("covariance")

  cov.add_variable(bin_i)
  cov.add_variable(bin_j)
  cov.add_variable(Covariance)

  return xsTable,cov

submission = Submission()

xsTable, cov = analysis_pe_costheta_ppi()
submission.add_table(xsTable)
submission.add_table(cov)

submission.add_table(flux_CV())

submission.add_additional_resource(description="Python conversion script used to build this submisson. Part of NUISANCE.",
    location="ToHepData.py",
    copy_file=True)

submission.add_additional_resource(description="Selection and projection function examples. Can be executued in the ProSelecta environment v1.0.",
    location="analysis.cxx",
    copy_file=True,
    file_type="ProSelecta")

submission.add_link(description="Use with NUISANCE3", location="https://github.com/NUISANCEMC/nuisance3")
submission.add_link(description="Adheres to the NUISANCE HEPData Conventions", location="https://github.com/NUISANCEMC/HEPData/tree/main")

submission.create_files(f"submission-{INSPIRE_id}", remove_old=True)
