#!/usr/bin/env python3

import yaml
import os
import ROOT
import re
from math import sqrt
import requests
import numpy as np

from hepdata_lib import Submission, Table, Variable, Uncertainty, RootFileReader

ref = "T2K-TN-496"
INSPIRE_id=999999

def GetTH2PolyTable(infilename, polyname, selfunc, xvar, yvar):
  infile = ROOT.TFile.Open(infilename)
  th2poly = infile.Get(polyname)

  XVar = Variable(xvar["name"], is_independent=True, is_binned=True, units=xvar["units"])
  XVar.values = []

  YVar = Variable(yvar["name"], is_independent=True, is_binned=True, units=yvar["units"])
  YVar.values = []

  CrossSection = Variable("cross_section", is_independent=False, is_binned=False, 
                          units=r"$\mathrm{cm}^{2} c/\mathrm{GeV} /\mathrm{degrees}$")
  CrossSection.values = []

  for bin in th2poly.GetBins():
    XVar.values.append([bin.GetXMin(),bin.GetXMax()])
    YVar.values.append([bin.GetYMin(),bin.GetYMax()])
    CrossSection.values.append(bin.GetContent())

  # print(f"{polyname}\n  x: {XVar.values}\n  y: {YVar.values}")

  CrossSection.add_qualifier("variable_type", "cross_section_measurement")
  CrossSection.add_qualifier("selectfunc", f"analysis.cxx:{selfunc}")
  CrossSection.add_qualifier(f"{xvar['name']}:projectfunc", f"analysis.cxx:{xvar['projfunc']}")
  CrossSection.add_qualifier(f"{xvar['name']}:prettyname", xvar['prettyname'])
  CrossSection.add_qualifier(f"{yvar['name']}:projectfunc", f"analysis.cxx:{yvar['projfunc']}")
  CrossSection.add_qualifier(f"{yvar['name']}:prettyname", yvar['prettyname'])
  CrossSection.add_qualifier("cross_section_units", "cm2|PerTargetNucleon|per_bin_width")
  CrossSection.add_qualifier("prettyname", r"$\mathrm{d}^{2}\sigma/\mathrm{d}\delta\alpha_{T}\mathrm{d}\delta{}p_{T}$")
  CrossSection.add_qualifier("probe_flux", "hepdata-sandbox:1722947187v1/t2kflux_2020_minus250ka_nd280_runcond:flux_numu")

  return XVar,YVar,CrossSection

def add_tables(sub):

  xvar = { 'name': 'dat', 
            'prettyname': r'$\delta\alpha_T$', 
            'units': "degrees",
            'projfunc' : "T2K_CC0Pi_CO_TKI_nu_2024_Project_DeltaAlphaT_deg" }
  yvar = { 'name': 'dpt', 
            'prettyname': r'$\delta{}p_T$', 
            'units': r"$\mathrm{GeV}/c$",
            'projfunc' : "T2K_CC0Pi_CO_TKI_nu_2024_Project_DeltaPT_GeV_c" }

  C_XVar,C_YVar,C_CrossSection = GetTH2PolyTable("check_calcxsec_1.root", 
    "data_result_C_2d", "T2K_CC0Pi_CO_TKI_nu_2024_Select_C", 
    xvar, yvar)
  O_XVar,O_YVar,O_CrossSection = GetTH2PolyTable("check_calcxsec_1.root", 
    "data_result_O_2d", "T2K_CC0Pi_CO_TKI_nu_2024_Select_O", 
    xvar, yvar)

  C_CrossSection.add_qualifier("target", "C")
  C_CrossSection.add_qualifier("errors", "covariance_c:covariance")
  O_CrossSection.add_qualifier("target", "O")
  O_CrossSection.add_qualifier("errors", "covariance_o:covariance")

  reader = RootFileReader("check_calcxsec_1.root")
  CO_cov_matrix = reader.read_hist_2d("cov_matrix")

  C_nbins = len(C_XVar.values)
  O_nbins = len(O_XVar.values)
  nbins = C_nbins + O_nbins
  covmat = np.array(CO_cov_matrix["z"]).reshape(nbins,nbins)

  C_TotalUncertainty = Uncertainty("total", is_symmetric=True)
  C_TotalUncertainty.values = np.sqrt(np.diagonal(covmat[:C_nbins,:C_nbins]))
  C_CrossSection.add_uncertainty(C_TotalUncertainty)

  O_TotalUncertainty = Uncertainty("total", is_symmetric=True)
  O_TotalUncertainty.values = np.sqrt(np.diagonal(covmat[C_nbins:,C_nbins:]))
  O_CrossSection.add_uncertainty(O_TotalUncertainty)

  C_xsTable = Table("cross_section_c")
  O_xsTable = Table("cross_section_o")

  C_xsTable.add_variable(C_XVar)
  C_xsTable.add_variable(C_YVar)
  C_xsTable.add_variable(C_CrossSection)

  C_xsTable.keywords["observables"] = ["D2SIG/DDPT/DALPHAT"]
  C_xsTable.keywords["reactions"] = ["NUMU C --> MU- P"]
  C_xsTable.keywords["phrases"] = ["Neutrino CC0Pi", "Cross Section"]

  O_xsTable.add_variable(O_XVar)
  O_xsTable.add_variable(O_YVar)
  O_xsTable.add_variable(O_CrossSection)

  O_xsTable.keywords["observables"] = ["D2SIG/DDPT/DALPHAT"]
  O_xsTable.keywords["reactions"] = ["NUMU C --> MU- P"]
  O_xsTable.keywords["phrases"] = ["Neutrino CC0Pi", "Cross Section"]

  C_bin_i = Variable("bin_i", is_independent=True, is_binned=False, units="")
  C_bin_i.values = []

  C_bin_j = Variable("bin_j", is_independent=True, is_binned=False, units="")
  C_bin_j.values = []

  #read_hist_2d loops over x then y (column-major) vs, TH2D which is row-major
  for i in range(C_nbins):
    for j in range(C_nbins):
      C_bin_i.values.append(i)
      C_bin_j.values.append(j)

  C_Covariance = Variable("covariance", is_independent=False, is_binned=False, units="")
  C_Covariance.values = covmat[:C_nbins,:C_nbins].flatten()

  C_Covariance.add_qualifier("variable_type", "error_table")
  C_Covariance.add_qualifier("error_type", "covariance")

  C_covTable = Table("covariance_c")
  C_covTable.add_variable(C_bin_i)
  C_covTable.add_variable(C_bin_j)
  C_covTable.add_variable(C_Covariance)


  O_bin_i = Variable("bin_i", is_independent=True, is_binned=False, units="")
  O_bin_i.values = []

  O_bin_j = Variable("bin_j", is_independent=True, is_binned=False, units="")
  O_bin_j.values = []

  #read_hist_2d loops over x then y (column-major) vs, TH2D which is row-major
  for i in range(O_nbins):
    for j in range(O_nbins):
      O_bin_i.values.append(i)
      O_bin_j.values.append(j)

  O_Covariance = Variable("covariance", is_independent=False, is_binned=False, units="")
  O_Covariance.values = covmat[C_nbins:,C_nbins:].flatten()

  O_Covariance.add_qualifier("variable_type", "error_table")
  O_Covariance.add_qualifier("error_type", "covariance")

  O_covTable = Table("covariance_o")
  O_covTable.add_variable(O_bin_i)
  O_covTable.add_variable(O_bin_j)
  O_covTable.add_variable(O_Covariance)

### Joint Measurement

  CO_bin_i = Variable("bin_i", is_independent=True, is_binned=False, units="")
  CO_bin_i.values = []

  CO_bin_j = Variable("bin_j", is_independent=True, is_binned=False, units="")
  CO_bin_j.values = []

  #read_hist_2d loops over x then y (column-major) vs, TH2D which is row-major
  for i in range(nbins):
    for j in range(nbins):
      CO_bin_i.values.append(i)
      CO_bin_j.values.append(j)

  CO_Covariance = Variable("covariance", is_independent=False, is_binned=False, units="")
  CO_Covariance.values = CO_cov_matrix["z"]

  CO_Covariance.add_qualifier("variable_type", "error_table")
  CO_Covariance.add_qualifier("error_type", "covariance")

  CO_covTable = Table("covariance_co")
  CO_covTable.add_variable(CO_bin_i)
  CO_covTable.add_variable(CO_bin_j)
  CO_covTable.add_variable(CO_Covariance)

  CO_xsTable = Table("cross_section_co")
  CO_CrossSection = Variable("cross_section_co", is_independent=False)
  CO_CrossSection.add_qualifier("variable_type", "composite_cross_section_measurement")
  CO_CrossSection.add_qualifier("sub_measurements", "cross_section_c,cross_section_o")
  CO_CrossSection.add_qualifier("errors", "covariance_co:covariance")
  CO_xsTable.add_variable(CO_CrossSection)

  sub.add_table(C_xsTable)
  sub.add_table(C_covTable)
  sub.add_table(O_xsTable)
  sub.add_table(O_covTable)
  sub.add_table(CO_xsTable)
  sub.add_table(CO_covTable)

submission = Submission()

add_tables(submission)

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
