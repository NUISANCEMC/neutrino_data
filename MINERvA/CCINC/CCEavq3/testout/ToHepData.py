#!/usr/bin/env python3

import yaml
import csv
import os
import ROOT
import re
import numpy as np

from hepdata_lib import Submission, Table, Variable, Uncertainty

probe = "nu"
expt = "MINERvA"
target = "CH"
species = "numu"
ref = "PRD.26.537"

INSPIRE_id=182176

binx = np.array([0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8])
biny = np.array([0.0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.20, 0.25, 0.30, 0.35, 0.40, 0.50, 0.60, 0.80])


  # CreateDataHistogram(7, binx, 17, biny);

  # SetDataValuesFromTextFile( fSettings.GetDataInput() );
  # ScaleData(1E-42);

  # SetMapValuesFromText( fSettings.GetMapInput() );

  # SetCholDecompFromTextFile( fSettings.GetCovarInput(), 67);
  # ScaleCovar(1E-16);


h = ROOT.TH2D("h","h",6, binx, 16, biny)
bins = []
for i in range(h.GetNbinsX()):
  for j in range(h.GetNbinsY()):
    bins.append( [h.GetXaxis().GetBinLowEdge(i+1),  h.GetXaxis().GetBinLowEdge(i+2),  
                  h.GetYaxis().GetBinLowEdge(i+1),  h.GetYaxis().GetBinLowEdge(i+2)] )
    

def add_default_qualifiers(var):
  var.add_qualifier("likelihood", "CorrelatedChi2")
  var.add_qualifier("finalize", "FATXNormalized")
  var.add_qualifier("weight", "DefaultWeight")
  var.add_qualifier("targets", "CH")
  var.add_qualifier("species", "numu")
  var.add_qualifier("flux", "inspire:182176")

def data_2D_load(path):
  f = open(path)
  vals = []
  for l in f:
    vals.append(l.strip().split())
  return vals

def data_table():

  vals = data_2D_load("data_2D.txt")
  maps = data_2D_load("map_2D.txt")
  cov = data_2D_load("covar_2D.txt")
  # bins = data_2D_load("binning.txt")

  # print(cov)
  dataset = {}
  dataerr = {}
  maxid = -1
  for i in range(len(vals)):
    for j in range(len(vals[i])):
      id = int(maps[i][j])
      if id > maxid: maxid = id
      dataset[id]  = float(vals[i][j])
      dataerr[id]  = np.sqrt(1.0) #float(cov[id][id]))

  data = {"values": [], "errors": [], "xbins": [], "ybins": []}
  for i in range(maxid):
      if i == 0: continue
      data["values"].append( dataset[i]*1E-42 )
      data["errors"].append( dataerr[i]*1E-42 )
      data["xbins"].append( [bins[i][0], bins[i][1]] )
      data["ybins"].append( [bins[i][2], bins[i][3]] )

  #### Build Submission
  q3Var = Variable("q3", is_independent=True, is_binned=True, units="[GeV/c]")
  q3Var.values = data["xbins"]

  EavVar = Variable("Eav", is_independent=True, is_binned=True, units="[GeV/c]")
  EavVar.values = data["ybins"]

  dSigma = Variable("dSigma", is_independent=False, is_binned=False, units="d^{2}#sigma/dEav dq_{3} [cm^{2}/nucleon/GeV^{2}]")
  dSigma.values = data["values"]

  unc_data = Uncertainty("Total errors", is_symmetric=True)
  unc_data.values = data["errors"]

  # dSigma.add_uncertainty(unc_data)

  add_default_qualifiers(dSigma)
  dSigma.add_qualifier("Filter", "MINERvA_CCINC_CCEavq3_Filter")
  dSigma.add_qualifier("q3", "MINERvA_CCINC_CCEavq3_Project_q3")
  dSigma.add_qualifier("Eav", "MINERvA_CCINC_CCEavq3_Project_Eav")


  xsTable = Table("dSigma-q3-Eav")
  xsTable.description = """
  """

  xsTable.add_variable(q3Var)
  xsTable.add_variable(EavVar)
  xsTable.add_variable(dSigma)

  xsTable.keywords["observables"] = ["DSIGMA"]
  xsTable.keywords["reactions"] = ["NUMU CH --> MU- CH"]
  xsTable.keywords["phrases"] = ["Neutrino", "CCINC", "Differential"]

  return xsTable

submission = Submission()
submission.read_abstract("abstract.txt")
submission.add_table(data_table())

# submission.add_additional_resource(description="Original digitised file used in NUISANCE.",
    # location="nuisance_digitized_file.root",
    # copy_file=True)

submission.add_additional_resource(description="Python conversion script used to build this submisson. Part of NUISANCE.",
    location="ToHepData.py",
    copy_file=True)

submission.add_additional_resource(description="Filter and Projection function examples. Can be run in the NUISANCE HepData Env 1.0.",
    location="analysis.cxx",
    copy_file=True)

submission.create_files("testout", remove_old=True)
