#!/usr/bin/env python3

import yaml
import csv
import os
import ROOT
import re
import numpy as np

from hepdata_lib import Submission, Table, Variable, Uncertainty

probe = "nu"
expt = "ANL"
target = "D2"
species = "numu"
ref = "PRD.26.537"

INSPIRE_id=182176

def add_default_qualifiers(var):
  var.add_qualifier("default_types", "SHAPE/DIAG")
  var.add_qualifier("allowed_types", "SHAPE/DIAG")
  var.add_qualifier("allowed_targets", "D2")
  var.add_qualifier("allowed_species", "numu")
  var.add_qualifier("suggested_flux", "inspire:182176")

def eventrate_q2():
  inFile ="nuisance_digitized_file.root"

  f = ROOT.TFile(inFile)
  h = f.Get("ANL_1DQ2_Data")
  
  data = {
    "values": [],
    "xbins": [],
    "errors": []
  }

  for i in range(h.GetNbinsX()):
    if h.GetBinContent(i+1) <= 0: continue
    data["values"].append( h.GetBinContent(i+1) )
    data["errors"].append( np.sqrt( h.GetBinContent(i+1) ) )
    data["xbins"].append( [h.GetXaxis().GetBinLowEdge(i+1), h.GetXaxis().GetBinLowEdge(i+2)] )
   
  #### Build Submission
  Q2Var = Variable("Q2", is_independent=True, is_binned=True, units="[GeV^{2}/C^{2}]")
  Q2Var.values = data["xbins"]

  EventCounts = Variable("EventCounts", is_independent=False, is_binned=False, units="Event Counts")
  EventCounts.values = data["values"]

  unc_data = Uncertainty("Statistical errors", is_symmetric=True)
  unc_data.values = data["errors"]

  # EventCounts.add_uncertainty(unc_data)

  add_default_qualifiers(EventCounts)
  EventCounts.add_qualifier("Filter", "ANL_CCQE_182176_Filter")
  EventCounts.add_qualifier("Q2", "ANL_CCQE_182176_Project_Q2")

  xsTable = Table("EventCounts-Q2")
  xsTable.description = """
  Event Counts for ANL Measurement binned in Q2.
  Comparisons should be shape only comparing to
  statistical errors.
  """

  xsTable.add_variable(Q2Var)
  xsTable.add_variable(EventCounts)

  xsTable.keywords["observables"] = ["EVENTS"]
  xsTable.keywords["reactions"] = ["NUMU D2 --> MU- D2"]
  xsTable.keywords["phrases"] = ["Neutrino", "CCQE", "Event Rate"]

  return xsTable

submission = Submission()
submission.read_abstract("abstract.txt")
submission.add_table(eventrate_q2())

submission.add_additional_resource(description="Original digitised file used in NUISANCE.",
    location="nuisance_digitized_file.root",
    copy_file=True)

submission.add_additional_resource(description="Python conversion script used to build this submisson. Part of NUISANCE.",
    location="ToHepData.py",
    copy_file=True)

submission.add_additional_resource(description="Filter and Projection function examples. Can be run in the NUISANCE HepData Env 1.0.",
    location="analysis.cxx",
    copy_file=True)

submission.create_files("testout", remove_old=True)
