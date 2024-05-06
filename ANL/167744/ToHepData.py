#!/usr/bin/env python3

import yaml
import csv
import os
import re
import numpy as np

from hepdata_lib import Submission, Table, Variable, Uncertainty
from subprocess import call

probe = "nu"
expt = "ANL"
target = "D2"
species = "numu"
ref = "Phys.Rev.D 25 (1982) 1161-1173"

INSPIRE_id=167744

def add_default_qualifiers(var):
  var.add_qualifier("likelihood", "Poisson")
  var.add_qualifier("scaling", "EventRateScaleToData")
  var.add_qualifier("weighting", "Default")
  var.add_qualifier("targets", "D2")
  var.add_qualifier("species", "numu")
  var.add_qualifier("suggested_flux", "inspire:unknown")


def eventrate_from_text(
    filename,
    channel,
    reaction,
    signal,
    project,
    independent,
    dependent,
    independent_unit,
    dependent_unit
):

  datafile = "nuisance_files/data/" + filename

  bins = []
  values = []
  for line in open(datafile):
    line = line.strip()
    print(line)
    try:
      bins.append(float(line.split(" ")[0]))
      values.append(float(line.split(" ")[1]))
    except:
      bins.append(float(line.split("\t")[0]))
      values.append(float(line.split("\t")[1]))

  data = {
    "values": [],
    "xbins": [],
    "errors": []
  }

  for i in range(len(values)-1):
    if values[i] <= 0.0: continue
    data["values"].append( values[i] )
    data["errors"].append( np.sqrt( values[i] ) )
    data["xbins"].append( [bins[i], bins[i+1]] )
  
  #### Build Submission
  IndepVar = Variable(independent, is_independent=True, is_binned=True, units=independent_unit)
  IndepVar.values = data["xbins"]

  DepVar = Variable("DepVar", is_independent=False, is_binned=False, units=dependent_unit)
  DepVar.values = data["values"]

  unc_data = Uncertainty("Statistical errors", is_symmetric=True)
  unc_data.values = data["errors"]

  # DepVar.add_uncertainty(unc_data)

  add_default_qualifiers(DepVar)
  DepVar.add_qualifier("Filter", signal)
  DepVar.add_qualifier(independent, project)

  xsTable = Table(f"{channel}-{dependent}-{independent}".lower())
  xsTable.description = f"""
  {dependent} for {channel} Measurement binned in {independent}
  Comparisons should be shape only comparing to
  statistical errors.
  """

  xsTable.add_variable(IndepVar)
  xsTable.add_variable(DepVar)

  xsTable.keywords["observables"] = ["EVENTS"]
  xsTable.keywords["reactions"] = [reaction]
  xsTable.keywords["phrases"] = ["Neutrino", channel, independent, dependent]

  return xsTable


def xsec_from_text(
    filename,
    channel,
    reaction,
    signal,
    project,
    independent,
    dependent,
    independent_unit,
    dependent_unit
):

  datafile = "nuisance_files/data/" + filename

  bins = []
  values = []
  errors = []
  for line in open(datafile):
    line = line.strip()
    line = line.replace("  ", " ")
    try:
      bins.append(float(line.split(" ")[0]))
      values.append(float(line.split(" ")[1]))
      errors.append(float(line.split(" ")[2]))
    except:
      bins.append(float(line.split("\t")[0]))
      values.append(float(line.split("\t")[1]))
      errors.append(float(line.split("\t")[2]))


  data = {
    "values": [],
    "xbins": [],
    "errors": []
  }

  for i in range(len(values)-1):
    if values[i] <= 0.0: continue
    data["values"].append( values[i] )
    data["errors"].append( errors[i] ) 
    data["xbins"].append( [bins[i], bins[i+1]] )
  
  #### Build Submission
  IndepVar = Variable(independent, is_independent=True, is_binned=True, units=independent_unit)
  IndepVar.values = data["xbins"]

  DepVar = Variable("DepVar", is_independent=False, is_binned=False, units=dependent_unit)
  DepVar.values = data["values"]

  unc_data = Uncertainty("Statistical errors", is_symmetric=True)
  unc_data.values = data["errors"]

  # DepVar.add_uncertainty(unc_data)

  DepVar.add_qualifier("likelihood", "SimpleLikelihood")
  DepVar.add_qualifier("scaling", "FATXNormalized")
  DepVar.add_qualifier("weighting", "Default")
  DepVar.add_qualifier("targets", "D2")
  DepVar.add_qualifier("species", "numu")
  DepVar.add_qualifier("suggested_flux", "inspire:unknown")

  DepVar.add_qualifier("Filter", signal)
  DepVar.add_qualifier(independent, project)

  xsTable = Table(f"{channel}-{dependent}-{independent}".lower())
  xsTable.description = f"""
  {dependent} for {channel} Measurement binned in {independent}
  Comparisons should be shape only comparing to
  statistical errors.
  """

  xsTable.add_variable(IndepVar)
  xsTable.add_variable(DepVar)

  xsTable.keywords["observables"] = ["dSIGMA"]
  xsTable.keywords["reactions"] = [reaction]
  xsTable.keywords["phrases"] = ["Neutrino", channel, independent, dependent]

  return xsTable

submission = Submission()
submission.read_abstract("abstract.txt")

# CC1npip
submission.add_table(eventrate_from_text(
  filename="CC1pip_on_n/ANL_CC1npip_cosmuStar.csv",
  channel="CC1npip",
  reaction="MU- + n --> MU- + pip",
  signal="ANL_167744_CC1npip_Selection",
  project="ANL_167744_CC1npip_Project_cosmuStar",
  independent="cosmuStar",
  independent_unit="cos(#theta*)",
  dependent="EventRate",
  dependent_unit="Event Counts")
)

submission.add_table(eventrate_from_text(
  filename="CC1pip_on_n/ANL_ppi_CC1npip.csv",
  channel="CC1npip",
  reaction="MU- + n --> MU- + pip",
  signal="ANL_167744_CC1npip_Selection",
  project="ANL_167744_CC1npip_Project_ppi",
  independent="ppi",
  independent_unit="p_{#pi} (GeV/c)",
  dependent="EventRate",
  dependent_unit="Event Counts")
)

submission.add_table(eventrate_from_text(
  filename="CC1pip_on_n/ANL_CC1pip_on_n_noEvents_Q2_W14GeV_rebin_firstQ2rem.txt",
  channel="CC1npip",
  reaction="MU- + n --> MU- + pip",
  signal="ANL_167744_CC1npip_Selection_W14Cut",
  project="ANL_167744_CC1npip_Project_Q2_W14",
  independent="Q2-lowW",
  independent_unit="Q^{2}_{rec} (GeV/c^{2})",
  dependent="EventRate",
  dependent_unit="Event Counts")
)

submission.add_table(eventrate_from_text(
  filename="CC1pip_on_n/ANL_CC1pip_on_n_noEvents_Q2_noWcut_HighQ2Gone.txt",
  channel="CC1npip",
  reaction="MU- + n --> MU- + pip",
  signal="ANL_167744_CC1npip_Selection",
  project="ANL_167744_CC1npip_Project_Q2_HighW",
  independent="Q2",
  independent_unit="Q^{2}_{rec} (GeV/c^{2})",
  dependent="EventRate",
  dependent_unit="Event Counts")
)


submission.add_table(eventrate_from_text(
  filename="CC1pip_on_n/ANL_CC1npip_Wmupi_per_0.04GeV.txt",
  channel="CC1npip",
  reaction="MU- + n --> MU- + pip",
  signal="ANL_167744_CC1npip_Selection",
  project="ANL_167744_CC1npip_Project_Wmupi",
  independent="Wmupi",
  independent_unit="GeV/c^{2}",
  dependent="EventRate",
  dependent_unit="Event Counts")
)


submission.add_table(eventrate_from_text(
  filename="CC1pip_on_n/ANL_CC1npip_WNmu_per_0.04GeV.txt",
  channel="CC1npip",
  reaction="MU- + n --> MU- + pip",
  signal="ANL_167744_CC1npip_Selection",
  project="ANL_167744_CC1npip_Project_WNmu",
  independent="WNmu",
  independent_unit="WNmu GeV/c^{2}",
  dependent="EventRate",
  dependent_unit="Event Counts")
)

submission.add_table(eventrate_from_text(
  filename="CC1pip_on_n/ANL_CC1npip_WNpi_per_0.02GeV.txt",
  channel="CC1npip",
  reaction="MU- + n --> MU- + pip",
  signal="ANL_167744_CC1npip_Selection",
  project="ANL_167744_CC1npip_Project_WNpi",
  independent="WNpi",
  independent_unit="GeV/c^{2}",
  dependent="EventRate",
  dependent_unit="Event Counts")
)


#### CC1pi0
submission.add_table(eventrate_from_text(
  filename="CC1pi0_on_n/ANL_CC1pi0_cosmuStar.csv",
  channel="CC1npi0",
  reaction="MU- + n --> MU- + pi0",
  signal= "ANL_167744_CC1npi0_Selection",
  project="ANL_167744_CC1npi0_Project_cosmuStar",
  independent="cosmustar",
  independent_unit="",
  dependent="EventRate",
  dependent_unit="Event Counts")
)

submission.add_table(eventrate_from_text(
  filename="CC1pi0_on_n/ANL_CC1pi0_on_n_noEvents_Q2_W14GeV_rebin_firstQ2rem.txt",
  channel="CC1npi0",
  reaction="MU- + n --> MU- + pi0",
  signal="ANL_167744_CC1npi0_Selection_lowW",
  project="ANL_167744_CC1npi0_Project_Q2_lowW",
  independent="Q2-lowW",
  independent_unit="GeV^{2}",
  dependent="EventRate",
  dependent_unit="Event Counts")
)

submission.add_table(eventrate_from_text(
  filename="CC1pi0_on_n/ANL_CC1pi0_on_n_noEvents_Q2_noWcut_HighQ2Gone.txt",
  channel="CC1npi0",
  reaction="MU- + n --> MU- + pi0",
  signal="ANL_167744_CC1npi0_Selection",
  project="ANL_167744_CC1npi0_Project_Q2",
  independent="Q2",
  independent_unit="GeV^{2}",
  dependent="EventRate",
  dependent_unit="Event Counts")
)

submission.add_table(eventrate_from_text(
  filename="CC1pi0_on_n/ANL_CC1pi0_Wmupi_per_0.04GeV.txt",
  channel="CC1npi0",
  reaction="MU- + n --> MU- + pi0",
  signal="ANL_167744_CC1npi0_Selection",
  project="ANL_167744_CC1npi0_Project_Wmupi",
  independent="Wmupi",
  independent_unit="GeV",
  dependent="EventRate",
  dependent_unit="Event Counts")
)

submission.add_table(eventrate_from_text(
  filename="CC1pi0_on_n/ANL_CC1pi0_WNmu_per_0.04GeV.txt",
  channel="CC1npi0",
  reaction="MU- + n --> MU- + pi0",
  signal="ANL_167744_CC1npi0_Selection",
  project="ANL_167744_CC1npi0_Project_WNmu",
  independent="WNmu",
  independent_unit="GeV",
  dependent="EventRate",
  dependent_unit="Event Counts")
)

submission.add_table(eventrate_from_text(
  filename="CC1pi0_on_n/ANL_CC1pi0_WNpi_per_0.03GeV.txt",
  channel="CC1npi0",
  reaction="MU- + n --> MU- + pi0",
  signal="ANL_167744_CC1npi0_Selection",
  project="ANL_167744_CC1npi0_Project_WNpi",
  independent="WNpi",
  independent_unit="GeV",
  dependent="EventRate",
  dependent_unit="Event Counts")
)

#### CC1pip
submission.add_table(eventrate_from_text(
  filename="CC1pip_on_p/ANL_CC1pip_on_p_noEvents_cosmuStar_1982.csv",
  channel="CC1ppip",
  reaction="MU- + p --> MU- + pip",
  signal= "ANL_167744_CC1ppip_Selection",
  project="ANL_167744_CC1ppip_Project_cosmuStar",
  independent="cosmustar",
  independent_unit="",
  dependent="EventRate",
  dependent_unit="Event Counts")
)

submission.add_table(eventrate_from_text(
  filename="CC1pip_on_p/ANL_CC1pip_on_p_noEvents_costhAdler_1982.csv",
  channel="CC1ppip",
  reaction="MU- + p --> MU- + pip",
  signal= "ANL_167744_CC1ppip_Selection_lowW",
  project="ANL_167744_CC1ppip_Project_costhAdler",
  independent="costhAdler",
  independent_unit="",
  dependent="EventRate",
  dependent_unit="Event Counts")
)

submission.add_table(eventrate_from_text(
  filename="CC1pip_on_p/ANL_CC1pip_on_p_noEvents_phiAdler_1982.csv",
  channel="CC1ppip",
  reaction="MU- + p --> MU- + pip",
  signal= "ANL_167744_CC1ppip_Selection_lowW",
  project="ANL_167744_CC1ppip_Project_phi",
  independent="phiAdler",
  independent_unit="",
  dependent="EventRate",
  dependent_unit="Event Counts")
)

submission.add_table(eventrate_from_text(
  filename="CC1pip_on_p/ANL_CC1pip_on_p_noEvents_ppi.csv",
  channel="CC1ppip",
  reaction="MU- + p --> MU- + pip",
  signal= "ANL_167744_CC1ppip_Selection_lowW",
  project="ANL_167744_CC1ppip_Project_ppi",
  independent="ppi",
  independent_unit="",
  dependent="EventRate",
  dependent_unit="Event Counts")
)


submission.add_table(eventrate_from_text(
  filename="CC1pip_on_p/ANL_CC1pip_on_p_noEvents_Q2_W14GeV_rebin_firstQ2rem.txt",
  channel="CC1ppip",
  reaction="MU- + p --> MU- + pip",
  signal= "ANL_167744_CC1ppip_Selection_lowW",
  project="ANL_167744_CC1ppip_Project_Q2_lowW",
  independent="Q2-lowW",
  independent_unit="",
  dependent="EventRate",
  dependent_unit="Event Counts")
)

submission.add_table(eventrate_from_text(
  filename="CC1pip_on_p/ANL_CC1pip_on_p_noEvents_Q2_noW_HighQ2Gone.txt",
  channel="CC1ppip",
  reaction="MU- + p --> MU- + pip",
  signal= "ANL_167744_CC1ppip_Selection",
  project="ANL_167744_CC1ppip_Project_Q2",
  independent="Q2",
  independent_unit="",
  dependent="EventRate",
  dependent_unit="Event Counts")
)

submission.add_table(eventrate_from_text(
  filename="CC1pip_on_p/ANL_CC1pip_on_p_noEvents_thProt.csv",
  channel="CC1ppip",
  reaction="MU- + p --> MU- + pip",
  signal= "ANL_167744_CC1ppip_Selection_lowW",
  project="ANL_167744_CC1ppip_Project_thpr",
  independent="thProt",
  independent_unit="",
  dependent="EventRate",
  dependent_unit="Event Counts")
)

submission.add_table(eventrate_from_text(
  filename="CC1pip_on_p/ANL_CC1ppip_Wmupi_per_0.02GeV.txt",
  channel="CC1ppip",
  reaction="MU- + p --> MU- + pip",
  signal= "ANL_167744_CC1ppip_Selection",
  project="ANL_167744_CC1ppip_Project_Wmupi",
  independent="Wmupi",
  independent_unit="",
  dependent="EventRate",
  dependent_unit="Event Counts")
)

submission.add_table(eventrate_from_text(
  filename="CC1pip_on_p/ANL_CC1ppip_WNmu_per_0.04GeV.txt",
  channel="CC1ppip",
  reaction="MU- + p --> MU- + pip",
  signal= "ANL_167744_CC1ppip_Selection",
  project="ANL_167744_CC1ppip_Project_WNmu",
  independent="WNmu",
  independent_unit="",
  dependent="EventRate",
  dependent_unit="Event Counts")
)

submission.add_table(eventrate_from_text(
  filename="CC1pip_on_p/ANL_CC1ppip_WNpi_per_0.02GeV.txt",
  channel="CC1ppip",
  reaction="MU- + p --> MU- + pip",
  signal= "ANL_167744_CC1ppip_Selection",
  project="ANL_167744_CC1ppip_Project_WNpi",
  independent="WNpi",
  independent_unit="",
  dependent="EventRate",
  dependent_unit="Event Counts")
)

submission.add_table(xsec_from_text(
  filename="CC1pip_on_p/ANL_CC1pip_on_p_dSigdQ2_W14_1982.txt",
  channel="CC1ppip",
  reaction="MU- + p --> MU- + pip",
  signal= "ANL_167744_CC1ppip_Selection_lowW",
  project="ANL_167744_CC1ppip_Project_Q2",
  independent="Q2",
  independent_unit="",
  dependent="dSigma",
  dependent_unit="XSec")
)

call(["tar","-czvf","nuisance_files.tar.gz","nuisance_files"])
submission.add_additional_resource(description="Original digitised file used in NUISANCE.",
    location="nuisance_files.tar.gz",
    copy_file=True)

submission.add_additional_resource(description="Filter and Projection function examples. Can be run in nuisance3.",
    location="analysis.cxx",
    copy_file=True)

#submission.add_additional_resource(description="Python conversion script used to build this submisson. Part of NUISANCE.",
#     location="ToHepData.py",
#     copy_file=True)

submission.create_files("entry", remove_old=True)
