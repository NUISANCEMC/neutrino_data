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
  var.add_qualifier("suggested_flux", "inspire:167744")

def eventrate_from_text(
    filename,
    channel,
    reaction,
    signal,
    project,
    independent,
    dependent
):

  datafile = "nuisance_digitized_text/" + filename

  bins = []
  values = []
  for line in open(datafile):
    line = line.strip()
    print(line)
    bins.append(float(line.split(" ")[0]))
    values.append(float(line.split(" ")[1]))

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
  IndepVar = Variable(independent, is_independent=True, is_binned=True, units="[GeV^{2}/C^{2}]")
  IndepVar.values = data["xbins"]

  DepVar = Variable("DepVar", is_independent=False, is_binned=False, units="Event Counts")
  DepVar.values = data["values"]

  unc_data = Uncertainty("Statistical errors", is_symmetric=True)
  unc_data.values = data["errors"]

  # DepVar.add_uncertainty(unc_data)

  add_default_qualifiers(DepVar)
  DepVar.add_qualifier("Filter", signal)
  DepVar.add_qualifier(independent, project)

  xsTable = Table(f"{dependent}-{independent}")
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

submission = Submission()
submission.read_abstract("abstract.txt")


# // ANL CC1ppip
#include "ANL_CC1ppip_Evt_1DQ2_nu.h"
#include "ANL_CC1ppip_Evt_1DWNmu_nu.h"
#include "ANL_CC1ppip_Evt_1DWNpi_nu.h"
#include "ANL_CC1ppip_Evt_1DWmupi_nu.h"
#include "ANL_CC1ppip_Evt_1DcosmuStar_nu.h"
#include "ANL_CC1ppip_Evt_1DcosthAdler_nu.h"
#include "ANL_CC1ppip_Evt_1Dphi_nu.h"
#include "ANL_CC1ppip_Evt_1Dppi_nu.h"
#include "ANL_CC1ppip_Evt_1Dthpr_nu.h"
#include "ANL_CC1ppip_XSec_1DEnu_nu.h"
#include "ANL_CC1ppip_XSec_1DQ2_nu.h"

# // ANL CC1npip
#include "ANL_CC1npip_Evt_1DQ2_nu.h"
#include "ANL_CC1npip_Evt_1DWNmu_nu.h"
#include "ANL_CC1npip_Evt_1DWNpi_nu.h"
#include "ANL_CC1npip_Evt_1DWmupi_nu.h"
#include "ANL_CC1npip_Evt_1DcosmuStar_nu.h"
#include "ANL_CC1npip_Evt_1Dppi_nu.h"
#include "ANL_CC1npip_XSec_1DEnu_nu.h"

# // ANL CC1pi0
#include "ANL_CC1pi0_Evt_1DQ2_nu.h"
#include "ANL_CC1pi0_Evt_1DWNmu_nu.h"
#include "ANL_CC1pi0_Evt_1DWNpi_nu.h"
#include "ANL_CC1pi0_Evt_1DWmupi_nu.h"
#include "ANL_CC1pi0_Evt_1DcosmuStar_nu.h"
#include "ANL_CC1pi0_XSec_1DEnu_nu.h"


submission.add_table(eventrate_from_text(
  filename="CC1pip_on_n/ANL_CC1npip_WNmu_per_0.04GeV.txt",
  channel="CC1npip",
  reaction="MU- + n --> MU- + pip",
  signal="ANL_167744_CC1npip_Filter_CC1pi3Prong",
  project="ANL_167744_CC1npip_Project_WNMu",
  independent="WNMu",
  dependent="EventRate")
)

submission.add_table(eventrate_from_text(
  filename="CC1pip_on_n/ANL_CC1npip_WNmu_per_0.04GeV.txt",
  channel="CC1npip",
  reaction="MU- + n --> MU- + pip",
  signal="ANL_167744_CC1npip_Filter_CC1pi3Prong",
  project="ANL_167744_CC1npip_Project_WNMu",
  independent="WNMu",
  dependent="EventRate")
)

submission.add_table(eventrate_from_text(
  filename="CC1pip_on_n/ANL_CC1npip_Wmupi_per_0.02GeV.txt",
  channel="CC1npip",
  reaction="MU- + n --> MU- + pip",
  signal="ANL_167744_CC1npip_Filter_CC1pi3Prong",
  project="ANL_167744_CC1npip_Project_Wmupi",
  independent="WNpi",
  dependent="EventRate")
)

submission.add_table(eventrate_from_text(
  filename="CC1pip_on_n/ANL_CC1pip_on_n_noEvents_Q2_W14GeV.txt",
  channel="CC1npip",
  reaction="MU- + n --> MU- + pip",
  signal="ANL_167744_CC1npip_Filter_CC1pi3Prong",
  project="ANL_167744_CC1npip_Project_Q2_W14GeV",
  independent="WNpi",
  dependent="EventRate")
)


  nuisance_digitized_text/CC1pip_on_n/ANL_CC1npip_WNmu_per_0.04GeV.txt				nuisance_digitized_text/CC1pip_on_n/ANL_CC1pip_on_n_noEvents_Q2_noWcut_HighQ2Gone.txt
nuisance_digitized_text/CC1pip_on_n/ANL_CC1npip_WNpi_per_0.02GeV.txt				nuisance_digitized_text/CC1pip_on_n/ANL_CC1pip_on_n_noEvents_Q2_noWcut_firstQ2rem.txt
nuisance_digitized_text/CC1pip_on_n/ANL_CC1npip_Wmupi_per_0.04GeV.txt				nuisance_digitized_text/CC1pip_on_n/anl82-numu-cc1npip-14Wcut.txt
nuisance_digitized_text/CC1pip_on_n/ANL_CC1pip_on_n_noEvents_Q2_W14GeV.txt			nuisance_digitized_text/CC1pip_on_n/anl82-numu-cc1npip-16Wcut.txt
nuisance_digitized_text/CC1pip_on_n/ANL_CC1pip_on_n_noEvents_Q2_W14GeV_firstQ2rem.txt		nuisance_digitized_text/CC1pip_on_n/anl82-numu-cc1npip-noWcut.txt
nuisance_digitized_text/CC1pip_on_n/ANL_CC1pip_on_n_noEvents_Q2_W14GeV_rebin.txt		nuisance_digitized_text/CC1pip_on_n/anl82corr-numu-n-to-mu-n-piplus-lowW_edges.txt
nuisance_digitized_text/CC1pip_on_n/ANL_CC1pip_on_n_noEvents_Q2_W14GeV_rebin_firstQ2rem.txt	nuisance_digitized_text/CC1pip_on_n/anl82corr-numu-n-to-mu-n-piplus-noW_edges.txt
nuisance_digitized_text/CC1pip_on_n/ANL_CC1pip_on_n_noEvents_Q2_noWcut.txt



submission.add_additional_resource(description="Original digitised file used in NUISANCE.",
    location="nuisance_digitized_text.tar.gz",
    copy_file=True)

submission.add_additional_resource(description="Python conversion script used to build this submisson. Part of NUISANCE.",
    location="ToHepData.py",
    copy_file=True)

submission.add_additional_resource(description="Filter and Projection function examples. Can be run in the NUISANCE HepData Env 1.0.",
    location="analysis.cxx",
    copy_file=True)

submission.create_files("testout", remove_old=True)
