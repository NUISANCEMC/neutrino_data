#!/usr/bin/env python3

import yaml
import csv
import os
import ROOT
import re

from hepdata_lib import Submission, Table, Variable, Uncertainty

NUISANCE_DATA_ROOT = os.environ.get("NUISANCE_DATA_ROOT")

if not NUISANCE_DATA_ROOT:
  print("[ERROR]: NUISANCE_DATA_ROOT is not set.")
  exit(1)

probe = "nu"
expt = "T2K"
target = "CH"
species = "numu"
ref = "PRD.93.112012"

release_dir = "/".join([NUISANCE_DATA_ROOT,"released", probe, expt, target, species, ref])
nuisance_dir = "/".join([NUISANCE_DATA_ROOT,"nuisance", probe, expt, target, species, ref])

INSPIRE_id=1421157

def add_default_qualifiers(var):
  var.add_qualifier("default_types", "FIX/FULL")
  var.add_qualifier("allowed_types", "FULL,DIAG/FREE,SHAPE,FIX/SYSTCOV/STATCOV")
  var.add_qualifier("allowed_targets", "C,H")
  var.add_qualifier("allowed_species", "numu")
  var.add_qualifier("suggested_flux", "inspire:xxxxxx")

def xsec_pmu_costhetamu_I():
  inFileName = "nd280data-numu-cc0pi-xs-on-c-2015/cross-section_analysisI.txt"
  inFile ="/".join([release_dir, inFileName])

  data = {
    "values": [],
  }
  CosThetaAxis = []
  PmuAxis = []
  line_num = 0
  first_line = 3

  with open(inFile) as file:
    for line in file:
        if (line_num >= first_line):
          splits = [ x.strip() for x in line.split("|") ]
          if (len(splits) != 4):
            continue
          CosThetaAxis.append(tuple([ float(x.strip()) for x in re.split(" - ", splits[1]) ]))
          PmuAxis.append(tuple([ float(x.strip()) for x in re.split(" - ", splits[2]) ]))
          data["values"].append(float(splits[3].strip()))

        line_num = line_num + 1


  #### Build Submission
  CosThetaVar = Variable("CosThetaMu", is_independent=True, is_binned=True, units="")
  CosThetaVar.values = CosThetaAxis

  PVar = Variable("PMu", is_independent=True, is_binned=True, units="GeV/c")
  PVar.values = PmuAxis

  CrossSection = Variable("CrossSection", is_independent=False, is_binned=False, units="$10^{-38}$ cm${}^{2}/GeV /Nucleon$")
  CrossSection.values = data["values"]

  add_default_qualifiers(CrossSection)
  CrossSection.add_qualifier("Filter", "T2K_CC0pi_XSec_2DPcos_nu_Filter")
  CrossSection.add_qualifier("CosThetaMu", "T2K_CC0pi_XSec_2DPcos_nu_Project_CosThetaMu")
  CrossSection.add_qualifier("PMu", "T2K_CC0pi_XSec_2DPcos_nu_I_Project_PMu")

  xsTable = Table("CrossSection-CosThetaMuPMu_AnalysisI")
  xsTable.description = """
  Something"""

  xsTable.add_variable(CosThetaVar)
  xsTable.add_variable(PVar)
  xsTable.add_variable(CrossSection)

  xsTable.keywords["observables"] = ["D2SIG/DCOSTHETA/DP"]
  xsTable.keywords["reactions"] = ["NUMU C --> MU- P"]
  xsTable.keywords["phrases"] = ["Neutrino CC0Pi", "Cross Section"]

  return xsTable

def covars_pmu_costhetamu_I():
  inFileName_flux = "nd280data-numu-cc0pi-xs-on-c-2015/covariance_fluxNormalizationSystematics_analysisI.txt"
  inFile_flux ="/".join([release_dir, inFileName_flux])

  inFileName_shape = "nd280data-numu-cc0pi-xs-on-c-2015/covariance_shapeSystematics_analysisI.txt"
  inFile_shape ="/".join([release_dir, inFileName_shape])

  inFileName_stats = "nd280data-numu-cc0pi-xs-on-c-2015/covariance_statisticUncertainty_analysisI.txt"
  inFile_stats ="/".join([release_dir, inFileName_stats])

  covars = {
    "flux": [],
    "shape": [],
    "stats": [],
    "total": [],
  }
  bini = []
  binj = []
  line_num = 0
  first_line = 3

  with open(inFile_flux) as file:
    for line in file:
        if (line_num >= first_line):
          splits = [ x.strip() for x in line.split("|") ]
          if (len(splits) != 2):
            continue
          binij = re.split(" - ", splits[0])
          bini.append(int(binij[0].strip()))
          binj.append(int(binij[1].strip()))
          covars["flux"].append(float(splits[1].strip()))

        line_num = line_num + 1

  line_num = 0
  with open(inFile_shape) as file:
    for line in file:
        if (line_num >= first_line):
          splits = [ x.strip() for x in line.split("|") ]
          if (len(splits) != 2):
            continue
          covars["shape"].append(float(splits[1].strip()))

        line_num = line_num + 1

  line_num = 0
  with open(inFile_stats) as file:
    for line in file:
        if (line_num >= first_line):
          splits = [ x.strip() for x in line.split("|") ]
          if (len(splits) != 2):
            continue
          covars["stats"].append(float(splits[1].strip()))

        line_num = line_num + 1

  #### Build Submission
  BiniVar = Variable("Bini", is_independent=True, is_binned=False, units="")
  BiniVar.values = bini

  BinjVar = Variable("Binj", is_independent=True, is_binned=False, units="")
  BinjVar.values = binj

  FluxNormCovariance = Variable("FluxNormCovariance", is_independent=False, is_binned=False, units="$10^{-38}$ cm${}^{2}/GeV /Nucleon$")
  FluxNormCovariance.values = covars["flux"]

  ShapeCovariance = Variable("ShapeCovariance", is_independent=False, is_binned=False, units="$10^{-38}$ cm${}^{2}/GeV /Nucleon$")
  ShapeCovariance.values = covars["shape"]

  StatsCovariance = Variable("StatsCovariance", is_independent=False, is_binned=False, units="$10^{-38}$ cm${}^{2}/GeV /Nucleon$")
  StatsCovariance.values = covars["stats"]

  if not ( (len(covars["flux"]) == len(covars["shape"])) 
      and (len(covars["flux"]) == len(covars["stats"]))):
    print("[ERROR]: found flux(%s), shape(%s), stats(%s) covariance entries but expect them to match." 
      % (len(covars["flux"]),len(covars["shape"]),len(covars["stats"])))
    exit(1)

  for i in range(len(covars["flux"])):
    covars["total"].append(covars["flux"][i]+covars["shape"][i]+covars["stats"][i])

  TotalCovariance = Variable("TotalCovariance", is_independent=False, is_binned=False, units="$10^{-38}$ cm${}^{2}/GeV /Nucleon$")
  TotalCovariance.values = covars["total"]

  covarTable = Table("Covar-CosThetaMuPMu_AnalysisI")
  covarTable.description = """
  Something"""

  covarTable.add_variable(BiniVar)
  covarTable.add_variable(BinjVar)
  covarTable.add_variable(FluxNormCovariance)
  covarTable.add_variable(ShapeCovariance)
  covarTable.add_variable(StatsCovariance)
  covarTable.add_variable(TotalCovariance)

  return covarTable

def xsec_pmu_costhetamu_II():
  inFileName = "nd280data-numu-cc0pi-xs-on-c-2015/rps_crossSection_analysis2.txt"
  inFile ="/".join([release_dir, inFileName])

  data = {
    "values": [],
  }
  CosThetaAxis = []
  PmuAxis = []
  line_num = 0
  first_line = 3

  with open(inFile) as file:
    for line in file:
        if (line_num >= first_line):
          splits = [ x.strip() for x in line.split("|") ]
          if (len(splits) != 4):
            continue
          splits = [ x.strip() for x in line.split("|") ]
          CosThetaAxis.append(tuple([ float(x.strip()) for x in re.split(" - ", splits[1]) ]))
          PmuAxis.append(tuple([ float(x.strip()) for x in re.split(" - ", splits[2]) ]))
          data["values"].append(float(splits[3].strip()))

        line_num = line_num + 1


  #### Build Submission
  CosThetaVar = Variable("CosThetaMu", is_independent=True, is_binned=True, units="")
  CosThetaVar.values = CosThetaAxis

  PVar = Variable("PMu", is_independent=True, is_binned=True, units="GeV/c")
  PVar.values = PmuAxis

  CrossSection = Variable("CrossSection", is_independent=False, is_binned=False, units="$10^{-38}$ cm${}^{2}/GeV /Nucleon$")
  CrossSection.values = data["values"]

  add_default_qualifiers(CrossSection)
  CrossSection.add_qualifier("Filter", "T2K_CC0pi_XSec_2DPcos_nu_Filter")
  CrossSection.add_qualifier("CosThetaMu", "T2K_CC0pi_XSec_2DPcos_nu_Project_CosThetaMu")
  CrossSection.add_qualifier("PMu", "T2K_CC0pi_XSec_2DPcos_nu_II_Project_PMu")

  xsTable = Table("CrossSection-CosThetaMuPMu_AnalysisII")
  xsTable.description = """
  Something"""

  xsTable.add_variable(CosThetaVar)
  xsTable.add_variable(PVar)
  xsTable.add_variable(CrossSection)

  xsTable.keywords["observables"] = ["D2SIG/DCOSTHETA/DP"]
  xsTable.keywords["reactions"] = ["NUMU C --> MU- P"]
  xsTable.keywords["phrases"] = ["Neutrino CC0Pi", "Cross Section"]

  return xsTable

submission = Submission()
submission.read_abstract("/".join([release_dir,"abstract.txt"]))
submission.add_table(xsec_pmu_costhetamu_I())
submission.add_table(covars_pmu_costhetamu_I())
submission.add_table(xsec_pmu_costhetamu_II())
submission.add_additional_resource(description="Python conversion script used to build this submisson. Part of NUISANCE.",
    location="/".join([nuisance_dir,"ToHepData.py"]),
    copy_file=True)
submission.add_additional_resource(description="Filter and Projection function examples. Can be run in the NUISANCE ProFilt Env 1.0.",
    location="/".join([nuisance_dir,"snippet.cxx"]),
    copy_file=True)
submission.add_additional_resource(description="Official data release from T2K",
    location="/".join([release_dir,"nd280data-numu-cc0pi-xs-on-c-2015.tar"]),
    copy_file=True)
submission.create_files("testout", remove_old=True)
