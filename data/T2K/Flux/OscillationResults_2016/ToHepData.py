#!/usr/bin/env python3

import os

from hepdata_lib import Submission, Table, Variable, Uncertainty, RootFileReader

ref = "PRD.91.072010"
INSPIRE_id=1343108

def flux_tables(inFileName, detname, tname):

  reader = RootFileReader(inFileName)
  
  fluxes =  { "nue": reader.read_hist_1d(f"enu_{detname}_nue"),
             "numu": reader.read_hist_1d(f"enu_{detname}_numu"),
             "nueb": reader.read_hist_1d(f"enu_{detname}_nueb"),
             "numub": reader.read_hist_1d(f"enu_{detname}_numub") }

  #### Build Submission
  EnuVar = Variable("e_nu", is_independent=True, is_binned=True, units="GeV")
  EnuVar.values = fluxes["nue"]["x_edges"]

  FluxTable = Table(tname)

  FluxTable.add_variable(EnuVar)

  for nuspec in fluxes.keys():

    FluxVar = Variable(f"flux_{nuspec}", is_independent=False, is_binned=False, units="/cm$^{2}$ /10$^{21}$ p.o.t /50 MeV$")
    FluxVar.values = fluxes[nuspec]["y"]

    FluxVar.add_qualifier("probe_particle", nuspec)
    FluxVar.add_qualifier("bin_content_type", "count_density")

    FluxTable.add_variable(FluxVar)

  return FluxTable

submission = Submission()
submission.read_abstract("abstract.txt")
submission.add_table(flux_tables("T2Kflux2016/t2kflux_2016_minus250kA.root", "nd280", "t2kflux_2016_minus250kA_nd280"))
submission.add_table(flux_tables("T2Kflux2016/t2kflux_2016_minus250kA.root", "sk", "t2kflux_2016_minus250kA_sk"))
submission.add_table(flux_tables("T2Kflux2016/t2kflux_2016_plus250kA.root", "nd280", "t2kflux_2016_plus250kA_nd280"))
submission.add_table(flux_tables("T2Kflux2016/t2kflux_2016_plus250kA.root", "sk", "t2kflux_2016_plus250kA_sk"))

submission.add_additional_resource(description="Python conversion script used to build this submisson. Part of NUISANCE.",
    location="ToHepData.py",
    copy_file=True)

submission.add_link(description="Zenodo data release", location="https://zenodo.org/records/5734267")

submission.create_files(f"submission-{INSPIRE_id}", remove_old=True)
