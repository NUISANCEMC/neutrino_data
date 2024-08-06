#!/usr/bin/env python3

import os, pprint

from hepdata_lib import Submission, Table, Variable, Uncertainty, RootFileReader

ref = "EPJC.83.782"
INSPIRE_id=2638628

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

    FluxVar.add_qualifier("probe_species", nuspec)
    FluxVar.add_qualifier("bin_content_type", "count_per_bin_width")

    FluxTable.add_variable(FluxVar)

  return FluxTable

def oscillation_tables(sub):
  bayes_reader = RootFileReader("T2K_arxiv2303/Bayesian_DataRelease.root")

  gr2D_th23_dCP_woRC_both_cred68_0 = bayes_reader.read_graph("gr2D_th23_dCP_woRC_both_cred68_0")

  th23 = Variable("th23", is_independent=True, is_binned=False, units="rad")
  th23.values = gr2D_th23_dCP_woRC_both_cred68_0["x"]
  dcp = Variable("dcp", is_independent=True, is_binned=False, units="rad")
  dcp.values = gr2D_th23_dCP_woRC_both_cred68_0["y"]

  th23_dCP_woRC_both_cred68_Table = Table("th23_dCP_woRC_both_cred68")
  th23_dCP_woRC_both_cred68_Table.add_variable(th23)
  th23_dCP_woRC_both_cred68_Table.add_variable(dcp)

  contour = Variable("contour_68", is_independent=False, is_binned=False)
  contour.values = [68,] * len(th23.values)
  th23_dCP_woRC_both_cred68_Table.add_variable(contour)

  sub.add_table(th23_dCP_woRC_both_cred68_Table)

submission = Submission()

# submission.read_abstract("abstract.txt")

submission.add_table(flux_tables("t2kflux_2020_public_release/t2kflux_2020_plus250kA_nominal.root", "nd280", "t2kflux_2020_plus250kA_nd280_nominal"))
submission.add_table(flux_tables("t2kflux_2020_public_release/t2kflux_2020_plus250kA_nominal.root", "sk", "t2kflux_2020_plus250kA_sk_nominal"))
submission.add_table(flux_tables("t2kflux_2020_public_release/t2kflux_2020_plus250kA_runcond.root", "nd280", "t2kflux_2020_plus250kA_nd280_runcond"))
submission.add_table(flux_tables("t2kflux_2020_public_release/t2kflux_2020_plus250kA_runcond.root", "sk", "t2kflux_2020_plus250kA_sk_runcond"))
submission.add_table(flux_tables("t2kflux_2020_public_release/t2kflux_2020_minus250kA_nominal.root", "nd280", "t2kflux_2020_minus250kA_nd280_nominal"))
submission.add_table(flux_tables("t2kflux_2020_public_release/t2kflux_2020_minus250kA_nominal.root", "sk", "t2kflux_2020_minus250kA_sk_nominal"))
submission.add_table(flux_tables("t2kflux_2020_public_release/t2kflux_2020_minus250kA_runcond.root", "nd280", "t2kflux_2020_minus250kA_nd280_runcond"))
submission.add_table(flux_tables("t2kflux_2020_public_release/t2kflux_2020_minus250kA_runcond.root", "sk", "t2kflux_2020_minus250kA_sk_runcond"))

oscillation_tables(submission)

submission.add_additional_resource(description="Python conversion script used to build this submisson. Part of NUISANCE.",
    location="ToHepData.py",
    copy_file=True)

submission.add_link(description="Zenodo data release (Flux)", location="https://zenodo.org/records/5734307")
submission.add_link(description="Zenodo data release (Oscillation)", location="https://zenodo.org/records/7741399")

submission.create_files(f"submission-{INSPIRE_id}", remove_old=True)
