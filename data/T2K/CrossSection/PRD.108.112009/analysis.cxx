int T2K_CC0Pi_onoffaxis_nu_SelectSignal(HepMC3::GenEvent const &ev) {
  auto [numu, muon] = ps::sel::PrimaryLeptonsForNuCC(ev, ps::pdg::kNuMu);
  if (!muon) {
    return 0;
  }

  auto nleps =
      ps::sel::OutPartsAny(ev, {ps::pdg::kMuon, -ps::pdg::kAMuon}).size();

  if ((nleps != 1)) {
    return 0;
  }

  auto notherparts = ps::sel::OutPartsExceptAny(
      ev, {ps::pdg::kMuon, ps::pdg::kProton, ps::pdg::kNeutron}).size();

  if ((notherparts != 0)) {
    return 0;
  }

  return 1; // 0pi
}

double T2K_CC0Pi_onoffaxis_nu_Project_CosThetaMu(HepMC3::GenEvent const &ev) {
  auto [numu, muon] = ps::sel::PrimaryLeptonsForNuCC(ev, ps::pdg::kNuMu);
  if (!muon) {
    return ps::kMissingDatum<double>;
  }

  return std::cos(muon->momentum().theta());
}

double T2K_CC0Pi_onoffaxis_nu_Project_PMu(HepMC3::GenEvent const &ev) {
  auto [numu, muon] = ps::sel::PrimaryLeptonsForNuCC(ev, ps::pdg::kNuMu);
  if (!muon) {
    return ps::kMissingDatum<double>;
  }

  return muon->momentum().p3mod() / ps::GeV;
}
