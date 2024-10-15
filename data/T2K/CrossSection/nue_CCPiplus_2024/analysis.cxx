int T2K_nueCCPiplus_2024_Select(HepMC3::GenEvent const &ev) {
  using namespace ps::pdg;

  if (!ps::event::has_beam_part(ev, kNuE) ||
      !ps::event::has_at_least_out_part(ev, ps::pids(kElectron, kPiPlus),
                                        {1, 1})) {
    return false;
  }

  return true;
}

double T2K_nueCCPiplus_2024_ElectronMomentum(HepMC3::GenEvent const &ev) {
  using namespace ps::pdg;

  if (!ps::event::has_beam_part(ev, kNuE) ||
      !ps::event::has_out_part(ev, kElectron)) {
    return ps::kMissingDatum<double>;
  }

  return ps::event::hm_out_part(ev, kElectron)->momentum().p3mod() /
         ps::unit::GeV;
}

double T2K_nueCCPiplus_2024_CosThetaElectron(HepMC3::GenEvent const &ev) {
  using namespace ps::pdg;

  if (!ps::event::has_beam_part(ev, kNuE) ||
      !ps::event::has_out_part(ev, kElectron)) {
    return ps::kMissingDatum<double>;
  }

  return ps::costheta(ps::event::hm_out_part(ev, kElectron));
}

double T2K_nueCCPiplus_2024_PiPlusMomentum(HepMC3::GenEvent const &ev) {
  using namespace ps::pdg;

  if (!ps::event::has_beam_part(ev, kNuE) ||
      !ps::event::has_out_part(ev, kPiPlus)) {
    return ps::kMissingDatum<double>;
  }

  return ps::event::hm_out_part(ev, kPiPlus)->momentum().p3mod() /
         ps::unit::GeV;
}
