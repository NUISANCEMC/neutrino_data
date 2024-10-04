#include <optional>

int T2K_CC0Pi_CO_TKI_nu_2024_Select(HepMC3::GenEvent const &ev) {

  using namespace ps::pdg;
  using namespace ps::unit;

  if (!ps::event::has_beam_part(ev, kNuMu) ||
      !ps::event::has_exact_out_part(ev, kMuon, 1) ||
      !ps::event::has_out_part(ev, kProton)) {
    return false;
  }

  if (ps::event::num_out_part(ev, ps::pids(kPiPlus, kPiZero, kPiMinus),
                              ps::flatten)) {
    return false;
  }

  // a projector object that calculates a particles costheta WRT neutrino
  // direction
  auto costheta_nu = ps::costheta(ps::event::beam_part(ev, kNuMu)->momentum());

  auto muon_phase_space =
      (ps::p3mod > 0.225_GeV) && (ps::p3mod < 10_GeV) && (costheta_nu > -0.6);

  if (!muon_phase_space(ps::event::hm_out_part(ev, kMuon))) {
    return false;
  }

  auto proton_phase_space =
      (ps::p3mod > 0.525_GeV) && (ps::p3mod < 1.1_GeV) && (costheta_nu > 0.3);

  auto passing_protons = ps::part::filter(proton_phase_space,
                                          ps::event::all_out_part(ev, kProton));

  return passing_protons.size();
}

int T2K_CC0Pi_CO_TKI_nu_2024_Select_C(HepMC3::GenEvent const &ev) {
  if (!T2K_CC0Pi_CO_TKI_nu_2024_Select(ev)) {
    return false;
  }

  return ps::event::target_part(ev)->pid() == ps::pdg::Carbon12;
}

int T2K_CC0Pi_CO_TKI_nu_2024_Select_O(HepMC3::GenEvent const &ev) {
  if (!T2K_CC0Pi_CO_TKI_nu_2024_Select(ev)) {
    return false;
  }

  return ps::event::target_part(ev)->pid() == ps::pdg::Oxygen16;
}

std::pair<HepMC3::FourVector, HepMC3::FourVector>
T2K_CC0Pi_CO_TKI_nu_2024_GetMuPPairPT(HepMC3::GenEvent const &ev) {
  using namespace ps::pdg;

  auto nu_dir = ps::event::beam_part(ev, kNuMu)->momentum();

  return {ps::vect::transverse(ps::event::hm_out_part(ev, kMuon)->momentum(),
                               nu_dir),
          ps::vect::transverse(ps::event::hm_out_part(ev, kProton)->momentum(),
                               nu_dir)};
}

double
T2K_CC0Pi_CO_TKI_nu_2024_Project_DeltaPT_GeV_c(HepMC3::GenEvent const &ev) {

  using namespace ps::pdg;
  using namespace ps::unit;

  // need at least these to form kinematics
  if (!ps::event::has_beam_part(ev, kNuMu) ||
      !ps::event::has_out_part(ev, kMuon) ||
      !ps::event::has_out_part(ev, kProton)) {
    return ps::kMissingDatum<double>;
  }

  auto const &mu_p_pt = T2K_CC0Pi_CO_TKI_nu_2024_GetMuPPairPT(ev);

  return (mu_p_pt.first + mu_p_pt.second).p3mod() / ps::unit::GeV_c;
}

double
T2K_CC0Pi_CO_TKI_nu_2024_Project_DeltaAlphaT_deg(HepMC3::GenEvent const &ev) {
  using namespace ps::pdg;
  using namespace ps::unit;

  // need at least these to form kinematics
  if (!ps::event::has_beam_part(ev, kNuMu) ||
      !ps::event::has_out_part(ev, kMuon) ||
      !ps::event::has_out_part(ev, kProton)) {
    return ps::kMissingDatum<double>;
  }

  auto const &mu_p_pt = T2K_CC0Pi_CO_TKI_nu_2024_GetMuPPairPT(ev);
  auto const &mu_pt = mu_p_pt.first;

  auto dpt = (mu_pt + mu_p_pt.second);

  return std::acos(-ps::vect::dot(mu_pt, dpt) / (mu_pt.p3mod() * dpt.p3mod())) /
         ps::unit::deg;
}