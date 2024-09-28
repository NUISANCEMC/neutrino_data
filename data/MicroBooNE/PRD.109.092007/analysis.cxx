int MicroBooNE_CC0Pi_GKI_nu_SelectSignal(HepMC3::GenEvent const &ev) {
  if (!ps::event::has_beam_part(ev, ps::pdg::kNuMu) ||
      !ps::event::has_exact_out_part(ev, ps::pdg::kMuon, 1) ||
      ps::event::num_out_part(ev, ps::pdg::kAMuon) ||
      ps::event::num_out_part(ev, ps::pdg::kPiZero) ||
      ps::event::num_out_part_except(
          ev, ps::pids(ps::pdg::kMuon, ps::pdg::kNeutron, ps::pdg::kProton,
                       ps::pdg::kPiPlus, ps::pdg::kPiMinus))) {
    return false;
  }

  auto nu = ps::event::beam_part(ev, ps::pdg::kNuMu);
  auto mu = ps::event::hm_out_part(ev, 13);

  auto p_mu_GeV = mu->momentum().p3mod() / ps::unit::GeV_c;

  if ((p_mu_GeV < 0.1) || (1.2 < p_mu_GeV)) {
    return false;
  }

  auto fs_protons = ps::event::all_out_part(ev, ps::pdg::kProton);
  auto fs_protons_in_ps = ps::part::filter(
      (ps::p3mod >= 300 * ps::unit::MeV_c) && (ps::p3mod < ps::unit::GeV_c),
      fs_protons);

  if (fs_protons_in_ps.size() != 1) {
    return false;
  }

  auto fs_cpi = ps::part::cat(ps::event::all_out_part(
      ev, ps::pids(ps::pdg::kPiPlus, ps::pdg::kPiMinus)));
  auto fs_cpi_observable =
      ps::part::filter(ps::p3mod >= 70 * ps::unit::MeV_c, fs_cpi);
  if (fs_cpi_observable.size()) {
    return false;
  }

  return true;
}

static int const kMicroBooNE_CC0Pi_GKI_nu_pn = 0;
static int const kMicroBooNE_CC0Pi_GKI_nu_alpha3d = 1;
static int const kMicroBooNE_CC0Pi_GKI_nu_phi3d = 2;
static int const kMicroBooNE_CC0Pi_GKI_nu_pn_para = 3;
static int const kMicroBooNE_CC0Pi_GKI_nu_pn_perp = 4;
static int const kMicroBooNE_CC0Pi_GKI_nu_pn_perp_x = 5;
static int const kMicroBooNE_CC0Pi_GKI_nu_pn_perp_y = 6;
static int const kMicroBooNE_CC0Pi_nu_dpt = 7;
static int const kMicroBooNE_CC0Pi_nu_dphit = 8;
static int const kMicroBooNE_CC0Pi_nu_dat = 9;
static int const kMicroBooNE_CC0Pi_GKI_nu_pn_perp_nkinematics = 10;

std::vector<double>
MicroBooNE_CC0Pi_GKI_nu_Kinematics(HepMC3::GenEvent const &ev) {

  std::vector<double> kinematics(kMicroBooNE_CC0Pi_GKI_nu_pn_perp_nkinematics,
                                 ps::kMissingDatum<double>);

  // can't construct kinematics without at least these particles
  if (!ps::event::has_beam_part(ev, ps::pdg::kNuMu) ||
      !ps::event::has_exact_out_part(ev, ps::pdg::kMuon, 1)) {
    return kinematics;
  }

  auto nu = ps::event::beam_part(ev, ps::pdg::kNuMu);
  auto nu_dir = ps::vect::direction(nu->momentum());

  auto mu = ps::event::hm_out_part(ev, 13);

  auto p_mu = mu->momentum();
  auto pt_mu = ps::vect::transverse(p_mu, nu_dir);

  auto fs_protons = ps::event::all_out_part(ev, ps::pdg::kProton);
  auto fs_protons_in_ps = ps::part::filter(
      (ps::p3mod >= 300 * ps::unit::MeV_c) && (ps::p3mod < ps::unit::GeV_c),
      fs_protons);

  if (!fs_protons_in_ps.size()) {
    return kinematics;
  }

  auto prot = ps::part::highest(ps::p3mod, fs_protons_in_ps);
  auto p_prot = prot->momentum();
  auto pt_prot = ps::vect::transverse(p_prot, nu_dir);

  auto dpt = (pt_prot + pt_mu);

  static double const Ar40EB = 30.9 * ps::unit::MeV;

  auto E_Cal = p_mu.e() + (p_prot.e() - p_prot.m()) + Ar40EB;
  auto P_L = nu_dir * (ps::vect::dot(p_mu, nu_dir) +
                       ps::vect::dot(p_prot, nu_dir) - E_Cal);
  auto qvec = (nu_dir * E_Cal) - p_mu;
  auto qvec_t = ps::vect::transverse(qvec, nu_dir);

  auto pn = P_L + dpt;
  auto phi3d =
      std::acos(ps::vect::dot(qvec, p_prot) / (qvec.p3mod() * p_prot.p3mod()));

  auto alpha3d =
      std::acos(ps::vect::dot(qvec, pn) / (qvec.p3mod() * pn.p3mod()));

  auto pn_perp_x =
      ps::vect::dot(ps::vect::cross(ps::vect::direction(qvec_t), nu_dir), pn);
  auto pn_perp_y = ps::vect::dot(
      ps::vect::cross(ps::vect::direction(qvec),
                      ps::vect::cross(ps::vect::direction(qvec_t), nu_dir)),
      pn);
  auto pn_perp = pn.p3mod() * std::sin(alpha3d);
  auto pn_para = pn.p3mod() * std::cos(alpha3d);

  kinematics[kMicroBooNE_CC0Pi_GKI_nu_pn] = pn.p3mod() / ps::unit::GeV_c;
  kinematics[kMicroBooNE_CC0Pi_GKI_nu_alpha3d] = alpha3d / ps::unit::deg;
  kinematics[kMicroBooNE_CC0Pi_GKI_nu_phi3d] = phi3d / ps::unit::deg;
  kinematics[kMicroBooNE_CC0Pi_GKI_nu_pn_para] = pn_para / ps::unit::GeV_c;
  kinematics[kMicroBooNE_CC0Pi_GKI_nu_pn_perp] = pn_perp / ps::unit::GeV_c;
  kinematics[kMicroBooNE_CC0Pi_GKI_nu_pn_perp_x] = pn_perp_x / ps::unit::GeV_c;
  kinematics[kMicroBooNE_CC0Pi_GKI_nu_pn_perp_y] = pn_perp_y / ps::unit::GeV_c;

  kinematics[kMicroBooNE_CC0Pi_nu_dpt] = dpt.p3mod() / ps::unit::GeV_c;

  auto phi = std::acos(ps::vect::dot(pt_mu, pt_prot) /
                       (pt_mu.p3mod() * pt_prot.p3mod()));
  auto alpha =
      std::acos(ps::vect::dot(pt_mu, dpt) / (pt_mu.p3mod() * dpt.p3mod()));

  kinematics[kMicroBooNE_CC0Pi_nu_dphit] = phi / ps::unit::deg;
  kinematics[kMicroBooNE_CC0Pi_nu_dat] = alpha / ps::unit::deg;

  return kinematics;
}

double MicroBooNE_CC0Pi_GKI_nu_pn(HepMC3::GenEvent const &ev) {
  return MicroBooNE_CC0Pi_GKI_nu_Kinematics(ev)[kMicroBooNE_CC0Pi_GKI_nu_pn];
}
double MicroBooNE_CC0Pi_GKI_nu_alpha3d(HepMC3::GenEvent const &ev) {
  return MicroBooNE_CC0Pi_GKI_nu_Kinematics(
      ev)[kMicroBooNE_CC0Pi_GKI_nu_alpha3d];
}
double MicroBooNE_CC0Pi_GKI_nu_phi3d(HepMC3::GenEvent const &ev) {
  return MicroBooNE_CC0Pi_GKI_nu_Kinematics(ev)[kMicroBooNE_CC0Pi_GKI_nu_phi3d];
}
double MicroBooNE_CC0Pi_GKI_nu_pn_para(HepMC3::GenEvent const &ev) {
  return MicroBooNE_CC0Pi_GKI_nu_Kinematics(
      ev)[kMicroBooNE_CC0Pi_GKI_nu_pn_para];
}
double MicroBooNE_CC0Pi_GKI_nu_pn_perp(HepMC3::GenEvent const &ev) {
  return MicroBooNE_CC0Pi_GKI_nu_Kinematics(
      ev)[kMicroBooNE_CC0Pi_GKI_nu_pn_perp];
}
double MicroBooNE_CC0Pi_GKI_nu_pn_perp_x(HepMC3::GenEvent const &ev) {
  return MicroBooNE_CC0Pi_GKI_nu_Kinematics(
      ev)[kMicroBooNE_CC0Pi_GKI_nu_pn_perp_x];
}
double MicroBooNE_CC0Pi_GKI_nu_pn_perp_y(HepMC3::GenEvent const &ev) {
  return MicroBooNE_CC0Pi_GKI_nu_Kinematics(
      ev)[kMicroBooNE_CC0Pi_GKI_nu_pn_perp_y];
}

double MicroBooNE_CC0Pi_nu_dpt(HepMC3::GenEvent const &ev) {
  return MicroBooNE_CC0Pi_GKI_nu_Kinematics(ev)[kMicroBooNE_CC0Pi_nu_dpt];
}

double MicroBooNE_CC0Pi_nu_dphit(HepMC3::GenEvent const &ev) {
  return MicroBooNE_CC0Pi_GKI_nu_Kinematics(ev)[kMicroBooNE_CC0Pi_nu_dphit];
}

double MicroBooNE_CC0Pi_nu_dat(HepMC3::GenEvent const &ev) {
  return MicroBooNE_CC0Pi_GKI_nu_Kinematics(ev)[kMicroBooNE_CC0Pi_nu_dat];
}