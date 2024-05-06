#include "HepMC3/GenEvent.h"

extern "C" {

double T2K_CC0pi_XSec_2DPcos_nu_Project_CosThetaMu(HepMC3::GenEvent const &ev) {
  auto nu = ps::sel::Beam(ev, ps::pdg::kNuMu);
  auto mu = ps::sel::OutPartHM(ev, ps::pdg::kMuon);
  if (!nu || !mu) {
    return 0xdeadbeef;
  }
  
  return ps::proj::parts::CosTheta(nu, mu);
}

double T2K_CC0pi_XSec_2DPcos_nu_I_Project_PMu(HepMC3::GenEvent const &ev) {  
  auto mu = ps::sel::OutPartHM(ev, ps::pdg::kMuon);
  if (!mu) {
    return 0xdeadbeef;
  }
  return mu->momentum().length2() * ps::GeV;
}

bool T2K_CC0pi_XSec_2DPcos_nu_Filter(HepMC3::GenEvent const &ev) {
  auto nu = ps::sel::Beam(ev, ps::pdg::kNuMu);
  auto mu = ps::sel::OutPartHM(ev, ps::pdg::kMuon);
  if (!nu || !mu) {
    return false;
  }
  return true;
}

}
