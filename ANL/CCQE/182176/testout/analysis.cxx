#include "HepMC3/GenEvent.h"

extern "C" {

double ANL_CCQE_182176_Project_Q2(HepMC3::GenEvent const &ev) {
  auto nu = ps::sel::Beam(ev, ps::pdg::kNuMu);
  auto mu = ps::sel::OutPartHM(ev, ps::pdg::kMuon);
  if (!nu || !mu) {
    return 0xdeadbeef;
  }
  
  return (mu->momentum() - nu->momentum()).length2();
}

bool ANL_CCQE_182176_Filter(HepMC3::GenEvent const &ev) {
  auto nu = ps::sel::Beam(ev, ps::pdg::kNuMu);
  auto mu = ps::sel::OutPartHM(ev, ps::pdg::kMuon);
  if (!nu || !mu) {
    return false;
  }
  // Check CCQE
  return NuHepMC::ER3::GetProcessID(et) == 200;
}

}
