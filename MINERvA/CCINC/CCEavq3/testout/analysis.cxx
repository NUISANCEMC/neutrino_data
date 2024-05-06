#include "HepMC3/GenEvent.h"
#include "ProSelecta/env/env.h"

extern "C" {

double MINERvA_CCINC_CCEavq3_Project_q3(HepMC3::GenEvent const &ev) {
  auto nu = ps::sel::Beam(ev, ps::pdg::kNuMu);
  if (!nu) return 0xdeadbeef;

  auto mu = ps::sel::OutPartHM(ev, ps::pdg::kMuon);
  if (!mu) return 0xdeadbeef;

  auto q = nu.momentum() - mu.momentum();
  double q0 = (q.E()) / 1.E3;
  double thmu = muon->momentum().Vect().Angle(neutrino->momentum().Vect());
  double pmu  = muon->momentum().Vect().Mag() / 1.E3;
  double emu  = muon->momentum().e() / 1.E3;
  double mmu  = muon->momentum().mag() / 1.E3;

  // Get Enu Rec
  double enu_rec = emu + q0;

  // Set Q2 QE
  double q2qe = 2 * enu_rec * (emu - pmu * cos(thmu)) - mmu * mmu;

  // Calc Q3 from Q2QE and EnuTree
  q3 = sqrt(q2qe + q0 * q0);

  return q3 / ps::units::GeV;
}

double MINERvA_CCINC_CCEavq3_Project_Eav(HepMC3::GenEvent const &ev) {
  auto nu = ps::sel::Beam(ev, ps::pdg::kNuMu);
  if (!nu) return 0xdeadbeef;

  auto mu = ps::sel::OutPartHM(ev, ps::pdg::kMuon);
  if (!mu) return 0xdeadbeef;

  auto parts_ke = ps::sel::OutPartAny(ev,
    [ps::pdg::kProton, ps::pdg::kPiPlus, ps::pdg::kPiMinus]);
  auto parts_te = ps::sel::OutPartAny(ev,
    [ps::pdg::kPiZero, ps::pdg::kElectron, ps::pdg::kAElectron, ps::pdg::kGamma]);

  double Eav = 0.0;
  for (auto const & part : parts_ke) {
    Eav += part.momentum().E();
  }
  for (auto const & part : parts_te) {
    Eav += part.momentum().E();
  }

  return Eav / ps::units::GeV;
}

int MINERvA_CCINC_CCEavq3_Filter(HepMC3::GenEvent const &ev) {
  auto nu = ps::sel::Beam(ev, ps::pdg::kNuMu);
  if (!nu) return false;

  auto mu = ps::sel::OutPartHM(ev, ps::pdg::kMuon);
  if (!mu) return false;

  double angle = ps::proj::event::CosLep(ev);
  if (cos(angle) <  0.93969262078) return false;
  
  if (ps::proj::event::ELep(ev) < 1.5 * ps::units::GeV) return false;
  return true;
}

}

