#include "FitEvent.h"

#include <cmath>
#include <iostream>
#include <string>

extern "C" {

bool T2K_CC0pi_XSec_2DPcos_nu_Filter(FitEvent *event) {

  int const nuPDG = 14;

  // Check for the desired PDG code
  if (!event->HasISParticle(nuPDG)) {
    return false;
  }
  // Check that the charged lepton we expect has been produced
  if (!event->HasFSParticle(nuPDG > 0 ? nuPDG - 1 : nuPDG + 1)) {
    return false;
  }

  // Veto event with mesons
  if (event->NumFSMesons() != 0) {
    return false;
  }

  // Veto events which don't have exactly 1 outgoing charged lepton
  if (event->NumFSLeptons() != 1) {
    return false;
  }

  return true;
}

bool T2K_CC0pi_XSec_2DPcos_nu_II_Filter(FitEvent *event) {

  if (!T2K_CC0pi_XSec_2DPcos_nu_Filter(event)) {
    return false;
  }
  

  return true;
}

double T2K_CC0pi_XSec_2DPcos_nu_Project_CosThetaMu(FitEvent *event) {

  auto nu = event->GetNeutrinoIn();
  auto mu = event->GetHMFSParticle(13);

  if (!nu || !mu) {
    return 0xdeadbeef;
  }

  return std::cos(nu->P3().Angle(mu->P3()));
}

double T2K_CC0pi_XSec_2DPcos_nu_I_Project_PMu(FitEvent *event) {
  auto mu = event->GetHMFSParticle(13);

  if (!mu) {
    return 0xdeadbeef;
  }

  const double to_GeV = 1E-3;
  return mu->p() * to_GeV;
}
}
