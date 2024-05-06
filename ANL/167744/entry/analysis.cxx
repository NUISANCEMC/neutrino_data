#include "HepMC3/GenEvent.h"
#include "HepMC3/FourVector.h"
extern "C" {

// HELPERS
FourVector BoostVector(const FourVector& vecin, const FourVector& boost) {
    return vecin;
}

double CosAngle(const FourVector& v1, const FourVector &v2) {
    return 0.0;
}

double HadronicMass_NucleonPi(const FourVector&  Pn, const FourVector& Ppip) {
    return 0.0;
}

bool isCC1pi3Prong(HepMC3::GenEvent const &ev,
    int nuPDG, int piPDG, int thirdPDG,
    double EnuMin, double EnuMax) {

    auto nu = ps::sel::Beam(ev, nuPDG);
    auto in = ps::sel::Target(ev, thirdPDG);

    auto mu = ps::sel::OutPartHM(ev, nuPDG > 0 ? nuPDG-1 : nuPDG+1);
    auto pip = ps::sel::OutPartHM(ev, piPDG);
    auto n = ps::sel::OutPartHM(ev, thirdPDG);

    if (!nu || !in || !mu || !pip || !n) {
        return 0xdeadbeef;
    }

    if (EnuMin != EnuMax){
        if (nu.momentum().e() < EnuMin) return false;
        if (nu.momentum().e() < EnuMax) return false;
    }
    
    int nMesons  = ps::sel::OutPartAny(ev, ps::pdg::kPiPlus).size();
    int nLeptons = ps::sel::OutPartAny(ev, ps::pdg::kPiPlus).size();
    int nPion    = ps::sel::OutPartAny(ev, ps::pdg::kPiPlus).size();

    // Check that the desired pion exists and is the only meson
    if (nPion != 1 || nMesons != 1) return false;

    // Check that there is only one final state lepton
    if (nLeptons != 1) return false;

    // Check that there are only three FS particles
    // if (event->NumFSParticle() != 3) return false; // Unsure how to do this.

    return true;
}


int ANL_167744_CC1npip_Filter_CC1pi3Prong(HepMC3::GenEvent const &ev) {
    return isCC1pi3Prong(ev, 14, 211, 2112, 0.0 * ps::GeV, 1.5 * ps::GeV);
}

int ANL_167744_CC1npip_Filter_CC1pi3Prong_HadMass() {
    if (!ANL_167744_CC1npip_Filter_CC1pi3Prong()) return 0;
    // return HadMassCut();
    return 1;
}


double ANL_167744_CC1npip_Project_cosmuStar(HepMC3::GenEvent const &ev) {
    auto nu = ps::sel::Beam(ev, ps::pdg::kNuMu);
    auto in = ps::sel::Target(ev, ps::pdg::kNeutron);
    auto mu = ps::sel::OutPartHM(ev, ps::pdg::kMuon);
    auto pip = ps::sel::OutPartHM(ev, ps::pdg::kPiPlus);
    auto n = ps::sel::OutPartHM(ev, ps::pdg::kNeutron);

    if (!nu || !in || !mu || !pip || !n) {
        return 0xdeadbeef;
    }

    double hadMass = HadronicMass_NucleonPi(n, pip);
    if (hadMass > 1400 * ps::MeV) {
        return 0xdeadbeef;
    }

    // Now need to boost into center-of-mass frame
    auto cms = (nu->momentum() + in->momentum());
    auto pmub = BoostVector(mu->momentum(), -cms);
    auto pnub = BoostVector(nu->momentum(), cms);

    return CosAngle(pmub, pnub);
}

}




// oid ANL_167744_CC1npip_Project_ppi(FitEvent *event) {

//   if (event->NumFSParticle(2112) == 0 ||
//       event->NumFSParticle(211) == 0 ||
//       event->NumFSParticle(13) == 0)
//     return;

//   TLorentzVector Pnu  = event->GetNeutrinoIn()->fP;
//   TLorentzVector Pn   = event->GetHMFSParticle(2112)->fP;
//   TLorentzVector Ppip = event->GetHMFSParticle(211)->fP;
//   TLorentzVector Pmu  = event->GetHMFSParticle(13)->fP;

//   double hadMass = FitUtils::MpPi(Pn, Ppip);
//   double ppip;

//   // This measurement has a 1.4 GeV M(Npi) constraint
//   if (hadMass < 1400) ppip = FitUtils::p(Ppip) * 1000.;
//   else ppip = -1.0;

//   fXVar = ppip;

//   return;
// }



// //********************************************************************
// void ANL_CC1npip_Evt_1DQ2_nu::FillEventVariables(FitEvent *event) {
// //********************************************************************

//   if (event->NumFSParticle(2112) == 0 || event->NumFSParticle(211) == 0 || event->NumFSParticle(13) == 0) {
//     return;
//   }

//   TLorentzVector Pnu  = event->GetNeutrinoIn()->fP;
//   TLorentzVector Pn   = event->GetHMFSParticle(2112)->fP;
//   TLorentzVector Ppip = event->GetHMFSParticle(211)->fP;
//   TLorentzVector Pmu  = event->GetHMFSParticle(13)->fP;

//   double hadMass = FitUtils::MpPi(Pn, Ppip);
//   double q2CCpip;

//   // ANL has a M(pi, p) < 1.4 GeV cut imposed (also no cut measurement but not useful for delta tuning)
//   if (hadMass < HadCut * 1000.) {
//     q2CCpip = -1.0 * (Pnu - Pmu).Mag2() / 1.E6;
//   } else {
//     q2CCpip = -1.0;
//   }

//   fXVar = q2CCpip;

//   return;
// };

// void ANL_CC1npip_Evt_1DWmupi_nu::FillEventVariables(FitEvent *event) {

//   if (event->NumFSParticle(2112) == 0 ||
//       event->NumFSParticle(211) == 0 ||
//       event->NumFSParticle(13) == 0)
//     return;

//   TLorentzVector Pmu  = event->GetHMFSParticle(13)->fP;
//   TLorentzVector Ppip = event->GetHMFSParticle(211)->fP;

//   double hadMass = (Pmu+Ppip).Mag()/1000.;

//   fXVar = hadMass;

//   return;
// };


// void ANL_CC1npip_Evt_1DWNmu_nu::FillEventVariables(FitEvent *event) {

//   if (event->NumFSParticle(2112) == 0 ||
//       event->NumFSParticle(211) == 0 ||
//       event->NumFSParticle(13) == 0)
//     return;

//   TLorentzVector Pn   = event->GetHMFSParticle(2112)->fP;
//   TLorentzVector Pmu  = event->GetHMFSParticle(13)->fP;

//   double hadMass = (Pn+Pmu).Mag()/1000.;

//   fXVar = hadMass;

//   return;
// };


// void ANL_CC1npip_Evt_1DWNpi_nu::FillEventVariables(FitEvent *event) {

//   if (event->NumFSParticle(2112) == 0 ||
//       event->NumFSParticle(211) == 0 ||
//       event->NumFSParticle(13) == 0)
//     return;

//   TLorentzVector Pn   = event->GetHMFSParticle(2112)->fP;
//   TLorentzVector Ppip = event->GetHMFSParticle(211)->fP;

//   double hadMass = (Pn+Ppip).Mag()/1000.;

//   fXVar = hadMass;

//   return;
// };



// int HadMassCut(){

// }

// int ANL_167744_CC1npip_Filter_CC1pi3Prong(){
//     return SignalDef::isCC1pi3Prong(event, 14, 211, 2112, 0.0*GeV, 1.5*GeV);
// }

// int ANL_167744_CC1npip_Filter_CC1pi3Prong_HadMass(){
//     if (! ANL_167744_CC1npip_Filter_CC1pi3Prong() ) return 0;
//     return HadMassCut();
// }


