#include "HepMC3/GenEvent.h"
#include "HepMC3/FourVector.h"


namespace ps {
    namespace sel {
        int NumOutPart(HepMC3::GenEvent const &ev, int pdg) { 
            return ps::sel::OutPartsAny(ev, {pdg}).size(); 
        }

        auto NeutrinoIn(HepMC3::GenEvent const &ev) {
            return ps::sel::BeamAny(ev, {-12,-14,-16,12,14,16});
        }
    }
}

// HELPERS
HepMC3::FourVector BoostVector(const HepMC3::FourVector& vecin, const HepMC3::FourVector& boost) {
    return vecin;
}

double CosAngle(const HepMC3::FourVector& v1, const HepMC3::FourVector &v2) {
    return 0.0;
}

double HadronicMass_NucleonPi(const HepMC3::FourVector&  Pn, const HepMC3::FourVector& Ppip) {
    return 0.0;
}

extern "C" {

bool isCC1pi3Prong(HepMC3::GenEvent const &ev,
    int nuPDG, int piPDG, int thirdPDG,
    double EnuMin, double EnuMax) {

    auto nu = ps::sel::Beam(ev, nuPDG);
    auto in = ps::sel::Target(ev);

    auto mu = ps::sel::OutPartHM(ev, nuPDG > 0 ? nuPDG-1 : nuPDG+1);
    auto pip = ps::sel::OutPartHM(ev, piPDG);
    auto n = ps::sel::OutPartHM(ev, thirdPDG);

    if (!nu || !in || !mu || !pip || !n) {
        return 0xdeadbeef;
    }

    if (EnuMin != EnuMax) {
        if (nu->momentum().e() < EnuMin) return false;
        if (nu->momentum().e() < EnuMax) return false;
    }

    int nMesons  = ps::sel::OutPartsAny(ev, {211, -211, 111}).size();

    int nLeptons = ps::sel::OutPartsAny(ev,
        ps::pdg::groups::kChargedLeptons).size();

    int nPion    = ps::sel::OutPartsAny(ev, {piPDG}).size();

    // Check that the desired pion exists and is the only meson
    if (nPion != 1 || nMesons != 1) return false;

    // Check that there is only one final state lepton
    if (nLeptons != 1) return false;

    // Check that there are only three FS particles
    // if (ps::sel::NumOutPart(ev,) != 3) return false; // Unsure how to do this.

    return true;
}


int ANL_167744_CC1npip_Filter_CC1pi3Prong(HepMC3::GenEvent const &ev) {
    return isCC1pi3Prong(ev, 14, 211, 2112, 0.0 * ps::GeV, 1.5 * ps::GeV);
}

int ANL_167744_CC1npip_Filter_CC1pi3Prong_HadMass(HepMC3::GenEvent const &ev) {
    if (!ANL_167744_CC1npip_Filter_CC1pi3Prong(ev)) return 0;
    // return HadMassCut();
    return 1;
}


double ANL_167744_CC1npip_Project_cosmuStar(HepMC3::GenEvent const &ev) {
    auto nu = ps::sel::Beam(ev, ps::pdg::kNuMu);
    auto in = ps::sel::Target(ev); // Need explicit initial state
    auto mu = ps::sel::OutPartHM(ev, ps::pdg::kMuon);
    auto pip = ps::sel::OutPartHM(ev, ps::pdg::kPiPlus);
    auto n = ps::sel::OutPartHM(ev, ps::pdg::kNeutron);

    if (!nu || !in || !mu || !pip || !n) {
        return 0xdeadbeef;
    }

    double hadMass = HadronicMass_NucleonPi(n->momentum(), pip->momentum());
    if (hadMass > 1400 * ps::MeV) {
        return 0xdeadbeef;
    }

    // Now need to boost into center-of-mass frame
    auto cms = (nu->momentum() + in->momentum());
    auto pmub = BoostVector(mu->momentum(), cms * -1.0);
    auto pnub = BoostVector(nu->momentum(), cms);

    return CosAngle(pmub, pnub);
}


double ANL_167744_CC1npip_Project_ppi(HepMC3::GenEvent const &ev) {

    if (ps::sel::NumOutPart(ev, 2112) == 0 ||
        ps::sel::NumOutPart(ev, 211) == 0 ||
        ps::sel::NumOutPart(ev, 13) == 0)
        return 0xdeadbeef;

    auto Pnu  = ps::sel::NeutrinoIn(ev)->momentum();
    auto Pn   = ps::sel::OutPartHM(ev, 2112)->momentum();
    auto Ppip = ps::sel::OutPartHM(ev, 211)->momentum();
    auto Pmu  = ps::sel::OutPartHM(ev, 13)->momentum();

    double hadMass = 0.0; //FitUtils::MpPi(Pn, Ppip);
    double mom_pip =  0.0;

    // This measurement has a 1.4 GeV M(Npi) constraint
    if (hadMass < 1.4 * ps::GeV)
        mom_pip = Ppip.length() * ps::GeV;
    else
        mom_pip = -1.0;

    return mom_pip;
}



double ANL_167744_CC1npip_Project_Q2(HepMC3::GenEvent const &ev) {

    if (ps::sel::NumOutPart(ev, 2112) == 0 ||
        ps::sel::NumOutPart(ev, 211) == 0 ||
        ps::sel::NumOutPart(ev, 13) == 0) {
        return 0xdeadbeef;
    }

    auto Pnu  = ps::sel::NeutrinoIn(ev)->momentum();
    auto Pn   = ps::sel::OutPartHM(ev, 2112)->momentum();
    auto Ppip = ps::sel::OutPartHM(ev, 211)->momentum();
    auto Pmu  = ps::sel::OutPartHM(ev, 13)->momentum();

    double q2CCpip = -1.0 * (Pnu - Pmu).length2() / 1.E6;

    return q2CCpip;
};

double ANL_167744_CC1npip_Project_Q2_LowW(HepMC3::GenEvent const &ev) {

    double hadMass = 0.0; //FitUtils::MpPi(Pn, Ppip);
    if (hadMass > 1.4*ps::GeV) {
        return 0xdeadbeef;
    }

    double q2CCpip = ANL_167744_CC1npip_Project_Q2(ev);
    return q2CCpip;
}


double ANL_167744_CC1npip_Project_Wmupi(HepMC3::GenEvent const &ev) {

    if (ps::sel::NumOutPart(ev, 2112) == 0 ||
        ps::sel::NumOutPart(ev, 211) == 0 ||
        ps::sel::NumOutPart(ev, 13) == 0)
        return 0xdeadbeef;

    auto Pmu  = ps::sel::OutPartHM(ev, 13)->momentum();
    auto Ppip = ps::sel::OutPartHM(ev, 211)->momentum();

    double hadMass = (Pmu+Ppip).length()/1000.;

    return hadMass;
}


double ANL_167744_CC1npip_Project_WNmu(HepMC3::GenEvent const &ev) {

    if (ps::sel::NumOutPart(ev, 2112) == 0 ||
        ps::sel::NumOutPart(ev, 211) == 0 ||
        ps::sel::NumOutPart(ev, 13) == 0)
        return 0xdeadbeef;

    auto Pn   = ps::sel::OutPartHM(ev, 2112)->momentum();
    auto Pmu  = ps::sel::OutPartHM(ev, 13)->momentum();

    return (Pn+Pmu).length()/1000.;
}


double ANL_167744_CC1npip_Project_Npi(HepMC3::GenEvent const &ev) {

    if (ps::sel::NumOutPart(ev, 2112) == 0 ||
        ps::sel::NumOutPart(ev, 211) == 0 ||
        ps::sel::NumOutPart(ev, 13) == 0)
        return 0xdeadbeef;

    auto Pn   = ps::sel::OutPartHM(ev, 2112)->momentum();
    auto Ppip = ps::sel::OutPartHM(ev, 211)->momentum();

    return (Pn+Ppip).length()/1000.;
}

}
