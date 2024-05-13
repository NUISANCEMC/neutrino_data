#include "HepMC3/GenEvent.h"
#include "HepMC3/FourVector.h"

#define Dot(u,v)  (u.x()*v.x() + u.y()*v.y() + u.z()*v.z())

HepMC3::FourVector Cross(const HepMC3::FourVector& v1,
    const HepMC3::FourVector& v2) {
    return HepMC3::FourVector(
    v1.y()*v2.z()-v2.y()*v1.z(),
    v1.z()*v2.x()-v2.z()*v1.x(),
    v1.x()*v2.y()-v2.x()*v1.y(),
    0.0);
}

// HELPERS
HepMC3::FourVector BoostVector(
    const HepMC3::FourVector& v1,
    const HepMC3::FourVector& boost) {

    HepMC3::FourVector vo;

    // Boost this Lorentz vector
    double bx = boost.x();
    double by = boost.y();
    double bz = boost.z();

    double b2 = bx*bx + by*by + bz*bz;
    double gamma = 1.0 / sqrt(1.0 - b2);
    double bp = bx*v1.x() + by*v1.y() + bz*v1.z();
    double gamma2 = b2 > 0 ? (gamma - 1.0)/b2 : 0.0;
    
    vo.set_x(v1.x() + gamma2*bp*bx + gamma*bx*v1.e());
    vo.set_y(v1.y() + gamma2*bp*by + gamma*by*v1.e());
    vo.set_z(v1.z() + gamma2*bp*bz + gamma*bz*v1.e());
    vo.set_e(gamma*(v1.e() + bp));

    return vo;
}

double Angle(const HepMC3::FourVector& v1, const HepMC3::FourVector &q) {
    double ptot2 = v1.length2()*q.length2();
    if (ptot2 <= 0) {
        return 0.0;
    } else {
        double arg = Dot(v1, q)/sqrt(ptot2);
        if (arg >  1.0) arg =  1.0;
        if (arg < -1.0) arg = -1.0;
        return acos(arg);
    }
}

double CosAngle(const HepMC3::FourVector& v1, const HepMC3::FourVector &v2) {
    return cos(Angle(v1, v2));
}

double InvariantMass(const HepMC3::FourVector& pp, const HepMC3::FourVector& ppi) {

    double E_p = pp.e();
    double p_p = pp.length();
    double m_p = sqrt(E_p * E_p - p_p * p_p);

    double E_pi = ppi.e();
    double p_pi = ppi.length();
    double m_pi = sqrt(E_pi * E_pi - p_pi * p_pi);

    double th_p_pi = Angle(pp, ppi);

    // fairly easy thing to derive since bubble chambers measure the proton!
    double invMass = sqrt(m_p * m_p + m_pi * m_pi + 2 * E_p * E_pi -
                            2 * p_pi * p_p * cos(th_p_pi));

    return invMass;

};

double CosThAdler(const HepMC3::FourVector& Pnu, const HepMC3::FourVector& Pmu,
    const HepMC3::FourVector& Ppi, const HepMC3::FourVector& Pprot) {

    // Get the "resonance" lorentz vector (pion proton system)
    HepMC3::FourVector Pres = Pprot + Ppi;

    // Boost the particles into the resonance rest frame so we can define the
    // x,y,z axis
    auto Pnub = BoostVector(Pnu, Pres * -1.0);
    auto Pmub = BoostVector(Pmu, Pres * -1.0);
    auto Ppib = BoostVector(Ppi, Pres * -1.0);

    // The z-axis is defined as the axis of three-momentum transfer, \vec{k}
    // Also unit normalise the axis
    auto zAxis = (Pnu - Pmu);
    zAxis *= 1.0 / zAxis.length();

    // Then the angle between the pion in the "resonance" rest-frame and the
    // z-axis is the theta Adler angle
    double costhAdler = CosAngle(Ppi, zAxis);

    return costhAdler;
}

double PhiAdler(const HepMC3::FourVector& Pnu, const HepMC3::FourVector& Pmu,
    const HepMC3::FourVector& Ppi, const HepMC3::FourVector& Pprot) {

    // Get the "resonance" lorentz vector (pion proton system)
    HepMC3::FourVector Pres = Pprot + Ppi;

    // Boost the particles into the resonance rest frame so we can define the
    // x,y,z axis
    auto Pnub = BoostVector(Pnu, Pres * -1.0);
    auto Pmub = BoostVector(Pmu, Pres * -1.0);
    auto Ppib = BoostVector(Ppi, Pres * -1.0);

    // The z-axis is defined as the axis of three-momentum transfer, \vec{k}
    // Also unit normalise the axis
    auto zAxis = (Pnu - Pmu);
    zAxis *= 1.0 / zAxis.length();

    // The y-axis is then defined perpendicular to z and muon momentum in the
    // resonance frame
    auto yAxis = Cross(Pnu, Pmu);
    yAxis /= yAxis.length();

    // And the x-axis is then simply perpendicular to z and x
    auto xAxis = Cross(yAxis, zAxis);
    xAxis /= xAxis.length();

    // Project the pion on to x and y axes
    double x = Dot(Ppi, xAxis);
    double y = Dot(Ppi, yAxis);

    double newphi = atan2(y, x) * (180. / M_PI);
    // Convert negative angles to positive
    if (newphi < 0.0)
        newphi += 360.0;

    return newphi;
}

namespace ANL {
namespace inspire167744 {

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
        if (nu->momentum().e() > EnuMax) return false;
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
    if (ps::sel::OutPartsAny(ev,{}).size() != 3) return false; // Unsure how to do this.

    return true;
}

int HadMassCut(HepMC3::GenEvent const &ev,
    int pion_pdg,
    int nucleon_pdg,
    double hadronic_mass_cut) {

    if (!ps::sel::OutPartHM(ev, nucleon_pdg)) return 0;
    if (!ps::sel::OutPartHM(ev, pion_pdg)) return 0;

    auto pp = ps::sel::OutPartHM(ev, nucleon_pdg)->momentum();
    auto ppi = ps::sel::OutPartHM(ev, pion_pdg)->momentum();

    double mm = InvariantMass(pp, ppi);
    if (mm > hadronic_mass_cut) return 0;
    return 1;
}

}
}

using namespace ANL;
using namespace inspire167744;

extern "C" {
// ANL_167744_CC1npip Analysis
int ANL_167744_CC1npip_Selection(HepMC3::GenEvent const &ev) {
    return isCC1pi3Prong(ev, 14, 211, 2112, 0.0 * ps::GeV, 1.5 * ps::GeV);
}

int ANL_167744_CC1npip_Selection_lowW(HepMC3::GenEvent const &ev) {
    if (!ANL_167744_CC1npip_Selection(ev)) return 0;
    if (!HadMassCut(ev, 211, 2112, 1.4 * ps::GeV)) return 0;
    return 1;
}

double ANL_167744_CC1npip_Project_cosmuStar(HepMC3::GenEvent const &ev) {

    if (!ANL_167744_CC1npip_Selection_lowW(ev)) {
        return 0xdeadbeef;
    }

    auto nu = ps::sel::Beam(ev, ps::pdg::kNuMu);
    auto in = ps::sel::Target(ev); // Need explicit initial state
    auto mu = ps::sel::OutPartHM(ev, ps::pdg::kMuon);
    auto pip = ps::sel::OutPartHM(ev, ps::pdg::kPiPlus);
    auto n = ps::sel::OutPartHM(ev, ps::pdg::kNeutron);

    if (!nu || !in || !mu || !pip || !n) return 0xdeadbeef;

    // Now need to boost into center-of-mass frame
    auto cms = (nu->momentum() + in->momentum());
    auto pmub = BoostVector(mu->momentum(), cms * -1.0);
    auto pnub = BoostVector(nu->momentum(), cms);

    return CosAngle(pmub, pnub);
}


double ANL_167744_CC1npip_Project_ppi(HepMC3::GenEvent const &ev) {

    if (!ANL_167744_CC1npip_Selection_lowW(ev)) {
        return 0xdeadbeef;
    }

    if (!ps::sel::Beam(ev, ps::pdg::kNuMu)) return 0xdeadbeef;
    if (!ps::sel::OutPartHM(ev, 2112)) return 0xdeadbeef;
    if (!ps::sel::OutPartHM(ev, 211)) return 0xdeadbeef;
    if (!ps::sel::OutPartHM(ev, 13)) return 0xdeadbeef;

    auto Pnu  = ps::sel::Beam(ev, ps::pdg::kNuMu)->momentum();
    auto Pn   = ps::sel::OutPartHM(ev, 2112)->momentum();
    auto Ppip = ps::sel::OutPartHM(ev, 211)->momentum();
    auto Pmu  = ps::sel::OutPartHM(ev, 13)->momentum();

    return Ppip.length() * ps::MeV;
}



double ANL_167744_CC1npip_Project_Q2_highW(HepMC3::GenEvent const &ev) {

    if (!ANL_167744_CC1npip_Selection(ev)) {
        return 0xdeadbeef;
    }

    if (!ps::sel::Beam(ev, ps::pdg::kNuMu)) return 0xdeadbeef;
    if (!ps::sel::OutPartHM(ev, 2112)) return 0xdeadbeef;
    if (!ps::sel::OutPartHM(ev, 211)) return 0xdeadbeef;
    if (!ps::sel::OutPartHM(ev, 13)) return 0xdeadbeef;

    auto Pnu  = ps::sel::Beam(ev, ps::pdg::kNuMu)->momentum();
    auto Pn   = ps::sel::OutPartHM(ev, 2112)->momentum();
    auto Ppip = ps::sel::OutPartHM(ev, 211)->momentum();
    auto Pmu  = ps::sel::OutPartHM(ev, 13)->momentum();

    double q2CCpip = -(Pnu - Pmu).m2() / ps::GeV / ps::GeV;

    return q2CCpip;
};

double ANL_167744_CC1npip_Project_Q2_lowW(HepMC3::GenEvent const &ev) {
    if (!ANL_167744_CC1npip_Selection_lowW(ev)) {
        return 0xdeadbeef;
    }
    return ANL_167744_CC1npip_Project_Q2_highW(ev);
}


double ANL_167744_CC1npip_Project_Wmupi(HepMC3::GenEvent const &ev) {
    if (!ANL_167744_CC1npip_Selection(ev)) {
        return 0xdeadbeef;
    }

    if (!ps::sel::Beam(ev, ps::pdg::kNuMu)) return 0xdeadbeef;
    if (!ps::sel::OutPartHM(ev, 211)) return 0xdeadbeef;
    if (!ps::sel::OutPartHM(ev, 13)) return 0xdeadbeef;

    auto Pmu  = ps::sel::OutPartHM(ev, 13)->momentum();
    auto Ppip = ps::sel::OutPartHM(ev, 211)->momentum();

    return (Pmu+Ppip).m()/ps::GeV;
}


double ANL_167744_CC1npip_Project_WNmu(HepMC3::GenEvent const &ev) {
    if (!ANL_167744_CC1npip_Selection(ev)) {
        return 0xdeadbeef;
    }

    if (!ps::sel::Beam(ev, ps::pdg::kNuMu)) return 0xdeadbeef;
    if (!ps::sel::OutPartHM(ev, 2112)) return 0xdeadbeef;
    if (!ps::sel::OutPartHM(ev, 13)) return 0xdeadbeef;

    auto Pn   = ps::sel::OutPartHM(ev, 2112)->momentum();
    auto Pmu  = ps::sel::OutPartHM(ev, 13)->momentum();

    return (Pn+Pmu).m()/ps::GeV;
}


double ANL_167744_CC1npip_Project_WNpi(HepMC3::GenEvent const &ev) {
    if (!ANL_167744_CC1npip_Selection(ev)) {
        return 0xdeadbeef;
    }

    if (!ps::sel::Beam(ev, ps::pdg::kNuMu)) return 0xdeadbeef;
    if (!ps::sel::OutPartHM(ev, 2112)) return 0xdeadbeef;
    if (!ps::sel::OutPartHM(ev, 211)) return 0xdeadbeef;

    auto Pn   = ps::sel::OutPartHM(ev, 2112)->momentum();
    auto Ppip = ps::sel::OutPartHM(ev, 211)->momentum();

    return (Pn+Ppip).m()/ps::GeV;
}

// ANL_167744_CC1npi0 Analyses
int ANL_167744_CC1npi0_Selection(HepMC3::GenEvent const &ev) {
    return isCC1pi3Prong(ev, 14, 111, 2212, 0.0 * ps::GeV, 1.5 * ps::GeV);
}

int ANL_167744_CC1npi0_Selection_lowW(HepMC3::GenEvent const &ev) {
    if (!ANL_167744_CC1npi0_Selection(ev)) return 0;
    if (!HadMassCut(ev, 111, 2212, 1.4 * ps::GeV)) return 0;
    return 1;
}

double ANL_167744_CC1npi0_Project_cosmuStar(HepMC3::GenEvent const &ev) {
    if (!ANL_167744_CC1npi0_Selection(ev)) {
        return 0xdeadbeef;
    }

    if (!ps::sel::Beam(ev, ps::pdg::kNuMu) || !ps::sel::Target(ev)) return 0xdeadbeef;
    if (!ps::sel::OutPartHM(ev, 2212)) return 0xdeadbeef;
    if (!ps::sel::OutPartHM(ev, 111)) return 0xdeadbeef;
    if (!ps::sel::OutPartHM(ev, 13)) return 0xdeadbeef;


    auto Pnu  = ps::sel::Beam(ev, ps::pdg::kNuMu)->momentum();
    auto Pin  = ps::sel::Target(ev)->momentum();
    auto Pp   = ps::sel::OutPartHM(ev, 2212)->momentum();
    auto Ppi0 = ps::sel::OutPartHM(ev, 111)->momentum();
    auto Pmu  = ps::sel::OutPartHM(ev, 13)->momentum();

    double cosmuStar = -999;

    // Now need to boost into center-of-mass frame
    auto CMS = Pnu + Pin;
    Pmu = BoostVector(Pmu, CMS * -1.0);
    Pnu = BoostVector(Pnu, CMS * +1.0);

    return CosAngle(Pmu, Pnu);
}

double ANL_167744_CC1npi0_Project_Q2(HepMC3::GenEvent const &ev) {
    if (!ANL_167744_CC1npi0_Selection(ev)) {
        return 0xdeadbeef;
    }

    if (!ps::sel::Beam(ev, ps::pdg::kNuMu)) return 0xdeadbeef;
    if (!ps::sel::OutPartHM(ev, 2212)) return 0xdeadbeef;
    if (!ps::sel::OutPartHM(ev, 111)) return 0xdeadbeef;
    if (!ps::sel::OutPartHM(ev, 13)) return 0xdeadbeef;


    auto Pnu  = ps::sel::Beam(ev, ps::pdg::kNuMu)->momentum();
    auto Pp   = ps::sel::OutPartHM(ev, 2212)->momentum();
    auto Ppi0 = ps::sel::OutPartHM(ev, 111)->momentum();
    auto Pmu  = ps::sel::OutPartHM(ev, 13)->momentum();

    return -(Pnu - Pmu).m2() / ps::GeV / ps::GeV;
};

double ANL_167744_CC1npi0_Project_Q2_lowW(HepMC3::GenEvent const &ev) {
    if (!ANL_167744_CC1npi0_Selection_lowW(ev)) {
        return 0xdeadbeef;
    }
    return ANL_167744_CC1npi0_Project_Q2(ev);
};

double ANL_167744_CC1npi0_Project_Wmupi(HepMC3::GenEvent const &ev) {
    if (!ANL_167744_CC1npi0_Selection(ev)) {
        return 0xdeadbeef;
    }

    if (!ps::sel::OutPartHM(ev, 111)) return 0xdeadbeef;
    if (!ps::sel::OutPartHM(ev, 13)) return 0xdeadbeef;

    auto Pmu  = ps::sel::OutPartHM(ev, 13)->momentum();
    auto Ppip = ps::sel::OutPartHM(ev, 111)->momentum();
    return (Pmu+Ppip).m()/ps::GeV;
};

double ANL_167744_CC1npi0_Project_WNpi(HepMC3::GenEvent const &ev) {
    if (!ANL_167744_CC1npi0_Selection(ev)) {
        return 0xdeadbeef;
    }

    if (!ps::sel::OutPartHM(ev, 111)) return 0xdeadbeef;
    if (!ps::sel::OutPartHM(ev, 2212)) return 0xdeadbeef;

    auto Pp   = ps::sel::OutPartHM(ev, 2212)->momentum();
    auto Ppip = ps::sel::OutPartHM(ev, 111)->momentum();
    return (Pp+Ppip).m()/ps::GeV;
};

double ANL_167744_CC1npi0_Project_WNmu(HepMC3::GenEvent const &ev) {
    if (!ANL_167744_CC1npi0_Selection(ev)) {
        return 0xdeadbeef;
    }

    if (!ps::sel::OutPartHM(ev, 2212)) return 0xdeadbeef;
    if (!ps::sel::OutPartHM(ev, 13)) return 0xdeadbeef;

    auto Pp   = ps::sel::OutPartHM(ev, 2212)->momentum();
    auto Ppip = ps::sel::OutPartHM(ev, 13)->momentum();
    return (Pp+Ppip).m()/ps::GeV;
};


int ANL_167744_CC1ppip_Selection(HepMC3::GenEvent const &ev) {
    if (!HadMassCut(ev, 211, 2212, 10.0 * ps::GeV)) return 0;
    return isCC1pi3Prong(ev, 14, 211, 2212, 0.0 * ps::GeV, 6.0 * ps::GeV);
}

int ANL_167744_CC1ppip_Selection_lowW(HepMC3::GenEvent const &ev) {
    if (!ANL_167744_CC1ppip_Selection(ev)) return 0;
    if (!HadMassCut(ev, 211, 2212, 1.4 * ps::GeV)) return 0;
    return 1;
}

double ANL_167744_CC1ppip_Project_cosmuStar(HepMC3::GenEvent const &ev) {
    if (!ANL_167744_CC1ppip_Selection_lowW(ev)) {
        return 0xdeadbeef;
    }

    if (!ps::sel::Beam(ev, ps::pdg::kNuMu)) return 0xdeadbeef;
    if (!ps::sel::OutPartHM(ev, 2212)) return 0xdeadbeef;
    if (!ps::sel::OutPartHM(ev, 211)) return 0xdeadbeef;
    if (!ps::sel::OutPartHM(ev, 13)) return 0xdeadbeef;

    auto Pnu  = ps::sel::Beam(ev, ps::pdg::kNuMu)->momentum();
    auto Pin  = ps::sel::Target(ev)->momentum();
    auto Pp   = ps::sel::OutPartHM(ev, 2212)->momentum();
    auto Ppip = ps::sel::OutPartHM(ev, 211)->momentum();
    auto Pmu  = ps::sel::OutPartHM(ev, 13)->momentum();

    // Now need to boost into center-of-mass frame
    auto cms = (Pnu + Pin);
    auto pmub = BoostVector(Pmu, cms * -1.0);
    auto pnub = BoostVector(Pnu, cms);

    return CosAngle(pmub, pnub);
};


double ANL_167744_CC1ppip_Project_costhAdler(HepMC3::GenEvent const &ev) {
    if (!ANL_167744_CC1ppip_Selection_lowW(ev)) {
        return 0xdeadbeef;
    }

    if (!ps::sel::Beam(ev, ps::pdg::kNuMu)) return 0xdeadbeef;
    if (!ps::sel::OutPartHM(ev, 2212)) return 0xdeadbeef;
    if (!ps::sel::OutPartHM(ev, 211)) return 0xdeadbeef;
    if (!ps::sel::OutPartHM(ev, 13)) return 0xdeadbeef;

    auto Pnu  = ps::sel::Beam(ev, ps::pdg::kNuMu)->momentum();
    auto Pp   = ps::sel::OutPartHM(ev, ps::pdg::kProton)->momentum();
    auto Ppip = ps::sel::OutPartHM(ev, ps::pdg::kPiPlus)->momentum();
    auto Pmu  = ps::sel::OutPartHM(ev, ps::pdg::kMuon)->momentum();

    // Get Adler cos theta
    return CosThAdler(Pnu, Pmu, Ppip, Pp);
};


double ANL_167744_CC1ppip_Project_phi(HepMC3::GenEvent const &ev) {
    if (!ANL_167744_CC1ppip_Selection_lowW(ev)) {
        return 0xdeadbeef;
    }

    if (!ps::sel::Beam(ev, ps::pdg::kNuMu)) return 0xdeadbeef;
    if (!ps::sel::OutPartHM(ev, 2212)) return 0xdeadbeef;
    if (!ps::sel::OutPartHM(ev, 211)) return 0xdeadbeef;
    if (!ps::sel::OutPartHM(ev, 13)) return 0xdeadbeef;

    auto Pnu  = ps::sel::Beam(ev, ps::pdg::kNuMu)->momentum();
    auto Pp   = ps::sel::OutPartHM(ev, ps::pdg::kProton)->momentum();
    auto Ppip = ps::sel::OutPartHM(ev, ps::pdg::kPiPlus)->momentum();
    auto Pmu  = ps::sel::OutPartHM(ev, ps::pdg::kMuon)->momentum();

    return PhiAdler(Pnu, Pmu, Ppip, Pp);
};


double ANL_167744_CC1ppip_Project_Q2(HepMC3::GenEvent const &ev) {
    if (!ANL_167744_CC1ppip_Selection(ev)) {
        return 0xdeadbeef;
    }

    if (!ps::sel::Beam(ev, ps::pdg::kNuMu)) return 0xdeadbeef;
    if (!ps::sel::OutPartHM(ev, 2212)) return 0xdeadbeef;
    if (!ps::sel::OutPartHM(ev, 211)) return 0xdeadbeef;
    if (!ps::sel::OutPartHM(ev, 13)) return 0xdeadbeef;

    auto Pnu  = ps::sel::Beam(ev, ps::pdg::kNuMu)->momentum();
    auto Pp   = ps::sel::OutPartHM(ev, ps::pdg::kProton)->momentum();
    auto Ppip = ps::sel::OutPartHM(ev, ps::pdg::kPiPlus)->momentum();
    auto Pmu  = ps::sel::OutPartHM(ev, ps::pdg::kMuon)->momentum();

    return -(Pnu - Pmu).m2() / ps::GeV / ps::GeV;
};

double ANL_167744_CC1ppip_Project_Q2_lowW(HepMC3::GenEvent const &ev) {
    if (!ANL_167744_CC1ppip_Selection_lowW(ev)) {
        return 0xdeadbeef;
    }
    return ANL_167744_CC1ppip_Project_Q2(ev);
};


double ANL_167744_CC1ppip_Project_ppi(HepMC3::GenEvent const &ev) {
    if (!ANL_167744_CC1ppip_Selection_lowW(ev)) {
        return 0xdeadbeef;
    }

    if (!ps::sel::Beam(ev, ps::pdg::kNuMu)) return 0xdeadbeef;
    if (!ps::sel::OutPartHM(ev, ps::pdg::kProton)) return 0xdeadbeef;
    if (!ps::sel::OutPartHM(ev, ps::pdg::kPiPlus)) return 0xdeadbeef;
    if (!ps::sel::OutPartHM(ev, ps::pdg::kMuon)) return 0xdeadbeef;


    auto Pnu  = ps::sel::Beam(ev, ps::pdg::kNuMu)->momentum();
    auto Pp   = ps::sel::OutPartHM(ev, ps::pdg::kProton)->momentum();
    auto Ppip = ps::sel::OutPartHM(ev, ps::pdg::kPiPlus)->momentum();
    auto Pmu  = ps::sel::OutPartHM(ev, ps::pdg::kMuon)->momentum();

    return Ppip.length() / ps::MeV;
}

double ANL_167744_CC1ppip_Project_thpr(HepMC3::GenEvent const &ev) {
    if (!ANL_167744_CC1ppip_Selection_lowW(ev)) {
        return 0xdeadbeef;
    }

    if (!ps::sel::Beam(ev, ps::pdg::kNuMu)) return 0xdeadbeef;
    if (!ps::sel::OutPartHM(ev, ps::pdg::kProton)) return 0xdeadbeef;
    if (!ps::sel::OutPartHM(ev, ps::pdg::kPiPlus)) return 0xdeadbeef;
    if (!ps::sel::OutPartHM(ev, ps::pdg::kMuon)) return 0xdeadbeef;

    auto Pnu  = ps::sel::Beam(ev, ps::pdg::kNuMu)->momentum();
    auto Pp   = ps::sel::OutPartHM(ev, ps::pdg::kProton)->momentum();
    auto Ppip = ps::sel::OutPartHM(ev, ps::pdg::kPiPlus)->momentum();
    auto Pmu  = ps::sel::OutPartHM(ev, ps::pdg::kMuon)->momentum();

    return CosAngle(Pnu, Pp);
};



double ANL_167744_CC1ppip_Project_Wmupi(HepMC3::GenEvent const &ev) {
    if (!ANL_167744_CC1ppip_Selection(ev)) {
        return 0xdeadbeef;
    }

    if (!ps::sel::Beam(ev, ps::pdg::kNuMu)) return 0xdeadbeef;
    if (!ps::sel::OutPartHM(ev, ps::pdg::kProton)) return 0xdeadbeef;
    if (!ps::sel::OutPartHM(ev, ps::pdg::kPiPlus)) return 0xdeadbeef;
    if (!ps::sel::OutPartHM(ev, ps::pdg::kMuon)) return 0xdeadbeef;

    auto Pmu  = ps::sel::OutPartHM(ev, ps::pdg::kMuon)->momentum();
    auto Ppip = ps::sel::OutPartHM(ev, ps::pdg::kPiPlus)->momentum();
    return (Pmu+Ppip).m()/ps::GeV;
};


double ANL_167744_CC1ppip_Project_WNmu(HepMC3::GenEvent const &ev) {
    if (!ANL_167744_CC1ppip_Selection(ev)) {
        return 0xdeadbeef;
    }

    if (!ps::sel::Beam(ev, ps::pdg::kNuMu)) return 0xdeadbeef;
    if (!ps::sel::OutPartHM(ev, ps::pdg::kProton)) return 0xdeadbeef;
    if (!ps::sel::OutPartHM(ev, ps::pdg::kPiPlus)) return 0xdeadbeef;
    if (!ps::sel::OutPartHM(ev, ps::pdg::kMuon)) return 0xdeadbeef;

    auto Pmu  = ps::sel::OutPartHM(ev, ps::pdg::kMuon)->momentum();
    auto Pp   = ps::sel::OutPartHM(ev, ps::pdg::kProton)->momentum();
    return (Pmu+Pp).m()/ps::GeV;
};

double ANL_167744_CC1ppip_Project_WNpi(HepMC3::GenEvent const &ev) {
    if (!ANL_167744_CC1ppip_Selection(ev)) {
        return 0xdeadbeef;
    }

    if (!ps::sel::Beam(ev, ps::pdg::kNuMu)) return 0xdeadbeef;
    if (!ps::sel::OutPartHM(ev, ps::pdg::kProton)) return 0xdeadbeef;
    if (!ps::sel::OutPartHM(ev, ps::pdg::kPiPlus)) return 0xdeadbeef;
    if (!ps::sel::OutPartHM(ev, ps::pdg::kMuon)) return 0xdeadbeef;

    auto Pp   = ps::sel::OutPartHM(ev, ps::pdg::kProton)->momentum();
    auto Ppip = ps::sel::OutPartHM(ev, ps::pdg::kPiPlus)->momentum();
    return (Pp+Ppip).m()/ps::GeV;
};


double ANL_167744_CC1ppip_XSec_1DQ2(HepMC3::GenEvent const &ev) {
    if (!ANL_167744_CC1ppip_Selection_lowW(ev)) {
        return 0xdeadbeef;
    }

    if (!ps::sel::Beam(ev, ps::pdg::kNuMu)) return 0xdeadbeef;
    if (!ps::sel::OutPartHM(ev, ps::pdg::kProton)) return 0xdeadbeef;
    if (!ps::sel::OutPartHM(ev, ps::pdg::kPiPlus)) return 0xdeadbeef;
    if (!ps::sel::OutPartHM(ev, ps::pdg::kMuon)) return 0xdeadbeef;

    auto Pnu  = ps::sel::Beam(ev, ps::pdg::kNuMu)->momentum();
    auto Pp   = ps::sel::OutPartHM(ev, ps::pdg::kProton)->momentum();
    auto Ppip = ps::sel::OutPartHM(ev, ps::pdg::kPiPlus)->momentum();
    auto Pmu  = ps::sel::OutPartHM(ev, ps::pdg::kMuon)->momentum();

    return  -(Pnu - Pmu).m2() / ps::GeV / ps::GeV;
};

}
