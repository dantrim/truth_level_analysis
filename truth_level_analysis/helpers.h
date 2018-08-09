#ifndef TRUTHHELPERS_H
#define TRUTHHELPERS_H

//std/stl
#include <string>

//xAOD
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODJet/Jet.h"

//ROOT
#include "TLorentzVector.h"


namespace truth
{

static std::string computeMethodName(const std::string& function, const std::string& prettyFunction) 
{
    size_t locFunName = prettyFunction.find(function);
    size_t begin = prettyFunction.rfind(" ", locFunName) + 1;
    size_t end = prettyFunction.find("(",locFunName + function.length());
    if(prettyFunction[end+1] == ')') {
        return (prettyFunction.substr(begin, end - begin) + "()");
    }
    else {
        return (prettyFunction.substr(begin, end - begin) + "(...)");
    }
}
#define MYFUNC computeMethodName(__FUNCTION__,__PRETTY_FUNCTION__).c_str() // c_str() optional
//constants
const double mev2gev = 1.0e-3;

struct IsPdg {
    int pdg;
    IsPdg(int p) : pdg(p) {}
    bool operator() (const xAOD::TruthParticle* l) { return l && l->absPdgId()==pdg; }
};

static IsPdg IsElectron(11);
static IsPdg IsMuon(13);
static IsPdg IsTau(15);

static bool IsEE( const xAOD::TruthParticle* l0, const xAOD::TruthParticle* l1 )
{
    return ( IsElectron(l0) && IsElectron(l1) );
}
static bool IsMM( const xAOD::TruthParticle* l0, const xAOD::TruthParticle* l1 )
{
    return ( IsMuon(l0) && IsMuon(l1) );
}
static bool IsEM( const xAOD::TruthParticle* l0, const xAOD::TruthParticle* l1 )
{
    return ( IsElectron(l0) && IsMuon(l1) );
}
static bool IsME( const xAOD::TruthParticle* l0, const xAOD::TruthParticle* l1 )
{
    return ( IsMuon(l0) && IsElectron(l1) );
}
static bool IsSF( const xAOD::TruthParticle* l0, const xAOD::TruthParticle* l1 )
{
    return ( IsEE(l0,l1) || IsMM(l0,l1) );
}
static bool IsDF( const xAOD::TruthParticle* l0, const xAOD::TruthParticle* l1 )
{
    return ( IsEM(l0,l1) || IsME(l0,l1) );
}

struct pt_greater {
    bool operator() (const xAOD::TruthParticle* a, const xAOD::TruthParticle* b) { return a->pt() > b->pt(); }
};
static pt_greater PtGreater;

struct pt_greaterTLV {
    bool operator() (const TLorentzVector& a, const TLorentzVector& b) { return a.Pt() > b.Pt(); }
};
static pt_greaterTLV PtGreaterTLV;

struct pt_greaterJet {
    bool operator() (const xAOD::Jet* a, const xAOD::Jet* b) { return a->pt() > b->pt(); }
};
static pt_greaterJet PtGreaterJet;

struct Histo
{
    Histo() :
        name(""),
        title(""),
        nbinsX(0),
        lowX(0),
        highX(0),
        nbinsY(0),
        lowY(0),
        highY(0) {}
    Histo(std::string name = "", std::string title = "",
        int nbinsX = 0, double lowX = 0, double highX = 0,
        int nbinsY = 0, double lowY = 0, double highY = 0) :
        name(name),
        title(title),
        nbinsX(nbinsX),
        lowX(lowX),
        highX(highX),
        nbinsY(nbinsY),
        lowY(lowY),
        highY(highY) {}

    std::string name;
    std::string title;
    int nbinsX;
    double lowX;
    double highX;
    int nbinsY;
    double lowY;
    double highY;

    bool is2D() { return nbinsY != 0; }
};

enum DileptonFlavor {
    EE=0,
    MM,
    TT,
    EM,
    ME,
    TE,
    ET,
    TM,
    MT,
    SF,
    Invalid
};

static DileptonFlavor get_dilepton_flavor(const xAOD::TruthParticle* l0, const xAOD::TruthParticle* l1) {

    bool e0 = l0->isElectron();
    bool m0 = l0->isMuon();
    bool t0 = l0->isTau();
    bool e1 = l1->isElectron();
    bool m1 = l1->isMuon();
    bool t1 = l1->isTau();

    if(e0 && e1) return DileptonFlavor::EE;
    else if(e0 && m1) return DileptonFlavor::EM;
    else if(e0 && t1) return DileptonFlavor::ET;
    else if(m0 && m1) return DileptonFlavor::MM;
    else if(m0 && e1) return DileptonFlavor::ME;
    else if(m0 && t1) return DileptonFlavor::MT;
    else if(t0 && t1) return DileptonFlavor::TT;
    else if(t0 && e1) return DileptonFlavor::TE;
    else if(t0 && m1) return DileptonFlavor::TM;
    else
        return DileptonFlavor::Invalid;
}


} // namespace

#endif
