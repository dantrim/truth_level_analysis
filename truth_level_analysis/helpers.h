#ifndef TRUTHHELPERS_H
#define TRUTHHELPERS_H

//std/stl
#include <string>

//xAOD
#include "xAODTruth/TruthParticleContainer.h"

namespace xAOD {
    class TEvent;
    class TStore;
}

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


} // namespace

#endif