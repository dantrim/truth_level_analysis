#include "truth_level_analysis/hh_truth_validation.h"

//std/stl
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
using namespace std;

//analysis
#include "truth_level_analysis/helpers.h"

//ROOT
#include "TFile.h"

//xAOD

namespace truth
{

HHTruthValidation::HHTruthValidation() :
    m_mc_weight(1.0),
    m_suffix(""),
    m_outfilename("")
{
}

void HHTruthValidation::SlaveBegin(TTree* /*tree*/)
{
    timer()->Start();
}

void HHTruthValidation::Terminate()
{
    TruthSelectorBase::Terminate();
}

void HHTruthValidation::setup_output()
{
    stringstream ofn;
    ofn << "hh_validation";
    if(suffix()!="") ofn << "_" << suffix();
    ofn << ".root";
    cout << MYFUNC << " setting up output file (name = " << ofn.str() << ")" << endl;
    m_rfile = new TFile(ofn.str().c_str(), "RECREATE");
    if(m_rfile->IsZombie()) {
        cout << MYFUNC << " ERROR: could not create output file" << endl;
        exit(1);
    }
    m_rfile->cd();

}

Bool_t HHTruthValidation::Process(Long64_t entry)
{
    static Long64_t chain_entry = -1;
    chain_entry++;
    event()->getEntry(chain_entry);

    if(chain_entry == 0) {
        setup_output();
    }

    // let's just build up the collections here

    const xAOD::EventInfo* ei = 0;
    RETURN_CHECK(GetName(), event()->retrieve(ei, "EventInfo"));

    const xAOD::TruthParticleContainer* truthParticles = 0;
    RETURN_CHECK(GetName(), event()->retrieve(truthParticles, "TruthParticles"));

    const xAOD::TruthParticleContainer* truthElectrons = 0;
    RETURN_CHECK(GetName(), event()->retrieve(truthElectrons, "TruthElectrons"));

    const xAOD::TruthParticleContainer* truthMuons = 0;
    RETURN_CHECK(GetName(), event()->retrieve(truthMuons, "TruthMuons"));

    const xAOD::TruthVertexContainer* truthVertices = 0;
    RETURN_CHECK(GetName(), event()->retrieve(truthVertices, "TruthVertices"));

    return true;
}

} // namespace
