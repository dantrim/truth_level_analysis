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
#include "xAODEventInfo/EventInfo.h"

namespace truth
{

HHTruthValidation::HHTruthValidation() :
    m_mc_weight(1.0),
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
    save_histograms(m_rfile);
}

void HHTruthValidation::setup_output()
{
    stringstream ofn;
    // assume that the output directory has already been tested for existence
    ofn << output_dir() << "/" << "hh_validation";
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

void HHTruthValidation::setup_histograms()
{
    add_histogram(Histo("l0_pt", ";Lead lepton p_{T} [GeV]", 100, 0, -1));
}

Bool_t HHTruthValidation::Process(Long64_t entry)
{
    static Long64_t chain_entry = -1;
    chain_entry++;
    event()->getEntry(chain_entry);

    if(chain_entry == 0) {
        setup_output();
        setup_histograms();
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

    // set the MC weight
    mc_weight(ei->mcEventWeight());

    for(const auto & p : *truthElectrons) {
        fill_histo("l0_pt", p->pt() * mev2gev);
    }

    return true;
}

} // namespace
