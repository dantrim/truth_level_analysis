#include "truth_level_analysis/truth_ntupler.h"

//std/stl
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
using namespace std;

//analysis
#include "truth_level_analysis/helpers.h"

//xAOD
#include "xAODEventInfo/EventInfo.h"
#include "xAODJet/Jet.h"
#include "xAODJet/JetContainer.h"

//ROOT
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"

#define FBRANCH(name) \
    do { \
        b_##name = output_tree()->Branch(#name, &m_float_vm[#name]); \
    } while(0);

#define IBRANCH(name) \
    do { \
        b_##name = output_tree()->Branch(#name, &m_int_vm[#name]); \
    } while(0);

namespace truth
{

TruthNtupler::TruthNtupler() :
    m_do_sumw(false),
    m_dsid(0),
    m_xsec(-1.0),
    m_sumw(-1.0)
{
}

void TruthNtupler::SlaveBegin(TTree* /*tree*/)
{
    timer()->Start();
}

void TruthNtupler::Terminate()
{
    if(do_sumw()) {
        cout << MYFUNC << " ===================================" << endl;
        cout << MYFUNC << "  Total sumw: " << m_sumw << endl;
        cout << MYFUNC << " ===================================" << endl;
    }
    TruthSelectorBase::Terminate();

    if(!do_sumw()) {
        output_file()->cd();
        output_tree()->Write();
    }
}

void TruthNtupler::setup_output()
{
    cout << MYFUNC << " setting up output ntupler structure" << endl;

    // setup the TTree and TFile
    m_output_tree_file = new TFile("test.root", "RECREATE");
    m_output_tree = new TTree("truth", "truth");

    FBRANCH(l0_pt);
    IBRANCH(nBJets);
}

void TruthNtupler::reset_ttree_vars()
{
    for(auto m : m_float_vm) {
        m_float_vm.at(m.first) = -99.;
    }
    for(auto m : m_int_vm) {
        m_int_vm.at(m.first) = -1;
    }
}

void TruthNtupler::var(string key, int val)
{
    try {
        m_int_vm.at(key) = val;
    }
    catch(std::exception& e)
    {
        cout << MYFUNC << " could not fill variable \"" << key << "\": " << e.what() << endl;
        exit(1);
    }
}

void TruthNtupler::var(string key, float val)
{
    try {
        m_float_vm.at(key) = val;
    }
    catch(std::exception& e)
    {
        cout << MYFUNC << " could not fill variable \"" << key << "\": " << e.what() << endl;
        exit(1);
    }
}

Bool_t TruthNtupler::Process(Long64_t entry)
{
    static Long64_t chain_entry = -1;
    chain_entry++;
    event()->getEntry(chain_entry);

    if(chain_entry == 0) {
        if(!do_sumw()) {
            setup_output();
        }
    }

    reset_ttree_vars();

    if(chain_entry % 1000 == 0) {
        cout << MYFUNC << " Processing entry " << chain_entry << endl;
    }

    if(do_sumw()) {
        update_sumw();
    }
    else {
        float l0_pt = 33.2;
        var("l0_pt", l0_pt);
        int nBJets = 3;
        var("nBJets", nBJets);
        output_tree()->Fill();
    }

    return true;
}

bool TruthNtupler::update_sumw()
{
    const xAOD::EventInfo* ei = 0;
    RETURN_CHECK(GetName(), event()->retrieve(ei, "EventInfo"));
    m_sumw += ei->mcEventWeight();

    return true;
}


} // namespace truth

