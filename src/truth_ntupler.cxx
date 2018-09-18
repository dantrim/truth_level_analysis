#include "truth_level_analysis/truth_ntupler.h"

//std/stl
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
using namespace std;

//analysis
#include "truth_level_analysis/helpers.h"
#include "truth_level_analysis/utility/string_utils.h"

//xAOD
#include "xAODEventInfo/EventInfo.h"
#include "xAODJet/Jet.h"
#include "xAODJet/JetContainer.h"

//ROOT
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"

//xAOD
#include "PathResolver/PathResolver.h"

#define FBRANCH(name) \
    do { \
        output_tree()->Branch(#name, &m_float_vm[#name]); \
    } while(0);

#define IBRANCH(name) \
    do { \
        output_tree()->Branch(#name, &m_int_vm[#name]); \
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

    IBRANCH(dsid)
    FBRANCH(eventweight)
    FBRANCH(w)
    FBRANCH(sumw)

    // leptons
    IBRANCH(l_flav)
    FBRANCH(l0_pt)
    FBRANCH(l1_pt)
    FBRANCH(l0_eta)
    FBRANCH(l1_eta)
    FBRANCH(l0_phi)
    FBRANCH(l1_phi)
    IBRANCH(isEE)
    IBRANCH(isMM)
    IBRANCH(isSF)
    IBRANCH(isDF)
    
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

void TruthNtupler::load_sumw_and_xsec()
{
    string xsec_dir_search = "truth_level_analysis/cross_section/";
    string sumw_dir_search = "truth_level_analysis/sumw/";

    string xsec_dir = PathResolverFindCalibDirectory(xsec_dir_search);
    string sumw_dir = PathResolverFindCalibDirectory(sumw_dir_search);

    if(xsec_dir == "") {
        cout << MYFUNC << " ERROR: could not locate xsec directory (=" << xsec_dir_search<<")" << endl;
        exit(1);
    }
    if(sumw_dir == "") {
        cout << MYFUNC << " ERROR: could not locate sumw directory (=" << sumw_dir_search<<")" << endl;
        exit(1);
    }

    cout << MYFUNC << " Found xsec directory: " << xsec_dir << endl;
    cout << MYFUNC << " Found sumw directory: " << sumw_dir << endl;

    // now load the files
    string xsec_filename = xsec_dir + "/xsec_data.txt";
    string sumw_filename = sumw_dir + "/sumw_data.txt";

    std::ifstream xsec_file(xsec_filename.c_str());
    if(!xsec_file.good()) {
        cout << MYFUNC << " ERROR: could not open or find xsec file (=" << xsec_filename << ")" << endl;
        exit(1);
    }
    std::ifstream sumw_file(sumw_filename.c_str());
    if(!sumw_file.good()) {
        cout << MYFUNC << " ERROR: could not open or find sumw file (=" << sumw_filename << ")" << endl;
        exit(1);
    }

    m_sumw_map.clear();
    m_xsec_map.clear();

    // read in xsec data
    string xsec_line;
    bool xsec_ok = true;
    while(std::getline(xsec_file, xsec_line)) {
        if(Susy::utils::startswith(xsec_line, "#")) continue;
        std::vector< std::string > tokens = Susy::utils::tokenizeString(xsec_line, ' ');
        int dsid = -1;
        double xsec = -1;
        try {
            dsid = std::stoi(tokens.at(0));
            xsec = std::stof(tokens.at(1));
        }
        catch(std::exception& e) {
            cout << MYFUNC << " ERROR: unable to interpret xsec data (xsec data = " << xsec_line << ")" << endl;
            xsec_ok = false;
        }

        // update the map
        m_xsec_map[dsid] = xsec;
    }
    if(!xsec_ok) exit(1);

    // read in sumw data
    string sumw_line;
    bool sumw_ok = true;
    while(std::getline(sumw_file, sumw_line)) {
        if(Susy::utils::startswith(sumw_line, "#")) continue;
        std::vector< std::string > tokens = Susy::utils::tokenizeString(sumw_line, ' ');
        int dsid = -1;
        double sumw = -1;
        try {
            dsid = std::stoi(tokens.at(0));
            sumw = std::stof(tokens.at(1));
        }
        catch(std::exception& e) {
            cout << MYFUNC << " ERROR: unable to interpret the sumw data (sumw data = " << sumw_line << ")" << endl;
            sumw_ok = false;
        }

        // update the map
        m_sumw_map[dsid] = sumw;
    }
    if(!sumw_ok) exit(1);

    bool data_ok = true;
    if(!(m_sumw_map.count(dsid()) > 0)) {
        cout << MYFUNC << " ERROR: sumw map does not have data for current DSID (=" << dsid() << ")" << endl;
        data_ok = false;
    }
    if(!(m_xsec_map.count(dsid()) > 0)) {
        cout << MYFUNC << " ERROR: xsec map does not have data for current DSID (=" << dsid() << ")" << endl;
        data_ok = false;
    }
    if(!data_ok) exit(1);

    cout << MYFUNC << " Setting sample cross-section to: " << m_xsec_map[dsid()] << " pb" << endl;
    cout << MYFUNC << " Setting sample sumw to         : " << m_sumw_map[dsid()] << endl;

    m_xsec = m_xsec_map[dsid()];
    m_sumw = m_sumw_map[dsid()];

}

Bool_t TruthNtupler::Process(Long64_t entry)
{
    static Long64_t chain_entry = -1;
    chain_entry++;
    event()->getEntry(chain_entry);

    if(chain_entry == 0) {
        if(!do_sumw()) {
            load_sumw_and_xsec();
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
        process_event();
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

bool TruthNtupler::process_event()
{
    // get the jets
    vector<const xAOD::Jet*> jets;
    const xAOD::JetContainer* xjets = 0;
    string jetkey = "AntiKt4TruthJets";
    if(!event()->contains<xAOD::JetContainer>(jetkey)) jetkey = "AntiKt4TruthDressedWZJets";
    if(!event()->retrieve(xjets, jetkey).isSuccess()) {
        cout << MYFUNC << " ERROR: could not retrieve truth jets with key " << jetkey << endl;
        exit(1);
    }
    for(const auto & jet : *xjets) {
        if(!(jet->pt() > 20. * mev2gev)) continue;
        jets.push_back(jet);
    }
    std::sort(jets.begin(), jets.end(), PtGreaterJet); 

    // get the electrons
    vector<const xAOD::TruthParticle*> electrons;
    const xAOD::TruthParticleContainer* xelectrons = 0;
    if(!event()->retrieve(xelectrions, "TruthElectrons").isSuccess()) {
        cout << MYFUNC << " ERROR: could not retrieve truth electrons" << endl;
        exit(1);
    }
    for(const auto & lep : *xelectrons) {
        if(lep->status() != 1) continue;
        if(std::abs(lep->eta())>2.47) continue;
        if(!(lep->pt() * mev2gev > 10.)) continue;
        electrons.push_back(lep);
    }

    // get the muons

    return true;
}


} // namespace truth

