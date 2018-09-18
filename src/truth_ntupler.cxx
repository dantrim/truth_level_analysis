#include "truth_level_analysis/truth_ntupler.h"

//std/stl
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <math.h>
using namespace std;

//analysis
#include "truth_level_analysis/helpers.h"
#include "truth_level_analysis/utility/string_utils.h"
#include "truth_level_analysis/MT2_ROOT.h"

//xAOD
#include "xAODEventInfo/EventInfo.h"
#include "xAODJet/Jet.h"
#include "xAODJet/JetContainer.h"
#include "xAODMissingET/MissingETContainer.h"

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
        cout << "\n" << MYFUNC << " Output ntuple stored: " << output_filename() << "\n" << endl;
    }
}

void TruthNtupler::setup_output()
{
    cout << MYFUNC << " setting up output ntupler structure" << endl;

    // setup the TTree and TFile
    stringstream ofn;
    ofn << output_dir();
    ofn << "truth_ntup_" << dsid();
    if(suffix() != "") ofn << "_" << suffix();
    ofn << ".root";
    m_output_tree_file = new TFile(ofn.str().c_str(), "RECREATE");
    m_output_filename = ofn.str();

    m_output_tree = new TTree("truth", "truth");

    IBRANCH(dsid)
    FBRANCH(eventweight)
    FBRANCH(w)
    FBRANCH(sumw)

    // leptons
    FBRANCH(l0_pt)
    FBRANCH(l1_pt)
    FBRANCH(l0_eta)
    FBRANCH(l1_eta)
    FBRANCH(l0_phi)
    FBRANCH(l1_phi)
    IBRANCH(isEE)
    IBRANCH(isMM)
    IBRANCH(isEM)
    IBRANCH(isME)
    IBRANCH(isSF)
    IBRANCH(isDF)
    IBRANCH(l0_q)
    IBRANCH(l1_q)
    FBRANCH(pTll)
    FBRANCH(dphi_ll)
    FBRANCH(dRll)
    FBRANCH(mll)

    // jets
    IBRANCH(nJets)
    IBRANCH(nBJets)
    IBRANCH(nSJets)
    FBRANCH(j0_pt)
    FBRANCH(j0_phi)
    FBRANCH(j0_eta)
    FBRANCH(j1_pt)
    FBRANCH(j1_phi)
    FBRANCH(j1_eta)
    FBRANCH(j2_pt)
    FBRANCH(j2_phi)
    FBRANCH(j2_eta)
    FBRANCH(sj0_pt)
    FBRANCH(sj0_phi)
    FBRANCH(sj0_eta)
    FBRANCH(sj1_pt)
    FBRANCH(sj1_phi)
    FBRANCH(sj1_eta)
    FBRANCH(sj2_pt)
    FBRANCH(sj2_phi)
    FBRANCH(sj2_eta)
    FBRANCH(bj0_pt)
    FBRANCH(bj0_phi)
    FBRANCH(bj0_eta)
    FBRANCH(bj1_pt)
    FBRANCH(bj1_phi)
    FBRANCH(bj1_eta)
    FBRANCH(bj2_pt)
    FBRANCH(bj2_phi)
    FBRANCH(bj2_eta)

    // met
    FBRANCH(met)
    FBRANCH(metPhi)

    // leptons + met
    FBRANCH(dr_llmet)
    FBRANCH(dphi_met_ll)
    FBRANCH(met_pTll)

    // di-bjet system
    FBRANCH(mbb)
    FBRANCH(dRbb)
    FBRANCH(dphi_bb)
    FBRANCH(pTbb)

    // leptons + di-bjet system
    FBRANCH(dr_llbb)
    FBRANCH(dphi_llbb)
    FBRANCH(dphi_bj0_ll)
    FBRANCH(dphi_bj0_l0)
    FBRANCH(dr_bj0_l0)
    FBRANCH(dphi_llmet_bb)

    // WWbb variables
    FBRANCH(HT2)
    FBRANCH(HT2Ratio)
    FBRANCH(sumpt)
    FBRANCH(mt2_llbb)
    FBRANCH(mt2_bb)
    FBRANCH(MT_1)
    FBRANCH(MT_1_scaled)


    
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

void TruthNtupler::ivar(string key, int val)
{
    try {
        m_int_vm.at(key) = val;
    }
    catch(std::exception& e)
    {
        cout << MYFUNC << " ERROR: could not fill variable \"" << key << "\": " << e.what() << endl;
        exit(1);
    }
}

void TruthNtupler::fvar(string key, float val)
{
    try {
        m_float_vm.at(key) = val;
    }
    catch(std::exception& e)
    {
        cout << MYFUNC << " ERROR: could not fill variable \"" << key << "\": " << e.what() << endl;
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

bool TruthNtupler::is_bjet(const xAOD::Jet* jet)
{
    int flavor = 0;
    if(jet->isAvailable<int>("HadronConeExclTruthLabelID")) {
        flavor = jet->auxdata<int>("HadronConeExclTruthLabelID");
    }
    else if(jet->isAvailable<int>("ConeTruthLabelID")) {
        flavor = jet->auxdata<int>("ConeTruthLabelID");
    }
    else if(jet->isAvailable<int>("PartonTruthLabelID")) {
        flavor = abs(jet->auxdata<int>("PartonTruthLabelID"));
    }
    else if(jet->isAvailable<int>("GhostBHadronsFinalCount")) {
        if(jet->auxdata<int>("GhostBHadronsFinalCount")) {
            flavor = 5;
        }
        else if(jet->auxdata<int>("GhostCHadronsFinalCount")) {
            flavor = 4;
        }
        else {
            flavor = 1;
        }
    }
    return flavor == 5;
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

    n_evt_processed++;

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
    if(!event()->retrieve(xelectrons, "TruthElectrons").isSuccess()) {
        cout << MYFUNC << " ERROR: could not retrieve truth electrons" << endl;
        exit(1);
    }
    for(const auto & lep : *xelectrons) {
        bool is_prompt = (lep->auxdata<unsigned int>("classifierParticleType" ) == 2);
        if(!is_prompt) continue;
        if(std::abs(lep->eta())>2.47) continue;
        if(!(lep->pt() * mev2gev > 10.)) continue;
        electrons.push_back(lep);
    }

    // get the muons
    vector<const xAOD::TruthParticle*> muons;
    const xAOD::TruthParticleContainer* xmuons = 0;
    if(!event()->retrieve(xmuons, "TruthMuons").isSuccess()) {
        cout << MYFUNC << " ERROR: could not retrieve truth muons" << endl;
        exit(1);
    }
    for(const auto & lep : *xmuons) {
        bool is_prompt = (lep->auxdata<unsigned int>("classifierParticleType" ) == 6);
        if(!is_prompt) continue;
        if(!(lep->pt() * mev2gev > 10.)) continue;
        if(std::abs(lep->eta())>2.4) continue;
        muons.push_back(lep);
    }
    std::sort(electrons.begin(), electrons.end(), PtGreater);
    std::sort(muons.begin(), muons.end(), PtGreater);

    // overlap step between leptons and jets
    for(int ijet = 0; ijet < (int) jets.size(); ijet++) {
        for(int iel = 0; iel < (int) electrons.size(); iel++) {
            if(jets.at(ijet)->p4().DeltaR(electrons.at(iel)->p4()) < 0.2) {
                jets.erase(jets.begin()+ijet);
                ijet--;
                break;
            }
        } // iel
    } // ijet

    for(int iel = 0; iel < (int) electrons.size(); iel++) {
        for(int ijet = 0; ijet < (int) jets.size(); ijet++) {
            if(electrons.at(iel)->p4().DeltaR(jets.at(ijet)->p4()) < 0.4) {
                electrons.erase(electrons.begin() + iel);
                iel--;
                break;
            }
        } // ijet
    } // iel

    for(int imu = 0; imu < (int) muons.size(); imu++) {
        for(int ijet = 0; ijet < (int) jets.size(); ijet++) {
            if(muons.at(imu)->p4().DeltaR(jets.at(ijet)->p4()) < 0.4) {
                muons.erase(muons.begin()+imu);
                imu--;
                break;
            }
        } // ijet
    } // imu

    vector<const xAOD::TruthParticle*> leptons;
    for(auto & el : electrons) leptons.push_back(el);
    for(auto & mu : muons) leptons.push_back(mu);
    std::sort(leptons.begin(), leptons.end(), PtGreater);

    // require at least 2 leptons
    if(!(leptons.size()>=2)) return true;

    vector<const xAOD::Jet*> bjets;
    vector<const xAOD::Jet*> sjets;

    for(const auto j : jets) {
        if(is_bjet(j)) {
            bjets.push_back(j);
        }
        else {
            sjets.push_back(j);
        }
    }
    std::sort(bjets.begin(), bjets.end(), PtGreaterJet);
    std::sort(sjets.begin(), sjets.end(), PtGreaterJet);

    // MET
    const xAOD::MissingETContainer* met = 0;
    RETURN_CHECK(GetName(), event()->retrieve(met, "MET_Truth"));
    TLorentzVector metTLV;
    metTLV.SetPxPyPzE((*met)["NonInt"]->mpx(), (*met)["NonInt"]->mpy(), 0., (*met)["NonInt"]->met());

    fill_ntuple(leptons, jets, sjets, bjets, metTLV);
    output_tree()->Fill();
    n_evt_stored++;
    

    return true;
}

void TruthNtupler::fill_ntuple(vector<const xAOD::TruthParticle*> leptons,
    vector<const xAOD::Jet*> jets, vector<const xAOD::Jet*> sjets, vector<const xAOD::Jet*> bjets,
    TLorentzVector met)
{
    // assume that we have >=2 leptons being provided
    fvar("l0_pt", leptons.at(0)->pt() * mev2gev);
    fvar("l1_pt", leptons.at(1)->pt() * mev2gev);
    fvar("l0_eta", leptons.at(0)->eta());
    fvar("l1_eta", leptons.at(1)->eta());
    fvar("l0_phi", leptons.at(0)->phi());
    fvar("l1_phi", leptons.at(1)->phi());

    int is_ee = 0;
    int is_mm = 0;
    int is_em = 0;
    int is_me = 0;
    int is_sf = 0;
    int is_df = 0;

    is_ee = ( (leptons.at(0)->absPdgId()==11 && leptons.at(1)->absPdgId()==11) ? 1 : 0 );
    is_mm = ( (leptons.at(0)->absPdgId()==13 && leptons.at(1)->absPdgId()==13) ? 1 : 0 );
    is_em = ( (leptons.at(0)->absPdgId()==11 && leptons.at(1)->absPdgId()==13) ? 1 : 0 );
    is_me = ( (leptons.at(0)->absPdgId()==13 && leptons.at(1)->absPdgId()==11) ? 1 : 0 );
    is_sf = (is_ee || is_mm);
    is_df = (is_em || is_me);
    ivar("isEE", is_ee);
    ivar("isMM", is_mm);
    ivar("isEM", is_em);
    ivar("isME", is_me);
    ivar("isSF", is_sf);
    ivar("isDF", is_df);
    ivar("l0_q", (leptons.at(0)->charge()));
    ivar("l1_q", (leptons.at(1)->charge()));

    float dphi_ll = (leptons.at(0)->p4().DeltaPhi(leptons.at(1)->p4()));
    float dRll = (leptons.at(0)->p4().DeltaR(leptons.at(1)->p4()));
    float pTll = (leptons.at(0)->p4() + leptons.at(1)->p4()).Pt() * mev2gev;
    float mll = (leptons.at(0)->p4() + leptons.at(1)->p4()).M() * mev2gev;
    fvar("dphi_ll", dphi_ll);
    fvar("dRll", dRll);
    fvar("pTll", pTll);
    fvar("mll", mll);

    // jets
    int nJets = jets.size();
    int nBJets = bjets.size();
    int nSJets = sjets.size();
    ivar("nJets", nJets);
    ivar("nBJets", nBJets);
    ivar("nSJets", nSJets);

    if(nJets>0) {
        fvar("j0_pt", jets.at(0)->pt() * mev2gev);
        fvar("j0_eta", jets.at(0)->eta());
        fvar("j0_phi", jets.at(0)->phi());
    }
    if(nJets>1) {
        fvar("j1_pt", jets.at(1)->pt() * mev2gev);
        fvar("j1_eta", jets.at(1)->eta());
        fvar("j1_phi", jets.at(1)->phi());
    }
    if(nJets>2) {
        fvar("j2_pt", jets.at(2)->pt() * mev2gev);
        fvar("j2_eta", jets.at(2)->eta());
        fvar("j2_phi", jets.at(2)->phi());
    }
    if(nBJets>0) {
        fvar("bj0_pt", bjets.at(0)->pt() * mev2gev);
        fvar("bj0_eta", bjets.at(0)->eta());
        fvar("bj0_phi", bjets.at(0)->phi());
    }
    if(nBJets>1) {
        fvar("bj1_pt", bjets.at(1)->pt() * mev2gev);
        fvar("bj1_eta", bjets.at(1)->eta());
        fvar("bj1_phi", bjets.at(1)->phi());
    }
    if(nBJets>2) {
        fvar("bj2_pt", bjets.at(2)->pt() * mev2gev);
        fvar("bj2_eta", bjets.at(2)->eta());
        fvar("bj2_phi", bjets.at(2)->phi());
    }
    if(nSJets>0) {
        fvar("sj0_pt", sjets.at(0)->pt() * mev2gev);
        fvar("sj0_eta", sjets.at(0)->eta());
        fvar("sj0_phi", sjets.at(0)->phi());
    }
    if(nSJets>1) {
        fvar("sj1_pt", sjets.at(1)->pt() * mev2gev);
        fvar("sj1_eta", sjets.at(1)->eta());
        fvar("sj1_phi", sjets.at(1)->phi());
    }
    if(nSJets>2) {
        fvar("sj2_pt", sjets.at(2)->pt() * mev2gev);
        fvar("sj2_eta", sjets.at(2)->eta());
        fvar("sj2_phi", sjets.at(2)->phi());
    }

    // met
    fvar("met", met.Pt() * mev2gev);
    fvar("metPhi", met.Phi());

    // met + leptons
    fvar("dr_llmet", (leptons.at(0)->p4() + leptons.at(1)->p4()).DeltaR(met));
    fvar("dphi_met_ll", (leptons.at(0)->p4() + leptons.at(1)->p4()).DeltaPhi(met));
    fvar("met_pTll", (met + leptons.at(0)->p4() + leptons.at(1)->p4()).Pt() * mev2gev);


    // 2 bjets
    if(nBJets>1) {
        fvar("mbb", (bjets.at(0)->p4() + bjets.at(1)->p4()).M() * mev2gev);
        fvar("dRbb", (bjets.at(0)->p4().DeltaR(bjets.at(1)->p4())));
        fvar("dphi_bb", (bjets.at(0)->p4().DeltaPhi(bjets.at(1)->p4())));
        fvar("pTbb", (bjets.at(0)->p4() + bjets.at(1)->p4()).Pt() * mev2gev);

        fvar("dr_llbb", (leptons.at(0)->p4()+leptons.at(1)->p4()).DeltaR((bjets.at(0)->p4()+bjets.at(1)->p4())));
        fvar("dphi_llbb", (leptons.at(0)->p4()+leptons.at(1)->p4()).DeltaPhi((bjets.at(0)->p4()+bjets.at(1)->p4())));
        fvar("dphi_llmet_bb", (leptons.at(0)->p4() + leptons.at(1)->p4() + met).DeltaPhi(bjets.at(0)->p4() + bjets.at(1)->p4()));
        
    }
    if(nBJets>0) {
        fvar("dphi_bj0_ll", (bjets.at(0)->p4().DeltaPhi( (leptons.at(0)->p4()+leptons.at(1)->p4()))));
        fvar("dphi_bj0_l0", (bjets.at(0)->p4().DeltaPhi( leptons.at(0)->p4() ) ));
    }

    // WWbb variables
    if(nBJets>1) {
        float ht2 = (bjets.at(0)->p4()+bjets.at(1)->p4()).Pt() + (met + leptons.at(0)->p4() + leptons.at(1)->p4()).Pt();
        fvar("HT2", ht2 * mev2gev);
        float sumpt = bjets.at(0)->pt() + bjets.at(1)->pt() + met.Pt() + leptons.at(0)->pt() + leptons.at(1)->pt();
        fvar("sumpt", sumpt * mev2gev);
        fvar("HT2Ratio", ht2 / sumpt);

        TLorentzVector vll = (leptons.at(0)->p4() + leptons.at(1)->p4());
        TLorentzVector vbb = (bjets.at(0)->p4() + bjets.at(1)->p4());
        ComputeMT2 calc_mt2llbb = ComputeMT2(vll, vbb, met);
        fvar("mt2_llbb", calc_mt2llbb.Compute() * mev2gev);

        ComputeMT2 calc_mt2bb = ComputeMT2(bjets.at(0)->p4(), bjets.at(1)->p4(), met);
        fvar("mt2_bb", calc_mt2bb.Compute() * mev2gev);

        // MT_1
        TLorentzVector bjet_system = (bjets.at(0)->p4() + bjets.at(1)->p4());
        float mbb = bjet_system.M();
        float scaling = 125.09 / mbb;

        TLorentzVector vis_system = (leptons.at(0)->p4() + leptons.at(1)->p4() + bjet_system);
        float pt_vis = vis_system.Pt();
        float m_vis = vis_system.M();
        float et_vis = sqrt(pt_vis*pt_vis + m_vis*m_vis);
        float x = (et_vis + met.Pt());
        float y = (vis_system + met).Pt();
        fvar("MT_1", sqrt( (x*x) - (y*y) ) * mev2gev);

        // MT_1_scaled
        bjet_system.SetPtEtaPhiE(bjet_system.Pt() * scaling, bjet_system.Eta(),
            bjet_system.Phi(), bjet_system.E() * scaling);
        vis_system = (leptons.at(0)->p4() + leptons.at(1)->p4() + bjet_system);
        pt_vis = vis_system.Pt();
        m_vis = vis_system.M();
        et_vis = sqrt(pt_vis*pt_vis + m_vis*m_vis);
        x = (et_vis + met.Pt());
        y = (vis_system + met).Pt();
        fvar("MT_1_scaled", sqrt( (x*x) - (y*y) ) * mev2gev);
    }

}



} // namespace truth

