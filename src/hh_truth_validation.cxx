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
#include "xAODJet/Jet.h"
#include "xAODJet/JetContainer.h"

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

    add_histogram(Histo("h_has_hbb", ";Event has h#rightarrow bb", 2, 0, 2));
    add_histogram(Histo("h_has_hWW", ";Event has h#rightarrow WW", 2, 0, 2));
    add_histogram(Histo("h_has_hh", ";Event has hh", 2, 0, 2));
    add_histogram(Histo("h_vx_mhh", ";Invariant mass of hh pair [GeV]", 100, 0, -1));
    add_histogram(Histo("h_vx_pThh", ";hh pair p_{T} [GeV]", 100, 0, -1));
    add_histogram(Histo("h_vx_dphi_hh", ";#Delta #phi_{hh}", 100, 0, -1));

    // histograms related to the H->WW or H->bb vertex (i.e. particles directly at the vertices)
    add_histogram(Histo("h_vx_dilep_flavor", ";Dilepton Flavor", DileptonFlavor::Invalid, 0, DileptonFlavor::Invalid));
    add_histogram(Histo("h_vx_h0_pt", ";Leading Higgs p_{T} [GeV]", 100, 0, -1));
    add_histogram(Histo("h_vx_h1_pt", ";Sub-leading Higgs p_{T} [GeV]", 100, 0, -1));
    add_histogram(Histo("h_vx_hWW_pt", ";p_{T} of h->WW Higgs [GeV]", 100, 0, -1));
    add_histogram(Histo("h_vx_hBB_pt", ";p_{T} of h->bb Higgs [GeV]", 100, 0, -1));
    add_histogram(Histo("h_vx_dphi_b0w0", ";#Delta #phi_{W0,b0}", 100, 0, -1));
    add_histogram(Histo("h_vx_mBB", ";Invariant mass of b-quark pair [GeV]", 100, 0, -1));
    add_histogram(Histo("h_vx_mWW", ";Invariant mass of WW pair [GeV]", 100, 0, -1));
    add_histogram(Histo("h_vx_mW0pt", ";Invariant mass of p_{T}-leading W [GeV]", 100, 0, -1));
    add_histogram(Histo("h_vx_mW1pt", ";Invariant mass of p_{T}-subleading W [GeV]", 100, 0, -1));
    add_histogram(Histo("h_vx_mW0p", ";Invariant mass of p-leading W [GeV]", 100, 0, -1));
    add_histogram(Histo("h_vx_mW1p", ";Invariant mass of p-subleading W [GeV]", 100, 0, -1));
    add_histogram(Histo("h2_vx_mW_Wpt", ";Invariant mass of W [GeV]; p_{T} of W [GeV]", 300, 0, -1, 100, 0, -1));
    add_histogram(Histo("h2_vx_eta_mW", ";#eta of W; Invariant mass of W [GeV]", 320, -3.2, 3.2, 300, 0, -1));
    add_histogram(Histo("h2_vx_mW_dRChildren", ";Invariant mass of W [GeV]; #Delta R between W children", 100, 0, -1, 100, 0, -1));
    add_histogram(Histo("h2_vx_mW_dPhiChildren", ";Invariant mass of W [GeV]; #Delta #phi between W children", 100, 0, -1, 100, 0, -1));
    add_histogram(Histo("h_vx_lead_lep_parent", ";p_{T} lead lepton is p_{T} lead W", 2, 0, 2));
    add_histogram(Histo("h_vx_dPhiWW", ";#Delta #phi_{WW}", 60, -3.2,3.2));
    add_histogram(Histo("h_vx_dPhiBB", ";#Delta #phi_{bb}", 60, -3.2,3.2));
    add_histogram(Histo("h_vx_dRWW", ";#Delta R_{WW}", 100, 0, 6));
    add_histogram(Histo("h_vx_dRBB", ";#Delta R_{bb}", 100, 0, 6));

    // histograms from final sate objects
    add_histogram(Histo("h_fs_raw_nLeptons", ";(Raw) Number of leptons (e/#mu)", 6, 0, 6));
    add_histogram(Histo("h_fs_raw_l0_pt", ";(Raw) Leading lepton p_{T} [GeV]", 100, 0, 200));
    add_histogram(Histo("h_fs_raw_l1_pt", ";(Raw) Sub-leading lepton p_{T} [GeV]", 50, 0, 100));
    add_histogram(Histo("h_fs_raw_l0_eta", ";(Raw) Leading lepton #eta", 320, -3.0, 3.0));
    add_histogram(Histo("h_fs_raw_l1_eta", ";(Raw) Sub-leading lepton #eta", 320, -3, 3));
    add_histogram(Histo("h_fs_raw_l0_parent", ";(Raw) Leading lepton parent pdgId", 100, 0, 100));
    add_histogram(Histo("h_fs_raw_l1_parent", ";(Raw) Sub-leading lepton parent pdgId", 100, 0, -1));
}

Bool_t HHTruthValidation::Process(Long64_t entry)
{
    static Long64_t chain_entry = -1;
    chain_entry++;
    event()->getEntry(chain_entry);
    n_evt_processed++;

    if(chain_entry % 1000 == 0) {
        cout << MYFUNC << " Processing entry " << chain_entry << endl;
    }

    if(chain_entry == 0) {
        setup_output();
        setup_histograms();
    }

    // set the MC weight
    const xAOD::EventInfo* ei = 0;
    RETURN_CHECK(GetName(), event()->retrieve(ei, "EventInfo"));
    mc_weight(ei->mcEventWeight());


    // let's just build up the collections here
    bool event_ok = process_strict_truth();
    event_ok = process_final_state();


    //for(const auto & p : *truthElectrons) {
    //    fill_histo("l0_pt", p->pt() * mev2gev);
    //}

    return true;
}

bool HHTruthValidation::process_strict_truth()
{
    // here we directly access the higgs bosons and their immediate children
    // at their production vertex

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

    const xAOD::TruthParticle* higgs_bosons[2] = { nullptr };

    HiggsCandidate hWW;
    HiggsCandidate hBB;

    vector<int> children;
    int n_higgs_candidates = 0;
    int n_hww = 0;
    int n_hbb = 0;

    for(const auto & vertex : *truthVertices) {
        if(vertex->nIncomingParticles() == 1 && vertex->nOutgoingParticles() >= 2) {
            const xAOD::TruthParticle* mother = vertex->incomingParticle(0);

            if(mother && mother->absPdgId() == 25) {
                if(debug_level() >= 5) {
                    cout << MYFUNC << " mother higgs nChildren = " << mother->nChildren() << " : child pdgs [ ";                
                    for(int ichild = 0; ichild < (int)mother->nChildren(); ichild++) {
                        cout << mother->child(ichild)->pdgId() << " ";
                    } // ichild
                    cout << "]" << endl;
                }

                n_higgs_candidates++;
                if(mother->nChildren() >= 2) {

                    // locate the children of the higgs that are either W's or b's
                    vector<int> child_idx;
                    const xAOD::TruthParticle* child0 = nullptr;
                    const xAOD::TruthParticle* child1 = nullptr;

                    for(int ichild = 0; ichild < (int)mother->nChildren(); ichild++) {
                        if(mother->child(ichild)->absPdgId() == 24 ||
                            mother->child(ichild)->absPdgId() == 5) {
                            child_idx.push_back(ichild);
                        }
                    } // ichild

                    // order the children by pT
                    int idx0 = child_idx.at(0);
                    int idx1 = child_idx.at(1);
                    if(mother->child(idx0)->pt() > mother->child(idx1)->pt()) {
                        child0 = mother->child(idx0);
                        child1 = mother->child(idx1);
                    }
                    else {
                        child0 = mother->child(idx1);
                        child1 = mother->child(idx0);
                    }

                    bool is_bb = (child0->absPdgId() == 5 && child1->absPdgId() == 5);
                    bool is_ww = (child0->absPdgId() == 24 && child1->absPdgId() == 24);

                    if(is_bb) {
                        n_hbb++;
                        if(hBB.parent != nullptr) {
                            cout << MYFUNC << " WARNING Found multiple H -> bb candidates in the event" << endl;
                        }

                        hBB.parent = mother;
                        hBB.child0 = child0;
                        hBB.child1 = child1;
                    }
                    else if(is_ww) {
                        n_hww++;
                        if(hWW.parent != nullptr) {
                            cout << MYFUNC << " WARNING Found multiple H -> WW candidates in the event" << endl;
                        }

                        hWW.parent = mother;
                        hWW.child0 = child0;
                        hWW.child1 = child1;
                    }
                    else {
                        cout << MYFUNC << " WARNING Found a non H -> WW or H -> bb decay (first 2 children are: " 
                            << child0->pdgId() << ", " << child1->pdgId() << ")" << endl;
                    }
                } // >=2 children
                else {
                    cout << MYFUNC << " WARNING Found a Higgs that does NOT have 2 children (number of children: "
                            << mother->nChildren() << ")" << endl;
                }
            } // mother is higgs
        }
    } // vertex

    bool n_higgs_ok = (n_higgs_candidates==2);
    bool higgs_decays_ok = (n_hww==1 && n_hbb==1);

    // from the H->WW candidate, fill get the lepton-neutrino system for the each of the W's
    WCandidate w0;
    WCandidate w1;
    bool w_ok = true;
    if(hWW.both_children()) {
        if(!get_w_decay(hWW.child0, w0)) { w_ok = false; cerr << MYFUNC << " WARNING Problem while getting W decay" << endl; }
        if(!get_w_decay(hWW.child1, w1)) { w_ok = false; cerr << MYFUNC << " WARNING Problem while getting W decay" << endl; }
    }

    fill_histo("h_has_hbb", (n_hbb == 1 ? 1 : 0));
    fill_histo("h_has_hWW", (n_hww == 1 ? 1 : 0));
    fill_histo("h_has_hh", (n_higgs_candidates == 2 ? 1 : 0));

    bool event_ok = (n_higgs_ok && higgs_decays_ok && w_ok);

    if(!event_ok) return true;


    fill_histo("h_vx_mhh", (hBB.parent->p4() + hWW.parent->p4()).M() * mev2gev);
    fill_histo("h_vx_pThh", (hBB.parent->p4() + hWW.parent->p4()).Pt() * mev2gev);

    // fill histograms related to the particles directly at the higgs decay vertices
    if(hWW.parent->pt() > hBB.parent->pt()) {
        fill_histo("h_vx_h0_pt", hWW.parent->pt() * mev2gev);
        fill_histo("h_vx_h1_pt", hBB.parent->pt() * mev2gev);
    }
    else {
        fill_histo("h_vx_h0_pt", hBB.parent->pt() * mev2gev);
        fill_histo("h_vx_h1_pt", hWW.parent->pt() * mev2gev);
    }

    fill_histo("h_vx_hWW_pt", hWW.parent->pt() * mev2gev);
    fill_histo("h_vx_hBB_pt", hBB.parent->pt() * mev2gev);
    fill_histo("h_vx_dphi_hh", hWW.parent->p4().DeltaPhi(hBB.parent->p4()));
    fill_histo("h_vx_dphi_b0w0", hWW.child0->p4().DeltaPhi(hBB.child0->p4()));
    fill_histo("h_vx_mBB", (hBB.child0->p4() + hBB.child1->p4()).M() * mev2gev);
    fill_histo("h_vx_mWW", (hWW.child0->p4() + hWW.child1->p4()).M() * mev2gev);
    fill_histo("h_vx_dPhiWW", (hWW.child0->p4().DeltaPhi(hWW.child1->p4())));
    fill_histo("h_vx_dPhiBB", (hBB.child0->p4().DeltaPhi(hBB.child1->p4())));
    fill_histo("h_vx_dRWW", (hWW.child0->p4().DeltaR(hWW.child1->p4())));
    fill_histo("h_vx_dRBB", (hBB.child0->p4().DeltaR(hBB.child1->p4())));
    fill_histo("h_vx_mW0pt", hWW.child0->p4().M() * mev2gev);
    fill_histo("h_vx_mW1pt", hWW.child1->p4().M() * mev2gev);
    if(hWW.child0->p4().P() > hWW.child1->p4().P()) {
        fill_histo("h_vx_mW0p", hWW.child0->p4().M() * mev2gev);
        fill_histo("h_vx_mW1p", hWW.child1->p4().M() * mev2gev);
    }
    else {
        fill_histo("h_vx_mW0p", hWW.child1->p4().M() * mev2gev);
        fill_histo("h_vx_mW1p", hWW.child0->p4().M() * mev2gev);
    }
    fill_histo("h2_vx_mW_Wpt", hWW.child0->pt() * mev2gev, hWW.child0->p4().M() * mev2gev);
    fill_histo("h2_vx_mW_Wpt", hWW.child1->pt() * mev2gev, hWW.child1->p4().M() * mev2gev);
    fill_histo("h2_vx_eta_mW", hWW.child0->eta(), hWW.child0->p4().M() * mev2gev);
    fill_histo("h2_vx_eta_mW", hWW.child1->eta(), hWW.child1->p4().M() * mev2gev);
    fill_histo("h2_vx_mW_dRChildren", w0.w->p4().M() * mev2gev, w0.lepton->p4().DeltaR(w0.neutrino->p4()));
    fill_histo("h2_vx_mW_dRChildren", w1.w->p4().M() * mev2gev, w1.lepton->p4().DeltaR(w1.neutrino->p4()));
    fill_histo("h2_vx_mW_dPhiChildren", w0.w->p4().M() * mev2gev, w0.lepton->p4().DeltaPhi(w0.neutrino->p4()));
    fill_histo("h2_vx_mW_dPhiChildren", w1.w->p4().M() * mev2gev, w1.lepton->p4().DeltaPhi(w1.neutrino->p4()));

    // let's fill these histograms unweighted, so set the weight to 1
    mc_weight(1);
    if(w0.lepton->pt() > w1.lepton->pt()) {
        fill_histo("h_vx_lead_lep_parent", 1);
        fill_histo("h_vx_dilep_flavor", (int)get_dilepton_flavor(w0.lepton, w1.lepton));
    }
    else {
        fill_histo("h_vx_lead_lep_parent", 0);
        fill_histo("h_vx_dilep_flavor", (int)get_dilepton_flavor(w1.lepton, w0.lepton));
    }
    // don't forget to reset the weight back to the MC gen weight
    mc_weight(ei->mcEventWeight());

    return true;
}

bool HHTruthValidation::get_w_decay(const xAOD::TruthParticle* wboson, WCandidate& w)
{
    if(wboson == nullptr) {
        if(debug_level() >= 1) { cerr << MYFUNC << " WARNING Input parent W is null" << endl; }
        return false;
    }
    if(!wboson->isW()) {
        if(debug_level() >= 1) { cerr << MYFUNC << " WARNING Input parent W is not a W" << endl; }
        return false;
    }

    // walk down the parent W boson's decay until you get to the lepton-neutrino vertex
    const xAOD::TruthParticle* wint = wboson;
    const xAOD::TruthParticle* neutrino = nullptr;
    const xAOD::TruthParticle* lepton = nullptr;

    bool self_decay = false;
    do {
        for(size_t i = 0; i < wint->nChildren(); i++) {
            self_decay = false;
            const xAOD::TruthParticle* child = wint->child(i);
            if(child->pdgId() == wint->pdgId()) {
                self_decay = true;
                wint = child;
                break;
            } // child has same pdgID
            else if(child->isChLepton()) {
                lepton = child;
            } // child is a charged lepton
            else if(child->isNeutrino()) {
                neutrino = child;
            } // child is neutrino
        } // i
    }
    while (self_decay);

    w.w = wint;
    w.lepton = lepton;
    w.neutrino = neutrino;

    return true;
}

bool HHTruthValidation::process_final_state()
{
    const xAOD::EventInfo* ei = 0;
    RETURN_CHECK(GetName(), event()->retrieve(ei, "EventInfo"));

    const xAOD::TruthParticleContainer* truthParticles = 0;
    RETURN_CHECK(GetName(), event()->retrieve(truthParticles, "TruthParticles"));

    const xAOD::TruthParticleContainer* truthElectrons = 0;
    RETURN_CHECK(GetName(), event()->retrieve(truthElectrons, "TruthElectrons"));

    const xAOD::TruthParticleContainer* truthMuons = 0;
    RETURN_CHECK(GetName(), event()->retrieve(truthMuons, "TruthMuons"));

    const xAOD::JetContainer* truthJets = 0;
    //RETURN_CHECK(GetName(), event()->retrieve(truthJets, "AntiKt4TruthDressedWZJets"));
    RETURN_CHECK(GetName(), event()->retrieve(truthJets, "AntiKt4TruthWZJets"));

    const xAOD::MissingETContainer* met = 0;
    RETURN_CHECK(GetName(), event()->retrieve(met, "MET_Truth" ));

    // first don't apply any selection to just see how the filter is behaving on these "raw" objects
    vector<const xAOD::TruthParticle*> leptons;
    for(const auto & lep : *truthElectrons) {
        if(lep->status() != 1) continue;
        leptons.push_back(lep);
    }
    for(const auto & lep : *truthMuons) {
        if(lep->status() != 1) continue;
        leptons.push_back(lep);
    }
    std::sort(leptons.begin(), leptons.end(), PtGreater);

    fill_histo("h_fs_raw_nLeptons", leptons.size());
    if(leptons.size()>0) {
        fill_histo("h_fs_raw_l0_pt", leptons.at(0)->pt() * mev2gev); 
        fill_histo("h_fs_raw_l0_eta", leptons.at(0)->eta());
        if(leptons.at(0)->parent(0)) fill_histo("h_fs_raw_l0_parent", (int)leptons.at(0)->parent(0)->absPdgId());
    }
    if(leptons.size()>1) {
        fill_histo("h_fs_raw_l1_pt", leptons.at(1)->pt() * mev2gev);
        fill_histo("h_fs_raw_l1_eta", leptons.at(1)->eta());
        if(leptons.at(1)->parent(0)) fill_histo("h_fs_raw_l1_parent", (int)leptons.at(1)->parent(0)->absPdgId());
    }

    return true;
}

} // namespace
