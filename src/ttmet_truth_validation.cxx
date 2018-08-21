#include "truth_level_analysis/ttmet_truth_validation.h"

//std/stl
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
using namespace std;

//math
#include <cmath>

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

TtMetTruthValidation::TtMetTruthValidation() :
    m_mc_weight(1.0),
    m_outfilename(""),
    m_rfile(nullptr)
{
}

void TtMetTruthValidation::SlaveBegin(TTree* /*tree*/)
{
    timer()->Start();
}

void TtMetTruthValidation::Terminate()
{
    TruthSelectorBase::Terminate();
    save_histograms(m_rfile);
}

void TtMetTruthValidation::setup_output()
{
    stringstream ofn;
    ofn << output_dir() << "/" << "ttmet_validation";
    if(suffix() != "") ofn << "_" << suffix();
    ofn << ".root";
    cout << MYFUNC << " setting up output file (name = " << ofn.str() << ")" << endl;
    m_rfile = new TFile(ofn.str().c_str(), "RECREATE");
    if(m_rfile->IsZombie()) {
        cout << MYFUNC << " ERROR: could not create the output file" << endl;
        exit(1);
    }
    m_rfile->cd();
}

void TtMetTruthValidation::setup_histograms()
{
    add_histogram(Histo("h_wPt", ";W boson p_{T} [GeV]", 100, 0, -1));
    add_histogram(Histo("h_wPt_children", ";(l+#nu) p_{T} [GeV]", 100, 0, -1));
    add_histogram(Histo("h_wMass", ";W boson mass [GeV]", 100, 0, -1));
    add_histogram(Histo("h_wMass_children", ";(l+#nu) mass [GeV]", 100, 0, -1));

    // stop pair stuff
    add_histogram(Histo("h_st0_pT", ";#tilde{t} p_{T} 0 [GeV]", 100, 0, -1));
    add_histogram(Histo("h_st1_pT", ";#tilde{t} p_{T} 1 [GeV]", 100, 0, -1));
    add_histogram(Histo("h_st0_m", ";#tilde{t} 0 mass [GeV]", 100, 0, -1));
    add_histogram(Histo("h_st1_m", ";#tilde{t} 1 mass [GeV]", 100, 0, -1));

    add_histogram(Histo("h_n10_pT", ";#tilde{#chi}_{1}^{0} 0 p_{T} [GeV]", 100, 0, -1));
    add_histogram(Histo("h_n11_pT", ";#tilde{#chi}_{1}^{0} 1 p_{T} [GeV]", 100, 0, -1));
    add_histogram(Histo("h_n10_m", ";#tilde{#chi}_{1}^{0} 0 mass [GeV]", 100, 0, -1));
    add_histogram(Histo("h_n11_m", ";#tilde{#chi}_{1}^{0} 1 mass [GeV]", 100, 0, -1));

    add_histogram(Histo("h_stst_pT", ";p_{T}^{#tilde{t} #tilde{t}} [GeV]", 100, 0, -1));
    add_histogram(Histo("h_stst_m", ";m_{#tilde{t} #tilde{t}} [GeV]", 100, 0, -1));
    add_histogram(Histo("h_n1n1_pT", ";p_{T}^{#tilde{#chi} #tilde{#chi}} [GeV]", 100, 0, -1));
    add_histogram(Histo("h_n1n1_m", ";m_{#tilde{#chi} #tilde{#chi}} [GeV]", 100, 0, -1));
    add_histogram(Histo("h_dphi_stst", ";#Delta #phi_{#tilde{t} #tilde{t}}", 100, 0, -1));
    add_histogram(Histo("h_dphi_n1n1", ";#Delta #phi_{#tilde{#chi} #tilde{#chi}}", 100, 0, -1));

    add_histogram(Histo("h_top0_m", ";t-quark 0 mass [GeV]", 100, 0, -1));
    add_histogram(Histo("h_top1_m", ";t-quark 1 mass [GeV]", 100, 0, -1));
    add_histogram(Histo("h_dphi_tt", ";#Delta #phi_{tt}", 100, 0, -1));
    add_histogram(Histo("h_b0_pT", ";b-quark 0 p_{T} [GeV]", 100, 0, -1));
    add_histogram(Histo("h_b1_pT", ";b-quark 1 p_{T} [GeV]", 100, 0, -1));

    add_histogram(Histo("h_wb0_pT", ";W+b 0 (l+#nu+b) p_{T} [GeV]", 100, 0, -1));
    add_histogram(Histo("h_wb0_m", ";W+b 0 mass (l+#nu+b) mass [GeV]", 100, 0, -1));
    add_histogram(Histo("h_cosThetal0", ";cos #theta_l 0", 100, 0, -1));
    add_histogram(Histo("h_wb1_pT", ";W+b 1 (l+#nu+b) p_{T} [GeV]", 100, 0, -1));
    add_histogram(Histo("h_wb1_m", ";W+b 1 mass (l+#nu+b) mass [GeV]", 100, 0, -1));
    add_histogram(Histo("h_cosThetal1", ";cos #theta_l 1", 100, 0, -1));

    // "final state"
    add_histogram(Histo("h_fs_met", ";Missing Transverse Momentum (Non-Int) [GeV]", 100, 0, -1));
    add_histogram(Histo("h_fs_l0_pT", ";Lead lepton p_{T} [GeV]", 100, 0, -1));
    add_histogram(Histo("h_fs_l1_pT", ";Sub-lead lepton p_{T} [GeV]", 100, 0, -1));
    add_histogram(Histo("h_fs_pTll", ";p_{T}^{ll} [GeV]", 100, 0, -1));
    add_histogram(Histo("h_fs_dphi_ll", ";#Delta #phi_{ll}", 100, 0, -1));
    add_histogram(Histo("h_fs_dphi_met_ll", ";#Delta #phi(MET,ll)", 100, 0, -1));

    add_histogram(Histo("h_fs_nJets", ";Jet multiplicity", 20, 0, 20));
    add_histogram(Histo("h_fs_nBJets", ";Truth b-Jet multiplicity", 10, 0, 10));
    add_histogram(Histo("h_fs_nSJets", ";Truth non-b-Jet multiplicity", 15, 0, 15));
    add_histogram(Histo("h_fs_j0_pT", ";Lead jet p_{T} [GeV]", 100, 0, -1));
    add_histogram(Histo("h_fs_j1_pT", ";Sub-lead jet p_{T} [GeV]", 100, 0, -1));
    add_histogram(Histo("h_fs_bj0_pT", ";Lead b-jet p_{T} [GeV]", 100, 0, -1));
    add_histogram(Histo("h_fs_bj1_pT", ";Sub-lead b-jet p_{T} [GeV]", 100, 0, -1));

    add_histogram(Histo("h_fs_dphi_l0_j0", ";#Delta #phi(lead-lep, lead-jet)", 100, 0, -1));
    add_histogram(Histo("h_fs_dphi_l0_bj0", ";#Delta #phi(lead-lep, lead b-jet)", 100, 0, -1));
    add_histogram(Histo("h_fs_dphi_ll_j0", ";#Delta #phi(ll, lead-jet)", 100, 0, -1));
    add_histogram(Histo("h_fs_dphi_ll_bj0", ";#Delta #phi(ll, lead b-jet)", 100, 0, -1));
}

Bool_t TtMetTruthValidation::Process(Long64_t entry)
{
    static Long64_t chain_entry = -1;
    chain_entry++;
    event()->getEntry(chain_entry);
    n_evt_processed++;

    if(chain_entry == 0) {
        setup_output();
        setup_histograms();
    }

    if(chain_entry % 500  == 0) {
        cout << MYFUNC << " Processing entry " << chain_entry << endl;
    }

    const xAOD::EventInfo* ei = 0;
    RETURN_CHECK(GetName(), event()->retrieve(ei, "EventInfo"));
    mc_weight(ei->mcEventWeight()); // this propagates the event weight value to later calls

    bool ok = process_stop_decays();
    if(!ok) return true;
    ok = process_final_state();

    return true;
}

bool TtMetTruthValidation::process_stop_decays()
{
    const xAOD::TruthVertexContainer* truthVertices = 0;
    RETURN_CHECK(GetName(), event()->retrieve(truthVertices, "TruthVertices"));

    //_________________________________________________________________________
    // find the W bosons and their decays

    const xAOD::TruthParticle *wboson = nullptr,
                                *wlepton = nullptr,
                                *wneutrino = nullptr;
    const xAOD::TruthParticle *wbosons[2] = { nullptr },
                                *wleptons[2] = { nullptr },
                                *wneutrinos[2] = { nullptr };
    bool self_decay = false;
    for(const auto & vertex : *truthVertices) {
        if(vertex->nIncomingParticles() == 1 && vertex->nOutgoingParticles() == 2) {
            const xAOD::TruthParticle* mother = vertex->incomingParticle(0);
            if(!mother->isW()) continue;
            unsigned int fill_idx = mother->pdgId() > 0 ? 0 : 1;
            wbosons[fill_idx] = mother;
        } // 1 --> 2
    } // vertex

    if(wbosons[0] == nullptr || wbosons[1] == nullptr) {
        return false;
    }

    // find W decay products
    for(unsigned int iw = 0; iw < 2; iw++) {
        if(wbosons[iw] != nullptr) {
            self_decay = false;
            wboson = wbosons[iw];
            do {
                for(unsigned int i = 0; i < wboson->nChildren(); i++) {
                    self_decay = false;
                    const xAOD::TruthParticle* child = wboson->child(i);
                    if(child->pdgId() == wboson->pdgId()) {
                        wboson = child;
                        self_decay = true;
                        break;
                    }
                    else if(child->isChLepton()) {
                        wlepton = child;
                        wleptons[iw] = wlepton;
                    }
                    else if(child->isNeutrino()) {
                        wneutrino = child;
                        wneutrinos[iw] = wneutrino;
                    }
                } // i
            } while(self_decay);
        } // wboson is not null
    } // iw

    if(wleptons[0] == nullptr || wleptons[1] == nullptr || wneutrinos[0] == nullptr || wneutrinos[1] == nullptr) {
        return false;
    }

    for(unsigned int i = 0; i < 2; i++) {
        fill_histo("h_wPt", wbosons[i]->pt() * mev2gev);
        fill_histo("h_wMass", wbosons[i]->m() * mev2gev);
        fill_histo("h_wPt_children", (wleptons[i]->p4() + wneutrinos[i]->p4()).Pt() * mev2gev);
        fill_histo("h_wMass_children", (wleptons[i]->p4() + wneutrinos[i]->p4()).M() * mev2gev);
    }

    //_________________________________________________________________________
    // find the stops and their decays
    // TODO: handle cases where this IS NOT stop signal
    const xAOD::TruthParticle *stop = nullptr,
                                *bquark = nullptr,
                                *neutralino = nullptr;
                                //*wboson = nullptr,
                                //*wlepton = nullptr,
                                //*wneutrino = nullptr;
    wboson = nullptr;
    const xAOD::TruthParticle *stops[2] = { nullptr },
                                *bquarks[2] =  { nullptr },
                                *neutralinos[2] = { nullptr };
                                //*wbosons[2] = { nullptr },
                                //*wleptons[2] = { nullptr },
                                //*wneutrinos[2] = { nullptr };
    wbosons[2] = { nullptr };
    wleptons[2] = { nullptr };
    wneutrinos[2] = { nullptr };
    self_decay = false;

    for(const auto & vertex : * truthVertices) {
        if(vertex->nIncomingParticles() == 1 && vertex->nOutgoingParticles() == 3) {
            const xAOD::TruthParticle* mother = vertex->incomingParticle(0);
            if(mother->absPdgId() != 1000006) continue; // stops only
            unsigned int fill_idx = mother->pdgId() > 0 ? 0 : 1;
            stops[fill_idx] = mother;
            for(unsigned int i = 0; i < 3; i++) {
                const xAOD::TruthParticle* child = vertex->outgoingParticle(i);
                if(child->absPdgId() == 1000022) { neutralinos[fill_idx] = child; }
                else if(child->absPdgId() == 5) { bquarks[fill_idx] = child; }
                else if(child->isW()) { wbosons[fill_idx] = child; }
            } // i
        } // 1 --> 3
    } // vertex

    if(stops[0] == nullptr || stops[1] == nullptr ||
        wbosons[0] == nullptr || wbosons[1] == nullptr ||
        neutralinos[0] == nullptr || neutralinos[1] == nullptr ||
        bquarks[0] == nullptr || bquarks[1] == nullptr) {

        cout << MYFUNC << " WARNING did not find all stop decay products" << endl;
        return false;
    }

    // find W decay products
    for(unsigned int iw = 0; iw < 2; iw++) {
        if(wbosons[iw] != nullptr) {
            self_decay = false;
            wboson = wbosons[iw];
            do {
                for(unsigned int i = 0; i < wboson->nChildren(); i++) {
                    self_decay = false;
                    const xAOD::TruthParticle* child = wboson->child(i);
                    if(child->pdgId() == wboson->pdgId()) {
                        wboson = child;
                        self_decay = true;
                        break;
                    }
                    else if(child->isChLepton()) {
                        wlepton = child;
                        wleptons[iw] = wlepton;
                    }
                    else if(child->isNeutrino()) {
                        wneutrino = child;
                        wneutrinos[iw] = wneutrino;
                    }
                } // i
            } while(self_decay);
        } // wboson is not null
    } // iw

    TLorentzVector top_system0 = wleptons[0]->p4() + wneutrinos[0]->p4() + bquarks[0]->p4();
    TLorentzVector top_system1 = wleptons[1]->p4() + wneutrinos[1]->p4() + bquarks[1]->p4();
    TLorentzVector lep0 = wleptons[0]->p4();
    TLorentzVector lep1 = wleptons[1]->p4();
    TLorentzVector chi0 = neutralinos[0]->p4();
    TLorentzVector chi1 = neutralinos[1]->p4();

    fill_histo("h_wb0_pT", top_system0.Pt() * mev2gev);
    fill_histo("h_wb1_pT", top_system1.Pt() * mev2gev);
    fill_histo("h_wb0_m", top_system0.M() * mev2gev);
    fill_histo("h_wb1_m", top_system1.M() * mev2gev);
    
    TVector3 boost0 = top_system0.BoostVector();
    TVector3 boost1 = top_system1.BoostVector();

    lep0.Boost(-boost0);
    chi0.Boost(-boost0);
    lep1.Boost(-boost1);
    chi1.Boost(-boost1);

    fill_histo("h_cosThetal0", cos( chi0.Angle(lep0.Vect()) ));
    fill_histo("h_cosThetal1", cos( chi1.Angle(lep1.Vect()) ));

    fill_histo("h_st0_pT", stops[0]->p4().Pt() * mev2gev);
    fill_histo("h_st1_pT", stops[1]->p4().Pt() * mev2gev);
    fill_histo("h_st0_m", stops[0]->p4().M() * mev2gev);
    fill_histo("h_st1_m", stops[1]->p4().M() * mev2gev);
    fill_histo("h_n10_pT", neutralinos[0]->p4().Pt() * mev2gev);
    fill_histo("h_n11_pT", neutralinos[1]->p4().Pt() * mev2gev);
    fill_histo("h_n10_m", neutralinos[0]->p4().M() * mev2gev);
    fill_histo("h_n11_m", neutralinos[1]->p4().M() * mev2gev);

    fill_histo("h_stst_pT", (stops[0]->p4() + stops[1]->p4()).Pt() * mev2gev);
    fill_histo("h_stst_m", (stops[0]->p4() + stops[1]->p4()).M() * mev2gev);
    fill_histo("h_n1n1_pT", (neutralinos[0]->p4() + neutralinos[1]->p4()).Pt() * mev2gev);
    fill_histo("h_n1n1_m", (neutralinos[0]->p4() + neutralinos[1]->p4()).M() * mev2gev);
    fill_histo("h_dphi_stst", stops[0]->p4().DeltaPhi(stops[1]->p4()));
    fill_histo("h_dphi_n1n1", neutralinos[0]->p4().DeltaPhi(neutralinos[1]->p4()));

    return true;
}

bool TtMetTruthValidation::process_final_state()
{

    const xAOD::TruthParticleContainer* truthElectrons = 0;
    RETURN_CHECK(GetName(), event()->retrieve(truthElectrons, "TruthElectrons"));
    const xAOD::TruthParticleContainer* truthMuons = 0;
    RETURN_CHECK(GetName(), event()->retrieve(truthMuons, "TruthMuons"));
    const xAOD::JetContainer* truthJets = 0;
    RETURN_CHECK(GetName(), event()->retrieve(truthJets, "AntiKt4TruthDressedWZJets"));
    const xAOD::MissingETContainer* met = 0;
    RETURN_CHECK(GetName(), event()->retrieve(met, "MET_Truth"));
    TLorentzVector metTLV;
    metTLV.SetPxPyPzE((*met)["NonInt"]->mpx() * mev2gev, (*met)["NonInt"]->mpy() * mev2gev, 0., (*met)["NonInt"]->met() * mev2gev);

    vector<const xAOD::Jet*> jets;
    for(const auto & j : *truthJets) {
        if(!(j->p4().Pt() > 20. * mev2gev)) continue;
        jets.push_back(j);
    }
    std::sort(jets.begin(), jets.end(), PtGreaterJet);

    vector<const xAOD::TruthParticle*> electrons;
    for(const auto & lep : *truthElectrons) {
        if(lep->status() != 1) continue;
        if(std::abs(lep->eta())>2.47) continue;
        if(!(lep->pt()*mev2gev > 10.)) continue;
        electrons.push_back(lep);
    } // lep
    vector<const xAOD::TruthParticle*> muons;
    for(const auto & lep : *truthMuons) {
        if(lep->status() != 1) continue;
        if(!(lep->pt()*mev2gev>10.)) continue;
        if(std::abs(lep->eta())>2.4) continue;
        muons.push_back(lep);
    }
    std::sort(electrons.begin(), electrons.end(), PtGreater);
    std::sort(muons.begin(), muons.end(), PtGreater);

    // overlap removal to get rid of duplicate objects
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

    // require 2 leptons
    if(!(leptons.size()>=2)) return true;

    vector<const xAOD::Jet*> bjets;
    vector<const xAOD::Jet*> sjets;
    for(const auto j : jets) {
        int flavor = std::abs(j->auxdata<int>("PartonTruthLabelID"));
        if(flavor==5) {
            bjets.push_back(j);
        }
        else {
            sjets.push_back(j);
        }
    }
    std::sort(bjets.begin(), bjets.end(), PtGreaterJet);
    std::sort(sjets.begin(), sjets.end(), PtGreaterJet);

    fill_histo("h_fs_met", metTLV.Pt());
    fill_histo("h_fs_l0_pT", leptons.at(0)->p4().Pt() * mev2gev);
    fill_histo("h_fs_l1_pT", leptons.at(1)->p4().Pt() * mev2gev);
    fill_histo("h_fs_pTll", (leptons.at(0)->p4() + leptons.at(1)->p4()).Pt() * mev2gev);
    fill_histo("h_fs_dphi_ll", leptons.at(0)->p4().DeltaPhi(leptons.at(1)->p4()));
    fill_histo("h_fs_dphi_met_ll", metTLV.DeltaPhi( (leptons.at(0)->p4() + leptons.at(1)->p4() )));

    fill_histo("h_fs_nJets", jets.size());
    fill_histo("h_fs_nBJets", bjets.size());
    fill_histo("h_fs_nSJets", sjets.size());
    if(jets.size()>0) {
        fill_histo("h_fs_j0_pT", jets.at(0)->pt() * mev2gev);
        fill_histo("h_fs_dphi_l0_j0", leptons.at(0)->p4().DeltaPhi(jets.at(0)->p4()));
        fill_histo("h_fs_dphi_ll_j0", (leptons.at(0)->p4() + leptons.at(1)->p4()).DeltaPhi(jets.at(0)->p4()));
    }
    if(jets.size()>1) {
        fill_histo("h_fs_j1_pT", jets.at(1)->pt() * mev2gev);
    }

    if(bjets.size()>0) {
        fill_histo("h_fs_bj0_pT", bjets.at(0)->pt() * mev2gev);
        fill_histo("h_fs_dphi_l0_bj0", leptons.at(0)->p4().DeltaPhi(bjets.at(0)->p4()));
        fill_histo("h_fs_dphi_ll_bj0", (leptons.at(0)->p4()+leptons.at(1)->p4()).DeltaPhi(bjets.at(0)->p4()));
    }
    if(bjets.size()>1) {
        fill_histo("h_fs_bj1_pT", bjets.at(1)->pt() * mev2gev);
    }

    return true;
}


} // namespace truth
