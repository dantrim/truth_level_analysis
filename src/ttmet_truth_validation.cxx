#include "truth_level_analysis/ttmet_truth_validation.h"

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
    add_histogram(Histo("h_st_pT", ";#tilde{t} p_{T} [GeV]", 100, 0, -1));
    add_histogram(Histo("h_st_m", ";#tilde{t} mass [GeV]", 100, 0, -1));
    add_histogram(Histo("h_n1_pT", ";#tilde{#chi}_{1}^{0} p_{T} [GeV]", 100, 0, -1));
    add_histogram(Histo("h_n1_m", ";#tilde{#chi}_{1}^{0} p_{T} [GeV]", 100, 0, -1));
    add_histogram(Histo("h_stst_pT", ";p_{T}^{#tilde{t} #tilde{t}} [GeV]", 100, 0, -1));
    add_histogram(Histo("h_stst_m", ";m_{#tilde{t} #tilde{t}} [GeV]", 100, 0, -1));
    add_histogram(Histo("h_n1n1_pT", ";p_{T}^{#tilde{#chi} #tilde{#chi}} [GeV]", 100, 0, -1));
    add_histogram(Histo("h_n1n1_m", ";m_{#tilde{#chi} #tilde{#chi}} [GeV]", 100, 0, -1));
    add_histogram(Histo("h_dphi_stst", ";#Delta #phi_{#tilde{t} #tilde{t}}", 100, 0, -1));
    add_histogram(Histo("h_dphi_n1n1", ";#Delta #phi_{#tilde{#chi} #tilde{#chi}}", 100, 0, -1));
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

    for(int i = 0; i < 2; i++) {
        fill_histo("h_st_pT", stops[i]->pt() * mev2gev);
        fill_histo("h_st_m", stops[i]->m() * mev2gev);
        fill_histo("h_n1_pT", neutralinos[i]->pt() * mev2gev);
        fill_histo("h_n1_m", neutralinos[i]->m() * mev2gev);
    }
    fill_histo("h_stst_pT", (stops[0]->p4() + stops[1]->p4()).Pt() * mev2gev);
    fill_histo("h_stst_m", (stops[0]->p4() + stops[1]->p4()).M() * mev2gev);
    fill_histo("h_n1n1_pT", (neutralinos[0]->p4() + neutralinos[1]->p4()).Pt() * mev2gev);
    fill_histo("h_n1n1_m", (neutralinos[0]->p4() + neutralinos[1]->p4()).M() * mev2gev);
    fill_histo("h_dphi_stst", stops[0]->p4().DeltaPhi(stops[1]->p4()));
    fill_histo("h_dphi_n1n1", neutralinos[0]->p4().DeltaPhi(neutralinos[1]->p4()));

    return true;
}


} // namespace truth
