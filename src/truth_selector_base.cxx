#include "truth_level_analysis/truth_selector_base.h"

//std/stl
#include <iostream>
#include <sstream>
using namespace std;

//xAOD
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/tools/ReturnCheck.h"
#include "xAODEventInfo/EventInfo.h"
#include "xAODRootAccess/TStore.h"

//ROOT
#include "TROOT.h"
#include "TChain.h"
#include "TTree.h"

//ClassImp(truth::TruthSelectorBase)

namespace truth
{

TruthSelectorBase::TruthSelectorBase() :
    m_dbg(0),
    n_evt_processed(0),
    n_evt_stored(0),
    m_output_tree_file(nullptr),
    m_output_tree(nullptr),
    m_tree(nullptr)
//    m_event(new xAOD::TEvent(xAOD::TEvent::kClassAccess)
    //m_tstore(new xAOD::TStore),
{
    cout << "TruthSelectorBase()" << endl;
    cout << MYFUNC << " TruthSelector wakes up" << endl;
    m_event = new xAOD::TEvent(xAOD::TEvent::kClassAccess);
    m_tstore = new xAOD::TStore();
}

TruthSelectorBase::~TruthSelectorBase()
{
    delete m_event;
    delete m_tstore;
}

void TruthSelectorBase::Init(TTree* tree)
{
    string fn;
    cout << MYFUNC << " " << m_input_samplename << endl;
    if(tree) {
        cout << MYFUNC << " Found tree with " << tree->GetEntries() << " entries" << endl;
    }
    else {
        cout << MYFUNC << " ERROR input tree is null" << endl;
        exit(1);
    }

    m_tree = tree;
    m_event->readFrom(tree);
    return;
}

void TruthSelectorBase::SlaveBegin(TTree* /*tree*/)
{
    cout << MYFUNC << endl;
    m_timer.Start();
}

void TruthSelectorBase::Begin(TTree* /*tree*/)
{
    cout << MYFUNC << endl;
}

Bool_t TruthSelectorBase::Process(Long64_t entry)
{
    cout << MYFUNC << endl;
    return kTRUE;
}

void TruthSelectorBase::Terminate()
{
    cout << MYFUNC << endl;
    cout << timer_summary() << endl;
}

void TruthSelectorBase::save_output()
{
    cout << MYFUNC << endl;

    m_output_tree_file = m_output_tree->GetCurrentFile();
    if(m_output_tree_file) {
        TDirectory* cd = gROOT->CurrentDirectory();
        m_output_tree_file->cd();
        TNamed n("input_sample", m_input_samplename.c_str());
        n.Write();
        cd->cd();
    }
    m_output_tree->Write(0, TObject::kOverwrite);
    cout << MYFUNC << " Output tree saved to : " << m_output_tree_file->GetName() << endl;
    m_output_tree_file->Close();
}

std::string TruthSelectorBase::timer_summary()
{
    double real_time = m_timer.RealTime();
    double cpu_time = m_timer.CpuTime();
    int hours = int(real_time / 3600);
    real_time -= hours * 3600;
    int min = int(real_time / 60);
    real_time -= min * 60;
    int sec = int(real_time);

    float speed = n_evt_processed / m_timer.RealTime() / 1000;
    TString line1; line1.Form("Real %d:%02d:%02d, CPU %.3f", hours, min, sec, cpu_time);
    TString line2; line2.Form("[kHz] : %2.3f", speed);
    ostringstream oss;
    oss << MYFUNC << " -------------------------------------------------------\n"
        << MYFUNC << "  Number of events processed : " << n_evt_processed << endl
        << MYFUNC << "  Number of events saved     : " << n_evt_stored << endl
        << MYFUNC << "  \t Analysis time        : " << line1 << endl
        << MYFUNC << "  \t Analysis speed       : " << line2 << endl
        << MYFUNC << " -------------------------------------------------------" << endl;
    return oss.str();
}

} // namespace truth
