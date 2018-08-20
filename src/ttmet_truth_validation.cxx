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
    m_outfilename("")
{
}

void TtMetTruthValidation::SlaveBegin(TTree* /*tree*/)
{
    timer()->Start();
}

void TtMetTruthValidation::Terminate()
{
    TruthSelectorBase::Terminate();
}

Bool_t TtMetTruthValidation::Process(Long64_t entry)
{
    static Long64_t chain_entry = -1;
    chain_entry++;
    event()->getEntry(chain_entry);
    n_evt_processed++;

    if(chain_entry % 500  == 0) {
        cout << MYFUNC << " Processing entry " << chain_entry << endl;
    }

    const xAOD::EventInfo* ei = 0;
    RETURN_CHECK(GetName(), event()->retrieve(ei, "EventInfo"));
    mc_weight(ei->mcEventWeight()); // this propagates the event weight value to later calls


    return true;
}


} // namespace truth
