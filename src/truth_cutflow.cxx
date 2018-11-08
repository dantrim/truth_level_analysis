//analysis
#include "truth_level_analysis/truth_cutflow.h"
#include "truth_level_analysis/helpers.h"

//std/stl
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
using namespace std;

//xAOD
#include "xAODEventInfo/EventInfo.h"
#include "xAODJet/Jet.h"
#include "xAODJet/JetContainer.h"


//EDM/tools
#include "PathResolver/PathResolver.h"

namespace truth
{
//_____________________________________________________________________________
TruthCutflow::TruthCutflow() :
    m_mc_weight(1.0),
    m_event_weight(1.0),
    m_dsid(0),
    m_lumi(1.),
    m_suffix("")
{
    cout << MYFUNC << endl;
}
//_____________________________________________________________________________
void TruthCutflow::SlaveBegin(TTree* /*tree*/)
{
    timer()->Start();
}
//_____________________________________________________________________________
void TruthCutflow::Terminate()
{
    TruthSelectorBase::Terminate();
}
//_____________________________________________________________________________
void TruthCutflow::load_maps()
{
    string xsec_dir_search = "truth_level_analysis/cross_section/";
    string sumw_dir_search = "truth_level_analysis/sumw/";

    string xsec_dir = PathResolverFindCalibDirectory(xsec_dir_search);
    string sumw_dir = PathResolverFindCalibDirectory(sumw_dir_search);

    if(xsec_dir == "") {
        stringstream err;
        err << MYFUNC << " ERROR: unable to locate xsec directory (=" << xsec_dir_search << ")";
        throw  std::runtime_error(err.str());
    }

    if(sumw_dir == "") {
        stringstream err;
        err << MYFUNC << " ERROR: unable to locate sumw directory (=" << sumw_dir_search << ")";
        throw std::runtime_error(err.str());
    }

}
//_____________________________________________________________________________
Bool_t TruthCutflow::Process(Long64_t /*entry*/)
{
    static Long64_t chain_entry = -1;
    chain_entry++;
    event()->getEntry(chain_entry);
    n_evt_processed++;

    if(chain_entry == 0) {
        load_maps();
    }

    if(chain_entry % 1000 == 0) {
        cout << MYFUNC << " Processing entry " << chain_entry << endl;
    }

    const xAOD::EventInfo* ei = 0;
    RETURN_CHECK(GetName(), event()->retrieve(ei, "EventInfo"));
    mc_weight(ei->mcEventWeight());

    return true;
}

} // namespace truth

