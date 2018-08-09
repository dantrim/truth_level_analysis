#ifndef TRUTH_SELECTOR_BASE_H
#define TRUTH_SELECTOR_BASE_H

//ROOT
#include "TSelector.h"
#include "TFile.h"
#include "TStopwatch.h"

//std/stl
#include <string>
#include <vector>

//analysis
#include "truth_level_analysis/helpers.h"

//xAOD
#include "xAODRootAccess/tools/ReturnCheck.h"
#include "xAODEventInfo/EventInfo.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODTruth/TruthEventContainer.h"

#undef CHECK
#define CHECK( ARG ) \
    do { \
        const bool result = ARG; \
        if( ! result ) { \
            ::Error("TruthSelector", "Failed to execute: \"%s\"", \
                #ARG ); \
            exit(1); \
        } \
    } while (false)

namespace truth
{

class TruthSelectorBase : public TSelector
{
    public :

        TruthSelectorBase();
        virtual ~TruthSelectorBase();

        // inherited TSelector methods
        virtual Int_t Version() const { return 2; } // ROOT magic numbers
        virtual void Init(TTree* tree);
        virtual Bool_t Notify() { return kTRUE; }
        virtual void Begin(TTree* tree);
        virtual void SlaveBegin(TTree* tree);
        virtual void Terminate();
        virtual Bool_t Process(Long64_t entry);

        void set_debug_level(int dbg_level) { m_dbg = dbg_level; }
        int debug_level() { return m_dbg; }

        void set_input_samplename(std::string name) { m_input_samplename = name; }
        void set_output_filename(std::string name) { m_output_filename = name; }
        std::string input_samplename() { return m_input_samplename; }
        std::string output_filename() { return m_output_filename; } 

        void save_output();

        TStopwatch* timer() { return &m_timer; }
        std::string timer_summary();

        //xAOD
        xAOD::TEvent* event() { return m_event; }
        xAOD::TStore* store() { return m_tstore; }

        //ROOT
        TFile* output_file() { return m_output_tree_file; }
        TTree* output_tree() { return m_output_tree; }

    protected :

        int m_dbg;

        int n_evt_processed;
        int n_evt_stored;

        TFile* m_output_tree_file;
        TTree* m_output_tree;

        TTree* m_tree;
        xAOD::TEvent* m_event;
        xAOD::TStore* m_tstore;
        TStopwatch m_timer;

        std::string m_input_samplename;
        std::string m_output_filename;

    //ClassDef(TruthSelectorBase, 0);

}; // class

} // namespace truth

#endif
