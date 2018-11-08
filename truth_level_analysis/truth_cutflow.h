#ifndef TRUTH_CUTFLOW_H
#define TRUTH_CUTFLOW_H

//analysis
#include "truth_level_analysis/truth_selector_base.h"

//std/stl
#include <map>

//xAOD

//ROOT

namespace truth
{
class TruthCutflow : public TruthSelectorBase
{
    public :
        TruthCutflow();
        virtual ~TruthCutflow(){};

        void set_suffix(std::string val) { m_suffix = val; }
        std::string suffix() { return m_suffix; }

        void set_dsid(int dsid) { m_dsid = dsid; }
        int dsid() { return m_dsid; }

        void set_lumi(float lumi) { m_lumi = lumi; }
        float lumi() { return m_lumi; }

        //misc
        void load_maps();
        std::map<int, double> & sumw_map() { return m_sumw_map; }
        std::map<int, double> & xsec_map() { return m_xsec_map; }

        // TSelector overrides
        void SlaveBegin(TTree* tree);
        void Terminate();
        virtual Bool_t Process(Long64_t entry);

    private :
        double m_mc_weight;
        double m_event_weight;
        int m_dsid;
        double m_lumi;
        std::string m_suffix;

        std::map<int, double> m_sumw_map;
        std::map<int, double> m_xsec_map;
        

}; // class TruthCutflow
} // namespace truth

#endif // TRUTH_CUTFLOW_H
