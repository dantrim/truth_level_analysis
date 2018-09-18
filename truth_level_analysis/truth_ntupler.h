#ifndef TRUTH_NTUPLER_H
#define TRUTH_NTUPLER_H

//analysis
#include "truth_level_analysis/truth_selector_base.h"

//std/stl
#include <map>

//ROOT
class TBranch;

namespace truth
{

class TruthNtupler : public TruthSelectorBase
{
    public :
        TruthNtupler();
        virtual ~TruthNtupler(){};

        void set_do_sumw(bool doit) { m_do_sumw = doit; }
        bool do_sumw() { return m_do_sumw; }

        void set_dsid(int dsid) { m_dsid = dsid; }
        int dsid() { return m_dsid; }

        void set_lumi(float lumi) { m_lumi = lumi; }
        float lumi() { return m_lumi; }

        // overrides of TSelector
        void SlaveBegin(TTree* tree);
        void Terminate();
        virtual Bool_t Process(Long64_t entry);

        // misc
        void setup_output();
        std::map<int, double> & sumw_map() { return m_sumw_map; }
        std::map<int, double> & xsec_map() { return m_xsec_map; }

        // flow
        void load_sumw_and_xsec();
        bool update_sumw();
        bool process_event();
        void reset_ttree_vars();
        void var(std::string, int val);
        void var(std::string, float val);

        // analysis
        double w() { return ((m_mc_weight * m_xsec * (m_lumi * 1000.)) / m_sumw); }

    private :

        bool m_do_sumw;

        int m_dsid;
        double m_xsec;
        double m_lumi; // in [1/fb]
        double m_sumw;

        std::map<int, double> m_sumw_map;
        std::map<int, double> m_xsec_map;

        // OUTPUT NTUPLE
        std::map<std::string, float> m_float_vm;
        std::map<std::string, int> m_int_vm;

}; // class TruthNtupler

} // namespace truth

#endif
