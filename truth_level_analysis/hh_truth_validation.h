#ifndef HH_TRUTH_VALIDATION_H
#define HH_TRUTH_VALIDATION_H

//analysis
#include "truth_level_analysis/truth_selector_base.h"

//xAOD
#include "xAODMissingET/MissingETContainer.h"
#include "xAODTruth/TruthParticleContainer.h"

//std/stl

//ROOT
class TH1F;
class TFile;


namespace truth
{

class HHTruthValidation : public TruthSelectorBase
{
    public :

        struct HiggsCandidate
        {
            HiggsCandidate() :
                parent(nullptr),
                child0(nullptr),
                child1(nullptr) {}
            const xAOD::TruthParticle* parent;
            const xAOD::TruthParticle* child0;
            const xAOD::TruthParticle* child1;
            bool both_children() { return (child0 && child1); }
        };

        struct WCandidate
        {
            WCandidate() :
                w(nullptr),
                lepton(nullptr),
                neutrino(nullptr) {}
            const xAOD::TruthParticle* w;
            const xAOD::TruthParticle* lepton;
            const xAOD::TruthParticle* neutrino;
            bool both_children() { return (lepton && neutrino); }
        };

        HHTruthValidation();
        virtual ~HHTruthValidation(){};

        // set the suffix to any output
        void set_suffix(std::string val) { m_suffix = val; }
        std::string suffix() { return m_suffix; }

        // setup our output file and histograms
        void setup_output();
        void setup_histograms();

        // override the TSelector methods 
        void SlaveBegin(TTree* tree);
        void Terminate();
        virtual Bool_t Process(Long64_t entry);

        bool process_strict_truth();
        bool get_w_decay(const xAOD::TruthParticle* wboson, WCandidate& w);

        bool process_final_state();

    private :

        double m_mc_weight;

        std::string m_outfilename;

        // ROOT
        TFile* m_rfile;

}; // class

} // namespace
#endif
