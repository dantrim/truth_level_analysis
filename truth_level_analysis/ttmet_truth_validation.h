#ifndef TTMET_TRUTH_VALIDATION_H
#define TTMET_TRUTH_VALIDATION_H

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

class TtMetTruthValidation : public TruthSelectorBase
{
    public :

        TtMetTruthValidation();
        virtual ~TtMetTruthValidation(){};

        // override the inherited TSelector methods
        void SlaveBegin(TTree* tree);
        void Terminate();
        virtual Bool_t Process(Long64_t entry);

        void setup_output();
        void setup_histograms();

        bool process_stop_decays();
        bool process_final_state();

    private :

        double m_mc_weight;
        std::string m_outfilename;
        TFile* m_rfile;


}; // class TtMetTruthValidation


} // namespace truth

#endif
