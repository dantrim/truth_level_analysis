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

    private :

        double m_mc_weight;

        std::string m_outfilename;

        // ROOT
        TFile* m_rfile;

}; // class

} // namespace
#endif
