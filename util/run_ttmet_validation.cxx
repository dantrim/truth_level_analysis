#include "truth_level_analysis/ttmet_truth_validation.h"

//analysis
#include "truth_level_analysis/utility/ChainHelper.h"

//xAOD
#include "xAODRootAccess/Init.h"

//boost
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

//std/stl
#include <string>
using namespace std;

//ROOT
#include "TChain.h"

int main(int argc, char** argv)
{
    string input_file = "";
    string output_dir = "./";
    int dbg_level = 0;
    string suffix = "";
    int n_entries = -1;

    namespace po = boost::program_options;
    string description = "Run MC tt+MET 2L generation validation";
    po::options_description desc(description);
    desc.add_options()
        ("help,h", "produce help message")
        ("input,i",
            po::value(&input_file),
                "input file [DAOD, text filelist, or directory of files]")
        ("output,o",
            po::value(&output_dir)->default_value("./"),
            "output directory to dump any produced files")
        ("nentries,n",
            po::value(&n_entries)->default_value(-1),
            "provide number of entries to process [default: all]")
        ("debug,d",
            po::value(&dbg_level)->default_value(0),
            "set the debug level [default: 0, lowest level of verbosity]")
        ("suffix",
            po::value(&suffix)->default_value(""),
            "provide an output suffix")
    ; // end input args
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
    po::notify(vm);

    if(vm.count("help")) {
        cerr << desc << endl;
        exit(0);
    }

    if(vm.count("input")) {
        cout << "input file : " << input_file << endl;
    }
    else {
        cerr << "ERROR: you did not provide an input file" << endl;
        exit(1);
    }

    boost::filesystem::path path(output_dir);
    if(boost::filesystem::exists(path)) {
        if(!boost::filesystem::is_directory(path)) {
            cerr << "ERROR: provided output directory (=" << output_dir << ") is not a directory" << endl;
            exit(1);
        }
    }
    else {
        cerr << "ERROR: provided output directory (=" << output_dir << ") does not exist" << endl;
        exit(1);
    }

    xAOD::Init("TtMetTruthValidation");
    TChain* chain = new TChain("CollectionTree");
    int file_err = ChainHelper::addInput(chain, input_file, true);
    if(file_err) return 1;
    Long64_t n_total = chain->GetEntries();
    chain->ls();

    if(n_entries < 0) n_entries = n_total;

    // setup the looper
    truth::TtMetTruthValidation* looper = new truth::TtMetTruthValidation();
    looper->set_debug_level(dbg_level);
    looper->set_output_dir(output_dir);
    looper->set_input_samplename(input_file);
    looper->set_suffix(suffix);
    chain->Process(looper, input_file.c_str(), n_entries);

    // clean up
    delete chain;
    delete looper;

    return 0;
}
