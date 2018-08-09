#include "truth_level_analysis/hh_truth_validation.h"

//analysis
#include "truth_level_analysis/utility/ChainHelper.h"

//xAOD
 #include "xAODRootAccess/Init.h"

//boost
#include <boost/program_options.hpp>

//std/stl
#include <string>
using namespace std;

//ROOT
#include "TChain.h"

void help()
{
    cout << "Options:" << endl;
    cout << " -i|--input            input file list, file, or dir" << endl;
    cout << " -n|--nentries         number of entries to process [default: all]" << endl;
    cout << " -d|--debug            set the debug level [default: 0]" << endl;
    cout << " --suffix              suffix to append to any output files" << endl;
    cout << " -h|--help             print this help message" << endl;
}

int main(int argc, char** argv)
{
    string input_file = "";
    int dbg_level = 0;
    string suffix = "";
    int n_entries = -1;
    
    namespace po = boost::program_options;
    string description = "Run MC HH generation validation";
    po::options_description desc(description);
    desc.add_options()
        ("help,h", "produce help message")
        ("input,i",
            po::value(&input_file),
            "input file [DAOD, text file, or directory of files]")
        ("nentries,n",
            po::value(&n_entries)->default_value(-1),
            "provide number of entries to process")
        ("debug,d",
            po::value(&dbg_level)->default_value(0),
            "set the debug level")
        ("suffix",
            po::value(&suffix)->default_value(""),
            "provide an output suffix")
    ;
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
    po::notify(vm);

    if(vm.count("help")) {
        cout << desc << endl;
//        help();
        exit(0);
    }

    if(vm.count("input")) {
        cout << "input file : " << input_file << endl;
    }
    else {
        cout << "ERROR: you did not provide an input file" << endl;
        exit(1);
    }

    xAOD::Init("HHTruthValidation");
    TChain* chain = new TChain("CollectionTree");
    int file_err = ChainHelper::addInput(chain, input_file, true);
    if(file_err) return 1;
    Long64_t n_total = chain->GetEntries();
    chain->ls();

    if(n_entries < 0) n_entries = n_total;

    // setup the looper
    truth::HHTruthValidation* looper = new truth::HHTruthValidation();
    looper->set_debug_level(dbg_level);
    looper->set_input_samplename(input_file);
    looper->set_suffix(suffix);
    chain->Process(looper, input_file.c_str(), n_entries);

    // clean up our old news
    delete chain;
    delete looper;
    return 0;
}
