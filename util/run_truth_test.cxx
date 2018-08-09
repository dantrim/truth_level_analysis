#include "truth_level_analysis/truth_selector_base.h"

//analysis
#include "truth_level_analysis/utility/ChainHelper.h"

//std/stl
#include <iostream>
using namespace std;

//boost
#include <boost/program_options.hpp>

//xAOD
#include "xAODRootAccess/Init.h"

int main(int argc, char** argv)
{
    string dummy = "";
    namespace po = boost::program_options;
    string description = "Test Executable";
    po::options_description desc(description);
    desc.add_options()
        ("help,h", "print out the help message")
        ("dummy", po::value(&dummy)->default_value("default"), "a dummy input variable")
    ;
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
    po::notify(vm);

    if(vm.count("help")) {
        cerr << desc << endl;
        exit(1);
    }

    // init xAOD
    xAOD::Init("TestTruth");

    truth::TruthSelectorBase* looper = new truth::TruthSelectorBase();

    delete looper;
    return 0;
}
