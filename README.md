# truth_level_analysis
Analysis code for processing the ATLAS TRUTHX data format

# Requirements
Be on a machine with access to `/cvmfs/`.

# Expectation
This code likes to run over **TRUTH1** format. Get out of here with that **TRUTH3** racket.

# Building and Running
```bash
mkdir analysis_dir; cd analysis_dir
mkdir source; cd source
asetup "AnalysisBase,21.2,latest,here"
git clone git@github.com:dantrim/truth_level_analysis.git
cd ..
mkdir build; cd build
cmake ../source
make -j
source x86_64*/setup.sh
run_truth_test
```
In reality, replace the executable `run_truth_test` with the one you want.

# Structure
The C++ header (re:**\*.h**) files are located in the directory *truth_level_analysis/*.
The C++ source (re:**\*.cxx**) files are located in *src/*. Main applications/executables
are locatedin *util/*. The `CMakeLists.txt` file will pick up any file in these directories
and will build an executable out of any main application in *util/*, there is no need to
explicitly add it.

# Current Goings On

Here is a listing of the current analysis loopers that can process input TRUTH DAODs:

|Executable| Looper Class (source) | Purpose|
|----------|-----------------------|--------|
|run_hh_validation | HHTruthValidation (hh_truth_validation) | Process HH signal samples and make generator/filter validation histograms |
|run_ttmet_validation| TtMetTruthValidation (ttmet_truth_vadliation) | Process 2L tt+MET signal samples and make generator/filter validation histograms |
