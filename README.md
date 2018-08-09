# truth_level_analysis
Analysis code for processing the ATLAS TRUTHX data format

# Requirements
Be on a machine with access to `/cvmfs/`.

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
