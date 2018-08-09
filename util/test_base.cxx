#include "truth_level_analysis/truth_selector_base.h"

//std/stl
#include <iostream>
using namespace std;

int main(int argc, char** argv)
{
    cout << "hello world" << endl;

    truth::TruthSelectorBase* base = new truth::TruthSelectorBase();
    delete base;

}
