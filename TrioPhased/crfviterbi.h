#ifndef CRFVITERBI
#define CRFVITERBI

#include "getdata.h"
#include <math.h>

using namespace std;

int** CrfViterbi(ProcessedInput* proessedInput, double ** rfProbs, int** windowCalls, double** marginals);

#endif
