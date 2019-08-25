#include "windowtosnp.h"

// Note: if I ever screen out SNPs (ones that have only one allele including adm) then keep separate windowLocs that stores where they are to save time

void WindowToSnpCalls(ProcessedInput* processedInput, int** windowCalls, int emIterationIndex, double** marginals, int emIteration){
  int numWindows = processedInput->numWindows;
  int numAdm = processedInput->numAdm;  // # of admixed haplotypes
  int numAncPops = processedInput->numAncPops;
  vector<int> * windowBeginIndexes = processedInput->windowBeginIndexes;
  vector<int> * windowEndIndexes = processedInput->windowEndIndexes;
  string outputName = processedInput->outputName;
  bool doForwardBackward = processedInput->doForwardBackward;
  //vector<int> * hapClasses = processedInput->hapClasses;
  
  ofstream outfile;
  string emIterationIndexString;
  stringstream emStream;
  emStream << emIterationIndex;
  emIterationIndexString = emStream.str();
  outfile.open((outputName + "." + emIterationIndexString + ".Viterbi.txt").c_str());
  
  ofstream margOutfile;
  if(doForwardBackward){
    margOutfile.open((outputName + "." + emIterationIndexString + ".ForwardBackward.txt").c_str());
  }
  
  if(processedInput->compressedOutput){
    ofstream numSnpsPerWindowFile;
    numSnpsPerWindowFile.open((outputName + "." + emIterationIndexString + ".SNPsPerWindow.txt").c_str());
    for(int windex = 0; windex < numWindows; windex++){
      numSnpsPerWindowFile << (*windowEndIndexes)[windex] - (*windowBeginIndexes)[windex] + 1 << endl;
      for(int admIndex = 0; admIndex < numAdm; admIndex++){
        int windowCall = windowCalls[admIndex][windex];
        outfile << windowCall+1 << " ";  // Add 1 so agrees with user input
        if(doForwardBackward){
          for(int ancIndex = 0; ancIndex < numAncPops; ancIndex++){
            margOutfile << marginals[admIndex][windex*numAncPops + ancIndex] << " ";
          }
        }
      }
      outfile << endl;
      if(doForwardBackward){
        margOutfile << endl;
      }
    }
  } else{
    // Write calls horizontally by individual, each column is individual haplotype
    // For marginals, with K populations, each group of K columns is for a haplotype
    for(int windex = 0; windex < numWindows; windex++){
      int numSnpsInWindow = (*windowEndIndexes)[windex] - (*windowBeginIndexes)[windex] + 1;
      for(int snpIndex = 0; snpIndex < numSnpsInWindow; snpIndex++){
        for(int admIndex = 0; admIndex < numAdm; admIndex++){
          int windowCall = windowCalls[admIndex][windex];
          outfile << windowCall+1 << " ";  // Add 1 so agrees with user input
          if(doForwardBackward){
            for(int ancIndex = 0; ancIndex < numAncPops; ancIndex++){
              margOutfile << marginals[admIndex][windex*numAncPops + ancIndex] << " ";
            }
          }
        }
        outfile << endl;
        if(doForwardBackward){
          margOutfile << endl;
        }
      }
    }
  }
  
  outfile.close();
  if(doForwardBackward){
    margOutfile.close();
  }
}
