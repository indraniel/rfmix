#include "crfviterbi.h"
#include <float.h>

// Now with forward-backward goodness!
int** CrfViterbi(ProcessedInput* processedInput, double ** rfProbs, int** windowCalls, double** marginals, int emIteration){
  int numWindows = processedInput->numWindows;
  int numAncPops = processedInput->numAncPops;  // # of ancestral populations
  int numAdm = processedInput->numAdm;  // # of admixed haplotypes
  int genSinceAdmixture = processedInput->genSinceAdmixture;
  vector<double> * windowBeginLocs = processedInput->windowBeginLocs;
  vector<double> * windowEndLocs = processedInput->windowEndLocs;
  vector<int> * windowBeginIndexes = processedInput->windowBeginIndexes;
  vector<int> * windowEndIndexes = processedInput->windowEndIndexes;
  int** numHetSitesPerInd = processedInput->numHetSitesPerInd;
  int* numAdmPhasings = processedInput->numAdmPhasings;
  vector<bool>* admWindows = processedInput->admWindows;
  vector<bool>* admWindowsPhasings = processedInput->admWindowsPhasings;
  vector<int> * hapClasses = processedInput->hapClasses;
  
  int numAdmInds = numAdm / 2;
  double nonLogSwitchProb = 0.07;
  double switchProb = log(nonLogSwitchProb);
  double noSwitchProb = log(1.0 - nonLogSwitchProb);
  
  // Calculate probabilites of recombination between each window
  double distBtWindowsMids[numWindows-1];
  for(int windex=0; windex < numWindows-1; windex++){
    double nextWindowHalf = ((*windowEndLocs)[windex+1] - (*windowBeginLocs)[windex+1])/2.0;
    double thisWindowHalf = ((*windowBeginLocs)[windex] - (*windowBeginLocs)[windex])/2.0;
    double distBtWindowSegs = ((*windowBeginLocs)[windex+1] - (*windowEndLocs)[windex]);
    distBtWindowsMids[windex] = nextWindowHalf + thisWindowHalf + distBtWindowSegs;
  }
  
  double transitionProbs[numWindows-1];
  double genSinceAdmixtureDouble = (double) genSinceAdmixture;
  for(int windex=0; windex < numWindows-1; windex++){
    transitionProbs[windex] = 1-exp(-distBtWindowsMids[windex]*(genSinceAdmixtureDouble)/100.0);  // divide by 100 bc converting from cM to M
  }
  
  // Perform Viterbi
  // The number of states differs between windows and individuals
  // # states in a window for a person = 2 x numWindowPhasings x numAncestryPops x numAncestryPops
  // We do Viterbi for each individual, not haplotype
  int numStatesPerWindow[numWindows];
  int numPhasingsSeen[numWindows];  // for indexing into things like rfProb
  for(int windex = 0; windex < numWindows; windex++){
    numPhasingsSeen[windex] = 0;
  }
  for(int admInd = 0; admInd < numAdmInds; admInd++){
    cout << "CRF on individual: " << admInd << endl;
    int maxStates = 0;
    for(int windex = 0; windex < numWindows; windex++){
      int indNumPhasings = numHetSitesPerInd[windex][admInd] + 1;
      int numStates = 2 * indNumPhasings * numAncPops * numAncPops;
      numStatesPerWindow[windex] = numStates;
      if(numStates > maxStates){
        maxStates = numStates;
      }
    }
    
    // points back to best previous state reaching the current state
    int ** stateBestPaths = new int*[numWindows - 1];
    for(int windex = 0; windex < numWindows - 1; windex++){
      stateBestPaths[windex] = new int[numStatesPerWindow[windex + 1]];    // make sure windows align
    }
    
    // stores sum of logProbs of best path that reaches this state at the current time
    // for simplicity just make it the size of the maximum number of states across all windows
    // ordering of states: same as ordering of rfProbs for one individual
    // i.e. phasing1chrom1anc1
    double cumulativeProbs[maxStates];
    double cumulativeProbsHelper[maxStates];
    
    // Handle initial window
    int window0NumHetSites = numHetSitesPerInd[0][admInd];
    int window0IndNumPhasings = window0NumHetSites + 1;
    // for phasing == 0
    double window0LogNoSwitchProb = noSwitchProb * ((double)window0NumHetSites);    // sums because logs
    // for phasing != 0
    double window0LogSwitchProb = switchProb + noSwitchProb * ((double)(window0NumHetSites - 1));
    
    for(int phasing = 0; phasing < window0IndNumPhasings; phasing++){
      for(int chrom = 0; chrom < 2; chrom++){  // if chrom==1 then flip the phasing so strands flipped
        for(int ancestryTop = 0; ancestryTop < numAncPops; ancestryTop++){
          for(int ancestryBottom = 0; ancestryBottom < numAncPops; ancestryBottom++){
            int state = phasing*2*numAncPops*numAncPops + chrom*numAncPops*numAncPops + ancestryTop*numAncPops + ancestryBottom;
                        
            double topRFprob = rfProbs[0][numPhasingsSeen[0]*2*numAncPops + phasing*2*numAncPops + chrom*numAncPops + ancestryTop];  // top strand for current orientation
            double bottomRFprob = rfProbs[0][numPhasingsSeen[0]*2*numAncPops + phasing*2*numAncPops + (1 - chrom)*numAncPops + ancestryBottom];  // i.e. chrom
            
            if(phasing == 0){
              cumulativeProbs[state] = log(topRFprob) + log(bottomRFprob) + window0LogNoSwitchProb;
            }else{
              cumulativeProbs[state] = log(topRFprob) + log(bottomRFprob) + window0LogSwitchProb;
            }
          }
        }
      }
    }
    
    // Do remaining windows
    int intLastNumPhasings = window0IndNumPhasings;
    for(int windex = 0; windex < numWindows-1; windex++){  // starts at 0 because recording state we came from
      int numHetSites = numHetSitesPerInd[windex + 1][admInd];
      int indNumPhasings = numHetSites + 1;
      //int numStates = 2 * indNumPhasings * numAncPops * numAncPops;
      double logNoSwitchProb = noSwitchProb * ((double)numHetSites);    // sums because logs
      double logSwitchProb = switchProb + noSwitchProb * ((double)(numHetSites - 1));
      
      // thisState translates to thisPhasingthisChromthisAncestry combo
      for(int thisPhasing = 0; thisPhasing < indNumPhasings; thisPhasing++){
        for(int thisChrom = 0; thisChrom < 2; thisChrom++){
          for(int thisAncestryTop = 0; thisAncestryTop < numAncPops; thisAncestryTop++){
            for(int thisAncestryBottom = 0; thisAncestryBottom < numAncPops; thisAncestryBottom++){
              int thisState = thisPhasing*2*numAncPops*numAncPops + thisChrom*numAncPops*numAncPops + thisAncestryTop*numAncPops + thisAncestryBottom;
              
              double thisTopRFprob = rfProbs[windex + 1][numPhasingsSeen[windex + 1]*2*numAncPops + thisPhasing*2*numAncPops + thisChrom*numAncPops + thisAncestryTop];  // top strand for current orientation
              double thisBottomRFprob = rfProbs[windex + 1][numPhasingsSeen[windex + 1]*2*numAncPops + thisPhasing*2*numAncPops + (1 - thisChrom)*numAncPops + thisAncestryBottom];  // i.e. chrom
              
              double maxScore = -DBL_MAX;
              int lastMaxState = 0;
              
              double thisStateProb = 0.0;
              if(thisPhasing == 0){
                thisStateProb = log(thisTopRFprob) + log(thisBottomRFprob) + logNoSwitchProb;
              }else{
                thisStateProb = log(thisTopRFprob) + log(thisBottomRFprob) + logSwitchProb;
              }
              
              bool thisLeftPurpleTop;
              if(thisChrom == 0){
                thisLeftPurpleTop = true;
              }else{
                thisLeftPurpleTop = false;
              }
              
              for(int lastPhasing = 0; lastPhasing < intLastNumPhasings; lastPhasing++){
                for(int lastChrom = 0; lastChrom < 2; lastChrom++){
                  for(int lastAncestryTop = 0; lastAncestryTop < numAncPops; lastAncestryTop++){
                    for(int lastAncestryBottom = 0; lastAncestryBottom < numAncPops; lastAncestryBottom++){
                      // only transition bw certain states
                      // match colors along top bc top matches iff bottom matches
                      bool lastRightPurpleTop;
                      if(lastChrom == 0 && lastPhasing == 0){
                        lastRightPurpleTop = true;
                      } else if (lastChrom == 0 && lastPhasing != 0){
                        lastRightPurpleTop = true;
                      } else{
                        lastRightPurpleTop = false;
                      }
                                            
                      int lastState = lastPhasing*2*numAncPops*numAncPops + lastChrom*numAncPops*numAncPops + lastAncestryTop*numAncPops + lastAncestryBottom;
                      double score = cumulativeProbs[lastState] + thisStateProb;
                      
                      if(lastRightPurpleTop != thisLeftPurpleTop){
                        // mismatch in phasings so skip
                        //continue;
                        // actually, need to consider a switch between windows
                        score = score - noSwitchProb + switchProb;
                      }
                      
                      // add into score the consistency of ancestries on BOTH strands
                      if(thisAncestryTop == lastAncestryTop){
                        score += log((1 - transitionProbs[windex])*(1.0/((double)numAncPops)) + transitionProbs[windex]*(1.0/((double)numAncPops))*(1.0/((double)numAncPops)));
                      } else{
                        score += log(transitionProbs[windex]*(1.0/((double)numAncPops))*(1.0/((double)numAncPops)));
                      }
                      
                      if(thisAncestryBottom == lastAncestryBottom){
                        score += log((1 - transitionProbs[windex])*(1.0/((double)numAncPops)) + transitionProbs[windex]*(1.0/((double)numAncPops))*(1.0/((double)numAncPops)));
                      }else{
                        score += log(transitionProbs[windex]*(1.0/((double)numAncPops))*(1.0/((double)numAncPops)));
                      }
                      
                      if (score > maxScore){
                        maxScore = score;
                        lastMaxState = lastState;
                      }

                    }
                  }
                }
              } // end of looking through last states
              
              stateBestPaths[windex][thisState] = lastMaxState;
              cumulativeProbsHelper[thisState] = maxScore;
              
            }
          }
        }
      }  // end of looking through this states
      
      for(int stateIndex = 0; stateIndex < maxStates; stateIndex++){
        cumulativeProbs[stateIndex] = cumulativeProbsHelper[stateIndex];
      }
      
      intLastNumPhasings = indNumPhasings;
    }
    
    // Get max cumulative probs and backtrace
    double maxScore = -DBL_MAX;
    int maxIndex = 0;
    int lastNumStates = numStatesPerWindow[numWindows - 1];
    for(int stateIndex = 0; stateIndex < lastNumStates; stateIndex++){
      double score = cumulativeProbs[stateIndex];
      if(score > maxScore){
        maxScore = score;
        maxIndex = stateIndex;
      }
    }
    
    // At each state get ancestryTop and ancestryBottom
    // ancestryTop goes to admIndex and ancestryBottom goes to admIndex+1
    int pointedTo = maxIndex;
    
    // Store the MAP path which is used in getting ancestry marginals with forward-backward
    int vitChrs[numWindows];
    int phasingChrs[numWindows];
    
    // Get ancestryTop and ancestryBottom from state for this window
    int ancBottom = pointedTo % numAncPops;
    int ancTop = (pointedTo % (numAncPops*numAncPops)) / numAncPops;
    windowCalls[2*admInd][numWindows - 1] = ancTop;
    windowCalls[2*admInd+1][numWindows - 1] = ancBottom;
    
    int chr = (pointedTo % (2*numAncPops*numAncPops)) / (numAncPops*numAncPops);
    vitChrs[numWindows-1] = chr;
    int phasing = (pointedTo % (numAdmPhasings[numWindows - 1]*2*numAncPops*numAncPops)) / (2*numAncPops*numAncPops);
    phasingChrs[numWindows-1] = phasing;
    
    int numSNPsInWindow = (*windowEndIndexes)[numWindows - 1] - (*windowBeginIndexes)[numWindows - 1] + 1;
    
    // Correct phasing
    if((*hapClasses)[2*admInd] == 0 || emIteration > 0){  // only rephase ORIGINAL admixed guys if emIteration==0
      for(int snpIndex = 0; snpIndex < numSNPsInWindow; snpIndex++){
        int spot1 = snpIndex*numAdm + admInd*2;
        int spot2 = spot1+1;
        
        admWindows[numWindows - 1][spot1] = admWindowsPhasings[numWindows - 1][numAdmPhasings[numWindows-1]*2*snpIndex + numPhasingsSeen[numWindows-1]*2 + phasing*2 + chr];
        admWindows[numWindows - 1][spot2] = admWindowsPhasings[numWindows - 1][numAdmPhasings[numWindows-1]*2*snpIndex + numPhasingsSeen[numWindows-1]*2 + phasing*2 + (1-chr)];
      }
    }
    
    for(int windex = numWindows-2; windex >=0; windex--){
      pointedTo = stateBestPaths[windex][pointedTo];    // sbp has length numwindows-1
      ancBottom = pointedTo % numAncPops;
      ancTop = (pointedTo % (numAncPops*numAncPops)) / numAncPops;
      windowCalls[2*admInd][windex] = ancTop;
      windowCalls[2*admInd+1][windex] = ancBottom;
      
      if((*hapClasses)[2*admInd] == 0 || emIteration > 0){  // only rephase ORIGINAL admixed guys
        chr = (pointedTo % (2*numAncPops*numAncPops)) / (numAncPops*numAncPops);
        vitChrs[windex] = chr;
        phasing = (pointedTo % (numAdmPhasings[windex]*2*numAncPops*numAncPops)) / (2*numAncPops*numAncPops);
        phasingChrs[windex] = phasing;
        
        numSNPsInWindow = (*windowEndIndexes)[windex] - (*windowBeginIndexes)[windex] + 1;
        
        for(int snpIndex = 0; snpIndex < numSNPsInWindow; snpIndex++){
          int spot1 = snpIndex*numAdm + admInd*2;
          int spot2 = spot1+1;
          
          admWindows[windex][spot1] = admWindowsPhasings[windex][numAdmPhasings[windex]*2*snpIndex + numPhasingsSeen[windex]*2 + phasing*2 + chr];
          admWindows[windex][spot2] = admWindowsPhasings[windex][numAdmPhasings[windex]*2*snpIndex + numPhasingsSeen[windex]*2 + phasing*2 + (1-chr)];
        }
      }
      
    }
    
    
    
    
    
    
    
    
    
    
    
    
    // ==== Forward Backward ====
    if(processedInput->doForwardBackward && ((*hapClasses)[2*admInd] == 0 || emIteration > 0)){
      for(int hapIndex = 0; hapIndex < 2; hapIndex++){
      
        // Get initial beliefs
        double initialProbs[numWindows-1][numAncPops][numAncPops]; // numWindows-1 clusters. Each is K x K matrix, row index is first variable values
        
        for(int windex = 0; windex < numWindows-1; windex++){
          // Enter transition probs
          for(int i = 0; i < numAncPops; i++){
            for(int j = 0; j < numAncPops; j++){
              if(i == j){
                initialProbs[windex][i][j] = (1.0 - transitionProbs[windex])*(1.0/((double)numAncPops)) + transitionProbs[windex]*(1.0/((double)numAncPops))*(1.0/((double)numAncPops));
              }else{
                initialProbs[windex][i][j] = transitionProbs[windex]*(1.0/((double)numAncPops))*(1.0/((double)numAncPops));
              }
              initialProbs[windex][i][j] *= rfProbs[windex][numPhasingsSeen[windex]*2*numAncPops + phasingChrs[windex]*2*numAncPops + abs(vitChrs[windex] - hapIndex)*numAncPops + i];   // rfProbs[admIndex][(windex)*numAncPops + i];
              if(windex == numWindows-2){    // last cluster so multiply in last rfProb variable
                initialProbs[windex][i][j] *= rfProbs[windex+1][numPhasingsSeen[windex+1]*2*numAncPops + phasingChrs[windex+1]*2*numAncPops + abs(vitChrs[windex+1] - hapIndex)*numAncPops + j]; // rfProbs[admIndex][(windex+1)*numAncPops + j];
              }
            }
          }
        }
        
        // Normalize
        //    for(int windex = 0; windex < numWindows-1; windex++){
        //      double sum = 0.0;
        //      for(int i = 0; i < numAncPops; i++){
        //        for(int j = 0; j < numAncPops; j++){
        //          sum += initialProbs[windex][i][j];
        //        }
        //      }
        //      for(int i = 0; i < numAncPops; i++){
        //        for(int j = 0; j < numAncPops; j++){
        //          initialProbs[windex][i][j] /= sum;
        //        }
        //      }
        //    }
        
        // Use initial beliefs to create messages. First is message passed to second cluster from first cluster.
        // Sum across rows to eliminate first variable
        // numWindows-1 clusters, so number of messages is numWindows-2
        // index of message is the index of cluster it multiplies
        double leftMessages[numWindows-2+1][numAncPops];  // added one so aligns
        for(int windex = 1; windex <= numWindows-2; windex++){
          for(int col = 0; col < numAncPops; col++){
            leftMessages[windex][col] = 0.0;
            for(int row = 0; row < numAncPops; row++){
              if(windex != 1){
                leftMessages[windex][col] += initialProbs[windex-1][row][col] * leftMessages[windex-1][row];
              }else{
                leftMessages[windex][col] += initialProbs[windex-1][row][col];
              }
            }
          }
          // Normalize
          double leftMessageSum = 0.0;
          for(int col = 0; col < numAncPops; col++){
            leftMessageSum += leftMessages[windex][col];
          }
          for(int col = 0; col < numAncPops; col++){
            leftMessages[windex][col] /= leftMessageSum;
          }
        }
        
        double rightMessages[numWindows-2][numAncPops];
        for(int windex = numWindows-3; windex >= 0; windex--){
          for(int row = 0; row < numAncPops; row++){
            rightMessages[windex][row] = 0.0;
            for(int col = 0; col < numAncPops; col++){
              if(windex != numWindows-3){
                rightMessages[windex][row] += initialProbs[windex+1][row][col] * rightMessages[windex+1][col];
              }else{
                rightMessages[windex][row] += initialProbs[windex+1][row][col];
              }
            }
          }
          // Normalize
          double rightMessageSum = 0.0;
          for(int col = 0; col < numAncPops; col++){
            rightMessageSum += rightMessages[windex][col];
          }
          for(int col = 0; col < numAncPops; col++){
            rightMessages[windex][col] /= rightMessageSum;
          }
        }
        
        // Calculate final beliefs over clusters
        double finalBeliefs[numWindows-1][numAncPops][numAncPops];
        for(int windex = 0; windex < numWindows-1; windex++){
          if(windex != 0 && windex != numWindows-2){
            for(int row = 0; row < numAncPops; row++){
              for(int col = 0; col < numAncPops; col++){
                finalBeliefs[windex][row][col] = initialProbs[windex][row][col] * leftMessages[windex][row] * rightMessages[windex][col];
              }
            }
          }else if(windex == 0){
            for(int row = 0; row < numAncPops; row++){
              for(int col = 0; col < numAncPops; col++){
                finalBeliefs[windex][row][col] = initialProbs[windex][row][col] * rightMessages[windex][col];
              }
            }
          }else{
            for(int row = 0; row < numAncPops; row++){
              for(int col = 0; col < numAncPops; col++){
                finalBeliefs[windex][row][col] = initialProbs[windex][row][col] * leftMessages[windex][row];
              }
            }
          }
        }
        
        // Normalize
        for(int windex = 0; windex < numWindows-1; windex++){
          double sum = 0.0;
          for(int row = 0; row < numAncPops; row++){
            for(int col = 0; col < numAncPops; col++){
              sum += finalBeliefs[windex][row][col];
            }
          }
          for(int row = 0; row < numAncPops; row++){
            for(int col = 0; col < numAncPops; col++){
              finalBeliefs[windex][row][col] /= sum;
            }
          }
        }
        
        // DEBUG
        //    cout << "Individual " << admIndex << endl;
        //    for(int windex = 0; windex < numWindows - 1; windex++){
        //      cout << "Window " << windex << endl;
        //      cout << "Var A Marginals ";
        //      for(int row = 0; row < numAncPops; row++){
        //        double varA = 0.0;
        //        for(int col = 0; col < numAncPops; col++){
        //          varA += finalBeliefs[windex][row][col];
        //        }
        //        cout << varA << " ";
        //      }
        //      cout << endl;
        //      cout << "Var B Marginals ";
        //      for(int col = 0; col < numAncPops; col++){
        //        double varB = 0.0;
        //        for(int row = 0; row < numAncPops; row++){
        //          varB += finalBeliefs[windex][row][col];
        //        }
        //        cout << varB << " ";
        //      }
        //      cout << endl;
        //      cout << "Left ";
        //      if(windex > 0 && windex <=numWindows-2){
        //        for(int col = 0; col < numAncPops; col++){
        //          cout << leftMessages[windex][col] << " ";
        //        }
        //      }
        //      cout << endl;
        //      cout << "Right ";
        //      if(windex >= 0 && windex <=numWindows-3){
        //        for(int col = 0; col < numAncPops; col++){
        //          cout << rightMessages[windex][col] << " ";
        //        }
        //      }
        //      cout << endl;
        //      cout << "rfProb ";
        //      for(int ancIndex = 0; ancIndex < numAncPops; ancIndex++){
        //       cout << rfProbs[admIndex][(windex)*numAncPops + ancIndex] << " ";
        //      }
        //      cout << endl << endl;
        //    }
        // END DEBUG
        
        // Use row variable in each cluster except for last window
        for(int windex = 0; windex < numWindows-1; windex++){
          for(int ancIndex = 0; ancIndex < numAncPops; ancIndex++){
            marginals[2*admInd + hapIndex][windex*numAncPops + ancIndex] = 0.0;
            for(int col = 0; col < numAncPops; col++){
              marginals[2*admInd + hapIndex][windex*numAncPops + ancIndex] += finalBeliefs[windex][ancIndex][col];
            }
          }
        }
        // Do last window which is from last cluster, shared with secondto last window
        int windex = numWindows-1;
        for(int ancIndex = 0; ancIndex < numAncPops; ancIndex++){
          marginals[2*admInd + hapIndex][windex*numAncPops + ancIndex] = 0.0;
          for(int row = 0; row < numAncPops; row++){
            marginals[2*admInd + hapIndex][windex*numAncPops + ancIndex] += finalBeliefs[windex-1][row][ancIndex];
          }
        }
      }
    }
    
    
    
    
    
    
    
    
    
    
    
    
    for(int windex = 0; windex < numWindows - 1; windex++){
      delete[] stateBestPaths[windex];
    }
    
    for(int windex = 0; windex < numWindows; windex++){
      numPhasingsSeen[windex] += numHetSitesPerInd[windex][admInd] + 1;
    }
  }
  
  // Write fixed alleles
  ofstream outfile;
  string emIterationString;
  stringstream emiterationStringStream;
  emiterationStringStream << emIteration;
  emIterationString = emiterationStringStream.str();
  string outputName = string(processedInput->outputName) + string(".allelesRephased") + emIterationString + string(".txt");
  outfile.open((outputName).c_str());
  for(int windex = 0; windex < numWindows; windex++){
    int numSNPsInWindow = (*windowEndIndexes)[windex] - (*windowBeginIndexes)[windex] + 1;
    for(int snpIndex = 0; snpIndex < numSNPsInWindow; snpIndex++){
      for(int spot = 0; spot < numAdm; spot++){
        outfile << (int) admWindows[windex][snpIndex*numAdm + spot];
      }
      outfile << endl;
    }
  }
  outfile.close();
  
  return windowCalls;
}
