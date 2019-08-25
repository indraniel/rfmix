#include "crfviterbi.h"
#include <float.h>

// Now with forward-backward goodness!
int** CrfViterbi(ProcessedInput* processedInput, double ** rfProbs, int** windowCalls, double** marginals){
  int numWindows = processedInput->numWindows;
  int numAncPops = processedInput->numAncPops;  // # of ancestral populations
  int numAdm = processedInput->numAdm;  // # of admixed haplotypes
  int genSinceAdmixture = processedInput->genSinceAdmixture;
  vector<double> * windowBeginLocs = processedInput->windowBeginLocs;
  vector<double> * windowEndLocs = processedInput->windowEndLocs;
  
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
//    cout << transitionProbs[windex] << endl;
  }
  
  // Perform Viterbi
  int stateBestPaths[numWindows - 1][numAncPops]; // doesn't need to be zero'd
  double cumulativeProbs[numAncPops];  // stores sum of logProbs of best path that reaches this state at the current time
//  int** windowCalls = new int*[numAdm];  EM
//  for(int i = 0; i < numAdm; i++){
//    windowCalls[i] = new int[numWindows];
//  }

  double cumulativeProbsHelper[numAncPops];
  
  for(int admIndex = 0; admIndex < numAdm; admIndex++){
    // Handle initial window
    for(int ancIndex = 0; ancIndex < numAncPops; ancIndex++){
      cumulativeProbs[ancIndex] = log(rfProbs[admIndex][ancIndex]);
    }
    
    // Do rest
    for(int windex = 0; windex < numWindows-1; windex++){  // starts at 0 because recording state we came from
      for(int thisAnc = 0; thisAnc < numAncPops; thisAnc++){
        double maxScore = -DBL_MAX;
        int lastMaxAnc = 0;
        double thisAncProb = log(rfProbs[admIndex][(windex+1)*numAncPops + thisAnc]);
        for(int lastAnc = 0; lastAnc < numAncPops; lastAnc++){
          double score;
          if(lastAnc == thisAnc){
            score = cumulativeProbs[lastAnc] + log((1 - transitionProbs[windex])*(1.0/((double)numAncPops)) + transitionProbs[windex]*(1.0/((double)numAncPops))*(1.0/((double)numAncPops))) + thisAncProb;
          }else{
            score = cumulativeProbs[lastAnc] + log(transitionProbs[windex]*(1.0/((double)numAncPops))*(1.0/((double)numAncPops))) + thisAncProb;
          }
          if (score > maxScore){
            maxScore = score;
            lastMaxAnc = lastAnc;
          }
        }
        stateBestPaths[windex][thisAnc] = lastMaxAnc;
        cumulativeProbsHelper[thisAnc] = maxScore;
      }
      for(int ancIndex = 0; ancIndex < numAncPops; ancIndex++){
        cumulativeProbs[ancIndex] = cumulativeProbsHelper[ancIndex];
      }
    }
    
    // Get max cumulative probs and backtrace
    double maxScore = -DBL_MAX;
    int maxIndex = 0;
    for(int ancIndex = 0; ancIndex < numAncPops; ancIndex++){
      double score = cumulativeProbs[ancIndex];
//      cout << score << endl;
      if(score > maxScore){
        maxScore = score;
        maxIndex = ancIndex;
      }
    }
    
    int pointedTo = maxIndex;
    windowCalls[admIndex][numWindows-1] = pointedTo;
    for(int windex = numWindows-2; windex >=0; windex--){
      pointedTo = stateBestPaths[windex][pointedTo];    // sbp has length numwindows-1
      windowCalls[admIndex][windex] = pointedTo;
    }
    
    // ==== Forward Backward ====
    if(processedInput->doForwardBackward){
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
            initialProbs[windex][i][j] *= rfProbs[admIndex][(windex)*numAncPops + i];
            if(windex == numWindows-2){    // last cluster so multiply in last rfProb variable
              initialProbs[windex][i][j] *= rfProbs[admIndex][(windex+1)*numAncPops + j];
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
          marginals[admIndex][windex*numAncPops + ancIndex] = 0.0;
          for(int col = 0; col < numAncPops; col++){
            marginals[admIndex][windex*numAncPops + ancIndex] += finalBeliefs[windex][ancIndex][col];
          }
        }
      }
      // Do last window which is from last cluster, shared with secondto last window
      int windex = numWindows-1;
      for(int ancIndex = 0; ancIndex < numAncPops; ancIndex++){
        marginals[admIndex][windex*numAncPops + ancIndex] = 0.0;
        for(int row = 0; row < numAncPops; row++){
          marginals[admIndex][windex*numAncPops + ancIndex] += finalBeliefs[windex-1][row][ancIndex];
        }
      }
    }
  }
  
//  for(int admIndex = 0; admIndex < numAdm; admIndex++){
//    cout << "Individual "<<admIndex<<endl;
//    for(int windex=0; windex < numWindows; windex++){
//      cout << windowCalls[admIndex][windex]<<endl;
//    }
//  }
  
  return windowCalls;
}
