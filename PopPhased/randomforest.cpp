#include "randomforest.h"
#include <math.h>
#include <vector>
#include <set>
#include <omp.h>

#define NOGOOD -1

void Tree(int numHomogSNPs, bool* snpHomogParam, int windex, int firstSNPIndex, int windowNumSNPs, unsigned int mtry, int numAncPops, int bootstrapSample[], int numAdm, int numAnc, int* ancHapClasses, vector<bool> * ancWindows, vector<bool> * admWindows, int admHere[], double** treeCounts, int depth, int minNodeSize);

double ** RandomForest(ProcessedInput* processedInput, double** admClassProbs, int emIteration){
  int numWindows = processedInput->numWindows;
  int numAncPops = processedInput->numAncPops;  // # of ancestral populations
  int numAdm = processedInput->numAdm;  // # of admixed haplotypes
  int numAnc = processedInput->numAnc;  // # of ancestral haplotypes
  int numTrees = processedInput->numTrees;
  vector<int> * windowBeginIndexes = processedInput->windowBeginIndexes;
  vector<int> * windowEndIndexes = processedInput->windowEndIndexes;
  int ** ancHapClassesPerWindow = processedInput->ancHapClassesPerWindow;
  vector<bool> * ancWindows = processedInput->ancWindows;
  //int ** admWindows = processedInput->admWindows;
  vector<bool> * admWindowsPhasings = processedInput->admWindowsPhasings;
  int bootstrapSampleSize = processedInput->bootstrapSampleSize;
  int bootstrappingMethod = processedInput->bootstrappingMethod;
  int** numHetSitesPerInd = processedInput->numHetSitesPerInd;
  int* numAdmPhasings = processedInput->numAdmPhasings;
  int minNodeSize = processedInput->minNodeSize;
  
  int numAdmInds = numAdm / 2;
  
  srand((unsigned int)time(0));
  
//  double ** treeCounts = new double*[numAdm];
//  for(int i = 0; i < numAdm; i++){
//    treeCounts[i] = new double[numWindows * numAncPops];
//    for(int j = 0; j < (numWindows * numAncPops); j++){
//      treeCounts[i][j] = 1.0;
//    }
//  }
  
  double ** treeCounts = new double*[numWindows];
  for(int windex = 0; windex < numWindows; windex++){
    treeCounts[windex] = new double[2*numAdmPhasings[windex] * numAncPops];
    for(int i = 0; i < (2*numAdmPhasings[windex] * numAncPops); i++){
      treeCounts[windex][i] = 1.0;
    }
  }
    
  #pragma omp parallel for
  for(int windex = 0; windex < numWindows; windex++){
    //cout << "WINDEX: " << windex << endl;
    //cout << omp_get_num_threads() << endl;
    // Phasing
    int thisWindowNumAdm = numAdmPhasings[windex]*2;
    
    //    cout << endl;
    //cout << "WINDOW: " << windex << endl;
    int * ancHapClasses = ancHapClassesPerWindow[windex];    // EM
    
    // Get number of haps in each ref class for stratified bootstrap sampling and EM
    // EM : Affects: numHapsInEachPop and popHapLocations, moved from getdata
    int * numHapsInEachPop = new int[numAncPops];
    for(int i = 0; i < numAncPops; i++){
      numHapsInEachPop[i] = 0;
    }
    for(int i = 0; i < numAnc; i++){
      numHapsInEachPop[ancHapClasses[i]]++;
    }
    
    // Get columns of each hap in each pop
    int ** popHapLocations = new int*[numAncPops];
    for(int i = 0; i < numAncPops; i++){
      popHapLocations[i] = new int[numHapsInEachPop[i]];
    }
    int currentSpotToAdd[numAncPops];  // helper
    for(int i = 0; i < numAncPops; i++){
      currentSpotToAdd[i] = 0;
    }
    for(int i = 0; i < numAnc; i++){
      int ancOfHap = ancHapClasses[i];
      popHapLocations[ancOfHap][currentSpotToAdd[ancOfHap]] = i;
      currentSpotToAdd[ancOfHap]++;
    }
    
    int windowNumSNPs = (*windowEndIndexes)[windex] - (*windowBeginIndexes)[windex] + 1;
    unsigned int mtry = (unsigned int) sqrt(windowNumSNPs);
    //    mtry = mtry * 2;
    //    mtry = (mtry>windowNumSNPs)?windowNumSNPs:mtry;
    int firstSNPIndex = (*windowBeginIndexes)[windex];
    int bootstrapSample[numAnc];  // holds number of times the ith haplotype got sampled
    
    if(emIteration != 0){
      int admHapStartIndex = 0;
      
      for(int admIndTesting = 0; admIndTesting < numAdmInds; admIndTesting++){
        //cout << "Admixed Individual: " << admIndTesting << endl;
        int thisAdmIndNumPhasings = numHetSitesPerInd[windex][admIndTesting] + 1;
        
        for(int treeIndex = 0; treeIndex < numTrees; treeIndex++){
          // cout << "TREE: " << treeIndex << endl;
          // Get bootstrapped sample
          // Initialize
          for(int i = 0; i < numAnc; i++){
            bootstrapSample[i] = 0;
          }
          
          if(bootstrappingMethod == 2){
            // do stratified bootstrap
            for(int ancestry = 0; ancestry < numAncPops; ancestry++){
              int numHapsInAnc = numHapsInEachPop[ancestry];
              int numHapsInAncInAdmInd = 0;
              if (ancHapClasses[admIndTesting*2] == ancestry){  // using adm = anc
                numHapsInAncInAdmInd++;
              }
              if (ancHapClasses[admIndTesting*2+1] == ancestry){
                numHapsInAncInAdmInd++;
              }
              if(numHapsInAnc == 0 || numHapsInAnc == numHapsInAncInAdmInd){
                continue;
              }
              for(int sampleIndex = 0; sampleIndex < numHapsInAnc; sampleIndex++){
                int ancHapSampled = rand() % numHapsInAnc;
                int hapColumn = popHapLocations[ancestry][ancHapSampled];
                while(hapColumn == admIndTesting*2 || hapColumn == admIndTesting*2+1){
                  ancHapSampled = rand() % numHapsInAnc;
                  hapColumn = popHapLocations[ancestry][ancHapSampled];
                }
                
                bootstrapSample[hapColumn] += 1;
              }
            }
          }
          else if(bootstrappingMethod == 1){
            for(int sampleIndex = 0; sampleIndex < bootstrapSampleSize; sampleIndex++){  // bootstrap the same number as there are ancestral haplotypes
              int ancSampled = rand() % numAncPops;
              int numHapsInAnc = numHapsInEachPop[ancSampled];
              int numHapsInAncInAdmInd = 0;
              if (ancHapClasses[admIndTesting*2] == ancSampled){  // using adm = anc
                numHapsInAncInAdmInd++;
              }
              if (ancHapClasses[admIndTesting*2+1] == ancSampled){
                numHapsInAncInAdmInd++;
              }
              while(numHapsInAnc == 0 || numHapsInAnc == numHapsInAncInAdmInd){
                // the second condition is there because will have infinite loop looking for a second hap to have as anc if the only one is in test set
                ancSampled = rand() % numAncPops;
                numHapsInAnc = numHapsInEachPop[ancSampled];
                numHapsInAncInAdmInd = 0;
                if (ancHapClasses[admIndTesting*2] == ancSampled){  // using adm = anc
                  numHapsInAncInAdmInd++;
                }
                if (ancHapClasses[admIndTesting*2+1] == ancSampled){
                  numHapsInAncInAdmInd++;
                }
              }
              int ancHapSampled = rand() % numHapsInAnc;
              int hapColumn = popHapLocations[ancSampled][ancHapSampled];
              while(hapColumn == admIndTesting*2 || hapColumn == admIndTesting*2+1){
                ancHapSampled = rand() % numHapsInAnc;
                hapColumn = popHapLocations[ancSampled][ancHapSampled];
              }
              bootstrapSample[hapColumn] += 1;
            }
          }
          else{
            for(int sampleIndex = 0; sampleIndex < bootstrapSampleSize; sampleIndex++){  // bootstrap the same number as there are ancestral haplotypes
              int hapColumn = rand() % numAnc;
              while(hapColumn == admIndTesting*2 || hapColumn == admIndTesting*2+1){
                hapColumn = rand() % numAnc;
              }
              bootstrapSample[hapColumn] += 1;
            }
          }

          
          // this is alligned with admWindowsPhasings, not admWindows
          int admHere[thisWindowNumAdm];  // tracks which admixed are in current branch
          for (int i = 0; i < thisWindowNumAdm; i++){
              admHere[i] = 0;
          }
          
          for(int phaseHapIndex = admHapStartIndex; phaseHapIndex < admHapStartIndex+thisAdmIndNumPhasings*2; phaseHapIndex++){
            admHere[phaseHapIndex] = 1;
          }
          
          bool snpHomog[windowNumSNPs];  // tracks which snps are homogeneous in current branch
          for (int i = 0; i < windowNumSNPs; i++){
            snpHomog[i] = false;
          }
          
          // Debug
//          if(treeIndex == 0){
//            cout << "====" << endl;
//            cout << admIndTesting << endl;
//            cout << admHapStartIndex <<endl;
//            cout << thisAdmIndNumPhasings << endl;
//            for (int i = 0; i < thisWindowNumAdm; i++){
//              cout << admHere[i] << " ";
//            }
//            cout << endl;
//            for(int i = 0; i < numAnc; i++){
//              cout << bootstrapSample[i] << " ";
//            }
//            cout << endl;
//          }
          
          Tree(0, snpHomog, windex, firstSNPIndex, windowNumSNPs, mtry, numAncPops, bootstrapSample, thisWindowNumAdm, numAnc, ancHapClasses, ancWindows, admWindowsPhasings, admHere, treeCounts, 0, minNodeSize);
        }
        
        admHapStartIndex += thisAdmIndNumPhasings * 2;
      }
    } else{
      for(int treeIndex = 0; treeIndex < numTrees; treeIndex++){
        //if(windex == 517){
        //  cout << "\nTREE: " << treeIndex << endl;
        //}

        // Get bootstrapped sample
        // Initialize
        for(int i = 0; i < numAnc; i++){
          bootstrapSample[i] = 0;
        }
        
        if(bootstrappingMethod == 2){
          // do stratified bootstrap
          for(int ancestry = 0; ancestry < numAncPops; ancestry++){
            int numHapsInAnc = numHapsInEachPop[ancestry];
            for(int sampleIndex = 0; sampleIndex < numHapsInAnc; sampleIndex++){
              int ancHapSampled = rand() % numHapsInAnc;
              int hapColumn = popHapLocations[ancestry][ancHapSampled];
              bootstrapSample[hapColumn] += 1;
            }
          }
        }
        else if(bootstrappingMethod == 1){
          for(int sampleIndex = 0; sampleIndex < bootstrapSampleSize; sampleIndex++){  // bootstrap the same number as there are ancestral haplotypes
            int ancSampled = rand() % numAncPops;
            int numHapsInAnc = numHapsInEachPop[ancSampled];
            int ancHapSampled = rand() % numHapsInAnc;
            int hapColumn = popHapLocations[ancSampled][ancHapSampled];
            bootstrapSample[hapColumn] += 1;
          }
        }
        else{
          for(int sampleIndex = 0; sampleIndex < bootstrapSampleSize; sampleIndex++){  // bootstrap the same number as there are ancestral haplotypes
            bootstrapSample[rand() % numAnc] += 1;
          }
        }
        
        int admHere[thisWindowNumAdm];  // tracks which admixed are in current branch
        for (int i = 0; i < thisWindowNumAdm; i++){
          admHere[i] = 1;
        }
        
        bool snpHomog[windowNumSNPs];  // tracks which snps are homogenous in current branch
        for (int i = 0; i < windowNumSNPs; i++){
          snpHomog[i] = false;
        }
        
        Tree(0, snpHomog, windex, firstSNPIndex, windowNumSNPs, mtry, numAncPops, bootstrapSample, thisWindowNumAdm, numAnc, ancHapClasses, ancWindows, admWindowsPhasings, admHere, treeCounts, 0, minNodeSize);
      }
    }
    delete[] numHapsInEachPop;
    for(int i = 0; i < numAncPops; i++){
      delete[] popHapLocations[i];
    }
    delete[] popHapLocations;
  }
  
//  cout << "Tree Counts\n";
//  for(int i = 0; i < numAdm; i++){
//    for(int j = 0; j < (numWindows * numAncPops); j++){
//      cout << treeCounts[i][j] << " " << endl;
//    }
//    cout << endl;
//  }
  
  //Convert counts to probabilities
//  double** admClassProbs = new double*[numAdm];    EM
  double numTreesDouble = (double) numTrees;
  double numAncPopsDouble = (double) numAncPops;
  for(int windex = 0; windex < numWindows; windex++){
    for(int admHap = 0; admHap < 2*numAdmPhasings[windex]; admHap++){
      for(int ancestry = 0; ancestry < numAncPops; ancestry++){
        admClassProbs[windex][admHap*numAncPops + ancestry] = treeCounts[windex][admHap*numAncPops + ancestry] / (numTreesDouble+numAncPopsDouble);
      }
    }
  }
  
  //Clean up
  for(int windex = 0; windex < numWindows; windex++){
    delete[] treeCounts[windex];
  }
  delete[] treeCounts;
  
  return admClassProbs;
}

void Tree(int numHomogSNPs, bool* snpHomogParam, int windex, int firstSNPIndex, int windowNumSNPs, unsigned int mtry, int numAncPops, int bootstrapSample[], int numAdm, int numAnc, int * ancHapClasses, vector<bool>* ancWindows, vector<bool>* admWindows, int admHere[], double** treeCounts, int depth, int minNodeSize){
  bool * snpHomog = new bool[windowNumSNPs];
  for(int i = 0; i < windowNumSNPs; i++){
    snpHomog[i] = snpHomogParam[i];
  }
  
  // check if done due to only one ref class left
  set<int> classesLeft;
  for(int i = 0; i < numAnc; i++){
    if(bootstrapSample[i] > 0){
      classesLeft.insert(ancHapClasses[i]);
    }
  }
  //  cout << "Classes Left: " << classesLeft.size() << endl;
  if(classesLeft.size() == 1){
    // only one type of class maps here, call all adm mapping here as that class
    set<int>::iterator clp = classesLeft.begin();
    int calledClass = *clp;
    
    for(int i = 0; i < numAdm; i++){
      //treeCounts[i][windex*numAncPops + calledClass] += (double)admHere[i];
      treeCounts[windex][i*numAncPops + calledClass] += (double)admHere[i];
    }
  } 
  //  else if(classesLeft.size() == 0){
  //    cout << "Error: node with no classes" << endl;
  //    exit(1);
  //  }
  else{
    int indexOfMinGini = NOGOOD;
    int * classCounts0 = new int[numAncPops];
    int * classCounts1 = new int[numAncPops];
    while(indexOfMinGini == NOGOOD){
      // sample without replacement mtry snps from the number of snps in the window
      set<int> sampledSNPs;
      //    cout << "Sampled SNPs: ";
      while(sampledSNPs.size() < mtry){
        int sample = rand() % windowNumSNPs;  // + firstSNPIndex
        sampledSNPs.insert(sample);
        //      cout << sample << " ";
      }
      //    cout << endl;
      
      // calculate Gini for each SNP
      double minGini = 1.0;
      int sumClassCounts0 = 0;
      int sumClassCounts1 = 0;
      for(int i = 0; i < numAncPops; i++){
        classCounts0[i] = 0;
        classCounts1[i] = 0;
      }
      set<int>::iterator it;
      for (it=sampledSNPs.begin(); it != sampledSNPs.end(); it++){
        int snp = *it;
        
        if(snpHomog[snp]){
          // only one allele so does not divide. This skips homog snps that were found to be homog previously.
          continue;
        }
        
        // Get number of haps from each class going along 0 or 1 branches
        for(int ancIndex = 0; ancIndex < numAnc; ancIndex++){
          int pop = ancHapClasses[ancIndex];
          int addTo1 = bootstrapSample[ancIndex]*ancWindows[windex][snp*numAnc + ancIndex];
          int addTo0 = bootstrapSample[ancIndex] - addTo1;
          classCounts0[pop] += addTo0;
          classCounts1[pop] += addTo1;
          sumClassCounts0 += addTo0;
          sumClassCounts1 += addTo1;
        }
        
        
        //DEBUG
        //      cout << "SNP " << snp<< " Class Counts" << endl;
        //      for(int i = 0; i < numAncPops; i++){
        //        cout << "Going 0, Class " << i << " Count " << classCounts0[i] << endl;
        //        cout << "Going 1, Class " << i << " Count " << classCounts1[i] << endl;
        //      }
        //ENDDEBUG
        
        if (sumClassCounts0 == 0 || sumClassCounts1 == 0){
          // new homog snp. this SNP sends all anc along one branch only, will also give nan when calculating gini below        
          numHomogSNPs++;
          snpHomog[snp] = true;
          if(numHomogSNPs == windowNumSNPs){
            // To prevent a branch with no ancestral going both 0 and 1 stop here and call adm on this branch
            int classVotes[numAncPops];
            for(int i = 0; i < numAncPops; i++){
              classVotes[i] = 0;
            }
            
            int totalVotes = 0;
            
            for(int i = 0; i < numAnc; i++){
              int thisClass = ancHapClasses[i];
              classVotes[thisClass] += bootstrapSample[i];
              totalVotes += bootstrapSample[i];
            }
            
            for(int i = 0; i < numAdm; i++){
              if(admHere[i] == 1){
                for(int j = 0; j < numAncPops; j++){
                  treeCounts[windex][i*numAncPops + j] += (double)classVotes[j] / (double)totalVotes;
                  //treeCounts[i][windex*numAncPops + j] += (double)classVotes[j] / (double)totalVotes;
                }
              }
            }
            
            delete[] snpHomog;
            delete[] classCounts0;
            delete[] classCounts1;
            return;
          }
        }
        else{
          double gini0 = 1.0;
          double gini1 = 1.0;
          for(int i = 0; i < numAncPops; i++){
            gini0 -= pow((double)classCounts0[i]/(double)sumClassCounts0, 2);
            gini1 -= pow((double)classCounts1[i]/(double)sumClassCounts1, 2);
          }
          
          int sumAllClassCounts = sumClassCounts0 + sumClassCounts1;
          double gini = ((double)sumClassCounts0 / (double)sumAllClassCounts)*gini0 + ((double)sumClassCounts1 / (double)sumAllClassCounts)*gini1;
          
          //        cout << "Gini " << gini << endl;
          
          if(gini < minGini){
            minGini = gini;
            indexOfMinGini = snp;
          }
          
          //        cout << "MinIndex " << indexOfMinGini << endl;
        }
        
        for(int i = 0; i < numAncPops; i++){
          classCounts0[i] = 0;
          classCounts1[i] = 0;
        }
        sumClassCounts0 = 0;
        sumClassCounts1 = 0;
      }
      
      // Sample obtained
    }
    
    delete[] classCounts0;
    classCounts0 = 0;
    delete[] classCounts1;
    classCounts1 = 0;
    
    bool go0 = true;
    bool go1 = true;
    
    // Branch anc to 0 and 1
    int * bootstrapSample0 = new int[numAnc];
    int * bootstrapSample1 = new int[numAnc];
    for(int ancIndex = 0; ancIndex < numAnc; ancIndex++){
      bootstrapSample1[ancIndex] = bootstrapSample[ancIndex]*ancWindows[windex][indexOfMinGini*numAnc+ancIndex];
      bootstrapSample0[ancIndex] = bootstrapSample[ancIndex] - bootstrapSample1[ancIndex];
    }
    
    // Branch adm to 0 and 1
    int* admHere0 = new int[numAdm];
    int* admHere1 = new int[numAdm];
    int numAdmGoing0 = 0;
    int numAdmGoing1 = 0;
    for (int i = 0; i < numAdm; i++){
      admHere1[i] = admWindows[windex][indexOfMinGini*numAdm + i]*admHere[i];
      numAdmGoing1 += admHere1[i];
      admHere0[i] = admHere[i] - admHere1[i];    // speedup by using local variable to store admHere0[] and admHere1[]
      numAdmGoing0 += admHere0[i];
    }
    
    //    cout << "Num Adm going 0: " << numAdmGoing0 << endl;
    //    cout << "Num Adm going 1: " << numAdmGoing1 << endl;
    
    // Don't go down a branch if no admixed going that way: NEED to delete this if want to use the trees later
    if(numAdmGoing0 == 0){
      go0 = false;
    } else if(numAdmGoing1 == 0){
      go1 = false;
    }
    
    // Stop here if a child node size is below threshold
    int numAncGoing0 = 0;
    int numAncGoing1 = 0;
    for (int i = 0; i < numAnc; i++){
      numAncGoing0 += bootstrapSample0[i];
      numAncGoing1 += bootstrapSample1[i];
    }
    
    if(numAncGoing0 < minNodeSize){
      go0 = false;
      
      int classVotes[numAncPops];
      for(int i = 0; i < numAncPops; i++){
        classVotes[i] = 0;
      }
      
      int totalVotes = 0;
      
      for(int i = 0; i < numAnc; i++){
        int thisClass = ancHapClasses[i];
        classVotes[thisClass] += bootstrapSample[i];
        totalVotes += bootstrapSample[i];
      }
      
      for(int i = 0; i < numAdm; i++){
        if(admHere0[i] == 1){
          for(int j = 0; j < numAncPops; j++){
            treeCounts[windex][i*numAncPops + j] += (double)classVotes[j] / (double)totalVotes;
          }
        }
      }
    }
    if (numAncGoing1 < minNodeSize){
      go1 = false;
      
      int classVotes[numAncPops];
      for(int i = 0; i < numAncPops; i++){
        classVotes[i] = 0;
      }
      
      int totalVotes = 0;
      
      for(int i = 0; i < numAnc; i++){
        int thisClass = ancHapClasses[i];
        classVotes[thisClass] += bootstrapSample[i];
        totalVotes += bootstrapSample[i];
      }
      
      for(int i = 0; i < numAdm; i++){
        if(admHere1[i] == 1){
          for(int j = 0; j < numAncPops; j++){
            treeCounts[windex][i*numAncPops + j] += (double)classVotes[j] / (double)totalVotes;
          }
        }
      }
    }
    
    // Branch 0
    if(go0){
      //      cout << "Going 0\n";
      Tree(numHomogSNPs, snpHomog, windex, firstSNPIndex, windowNumSNPs, mtry, numAncPops, bootstrapSample0, numAdm, numAnc, ancHapClasses, ancWindows, admWindows, admHere0, treeCounts, depth+1, minNodeSize);
      //      cout << "\nDepth: " << depth << " Back from 0\n";
    }
    
    // Branch 1
    if(go1){
      //      cout << "Going 1\n";
      Tree(numHomogSNPs, snpHomog, windex, firstSNPIndex, windowNumSNPs, mtry, numAncPops, bootstrapSample1, numAdm, numAnc, ancHapClasses, ancWindows, admWindows, admHere1, treeCounts, depth+1, minNodeSize);
      //      cout << "\nDepth: " << depth << " Back from 1\n";
    }
    
    delete[] bootstrapSample0;
    bootstrapSample0 = 0;
    delete[] bootstrapSample1;
    bootstrapSample1 = 0;
    delete[] admHere0;
    admHere0 = 0;
    delete[] admHere1;
    admHere1 = 0;
  }
  delete[] snpHomog;
}
