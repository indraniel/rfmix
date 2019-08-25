#include "randomforest.h"
#include <math.h>
#include <vector>
#include <set>
#include <omp.h>

#define NOGOOD -1

void Tree(int numHomogSNPs, bool* snpHomogParam, int windex, int firstSNPIndex, int windowNumSNPs, unsigned int mtry, int numAncPops, int bootstrapSample[], int numAdm, int numAnc, int* ancHapClasses, char ** ancWindows, char ** admWindows, int admHere[], double** treeCounts, int depth, int minNodeSize);

double ** RandomForest(ProcessedInput* processedInput, double** admClassProbs, int emIteration){
  int numWindows = processedInput->numWindows;
  int numAncPops = processedInput->numAncPops;  // # of ancestral populations
  int numAdm = processedInput->numAdm;  // # of admixed haplotypes
  int numAnc = processedInput->numAnc;  // # of ancestral haplotypes
  int numTrees = processedInput->numTrees;
  vector<int> * windowBeginIndexes = processedInput->windowBeginIndexes;
  vector<int> * windowEndIndexes = processedInput->windowEndIndexes;
  int ** ancHapClassesPerWindow = processedInput->ancHapClassesPerWindow;
  char ** ancWindows = processedInput->ancWindows;
  char ** admWindows = processedInput->admWindows;
  int bootstrapSampleSize = processedInput->bootstrapSampleSize;
  int bootstrappingMethod = processedInput->bootstrappingMethod;
  unsigned int mtryFactor = processedInput->mtryFactor;
  int minNodeSize = processedInput->minNodeSize;
  
  int numDivisions = 10;
  
  srand((int) time(0));
  
  double ** treeCounts = new double*[numAdm];
  for(int i = 0; i < numAdm; i++){
    treeCounts[i] = new double[numWindows * numAncPops];
    for(int j = 0; j < (numWindows * numAncPops); j++){
      treeCounts[i][j] = 1.0;
    }
  }
  
  #pragma omp parallel for
  for(int windex = 0; windex < numWindows; windex++){
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
    int mtry = (int) sqrt(windowNumSNPs);
    mtry = mtry * mtryFactor;
    mtry = (mtry>windowNumSNPs)?windowNumSNPs:mtry;
    //cout << "Window: " << windex << " NumSNPS: " << windowNumSNPs << "mtry: " << mtry << endl;
    int firstSNPIndex = (*windowBeginIndexes)[windex];
    int bootstrapSample[numAnc];  // holds number of times the ith haplotype got sampled
    
    if(emIteration != 0){
      // Use this instead of numPerBucket because some buckets will have an extra hap added
      int popNumPerDivision[numAncPops][numDivisions];
      for(int ancestry = 0; ancestry < numAncPops; ancestry++){
        int numHapsInAnc = numHapsInEachPop[ancestry];
        for(int divisionIndex = 0; divisionIndex < numDivisions; divisionIndex++){
          int numPerBucket = numHapsInAnc / numDivisions;
          popNumPerDivision[ancestry][divisionIndex] = numPerBucket;
        }
      }
      
      // Just drawing numPerBucket = numHapsInAnc / numDivisions will leave out 6 if 16 haps in pop and numDivisions == 10
      // Thus, pick which divisions get to pick an extra hap, for each ancestry
      // All this info is put into popNumPerDivision initialized above
      bool extraHapDivisions[numDivisions][numAncPops];
      for(int divisionIndex = 0; divisionIndex < numDivisions; divisionIndex++){
        for(int ancestry = 0; ancestry < numAncPops; ancestry++){
          extraHapDivisions[divisionIndex][ancestry] = false;
        }
      }
      for(int ancestry = 0; ancestry < numAncPops; ancestry++){
        int numHapsInAnc = numHapsInEachPop[ancestry];
        int numPerBucket = numHapsInAnc / numDivisions;
        int numExtra = numHapsInAnc - numPerBucket * numDivisions;
        for(int extraIndex = 0; extraIndex < numExtra; extraIndex++){
          int drawnBucket = rand() % numDivisions;
          while(extraHapDivisions[drawnBucket][ancestry]){
            drawnBucket = rand() % numDivisions;
          }
          extraHapDivisions[drawnBucket][ancestry] = true;
          popNumPerDivision[ancestry][drawnBucket]++;
        }
      }
      
      int maxHapsInPop = 0;
      for(int i = 0; i < numAncPops; i++){
        int numHapsInAnc = numHapsInEachPop[i];
        if(numHapsInAnc > maxHapsInPop){
          maxHapsInPop = numHapsInAnc;
        }
      }
      
      bool alreadyTested[numAncPops][maxHapsInPop];
      for(int i = 0; i < numAncPops; i++){
        for(int j = 0; j < maxHapsInPop; j++){
          alreadyTested[i][j] = false;
        }
      }
      
//      int ** alreadyTested = new int*[numAncPops];  // 1 if haplotype already used in test set
//      for (int i = 0; i < numAncPops; i++){
//        int numHapsInAnc = numHapsInEachPop[i];
//        alreadyTested[i] = new int[numHapsInAnc];
//        for(int j = 0; j < numHapsInAnc; j++){
//          alreadyTested[i][j] = 0;
//        }
//      }
      
      int thisTestHapFormat[numAnc];  // stores which haps in alleles file are used for adm test
      
      for(int divisionIndex = 0; divisionIndex < numDivisions; divisionIndex++){
        // Split into training and test sets
        // Because not pre-EM, the Anc and Adm alleles are the same. Either original anc are thrown out or they are operated on too
        // Divide each population into numDividions equal sized pieces
        int ** thisTest = new int*[numAncPops]; // hapIndices of those used in adm "test" set
        for (int ancestry = 0; ancestry < numAncPops; ancestry++){
          thisTest[ancestry] = new int[popNumPerDivision[ancestry][divisionIndex]];
        }
        
        for(int i = 0; i < numAdm; i++){
          thisTestHapFormat[i] = 0;
        }
        
        // Get adm test set: will need to convert index in pop to index in alleles file: popHapLocations
        for(int ancestry = 0; ancestry < numAncPops; ancestry++){
          int numHapsInAnc = numHapsInEachPop[ancestry];
          int numPerBucket = popNumPerDivision[ancestry][divisionIndex];
          for(int drawIndex = 0; drawIndex < numPerBucket; drawIndex++){
            int drawn = rand() % numHapsInAnc;
            while(alreadyTested[ancestry][drawn]){
              drawn = rand() % numHapsInAnc;
            }
            alreadyTested[ancestry][drawn] = true;
            thisTest[ancestry][drawIndex] = drawn;
            thisTestHapFormat[popHapLocations[ancestry][drawn]] = 1;
          }
        }
        
        for(int treeIndex = 0; treeIndex < numTrees; treeIndex++){
          //      cout << "\nTREE: " << treeIndex << endl;
          // Get bootstrapped sample
          // Initialize
          for(int i = 0; i < numAnc; i++){
            bootstrapSample[i] = 0;
          }
          
          if(bootstrappingMethod == 2){
            // do stratified bootstrap
            for(int ancestry = 0; ancestry < numAncPops; ancestry++){
              int numHapsInAnc = numHapsInEachPop[ancestry];
              if(numHapsInAnc == 0 || (numHapsInAnc == 1 && thisTestHapFormat[popHapLocations[ancestry][0]])){
                continue;
              }
              for(int sampleIndex = 0; sampleIndex < numHapsInAnc; sampleIndex++){
                int ancHapSampled = rand() % numHapsInAnc;
                int hapColumn = popHapLocations[ancestry][ancHapSampled];
                while(thisTestHapFormat[hapColumn]){
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
              while(numHapsInAnc == 0 || (numHapsInAnc == 1 && thisTestHapFormat[popHapLocations[ancSampled][0]])){
                // the second condition is there because will have infinite loop looking for a second hap to have as anc if the only one is in test set
                ancSampled = rand() % numAncPops;
                numHapsInAnc = numHapsInEachPop[ancSampled];
              }
              int ancHapSampled = rand() % numHapsInAnc;
              int hapColumn = popHapLocations[ancSampled][ancHapSampled];
              while(thisTestHapFormat[hapColumn]){
                ancHapSampled = rand() % numHapsInAnc;
                hapColumn = popHapLocations[ancSampled][ancHapSampled];
              }
              bootstrapSample[hapColumn] += 1;
            }
          }
          else{
            for(int sampleIndex = 0; sampleIndex < bootstrapSampleSize; sampleIndex++){  // bootstrap the same number as there are ancestral haplotypes
              int hapColumn = rand() % numAnc;
              while(thisTestHapFormat[hapColumn]){
                hapColumn = rand() % numAnc;
              }
              bootstrapSample[hapColumn] += 1;
            }
          }
          
          int admHere[numAdm];  // tracks which admixed are in current branch
          for (int i = 0; i < numAdm; i++){
            admHere[i] = 0;
          }
          // Only include those in this adm test set
          for(int ancestry = 0; ancestry < numAncPops; ancestry++){
            int numPerBucket = popNumPerDivision[ancestry][divisionIndex];
            for(int drawIndex = 0; drawIndex < numPerBucket; drawIndex++){
              int drawn = thisTest[ancestry][drawIndex];
              admHere[popHapLocations[ancestry][drawn]] = 1;
            }
          }
          
          bool snpHomog[windowNumSNPs];  // tracks which snps are homogenous in current branch
          for (int i = 0; i < windowNumSNPs; i++){
            snpHomog[i] = false;
          }
          
          Tree(0, snpHomog, windex, firstSNPIndex, windowNumSNPs, mtry, numAncPops, bootstrapSample, numAdm, numAnc, ancHapClasses, ancWindows, admWindows, admHere, treeCounts, 0, minNodeSize);
        }
        
        for (int ancestry = 0; ancestry < numAncPops; ancestry++){
          delete[] thisTest[ancestry];
        }
        delete[] thisTest;
      }
    } else{      
      for(int treeIndex = 0; treeIndex < numTrees; treeIndex++){
        //      cout << "\nTREE: " << treeIndex << endl;
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
        
        int admHere[numAdm];  // tracks which admixed are in current branch
        for (int i = 0; i < numAdm; i++){
          admHere[i] = 1;
        }
        
        bool snpHomog[windowNumSNPs];  // tracks which snps are homogenous in current branch
        for (int i = 0; i < windowNumSNPs; i++){
          snpHomog[i] = false;
        }
        
        Tree(0, snpHomog, windex, firstSNPIndex, windowNumSNPs, mtry, numAncPops, bootstrapSample, numAdm, numAnc, ancHapClasses, ancWindows, admWindows, admHere, treeCounts, 0, minNodeSize);
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
  for(int i = 0; i < numAdm; i++){
    //    cout << "Individual: " << i << endl;
    //admClassProbs[i] = new double[numWindows * numAncPops];    EM
    for(int j = 0; j < numWindows; j++){
      for(int k = 0; k < numAncPops; k++){
        admClassProbs[i][j*numAncPops+k] = treeCounts[i][j*numAncPops+k] / (numTreesDouble+numAncPopsDouble);
        //        cout << admClassProbs[i][j*numAncPops+k] << " ";
      }
      //      cout << endl;
    }
    //    cout << endl;
  }
  
  //Clean up
  for(int i = 0; i < numAdm; i++){
    delete[] treeCounts[i];
  }
  delete[] treeCounts;
  
  return admClassProbs;
}

void Tree(int numHomogSNPs, bool* snpHomogParam, int windex, int firstSNPIndex, int windowNumSNPs, unsigned int mtry, int numAncPops, int bootstrapSample[], int numAdm, int numAnc, int * ancHapClasses, char ** ancWindows, char ** admWindows, int admHere[], double** treeCounts, int depth, int minNodeSize){
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
      treeCounts[i][windex*numAncPops + calledClass] += (double)admHere[i];
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
                  treeCounts[i][windex*numAncPops + j] += (double)classVotes[j] / (double)totalVotes;
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
            treeCounts[i][windex*numAncPops + j] += (double)classVotes[j] / (double)totalVotes;
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
            treeCounts[i][windex*numAncPops + j] += (double)classVotes[j] / (double)totalVotes;
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
