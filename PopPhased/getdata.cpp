/*
 *  Parameters:
 *    -w: window size
 *    -G: num generations
 *    -o: output file name
 *    -m: marker locations
 *    -a: file with phased adm and ancestral alleles, one column per haplotype
 *    -p: population of each, 0 being admixed
 *    -t: number of trees to grow in random forest
 *    -r: allow OpenMP (shared memory) parallelization. 1-true, 0-false
 *    -h: number of threads to use if doing parallelization
 *    -b: number of bootstrap samples used per tree
 *    -s: kind of bootstrap sampling. 0-sample with equal prob from each haplotype, 1-sample with equal prob from each class, 2-stratified by class
 *    -e: number of EM iterations. Default is 0. Non-EM round does inference on ref panels too, EM keeps ref panels as adm. Idea is if there's a lot of adm then the signal will override the prior (from anc)
 *    -u: 0-don't use anc as anc for iterations after the first, 1-otherwise
 *    -c: 0-don't make calls on anc as if they were adm, 1-do so after first iteration. UNUSED
 *    -x: file with indices of SNPs to prune. 0 indexed.
 *    -n: minimum number of ancestral haplotypes to have in a tree node
 *    -fb: 1-Calculate and output forward backward, 0-Don't
 *    -co: compressed output: output window-based calls with additional file listing # SNPs per window
 */

#include "getdata.h"
#include <omp.h>

ProcessedInput* processInput(int argc, char ** argv, ofstream * logfile){
  // ===== Step 1: Make sure an input file has been specified =====
  if(argc<2){
    cout<<"Error: No parameters\n";
    exit(1);
  }
    
  // ===== Step 2: Get the parameters =====
  map<string, string> params;
  vector<string> ancfiles;
  
  //Set default values
  params["-w"]="0.2";
  params["-G"]="8";
  params["-o"]="RFMix_results";
  params["-t"]="100";
  params["-r"]="1";
  params["-s"]="1";
  params["-e"]="0";
  params["-u"]="0";
  params["-c"]="0";
  params["-x"]="";
  params["-n"]="1";
  params["-fb"]="0";
  params["-co"]="0";
  
  // -h parameter requires some extra work. Initializes to default setting
  ostringstream threadSC;
  threadSC << omp_get_num_procs();
  params["-h"]=threadSC.str();
  
  //Parse commandline: assumes adm and ancestral all in one file
  for(int i=1; i<argc-1; i+=2){
    string s(argv[i]);
    params[argv[i]]=argv[i+1];
  }
  
  int numThreads;
  istringstream threadStringConverter(params["-h"]);
  if(!(threadStringConverter >> numThreads)){
    cout << "Invalid number of threads" << endl;
    exit(1);
  }
  omp_set_num_threads(numThreads); // Note: if multiple parallel regions, I think this only applies to the first
  
  // Note: check if person wants parallelization below. If not, numThreads set to 1
  
  cout << "Number of processors available: " << omp_get_num_procs() << endl;
  
  int numTrees;
  istringstream treeStringConverter(params["-t"]);
  if(!(treeStringConverter >> numTrees)){
    cout << "Invalid number of trees" << endl;
    exit(1);
  }
  
  int bootstrappingMethod;
  istringstream stratStringConverter(params["-s"]);
  if(!(stratStringConverter >> bootstrappingMethod)){
    cout << "Invalid bootstrap method flag, should be nonnegative integer" << endl;
    exit(1);
  }
  
  bool doForwardBackward;
  istringstream fbStringConverter(params["-fb"]);
  if(!(fbStringConverter >> doForwardBackward)){
    cout << "Invalid forward backward flag, should be 0 or 1" << endl;
    exit(1);
  }
  
  bool coFactor;
  istringstream coStringConverter(params["-co"]);
  if(!(coStringConverter >> coFactor)){
    cout << "Invalid co flag, should be 1 or 0" << endl;
    exit(1);
  }
  
  int genSinceAdmixture;
  istringstream genStringConverter(params["-G"]);
  if(!(genStringConverter >> genSinceAdmixture)){
    cout << "Invalid generations since admixture" << endl;
    exit(1);
  }
  
  int minNodeSize;
  istringstream mnsStringConverter(params["-n"]);
  if(!(mnsStringConverter >> minNodeSize)){
    cout << "Invalid minimum node size" << endl;
    exit(1);
  }
  
  cout << "Minimum node size: " << minNodeSize << endl;
  
  bool allowParallel;
  istringstream parallelStringConverter(params["-r"]);
  if(!(parallelStringConverter >> allowParallel)){
    cout << "Invalid parelleization flag, should be 0 or 1" << endl;
    exit(1);
  }
  
  if(!allowParallel){
    omp_set_num_threads(1);
  }
  
  int numEmIterations;
  istringstream emStringConverter(params["-e"]);
  if(!(emStringConverter >> numEmIterations)){
    cout << "Invalid number of iterations" << endl;
    exit(1);
  }
  
  int callAnc;
  istringstream callAncStringConverter(params["-c"]);
  if(!(callAncStringConverter >> callAnc)){
    cout << "Invalid -c argument" << endl;
    exit(1);
  }
  
  int useAnc;
  istringstream useAncStringConverter(params["-u"]);
  if(!(useAncStringConverter >> useAnc)){
    cout << "Invalid -u argument" << endl;
    exit(1);
  }
  
  (*logfile).open((params["-o"]+".log.txt").c_str());
  (*logfile) << "Parameters fetched" << endl;
  
  // ===== Read in SNPs to exclude. 0 indexed =====
  int numPruned = 0;
  vector<int> excludeIndices;
  if(params["-x"].length() != 0){
    ifstream excludes;  // assume one column
    excludes.open(params["-x"].c_str());
    
    if(!excludes.is_open()){
      cout<<"Error: Unable to open the pruned markers file: "<<endl;
      cout<<params["-x"]<<endl;
      
      (*logfile)<<"Error: Unable to open the pruned markers file: "<<endl;
      (*logfile)<<params["-x"]<<endl;
      
      exit(1);
    }
    
    string excludeLine;
    int excludeIndex;
    while(!excludes.eof()){
      getline(excludes, excludeLine);
      
      if(excludeLine.length() == 0){
        continue;
      }
      
      istringstream isse(excludeLine);  // slow
      if(!(isse >> excludeIndex)){
        cout<<"Error: Exclude index not an integer: "<<endl;
        cout<<excludeIndex<<endl;
        
        (*logfile)<<"Error: Exclude index not an integer: "<<endl;
        cout<<excludeIndex<<endl;
        
        exit(1);
      }
      
      if(excludeIndices.size() > 0 && excludeIndices.back() > excludeIndex){
        cout << "Exclude indices not in order" << endl;
        exit(1);
      }
      
      excludeIndices.push_back(excludeIndex);
    }
    
    numPruned = (int) excludeIndices.size();
    
    excludes.close();
  }
  
  // ===== Read in SNP positions in cM =====
  // ===== Use this to get index of beginning of each window and distances b/w windows =====  
  ifstream snpPos;  // assume one column
  snpPos.open(params["-m"].c_str());
    
  if (!snpPos.is_open() ) {
    cout<<"Error: Unable to open the marker locations file: "<<endl;
    cout<<params["-m"]<<endl;
    
    (*logfile)<<"Error: Unable to open the marker locations file: "<<endl;
    (*logfile)<<params["-m"]<<endl;
    
    exit(1);
  }
  
  //Read in file and get number of SNPs and windows
  string locLine;
  double markerLoc;
  vector<double> * windowBeginLocs = new vector<double>(); //ith element has location in cM of first snp in window i
  windowBeginLocs->reserve(1000);   // 1000 windows is probably max ever needed
  vector<double> * windowEndLocs = new vector<double>();  // used for determining distances between SNPs
  windowEndLocs->reserve(1000);
  vector<int> * windowBeginIndexes = new vector<int>();  // used for determining which SNP is in which window
  windowBeginIndexes->reserve(1000);
  vector<int> * windowEndIndexes = new vector<int>();
  windowEndIndexes->reserve(1000);
  bool firstSNP = true;
  double windowBegin = 0.0;
  double lastLoc = 0.0;
  int forPruneLineIndex = 0;
  int lineIndex = 0;
  double windowSize = atof((params["-w"]).c_str());
  unsigned int pruneIndex = 0;
  while(!snpPos.eof()){
    getline(snpPos, locLine);
    
    if(pruneIndex < excludeIndices.size() && excludeIndices[pruneIndex] == forPruneLineIndex){
      pruneIndex++;
      forPruneLineIndex++;
      continue;
    }
    
    if(locLine.length() == 0){
      continue;
    }
    
    istringstream iss(locLine);  // slow
    if(!(iss >> markerLoc)){
      cout<<"Error: Marker location not a decimal number: "<<endl;
      cout<<markerLoc<<endl;
      
      (*logfile)<<"Error: Marker location not a decimal number: "<<endl;
      cout<<markerLoc<<endl;
      
      exit(1);
    }
        
    if(firstSNP){
      firstSNP = false;
      windowBeginLocs->push_back(markerLoc);
      windowBeginIndexes->push_back(lineIndex);
      windowBegin = markerLoc;
    } else if((markerLoc - windowBegin) > windowSize){
      windowEndLocs->push_back(lastLoc);
      windowEndIndexes->push_back(lineIndex-1);
      windowBeginLocs->push_back(markerLoc);
      windowBeginIndexes->push_back(lineIndex);
      windowBegin = markerLoc;
    }
    
    lastLoc = markerLoc;
    lineIndex++;
    forPruneLineIndex++;
  }
  
  windowEndLocs->push_back(lastLoc);
  windowEndIndexes->push_back(lineIndex-1);
  
  int numSNPs = lineIndex;  // not -1 because lineIndex starts at 0
  cout << "Number of SNPs Read: " << forPruneLineIndex << endl;
  cout << "Number of SNPs to Exclude: " << numPruned << endl;
  cout << "Number left after pruning: " << numSNPs << endl;
  int numWindows = (int) windowBeginLocs->size();
  cout << "Number of windows: " << numWindows << endl;
    
  snpPos.close();
  
  // Read in classes, get number of haplotypes
  ifstream classes;  // assume one row
  classes.open(params["-p"].c_str());
  
  if (!classes.is_open() ) {
    cout<<"Error: Unable to open the classes file: "<<endl;
    cout<<params["-p"]<<endl;
    
    (*logfile)<<"Error: Unable to open the classes file: "<<endl;
    (*logfile)<<params["-p"]<<endl;
    
    exit(1);
  }
  
  string classLine;
  getline(classes, classLine);
  istringstream classSS(classLine);
  int hapClass;
  vector<int> * hapClasses = new vector<int>();  // holds classes of nonAdm haps
  int numAncPops = 0;  // Assumes 0 is admixed label and ancestries are 1,2,...,numAncPops labeled with no skips
  int numAdm = 0;
  while(classSS >> hapClass){
    hapClasses->push_back(hapClass);
    
    if(hapClass > numAncPops){
      numAncPops = hapClass;
    }
       
    if(hapClass == 0){
      numAdm++;
    }
  }
  
  classes.close();
  
  int numHaps = (int) hapClasses->size();
  int numAnc = numHaps - numAdm;  // for initial round we need this, later numHaps=numAnc=numAdm
  if(useAnc == 1){
    numAdm = numHaps;    // EM
  }
  
  map<string, string>::iterator bootIt = params.find("-b");
  if(bootIt == params.end()){
    // not specified by user
    stringstream ssboot;
    ssboot << numAnc;  // default bootstrap sample size
    params["-b"] = ssboot.str();
  }
  int bootstrapSampleSize;
  istringstream bootStringConverter(params["-b"]);
  if(!(bootStringConverter >> bootstrapSampleSize)){
    cout << "Invalid bootstrap sample size" << endl;
    exit(1);
  }
    
  //Read in alleles
  ifstream alleles;
  alleles.open(params["-a"].c_str());
  
  if (!alleles.is_open() ) {
    cout<<"Error: Unable to open the alleles file: "<<endl;
    cout<<params["-a"]<<endl;
    
    (*logfile)<<"Error: Unable to open the alleles file: "<<endl;
    (*logfile)<<params["-a"]<<endl;
    
    exit(1);
  }
  
  // each window is a single matrix to keep everything together in memory. Use multiplication to get different snp lines
  // Keep this the same for EM because no calls on Adm yet
  vector<bool> * ancWindows = new vector<bool>[numWindows];
  for(int windex = 0; windex < numWindows; windex++){
    int numSNPsInWindow = (*windowEndIndexes)[windex] - (*windowBeginIndexes)[windex] + 1;
    int windowBlockSize = numSNPsInWindow * numAnc;
    ancWindows[windex] = vector<bool>(windowBlockSize, false);
  }
  
  // EM: nothing changed because changed numAdm above
  vector<bool> * admWindows = new vector<bool>[numWindows];
  for(int windex = 0; windex < numWindows; windex++){
    int numSNPsInWindow = (*windowEndIndexes)[windex] - (*windowBeginIndexes)[windex] + 1;
    int windowBlockSize = numSNPsInWindow * numAdm;
    admWindows[windex] = vector<bool>(windowBlockSize, false);
  }
  
  int* ancHapClasses = new int[numAnc];
  
  string snpLineString;
  int snpIndex = 0;
  int windowIndex = 0;
  int numSNPsInWindow = (*windowEndIndexes)[windowIndex] - (*windowBeginIndexes)[windowIndex] + 1;
  int alleleLineIndex = 0;
  unsigned int allelePruneIndex = 0;
  while(!alleles.eof()){
    getline(alleles, snpLineString);
    
    if(allelePruneIndex < excludeIndices.size() && excludeIndices[allelePruneIndex] == alleleLineIndex){
      allelePruneIndex++;
      alleleLineIndex++;
      continue;
    }
    
    if(snpLineString.length() == 0){
      continue;
    }
    
    int admHapIndex = 0;
    int ancHapIndex = 0;
    
    for(int hapIndex = 0; hapIndex < numHaps; hapIndex++){
      if((*hapClasses)[hapIndex] != 0){
        ancWindows[windowIndex][snpIndex*numAnc + ancHapIndex] = snpLineString[hapIndex] - '0';
        ancHapClasses[ancHapIndex] = (*hapClasses)[hapIndex] - 1;  //make 0 based
        ancHapIndex++;
        if(useAnc == 1){
          // EM
          admWindows[windowIndex][snpIndex*numAdm + admHapIndex] = snpLineString[hapIndex] - '0';
          admHapIndex++;
        }
      }else{
        // EM
        admWindows[windowIndex][snpIndex*numAdm + admHapIndex] = snpLineString[hapIndex] - '0';
        admHapIndex++;
      }
    }
    
    snpIndex++;
    
    if(snpIndex == numSNPsInWindow){
      snpIndex = 0;
      windowIndex++;
      numSNPsInWindow = (*windowEndIndexes)[windowIndex] - (*windowBeginIndexes)[windowIndex] + 1;
    }
    
    alleleLineIndex++;
  }
  
  // EM
  int** ancHapClassesPerWindow = new int*[numWindows];
  for(int windex = 0; windex < numWindows; windex++){
    ancHapClassesPerWindow[windex] = new int[numAnc];
    for(int i = 0; i < numAnc; i++){
      ancHapClassesPerWindow[windex][i] = ancHapClasses[i];
    }
  }
  
  alleles.close();
  
  delete[] ancHapClasses;
  
  // Phasing
  // Stores number of het sites in each admixed individual per window
  cout << "Creating phasings" << endl;
  int numAdmInds = numAdm / 2;
  int** numHetSitesPerInd = new int*[numWindows];    // numWindows x numAdmInds
  for(int windex = 0; windex < numWindows; windex++){
    numHetSitesPerInd[windex] = new int[numAdmInds];
    int numSNPsInWindow = (*windowEndIndexes)[windex] - (*windowBeginIndexes)[windex] + 1;
    for(int admInd = 0; admInd < numAdmInds; admInd++){
      numHetSitesPerInd[windex][admInd] = 0;
    }
    for(int snpIndex = 0; snpIndex < numSNPsInWindow; snpIndex++){
      for(int admInd = 0; admInd < numAdmInds; admInd++){
        int allele1 = admWindows[windex][snpIndex*numAdm + 2*admInd];
        int allele2 = admWindows[windex][snpIndex*numAdm + 2*admInd+1];
        
        if(allele1 != allele2){
          numHetSitesPerInd[windex][admInd] += 1;
        }
      }
    }
  }
  
  // Figure out dimensions for array containing the different possible phasings of each admixed individual
  // Assume 0 or 1 switches per window
  int * numAdmPhasings = new int[numWindows];
  for(int windex = 0; windex < numWindows; windex++){
    numAdmPhasings[windex] = 0;
    for(int admInd = 0; admInd < numAdmInds; admInd++){
      numAdmPhasings[windex] += numHetSitesPerInd[windex][admInd] + 1;    // original phasing plus one switch per het site
    }
  }
  
  // Stores all phasings
  // For each person, phasing go: original, last switched, 2nd to last and after switched, etc.
  vector<bool> * admWindowsPhasings = new vector<bool>[numWindows];
  for(int windex = 0; windex < numWindows; windex++){
    int numSNPsInWindow = (*windowEndIndexes)[windex] - (*windowBeginIndexes)[windex] + 1;
    int windowBlockSize = numSNPsInWindow * numAdmPhasings[windex]*2;
    admWindowsPhasings[windex] = vector<bool>(windowBlockSize, false);
  }
  
  for(int windex = 0; windex < numWindows; windex++){
    int start = 0;
    int numSNPsInWindow = (*windowEndIndexes)[windex] - (*windowBeginIndexes)[windex] + 1;
    int numHapsInWindow = numAdmPhasings[windex]*2;
    for(int admInd = 0; admInd < numAdmInds; admInd++){
      int indNumHetSites = numHetSitesPerInd[windex][admInd];
      int indNumHaps = 2*(indNumHetSites+1);
      int indNumPhasings = indNumHetSites+1;
      
      int hetIndex = 0;
      // go bottom to top
      for(int snpIndex = numSNPsInWindow-1; snpIndex >= 0; snpIndex--){
        int allele1 = admWindows[windex][snpIndex*numAdm + 2*admInd];
        int allele2 = admWindows[windex][snpIndex*numAdm + 2*admInd + 1];
        if (allele1 != allele2){
          hetIndex += 1;
        }
        for(int phasing = 0; phasing < indNumPhasings; phasing++){
          // if phasing < hetIndex just copy; otherwise switch; do for hom sites too bc doesn't matter
          int col1 = start + 2*phasing;
          int col2 = start + 2*phasing + 1;
          
          if(phasing < hetIndex){
            admWindowsPhasings[windex][snpIndex*numHapsInWindow + col1] = allele1;
            admWindowsPhasings[windex][snpIndex*numHapsInWindow + col2] = allele2;
          } else{
            admWindowsPhasings[windex][snpIndex*numHapsInWindow + col1] = allele2;
            admWindowsPhasings[windex][snpIndex*numHapsInWindow + col2] = allele1;
          }
        }
      }
      
      start = start + indNumHaps;
    }
  }

  // DEBUG
//  for(int windex = 0; windex < 5; windex++){
//    int numSNPsInWindow = (*windowEndIndexes)[windex] - (*windowBeginIndexes)[windex] + 1;
//    int numHapsInWindow = numAdmPhasings[windex]*2;
//    for(int snpIndex = 0; snpIndex < numSNPsInWindow; snpIndex++){
//      for(int hapIndex = 0; hapIndex < numHapsInWindow; hapIndex++){
//        cout << admWindowsPhasings[windex][snpIndex*numHapsInWindow + hapIndex] << " ";
//      }
//      cout << endl;
//    }
//  }
  
  cout << "Done creating phasings" << endl;
  //DEBUG
//  for(int windex = 0; windex < 5; windex++){
//    for(int admInd = 0; admInd < numAdmInds; admInd++){
//      cout << numHetSitesPerInd[windex][admInd] << " ";
//    }
//    cout << "\n" << endl;
//  }

  
  //Step 9: Store the data
  ProcessedInput* processedInput = new ProcessedInput();
  processedInput->numSNPs = numSNPs;
  processedInput->numWindows = (int) numWindows;
  processedInput->numHaps = numHaps;
  processedInput->numAncPops = numAncPops;
  processedInput->numAdm = numAdm;
  processedInput->numAnc = numAnc;
  processedInput->windowSize = windowSize;
  processedInput->numTrees = numTrees;
  processedInput->genSinceAdmixture = genSinceAdmixture;
  processedInput->windowBeginLocs = windowBeginLocs;
  processedInput->windowEndLocs = windowEndLocs;
  processedInput->windowBeginIndexes = windowBeginIndexes;
  processedInput->windowEndIndexes = windowEndIndexes;
  processedInput->ancHapClassesPerWindow = ancHapClassesPerWindow;
  processedInput->ancWindows = ancWindows;
  processedInput->admWindows = admWindows;
  processedInput->outputName = params["-o"];
  processedInput->allowParallel = allowParallel;
  processedInput->bootstrapSampleSize = bootstrapSampleSize;
  processedInput->bootstrappingMethod = bootstrappingMethod;
  processedInput->numEmIterations = numEmIterations;
  processedInput->hapClasses = hapClasses;
  processedInput->callAnc = callAnc;
  processedInput->useAnc = useAnc;
  processedInput->admWindowsPhasings = admWindowsPhasings;
  processedInput->numHetSitesPerInd = numHetSitesPerInd;
  processedInput->numAdmPhasings = numAdmPhasings;
  processedInput->minNodeSize = minNodeSize;
  processedInput->doForwardBackward = doForwardBackward;
  processedInput->compressedOutput = coFactor;
  
  return processedInput;
}
