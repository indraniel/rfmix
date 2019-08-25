import argparse
import os
import sys

parser = argparse.ArgumentParser()

# Required arguments
parser.add_argument("admixed_phasing", help="Specifies if admixed samples were trio-phased (TrioPhased) or population-phased (PopPhased). PopPhased has RFMix simultaneously attempt to correct phasing errors. Use if admixed samples are population-phased. PopPhased significantly slows analysis. TrioPhased-assume trio-phased samples, PopPhased-assume population-phased samples", type=str)
parser.add_argument("alleles", help="file containing binary-coded alleles", type=file)
parser.add_argument("classes", help="file containing ancestry of haplotypes in alleles file", type=file)
parser.add_argument("snp_locations", help="file containing genetic coordinates of snps in alleles file", type=file)

# Optional arguments
parser.add_argument("-w", "--window-size", help="window size in genetic distance. Default=0.2cM", type=float)
parser.add_argument("-G", "--generations", help="generations since admixture event. Default=8", type=int)
parser.add_argument("-o", "--output-name", help="prefix for output file. Default=ancestry_output")
parser.add_argument("-t", "--trees", help="number of trees to generate per random forest. Default=100", type=int)
parser.add_argument("--disable-parallel", help="turn off OpenMP (shared memory) parallelization", action="store_true")
parser.add_argument("--num-threads", help="number of threads to use if doing parallelization. Default is number of processors on computer/node", type=int)
parser.add_argument("-b", "--bootstrap-sample-size", help="number of bootstrap samples used per tree. Default=<number of reference haplotypes>", type=int)
parser.add_argument("-s", "--bootstrap-type", help="kind of bootstrap sampling. 0-sample with equal prob from each haplotype, 1-(Default)-sample with equal probability from each class, 2-stratified by class", type=int, choices=[0,1,2])
parser.add_argument("-e", "--em-iterations", help="number of EM iterations. Default=0", type=int)
parser.add_argument("--use-reference-panels-in-EM", help="not using this flag means that reference panels are discarded after the initial inference step", action="store_true")
parser.add_argument("-x", "--discard-rare-variants", help="ignore SNPs with minor allele occurrence less than the given value e.g. an argument of 2 will ignore SNPs with minor alleles that occur 1 or less times across both the admixed and reference panels. Outputs file with suffic toExclude containing 0-based indices of removed SNPs. Default=0", type=int)
parser.add_argument("-f", "--mtry-factor", help="mtry is the number of SNPs to consider when creating a tree node in a random forest. The standard is to use the square root of the number of SNPs in the window. This multipiplies that number", type=float)
parser.add_argument("--forward-backward", help="output the forward-backward probabilities. Defaults to not doing so because of potentially large file sizes", action="store_true")
parser.add_argument("--skip-check-input-format", help="don't do initial check of input files", action="store_true")
parser.add_argument("-n", "--min-node-size", help="minimum number of reference haplotypes per tree node. Default=1", type=int)
parser.add_argument("--succinct-output", help="output calls per window along with file listing number of SNPs per window", action="store_true")

args = parser.parse_args()

if args.admixed_phasing not in ["TrioPhased", "PopPhased"]:
	print "First argument must specify TrioPhased or PopPhased"
	sys.exit()

# Input checking
if not args.skip_check_input_format:
  print "Checking input..."
  # Get number of SNPs in locations file
  locations_count = 0
  for lineIndex,line in enumerate(args.snp_locations):
    if len(line.split()) == 0:
      continue
    locations_count += 1
      
  # Check if alleles file is binary
  alleles_count = 0
  for lineIndex,line in enumerate(args.alleles):
    if len(line.split()) == 0:
      continue

    alleles_count += 1
    
    lineStrip = line.strip()
    
    if lineIndex == 0:
      num_haps_in_alleles_file = len(lineStrip)
    
    for allele in lineStrip:
      if int(allele) not in [0,1]:
        raise argparse.ArgumentTypeError("Allele " + allele + " on line " + str(lineIndex) + " in alleles file is not 0 or 1")

  # See if number of SNPs in locations and alleles files agree
  if locations_count != alleles_count:
    raise argparse.ArgumentTypeError("Number of snps in marker locations file (" + str(locations_count) + ") does not equal number of snps in alleles file (" + str(alleles_count) + ")")

  # See if classes file assigns ancestry indices appropriately
  maxAnc = 0
  classes = [int(x) for x in args.classes.readline().split()]
  for hapClass in classes:
    if hapClass > maxAnc:
      maxAnc = hapClass
  classUsed = [False]*(maxAnc+1)
  for ancClass in classes:
    classUsed[ancClass] = True
  for classIndex in range(maxAnc+1):
    if not classUsed[classIndex]:
      raise argparse.ArgumentTypeError("Class " + str(classIndex) + " not present in classes file")
  print "Done checking input"

# Look for non/rare variant sites to exclude
if args.discard_rare_variants:
  excludeFileName = "toExclude.txt"
  if args.output_name:
    excludeFileName = args.output_name + excludeFileName
  else:
    excludeFileName = "ancestry_output" + excludeFileName
  excludeFile = open(excludeFileName, 'w')
  minCount = args.discard_rare_variants
  for lineIndex,line in enumerate(args.alleles):
    if len(line.split()) == 0:
      continue
    if lineIndex == 0:
      numHaps = len(line.strip())
    zeroCount = 0
    for allele in line.strip():
      if allele == "0":
        zeroCount += 1
    if zeroCount < minCount or numHaps-zeroCount < minCount:
      excludeFile.write(str(lineIndex) + "\n")
  excludeFile.close()

# Run RFMix
argDict = {"window_size":"-w", "generations":"-G", "output_name": "-o", "trees": "-t", "disable_parallel":"-r", "num_threads":"-h", "bootstrap_sample_size":"-b", "bootstrap_type":"-s", "em_iterations":"-e", "use_reference_panels_in_EM":"-u", "discard_rare_variants":"-x", "mtry_factor":"-f", "forward_backward": "-fb", "min_node_size": "-n", "succinct_output": "-co"}

parameters = ["-a", args.alleles.name, "-p", args.classes.name, "-m", args.snp_locations.name]

opts = vars(args)
options = {k : opts[k] for k in opts if k in argDict.keys() and opts[k] not in [None, False]}
for option in options.keys():
  parameters.append(argDict[option])
  if option == "disable_parallel":
    parameters.append("0")
  elif option == "use_reference_panels_in_EM":
    parameters.append("1")
  elif option == "discard_rare_variants":
    parameters.append(excludeFileName)
  elif option == "forward_backward":
    parameters.append("1")
  elif option == "succinct_output":
    parameters.append("1")
  else:
    parameters.append(str(options[option]))

if args. admixed_phasing == "TrioPhased":
  command = " ".join(["./TrioPhased/RFMix_TrioPhased"] + parameters)
  print command
  os.system(command)
elif args.admixed_phasing == "PopPhased":
  command = " ".join(["./PopPhased/RFMix_PopPhased"] + parameters)
  print command
  os.system(command)
  

