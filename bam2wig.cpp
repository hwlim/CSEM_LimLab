#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<string>
#include<vector>
#include "wiggle.h"

using namespace std;

vector<string> args;

void printUsage() {
  printf("Usage: csem-bam2wig sorted_bam_input wig_output wiggle_name [--no-fractional-weight] [--extend-reads fragment_length] [--only-midpoint] [--help]\n");
  printf("sorted_bam_input\t\t: Input BAM format file, must be sorted\n");
  printf("wig_output\t\t\t: Output wiggle file's name, e.g. output.wig\n");
  printf("wiggle_name\t\t\t: The name of this wiggle plot\n");
  printf("--no-fractional-weight\t\t: If this is set, RSEM will not look for \"ZW\" tag and each alignment appeared in the BAM file has weight 1. Set this if your BAM file is not generated by CSEM\n");
  printf("--extend-reads fragment_length\t: Extend reads to their full fragment length. fragment_length is the average fragment length for the data set and should be positive\n");
  printf("--only-midpoint\t\t\t: represent each fragment by its midpoint. This can be set only if --extend-reads is set\n");
  printf("--help\t\t\t\t: Show help information\n");
  exit(-1);
}

int main(int argc, char* argv[]) {

  args.clear();
  // skip argv[0], which is the program's name
  for (int i = 1; i < argc; i++) {
    if (!strncmp(argv[i], "--", 2)) {
      if (!strcmp(argv[i], "--no-fractional-weight")) no_fractional_weight = true;
      else if (!strcmp(argv[i], "--extend-reads")) { 
	if (i + 1 == argc || (fragment_length = atoi(argv[i + 1])) <= 0) { printf("--extend-reads option is not set correctly!\n"); printUsage(); }
	++i; // change i, to skip fragment_length
      }
      else if (!strcmp(argv[i], "--only-midpoint")) only_midpoint = true;
      else if (!strcmp(argv[i], "--help")) printUsage();
      else { printf("Cannot recognize option \"%s\"!\n", argv[i]); printUsage(); }
    }
    else args.push_back(argv[i]);
  }

  if (only_midpoint && fragment_length <= 0) { printf("--only-midpoint cannot be set if --extend-reads is not set!\n"); printUsage(); }

  if (args.size() != 3) { printf("Number of arguments does not match!\n"); printUsage(); } 

  UCSCWiggleTrackWriter track_writer(args[1], args[2]);
  build_wiggles(args[0], track_writer);
  
  return 0;
}
