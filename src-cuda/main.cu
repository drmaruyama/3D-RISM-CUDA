#include <iostream>
#include <fstream>
#include <unistd.h>
#include "rism3d.h"

int main (int argc, char * argv[]) {
  RISM3D * system;
  int ch;
  int cu, dn;
  string input;
  string structure;
  string esp;
  bool centering = true;

  cu = dn = 0;
  system = new RISM3D;

  while ((ch = getopt(argc, argv, "c:d:i:s:e:f")) != -1) {
    switch (ch){
    case 'c':
      cu = atoi(optarg);
      break;
    case 'd':
      dn = atoi(optarg);
      break;
    case 'i':
      input = optarg;
      break;
    case 's':
      structure = optarg;
      break;
    case 'e':
      esp = optarg;
    case 'f':
      centering = false;
      break;
    }
  }

  if (input.empty() || structure.empty()) {
    if (argv[optind] == NULL) {
      cout << "No input file!" << endl;
      return (1);
    }
    input = argv[optind];
  }

  cout << "Set device " << dn << endl ;
  cudaSetDevice(dn);
  if (cu > 0) cout << "Charge up " << cu << endl;
  system -> initialize(input, structure, esp, centering);
  system -> iterate(cu);
  system -> output();    

  return(0);
}
