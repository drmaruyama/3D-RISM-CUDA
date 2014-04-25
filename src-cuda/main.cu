#include <iostream>
#include <fstream>
#include <unistd.h>
#include "rism3d.h"

int main (int argc, char * argv[]) {
  RISM3D * system;
  int ch;
  int cu, dn;

  cu = dn = 0;
  system = new RISM3D;

  while ((ch = getopt(argc, argv, "c:d:")) != -1) {
    switch (ch){
    case 'c':
      cu = atoi(optarg);
      break;
    case 'd':
      dn = atoi(optarg);
      break;
    }
  }
  if (argc == 1) {
    cout << "No parameter file!" << endl ;
    return (1) ;
  }

  cout << "Set device " << dn << endl ;
  cudaSetDevice(dn);
  if (cu > 0) cout << "Charge up " << cu << endl;
  system -> initialize(argv[optind]);
  system -> iterate(cu);
  system -> output();    

  return(0);
}
