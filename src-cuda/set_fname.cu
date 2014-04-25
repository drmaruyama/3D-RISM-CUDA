#include <iostream>
#include <fstream>
#include "rism3d.h"

void RISM3D :: set_fname (char inputfile[]) {
  fname.append(inputfile);
  fname = fname.substr(0, fname.rfind("."));
}
