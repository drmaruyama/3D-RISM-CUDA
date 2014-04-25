#include <iostream>
#include "rism3d.h"

void RISM3D :: initialize(char inputfile[]) {

  read_input(inputfile);
  set_cuda();
  set_fname(inputfile);
  initialize_g();
  set_solvent();
  cal_potential();
} 
