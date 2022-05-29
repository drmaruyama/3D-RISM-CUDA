#include <iostream>
#include "rism3d.h"

void RISM3D :: initialize(string control, string structure) {
  read_input(control, structure);
  set_cuda();
  set_fname(control, structure);
  initialize_g();
  set_solvent();
  cal_potential();
} 
