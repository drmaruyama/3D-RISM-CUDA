#include <algorithm>
#include <iostream>
#include <fstream>
#include "rism3d.h"

void RISM3D :: output() {

  transform(outlist.begin(), outlist.end(), outlist.begin(), ::tolower);

  int flag = 0;
  
  if (outlist.find("m") != string::npos) flag += 0x01;
  if (outlist.find("d") != string::npos) flag += 0x02;
  if (outlist.find("c") != string::npos) flag += 0x04;
  if (outlist.find("g") != string::npos) flag += 0x08;
  if (outlist.find("h") != string::npos) flag += 0x10;

  if ((flag & 0x01) == 0x01) {
    double pmv = cal_pmv();
    double * xmu = new double[sv -> natv * 2];
    cal_exchem(xmu);
    output_xmu(xmu, pmv);
    delete[] xmu;
  }

  if ((flag & 0x02) == 0x02) {
    double * du;
    du = new double[su -> num * 3];
    cal_grad(du);
    output_grad(du);
    delete[] du;
  }

  if ((flag & 0x04) == 0x04) {
    output_cuv();
  }

  if ((flag & 0x08) == 0x08) {
    output_guv();
  }

  if ((flag & 0x10) == 0x10) {
    output_huv();
  }
} 
