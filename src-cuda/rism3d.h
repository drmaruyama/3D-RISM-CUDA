#ifndef RISM3D_H
#define RISM3D_H

#include <complex>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <valarray>
#include <vector>
#include "physical.h"
#include "cell.h"
#include "control.h"
#include "solute.h"
#include "solvent.h"
#include "anderson.h"
#include "fft3d.h"
using namespace std;

class RISM3D {
public:
  RISM3D () {ce = new Cell; co = new Control; su = new Solute;
    sv = new Solvent; ma = new AN2; fft = new FFT3D;}
  ~RISM3D () {delete ce, co, su, sv;} 
  void initialize (char[]);
  void iterate (int);    
  void output ();
private:
  void add_tuv (double);
  void cal_Coulomb ();
  void cal_exchem (double * &);
  void cal_grad (double * &);
  void cal_LJ ();
  double cal_pmv ();
  void cal_potential ();
  double cal_rms ();
  void calculate (double);
  void initialize_g ();
  void initialize_tuv ();
  void output_cuv ();
  void output_grad (double * &);
  void output_guv ();
  void output_huv ();
  void output_xmu (double * &, double);
  void read_input (char[]);
  void read_tuv ();
  void set_fname (char[]);
  void set_cuda ();
  void set_solvent ();
  void write_tuv ();

  // Host
  vector <complex <double> *> guv;
  vector <complex <double> *> huv;
  vector <double *> tuv;
  vector <double> ga;
  double * siguv;
  double * epsuv;
  int * indga;
  int clos;
  int nga;
  string outlist;
  string fsolvent;
  string fname;
  dim3 g, b;

  // Device
  double2 * dguv;
  double2 * dhuv;
  double * dt;
  double * dtr;
  double * du;
  double * de;
  double4 * dgv;
  double * dsig;
  double * deps;
  double * dx;
  double * dfr;
  double2 * dfk;
  double * ds;

  Cell * ce;
  Control * co;
  Solute * su;
  Solvent * sv;
  AN2 * ma;
  FFT3D * fft;
};

#endif
