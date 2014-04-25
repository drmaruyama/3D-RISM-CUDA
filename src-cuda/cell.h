#ifndef Cell_H
#define Cell_H

class Cell {
public:
  Cell() {box = new double[3]; dr = new double[3]; grid = new int[3];}
  ~Cell() {delete[] box, dr, grid;}
  void setup();
  double * box;
  double * dr;
  int * grid;
  double volume, dv;
  int ngrid;
};

#endif // Cell_H
