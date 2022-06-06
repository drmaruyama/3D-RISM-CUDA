#include <stdlib.h>
#include <iostream>
#include "cell.h"
using namespace std;

void Cell :: setup() {
  volume = box[0] * box[1] * box[2];
  ngrid = grid[0] * grid[1] * grid[2];
  dv = volume / ngrid;
  dr[0] = box[0] / grid[0];
  dr[1] = box[1] / grid[1];
  dr[2] = box[2] / grid[2];
  shift[0] = shift[1] = shift[2] = 0.0;
}
