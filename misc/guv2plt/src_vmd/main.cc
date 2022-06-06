#include <iostream>
#include <fstream>
#include <string>
#include <unistd.h>
using namespace std;

#define MAGICNUMBER 3
#define TYPEOS 200

int main(int argc, char * argv[]) {
  string set_fname(char[]);

  int ch;
  float xs, ys, zs;
  
  while ((ch = getopt(argc, argv, "x:y:z:")) != -1) {
    switch (ch){
    case 'x':
      xs = atof(optarg);
      break;
    case 'y':
      ys = atof(optarg);
      break;
    case 'z':
      zs = atof(optarg);
      break;
    }
  }
  if (argc == 1) {
    cout << "No parameter file!" << endl ;
    return (1) ;
  }

  char* filename = argv[optind];
  std::string fname = set_fname(filename);
  ifstream fin(filename, ios::in | ios::binary);

  int nv;
  int g[3];
  double b[3];
  double s[3];
  fin.read((char *) &nv, sizeof(int));
  fin.read((char *) &g[0], sizeof(int)*3);
  fin.read((char *) &b[0], sizeof(double)*3);
  fin.read((char *) &s[0], sizeof(double)*3);

  int grid = g[0] * g[1] * g[2];
  double * guv = new double[grid * nv];
  fin.read((char *) &guv[0], sizeof(double) * grid * nv);
  fin.close();

  float xmin = - b[0] / 2 + xs - s[0];
  float ymin = - b[1] / 2 + ys - s[1];
  float zmin = - b[2] / 2 + zs - s[2];
  float xmax = b[0] / 2 - b[0] / g[0] + xs - s[0];
  float ymax = b[1] / 2 - b[1] / g[1] + ys - s[1];
  float zmax = b[2] / 2 - b[2] / g[2] + zs - s[2];

  int rank = MAGICNUMBER;
  int TypeOfSurface = TYPEOS;
  for (int iv = 0; iv < nv; ++iv) {
    char nplt[10];
    sprintf(nplt, "-%d.plt", iv);
    ofstream fout((fname + nplt).c_str(), ios::out | ios::binary);
    fout.write((char *) &rank, sizeof(int));
    fout.write((char *) &TypeOfSurface, sizeof(int));
    fout.write((char *) &g[2], sizeof(int));
    fout.write((char *) &g[1], sizeof(int));
    fout.write((char *) &g[0], sizeof(int));
    fout.write((char *) &zmin, sizeof(float));
    fout.write((char *) &zmax, sizeof(float));
    fout.write((char *) &ymin, sizeof(float));
    fout.write((char *) &ymax, sizeof(float));
    fout.write((char *) &xmin, sizeof(float));
    fout.write((char *) &xmax, sizeof(float));
    int ip = grid * iv;
    for (int ig = 0; ig < grid; ++ig) {
      float tmp = guv[ig + ip];
      fout.write((char *) &tmp, sizeof(float));      
    }
    fout.close();
  }
  return 0;
}
