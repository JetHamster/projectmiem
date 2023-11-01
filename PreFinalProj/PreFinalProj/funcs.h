#pragma once
#include <iostream>
#include <math.h>
using namespace std;
#include <limits>
#include <vector>
#include <fstream>
#include <iomanip>
#include <algorithm>
#define PI 3.14
typedef std::numeric_limits< double > dbl;

bool Check(vector<double> h);
double Hknorm(vector<double> h);

void write(const char* resfile, vector<double> rezvec);
void clear(const char* resfile);

void showvec(vector<double> v);
void showmatr(vector<vector<double>> A);

vector<double> msolve(vector<double> c, vector<double> d, vector<double> e, vector<double> f);

vector<double> Deriv(vector <double> vec, vector<double> grid);

vector<double> summ(vector<double> f1, vector<double> f2);
vector<double> diff(vector<double> f1, vector<double> f2);
vector<double> mult(vector<double> f1, double k);

vector<double> makesimplegrid(int size, double l, double r);
vector<double> makeevengrid(int size, double l, double r);
vector<double> makegrid(int size, double l, double r, double epsilon);


double l2norm(vector<double> vec);



double scalar(vector<double> v1, vector<double> v2);

vector<double> zeros(int n);
vector<double> Adotv(vector<vector<double>> A, vector<double> v);

vector<double> cg(vector<vector<double>> A, vector<double> b, vector<double> xprev);
vector<vector<double>> makematr(vector<double> c, vector<double> d, vector<double> e);
vector<double> linpsi(vector<double> grid);
vector<double> mypsi(vector<double> grid, double gamma);


double min(double a, double b);
double dtest(double gamma);

void dets(vector<double> c, vector<double> d, vector<double> e);