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

bool Check(vector<double> h) {
	auto it = max_element(std::begin(h), std::end(h));
	//cout <<  setprecision(20) << *it << " k" << endl;
	if (*it < 1.0e-9)
		return false;
	else
		return true;
}

double Hknorm(vector<double> h) {
	auto it = max_element(std::begin(h), std::end(h));
	//cout <<  setprecision(20) << *it << " k" << endl;
	return *it;
}

void write(const char* resfile, vector<double> rezvec) {
	ofstream f;

	f.open(resfile, ios::app);
	//вводим количество вещественных чисел
	for (int i = 0; i < rezvec.size(); i++)
	{
		f << rezvec[i] << " ";
	}
	f << endl;
	f.close();
}
void clear(const char* resfile) {
	ofstream f;

	f.open(resfile, ios::out);
	f << "";
	f.close();
}


void showvec(vector<double> v) {
	for (int i = 0; i < v.size(); i++)
		cout << v[i] << " ";
	cout << endl << endl;
}

void showmatr(vector<vector<double>> A) {
	for (int i = 0; i < A.size(); i++)
		showvec(A[i]);
}

vector<double> msolve(vector<double> c, vector<double> d, vector<double> e, vector<double> f) {
	int size = f.size();
	double len;
	double as = -e[0] / d[0], bs = f[0] / d[0];
	vector<double> a, b;
	a.push_back(as), b.push_back(bs);
	double atemp, btemp;
	for (int i = 1; i < size - 1; i++) {
		atemp = (-e[i] / (c[i] * a[i - 1] + d[i]));
		btemp = ((f[i] - c[i] * b[i - 1]) / (c[i] * a[i - 1] + d[i]));
		a.push_back(atemp), b.push_back(btemp);
	}
	len = a.size();
	vector<double> x;
	x.push_back((f[size - 1] - c[size - 1] * b[len - 1]) / (c[size - 1] * a[len - 1] + d[size - 1]));
	for (int i = 0; i < len; i++)
		x.push_back(x[i] * a[len - i - 1] + b[len - i - 1]);
	reverse(x.begin(), x.end());
	return x;
}

/*
vector<double> Deriv(vector <double> vec, vector<double> grid) {
	vector<double> rez;
	double dd;
	rez.push_back(vec[0]);
	for (int i = 1; i < vec.size() - 1; i++) {
		dd = grid[i] - grid[i - 1];
		rez.push_back((vec[i - 1] - 2 * vec[i] + vec[i + 1]) / pow(dd, 2));
	}
	//showvec(rez);
	rez.push_back(vec[vec.size() - 1]);
	return rez;
}
*/

vector<double> Deriv(vector <double> vec, vector<double> grid) {
	vector<double> rez;
	rez.push_back(vec[0]);
	for (int i = 1; i < vec.size() - 1; i++) {
		rez.push_back(2 * (((vec[i - 1] - vec[i]) / (grid[i - 1] - grid[i])) - ((vec[i] - vec[i + 1]) / (grid[i] - grid[i + 1]))) / (grid[i - 1] - grid[i + 1]));
	}
	//showvec(rez);
	rez.push_back(vec[vec.size() - 1]);
	return rez;
}

vector<double> summ(vector<double> f1, vector<double> f2) {
	vector<double> rez;
	for (int i = 0; i < f1.size(); i++)
		rez.push_back(f1[i] + f2[i]);
	return rez;
}

vector<double> diff(vector<double> f1, vector<double> f2) {
	vector<double> rez;
	for (int i = 0; i < f1.size(); i++)
		rez.push_back(f1[i] - f2[i]);
	return rez;
}

vector<double> mult(vector<double> f1, double k) {
	vector<double> rez;
	for (int i = 0; i < f1.size(); i++)
		rez.push_back(f1[i] * k);
	return rez;
}

vector<double> makesimplegrid(int size, double l, double r) {
	vector<double> grid;
	double delta = (r - l) / size;
	for (int i = 0; i <= size; i++) {
		grid.push_back(l + i * delta);
	}
	return grid;
}
vector<double> makeevengrid(int size, double l, double r) {
	double start = 0.0, step = 0.0, delta = (double) size * 0.25 / ((double)size*0.5+1.0)/12.5/size;
	cout << "d" << delta << endl;
	vector<double> grid;
	for (int i = 0; i < size / 2 + 1; i++) {
		start += step;
		step += delta;
		grid.push_back(start);
		if (start != 0.0)
			grid.insert(grid.begin(), - start);
	}
	return grid;
}
vector<double> makegrid(int size, double l, double r, double epsilon) {
	vector<double> grid;
	double delta = (r - l - 2 * epsilon) / (0.5 * (size)), smalldelta = (2 * epsilon) / (0.5 * size);
	//cout << delta;
	for (int i = 0; i <= size / 4; i++) {
		grid.push_back(l + i * delta);
	}
	for (int i = 1; i < size / 2; i++) {
		grid.push_back(-epsilon + i * smalldelta);
	}
	for (int i = size / 4; i >= 0; i--) {
		grid.push_back(r - i * delta);
	}
	cout << "delta: " << delta << " " << "smalldelta" << smalldelta << endl;
	return grid;
}

double l2norm(vector<double> vec) {
	double sum = 0;
	for (int i = 0; i < vec.size(); i++)
		sum += pow(vec[i], 2);
	return pow(sum, 0.5);
}





double scalar(vector<double> v1, vector<double> v2) {
	double sum = 0.0;
	for (int i = 0; i < v1.size(); i++)
		sum += v1[i] * v2[i];
	return sum;
}
vector<double> zeros(int n) {
	vector<double> rez = {};
	for (int i = 0; i < n + 1; i++)
		rez.push_back(0.0);
	return rez;
}
vector<double> Adotv(vector<vector<double>> A, vector<double> v) {
	double sum = 0;
	vector<double> x;
	for (int i = 0; i < A.size(); i++) {
		sum = 0;
		for (int j = 0; j < A[i].size(); j++) {
			sum += A[i][j] * v[j];
		}
		x.push_back(sum);
	}
	return x;
}

vector<double> cg(vector<vector<double>> A, vector<double> b, vector<double> x) {
	int len = b.size();
	vector<double> /*xprev,*/ rprev, zprev, xprev, r, z;
	double alpha, beta;
	for (int i = 0; i < len; i++) {
		xprev.push_back(0.0);
		//xprev.push_back(1.0);
	}
	rprev = diff(b, Adotv(A, xprev));
	zprev = rprev;
	int count = 0;
	while (Check(diff(xprev, x)) || count == 0) {
		cout << l2norm(diff(xprev, x)) << endl;
		count++;
		xprev = x;
		alpha = (scalar(rprev, rprev)) / scalar(Adotv(A, zprev), zprev);
		x = summ(xprev, mult(zprev, alpha));
		r = diff(rprev, mult(Adotv(A, zprev), alpha));
		beta = (scalar(r, r)) / scalar(rprev, rprev);
		z = summ(r, mult(zprev, beta));
		zprev = z;
		rprev = r;
	}
	cout << l2norm(diff(xprev, x)) << endl;
	return mult(x,0.5);
}

vector<vector<double>> makematr(vector<double> c, vector<double> d, vector<double> e) {
	int len = c.size();
	vector<vector<double>> A;
	vector<double> temp = {};
	temp.push_back(d[0]);
	temp.push_back(e[0]);
	for (int i = 0; i < len - 2; i++)
		temp.push_back(0.0);
	A.push_back(temp);
	for (int i = 1; i < len - 1; i++)
	{
		temp = {};
		for (int j = 0; j < i - 1; j++)
			temp.push_back(0.0);
		temp.push_back(c[i]);
		temp.push_back(d[i]);
		temp.push_back(e[i]);
		for (int j = i + 2; j < len; j++)
			temp.push_back(0.0);
		A.push_back(temp);
	}
	temp = {};
	for (int i = 0; i < len - 2; i++)
		temp.push_back(0.0);
	temp.push_back(c[len - 1]);
	temp.push_back(d[len - 1]);
	A.push_back(temp);
	return A;
}

vector<double> linpsi(vector<double> grid) {
	vector<double> rez;
	for (int i = 0; i < grid.size(); i++)
		rez.push_back(-grid[i]);
	return rez;
}
vector<double> mypsi(vector<double> grid, double gamma) {
	vector<double> rez;
	for (int i = 0; i < grid.size(); i++)
		rez.push_back(- tanh(gamma*grid[i])/tanh(gamma));
	return rez;

}

double min(double a, double b) {
	if (a > b)
		return b;
	else return a;
}

double dtest(double gamma) {
	return min( 3 / abs(1- gamma/tanh(gamma)), 0.1);
}

void dets(vector<double> c, vector<double> d, vector<double> e) {
	double f1 = 0, f2 = 1, f3 = d[0];
	for (int i = 1; i < c.size(); i++) {
		cout << "d" << f3 << endl;
		f1 = f2;
		f2 = f3;
		f3 = d[i] * f2 - c[i] * e[i] * f1;
		
	}
}


