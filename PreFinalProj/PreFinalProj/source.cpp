#include <iostream>
#include <math.h>
using namespace std;
#include <limits>
#include <vector>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include "funcs.h"
#define PI 3.14

typedef std::numeric_limits< double > dbl;


class LinAlg {
public:
	//Конструкторы
	LinAlg(vector<double> gr, double v, double e, double g, vector<double> ps) :grid(gr), V(v), eps(e), gamma(g), psi(ps) {
		d2psi = Deriv(ps, gr);
		linsolve();
	}
	//выводы
	vector<double> getgrid() {
		return grid;
	}

	vector<double> getfvec() {
		vector<double> fv;
		for (int i = 0; i < grid.size(); i++) {
			fv.push_back(f(grid[i], psi[i], d2psi[i]));
		}
		return fv;
	}

	vector<double> getgvec() {
		vector<double> gv;
		for (int i = 0; i < grid.size(); i++) {
			gv.push_back(g(grid[i], psi[i], d2psi[i]));
		}
		return gv;
	}

	vector<double> getrezvec() {

		return rezvec;
	}
private:
	double V, eps, gamma;
	vector<double> grid;
	vector<double> psi, d2psi;
	vector<double>  rezvec;
	

	vector<double> OneStep(vector<double> psi, vector<double> dpsi) {
		vector<double> itvec, F, newvec;
		vector<double> c, d, e;
		double s = psi.size();
		for (int i = 1; i < s - 1; i++) {
			F.push_back(FF(grid[i], psi[i], dpsi[i]));
		}
		//C
		c.push_back(1.0);
		for (int i = 2; i < s - 1; i++)
			c.push_back(2.0 / (grid[i - 1] - grid[i + 1]) / (grid[i - 1] - grid[i]));

		//D
		for (int i = 1; i < s - 1; i++) {
			d.push_back(DD(grid[i - 1], grid[i], grid[i + 1], psi[i], dpsi[i]));
		}
		//E`
		for (int i = 1; i < s - 2; i++)
			e.push_back(2.0 / (grid[i - 1] - grid[i + 1]) / (grid[i] - grid[i + 1]));
		e.push_back(1.0);
		//dets(c, d, e);
		//showmatr(makematr(c, d, e));
		//showvec(F);
		//showvec(grid);
		vector<double> rez = msolve(c, d, e, F);
		//vector<double> rez = cg(makematr(c, d, e), F, r);
		rez.insert(rez.begin(), 0.0);
		rez.push_back(0.0);
		return rez;
	}

	void linsolve() {
		rezvec = OneStep(psi, d2psi);
	}
	/*
	double q(double x) {
		if (x < 0)
			return -1;
		else if (x == 0)
			return 0;
		else
			return 1;
	}*/
	
	
	double q(double x) {
		if (x < 0)
			return -(-x * x - x + 1);
		else if (x == 0)
			return 0;
		else
			return -(x*x-x-1);
	}
	
	/*
	double q(double x) {
		if (x < 0)
			return -(0.1 * sin(2 * PI * x) + 1);
		else if (x == 0)
			return 0;
		else
			return -(-0.1 * sin(2 * PI * x)-1);
	}*/
	/*
	double q(double x) {
		return -V * pow(eps, 2) * 2.0 * pow(gamma, 2) * tanh(gamma * x) / tanh(gamma) / pow(cosh(gamma * x), 2) + sinh(V * tanh(gamma * x) / tanh(gamma)) / sinh(V);
	}*/
	/*
	double f(double x, double Psi, double D2psi) {
		return pow(eps, 2) * pow(gamma, 2) * (2.0) / tanh(gamma) / pow(cosh(gamma * x), 2) * tanh(gamma * x) - g(x, Psi, D2psi) * (x - tanh(gamma * x) / (tanh(gamma)));
	}*/
	
	
	double f(double x, double Psi, double D2psi) {
		return q(x) - G(x, Psi, D2psi) - pow(eps, 2)* D2psi;
	}
	
	double g(double x, double Psi, double D2psi) {
		return -cosh(Psi) / sinh(V);
	}

	double G(double x, double Psi, double D2psi) {
		return sinh(Psi) / sinh(V);
	}

	double FF(double x, double Psii, double D2psii) {
		return f(x, Psii, D2psii) / pow(eps, 2);

	}
	double DD(double x0, double x1, double x2, double Psii, double D2psii) {
		return -2.0 / (x0 - x1) / (x1 - x2) - g(x1, Psii, D2psii) / pow(eps, 2);
	}

};

class FullAlg {
public:
	FullAlg(vector<double> gr, double v, double e, double g) :grid(gr), V(v), eps(e), gamma(g) {
		rezvec = algsolve();
	}

	vector<double> getgrid() {
		return grid;
	}

	vector<double> getrezvec() {

		return rezvec;
	}


private:
	double V, eps, gamma;
	vector<double> grid;
	vector<double> psi0;
	vector<double>  rezvec;

	vector<double> algsolve() {
		vector<double> psik, Hk;
		int count = 0;
		for (int i = 0; i < grid.size(); i++) {
			Hk.push_back(1.0);
			psik.push_back(Psi0func(grid[i]));
		}

		do {
			count++;
			LinAlg ll(grid, V, eps, gamma, psik);
			Hk = ll.getrezvec();
			cout << Hknorm(Hk) << endl;
			psik = summ(psik, Hk);
		} while (Check(Hk));
		//cout << count << endl;
		return psik;

		//return psi0;
	}
	
	double Psi0func(double x) {
		return V * x;
	}
	
	/*
	double Psi0func(double x) {
		return - tanh(gamma*x)/tanh(gamma);
	}
	*/
};

void linplots(const char* resf, vector<double> grid, vector<double> vs, vector<double> epss, vector<double> gammas, vector<double> psi) {
	clear(resf);
	write(resf, vs);
	write(resf, epss);
	write(resf, gammas);
	write(resf, grid);
	for (int i = 0; i < vs.size(); i++) {
		for (int j = 0; j < epss.size(); j++) {
			for (int k = 0; k < gammas.size(); k++) {
				LinAlg lin(grid, vs[i], epss[j], gammas[k], psi);
				write(resf, lin.getrezvec());
				//showvec(lin.getrezvec());
			}
		}
	}
}

void linfgplots(const char* resf, vector<double> grid, vector<double> vs, vector<double> epss, vector<double> gammas, vector<double> psi) {
	clear(resf);
	write(resf, vs);
	write(resf, epss);
	write(resf, gammas);
	write(resf, grid);
	showvec(gammas);
	for (int i = 0; i < vs.size(); i++) {
		for (int j = 0; j < epss.size(); j++) {
			for (int k = 0; k < gammas.size(); k++) {
				LinAlg lin(grid, vs[i], epss[j], gammas[k], psi);
				write(resf, lin.getrezvec());
				write(resf, lin.getfvec());
				write(resf, lin.getgvec());
			}
		}
	}
}


void plots(const char* resf, vector<double> grid, vector<double> vs, vector<double> epss, vector<double> gammas) {
	clear(resf);
	write(resf, vs);
	write(resf, epss);
	write(resf, gammas);
	write(resf, grid);
	for (int i = 0; i < vs.size(); i++) {
		for (int j = 0; j < epss.size(); j++) {
			for (int k = 0; k < gammas.size(); k++) {
				FullAlg al(grid, vs[i], epss[j], gammas[k]);
				write(resf, al.getrezvec());
			}
		}
	}
}


int main() {
	vector<double> mygrid = makesimplegrid(2000, -1, 1), newgrid = makegrid(2000, -1, 1, 0.05);

	//LinAlg L(mygrid, -1.0, 0.01, 100.0, zeros(20000));
	//FullAlg Al(mygrid, -1.0, 0.01, 100.0);
	//showvec(L.getgrid());
	//showvec(Al.getrezvec());
	/*
	showvec(mygrid);
	showvec(newgrid);
	cout << newgrid.size();*/
	//linplots("test.txt", newgrid, { -1.0 }, { 0.001 }, {100.0 }, linpsi(newgrid));
	//cout << dtest(100) << setprecision(9);
	plots("test.txt", mygrid, { -5.0 }, { 0.001 }, { 100.0 });
	//showvec(linpsi(newgrid));

}
