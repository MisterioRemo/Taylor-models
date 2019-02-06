#pragma once
#include "series.h"
#include <functional>
#include <fstream>
#include <math.h>
using std::function;

template <typename T>
class equation {

	using mfunction = powerSeries<T>(equation<T>::*)(vector<powerSeries<T> > &);

private:
	multSerCoef *coef;
	vector<interval<T> > parameter;
	vector<powerSeries<T> > u;
	int sizeVar;		// сколько первых уравнений системы действительно считаем
	int sizeParam;		// сколько параметров представленно в виде ряда
	double h;
	const double EPS = 0.000001;
	std::ofstream fout;


	powerSeries<T> pFun1(vector<powerSeries<T> > &u);
	powerSeries<T> pFun2(vector<powerSeries<T> > &u);


	T startInterval(const T &begin, const T &end);
	void findFirstPositionInSerie(vector<int> *vec);	

	inline char nextState(char);
	void searchPoints(vector<char>, int, int);
	vector<T> statesToStatesT(vector<char>);
	void printPoints(vector<char>, int);
	void printPlot();


public:
	class notSimetricStartInterval {};

	vector<mfunction> pFun = { &equation<T>::pFun1, &equation<T>::pFun2 };

	equation(int nvar, int param, int order) {
		coef = new multSerCoef(nvar, param, order);		
		sizeVar = coef->realVariable();
		sizeParam = coef->realParameter();

		for (int i = 0; i < sizeVar + sizeParam; i++)
			u.push_back(powerSeries<T>(coef->serieSize(), coef));

	};

	~equation() {
		delete coef;
	};


	inline vector<powerSeries<T> > getODU() const;
	inline powerSeries<T> getODU(int i) const;

	void initialFlow(vector<interval<T> >*);
	void RungeKutta(double, double, double, bool = false, int = 0, std::string = "function.dat");
	void printPlot(std::string);
};

// для задания симметричного начального интервала на [-1; 1]
template <typename T>
T equation<T>::startInterval(const T &begin, const T &end) {
	T i = (end - begin) / 2.0;
	i = (i < 0) ? -i : i;

	if (end - i != begin + i)
		throw notSimetricStartInterval();

	return i;
}

// номера позиций членов в первой степени 
template <typename T>
void equation<T>::findFirstPositionInSerie(vector<int> *vec) {
	for (int i = 0; i < coef->serieSize(); i++) {
		if (coef->getMultOrder(i) == 1)
			vec->push_back(i);
	}
}

template <typename T>
void equation<T>::initialFlow(vector<interval<T> > *points) {
	int size = sizeVar + sizeParam;
	T p;
	vector<int> pos;
	findFirstPositionInSerie(&pos);

	if (points->size() < size)
		points->resize(size, interval<T>(0, 0));

	for (int i = 0; i < size; i++) {

		p = startInterval((*points)[i].begin(), (*points)[i].end());

		parameter.push_back(interval<T>(-p, p));
		u[i].serie(0, (*points)[i].begin() + p);
		u[i].serie(pos[i], 1);
	}
}

template <typename T> 
inline vector<powerSeries<T> > equation<T>::getODU() const {
	return u;
}

template <typename T> 
inline powerSeries<T> equation<T>::getODU(int i) const {
	if (u.size() > i) 
		return u[i];
}

template <typename T>
powerSeries<T> equation<T>::pFun1(vector<powerSeries<T> > &v) {
	return v[1];
}

template <typename T>
powerSeries<T> equation<T>::pFun2(vector<powerSeries<T> > &v) {
	return v[0] * v[0];
}

template <typename T>
void equation<T>::RungeKutta(double tStart, double tEnd, double h, bool plot = false, int plotStep = 0, std::string filename = "function.dat") {
	vector<powerSeries<T> > K1(sizeVar), K2(sizeVar), K3(sizeVar), K4(sizeVar), v(sizeVar);
	int j, i, k = 0,
		r = 1.0 / h / 2;
	if (plot) {
		fout.close();
		fout.open(filename);
		if (plotStep > 0) 
			r = plotStep;
	}
	

	while (tStart < tEnd + EPS) {
		if (plot && k % r == 0) {
			printPlot();
			k = 0;
		}
		k++;

		// runge-kutta
		for (i = 0; i < sizeVar; i++) //k1
			K1[i] = (this->*pFun[i])(u) * h;
		for (j = 0; j < sizeVar; j++) //v2
			v[j] = u[j] + K1[j] / 2.0;

		for (i = 0; i < sizeVar; i++) //k2
			K2[i] = (this->*pFun[i])(v) * h;
		for (j = 0; j < sizeVar; j++) //v3
			v[j] = u[j] + K2[j] / 2.0;

		for (i = 0; i < sizeVar; i++) //k3
			K3[i] = (this->*pFun[i])(v) * h;
		for (j = 0; j < sizeVar; j++) //v4
			v[j] = u[j] + K3[j];

		for (i = 0; i < sizeVar; i++) //k4
			K4[i] = (this->*pFun[i])(v) * h;

		for (i = 0; i < sizeVar; i++)
			u[i] = u[i] + (K1[i] + (K2[i] + K3[i]) * 2 + K4[i]) / 6;

		tStart += h;
	}

	if (fout) fout.close();
	return;
}

////////////////////////////////////////////////
//	print plot
////////////////////////////////////////////////
template <typename T> 
inline char equation<T>::nextState(char c) {
	return (c < 1) ? c + 1 : c;
}

template <typename T>
void equation<T>::printPlot(std::string filename) {
	fout.open(filename);
	printPlot();
	fout.close();
}

template <typename T>
void equation<T>::printPlot() {
	int pSize = parameter.size();

	for (int i = 0; i < pSize; i++) {
		vector<char> states(pSize, 0);
		printPoints(states, i);
		searchPoints(states, i, 0);
	}
	fout << "\n";
}

template <typename T>
void equation<T>::searchPoints(vector<char> states, int cur, int pos) {
	if (pos >= states.size())
		return;

	searchPoints(states, cur, pos + 1);
	states[pos] = nextState(states[pos]);

	searchPoints(states, cur, pos + 1);
	// for (char c : states)			std::cout << +c << " ";		std::cout << endl;
	printPoints(states, cur);
}

template <typename T>
vector<T> equation<T>::statesToStatesT(vector<char> states) {
	vector<T> statesT;

	for (int i = 0; i < states.size(); i++)
		statesT.push_back((states[i] == 0) ? parameter[i].begin() : parameter[i].end());

	return statesT;
}

template <typename T>
void equation<T>::printPoints(vector<char> states, int cur) {
	vector<T> statesT = statesToStatesT(states);
	T sum, p,
		h = (parameter[cur].end() - parameter[cur].begin()) / 30.0;
	int i = 0;

	if (h == 0) return;

	for (; statesT[cur] < parameter[cur].end() + EPS; statesT[cur] += h) {
		for (i = 0; i < sizeVar; i++) {
			sum = 0;

			for (int j = 0; j < coef->serieSize(); j++) {	// сумма ряда для u[i]
				p = 1;
				for (int k = 0; k < coef->realVariable() + coef->realParameter(); k++) {
					p *= pow(statesT[k], coef->orderTable[j][k]);
				}
				sum += u[i].serie(j) * p;
			}

			fout << sum << " ";
		}
		fout << "\n";
	}
	if (i != 0) fout << "\n";	// для нормального постороения в gnuplot
	return;
}