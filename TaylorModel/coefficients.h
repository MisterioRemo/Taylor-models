#pragma once
#include <vector>
#include <math.h>
#include <algorithm>
#include <iostream>
using std::vector;

// multiplication Series Coefficients
class multSerCoef {
	template <typename T> friend class equation;	// чтобы получить доступ к orderTable

private:
	vector<int> _sumOrder;	// sumOrder[i] = sum( orderTable[i][0..j] )
	vector<vector<int> > orderTable;
	vector< vector<int> > D;
	vector< vector<int> > C;
	int _order;         // порядок
	int _variable;      // кол-во параметров системы (т.е. сколько переменных учавствует в сосотаве ряда)
						// всегда приводится к чётному значению из-за алгоритма перемножения
	int _realParameter;     // кол-во параметров системы (alpha, beta...)
	int _realVariable;	// кол-во переменных (искомых), т.е. x, y, z...
	int _seriesSize;

	struct sortTableCoef {
		int _c1;
		int _c2;
		int _sumOrder;
		vector<int> _order;

		sortTableCoef() : _c1(0), _c2(0), _sumOrder(0) {};
		sortTableCoef(int c1, int c2, int sum) : _c1(c1), _c2(c2), _sumOrder(sum) {};
		inline void setToZero() { _c1 = _c2 = _sumOrder = 0; }	// правильно бы было добавить _order.clear(), 
																// но в данной реализации это не нужно
	};

	void findC();
	void findPowArray(vector<int>&);
	double findSeriesSize(double, double, double);
	static bool sortCFirstStep(const sortTableCoef&, const sortTableCoef&);
	void sortCSecondStep(vector<sortTableCoef>&);
	void removeByValue(vector<sortTableCoef>&, int);

	void findD();
	int findDElementC1(int);
	int findDElementC2(int);



public:
	multSerCoef() {};
	multSerCoef(int, int, int);
	~multSerCoef() {};

	inline int order() const { return _order; }
	inline int realVariable() const { return _realVariable; }
	inline int variableEven() const { return _variable; }
	inline int realParameter() const { return _realParameter; }
	inline int serieSize() const { return _seriesSize; }

	int getMultIndex(int, int) const;
	int getMultOrder(int) const;

	void printTableC() const;
	void printTableD() const;
};


