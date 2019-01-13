#include "coefficients.h"

multSerCoef::multSerCoef(int nvar, int param, int order) {
	_order = order;
	_realParameter = param;
	_realVariable = nvar;
	nvar += param;
	_variable = (nvar % 2) ? nvar + 1 : nvar;
	_seriesSize = findSeriesSize(_order + _variable, _variable, _order);

	C.resize(2);
	D.resize(2);

	findC();
	findD();
}

double multSerCoef::findSeriesSize(double np, double n, double p) {
	if (np < 1) return 1;
	n = (n < 1) ? 1 : n;
	p = (p < 1) ? 1 : p;
	return (np / (n*p)) * findSeriesSize(np - 1, n - 1, p - 1);
}

/*
Реализация алгоритма, рассмотренного в статье
Algorithms for higher order automatic differentiation in many variables with applications to beam physics
Martin Berz
*/
void multSerCoef::findC() {
	vector<int> indices(_variable, 0);
	vector<int> alphabet(_order + 1);
	for (int i = 0; i <= _order; i++)
		alphabet[i] = i;
	const int stopIndex = alphabet.size() - 1;


	int count = 0;  // счётчик найденных элементов ряда
	vector<int> powArray(_variable / 2); // 1 order+1 (order+1)^2 ... (order+1)^(variable/2)
	findPowArray(powArray);

	vector<sortTableCoef> sortVec;
	sortTableCoef elem(0, 0, 0);
	elem._order.resize(_variable);

	bool flag = true;
	while (flag) {
		if (count >= _seriesSize) break;

		elem.setToZero();
		for (int i = 0; i < _variable / 2; i++) {
			elem._c1 += alphabet[indices[i]] * powArray[i]; // вычисление коэф. по формуле
			elem._c2 += alphabet[indices[i + _variable / 2]] * powArray[i];
			elem._sumOrder += alphabet[indices[i]] + alphabet[indices[i + _variable / 2]];

			elem._order[i] = alphabet[indices[i]];
			elem._order[i + _variable / 2] = alphabet[indices[i + _variable / 2]];
		}

		if (elem._sumOrder <= _order) {
			sortVec.push_back(elem);
			count++;
		}

		// перебор строк длины _variable для алфавита = {0.._order}
		int pos = _variable;
		while (indices[--pos] == stopIndex) {   // если на данной позиции перебраны все символы из алфавита
			if (pos == 0) {
				flag = false;                   // заканчиваем работу, т.к. перебраны все варианты
				break;
			}
			indices[pos] = 0;
		}
		++indices[pos];                         // увеличиваем индекс левого соседа

	}


	std::stable_sort(sortVec.begin(), sortVec.end(), sortCFirstStep);
	sortCSecondStep(sortVec);	
}

void multSerCoef::findPowArray(vector<int> &arr) {
	arr[0] = 1;
	for (int i = 1; i < arr.size(); i++) {
		arr[i] = arr[i - 1] * (_order + 1);
	}
}

// устойчивая(!) сортировка в первую очередь по c2,
// затем по степени элемента ряда,
// затем в зависимости от c1
bool multSerCoef::sortCFirstStep(const sortTableCoef &a, const sortTableCoef &b) {
	return (a._c2 == b._c2) ?
		((a._sumOrder == b._sumOrder) ? a._c1 < b._c1 : a._sumOrder < b._sumOrder)
		: ((a._sumOrder == b._sumOrder) ? a._c2 < b._c2 : a._sumOrder < b._sumOrder);
}

void multSerCoef::sortCSecondStep(vector<sortTableCoef> &v) {
	int k = 0;
	int j;

	while (!v.empty()) {
		k = v[0]._c2;

		for (j = 0; j < v.size(); j++) {
			if (v[j]._c2 == k) {
				C[0].push_back(v[j]._c1);
				C[1].push_back(v[j]._c2);
				_sumOrder.push_back(v[j]._sumOrder);
				orderTable.push_back(v[j]._order);
			}
		}
		removeByValue(v, k);    // удалить все элементы, где _c2 = k
	}
}

void multSerCoef::removeByValue(vector<sortTableCoef> &v, int val) {
	v.erase(std::remove_if(v.begin(), v.end(), [&](sortTableCoef &stc){ return stc._c2 == val; }), v.end());
}

void multSerCoef::findD() {
	int size = pow(_order + 1, _variable / 2);
	for (int i = 0; i < size; i++) {
		D[0].push_back(findDElementC1(i));
		D[1].push_back(findDElementC2(i));
	}
}

int multSerCoef::findDElementC1(int elem) {
	for (int i = 0; i < _seriesSize; i++) {
		if (C[0][i] == elem)
			return i + 1;
	}
	return 0;
}

int multSerCoef::findDElementC2(int elem) {
	for (int i = 0; i < _seriesSize; i++) {
		if (C[1][i] == elem)
			return i;
	}
	return 0;
}

int multSerCoef::getMultIndex(int index1, int index2) const {
	if (getMultOrder(index1) + getMultOrder(index2) > _order)
		return -1;

	int c1 = C[0][index1] + C[0][index2];
	int c2 = C[1][index1] + C[1][index2];
	return D[0][c1] + D[1][c2] - 1;
}

int multSerCoef::getMultOrder(int index) const {
	return _sumOrder[index];
}

void multSerCoef::printTableC() const {
	std::cout << "I\t  C1\t|  C2\t|  sumOrder\t|  order" << std::endl;
	for (int i = 0; i < C[0].size(); i++) {
		std::cout << i + 1 << "\t  " << C[0][i] << "\t|  " << C[1][i] << "\t|  " << _sumOrder[i] << "\t|  ";
		for (int e : orderTable[i])
			std::cout << e << "  ";
		std::cout << std::endl;
	}
}

void multSerCoef::printTableD() const {
	std::cout << "I\t  D1\t|  D2" << std::endl;
	for (int i = 0; i < D[0].size(); i++) {
		std::cout << i  << "\t  " << D[0][i] << "\t|  " << D[1][i] << std::endl;
	}
}
