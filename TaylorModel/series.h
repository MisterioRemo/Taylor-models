/*
Ряды и операции над ними.
Реализованные арифметические операции описаны в
Taylor models and floating-point arithmetic: proof that arithmetic operations are validated in COSY
N. Revol, K. Makino, M.Berz
*/

#pragma once
#include "interval.h"
#include "coefficients.h"
#include <vector>
using std::vector;

const double Em = 1e-15;
const double Ec = 1e-20;
const double E = 2;


template <typename T>
inline T mabs(T t) { return (t > 0) ? t : -t; }

template <typename T>
class powerSeries {
private:
	vector<T> _series;
	interval<T> _error;

public:
	static multSerCoef *_coef;
	class notTheSameLength {};
	class outOfRange {};
	class divideByZero {};

	powerSeries() : _error(interval<T>(0)) {};
	powerSeries(int size, multSerCoef *coef) {
		_series.resize(size);
		_coef = coef;

		for (int i = 0; i < size; i++) _series[i] = 0;
		_error = interval<T>(0);
	}

	~powerSeries() {};

	inline vector<T> serie() const { return _series; }
	inline T serie(int index) const { return _series[index]; }
	inline void serie(int index, T t) { _series[index] = t; }

	inline interval<T> error() const { return _error; }
	inline void error(T begin, T end) { _error._begin = begin; _error._end = end; }


	powerSeries& operator=(const powerSeries &ps);


	inline T operator[](int index) const { return _series[index]; }
	inline T& operator[](int index) {
		if (index < 0 || index > _series.size())
			throw outOfRange();
		return _series[index];
	}


	powerSeries& operator+=(const powerSeries&);
	powerSeries operator+(const powerSeries&) const;

	powerSeries& operator-=(const powerSeries&);
	powerSeries operator-(const powerSeries&) const;

	powerSeries operator*(const T&) const;
	powerSeries operator*(const powerSeries &ps) const;


	powerSeries operator/(const T &a) const;

};


template <typename T> multSerCoef *powerSeries<T>::_coef;

template <typename T>
powerSeries<T>& powerSeries<T>::operator=(const powerSeries &ps) {
	if (this != &ps) {
		_error = ps._error;
		_series.resize(ps._series.size());
		_series = ps._series;
	}
	return *this;
}

template <typename T>
powerSeries<T>& powerSeries<T>::operator+=(const powerSeries &ps) {
	if (_series.size() != ps._series.size())
		throw notTheSameLength();

	T t = 0;
	T s = 0;
	for (int i = 0; i < _series.size(); i++) {
		t += (mabs(_series[i]) > mabs(ps._series[i])) ? mabs(_series[i]) : mabs(ps._series[i]);
		_series[i] += ps._series[i];

		if (mabs(_series[i]) < Ec) {
			s += mabs(_series[i]);
			_series[i] = 0;
		}
	}
	_error += ps._error + interval<T>(-t, t)*Em*E + interval<T>(-s, s)*E;
	return *this;
}

template <typename T>
powerSeries<T> powerSeries<T>::operator+(const powerSeries &ps) const {
	if (_series.size() != ps._series.size())
		throw notTheSameLength();

	T t = 0;
	T s = 0;
	powerSeries sum;
	for (int i = 0; i < _series.size(); i++) {
		t += (mabs(_series[i]) > mabs(ps._series[i])) ? mabs(_series[i]) : mabs(ps._series[i]);
		sum._series.push_back(_series[i] + ps._series[i]);

		if (mabs(sum._series[i]) < Ec) {
			s += mabs(_series[i]);
			sum._series[i] = 0;
		}
	}
	sum._error = _error + ps._error + interval<T>(-t, t)*Em*E + interval<T>(-s, s)*E;
	return sum;
}

template <typename T>
powerSeries<T>& powerSeries<T>::operator-=(const powerSeries &ps) {
	if (_series.size() != ps._series.size())
		throw notTheSameLength();

	T t = 0;
	T s = 0;
	for (int i = 0; i < _series.size(); i++) {
		t += (mabs(_series[i]) > mabs(ps._series[i])) ? mabs(_series[i]) : mabs(ps._series[i]);
		_series[i] -= ps._series[i];

		if (mabs(_series[i]) < Ec) {
			s += mabs(_series[i]);
			_series[i] = 0;
		}
	}

	_error -= ps._error + interval<T>(-t, t)*Em*E + interval<T>(-s, s)*E;
	return *this;
}

template <typename T>
powerSeries<T> powerSeries<T>::operator-(const powerSeries &ps) const {
	if (_series.size() != ps._series.size())
		throw notTheSameLength();

	T t = 0;
	T s = 0;
	powerSeries sub;
	for (int i = 0; i < _series.size(); i++) {
		t += (mabs(_series[i]) > mabs(ps._series[i])) ? mabs(_series[i]) : mabs(ps._series[i]);
		sub._series.push_back(_series[i] - ps._series[i]);

		if (mabs(sub._series[i]) < Ec) {
			s += mabs(sub._series[i]);
			sub._series[i] = 0;
		}
	}

	sub._error = _error - ps._error + interval<T>(-t, t)*Em*E + interval<T>(-s, s)*E;
	return sub;
}

template <typename T>
powerSeries<T> powerSeries<T>::operator*(const T &a) const {
	T t = 0;
	T s = 0;

	powerSeries ps;
	for (int i = 0; i < _series.size(); i++) {
		ps._series.push_back(_series[i] * a);
		t += ps._series[i];

		if (mabs(ps._series[i]) < Ec) {
			s += mabs(ps._series[i]);
			ps._series[i] = 0;
		}
	}
	ps._error = _error * a + interval<T>(-t, t)*Em*E + interval<T>(-s, s)*E;

	return ps;
}

template <typename T>
powerSeries<T> powerSeries<T>::operator*(const powerSeries &ps) const {
	powerSeries mul(_series.size(), _coef);
	T p = 0;
	T t = 0;
	int index;

	for (int i = 0; i < _series.size(); i++) {
		interval<T> J = interval<T>(0, 0);

		for (int j = 0; j < _series.size(); j++) {

			if ((index = _coef->getMultIndex(i, j)) != -1) {
				p = _series[i] * ps._series[j];
				t += mabs(p);
				t += (mabs(mul._series[index]) > mabs(p))	? mabs(mul._series[index]) : mabs(p);
				mul._series[index] += p;
			}
			else {
				J += interval<T>(-mabs(ps._series[j]), mabs(ps._series[j]));
			}
		}
		mul._error += interval<T>(-mabs(_series[i]), mabs(_series[i])) * (J + ps._error);
	}

	interval<T> temp(0, 0);
	for (int j = 0; j < _series.size(); j++) {
		temp += interval<T>(-mabs(ps._series[j]), mabs(ps._series[j]));
	}
	mul._error += _error * (ps._error + temp);

	T s = 0;
	for (int k = 0; k < _series.size(); k++) {
		if (mabs(mul._series[k]) < Ec) {
			s += mabs(mul._series[k]);
			mul._series[k] = 0;
		}
	}
	mul._error += interval<T>(-t, t)*E*Em + interval<T>(-s, s)*E;

	return mul;
}

template <typename T>
powerSeries<T> powerSeries<T>::operator/(const T &a) const {
	if (a == 0)
		throw divideByZero();

	T t = 0;
	T s = 0;
	powerSeries ps;

	for (int i = 0; i < _series.size(); i++) {
		ps._series.push_back(_series[i] / a);
		t += ps._series[i];

		if (mabs(ps._series[i]) < Ec) {
			s += mabs(ps._series[i]);
			ps._series[i] = 0;
		}
	}
	ps._error = _error / a + interval<T>(-t, t)*Em*E + interval<T>(-s, s)*E;

	return ps;
}



