#pragma once
#include <algorithm>

template <typename T>
class interval {
private:
	T _begin;
	T _end;

public:
	class divideByZero {};

	interval() : _begin(0), _end(0) {};
	interval(T t) : _begin(t), _end(t) {};
	interval(T begin, T end) : _begin(begin), _end(end) {};
	interval(const interval &i) : _begin(i._begin), _end(i._end) {};


	~interval() {};

	T begin() const { return _begin; }
	T end() const { return _end; }

	bool operator==(const interval &rhs) const;

	interval& operator+=(const interval &rhs);
	interval operator+(const interval &rhs) const;

	interval& operator-=(const interval &rhs);
	interval operator-(const interval &rhs) const;

	interval& operator*=(const interval &rhs);
	interval operator*(const interval &rhs) const;
	interval& operator*=(const T &t);
	interval operator*(const T &t) const;

	interval& operator/=(const interval &rhs);
	interval operator/(const interval &rhs) const;
	interval& operator/=(const T &t);
	interval operator/(const T &t) const;
};

template <typename T> inline
bool interval<T>::operator==(const interval &rhs) const {
	return (_begin == rhs._begin && _end == rhs._end);
}

template <typename T> inline
interval<T>& interval<T>::operator+=(const interval &rhs) {
	_begin += rhs._begin;
	_end += rhs._end;
	return *this;
}

template <typename T> inline
interval<T> interval<T>::operator+(const interval &rhs) const {
	return interval(_begin + rhs._begin, _end + rhs._end);
}

template <typename T> inline
interval<T>& interval<T>::operator-=(const interval &rhs) {
	_begin -= rhs._end;
	_end -= rhs._begin;
	return *this;
}

template <typename T> inline
interval<T> interval<T>::operator-(const interval &rhs) const {
	return interval(_begin - rhs._end, _end - rhs._begin);
}

template <typename T> inline
interval<T>& interval<T>::operator*=(const interval &rhs) {
	T arr[] = { _begin*rhs._begin, _end*rhs._end, _begin*rhs._end, _end*rhs._begin };
	_begin = *std::min_element(arr, arr + 4);
	_end = *std::max_element(arr, arr + 4);
	return *this;
}

template <typename T> inline
interval<T> interval<T>::operator*(const interval &rhs) const {
	T arr[] = { _begin*rhs._begin, _end*rhs._end, _begin*rhs._end, _end*rhs._begin };
	return interval(*std::min_element(arr, arr + 4), *std::max_element(arr, arr + 4));
}

template <typename T> inline
interval<T>& interval<T>::operator*=(const T &t) {
	_begin *= t;
	_end *= t;
	return *this;
}

template <typename T> inline
interval<T> interval<T>::operator*(const T &t) const {
	return interval(t * _begin, t * _end);
}


template <typename T> inline
interval<T>& interval<T>::operator/=(const interval &rhs) {
	if (rhs._begin == 0 || rhs._end == 0) throw divideByZero();
	T arr[] = { _begin / rhs._begin, _end / rhs._end, _begin / rhs._end, _end / rhs._begin };
	_begin = *std::min_element(arr, arr + 4);
	_end = *std::max_element(arr, arr + 4);
	return *this;
}

template <typename T> inline
interval<T> interval<T>::operator/(const interval &rhs) const {
	if (rhs._begin == 0 || rhs._end == 0) throw divideByZero();
	T arr[] = { _begin / rhs._begin, _end / rhs._end, _begin / rhs._end, _end / rhs._begin };
	return interval(*std::min_element(arr, arr + 4), *std::max_element(arr, arr + 4));
}

template <typename T> inline
interval<T>& interval<T>::operator/=(const T &t) {
	if (t == 0) throw divideByZero();
	_begin /= t;
	_end /= t;
	return *this;
}

template <typename T> inline
interval<T> interval<T>::operator/(const T &t) const {
	if (t == 0) throw divideByZero();
	return interval(_begin / t, _end / t);
}

