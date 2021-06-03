/*
 * Polynomial.h
 *
 *  Created on: 2021Äê4ÔÂ30ÈÕ
 *      Author: Administrator
 */

#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <algorithm>
#include <iostream>
#include <vector>
#include <utility>
#include <math.h>

using namespace std;
template <typename T>
class Polynomial {
public:
    vector<T> data;

    void normalize(vector<T> &coef) {
        for (auto iter = coef.rbegin(); iter != coef.rend();) {
            if (*iter == T(0)) {
				iter = typename vector<T>::reverse_iterator(coef.erase((++iter).base()));
            } else {
                break;
            }
        }
    }

//public:
    Polynomial<T>(const vector<T> &a):data(a) {
        normalize(data);
    }

    template <typename Iter>
    Polynomial<T>(Iter first, Iter last) :data(first, last) {
        normalize(data);
    }

    Polynomial<T>(const T &num = T()) {
        data.push_back(num);
        normalize(data);
    }

    T operator [] (size_t i) const {
        if (i >= data.size()) {
            return T(0);
        } else {
            return data[i];
        }
    }

    long long Degree() const {
        if (data.empty()) {
            return -1;
        } else {
            return static_cast<int>(data.size()) - 1;
        }
    }

    bool operator == (const Polynomial<T> &other) const {
        if (data.size() != other.data.size()) {
            return false;
        }
        for (size_t i = 0; i != data.size(); ++i) {
            if (data[i] != other[i]) {
                return false;
            }
        }
        return true;
    }

    bool operator != (const Polynomial<T> &other) const {
        return !(*this == other);
    }

    bool operator ==(const T &num) {
        return *this == Polynomial<T>(num);
    }

    bool operator !=(const T &num) {
        return *this != Polynomial<T>(num);
    }

    Polynomial<T>& operator += (const Polynomial<T> &other) {
        data.resize(max(other.data.size(), data.size()), T(0));
        for (size_t i = 0; i != min(data.size(), other.data.size()); ++i) {
            data[i] += other.data[i];
        }
        normalize(data);
        return *this;
    }

    Polynomial<T>& operator -= (const Polynomial<T> &other) {
        data.resize(max(other.data.size(), data.size()), T(0));
        for (size_t i = 0; i != min(data.size(), other.data.size()); ++i) {
            data[i] -= other.data[i];
        }
        normalize(data);
        return *this;
    }

    Polynomial<T>& operator +=(const T &num) {
        *this += Polynomial<T>(num);
        normalize(data);
        return *this;
    }

    Polynomial<T>& operator -=(const T &num) {
        *this -= Polynomial<T>(num);
        normalize(data);
        return *this;
    }

    Polynomial<T>& operator *=(const Polynomial<T> &other) {
        vector<T> temp(data.size() + other.data.size(), T(0));
        for (size_t i = 0; i != data.size(); ++i) {
            for (size_t j = 0; j != other.data.size(); ++j) {
                temp[i + j] += data[i] * other.data[j];
            }
        }
        normalize(temp);
        *this = Polynomial(temp);
        return *this;
    }

    Polynomial<T>& operator *=(const T &num) {
        for (size_t i = 0; i != data.size(); ++i) {
            data[i] *= num;
        }
        normalize(data);
        return *this;
    }

    T operator () (const T &point) const {
        T ans = T(0);
        for (auto iter = data.rbegin(); iter != data.rend(); ++iter) {
            ans += *iter;
            if ((iter + 1) != data.rend()) {
                ans *= point;
            }
        }
        return ans;
    }


    friend ostream& operator << (ostream& out, const Polynomial<T> &pol) {
        bool flag = false;
        unsigned long long degree = pol.data.size() - 1;
        for (auto iter = pol.data.rbegin(); iter != pol.data.rend(); ++iter, --degree) {
            T coef = *iter;
            if (coef != T(0)) {
                if (coef > T(0) && flag) {
                    out << '+';
                }
                flag = true;
                if (degree == 0) {
                    out << coef;
                } else if (coef == T(1)) {
                    out << 'x';
                } else if (coef == T(-1)) {
                    out << "-x";
                } else {
                    out << coef << "*x";
                }
                if (degree > 1) {
                    out << '^' << degree;
                }
            }
        }
        if (pol.data.size() == 0) {
            out << 0;
        }
        return out;
    }

    friend Polynomial<T> operator&(const Polynomial<T> &first, const Polynomial<T> &second) {
        Polynomial<T> comp(first.data.at(0));
        Polynomial<T> copy(second.data);
        size_t iter = 1;
        for (size_t degree = 1; degree != first.data.size(); ++degree) {
            for (; iter != degree; ++iter) {
                copy *= second;
            }
            comp += copy * first.data[degree];
        }
        return comp;
    }

    Polynomial<T> &operator /= (const Polynomial<T> &other) {
        Polynomial<T> priv(T(0));
        while (data.size() >= other.data.size()) {
            T coef = data.back() / other.data.back();
            size_t degree = data.size() - other.data.size();
            vector<T> div(degree + 1);
            div.back() = coef;
            Polynomial<T> temp(div);
            *this -= other * temp;
            priv += temp;
        }
        data = priv.data;
        return *this;
    }

    Polynomial<T> &operator %= (const Polynomial<T> &other) {
        Polynomial<T> quotient = *this / other;
        *this -= other * quotient;
        return *this;
    }

    friend Polynomial<T> operator,(const Polynomial<T> &first, const Polynomial<T> &second) {
        Polynomial<T> gcd = first;
        Polynomial<T> copy = second;
        while (copy.data.size() != 0) {
            gcd %= copy;
            swap(gcd, copy);
        }
        if (gcd.data.size() != 0) {
            Polynomial<T> temp(gcd[gcd.data.size() - 1]);
            gcd /= temp;
        }
        return gcd;
    }

	/*
    auto begin() const {
        return data.begin();
    }

    auto end() const {
        return data.end();
    }
	*/
};


//will change first, be careful!!!
template <typename T>
Polynomial<T> operator *(Polynomial<T> first, const Polynomial<T> &second) {
    return first *= second;
}

template <typename T>
Polynomial<T> operator +(Polynomial<T> first, const Polynomial<T> &second) {
    return first += second;
}

template <typename T>
Polynomial<T> operator -(Polynomial<T> first, const Polynomial<T> &second) {
    return first -= second;
}

template<typename T>
Polynomial<T> operator /(const Polynomial<T> &first, const Polynomial<T> &second) {
    auto copy = first;
    copy /= second;
    return copy;
}

template<typename T>
Polynomial<T> operator %(const Polynomial<T> &first, const Polynomial<T> &second) {
    auto copy = first;
    copy %= second;
    return copy;
}

template <typename T>
Polynomial<T> operator +(Polynomial<T> poly, const T &num) {
    return poly += Polynomial<T>(num);
}

template <typename T>
Polynomial<T> operator +(const T &num, Polynomial<T> poly) {
    return poly += Polynomial<T>(num);
}

template <typename T>
Polynomial<T> operator -(Polynomial<T> poly, const T &num) {
    return poly -= Polynomial<T>(num);
}

template <typename T>
Polynomial<T> operator -(const T &num, Polynomial<T> poly) {
    return Polynomial<T>(num) -= poly;
}

template <typename T>
Polynomial<T> operator *(Polynomial<T> poly, const T &num) {
    return  poly *= Polynomial<T>(num);
}

template <typename T>
Polynomial<T> operator *(const T &num, Polynomial<T> poly) {
    return  poly *= Polynomial<T>(num);
}

template <typename T>
//assume constant of  integral as 0
Polynomial<T> integral(const Polynomial<T>& other, const T &C = 0) {
	Polynomial<T> res;
	vector<T>& data = res.data;
	data.resize(other.data.size()+1, 0);
    for (size_t i = 1; i < data.size(); ++i) {
        data[i] = other.data[i-1]/(i);
    }
    res.normalize(data);
    return res;
}

//here i mean current depth
template <typename T>
Polynomial<T> multiIntegral_intern_recursive(Polynomial<T> **vecPoly, int i, const T &left) {
	if(vecPoly[i] == NULL) return Polynomial<T>(T(0)); // if  vecPoly[i] == NULL
	if (i == 0) {
		Polynomial<T> res = integral(*(vecPoly[0]));
		res = res - res(left);
		return res;
	}

	Polynomial<T> tmp = multiIntegral_intern_recursive(vecPoly, i-1, left) * (*(vecPoly[i]));
	Polynomial<T> res = integral(tmp);

    return (res - res(left));
}

//here i is len(vecPoly)-1, to keep same as recursive version
template <typename T>
Polynomial<T> multiIntegral_intern(Polynomial<T> **vecPoly, int i, const T &left) {
	Polynomial<T> res = Polynomial<T>(T(1));
	for (int j = 0; j <= i; j++) {
		if (vecPoly[j] == NULL) return Polynomial<T>(T(0));
		Polynomial<T> tmp = res * (*(vecPoly[j]));
		res = integral(tmp);
		res = res - res(left);
	}
	return res;
}

template <typename T>
//give vecPoly as the distribution of each variable ordered by x(0),x(1),x(n)..., and left as the lower bound, variable as upper bound, like |p(x(n-1))|{left, x(n)}
//we return the polynomial of last x(n)
T multiIntegral(Polynomial<T>** vecPoly, const size_t len, const T &left, const T &right) {
	Polynomial<T> tmp = multiIntegral_intern_recursive(vecPoly, len-1, left);
	return tmp(right);
}

#endif
