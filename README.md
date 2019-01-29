# Численный метод решения задач Коши систем ОДУ с интервальными параметрами
---

Программа предназначена для решения систем ОДУ, содержащих некоторую неопределённость или неоднозначность в начальных условиях или параметрах (т.е. в тех случаях, когда начальные условия задаются интервалами). Реализация метода основана на символьных вычислениях и использовании модели Тейлора.<br>
Работа базируется главным образом на статье Markus Neher [«From Interval Analysis to Taylor Models - An Overview»](http://na.math.kit.edu/neher/preprnts/neher_2005_taylor_models_IMACS05.pdf). Математические операции над рядами описаны в статье N. Revol, K. Makino, M. Berz [«Taylor Models and Floating-Point Arithmetic: Proof that Arithmetic Operations are Validated in COSY»](https://bt.pa.msu.edu/cgi-bin/display.pl?name=TMJLAP03). Особое внимание уделено алгоритму перемножения рядов, которое описано в статье Martin Berz [«Algorithms for Higher Order Automatic Differentiation in Many Variables with Applications to Beam Physics»](https://bt.pa.msu.edu/cgi-bin/display.pl?name=adalgo).


Программа является ресурсозатратной, так как она, во-первых, оперирует рядами и, во-вторых, дополнительно вычисляет вспомогательные таблицы коэффициентов для быстрого перемножения рядов. Потому программа писалась с изначальным расчётом на то, что она будет решать одну-единственную систему ОДУ за раз (смиритесь =3 ).


---

### Примеры

#### Пример № 1
Разберём, как записать систему, описанную в статье «From Interval Analysis to Taylor Models - An Overview».

![\begin{array}{l} \\ u' = v \\ v' = u^2 \\ u(0) = 1 + a \\ v(0) = -1 + b \\ a, b \in [-0.05; 0.05] \\ t \in [0; 6] \end{array}](https://latex.codecogs.com/svg.latex?%5Cbegin%7Barray%7D%7Bl%7D%20%5C%5C%20u%27%20%3D%20v%20%5C%5C%20v%27%20%3D%20u%5E2%20%5C%5C%20u%280%29%20%3D%201%20&plus;%20a%20%5C%5C%20v%280%29%20%3D%20-1%20&plus;%20b%20%5C%5C%20a%2C%20b%20%5Cin%20%5B-0.05%3B%200.05%5D%20%5C%5C%20t%20%5Cin%20%5B0%3B%206%5D%20%5Cend%7Barray%7D)

Прежде всего следует понять, сколько здесь переменных, а сколько параметров. Подставив *a* и *b* в начальные условия, можно увидеть, что у системы нет параметров, есть только переменные. Записываем исходные интервалы в вектор (порядок важен!).
```cpp
vector<interval<double> > initPoint;
initPoint.push_back(interval<double>(0.95, 1.05));
initPoint.push_back(interval<double>(-1.05, -0.95));
```

Для инициализации системы требуется указать количество переменных и параметров системы и порядок ряда, в котором будут производиться расчёты. А также начальные условия системы, которые предаются в виде вектора интервалов.
```cpp
equation<double> odu(2, 0, 18);
odu.initialFlow(&initPoint);
```

Затем необходимо определить сами рассчитываемые функции. В классе *equation* объявляем все функции системы и помещаем их в вектор функций.
```cpp
powerSeries<T> pFun1(vector<powerSeries<T> >&);
powerSeries<T> pFun2(vector<powerSeries<T> >&);
vector<mfunction> pFun = { &equation<T>::pFun1, &equation<T>::pFun2 };
```

Определяем функции.
```cpp
template <typename T>
powerSeries<T> equation<T>::pFun1(vector<powerSeries<T> > &v) {
	return v[1];
}

template <typename T>
powerSeries<T> equation<T>::pFun2(vector<powerSeries<T> > &v) {
	return v[0] * v[0];
}
```

Наконец, осуществляем расчёт методом Рунге-Кутты 4-го порядка.
```cpp
odu.RungeKutta(0, 6, 0.01);
```

![function1](https://github.com/MisterioRemo/misterioremo.github.io/blob/master/taylor-model-img/fun1.gif?raw=true)


#### Пример № 2
Все последующие примеры взяты из статьи А.Ю. Морозова и Д.Л. Ревизникова «Алгоритм адаптивной интерполяции на основе kd-дерева для численного интегрирования систем ОДУ с интервальными начальными условиями»

![\begin{array}{l} \\ x' = -0.9x+0.5xy\\ y' = \alpha y - 0.8xy \\ x(0) \in [0.9; 1.1] \\ y(0) \in [1.9; 2.1] \\ \alpha \in [0.7; 0.75] \\ t \in [0; 100] \end{array}](https://latex.codecogs.com/gif.latex?%5Cbegin%7Barray%7D%7Bl%7D%20%5C%5C%20x%27%20%3D%20-0.9x&plus;0.5xy%20%5C%5C%20y%27%20%3D%20%5Calpha%20y%20-%200.8xy%20%5C%5C%20x%280%29%20%5Cin%20%5B0.9%3B%201.1%5D%20%5C%5C%20y%280%29%20%5Cin%20%5B1.9%3B%202.1%5D%20%5C%5C%20%5Calpha%20%5Cin%20%5B0.7%3B%200.75%5D%20%5C%5C%20t%20%5Cin%20%5B0%3B%20100%5D%20%5Cend%7Barray%7D)

Запишем начальные условия, сначала идут переменные и только затем параметры. Это важно!
```cpp
vector<interval<double> > initPoint;
initPoint.push_back(interval<double>(0.9, 1.1));
initPoint.push_back(interval<double>(1.9, 2.1));
initPoint.push_back(interval<double>(0.7, 0.75));

equation<double> odu(2, 1, 4);
odu.initialFlow(&initPoint);
```

Затем определяем функции.<br>
*Обратите внимание, что к переменным мы обращаемся через агрумент функции v[i], а вот к параметру системы через u[j]. Как определить j? У него тот же индекс, что и в векторе initPoint.*
```cpp
template <typename T>
powerSeries<T> equation<T>::pFun1(vector<powerSeries<T> > &v) {
	return v[0] * (-0.9) + v[0] * v[1] * 0.5;
}

template <typename T>
powerSeries<T> equation<T>::pFun2(vector<powerSeries<T> > &v) {
	return u[2] * v[1] + v[0] * v[1] * (-0.8);
}
```

Рассчитываем систему с выводом графика с шагом печати 30*h.
```cpp
equation<double> odu(2, 1, 4);
odu.initialFlow(&initPoint);

odu.RungeKutta(0, 100, 0.01, true, 30);
```

![function2](https://github.com/MisterioRemo/misterioremo.github.io/blob/master/taylor-model-img/fun2.gif?raw=true)

*На графике отражена проекция множества решений из трехмерного пространства в двухмерное.*