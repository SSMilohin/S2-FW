# Математическая библиотека 
Итоговая работа по основам программирования за второй семестр  

Класс Matrix является шаблонным, где T - тип элементов матрицы, N - количество строк матрицы, M - количество столбцов матрицы
## Конструкторы класса Matrix:
Matrix() - матрица состоящая из нулей  
Matrix(int diag) - матрица, состоящая из нулей, кроме элементов главной диагонали, которые равны значению diag  
Matrix(const T mas[N][M]) - матрица созданная с помощью массива  
Matrix(const MasWrapper<T, N, M>& mas) - матрица созданная с помощью структуры MasWrapper  
Matrix(const Matrix<T, N, M>& mat) - матрица созданная с помощью другой матрицы

## Для класса Matrix перегружены следующие операторы:  
"-","+","=","+=","*", а также операторы ввода и вывода

## Методы класса Matrix:
.getN() - возвращает количество строк матрицы (int)  
.getM() - возвращает количество строк матрицы (int)  
.get(int i, int j) - возвращает элемент матрицы из i-ой строки и j-го столбца (T)  
.set(int i, int j, T data) - устанавливает значение элемента i-ой строки и j-го столбца равное "data"  
.transpose() - Возвращает транспонированную матрицу размерами MxN  
.attach(const Matrix<T, N, L>& mat) - Возвращает матрицу размерами Nx(M+L), к исходной матрице (NxM) "дописывается" матрица mat (NxL) для дальнейшего решения СЛАУ, либо нахождения обратной матрицы  
.toStepMatrix() - Возвращает матрицу такого же размера (NxM), приведённую к ступенчатому виду по строкам  
.toCanonicalStepMatrix() - Возвращает матрицу такого же размера (NxM), приведённую к каноническому виду по строкам  
.determinant() - Возвращает значение определителя матрицы (double)  
.inversion() - Возвращает обратную матрицу такого же размера (NxM)  
.SLAU(const Matrix<T, N, 1>& mat) - Возвращает матрицу размерами (Nx1)[Столбец], которая является решением СЛАУ, к которой "приписали" столбец свободных членов mat

## Для разработки  
Необходим Cmake 3.15, компилятор VS 2019  