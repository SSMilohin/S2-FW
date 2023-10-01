#pragma once
#include <iostream>
#include <cmath>

namespace sm::math
{
	template<typename T, int N, int M>
	struct MasWrapper
	{
		T mas[N][M];
	};

	template<typename T, int N, int M>
	class Matrix
	{
		int m_n, m_m;
		T m_mat[N][M];
	public:
		Matrix() //Пустой конструктор матрицы
		{
			m_n = N;
			m_m = M;
			for (int i = 0; i < m_n; i++)
				for (int j = 0; j < m_m; j++)
					m_mat[i][j] = 0;
		}

		Matrix(T diag) //Конструктор матрицы с элементами diag на главной диагонали
		{
			m_n = N;
			m_m = M;
			for (int i = 0; i < m_n; i++)
				for (int j = 0; j < m_m; j++)
				{
					if (i != j) m_mat[i][j] = 0;
					if (i == j) m_mat[i][j] = diag;
				}
		}

		Matrix(const T mas[N][M]) //Конструктор матрицы из массива
		{
			m_n = N;
			m_m = M;
			for (int i = 0; i < m_n; i++)
				for (int j = 0; j < m_m; j++)
					m_mat[i][j] = mas[i][j];
		}

		Matrix(const MasWrapper<T, N, M>& mas) //Конструктор матрицы из структуры типа "MasWrapper"
		{
			m_n = N;
			m_m = M;
			for (int i = 0; i < m_n; i++)
				for (int j = 0; j < m_m; j++)
					m_mat[i][j] = mas.mas[i][j];
		}

		Matrix(const Matrix<T, N, M>& mat) //Конструктор матрицы из матрицы
		{
			m_n = mat.m_n;
			m_m = mat.m_m;

			for (int i = 0; i < m_n; i++)
				for (int j = 0; j < m_m; j++)
					m_mat[i][j] = mat.m_mat[i][j];
		}

		int getN() const { return m_n; } //Получить количество строк матрицы
		int getM() const { return m_m; } //Получить количество столбцов матрицы
		T get(int i, int j) const { return m_mat[i][j]; } //Получить значение i-ой строки и j-ого столбца матрицы
		void set(int i, int j, T data) { m_mat[i][j] = data; }  //Установить значение "data" в i-ой строке и j-ом столбце матрицы

		Matrix<T, N, M>& operator=(const Matrix<T, N, M>& mat) //Перегрузка оператора присваивания матриц
		{
			m_n = mat.m_n;
			m_m = mat.m_m;

			for (int i = 0; i < m_n; i++)
				for (int j = 0; j < m_m; j++)
					m_mat[i][j] = mat.m_mat[i][j];

			return *this;
		}

		Matrix<T, N, M> operator+(const Matrix<T, N, M>& mat) //Перегрузка оператора сложения матриц
		{
			Matrix<T, N, M> tmp;
			for (int i = 0; i < m_n; i++)
				for (int j = 0; j < m_m; j++)
					tmp.m_mat[i][j] = m_mat[i][j] + mat.m_mat[i][j];
			return tmp;
		}

		Matrix<T, N, M> operator+=(const Matrix<T, N, M>& mat) //Перегрузка оператора сложения матриц с присваиванием 
		{
			m_n = mat.m_n;
			m_m = mat.m_m;

			for (int i = 0; i < m_n; i++)
				for (int j = 0; j < m_m; j++)
					m_mat[i][j] += mat.m_mat[i][j];

			return *this;
		}

		Matrix<T, N, M> operator-(const Matrix<T, N, M>& mat) //Перегрузка оператора вычитания матриц
		{
			Matrix<T, N, M> tmp;
			for (int i = 0; i < m_n; i++)
				for (int j = 0; j < m_m; j++)
					tmp.m_mat[i][j] = m_mat[i][j] - mat.m_mat[i][j];
			return tmp;
		}

		template<int L>
		Matrix<T, N, L> operator*(const Matrix<T, M, L>& mat) //Перегрузка оператора умножения матриц
		{
			Matrix<T, N, L> tmp;
			for (int i = 0; i < m_n; i++) {
				for (int j = 0; j < mat.getM(); j++) {
					T sum = 0;
					for (int k = 0; k < m_m; k++)
						sum += m_mat[i][k] * mat.get(k, j);
					tmp.set(i, j, sum);
				}
			}
			return tmp;
		}

		Matrix<T, M, N> transpose() //Нахождение транспонированной матрицы
		{
			Matrix<T, M, N> tr_mat;
			for (int i = 0; i < M; i++) {
				for (int j = 0; j < N; j++) {
					tr_mat.set(i, j, m_mat[j][i]);
				}
			}
			return tr_mat;
		}
		
		template<int L>
		Matrix<T, N, M + L> attach(const Matrix<T, N, L>& mat) //Прикрепление матрицы к исходной для дальнейшего решения СЛАУ, либо нахождения обратной матрицы
		{
			Matrix<T, N, M + L> tempMatrix;
			for (int i = 0; i < m_n; i++) {
				for (int j = 0; j < m_m; j++)
					tempMatrix.set(i, j, m_mat[i][j]);
				for (int j = m_m; j < m_m + mat.getM(); j++)
					tempMatrix.set(i, j, mat.get(i, j - m_m));
			}
			return tempMatrix;
		}

		Matrix<double, N, M> toStepMatrix() //Приведение матрицы к ступенчатому виду путём зануления нижнего левого угла
		{
			Matrix<double, N, M> stepMatrix = *this;
			for (int k = 0; k < m_n; k++)
			{
				for (int i = k + 1; i < m_n; i++)
				{
					double K = -stepMatrix.m_mat[i][k] / stepMatrix.m_mat[k][k];
					for (int j = 0; j < m_m; j++)
						stepMatrix.m_mat[i][j] += stepMatrix.m_mat[k][j] * K;
				}
			}
			return stepMatrix;
		}

		Matrix<double, N, M> toCanonicalStepMatrix() //Приведение матрицы к приведённому ступенчатому виду по строкам путём зануления правого верхнего угла, элементы главной диагонали превращаются в единицы
		{
			Matrix<double, N, M> diagonalMatrix = toStepMatrix();
			for (int k = m_n - 1; k > -1; k--)
			{
				T tmp = diagonalMatrix.m_mat[k][k];
				for (int i = m_m - 1; i > -1; i--)
					diagonalMatrix.m_mat[k][i] /= tmp;
				for (int i = k - 1; i > -1; i--)
				{
					double K = -diagonalMatrix.m_mat[i][k] / diagonalMatrix.m_mat[k][k];
					for (int j = m_m - 1; j > -1; j--)
						diagonalMatrix.m_mat[i][j] += diagonalMatrix.m_mat[k][j] * K;
				}
			}
			return diagonalMatrix;
		}

		double determinant() //Нахождение определителя матрицы
		{
			if (m_n != m_m) throw std::exception("This operation is only available for square matrixes.");
			Matrix<double, N, M> tempMatrix = toStepMatrix();
			double det = 1;
			for (int i = 0; i < m_n; i++) det *= tempMatrix.m_mat[i][i];
			return det;
		}

		Matrix<double, N, M> inversion() //Нахождение обратной матрицы
		{
			if (m_n != m_m) throw std::exception("This operation is only available for square matrixes.");
			double det = determinant();
			if (det == 0) throw std::exception("Matrix has no inverse.");
			Matrix<double, N, M> E(1);
			Matrix<double, N, M * 2> tempMatrix = attach(E).toCanonicalStepMatrix();
			Matrix<double, N, M> inverseMatrix;
			for (int i = 0; i < m_n; i++)
				for (int j = 0; j < m_m; j++)
					inverseMatrix.m_mat[i][j] = tempMatrix.get(i, j + m_m);
			return inverseMatrix;
		}

		Matrix<double, N, 1> SLAU(const Matrix<T, N, 1>& mat) //Нахождение решения СЛАУ
		{
			if (m_n != m_m) throw std::exception("This operation is only available for square matrixes.");
			if (determinant() == 0) throw std::exception("The system has no solutions");
			Matrix<double, N, 1> answer;
			Matrix<double, N, M + 1> tempMatrix = attach(mat).toCanonicalStepMatrix();

			for (int i = 0; i < m_n; i++)
				answer.set(i, 0, tempMatrix.get(i, m_m));
			return answer;
		}

		~Matrix() {} //Пустой деструктор

		template<typename T, int N, int M>
		friend std::istream& operator>>(std::istream& os, Matrix<T, N, M>& mat); //Объявление перегрузки оператора ввода матрицы
		template<typename T, int N, int M>
		friend std::ostream& operator<<(std::ostream& os, const Matrix<T, N, M>& mat); //Объявление перегрузки оператора вывода матрицы
	};

	template<typename T, int N, int M>
	std::istream& operator>>(std::istream& in, Matrix<T, N, M>& mat) //Перегрузка оператора ввода матрицы
	{
		for (int i = 0; i < mat.m_n; i++)
			for (int j = 0; j < mat.m_m; j++)
				in >> mat.m_mat[i][j];
		return in;
	}

	template<typename T, int N, int M>
	std::ostream& operator<<(std::ostream& out, const Matrix<T, N, M>& mat) //Перегрузка оператора вывода матрицы
	{
		out << "Matrix " << mat.m_n << "x" << mat.m_m << std::endl;
		for (int i = 0; i < mat.m_n; i++) {
			for (int j = 0; j < mat.m_m; j++)
			{
				out.width(5);
				out.precision(3);
				out << round(mat.m_mat[i][j] * 1000) / 1000 << " ";
			}
			out << std::endl;
		}
		return out;
	}

	using Vec2i = Matrix<int, 2, 1>;
	using Vec2d = Matrix<double, 2, 1>;
	using Mat22i = Matrix<int, 2, 2>;
	using Mat22d = Matrix<double, 2, 2>;
	using Mat33i = Matrix<int, 3, 3>;
	using Mat33d = Matrix<double, 3, 3>;
	using Mat44i = Matrix<int, 4, 4>;
	using Mat44d = Matrix<double, 4, 4>;
	using Mat55i = Matrix<int, 5, 5>;
	using Mat55d = Matrix<double, 5, 5>;
}