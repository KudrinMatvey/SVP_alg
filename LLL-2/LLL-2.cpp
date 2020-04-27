// LLL-2.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include<vector>
#include<map>
#include <numeric>
#include "matplotlibcpp.h"
#include "LLL-2.h"

using namespace std;
namespace plt = matplotlibcpp;
template<class T> using Matrix = vector<vector<T>>;
using _type = float;
constexpr double _c = 0.75;


double dotProduct(const vector<_type>& v1, const vector<_type>& v2) {
	double ret = 0;
	for (int i = 0; i < v1.size(); ++i)
		ret += (v1.at(i) * v2.at(i));
	return ret;
}

vector<_type> operator * (const double& c, const vector<_type>& v) {
	vector<_type> ret = v;
	for (auto& d : ret)
		d *= c;
	return ret;
}

vector<_type> operator - (const vector<_type>& v1, const vector<_type>& v2) {
	vector<_type> ret;
	for (int i = 0; i < v1.size(); ++i)
		ret.push_back(v1.at(i) - v2.at(i));
	return ret;
}
vector<_type> operator + (const vector<_type>& v1, const vector<_type>& v2) {
	vector<_type> ret;
	for (int i = 0; i < v1.size(); ++i)
		ret.push_back(v1.at(i) + v2.at(i));
	return ret;
}

vector<_type> proj(const vector<_type>& v1, const vector<_type>& v2) {
	vector<_type> proj = (dotProduct(v1, v2) / dotProduct(v1, v1)) * v1;
	return proj;
}


Matrix<_type> gSchmidt(const Matrix<_type>& vspace) {
	Matrix<_type> uspace;
	for (int i = 0; i < vspace.size(); ++i)
	{
		uspace.push_back(vspace.at(i));
		for (int j = 0; j < i; ++j)
		{
			uspace.at(i) = uspace.at(i) - proj(uspace.at(j), vspace.at(i));
		}

	}
	return uspace;
}
Matrix<_type> gSchmidt(const Matrix<_type>& vspace, int limit) {
	Matrix<_type> uspace;
	for (int i = 0; i <= limit; ++i)
	{
		uspace.push_back(vspace.at(i));
		for (int j = 0; j < i; ++j)
		{
			uspace.at(i) = uspace.at(i) - proj(uspace.at(j), vspace.at(i));
		}

	}
	return uspace;
}

void printVS(const Matrix<_type>& vspace) {
	for (auto& v : vspace) {
		cout << '|';
		for (auto& d : v)
			cout << d << ' ';
		cout << "|\n";
	}
}
void printVS(const vector<_type>& vspace) {
		cout << '|';
		for (auto& d : vspace)
			cout << d << ' ';
		cout << "|\n";
}
void printVec(const vector<_type>& v) {
	cout << endl;
	for (auto& d : v)
		cout << d << ' ';
}


float mu(int k1, int k2, Matrix<_type>* m, Matrix<_type>* gsm) {
	vector<_type>* a = &m->at(k1);
	vector<_type>* b = &gsm->at(k2);
	return std::inner_product(a->begin(), a->end(), b->begin(), 0.0) / std::inner_product(b->begin(), b->end(), b->begin(), 0.0);
}
float mu(vector<_type>& a, vector<_type>& b) {
	double d = std::inner_product(a.begin(), a.end(), b.begin(), 0.0);
	double dd = std::inner_product(b.begin(), b.end(), b.begin(), 0.0);
	return d / dd;
}
float squaredLengthOfVector(vector<_type>* v) {
	return std::inner_product(v->begin(), v->end(), v->begin(), 0.0);
}
float squaredLengthOfVector(vector<_type>& v) {
	return std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
}
vector<_type> reduce(vector<_type>& gsreducing, vector<_type>& target, vector<_type>& reducing) {
	vector<_type> t = target - round(mu(target, gsreducing)) * reducing;
	return t;
}

int KabatianskyLevenshteinC(Matrix<_type> &M) {
	float angle = 1;
	for (int i = 0; i < M.size(); i++)
	{
		auto v1 = M[i];
		for (int j = i + 1; j < M.size(); j++)
		{
			auto v2 = M[j];
			float t = dotProduct(v1, v2) / (sqrt(squaredLengthOfVector(&v1)) * sqrt(squaredLengthOfVector(&v2)));
			angle = fmin(angle, t);
		}
	}
	return fmax((int)(- 1.0 / 2.0 * log(1 - cos(angle)) - 0.099), 1);
}

void LLL(Matrix<_type>& B) {
	bool cont;
	Matrix<_type> gramSmidt;
	do {
		gramSmidt = gSchmidt(B);
		cont = false;
		for (int i = 1; i < B.size(); i++) {
			for (int j = 0; j < i; j++)
			{
				B.at(i) = reduce(gramSmidt.at(j), B.at(i), B.at(j));
			}
			float a1 = squaredLengthOfVector(&gramSmidt.at(i - 1));
			float a2 = squaredLengthOfVector(&gramSmidt.at(i));
			float mu1 = mu(i, i - 1, &B, &gramSmidt);
			if (a2 < (_c - mu1 * mu1) * a1) {
				B.at(i - 1).swap(B.at(i));
				cont = true;
			}
		}
	} while (cont);
}
void preparePlot(Matrix<_type> res, Matrix<_type> og) {
	map<_type, _type> m;
	vector<_type> x, y, z;
	int maxSize = 20;
	if (res.size() == 2) {
		vector<_type> bounds = 2 * og[0] + 2 * og[1];
		int boundX = bounds[0];
		int boundY = bounds[1];
		for (int i = -maxSize; i < maxSize; i++) {
			for (int j = -maxSize; j < maxSize; j++)
			{
				vector<_type> t = i * res.at(0) + j * res.at(1);
				//m.insert({ make_pair(t[0], t[1]) });
				if (-boundX < t[0] < boundX && -boundY < t[1] < boundY)
				{
					x.push_back(t.at(0));
					y.push_back(t.at(1));
				}
			}
		}
		plt::plot(x, y, ".");
	}
	if (res.size() == 3) {
		vector<_type> bounds = 2 * og[0] + 2 * og[1];
		int boundX = bounds[0];
		int boundY = bounds[1];
		int boundZ = bounds[2];
		for (int i = -maxSize; i < maxSize; i++) {
			for (int j = -maxSize; j < maxSize; j++)
				for (int k = -maxSize; k < maxSize; j++)
				{
					vector<_type> t = i * res.at(0) + j * res.at(1) + k * res.at(2);
					//m.insert({ make_pair(t[0], t[1]) });
					if (-boundX < t[0] < boundX && -boundY < t[1] < boundY)
					{
						x.push_back(t.at(0));
						y.push_back(t.at(1));
						z.push_back(t.at(2));
					}
				}
		}
		//plt::plot_surface(x, y, z);
	}

	//m.insert({ make_pair(t[0], t[1]) });
	//for (auto& v : m)
	//{
	//	x.push_back(v.first);
	//	y.push_back(v.second);
	//}
	plt::plot({ 0 }, { 0 }, "bo");
	plt::plot({ res.at(0).at(0) }, { res.at(0).at(1) }, "r+");
	plt::plot({ res.at(1).at(0) }, { res.at(1).at(1) }, "r+");
	plt::plot({ og.at(0).at(0) }, { og.at(0).at(1) }, "g*");
	plt::plot({ og.at(1).at(0) }, { og.at(1).at(1) }, "g*");
}
//todo discover
vector<_type> sample_gaussian(Matrix<_type>& B) {
	vector<_type> a;
	a = (5 - rand() % 10) * B.back();
	for (auto b : B)
		a = a + (5 - rand() % 10) * b;
	return a;
}
vector<_type> gauss_reduce(vector<_type>& v_new, Matrix<_type>& L, Matrix<_type>& S) {
	bool cont = true;
	while (cont) {
		auto vn = std::find_if(L.begin(), L.end(), [&v_new](vector<_type>& v) {
			return squaredLengthOfVector(&v) <= squaredLengthOfVector(&v_new)
				&& squaredLengthOfVector(&(v_new - v)) <= squaredLengthOfVector(&v_new);
			});
		if (vn != L.end()) {
			v_new = v_new - *vn;
		}
		else
		{
			cont = false;
		}
	}
	auto newEnd = std::remove_if(L.begin(), L.end(), [&v_new, &S](vector<_type>& v) {
		if (squaredLengthOfVector(&v) > squaredLengthOfVector(&v_new)
			&& squaredLengthOfVector(&(v - v_new)) <= squaredLengthOfVector(&v)) {
			S.push_back(v - v_new);
			return true;
		} return false;
		});
	L.erase(newEnd, L.end());

	return v_new;
}
//https://cseweb.ucsd.edu/~daniele/papers/Sieve.pdf
vector<_type> gaussSieve(Matrix<_type> &B) {
	Matrix<_type> L, S;
	vector<_type> v_new;
	int K = 0;
	int c = KabatianskyLevenshteinC(B);
	do
	{
		if (S.size() > 0) {
			v_new = S.back();
			S.pop_back();
		}
		else
		{
			v_new = sample_gaussian(B);
		}
		v_new = gauss_reduce(v_new, L, S);
		if (!v_new.size() || std::count(v_new.begin(), v_new.end() + 1, 0) == v_new.size()) {
			K++;
		}
		else
		{
			L.push_back(v_new);
		}
	} while (K < c);
	vector<_type> shortestVector = S[0];
	for (int i = 1; i < S.size(); i++) {
		if (squaredLengthOfVector(shortestVector) > squaredLengthOfVector(S[i])) {
			shortestVector = S[i];
		}
	}
	return shortestVector;
}
int main()
{
	vector<int> v_new;
	cout << (bool)(v_new.begin() == v_new.end());

	int currReduced = 1, currRedusing = 0; // текущий вектор
	 /*{ {201 , 37},
					  {1648,  297 } };*/
	Matrix<_type> matrix = { {15 , 23, 11},
							{32,1,1},
							{46,  15, 3} };
	Matrix<_type> matrix2 = matrix;
	Matrix<_type> original(matrix2);
	Matrix<_type> gramSmidt;
	//for (; currReduced < matrix.size(); currReduced++) {
	//	bool cont = true;
	//	for (; currRedusing < currReduced; currRedusing++) {
	//		while (cont)
	//		{
	//			gramSmidt = gSchmidt(matrix, currReduced);
	//			//printVS(gramSmidt);
	//			matrix.at(currReduced) = reduce(gramSmidt.at(currRedusing), matrix.at(currReduced), matrix.at(currRedusing));
	//			float a1 = squaredLengthOfVector(&gramSmidt.at(currRedusing));
	//			float a2 = squaredLengthOfVector(&gramSmidt.at(currReduced));
	//			float mu1 = mu(1, 0, &matrix, &gramSmidt);
	//			if (a2 < (_c - mu1 * mu1) * a1) {
	//				matrix.at(0).swap(matrix.at(1));
	//			}
	//			else {
	//				cont = false;
	//			}
	//		}
	//	}
	//}
	LLL(matrix);
	printVS(matrix);
	printVec(gaussSieve(matrix2));
	//preparePlot(matrix2, original);
	//plt::show();
	//cout << mu(1, 0, &matrix, &matrix);
}

// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"

// Советы по началу работы 
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.
