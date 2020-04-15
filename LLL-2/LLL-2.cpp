// LLL-2.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include<vector>
#include <numeric>
#include "LLL-2.h"

using namespace std;
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

void printVS(const Matrix<_type>& vspace) {
	for (auto& v : vspace) {
		cout << '|';
		for (auto& d : v)
			cout << d << ' ';
		cout << "|\n";
	}
}
void printVS(const vector<_type>& v) {
		cout << endl;
		for (auto& d : v)
			cout << d << ' ';
}


float mu(int k1, int k2, Matrix<_type>* m, Matrix<_type>* gsm) {
    vector<_type> *a = &m->at(k1);
    vector<_type> *b = &gsm->at(k2);
    return std::inner_product(a->begin(), a->end(), b->begin(), 0.0) / std::inner_product(b->begin(), b->end(), b->begin(), 0.0);
}
float mu(vector<_type>& a, vector<_type>& b) {
	double d = std::inner_product(a.begin(), a.end(), b.begin(), 0.0);
	double dd = std::inner_product(b.begin(), b.end(), b.begin(), 0.0);
	return d / dd;
}
float squaredLengthOfVector(vector<_type> *v) {
    return std::inner_product(v->begin(), v->end(), v->begin(), 0.0);
}
vector<_type> reduce(vector<_type>& gsreducing, vector<_type>& target, vector<_type>& reducing) {
	vector<_type> t =  target - round(mu(target, gsreducing)) * reducing;
	return t;
}
int main()
{
    int currReduced = 1, currRedusing = 0;// текущий вектор
    Matrix<_type> matrix = { {201 , 37},
                      {1648,  297 } };
	bool cont = true;
	Matrix<_type> gramSmidt;

	while (cont)
	{
		gramSmidt = gSchmidt(matrix);
		printVS(gramSmidt);
		matrix.at(currReduced) = reduce(gramSmidt.at(currRedusing), matrix.at(currReduced), matrix.at(currRedusing));
		float a1 = squaredLengthOfVector(&gramSmidt.at(0));
		float a2 = squaredLengthOfVector(&gramSmidt.at(1));
		float mu1 = mu(1, 0, &matrix, &gramSmidt);
		if (a2 < (_c - mu1 * mu1) * a1) {
			matrix.at(0).swap(matrix.at(1));
		}
		else {
			cont = false;
		}
	}
	printVS(matrix);
    cout << mu(1, 0, &matrix, &matrix);
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
