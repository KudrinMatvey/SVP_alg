// LLL-2.cpp : Этот файл содержит функцию "main". Здесь начинается и
// заканчивается выполнение программы.
//
#include <random>

#include <iostream>

#include <vector>

#include <map>

#include <numeric>

#include "matplotlibcpp.h"

#include "LLL-2.h"

using namespace std;
namespace plt = matplotlibcpp;
template <class T> using Matrix = vector<vector<T>>;
using _type = float;
constexpr double _c = 0.75;


void getCofactor(Matrix<_type>& A, Matrix<_type>& temp, int p, int q, int n)
{
    int i = 0, j = 0;

    // Looping for each element of the matrix 
    for (int row = 0; row < n; row++)
    {
        for (int col = 0; col < n; col++)
        {
            //  Copying into temporary matrix only those element 
            //  which are not in given row and column 
            if (row != p && col != q)
            {
                temp[i][j++] = A[row][col];

                // Row is filled, so increase row index and 
                // reset col index 
                if (j == n - 1)
                {
                    j = 0;
                    i++;
                }
            }
        }
    }
}

/* Recursive function for finding determinant of matrix.
   n is current dimension of A[][]. */
_type determinant(Matrix<_type>& A,  int n)
{

    //  Base case : if matrix contains single element 
    if (n == 1)
        return A[0][0];
    _type D = 0; // Initialize result 
    Matrix<_type> temp(n, vector<_type> (n , 0)); // To store cofactors 

    float sign = 1;  // To store sign multiplier 

     // Iterate for each element of first row 
    for (int f = 0; f < n; f++)
    {
        // Getting Cofactor of A[0][f] 
        getCofactor(A, temp, 0, f, n);
        D += sign * A[0][f] * determinant(temp, n - 1);

        // terms are to be added with alternate sign 
        sign = -sign;
    }

    return D;
}

// Function to get adjoint of A[N][N] in adj[N][N]. 
void adjoint(Matrix<_type>& A, Matrix<_type>& adj)
{
    int n = A.size();
    if (n== 1)
    {
        adj[0][0] = 1;
        return;
    }

    // temp is used to store cofactors of A[][] 
    _type sign = 1;
    Matrix<_type> temp(n, vector<_type>(n, 0));

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            // Get cofactor of A[i][j] 
            getCofactor(A, temp, i, j, n);

            // sign of adj[j][i] positive if sum of row 
            // and column indexes is even. 
            sign = ((i + j) % 2 == 0) ? 1 : -1;

            // Interchanging rows and columns to get the 
            // transpose of the cofactor matrix 
            adj[j][i] = (sign) * (determinant(temp, n - 1));
        }
    }
}

// Function to calculate and store inverse, returns false if 
// matrix is singular 
bool inverse(Matrix<_type>& A, Matrix<_type>& inverse)
{
    // Find determinant of A[][] 
    int det = determinant(A, A.size());
    //_type det = 6;
    if (det == 0)
    {
        cout << "Singular matrix, can't find its inverse";
        return false;
    }

    // Find adjoint
    Matrix<_type> adj;
    for (int i = 0; i < A.size(); i++) {
        adj.push_back(vector<_type>());
        for (int j = 0; j < A.size(); j++)
            adj[i].push_back(0);
    }
    adjoint(A, adj);

    // Find Inverse using formula "inverse(A) = adj(A)/det(A)" 
    for (int i = 0; i < A.size(); i++)
        for (int j = 0; j < A.size(); j++)
            inverse[i][j] = adj[i][j] / float(det);

    return true;
}


double dotProduct(const vector<_type>& v1, const vector<_type>& v2) {
    double ret = 0;
    for (int i = 0; i < v1.size(); ++i)
        ret += (v1.at(i) * v2.at(i));
    return ret;
}

vector<_type> operator*(const double& c, const vector<_type>& v) {
    vector<_type> ret = v;
    for (auto& d : ret)
        d *= c;
    return ret;
}
vector<_type> operator*(Matrix<_type>& M, const vector<_type>& v) {
    int n = v.size();
    vector<_type> res(n, 0);
    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < n; j++)
        {
            res[i] += M[j][i] * v[j];
        }
    }
    return res;
}

vector<_type> operator-(const vector<_type>& v1, const vector<_type>& v2) {
    vector<_type> ret;
    for (int i = 0; i < v1.size(); ++i)
        ret.push_back(v1.at(i) - v2.at(i));
    return ret;
}
bool operator==(const vector<_type>& v1, const vector<_type>& v2) {
    if (v1.size() == v2.size()) {
        if (v1.size()) {
            if(v1[0] == v2[0])
            { 
                for (int i = 1; i < v1.size(); i++)
                {
                    if (v1[i] != v2[i])
                    {
                        return false;
                    }
                }
                return true;
            }
            if(v1[0] == -v2[0])
            { 
                for (int i = 1; i < v1.size(); i++)
                {
                    if (v1[i] != -v2[i])
                    {
                        return false;
                    }
                }
                return true;
            }
            return false;
        }
        else return true;
    }
    else return false;
}
vector<_type> operator+(const vector<_type>& v1, const vector<_type>& v2) {
    vector<_type> ret;
    for (int i = 0; i < v1.size(); ++i)
        ret.push_back(v1.at(i) + v2.at(i));
    return ret;
}
vector<_type> roundV(vector<_type> a) {
    for (size_t i = 0; i < a.size(); i++)
    {
        a[i] = round(a[i]);
    }
    return a;
}

vector<_type> proj(const vector<_type>& v1, const vector<_type>& v2) {
    vector<_type> proj = (dotProduct(v1, v2) / dotProduct(v1, v1)) * v1;
    return proj;
}

Matrix<_type> gSchmidt(const Matrix<_type>& vspace) {
    Matrix<_type> uspace;
    for (int i = 0; i < vspace.size(); ++i) {
        uspace.push_back(vspace.at(i));
        for (int j = 0; j < i; ++j) {
            uspace.at(i) = uspace.at(i) - proj(uspace.at(j), vspace.at(i));
        }
    }
    return uspace;
}
Matrix<_type> gSchmidt(const Matrix<_type>& vspace, int limit) {
    Matrix<_type> uspace;
    for (int i = 0; i <= limit; ++i) {
        uspace.push_back(vspace.at(i));
        for (int j = 0; j < i; ++j) {
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
    cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<endl;
    for (auto& d : v)
        cout << d << ' ';
    cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;

}

float mu(int k1, int k2, Matrix<_type>* m, Matrix<_type>* gsm) {
    vector<_type>* a = &m -> at(k1);
    vector<_type>* b = &gsm -> at(k2);
    return std::inner_product(a -> begin(), a -> end(), b -> begin(), 0.0) /
        std::inner_product(b -> begin(), b -> end(), b -> begin(), 0.0);
}
float mu(vector<_type>& a, vector<_type>& b) {
    double d = std::inner_product(a.begin(), a.end(), b.begin(), 0.0);
    double dd = std::inner_product(b.begin(), b.end(), b.begin(), 0.0);
    return d / dd;
}
float squaredLengthOfVector(vector<_type>* v) {
    return std::inner_product(v -> begin(), v -> end(), v -> begin(), 0.0);
}
float squaredLengthOfVector(vector<_type>& v) {
    return std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
}
vector<_type> reduce(vector<_type>& gsreducing, vector<_type>& target,
    vector<_type>& reducing) {
    vector<_type> t = target - round(mu(target, gsreducing)) * reducing;
    return t;
}

int KabatianskyLevenshteinC(Matrix<_type>& M) {
    float angle = 1;
    for (int i = 0; i < M.size(); i++) {
        auto v1 = M[i];
        for (int j = i + 1; j < M.size(); j++) {
            auto v2 = M[j];
            float t = dotProduct(v1, v2) / (sqrt(squaredLengthOfVector(&v1)) *
                sqrt(squaredLengthOfVector(&v2)));
            angle = fmin(angle, t);
        }
    }
    return fmax((int)(-1.0 / 2.0 * log(1 - cos(angle)) - 0.099), 1);
}

void LLL(Matrix<_type>& B) {
    bool cont;
    Matrix<_type> gramSmidt;
    do {
        gramSmidt = gSchmidt(B);
        cont = false;
        for (int i = 1; i < B.size(); i++) {
            for (int j = 0; j < i; j++) {
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
            for (int j = -maxSize; j < maxSize; j++) {
                vector<_type> t = i * res.at(0) + j * res.at(1);
                // m.insert({ make_pair(t[0], t[1]) });
                if (-boundX < t[0] < boundX && -boundY < t[1] < boundY) {
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
                for (int k = -maxSize; k < maxSize; j++) {
                    vector<_type> t = i * res.at(0) + j * res.at(1) + k * res.at(2);
                    // m.insert({ make_pair(t[0], t[1]) });
                    if (-boundX < t[0] < boundX && -boundY < t[1] < boundY) {
                        x.push_back(t.at(0));
                        y.push_back(t.at(1));
                        z.push_back(t.at(2));
                    }
                }
        }
     }
    plt::plot({ 0 }, { 0 }, "bo");
    plt::plot({ res.at(0).at(0) }, { res.at(0).at(1) }, "r+");
    plt::plot({ res.at(1).at(0) }, { res.at(1).at(1) }, "r+");
    plt::plot({ og.at(0).at(0) }, { og.at(0).at(1) }, "g*");
    plt::plot({ og.at(1).at(0) }, { og.at(1).at(1) }, "g*");
}
//https://tel.archives-ouvertes.fr/tel-01245066v2/document
vector<_type> kleinSample(Matrix<_type>& B, Matrix<_type>& GSO, float deviation, vector<_type>& c) {
    vector<_type> v(B.size(), 0);
    std::random_device mch;
    std::default_random_engine generator(mch());

    //vector<_type> z(B.size());
    double d, z;
    for (int i = B.size() - 1; i >= 0; i--) {
        double length = squaredLengthOfVector(GSO.at(i));
        d = dotProduct(c, GSO.at(i)) / length ;
        double sigma = deviation / sqrt(length);
        std::normal_distribution<double> distribution(d, 1);
        //randomised rounding
        //z =  rand() / RAND_MAX > sigma ? d : 0;
        z = round(distribution(generator));
        std::cout << z << " " ;
        c = c - z * B[i];
        v = v + z * B[i];
    }
    std::cout << std::endl;
    return v;
}

// todo discover
vector<_type> sample_gaussian(Matrix<_type>& B) {
    vector<_type> a;
    a = (5 - rand() % 10) * B.back();
    for (auto b : B)
        a = a + (5 - rand() % 10) * b;
    return a;
}
vector<_type> gauss_reduce(vector<_type>& v_new, Matrix<_type>& L,
    Matrix<_type>& S) {
    bool cont = true;
    while (cont) {
        auto vn = std::find_if(L.begin(), L.end(), [&v_new](vector<_type>& v) {
            return squaredLengthOfVector(&v) <= squaredLengthOfVector(&v_new) &&
                squaredLengthOfVector(&(v_new - v)) <=
                squaredLengthOfVector(&v_new);
            });
        if (vn != L.end()) {
            v_new = v_new - *vn;
        }
        else {
            cont = false;
        }
    }
    auto newEnd =
        std::remove_if(L.begin(), L.end(), [&v_new, &S](vector<_type>& v) {
        if (squaredLengthOfVector(&v) > squaredLengthOfVector(&v_new) &&
            squaredLengthOfVector(&(v - v_new)) <= squaredLengthOfVector(&v)) {
            S.push_back(v - v_new);
            return true;
        }
        return false;
            });
    L.erase(newEnd, L.end());

    return v_new;
}
// https://cseweb.ucsd.edu/~daniele/papers/Sieve.pdf
vector<_type> gaussSieve(Matrix<_type>& B) {
    Matrix<_type> L, S;
    vector<_type> v_new;
    int K = 0;
    int c = KabatianskyLevenshteinC(B);
    do {
        if (S.size() > 0) {
            v_new = S.back();
            S.pop_back();
        }
        else {
            //v_new = sample_gaussian(B);
            Matrix<_type> gramSmidt = gSchmidt(B);
            vector<_type> a(3, 0);
            v_new = kleinSample(B, gramSmidt, 1, a);
        }
        v_new = gauss_reduce(v_new, L, S);
        if (!v_new.size() ||
            std::count(v_new.begin(), v_new.end(), 0) == v_new.size()) {
            K++;
        }
        else {
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

vector<_type> randomInBall(const int dimm, float r) {
    vector<_type> randVector(dimm);
    for (int i = 0; i < dimm; i++) {
        randVector[i] = (rand() % (int)(r + 1));
    }
    return (r * (rand()/double(RAND_MAX)) / sqrt(squaredLengthOfVector(randVector))) * randVector;
}
//https://crypto.stackexchange.com/questions/29661/how-to-find-the-value-of-a-vector-modulo-a-basis-in-lattice-based-cryptography/29701#29701
vector<_type> mod(vector<_type> vec, Matrix<_type> m) {
    Matrix<_type> r(m.size(), vector<_type>(m.size(), 0));
    inverse(m, r);
    return m * roundV(r * vec);
}
void ListReduce(vector<_type> vec, Matrix<_type> L, float delta) {
    bool cont = true;
    while (cont) {
        auto f = std::find_if(L.begin(), L.end(), [&vec, &delta](vector<_type>& l) {
            return squaredLengthOfVector(&(vec - l)) <= delta * delta * squaredLengthOfVector(&vec);
            });
        if (f == L.end()) return;
        vec = vec - *f;
    }
}
void ListReduce2(vector<_type> p, Matrix<_type> L) {
    bool cont = true;
    while (cont) {
        auto f = std::find_if(L.begin(), L.end(), [&p](vector<_type>& v) {
            return squaredLengthOfVector(&(p - v)) <= squaredLengthOfVector(p) &&
                squaredLengthOfVector(p ) >= squaredLengthOfVector(v);
            });
        if (f == L.end()) return;
        p = p - *f;
    }
}
//https://math.stackexchange.com/questions/396382/what-does-this-dollar-sign-over-arrow-in-function-mapping-mean
//vector<_type> ListSieve(Matrix<_type>& B, float mu) {
//    vector<vector<_type>> L;
//    float delta = 1 - (1 / B.size());
//    float eps = 0.685;
//    int K = pow(2, KabatianskyLevenshteinC(B) * B.size());
//    for (int i = 0; i < K;++i) {
////        vector<_type> e = randomInBall(B.size(), eps * mu);
////        vector<_type> p = mod(e, B);
//        vector<_type> v = sample_gaussian(B);
//        ListReduce(v, L, delta);
////        vector<_type> v = p - e;
//        auto k = std::find(L.begin(), L.end(), v);
//        if (k != L.end() || L.size() == 0) {
//            k = std::find_if(L.begin(), L.end(), [&v, &mu](vector<_type>& vj) {return squaredLengthOfVector(&(v - vj)) < mu * mu; });
//            if (k != L.end()) {
//                return *k;
//            }
//            L.push_back(v);
//        }   
//    }
//    vector<_type> shortestVector = L[0];
//    for (int i = 1; i < L.size(); i++) {
//        if (squaredLengthOfVector(shortestVector) > squaredLengthOfVector(L[i])) {
//            shortestVector = L[i];
//        }
//    }
//    return shortestVector;
//}
//https://eprint.iacr.org/2014/714.pdf
vector<_type> ListSieve2(Matrix<_type>& B, float mu) {
    vector<vector<_type>> L;
    float delta = 1 - (1 / B.size());
    float eps = 0.685;
    Matrix<_type> gramSmidt = gSchmidt(B);
    vector<_type> a(3, 0);
    // c?
    //int K = pow(2, KabatianskyLevenshteinC(B) * B.size());
    int K = 3;
    for (int i = 0; i < K;) {
        vector<_type> v = kleinSample(B, gramSmidt, 1, a);
        ListReduce2(v, L);
        bool isZero = 1;
        for (auto i : v) {
            if (i != 0) {
                isZero = false;
                break;
            }
        }
        if (isZero) {
            i++;
        }
        else {
            L.push_back(v);
        }
    }
    vector<_type> shortestVector = L[0];
    for (int i = 1; i < L.size(); i++) {
        if (squaredLengthOfVector(shortestVector) > squaredLengthOfVector(L[i])) {
            shortestVector = L[i];
        }
    }
    return shortestVector;
}

int main() {
    vector<int> v_new;
    cout << (bool)(v_new.begin() == v_new.end());

    int currReduced = 1, currRedusing = 0; // текущий вектор
                                           /*{ {201 , 37},
                                                                                   {1648,  297 } };*/
    Matrix<_type> matrix = { {15, 23, 11},
                             {32, 1, 1},
                             {46, 15, 3} };
    randomInBall(15, 5);
    Matrix<_type> matrix2 = matrix;
    Matrix<_type> matrix3 = matrix;
    Matrix<_type> original(matrix2);
    Matrix<_type> gramSmidt = gSchmidt(matrix);
    //LLL(matrix);
    printVS(matrix);
    vector<_type> sv = gaussSieve(matrix2);
    float shortestVectorLength = sqrt(squaredLengthOfVector(sv));
    printVec(sv);
    float mu = shortestVectorLength;
    printVec(ListSieve2(matrix3, mu));
   /* for (int i = 0; i < 10; i++) {
        vector<_type> a(3, 0);
        printVec(kleinSample(matrix, gramSmidt, 1, a));
    }*/
    }

// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"

// Советы по началу работы
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и
//   другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый
//   элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий
//   элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" >
//   "Открыть" > "Проект" и выберите SLN-файл.