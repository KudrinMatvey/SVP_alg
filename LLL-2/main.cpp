#include <Eigen/Dense>
#include <random>
#include <iostream>
#include <vector>
#include <map>
#include <numeric>
#include <fstream>
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <cmath>

using namespace std;
//namespace plt = matplotlibcpp;
template <class T>
using Matrix = vector<vector<T>>;
using _type = float;
constexpr double _c = 0.75;

template <typename T>
void getCofactor(Matrix<T>& A, Matrix<T>& temp, int p, int q, int n)
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
template <typename T>
T determinant(Matrix<T>& A, int n)
{

    //  Base case : if matrix contains single element
    if (n == 1)
        return A[0][0];
    T D = 0;                            // Initialize result
    Matrix<T> temp(n, vector<T>(n, 0)); // To store cofactors

    float sign = 1; // To store sign multiplier

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
    if (n == 1)
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
    for (int i = 0; i < A.size(); i++)
    {
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

template <typename T, typename T2>
double dotProduct(const vector<T>& v1, const vector<T2>& v2)
{
    double ret = 0;
    for (int i = 0; i < v1.size(); ++i)
        ret += (v1.at(i) * v2.at(i));
    return ret;
}

template <typename T>
vector<T> operator*(const double& c, const vector<T>& v)
{
    vector<T> ret = v;
    for (auto& d : ret)
        d *= c;
    return ret;
}

vector<_type> operator*(Matrix<_type>& M, const vector<_type>& v)
{
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
template <typename T, typename T2>
vector<T2> subtract(const vector<T>& v1, const vector<T2>& v2)
{
    vector<T2> ret(v1.size());
    for (int i = 0; i < v1.size(); ++i)
        ret[i] = (v1[i] - v2[i]);
    return ret;
}
template <typename T, typename T2>
vector<T> operator-(const vector<T>& v1, const vector<T2>& v2)
{
    vector<T> ret(v1.size());
    for (int i = 0; i < v1.size(); ++i)
        ret[i] = v1[i] - v2[i];
    return ret;
}
template <typename T>
vector<T> inverse(const vector<T>& v1)
{
    vector<T> ret(v1.size());
    for (int i = 0; i < v1.size(); ++i)
        ret[i] = -v1[i];
    return ret;
}
template <typename T>
bool operator==(const vector<T>& v1, const vector<T>& v2)
{
    if (v1.size() == v2.size())
    {
        if (v1.size())
        {
            if (v1[0] == v2[0])
            {
                if (v1[0] == 0)
                {
                    int sign = 0;
                    for (int i = 1; i < v1.size(); i++)
                    {
                        if (v1[i] == 0 && v2[i] == 0)
                            continue;
                        if (sign == 0)
                        {
                            if (v1[i] == v2[i])
                                sign = 1;
                            else if (v1[i] == -v2[i])
                                sign = -1;
                            else
                                return false;
                        }
                        else if (v1[i] != sign * v2[i])
                        {
                            return false;
                        }
                    }
                }
                else
                {
                    for (int i = 1; i < v1.size(); i++)
                    {
                        if (v1[i] != v2[i])
                        {
                            return false;
                        }
                    }
                }
                return true;
            }
            if (v1[0] == -v2[0])
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
        else
            return true;
    }
    else
        return false;
}
template <typename T>
vector<T> operator+(const vector<T>& v1, const vector<T>& v2)
{
    vector<T> ret(v1.size());
    for (int i = 0; i < v1.size(); ++i)
        ret[i] = (v1.at(i) + v2.at(i));
    return ret;
}
vector<_type> roundV(vector<_type> a)
{
    for (size_t i = 0; i < a.size(); i++)
    {
        a[i] = round(a[i]);
    }
    return a;
}

vector<_type> proj(const vector<_type>& v1, const vector<_type>& v2)
{
    vector<_type> proj = (dotProduct(v1, v2) / dotProduct(v1, v1)) * v1;
    return proj;
}

vector<_type> proj(const vector<_type>& v1, const vector<int>& v2)
{
    vector<_type> proj = (dotProduct(v1, v2) / dotProduct(v1, v1)) * v1;
    return proj;
}

Matrix<_type> gSchmidt(const Matrix<int>& vspace)
{
    Matrix<_type> uspace(vspace.size(), vector<_type>(vspace[0].size()));
    for (int i = 0; i < vspace[0].size(); i++)
    {
        uspace[0][i] = vspace[0][i];
    }
    for (int i = 1; i < vspace.size(); ++i)
    {
        for (int j = 0; j < i; ++j)
        {
            uspace.at(i) = subtract(vspace.at(i), proj(uspace.at(j), vspace.at(i)));
        }
    }
    return uspace;
}
Matrix<_type> gSchmidt(const Matrix<_type>& vspace, int limit)
{
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
template <typename T>
void printVS(const Matrix<T>& vspace)
{
    for (auto& v : vspace)
    {
        cout << '|';
        for (auto& d : v)
            cout << d << ' ';
        cout << "|\n";
    }
}
template <typename T>
void printVS(const Matrix<T>& vspace, ostream& file)
{
    file << '[';
    for (auto& v : vspace)
    {
        file << '[';
        for (auto& d : v)
            file << d << ' ';
        file << "]\n";
    }
    file << "]\n";
}
template <typename T>
void printVec(const vector<T>& v)
{
    cout << endl;
    cout << "$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
    for (auto& d : v)
        cout << d << ' ';
    cout << "$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
}
template <typename T>
void printVec(const vector<T>& v, ostream& file)
{
    file << endl;
    file << "$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
    for (auto& d : v)
        file << d << ' ';
    file << endl
        << "$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
}
template <typename T>
float mu(int k1, int k2, Matrix<T>* m, Matrix<_type>* gsm)
{
    T init = 0;
    vector<T>* a = &m->at(k1);
    vector<_type>* b = &gsm->at(k2);
    return std::inner_product(a->begin(), a->end(), b->begin(), init) /
        std::inner_product(b->begin(), b->end(), b->begin(), init);
}
template <typename T, typename T2>
float mu(vector<T>& a, vector<T2>& b)
{
    float init = 0;
    double d = std::inner_product(a.begin(), a.end(), b.begin(), init);
    double dd = std::inner_product(b.begin(), b.end(), b.begin(), init);
    return d / dd;
}
template <typename T>
double squaredLengthOfVector(vector<T>* v)
{
    return std::inner_product(v->begin(), v->end(), v->begin(), 0.0);
}
template <typename T>
double squaredLengthOfVector(vector<T>& v)
{
    return std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
}
vector<int> lll_reduce(vector<_type>& gsreducing, vector<int>& target,
    vector<int>& reducing)
{
    vector<int> t = target - round(mu(target, gsreducing)) * reducing;
    return t;
}

int KabatianskyLevenshteinC(Matrix<int>& M)
{
    float angle = 1;
    for (int i = 0; i < M.size(); i++)
    {
        auto v1 = M[i];
        for (int j = i + 1; j < M.size(); j++)
        {
            auto v2 = M[j];
            float t = dotProduct(v1, v2) / (sqrt(squaredLengthOfVector(&v1)) *
                sqrt(squaredLengthOfVector(&v2)));
            angle = fmin(angle, t);
        }
    }
    return fmax((int)(-1.0 / 2.0 * log(1 - cos(angle)) - 0.099), 1);
}

void LLL(Matrix<int>& B, ostream& o)
{
    bool cont;
    Matrix<_type> gramSmidt;
    int iter = 0;
    do
    {
        iter++;
        gramSmidt = gSchmidt(B);
        cont = false;
        for (int i = 1; i < B.size(); i++)
        {
            for (int j = 0; j < i; j++)
            {
                B.at(i) = lll_reduce(gramSmidt.at(j), B.at(i), B.at(j));
            }
            float a1 = squaredLengthOfVector(&gramSmidt.at(i - 1));
            float a2 = squaredLengthOfVector(&gramSmidt.at(i));
            float mu1 = mu(i, i - 1, &B, &gramSmidt);
            if (a2 < (_c - mu1 * mu1) * a1)
            {
                B.at(i - 1).swap(B.at(i));
                cont = true;
            }
        }
    } while (cont);
    o << "lll iter : " << iter << "|\n";
}

void LLL(Matrix<int>& B)
{
    bool cont;
    Matrix<_type> gramSmidt;
    int iter = 0;
    do
    {
        iter++;
        gramSmidt = gSchmidt(B);
        cont = false;
        for (int i = 1; i < B.size(); i++)
        {
            for (int j = 0; j < i; j++)
            {
                B.at(i) = lll_reduce(gramSmidt.at(j), B.at(i), B.at(j));
            }
            float a1 = squaredLengthOfVector(&gramSmidt.at(i - 1));
            float a2 = squaredLengthOfVector(&gramSmidt.at(i));
            float mu1 = mu(i, i - 1, &B, &gramSmidt);
            if (a2 < (_c - mu1 * mu1) * a1)
            {
                B.at(i - 1).swap(B.at(i));
                cont = true;
            }
        }
    } while (cont);
}
void preparePlot(Matrix<_type> res, Matrix<_type> og)
{
    map<_type, _type> m;
    vector<_type> x, y, z;
    int maxSize = 20;
    if (res.size() == 2)
    {
        vector<_type> bounds = 2 * og[0] + 2 * og[1];
        int boundX = bounds[0];
        int boundY = bounds[1];
        for (int i = -maxSize; i < maxSize; i++)
        {
            for (int j = -maxSize; j < maxSize; j++)
            {
                vector<_type> t = i * res.at(0) + j * res.at(1);
                if (-boundX < t[0] < boundX && -boundY < t[1] < boundY)
                {
                    x.push_back(t.at(0));
                    y.push_back(t.at(1));
                }
            }
        }
        //plt::plot(x, y, ".");
    }
    if (res.size() == 3)
    {
        vector<_type> bounds = 2 * og[0] + 2 * og[1];
        int boundX = bounds[0];
        int boundY = bounds[1];
        int boundZ = bounds[2];
        for (int i = -maxSize; i < maxSize; i++)
        {
            for (int j = -maxSize; j < maxSize; j++)
                for (int k = -maxSize; k < maxSize; j++)
                {
                    vector<_type> t = i * res.at(0) + j * res.at(1) + k * res.at(2);
                    if (-boundX < t[0] < boundX && -boundY < t[1] < boundY)
                    {
                        x.push_back(t.at(0));
                        y.push_back(t.at(1));
                        z.push_back(t.at(2));
                    }
                }
        }
    }
    /*plt::plot({ 0 }, { 0 }, "bo");
    plt::plot({ res.at(0).at(0) }, { res.at(0).at(1) }, "r+");
    plt::plot({ res.at(1).at(0) }, { res.at(1).at(1) }, "r+");
    plt::plot({ og.at(0).at(0) }, { og.at(0).at(1) }, "g*");
    plt::plot({ og.at(1).at(0) }, { og.at(1).at(1) }, "g*");*/
}
//https://tel.archives-ouvertes.fr/tel-01245066v2/document
vector<int> kleinSample(Matrix<int>& B, Matrix<_type>& GSO, float deviation, vector<int>& c)
{
    vector<int> v(c.size(), 0);
    std::random_device mch;
    std::default_random_engine generator(mch());

    double d, z;
    for (int i = B.size() - 1; i >= 0; i--)
    {
        double length = squaredLengthOfVector(GSO.at(i));
        d = dotProduct(c, GSO.at(i)) / length;
        double sigma = deviation / sqrt(length);
        std::normal_distribution<double> distribution(d, 1);
        z = round(distribution(generator));
        //std::cout << z << " " ;
        c = c - z * B[i];
        v = v + z * B[i];
    }
    //std::cout << std::endl;
    return v;
}

// todo discover
vector<_type> sample_gaussian(Matrix<_type>& B)
{
    vector<_type> a;
    a = (5 - rand() % 10) * B.back();
    for (auto b : B)
        a = a + (5 - rand() % 10) * b;
    return a;
}

bool reduce(vector<int>& p, Matrix<int>& L, vector<double>& lengths)
{
    double s = squaredLengthOfVector(p);
    for (size_t i = 0; i < L.size(); i++)
    {
        if (s >= lengths[i])
        {
            auto diff = p - L[i];
            double s2 = squaredLengthOfVector(diff);
            if (s2 < s)
            {
                p = diff;
                return true;
            }
        }
    }
    return false;
}

vector<int> gauss_reduce(vector<int>& v_new, Matrix<int>& L, Matrix<int>& S, vector<double>& lengths)
{
    bool cont = true;
    while (cont)
    {
        auto same = std::find_if(L.begin(), L.end(), [&v_new](vector<int>& v) { return v_new == v; });
        if (same != L.end())
        {
            return vector<int>(v_new.size());
        }
        // 2 раза считается разница, можно сохранить и сравнить с минимумом
        cont = reduce(v_new, L, lengths);
    }
    auto it = L.begin();
    auto it2 = lengths.begin();
    int i = 0;
    vector<int> a(v_new.size());
    while (it != L.end())
    {
        double v_new_length = squaredLengthOfVector(&v_new);
        vector<int> v = *it;
        if (lengths[i] > v_new_length)
        {
            a = v - v_new;
            if (squaredLengthOfVector(a) <= lengths[i])
            {
                S.push_back(a);
                it = L.erase(it);
                it2 = lengths.erase(it2);
            }
            else
            {
                ++it;
                ++it2, ++i;
            }
        }
        else
        {
            ++it;
            ++it2, ++i;
        }
    }
    return v_new;
}
// https://cseweb.ucsd.edu/~daniele/papers/Sieve.pdf
// c= 1000
vector<int> gaussSieve(Matrix<int>& B, float bound, ostream& o)
{
    Matrix<int> L, S;
    vector<int> v_new;
    Matrix<_type> gramSmidt = gSchmidt(B);
    vector<double> squarelengthsVector;
    int K = 0, iter = 0;
    //int c = KabatianskyLevenshteinC(B);
    map<int, int> _map;

    int c = 1000;
    do
    {
        _map.insert(make_pair(iter, K));

        iter++;
        if (S.size() > 0)
        {
            v_new = S.back();
            S.pop_back();
        }
        else
        {
            vector<int> a(B[0].size(), 0);
            v_new = kleinSample(B, gramSmidt, 1, a);
        }
        v_new = gauss_reduce(v_new, L, S, squarelengthsVector);
        if (!v_new.size() ||
            std::count(v_new.begin(), v_new.end(), 0) == v_new.size())
        {
            K++;
        }
        else if (sqrt(squaredLengthOfVector(v_new)) < bound)
        {
            o << "early exit gauss iter : " << iter << "|\n";
            return v_new;
        }
        else
        {
            L.push_back(v_new);
            squarelengthsVector.push_back(squaredLengthOfVector(v_new));
        }
    } while (K < c || L.size() == 0);
    vector<int> shortestVector = L[0];
    int _a = 0, count = 0;
    for (int i = 1; i < L.size(); i++)
    {
        if (sqrt(squaredLengthOfVector(L[i])) <= bound)
        {
            count++;
        }
        if (squaredLengthOfVector(shortestVector) > squaredLengthOfVector(L[i]))
        {
            shortestVector = L[i];
            _a = i;
        }
    }
    cout << "count " << count << endl;
    cout << "gauss at : " << _map.at(_a) << endl;
    o << "gauss iter : " << iter << "|\n";
    return shortestVector;
}

vector<_type> randomInBall(const int dimm, float r)
{
    vector<_type> randVector(dimm);
    for (int i = 0; i < dimm; i++)
    {
        randVector[i] = (rand() % (int)(r + 1));
    }
    return (r * (rand() / double(RAND_MAX)) / sqrt(squaredLengthOfVector(randVector))) * randVector;
}
//https://crypto.stackexchange.com/questions/29661/how-to-find-the-value-of-a-vector-modulo-a-basis-in-lattice-based-cryptography/29701#29701
vector<_type> mod(vector<_type> vec, Matrix<_type> m)
{
    Matrix<_type> r(m.size(), vector<_type>(m.size(), 0));
    inverse(m, r);
    return m * roundV(r * vec);
}

void ListReduce(vector<int>& p, Matrix<int>& L, vector<double>& lengths)
{
    bool cont = true;
    while (cont)
    {
        auto same = std::find(L.begin(), L.end(), p);
        if (same != L.end())
        {
            p = vector<int>(p.size());
            return;
        }
        if (!reduce(p, L, lengths))
            return;
        /*auto f = std::find_if(L.begin(), L.end(), [&p](vector<int>& v) {
            return squaredLengthOfVector(&p) >= squaredLengthOfVector(&v) && squaredLengthOfVector(&(p - v)) < squaredLengthOfVector(&p);
            });
        if (f == L.end()) return;
        p = p - *f;*/
    }
}
//https://math.stackexchange.com/questions/396382/what-does-this-dollar-sign-over-arrow-in-function-mapping-mean
//https://eprint.iacr.org/2014/714.pdf
vector<int> ListSieve(Matrix<int>& B, float bound, ostream& o)
{
    Matrix<int> L;
    Matrix<_type> gramSmidt = gSchmidt(B);
    vector<int> a(B[0].size(), 0);
    vector<double> lengths;
    map<int, int> _map;
    int K = 1000;
    int iter = 0;
    for (int i = 0, k = 0; i < K || L.size() == 0; k++)
    {
        iter++;
        vector<int> v = kleinSample(B, gramSmidt, 1, a);
        ListReduce(v, L, lengths);
        bool isZero = true;
        for (auto in : v)
        {
            if (in != 0)
            {
                isZero = false;
                break;
            }
        }
        _map.insert(make_pair(k, i));
        if (isZero)
        {
            i++;
        }
        else if (sqrt(squaredLengthOfVector(v)) < bound)
        {
            o << "early list iter : " << iter << "|\n";
            return v;
        }
        else
        {
            L.push_back(v);
            lengths.push_back(squaredLengthOfVector(v));
        }
    }
    vector<int> shortestVector = L[0];
    int _a = 0, count = 0;
    for (int i = 1; i < L.size(); i++)
    {
        if (sqrt(squaredLengthOfVector(L[i])) <= bound)
        {
            count++;
        }
        if (squaredLengthOfVector(shortestVector) > squaredLengthOfVector(L[i]))
        {
            shortestVector = L[i];
            _a = i;
        }
    }
    cout << "count " << count << endl;
    cout << "list at : " << _map.at(_a) << endl;
    o << "list iter : " << iter << "|\n";
    return shortestVector;
}

vector<int> sumVectorsWithCoeff(Matrix<int>& matrix, vector<int>& c)
{
    vector<int> a(c.size(), 0);
    for (size_t i = 0; i < c.size(); i++)
    {
        a = a + c[i] * matrix[i];
    }
    return a;
}

bool checkifBiggerByEachCoordinate(vector<int>& a, vector<int>& b)
{
    for (size_t i = 0; i < a.size(); i++)
    {
        if (abs(a[i]) > abs(b[i]))
            return false;
    }
    return true;
}

vector<int> LLLenumerate(Matrix<int>& matrix, ofstream& myfile)
{
    Matrix<int> enumeration;
    int n = matrix.size();
    int lim = pow(2, n * (n - 1) / 4.0);
    vector<int> c(n, 0);
    c[0] = 1;
    vector<int> shortest = sumVectorsWithCoeff(matrix, c);
    auto shortestL = squaredLengthOfVector(shortest);
    long long tmpL;
    vector<int> tmp;
    tmp = sumVectorsWithCoeff(matrix, c);
    if (!checkifBiggerByEachCoordinate(shortest, tmp))
    {
        tmpL = squaredLengthOfVector(tmp);
        if (tmpL > 0 && shortestL > tmpL)
        {
            shortestL = tmpL;
            shortest = tmp;
        }
    }
    for (int i = 0; i < n;)
    {
        int& a = c[i];
        if (a == -lim)
        {
            for (int j = 0; j < i; j++)
                c[j] = -lim;
            a++;
            i = 0;
            tmp = sumVectorsWithCoeff(matrix, c);
            if (!checkifBiggerByEachCoordinate(shortest, tmp))
            {
                tmpL = squaredLengthOfVector(tmp);
                if (tmpL > 0 && shortestL > tmpL)
                {
                    shortestL = tmpL;
                    shortest = tmp;
                }
            }
            // for (auto &bb : c)
            //     cout << bb;
            // cout << endl;
            // // enumeration.push_back(sumVectorsWithCoeff(matrix, c));
        }
        else if (a < lim)
        {
            a++;
            if (a == lim)
                i++;
            else
                i = 0;
            tmp = sumVectorsWithCoeff(matrix, c);
            if (!checkifBiggerByEachCoordinate(shortest, tmp))
            {
                tmpL = squaredLengthOfVector(tmp);
                if (tmpL > 0 && shortestL > tmpL)
                {
                    shortestL = tmpL;
                    shortest = tmp;
                }
            }
            // for (auto &bb : c)
            //     cout << bb;
            // cout << endl;
            // // enumeration.push_back(sumVectorsWithCoeff(matrix, c));
        }
        else if (a == lim)
        {
            i++;
        }
    }
    return shortest;
}

//https://web.eecs.umich.edu/~cpeikert/pubs/lattice-survey.pdf
Matrix<int> generateLattice(int dimm, int diff, int iter, bool print, ofstream& file)
{
    Matrix<int> m(dimm, vector<int>(dimm, 0));

    file << "---------------------------------------------\n\n"
        << "dimm: " << dimm << " iter: " << iter << " diff: " << diff << "|\n";
    int summ = 0;
    //srand(time(0));
    file << '[';
    for (auto& v : m)
    {
        file << '[';
        for (auto& el : v)
        {
            float f = rand();
            el = f / RAND_MAX * (float)diff * 20;
            summ += el;
            file << el << ' ';
        }
        file << "]\n";
    }
    file << "]\n";
    file << " summ: " << summ << "\n\n";
    return m;
}
void preprocessVoronoi(Matrix<int>& B)
{
    LLL(B);
}

vector<int> findRelevantHelper(Matrix<int>& B, int p, int n)
{
    vector<int> res(B.size(), 0);
    for (size_t i = 0; i < n; i++)
    {
        if (p & 1)
        {
            res = res - B[i];
        }
        p = p >> 1;
    }
    for (auto& b : res)
        b = b / 2;
    return res;
}

vector<int> CVPP2V(vector<int>& t, Matrix<vector<int>>& Voronoi) {
    int i = 0;
    vector<int> t0(t.size()), ti(t.size());
    while (find_if(Voronoi.begin(), Voronoi.end(), [&ti](auto& v) {
        return v[0] == ti
    }) == Voronoi.end()) {
        float tmp = mu(ti, Voronoi[0][0]);
        int maxi = 0;
        for (size_t j = 1; j < Voronoi.size(); j++)
        {
            if (tmp < mu(ti,Voronoi[j][0])) {
                maxi = j;
            }
        }
        ti = ti - ti;
        i++;
    }

    return t0 - ti;
}
vector<int> CVPP(vector<int>& t, Matrix<int>& B, Matrix<vector<int>>& Voronoi) {
    vector<int> res(t.size(), 0), tp;
    for (auto& v : B) {
        res = res + mu(t, v) * v;
    }
    // not sure about this 
    int p = 0;
    while(true) {
        p++;
        if (squaredLengthOfVector(t) < (squaredLengthOfVector(pow(2, p) * Voronoi[0][0]))) {
            break;
        }
    }

    for (int i = p - 1; i > 0; i--) {
        tp = tp - CVPP2V(tp, pow(2, i - 1) * Voronoi);
    }
    return res - tp;
}

vector<int> rankReduceCVP(const vector<int>& t, Matrix<int>& B, Matrix<vector<int>>& Voronoi, int H) {
    vector<int> res(t.size(), 0), v, closest, tmp;
    int h = mu(t, B.at(Voronoi.size()));
    float length, tmplength;
    int i = h - H + 1;
    res = CVPP(t - h * B.at(Voronoi.size()), Voronoi.back());
    length = squaredLengthOfVector(res);
    i++;
    // todo is it correct?
    for (; abs(i - h) < H; i++) {
        v = CVPP(t - i * B.at(Voronoi.size()), Voronoi.back());
        tmp = v - t;
        tmplength = squaredLengthOfVector(tmp);
        if (tmplength < length) {
            length = tmplength;
            res = tmp;
        }
    }
    return res;
}

Matrix<vector<int>> findRelevant(Matrix<int>& B, Matrix<vector<int>>& Voronoi, int H)
{
    Matrix<vector<int>> res;
    vector<int> t, c;
    int l = pow(2, B.size());
    // p - вектор из нелей и 1
    for (int p = 1; p < l; p++)
    {
        t = findRelevantHelper(B, p, B.size());
        c = rankReduceCVP(t, B, Voronoi, H);
        t = 2 * (t - c);
        c = inverse(t);
        res.push_back({ t, c });
    }
    return res;
}
void removeNonRelevant(Matrix<vector<int>>& res)
{
    for (size_t i = 0; i < res.size(); i++)
    {
        for (size_t j = 0; j < res.size(); j++)
        {
            if (i != j && squaredLengthOfVector(res[j][0] - 0.5 * res[i][0]) <= squaredLengthOfVector(0.5 * res[i][0])) {
                // debug whats in it
                res.erase(res.begin() + i);
            }
        }

    }

}

Matrix<vector<int>> computeVcell(Matrix<int>& B, Matrix<vector<int>>& Voronoi, int H)
{
    Matrix<vector<int>> res = findRelevant(B, Voronoi, H);
    removeNonRelevant(res);
    return res;
}
// https://cseweb.ucsd.edu/~pvoulgar/files/voronoi_full.pdf
void basicVoronoiCell(Matrix<int>& B)
{
    preprocessVoronoi(B);
    Matrix<vector<int>> Voronoi(1);
    Matrix<int> newB;
    Voronoi[0] = { B[0], inverse(B[0]) };
    for (int i = 1; i < B.size(); i++)
    {
        for (size_t j = 0; j < i; j++)
        {
            newB.push_back(B[j]);
        }
        // nedd to remake to vec of vec of vec
        Voronoi = computeVcell(newB, Voronoi, pow(2, 0.5 * i));
    }
}

int main()
{
    ofstream myfile;
    myfile.open("example.txt");
    try
    {
        for (int dimm = 2; dimm < 15; dimm += 2)
        {
            for (int diff = 2; diff <= 10; diff += 2)
            {
                for (int iter = 0; iter < 3; iter++)
                {
                    myfile << dimm << " " << diff << " " << iter << endl;
                    Matrix<int> matrix = generateLattice(dimm, diff, iter, true, myfile);
                    Matrix<int> matrix2 = matrix;
                    Matrix<int> matrix3 = matrix;
                    Matrix<int> original(matrix2);
                    Matrix<float> gramSmidt = gSchmidt(matrix);
                    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> mm(dimm, dimm);
                    for (int i = 0; i < dimm; i++)
                    {
                        for (int j = 0; j < dimm; j++)
                        {
                            mm(i, j) = matrix[i][j];
                        }
                    }

                    cout << "starting lll" << endl;
                    auto t1 = std::chrono::high_resolution_clock::now();
                    LLL(matrix, myfile);
                    vector<int> lllshortest = LLLenumerate(matrix, myfile);
                    auto t2 = std::chrono::high_resolution_clock::now();
                    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
                    myfile << "lll duration " << duration;
                    myfile << "\n\n";

                    long long det = mm.determinant();
                    // double pp = abs(det);
                    double pp = 1;
                    float bound = 1.01 * (sqrt(dimm)) * pow(pp, 1.0 / dimm);
                    myfile << "det: " << (long long)det << "\nbound: " << bound << endl;

                    cout << "starting gauss" << endl;
                    t1 = std::chrono::high_resolution_clock::now();
                    vector<int> sv = gaussSieve(matrix2, bound, myfile);
                    t2 = std::chrono::high_resolution_clock::now();
                    duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
                    myfile << "gauss duration: " << duration;
                    myfile << "\n\n";
                    cout << "f\n\n";

                    cout << "starting list" << endl;
                    t1 = std::chrono::high_resolution_clock::now();
                    vector<int> v3 = ListSieve(matrix3, bound, myfile);
                    t2 = std::chrono::high_resolution_clock::now();
                    duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
                    myfile << "list duration: " << duration;

                    if (!(v3 == sv) || !(v3 == lllshortest))
                    {
                        myfile << "\n\n\n\n !!!!!! error \n\n\n\n";
                        myfile << "gauss" << endl;
                        printVec(sv, myfile);
                        myfile << "\n\n";
                        myfile << "list" << endl;
                        printVec(v3, myfile);
                        myfile << "\n\n";
                        myfile << "lll" << endl;
                        printVec(lllshortest, myfile);
                        myfile << "\n\n";
                        myfile << "gauss length: " << (long)squaredLengthOfVector(sv) << "\n\n list length: " << (long)squaredLengthOfVector(v3)
                            << "\n\n lll length: " << (long)squaredLengthOfVector(lllshortest) << " \n\n\n\n";
                    }
                    else
                    {
                        myfile << "shortest" << endl;
                        printVec(lllshortest, myfile);
                    }
                    myfile << "\n\n";
                }
            }
        }
    }
    catch (exception e)
    {
        cout << e.what();
    }
    myfile.close();
    return 0;
}