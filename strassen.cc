#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <chrono>
#include <stdexcept>

using std::vector;
using namespace std;

vector< vector<vector<int> > > read_matrices_from_file(const char* file_path, int dimension) {
    ifstream file(file_path);
    vector<vector<int> > matrix_a(dimension, vector<int>(dimension));
    vector<vector<int> > matrix_b(dimension, vector<int>(dimension));

    for (int i = 0; i < 2 * dimension * dimension; i++) {
        int x;
        file >> x;
        if (i < dimension * dimension) {
            matrix_a[i / dimension][i % dimension] = x;
        } else {
            matrix_b[(i - dimension * dimension) / dimension][(i - dimension * dimension) % dimension] = x;
        }
    }

    vector<vector<vector<int> > >result;
    result.push_back(matrix_a);
    result.push_back(matrix_b);
    return result;
}

vector<vector<int> > addition(const vector<vector<int> >& matrix_1, const vector<vector<int> >& matrix_2) {
    int n = matrix_1.size();
    vector<vector<int> > result(n, vector<int>(n));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            result[i][j] = matrix_1[i][j] + matrix_2[i][j];
        }
    }

    return result;
}

vector<vector<int> > subtraction(const vector<vector<int> >& matrix_1, const vector<vector<int> >& matrix_2) {
    int n = matrix_1.size();
    vector<vector<int> > result(n, vector<int>(n));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            result[i][j] = matrix_1[i][j] - matrix_2[i][j];
        }
    }

    return result;
}

vector<vector<int> > trad_multiplication(const vector<vector<int> >& matrix_1, const vector<vector<int> >& matrix_2) {
    int n = matrix_1.size();
    vector<vector<int> > result(n, vector<int>(n));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            result[i][j] = 0;
            for (int k = 0; k < n; k++) {
                result[i][j] += matrix_1[i][k] * matrix_2[k][j];
            }
        }
    }

    return result;
}

vector<vector<int>> optimized_multiplication(const vector<vector<int>>& matrix_1, const vector<vector<int>>& matrix_2) {
    int n1 = matrix_1.size();
    int m1 = matrix_1[0].size();
    int n2 = matrix_2.size();
    int m2 = matrix_2[0].size();

    if (m1 != n2) {
        throw std::invalid_argument("The dimensions of the input matrices are not compatible for multiplication.");
    }

    vector<vector<int>> result(n1, vector<int>(m2, 0));

    const int blockSize = 32; // Adjust this value based on your system's cache size

    for (int ii = 0; ii < n1; ii += blockSize) {
        for (int jj = 0; jj < m2; jj += blockSize) {
            for (int kk = 0; kk < m1; kk += blockSize) {
                for (int i = ii; i < std::min(ii + blockSize, n1); i++) {
                    for (int k = kk; k < std::min(kk + blockSize, m1); k++) {
                        for (int j = jj; j < std::min(jj + blockSize, m2); j++) {
                            result[i][j] += matrix_1[i][k] * matrix_2[k][j];
                        }
                    }
                }
            }
        }
    }

    return result;
}
vector<vector<int> > strassens_helper(const vector<vector<int> >& A, const vector<vector<int> >& B, int n0) {
    int n = A.size();

    if (n <= n0) {
        return trad_multiplication(A, B);
    }

    int mid = n / 2;
    vector<vector<int> > a(mid, vector<int>(mid));
    vector<vector<int> > b(mid, vector<int>(n - mid));
    vector<vector<int> > c(n - mid, vector<int>(mid));
    vector<vector<int> > d(n - mid, vector<int>(n - mid));
    vector<vector<int> > e(mid, vector<int>(mid));
    vector<vector<int> > f(mid, vector<int>(n - mid));
    vector<vector<int> > g(n - mid, vector<int>(mid));
    vector<vector<int> > h(n - mid, vector<int>(n - mid));

    for (int i = 0; i < mid; i++) {
        for (int j = 0; j < mid; j++) {
            a[i][j] = A[i][j];
            e[i][j] = B[i][j];
        }
    }
for (int i = 0; i < mid; i++) {
    for (int j = mid; j < n; j++) {
        b[i][j - mid] = A[i][j];
        f[i][j - mid] = B[i][j];
    }
}

for (int i = mid; i < n; i++) {
    for (int j = 0; j < mid; j++) {
        c[i - mid][j] = A[i][j];
        g[i - mid][j] = B[i][j];
    }
}

for (int i = mid; i < n; i++) {
    for (int j = mid; j < n; j++) {
        d[i - mid][j - mid] = A[i][j];
        h[i - mid][j - mid] = B[i][j];
    }
}

vector<vector<int> > p1 = strassens_helper(a, subtraction(f, h), n0);
vector<vector<int> > p2 = strassens_helper(addition(a, b), h, n0);
vector<vector<int> > p3 = strassens_helper(addition(c, d), e, n0);
vector<vector<int> > p4 = strassens_helper(d, subtraction(g, e), n0);
vector<vector<int> > p5 = strassens_helper(addition(a, d), addition(e, h), n0);
vector<vector<int> > p6 = strassens_helper(subtraction(b, d), addition(g, h), n0);
vector<vector<int> > p7 = strassens_helper(subtraction(c, a), addition(e, f), n0);

vector<vector<int> > C11 = addition(subtraction(addition(p4, p5), p2), p6);
vector<vector<int> > C12 = addition(p1, p2);
vector<vector<int> > C21 = addition(p3, p4);
vector<vector<int> > C22 = addition(subtraction(addition(p1, p5), p3), p7);

vector<vector<int> > C(n, vector<int>(n));
for (int i = 0; i < mid; i++) {
    for (int j = 0; j < mid; j++) {
        C[i][j] = C11[i][j];
    }
}

for (int i = 0; i < mid; i++) {
    for (int j = mid; j < n; j++) {
        C[i][j] = C12[i][j - mid];
    }
}

for (int i = mid; i < n; i++) {
    for (int j = 0; j < mid; j++) {
        C[i][j] = C21[i - mid][j];
    }
}

for (int i = mid; i < n; i++) {
    for (int j = mid; j < n; j++) {
        C[i][j] = C22[i - mid][j - mid];
    }
}

return C;
}

vector<vector<int> > strassens(const vector<vector<int> >& A, const vector<vector<int> >& B, int n0) {
int n = A.size();

if (n != A[0].size() || n != B.size() || n != B[0].size()) {
    return {};
}

int m = pow(2, ceil(log2(n)));
vector<vector<int> > A_padded(m, vector<int>(m));
vector<vector<int> > B_padded(m, vector<int>(m));
for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
        A_padded[i][j] = A[i][j];
        B_padded[i][j] = B[i][j];
    }
}

vector<vector<int> > C_padded = strassens_helper(A_padded, B_padded, n0);
vector<vector<int> > C(n, vector<int>(n));

for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
        C[i][j] = C_padded[i][j];
    }
}

return C;
}

void print_to_stdout(const char* str) {
cout << str;
}

int main(int argc, char* argv[]) {
if (argc != 4) {
print_to_stdout("Usage: ./strassen flag dimension inputfile");
return 1;
}
int flag = atoi(argv[1]);
int dimension = atoi(argv[2]);
const char* input_file = argv[3];

//auto matrix_a = read_matrices_from_file(input_file, dimension)[0];
//auto matrix_b = read_matrices_from_file(input_file, dimension)[1];
dimension = 63;
std::vector<int> values = {-1,0, 1};

    // Seed the random number generator
    std::srand(static_cast<unsigned>(std::time(0)));
    
    //Create a 2D vector of size (rows, cols) with random elements from the 'values' vector
    std::vector<std::vector<int>> matrix_a(dimension, std::vector<int>(dimension));
    std::vector<std::vector<int>> matrix_b(dimension, std::vector<int>(dimension));

    for (int i = 0; i < dimension; ++i) {
        for (int j = 0; j < dimension; ++j) {
            int random_index = std::rand() % values.size();
            matrix_a[i][j] = values[random_index];
        }
    }
    for (int i = 0; i < dimension; ++i) {
        for (int j = 0; j < dimension; ++j) {
            int random_index = std::rand() % values.size();
            matrix_b[i][j] = values[random_index];
        }
    }

// for (int i = 1; i < dimension; ++i){
// auto start1 = std::chrono::high_resolution_clock::now();
// auto results1 = strassens(matrix_a, matrix_b, i);
// auto end1 = std::chrono::high_resolution_clock::now();
// auto duration1 = std::chrono::duration_cast<std::chrono::nanoseconds>(end1 - start1).count();
// std::cout << "Elapsed time: " << duration1 << " nanoseconds" << std::endl;
// auto start2 = std::chrono::high_resolution_clock::now();
// auto results2 = trad_multiplication(matrix_a, matrix_b);
// auto end2 = std::chrono::high_resolution_clock::now();
// auto duration2 = std::chrono::duration_cast<std::chrono::nanoseconds>(end2 - start2).count();
// std::cout << "Elapsed time trad: " << duration2 << " nanoseconds" << std::endl;
//}
auto results = strassens(matrix_a, matrix_b, 1);
for (int i = 0; i < dimension; i++){
    std::cout << (results[i][i]) << std::endl;
}


return 0;
}