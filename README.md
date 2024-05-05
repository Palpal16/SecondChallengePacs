
# Sparse matrix Class

This is a C++ implementation of a Matrix class that supports both row-major and column-major storage orders. The class is defined in the `algebra` namespace.

## Features

- The class supports both row-major (StorageOrder = 0) and column-major (StorageOrder = 1) storage orders.
- The class provides methods to resize the matrix, access elements, compress and uncompress the matrix, read the matrix from a file, and print the matrix.
- The class provides getters for the dimensions of the matrix.
- The class provides a method to check if the matrix is compressed.
- The class overloads the multiplication operator to support multiplication with a vector and with another matrix.

## Usage

```
cpp
#include "Matrix.hpp"

// Create a matrix
algebra::Matrix<int, 0> matrix(5, 5);

// Resize the matrix
matrix.resize(10, 10);

// Access elements
int value = matrix(2, 3);

// Compress and uncompress the matrix
matrix.compress();
matrix.uncompress();

// Read the matrix from a file
matrix.read("matrix.txt");

// Print the matrix
matrix.print();

// Get the dimensions of the matrix
std::size_t rows = matrix.getRows();
std::size_t cols = matrix.getCols();

// Check if the matrix is compressed
bool isCompressed = matrix.is_compressed();

// Multiply the matrix with a vector
std::vector<int> vec(10, 1);
std::vector<int> result = matrix * vec;

// Multiply the matrix with another matrix
algebra::Matrix<int, 0> otherMatrix(10, 10);
algebra::Matrix<int, 0> resultMatrix = matrix * otherMatrix;

```

Note
The multiplication operation is performed differently based on whether the matrices are compressed or not and the storage order. If the dimensions do not match, an exception is thrown. If the matrices are not in the same both compressed or not, a warning is printed.