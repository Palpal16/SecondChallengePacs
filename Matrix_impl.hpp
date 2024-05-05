#ifndef MATRIX_IMPL_H
#define MATRIX_IMPL_H

#include "Matrix.hpp"

using namespace algebra;

template <std::size_t StorageOrder>
bool
compare<StorageOrder>::operator()(const std::array<std::size_t, 2>& lhs, const std::array<std::size_t, 2>& rhs) const {
    if constexpr (StorageOrder == 0) {
        return lhs[0] < rhs[0] || (lhs[0] == rhs[0] && lhs[1] < rhs[1]);
    } else {
        return lhs[1] < rhs[1] || (lhs[1] == rhs[1] && lhs[0] < rhs[0]);
    }
}


template <typename T, std::size_t StorageOrder>
void
Matrix<T,StorageOrder>::resize(std::size_t R, std::size_t C){

    // If the matrix is compressed, uncompress it
    if(Compressed){
        uncompress();
    }

    // If the new dimensions are smaller, warn the user that some data might be erased
    if(R<rows || C<cols){
        std::cerr << "Warning: Resizing to smaller dimensions might erase some data" << std::endl;
        // Erase the data that is out of the new dimensions
        for(auto it = data.begin(); it != data.end();){
            if(it->first[0] >= R || it->first[1] >= C){
                it = data.erase(it);
            }
            else{
                ++it;
            }
        }
    }
    // Set the new dimensions
    rows = R;
    cols = C;
}



/*  -------------  CALL OPERATORS  -------------  */

// Non-const call operator
template <typename T, std::size_t StorageOrder>
T&
Matrix<T,StorageOrder>::operator()(const std::array<std::size_t, 2>& index){
    // If the matrix is compressed
    if(!Compressed){
        // If the index is outside matrix dimensions, update the dimensions
        if(index[0] >= rows){
            rows = index[0] + 1;
        }
        if(index[1] >= cols){
            cols = index[1] + 1;
        }
        // Return the value of the index
        return data[index];
    }

    // If the matrix is compressed
    // If the index is not found in the compressed data, return an error
    if(index[0] >= Inner_dir.size() - 1 || index[1] >= Inner_dir.size() - 1){
        throw std::out_of_range("Index outside the compressed matrix");
    }
    
    // Determine the row and column based on the storage order
    std::size_t row = index[StorageOrder == 0 ? 0 : 1];
    std::size_t col = index[StorageOrder == 0 ? 1 : 0];
    // Search for the value in the compressed data
    for (std::size_t j = Inner_dir[row]; j < Inner_dir[row + 1]; ++j) {
        if (Outer_dir[j] == col) {
            return Values[j];
        }
    }

    // If the index is not found in the compressed data, return an error
    throw std::out_of_range("Index not found in the compressed matrix");
}

// Const call operator
template <typename T, std::size_t StorageOrder>
T
Matrix<T,StorageOrder>::operator()(const std::array<std::size_t, 2>& index) const {
    // If the index is outside matrix dimensions, return an error
    if (index[0] >= rows || index[1] >= cols) {
        throw std::out_of_range("Index outside the sparse matrix");
    }

    // If the matrix is not compressed
    if (!Compressed) {
        auto it = data.find(index);
        // If the index is not found, return a default value
        if (it == data.end()) {
            return T();
        }
        // If the index is found, return the value
        return it->second;
    }
    
    // If the matrix is compressed
    // Call for compressed matrix
    std::size_t row = index[StorageOrder == 0 ? 0 : 1];
    std::size_t col = index[StorageOrder == 0 ? 1 : 0];
    for (std::size_t j = Inner_dir[row]; j < Inner_dir[row + 1]; ++j) {
        if (Outer_dir[j] == col) {
            return Values[j];
        }
    }

    // If the index is not found in the compressed data (returns 0)
    return T();
}



// Call operator that returns a row or a column based on the StorageOrder
template <typename T, std::size_t StorageOrder>
std::vector<T>
Matrix<T,StorageOrder>::operator()(std::size_t index) {
    // If the matrix is not compressed
    if(!Compressed){
        // If the storage order is row-major
        if constexpr (StorageOrder == 0) {
            // Create a vector to store the result
            std::vector<T> result(cols);
            for(auto it =data.lower_bound({index, 0}); it!=data.lower_bound({index + 1, 0}); ++it){
                result[it->first[1]] = it->second;
            }
            return result;
        } else {
            // If the storage order is column-major
            std::vector<T> result(rows);
            for(auto it =data.lower_bound({0, index}); it!=data.lower_bound({0, index + 1}); ++it){
                result[it->first[0]] = it->second;
            }
            return result;
        }
    }

    // If the matrix is compressed
    if constexpr (StorageOrder == 0) {
        std::vector<T> result(cols);
        for (std::size_t j = Inner_dir[index]; j < Inner_dir[index + 1]; ++j) {
            result[Outer_dir[j]] = Values[j];
        }
        return result;
    } else {
        std::vector<T> result(rows);
        for (std::size_t i = Inner_dir[index]; i < Inner_dir[index + 1]; ++i) {
            result[Outer_dir[i]] = Values[i];
        }
        return result;
    }
}


/*  -------------  COMPRESS AND UNCOMPRESS  -------------  */

// define the methods to compress the matrix using CSR or CSC based on StorageOrder
template <typename T, std::size_t StorageOrder>
void
Matrix<T,StorageOrder>::compress() {
    // If the matrix is already compressed, do nothing
    if (Compressed) {
        return;
    }

    // Initialize the current row
    std::size_t current_row = 0;
    // Add the first element to the inner directory
    Inner_dir.emplace_back(0);

    for (const auto& [key, value] : data) {
        // If the key is greater than the current row
        // Then the row is empty, so we add the size, to create an empty row in compressed form
        while (key[StorageOrder] > current_row) {
            Inner_dir.emplace_back(Values.size());
            ++current_row;
        }
        // Add the value to the values and the other key to the outer directory
        Values.emplace_back(value);
        Outer_dir.emplace_back(key[1 - StorageOrder]);
    }

    // Add last element or all the elements remaining to get to the right size
    while(Inner_dir.size() - 1 < (StorageOrder == 0 ? rows : cols)){
        Inner_dir.emplace_back(Values.size());
    }

    // Set the matrix as compressed
    Compressed = true;

    // Clear the uncompressed data
    data.clear();
}

// Method to uncompress the matrix
template <typename T, std::size_t StorageOrder>
void
Matrix<T,StorageOrder>::uncompress() {
    // If the matrix is not compressed, do nothing
    if (!Compressed) {
        return;
    }

    for (std::size_t i = 0; i < Inner_dir.size() - 1; ++i) {
        // For each element in the range defined by the inner directory
        for (std::size_t j = Inner_dir[i]; j < Inner_dir[i + 1]; ++j) {
            // Create the key based on the storage order
            std::array<std::size_t, 2> key = {StorageOrder == 0 ? i : Outer_dir[j], StorageOrder == 0 ? Outer_dir[j] : i};
            // Add the key and value to the data
            data.emplace(key, Values[j]);
        }
    }

    // Clear the compressed data
    Values.clear();
    Inner_dir.clear();
    Outer_dir.clear();

    // Set the matrix as uncompressed
    Compressed = false;
}



/*  -------------  READ AND PRINT  -------------  */

// Method to read from a file in matrix market format
template <typename T, std::size_t StorageOrder>
void
Matrix<T,StorageOrder>::read (const std::string &filename){
    // Open the file
    std::ifstream ifs(filename);
    // If the file cannot be opened, print a warning and return
    if(!ifs){
        std::cerr << "Warning: Impossibile aprire il file " << filename << std::endl;
        return;
    }

    // Read the first line that is not a comment or empty
    std::string line;
    do {
        std::getline(ifs, line);
    } while (line[0] == '%' || line.empty());
    
    // Parse the first line to get the number of rows, columns, and values
    std::istringstream first_line (line);
    std::size_t num_val;
    first_line >> rows >> cols >> num_val;
    // Clear the data and set the matrix as uncompressed
    data.clear();
    Compressed=false;

    // Read the rest of the lines
    while(std::getline(ifs, line)){
        // Parse each line to get the row, column, and value
        std::istringstream current_line (line);
        std::size_t i, j;
        T val;
        current_line >> i >> j >> val;
        // Add the value to the data
        data[{i-1,j-1}] = val;
    }

    // If the number of values read does not match the number of values in the header, print a warning
    if(num_val != data.size()){
        std::cerr << "Warning: Number of values read does not match the number of values in the header" << std::endl;
    }
    // Close the file
    ifs.close();
}


template <typename T, std::size_t StorageOrder>
void
Matrix<T,StorageOrder>::print() const {
    // Print the compression status and dimensions of the matrix
    std::cout << (Compressed==1 ? "Compressed" : "Uncompressed" ) << " matrix with dimensions: " << rows << "x" << cols << std::endl;
    if(Compressed){
        // If compressed simply print each value with its row and column (by parsing with the defined ordering)
        for (std::size_t i = 0; i < Inner_dir.size() - 1; ++i) {
            for (std::size_t j = Inner_dir[i]; j < Inner_dir[i + 1]; ++j) {
                std::cout << (StorageOrder == 0 ? i : Outer_dir[j]) << " " << (StorageOrder == 0 ? Outer_dir[j] : i) << " : " << Values[j] << std::endl;
            }
        }
        return;
    }
    // If the matrix is not compressed print each value in the data
    for (const auto& [key, value] : data) {
        std::cout << key[0] << " " << key[1] << " : " << value << std::endl;
    }
}



/*  -------------  PRODUCTS  -------------  */

namespace algebra {
    // Define the operator* to multiply a matrix by a vector
    template <typename T, std::size_t StorageOrder>
    std::vector<T>
    operator*(const Matrix<T, StorageOrder>& lhs, const std::vector<T>& rhs) {
        // If the number of columns of the matrix does not match the size of the vector, throw an exception
        if(lhs.cols != rhs.size()){
            throw std::invalid_argument("Matrix and vector dimensions do not match");
        }

        // Initialize the result vector with the number of rows of the matrix
        std::vector<T> result(lhs.rows);

        if( ! lhs.Compressed){
            // If not compressed, simply iterate over the data and multiply the values
            for (const auto& [key, value] : lhs.data) {
                result[key[0]] += value * rhs[key[1]];
            }
            return result;
        }

        // If compressed, iterate over the compressed data and multiply the values
        if constexpr (StorageOrder == 0){
            for (std::size_t i = 0; i < lhs.Inner_dir.size()-1; ++i) {
                for (std::size_t j = lhs.Inner_dir[i]; j < lhs.Inner_dir[i + 1]; ++j) {
                    result[i] += lhs.Values[j] * rhs[lhs.Outer_dir[j]];
                }
            }
            return result;
        }

        for (std::size_t j = 0; j < lhs.Inner_dir.size() -1; ++j) {
            for (std::size_t i = lhs.Inner_dir[j]; i < lhs.Inner_dir[j + 1]; ++i) {
                result[lhs.Outer_dir[i]] += lhs.Values[i] * rhs[j];
            }
        }
        return result;
    }


    // Define the operator* to multiply a matrix by a matrix
    // Implemented only when the return matrix is row-wise and the multiplied matrices are both compressed or both not
    template <typename T, std::size_t SO1, std::size_t SO2>
    Matrix<T,0>
    operator*(const Matrix<T, SO1>& lhs, const Matrix<T, SO2>& rhs) {
        // If the number of columns of the first matrix does not match the number of rows of the second matrix, throw an exception
        if(lhs.cols != rhs.rows){
            throw std::invalid_argument("Matrix dimensions do not match");
        }

        // Initialize the result matrix with the right dimensions
        Matrix<T,0> result(lhs.rows, rhs.cols);

        // If the second matrix is row-wise, then I go through the elements of the first matrix
        // and for each element I look for the corresponding row in the second matrix
        if constexpr (SO2 == 0){
            // Both matrices are not compressed
            if(! lhs.Compressed && ! rhs.Compressed){
                for (const auto& [key, value] : lhs.data) {
                    for(auto it =rhs.data.lower_bound({key[1], 0}); it!=rhs.data.lower_bound({key[1] + 1, 0}); ++it){
                        result(key[0], it->first[1]) += value * it->second;
                    }
                }
                return result;
            }
            // Both matrices are compressed
            if(lhs.Compressed && rhs.Compressed){
                for (std::size_t i = 0; i < lhs.Inner_dir.size()-1; ++i) {
                    for (std::size_t j = lhs.Inner_dir[i]; j < lhs.Inner_dir[i + 1]; ++j) {
                        // We change the way to find the row and column indexes of the first matrix based on the StorageOrder
                        for(std::size_t k = rhs.Inner_dir[SO1 == 0 ? lhs.Outer_dir[j] : i]; k < rhs.Inner_dir[(SO1 == 0 ? lhs.Outer_dir[j] : i) + 1]; ++k){
                            result(SO1 == 0 ? i : lhs.Outer_dir[j], rhs.Outer_dir[k]) += lhs.Values[j] * rhs.Values[k];
                        }
                    }
                }
                result.compress();
                return result;
            }
        }
        else{
            // If the second matrix is column-wise
            if(! lhs.Compressed && ! rhs.Compressed){
                for (const auto& [key, value] : lhs.data) {
                    for(const auto& [key2, value2] : rhs.data){
                        if(key[1] == key2[0]){
                            result(key[0], key2[1]) += value * value2;
                        }
                    }
                }
                return result;
            }
            if(lhs.Compressed && rhs.Compressed){
                for (std::size_t j=0; j<rhs.Inner_dir.size() -1; ++j) {
                    for (std::size_t i = rhs.Inner_dir[j]; i < rhs.Inner_dir[j + 1]; ++i) {
                        if constexpr (SO1 == 1){
                            for(std::size_t k = lhs.Inner_dir[rhs.Outer_dir[i]]; k < lhs.Inner_dir[rhs.Outer_dir[i] + 1]; ++k){
                                result(lhs.Outer_dir[k], j) += lhs.Values[k] * rhs.Values[i];
                            }
                        }
                        else{
                            for(std::size_t k = 0; k < lhs.Inner_dir.size()-1; ++k){
                                T left = lhs(k, rhs.Outer_dir[i]);
                                if(left != 0){
                                    result(k, j) += left * rhs.Values[i];
                                }
                            }
                        }
                    }
                }
                result.compress();
                return result;
            }
        }

        // If the method is not implemented, print a warning
        std::cerr << "Warning: Not implemented yet" << std::endl;
        return result;        
    }

}

#endif // MATRIX_IMPL_H