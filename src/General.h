#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <stdio.h>
#include <fstream>
#include <sstream>

typedef std::vector<double> vecDouble;
typedef std::vector<int> vecInt;
typedef std::vector<vecDouble> matDouble;

using std::cout;
using std::endl;

/// Modifies an argument with the product between two matrices (2x2)
/// @param matDouble mat1 @param matDouble mat2 @param matDouble result = mat1 * mat2
void productMat2x2(const matDouble &mat1, const matDouble &mat2, matDouble &result);

/// Modifies an argument with the multiplication of a matrices (2x2) by a scalar
/// @param matDouble mat1 @param double scalar @param matDouble result = scalar * mat1
void scalarMat2x2(const matDouble &mat1, const double &scalar, matDouble &result);

/// Modifies an argument with the transpose matrix
/// @param matDouble matrix @param matDouble transpose (transpose[i][j] = matrix[j][i])
void transposeMat2x2(const matDouble &matrix, matDouble &transpose);

/// Modifies an argument with the contraction between two matrices
/// @param matDouble mat1 @param matDouble mat2 @param matDouble result = mat1 : mat2
void contractionMat2x2(const matDouble &mat1, const matDouble &mat2, double &result);

/// Sets zero to a matrix
/// @param matDouble mat (mat[i][j] = 0.0)
void setZero2x2MatDouble(matDouble &mat);

/// Modifies an argument with the inverse of a matrix and the determinant
/// @param matDouble matrix @param matDouble inverseMatrix @param double det
void computeInverseAndDeterminant2x2(const matDouble &matrix, matDouble &inverseMatrix, double &det);

/// Sets a number of columns for a matrix with two lines
/// @param matDouble matrix @param int number of columns
void setColumnToMatDoubleWith2Lines(matDouble &mat, const int &nColumn);

/// Sets a number of columns for a matrix with three lines
/// @param matDouble matrix @param int number of columns
void setColumnToMatDoubleWith3Lines(matDouble &mat, const int &nColumn);

/// Sets a size for a matrix
/// @param matDouble matrix @param int number of lines @param int number of columns @param double value (optional)
void setSizeMatDouble(matDouble &mat, const int &nLine, const int &nColumn, const double &value = 0.0);

std::vector<std::string> split(std::string str, std::string delim);

/// Modifies an argument with the sum of two matrices
/// @param matDouble mat1 @param matDouble mat2 @param matDouble result = mat1 + mat2
void sumMat2x2(const matDouble &mat1, const matDouble &mat2, matDouble &result);


// void product(const matDouble &mat, const vecDouble &vec, vecDouble &result); //A*b

// void product(const matDouble &mat1, const matDouble &mat2, matDouble &result); //A*A

// void innerproduct(const vecDouble &vec1, const vecDouble &vec2, double &result); //b*b

// void contraction(const matDouble &mat1, const matDouble &mat2, double &result);

void setValueMatDouble(matDouble &mat, const double &value = 0.0);

void setValueVecDouble(vecDouble &mat, const double &value = 0.0);

// void computeInverseAndDeterminant(const matDouble &matrix, matDouble &inverseMatrix, double &det);

// void computeInverse(const matDouble &matrix, matDouble &inverseMatrix);

// void transposeMat(const matDouble &matrix, matDouble &transpose);

// void addition(const matDouble &mat1, const matDouble &mat2, matDouble &result);

// void setZero3x3MatDouble(matDouble &mat);

// matDouble getProduct(const matDouble &mat1, const matDouble &mat2);

// matDouble getAddition(const matDouble &mat1, const matDouble &mat2);

/// Returns the transpose matrix
/// @param matDouble matrix
/// @return transpose of matrix
matDouble getTransposeMat2x2(const matDouble &matrix);

/// Returns the product between two matrices
/// @param matDouble mat1 @param matDouble mat2
/// @return mat1 * mat2
matDouble getProductMat2x2(const matDouble &mat1, const matDouble &mat2);

/// Returns the sum of two matrices
/// @param matDouble mat1 @param matDouble mat2
/// @return mat1 + mat2
matDouble getSumMat2x2(const matDouble &mat1, const matDouble &mat2);

/// Returns the contraction between two matrices
/// @param matDouble mat1 @param matDouble mat2
/// @return mat1 : mat2
double getContractionMat2x2(const matDouble &mat1, const matDouble &mat2);


void jumpLine(const int &nLines, std::ifstream &file);


/// Modifies an argument with the subtraction of two matrices
/// @param matDouble mat1 @param matDouble mat2 @param matDouble result = mat1 - mat2
void subtractionMat3x3(const matDouble &mat1, const matDouble &mat2, matDouble &result);

/// Returns the sum of two matrices
/// @param matDouble mat1 @param matDouble mat2
/// @return mat1 + mat2
void sumMat3x3(const matDouble &mat1, const matDouble &mat2, matDouble &result);

/// Returns the sum of two matrices
/// @param matDouble mat1 @param matDouble mat2
/// @return mat1 + mat2
matDouble getSumMat3x3(const matDouble &mat1, const matDouble &mat2);

/// Returns the multiplication of a matrices (3x3) by a scalar
/// @param matDouble mat1 @param double scalar
matDouble getScalarMat3x3(const matDouble &mat1, const double &scalar);

/// Returns the contraction between two matrices
/// @param matDouble mat1 @param matDouble mat2
/// @return mat1 : mat2
double getContractionMat3x3(const matDouble &mat1, const matDouble &mat2);