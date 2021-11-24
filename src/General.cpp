#include "General.h"

void productMat2x2(const matDouble &mat1, const matDouble &mat2, matDouble &result) //matriz dois por dois
{
    result[0][0] = mat1[0][0] * mat2[0][0] + mat1[0][1] * mat2[1][0];
    result[0][1] = mat1[0][0] * mat2[0][1] + mat1[0][1] * mat2[1][1];
    result[1][0] = mat1[1][0] * mat2[0][0] + mat1[1][1] * mat2[1][0];
    result[1][1] = mat1[1][0] * mat2[0][1] + mat1[1][1] * mat2[1][1];
}

void scalarMat2x2(const matDouble &mat1, const double &scalar, matDouble &result)
{
    result[0][0] = scalar * mat1[0][0];
    result[0][1] = scalar * mat1[0][1];
    result[1][0] = scalar * mat1[1][0];
    result[1][1] = scalar * mat1[1][1];
}

void transposeMat2x2(const matDouble &matrix, matDouble &transpose)
{
    transpose[0][0] = matrix[0][0];
    transpose[0][1] = matrix[1][0];
    transpose[1][0] = matrix[0][1];
    transpose[1][1] = matrix[1][1];
}

void contractionMat2x2(const matDouble &mat1, const matDouble &mat2, double &result)
{
    result = mat1[0][0] * mat2[0][0] + mat1[0][1] * mat2[0][1] + mat1[1][0] * mat2[1][0] + mat1[1][1] * mat2[1][1];
}

void setZero2x2MatDouble(matDouble &mat)
{
    mat[0][0] = 0.0;
    mat[0][1] = 0.0;
    mat[1][0] = 0.0;
    mat[1][1] = 0.0;
}

void computeInverseAndDeterminant2x2(const matDouble &matrix, matDouble &inverseMatrix, double &det)
{
    det = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
    double detinv = 1.0 / det;
    inverseMatrix[0][0] = matrix[1][1] * detinv;
    inverseMatrix[1][0] = -matrix[1][0] * detinv;
    inverseMatrix[0][1] = -matrix[0][1] * detinv;
    inverseMatrix[1][1] = matrix[0][0] * detinv;
}

void setColumnToMatDoubleWith2Lines(matDouble &mat, const int &nColumn)
{
    mat.resize(2);
    mat[0].resize(nColumn, 0.0);
    mat[1].resize(nColumn, 0.0);
}

void setColumnToMatDoubleWith3Lines(matDouble &mat, const int &nColumn)
{
    mat.resize(3);
    mat[0].resize(nColumn, 0.0);
    mat[1].resize(nColumn, 0.0);
    mat[2].resize(nColumn, 0.0);
}

void setSizeMatDouble(matDouble &mat, const int &nLine, const int &nColumn, const double &value)
{
    mat.resize(nLine);
    for (int i = 0; i < nLine; i++)
    {
        mat[i].resize(nColumn, value);
    }
}

void setValueMatDouble(matDouble &mat, const double &value)
{
    for (int j = 0, nLine = mat.size(); j < nLine; j++)
    {
        for (int i = 0, nColumn = mat[0].size(); i < nColumn; i++)
        {
            mat[j][i] = value;
        }
    }
}

void setValueVecDouble(vecDouble &vec, const double &value)
{
    for (int j = 0, nLine = vec.size(); j < nLine; j++)
    {
        vec[j] = value;
    }
}

void sumMat2x2(const matDouble &mat1, const matDouble &mat2, matDouble &result)
{
    result[0][0] = mat1[0][0] + mat2[0][0];
    result[0][1] = mat1[0][1] + mat2[0][1];
    result[1][0] = mat1[1][0] + mat2[1][0];
    result[1][1] = mat1[1][1] + mat2[1][1];
}

std::vector<std::string> split(std::string str, std::string delim)
{
    std::istringstream is(str);
    std::vector<std::string> values;
    std::string token;
    while (getline(is, token, ' '))
        values.push_back(token);

    return values;
}

matDouble getTransposeMat2x2(const matDouble &matrix)
{
    matDouble transpose;
    setColumnToMatDoubleWith2Lines(transpose, 2);
    transpose[0][0] = matrix[0][0];
    transpose[0][1] = matrix[1][0];
    transpose[1][0] = matrix[0][1];
    transpose[1][1] = matrix[1][1];
    return transpose;
}

matDouble getProductMat2x2(const matDouble &mat1, const matDouble &mat2)
{
    matDouble result;
    setColumnToMatDoubleWith2Lines(result, 2);
    result[0][0] = mat1[0][0] * mat2[0][0] + mat1[0][1] * mat2[1][0];
    result[0][1] = mat1[0][0] * mat2[0][1] + mat1[0][1] * mat2[1][1];
    result[1][0] = mat1[1][0] * mat2[0][0] + mat1[1][1] * mat2[1][0];
    result[1][1] = mat1[1][0] * mat2[0][1] + mat1[1][1] * mat2[1][1];
    return result;
}

matDouble getSumMat2x2(const matDouble &mat1, const matDouble &mat2)
{
    matDouble result;
    setColumnToMatDoubleWith2Lines(result, 2);
    result[0][0] = mat1[0][0] + mat2[0][0];
    result[0][1] = mat1[0][1] + mat2[0][1];
    result[1][0] = mat1[1][0] + mat2[1][0];
    result[1][1] = mat1[1][1] + mat2[1][1];
    return result;
}

double getContractionMat2x2(const matDouble &mat1, const matDouble &mat2)
{
    double result = mat1[0][0] * mat2[0][0] + mat1[0][1] * mat2[0][1] + mat1[1][0] * mat2[1][0] + mat1[1][1] * mat2[1][1];
    return result;
}

void jumpLine(const int &nLines, std::ifstream &file)
{
    std::string line;
    for (int i = 0; i < nLines; i++)
    {
        std::getline(file, line);
    }
}

void subtractionMat3x3(const matDouble &mat1, const matDouble &mat2, matDouble &result)
{
    result[0][0] = mat1[0][0] - mat2[0][0];
    result[0][1] = mat1[0][1] - mat2[0][1];
    result[0][2] = mat1[0][2] - mat2[0][2];
    result[1][0] = mat1[1][0] - mat2[1][0];
    result[1][1] = mat1[1][1] - mat2[1][1];
    result[1][2] = mat1[1][2] - mat2[1][2];
    result[2][0] = mat1[2][0] - mat2[2][0];
    result[2][1] = mat1[2][1] - mat2[2][1];
    result[2][2] = mat1[2][2] - mat2[2][2];
}

void sumMat3x3(const matDouble &mat1, const matDouble &mat2, matDouble &result)
{
    result[0][0] = mat1[0][0] + mat2[0][0];
    result[0][1] = mat1[0][1] + mat2[0][1];
    result[0][2] = mat1[0][2] + mat2[0][2];
    result[1][0] = mat1[1][0] + mat2[1][0];
    result[1][1] = mat1[1][1] + mat2[1][1];
    result[1][2] = mat1[1][2] + mat2[1][2];
    result[2][0] = mat1[2][0] + mat2[2][0];
    result[2][1] = mat1[2][1] + mat2[2][1];
    result[2][2] = mat1[2][2] + mat2[2][2];
}

matDouble getSumMat3x3(const matDouble &mat1, const matDouble &mat2)
{
    matDouble result;
    setColumnToMatDoubleWith3Lines(result, 3);
    result[0][0] = mat1[0][0] + mat2[0][0];
    result[0][1] = mat1[0][1] + mat2[0][1];
    result[0][2] = mat1[0][2] + mat2[0][2];
    result[1][0] = mat1[1][0] + mat2[1][0];
    result[1][1] = mat1[1][1] + mat2[1][1];
    result[1][2] = mat1[1][2] + mat2[1][2];
    result[2][0] = mat1[2][0] + mat2[2][0];
    result[2][1] = mat1[2][1] + mat2[2][1];
    result[2][2] = mat1[2][2] + mat2[2][2];
    return result;
}

matDouble getScalarMat3x3(const matDouble &mat1, const double &scalar)
{
    matDouble result;
    setColumnToMatDoubleWith3Lines(result, 3);
    result[0][0] = scalar * mat1[0][0];
    result[0][1] = scalar * mat1[0][1];
    result[0][2] = scalar * mat1[0][2];
    result[1][0] = scalar * mat1[1][0];
    result[1][1] = scalar * mat1[1][1];
    result[1][2] = scalar * mat1[1][2];
    result[2][0] = scalar * mat1[2][0];
    result[2][1] = scalar * mat1[2][1];
    result[2][2] = scalar * mat1[2][2];
    return result;
}

double getContractionMat3x3(const matDouble &mat1, const matDouble &mat2)
{
    double result = mat1[0][0] * mat2[0][0] + mat1[0][1] * mat2[0][1] + mat1[0][2] * mat2[0][2] + mat1[1][0] * mat2[1][0] + mat1[1][1] * mat2[1][1] + mat1[1][2] * mat2[1][2] + mat1[2][0] * mat2[2][0] + mat1[2][1] * mat2[2][1] + mat1[2][2] * mat2[2][2];
    return result;
}