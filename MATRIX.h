#ifndef _My_Matrix_H
#define _My_Matrix_H

#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

class Matrix
{
  public:
    Matrix();
    Matrix( int rows, int cols );  // size rows x cols
 	~Matrix();   
    int row();  //number of rows
    int col();  //number of columns
    double dtm(); //determinant
    int display(ostream& outs=cout);
	int assign(int,int,double val); //assign value to matrix item
    int resize( int newRows, int newCols ); // resizes matrix to newRows x newCols
    vector<double>& operator[] (int row);  //assign value to matrix item
    Matrix operator+ (Matrix m2); //matrix addition
    Matrix operator- (Matrix m2); //matrix subtraction
    Matrix operator* (Matrix m2); //matrix multiplication
	Matrix operator* (double a);
    Matrix operator^ (int); //inverse of self-multiplication
    friend ostream& operator<< (ostream& outs, Matrix m2); //matrix output
    Matrix inverse();
	void rotate(Matrix &, Matrix &,double);
    int nRows;    // # of rows (capacity)
    int nCols;    // # of cols (capacity)
    vector< vector <double> > Mtrx; // the matrix of items 
	int check_range (int row); //range protection 
	int rank;     //rank of matrix
	double det;   //determinant   
	int sol(Matrix &, Matrix &,Matrix &);
	int func(Matrix &, Matrix &, Matrix &,Matrix &);
	int dfun(Matrix &, Matrix &, Matrix &,Matrix &);
	int ludcmp(Matrix &a,int *indx,double &d);  //LU decomposition
    int lubksb(Matrix &a,int *indx,double *b);  //LU back substitution
};  // end Matrix class specification

#endif
