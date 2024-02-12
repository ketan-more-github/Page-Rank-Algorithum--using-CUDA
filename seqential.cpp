#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;
#define N 6

float *mat_vect(float *m, float *v)
{

    float *r;
    r = (float *)malloc(N * sizeof(float));
    for (int i = 0; i < N; i++)
    {
        r[i] = 0;
        for (int j = 0; j < N; j++)
        {
            r[i] += m[i * N + j] * v[i];
        }
    }

    return r;
}

int main()
{

    float *matrix;
    float *vector;
    float *result;
    float *result1;
    float *result2;
    int i, j;
    float elem;

    matrix = (float *)malloc(N * N * sizeof(float));
    vector = (float *)malloc(N * sizeof(float));
    result1 = (float *)malloc(N * sizeof(float));
    result2 = (float *)malloc(N * sizeof(float));

    ifstream is("input.txt");

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            // cout<<"enter the matrix element"<<endl;
            // cin>>elem;
            float a;
            is >> a;
            matrix[i * N + j] = a;
        }
    }

    for (int i = 0; i < N; i++)
    {
        vector[i] = i + 1;
    }

    cout << "vector elem are " << endl;

    for (i = 0; i < N; i++)
    {
        cout << vector[i] << " ";
    }
    cout << endl;
    cout << "matrix elem are" << endl;
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            cout << i * N + j << ")" << matrix[i * N + j] << "   ";
        }
    }
    cout << endl;

        result = mat_vect(matrix, vector);
     

    for (i = 0; i < N; i++)
    {
        cout << result[i] << endl;
    }
}