#include <iostream>
#include <chrono>
#include <omp.h>
#include <vector>
#include <cmath>
#include <array>
#include <algorithm>
#include <fstream>

//use the GLS C++ library
/* #include "/usr/local/include/gsl/gsl_sf_bessel.h"
#include "../packages/gsl-2.6/specfunc/gsl_sf_bessel.h"
 */
using namespace std;

void getNumberOfThreads()
{
    //int omp_get_thread_num();
#pragma omp parallel
    {
        printf("Thread rank: %d\n", omp_get_thread_num());
    }
}

void checkParallelIfAvaiable()
{
    int count = 0;
#pragma omp parallel
    {
        int id = omp_get_thread_num();
        int data = id;
        int total = omp_get_num_threads();
        printf("Greetings from process %d out of %d with Data %d\n", id, total, data);
        count++;
    }
    printf("this is count %d \n", count);
}

ofstream outparallel("../output/outputParallel.dat");
ofstream outserial("../output/outputSerial.dat");

double myFunction(double a, double b, double c)
{
    return (double)a * log(b) + exp(a / b) + c;
}

double linearFunction(double a1, double a2, double a3)
{
    int n = 34;
    bool ok = 1;
    double sum1 = 0, sum2 = 0, sum3 = 0;
    while (ok)
    {
        for (int i = 1; i <= n && ok; ++i)
        {
            double temp = sqrt(a1 * sqrt(myFunction(a1 * i, a2, a3)) * log(2 + i) + a1 * a1 * i);
            if (isnan(temp))
            {
                ok = 0;
                break;
            }
            if (!isnan(temp))
                sum1 += temp;
        }
        for (int i = 1; i <= n - 5 && ok; ++i)
        {
            double temp = sqrt(a2 * sqrt(myFunction(a1, a2 * i, a3)) * i + a2 * i);
            if (isnan(temp))
            {
                ok = 0;
                break;
            }
            if (!isnan(temp))
                sum2 += temp;
        }
        for (int i = 1; i <= n + 7 && ok; ++i)
        {
            double temp = sqrt(a3 * sqrt(myFunction(a1, a2, a3 * i)) + a3 * a3 - i * log(i));
            if (isnan(temp))
            {
                ok = 0;
                break;
            }
            if (!isnan(temp))
                sum3 += temp;
        }
        if (ok)
            ok = 0;
    }
    if (sum1 && sum2 && sum3)
        return (double)abs(sum1 + sum2 + sum3);
    return 1;
}

void serialAlgorithm(int spaceDIM)
{
    auto sum = 0;
    auto timeStart = chrono::high_resolution_clock::now();
    for (int i = 1; i <= spaceDIM; ++i)
        for (int j = 1; j <= spaceDIM; ++j)
        {
            sum = i * j;
            //#pragma omp ordered
            outserial << i << " " << j << " " << sum << endl;
            //cout << i << " " << j << " " << sum << endl;
        }
    auto timeStop = chrono::high_resolution_clock::now();
    auto execTime = chrono::duration_cast<chrono::seconds>(timeStop - timeStart).count();
    cout << "Serial execution time = " << execTime << " seconds";
    cout << endl;
    outserial << "Execution time = " << execTime << " seconds";
    outserial << endl;
}

//comaprison segment 1
void serialAlgorithmTripleLoop(int spaceDIM)
{
    auto sum = 0;
    vector<double> fArray;
    auto timeStart = chrono::high_resolution_clock::now();
    for (int i = 1; i <= spaceDIM; ++i)
        for (int j = 1; j <= spaceDIM; ++j)
            for (int k = 1; k <= spaceDIM; ++k)
            {
                //sum += i * j * k;
                fArray.push_back(myFunction((double)sqrt(i), (double)sqrt(j), (double)sqrt(k)));
                //outserial << i << " " << j << " " << k << " " << myFunction((double)i, (double)j, (double)k) << "\n";
            }
    auto timeStop = chrono::high_resolution_clock::now();
    auto execTime = chrono::duration_cast<chrono::milliseconds>(timeStop - timeStart).count();
    cout << "Serial execution time = " << (double)execTime / 1000 << " seconds";
    cout << endl;
    outserial << "Execution time = " << (double)execTime / 1000 << " seconds";
    outserial << endl;
}

void parallelAlgorithm(int spaceDIM)
{
    //start of the actual algorithm
    auto sum = 0;
    auto timeStart = chrono::high_resolution_clock::now();
#pragma omp parallel
    {
#pragma omp for ordered
        for (int i = 1; i <= spaceDIM; ++i)
            for (int j = 1; j <= spaceDIM; ++j)
            {
                sum = i * j;
#pragma omp ordered
                outparallel << i << " " << j << " " << sum << endl;
                //cout << i << " " << j << " " << sum << endl;
            }
    }
    auto timeStop = chrono::high_resolution_clock::now();
    auto execTime = chrono::duration_cast<chrono::seconds>(timeStop - timeStart).count();
    cout << "Execution time = " << execTime << " seconds";
    cout << endl;
    outparallel << "Execution time = " << execTime << " seconds";
    outparallel << endl;
}
//comaprison segment 2
void parallelAlgorithmTripleLoop(int spaceDIM)
{
    //start of the actual algorithm
    auto sum = 0;
    auto timeStart = chrono::high_resolution_clock::now();
    vector<double> fArray;
#pragma omp parallel for ordered schedule(dynamic, 2)
    for (int i = 1; i <= spaceDIM; ++i)
        for (int j = 1; j <= spaceDIM; ++j)
            for (int k = 1; k <= spaceDIM; ++k)
            {
                double temp = myFunction((double)i, (double)j, (double)k);
//   sum += i * j * k;
// vector<double> fArray;
#pragma omp ordered
                fArray.push_back(temp);
                //              outparallel << i << " " << j << " " << myFunction((double)i, (double)j, (double)k) << "\n";
            }
    auto timeStop = chrono::high_resolution_clock::now();
    auto execTime = chrono::duration_cast<chrono::milliseconds>(timeStop - timeStart).count();
    cout << "Parallel execution time = " << (double)execTime / 1000 << " seconds";
    cout << endl;
    outparallel << "Execution time = " << (double)execTime / 1000 << " seconds";
    outparallel << endl;
}

typedef vector<double> myArray;
void linearAlgoritm(const int dim, const int innerDim, double *avgExecTime, double *avgIOTime)
{
    auto startTime = chrono::high_resolution_clock::now();
    myArray v1(innerDim);
    for (int a = 1; a <= innerDim; ++a)
    {
        double temp = linearFunction(a, a + 1.1, a + 0.7);
        if (temp)
            v1.push_back(temp);
    }
    auto endTime = chrono::high_resolution_clock::now();
    auto execTime = chrono::duration_cast<chrono::milliseconds>(endTime - startTime).count();
    startTime = chrono::high_resolution_clock::now();
    for (auto i = 0; i < v1.size(); ++i)
        outserial << i << " " << v1.at(i) << "\n";
    endTime = chrono::high_resolution_clock::now();
    auto ioTime = chrono::duration_cast<chrono::milliseconds>(endTime - startTime).count();
    cout << "Execution time of the algorithm is: " << (double)execTime / 1000.0 << " seconds\n";
    cout << "Execution time of the I/O process is: " << (double)ioTime / 1000.0 << " seconds\n";
    cout << endl;
    outserial << "Execution time of the algorithm is: " << (double)execTime / 1000.0 << " seconds\n";
    outserial << "Execution time of the I/O process is: " << (double)ioTime / 1000.0 << " seconds\n";
    outserial << endl;
    *avgIOTime = (double)ioTime / 1000.0;
    *avgExecTime = (double)execTime / 1000.0;
}

void multipleExecutions(int nExec, const int dim, const int innerDim)
{
    double avgExecTime;
    double avgIOTime;
    myArray avgExec(nExec), avgIO(nExec);
    int count = 0;
    const int executions = nExec;
    while (nExec)
    {
        linearAlgoritm(dim, innerDim, &avgExecTime, &avgIOTime);
        avgExec[count] = avgExecTime;
        avgIO[count] = avgIOTime;
        nExec--;
        count++;
    }
    auto sum = 0.0;
    for_each(avgIO.begin(), avgIO.end(), [&](auto d) {
        sum += d;
    });
    cout << "Average IO time for " << executions << " executions: ";
    cout << (double)sum / avgIO.size() << "\n";
    sum = 0.0;
    for_each(avgExec.begin(), avgExec.end(), [&](auto d) {
        sum += d;
    });
    cout << "Average execution time for " << executions << " executions: ";
    cout << (double)sum / avgExec.size();
    cout << endl;
}

int main()
{
    cout << "FUNCTION OPTIMIZATION" << endl;
    //serialAlgorithm();
    //parallelAlgorithm();
    //check for parallel processing available
    //getNumberOfThreads();
    //checkParallelIfAvaiable();
    //serialAlgorithmTripleLoop(dim);
    //parallelAlgorithmTripleLoop(dim);
    multipleExecutions(10, 1000, 123450);
    return 0;
}