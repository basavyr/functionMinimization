#include <iostream>
#include <chrono>
#include <omp.h>
#include <cmath>

using namespace std;

void getNumberOfThreads()
{
    int omp_get_thread_num();
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

int main()
{
    cout << "FUNCTION OPTIMIZATION" << endl;
    auto sum = 0;
    //check for parallel processing available
    //getNumberOfThreads();
    //checkParallelIfAvaiable();

    //start of the actual algorithm
    auto timeStart = chrono::high_resolution_clock::now();
#pragma omp parallel
    {
#pragma omp for
        for (int i = 1; i <= 200000000; ++i)
            sum += i * log(i);
    }
    auto timeStop = chrono::high_resolution_clock::now();
    auto execTime = chrono::duration_cast<chrono::seconds>(timeStop - timeStart).count();
    cout << "Execution time = " << execTime << " seconds";
    cout << endl;
    return 0;
}