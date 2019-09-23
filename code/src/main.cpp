#include <iostream>
#include <chrono>
#include <omp.h>
#include <cmath>

using namespace std;

int main()
{
    auto sum = 0;
    auto timeStart = chrono::high_resolution_clock::now();
    for (int i = 1; i <= 1000000000; ++i)
    {
        sum += i * log(i);
    }
    auto timeStop = chrono::high_resolution_clock::now();
    auto execTime = chrono::duration_cast<chrono::seconds>(timeStop - timeStart).count();
    cout << "FUNCTION OPTIMIZATION" << endl;
    cout << "Execution time = " << execTime << " seconds";
    cout << endl;
    return 0;
}