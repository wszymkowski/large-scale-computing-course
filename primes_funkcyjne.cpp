#include <iostream>
#include <cstdlib>
#include <omp.h>
#include <math.h>
#include <vector>
#include <time.h>
#include <algorithm>
#include <bits/stdc++.h>


using namespace std;

clock_t time_start, time_stop;
double wall_time_start, wall_time_stop;


void fill_tab_with_range(int tab[], int start, int end){
    int value = start;
    for(int i=0;i<=(end-start);i++){
        tab[i]=value;
        value++;
    }
}

void simpleSieve(int limit, vector<int>& primes){
    bool mark[limit+1];
    memset(mark, true, sizeof(mark));
    for(int i = 2; i <= limit; i++){
        if(mark[i] == true){
            primes.push_back(i);
            for(int j = i; j<= limit; j += i){
                mark[j] = false;
            }
        }
    }
}

void primesInRange(int low, int high, const vector<int> &primes, int tab[]){
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic)
        for(int prime : primes){
            int low_limit = floor(low / prime) * prime;
            if(low_limit < low) low_limit += prime; 
            for(int i = low_limit; i <= high; i += prime){
                if(tab[i - low] != -1){
                    tab[i - low] = -1;
                }
            }
        }
    }
}

void printVector(vector<int> vec){
    for(int i=0; i<vec.size(); i++){
        cout<< vec[i] << " ";
    }
    cout<<endl;
}

void printPrimes(int tab[], int elements){
    for(int i=0; i< elements; i++){
        if(tab[i] != -1){
            cout << tab[i] << " ";
        }
    }
    cout << endl;
}

int main(int argc, char *argv[]){
    int start = atoi(argv[1]);
    int end = atoi(argv[2]);
    int num_threads = atoi(argv[3]);
    omp_set_num_threads(num_threads);
    int limit = ceil(sqrt(end));
    int elements = (end-start) + 1;
    vector<int> primes;
    int* tab_numbers = new int[elements];
    

    fill_tab_with_range(tab_numbers, start, end);
    simpleSieve(limit, primes);

    wall_time_start = omp_get_wtime();
    time_start = clock();
    primesInRange(start, end, primes, tab_numbers);
    time_stop = clock();
    wall_time_stop = omp_get_wtime();

    //printPrimes(tab_numbers, elements);
    
    
    cout << "Czas przetwarzania zegarowego: " << (wall_time_stop - wall_time_start) << endl;
    cout << "Czas przetwarzania procesora: " << double(time_stop-time_start)/CLOCKS_PER_SEC << endl;
    cout << "Czas przetwarzania jednego procesora: " << (double(time_stop-time_start)/CLOCKS_PER_SEC) / num_threads << endl;

    return 0;
}