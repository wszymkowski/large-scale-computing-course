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

void printVector(vector<int> vec){
    for(int i=0; i<vec.size(); i++){
        cout<< vec[i] << " ";
    }
    cout<<endl;
}

vector<int> singleSiveBlock(const vector<int> &primes, const int from, const int to){
    int elements = (to - from) + 1;
    int* tab_numbers = new int [elements];
    fill_tab_with_range(tab_numbers, from, to);
    vector<int> results;

    for(int prime : primes){
        int low_limit = floor(from / prime) * prime;
        if(low_limit < from) low_limit += prime; 
        for(int i = low_limit; i <= to; i += prime){
            tab_numbers[i - from] = -1;
        }
    }
    for(int i = 0; i < elements; i++){
        if(tab_numbers[i] != -1){
            results.push_back(tab_numbers[i]);
        }
    }
    
    return results;
}

int main(int argc, char *argv[]){

    int start = atoi(argv[1]);
    int end = atoi(argv[2]);
    int sliceSize = atoi(argv[3]);
    int num_threads = atoi(argv[4]);
    omp_set_num_threads(num_threads);
    
    int limit = ceil(sqrt(end));
    int elements = (end-start) + 1;
    vector<int> primes;
    vector<int> result;
    
    simpleSieve(limit, primes);

    wall_time_start = omp_get_wtime();
    time_start = clock();
    
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic)
        for (int from = start; from <= end; from += sliceSize){
            int to = from + sliceSize;
            if (to > end) to = end;
            
            vector<int> slice_result = singleSiveBlock(primes, from, to);
            #pragma omp critical
            {
                result.insert(result.end(),slice_result.begin(), slice_result.end());
            }
        }
    }
  
    time_stop = clock();
    wall_time_stop = omp_get_wtime();
    
    cout << "Czas przetwarzania zegarowego: " << (wall_time_stop - wall_time_start) << endl;
    cout << "Czas przetwarzania procesora: " << double(time_stop-time_start)/CLOCKS_PER_SEC << endl;
    cout << "Czas przetwarzania jednego procesora: " << (double(time_stop-time_start)/CLOCKS_PER_SEC) / num_threads << endl;

    return 0;
}