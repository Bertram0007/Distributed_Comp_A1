#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <pthread.h>
#include <math.h>

struct thread_data
{
    double** matrix;
    double* result;
    double* blockingResult;
    int n;
    int m;
    int b;
    int startN;
    int startSubN;
    int number;
    int threadCounter;
    int blocks;
    int startNumA;
    int startIterateA;
};

void test(double **matrix, double *testResult, int n, int m){
    int testCounter = 0;
    for(int row=0; row<n; row++){
        for(int iterateRow = row; iterateRow < n; iterateRow++){
            double temp = 0;
            for(int column=0; column<m; column++){
                temp += matrix[row][column] * matrix[iterateRow][column];
            }
            testResult[testCounter++] = temp;
        }
    }
}

void *baseline(void *data){
    struct thread_data *p = (struct thread_data *)data;
    int row, iterateRow;
    int index = 0;
    int counter = 0;
    int numberTemp = p->number;
    for(row=p->startN; row<p->n; row++){
        for(iterateRow = row == p->startN ? p->startSubN : row ;iterateRow < p->n; iterateRow++){
            counter += 2;
            numberTemp--;
            double temp = 0;
            for(int column=0; column<p->m; column++){
                temp += p->matrix[row][column] * p->matrix[iterateRow][column];
            }
            index=0;
            int tempN = p->n-1;
            if(tempN>0){
                for(int i=0; i<row; i++){
                    index += tempN--;
                }
            }
            index+=iterateRow;
            p->result[index] = temp;
            if(numberTemp == 0){
                return NULL;
            }
        }
    }
    return NULL;
}

void *blockingAndLoopUnrolling(void *data){
    struct thread_data *p = (struct thread_data *)data;
    int blockCounter = 0, iterateA = p->startIterateA;
    for (int numA = p->startNumA; numA < p -> blocks; numA++){
        while (iterateA < p -> blocks){
            if (numA > iterateA){
                continue;
            }
            for (int row = numA * p->b; row < (numA+1) * p->b; row++){
                if(row >= p->n){
                    continue;
                }
                int tempN = p->n;
                int indexRow = 0;
                for(int i=0; i<row; i++){
                    indexRow += tempN--;
                }
                for (int transposeRow = iterateA * p->b; transposeRow < (iterateA+1) * p->b; transposeRow++){
                    if (row > transposeRow || transposeRow >= p->n){
                        continue;
                    }
                    double sum = 0;
                    for (int transposeColumn = 0; transposeColumn < p->m; transposeColumn++) {   //iterate the columns for transpose matrix
                        if (transposeColumn + 4 < p->m) {
                            sum += p->matrix[row][transposeColumn] * p->matrix[transposeRow][transposeColumn];
                            sum += p->matrix[row][transposeColumn + 1] * p->matrix[transposeRow][transposeColumn + 1];
                            sum += p->matrix[row][transposeColumn + 2] * p->matrix[transposeRow][transposeColumn + 2];
                            sum += p->matrix[row][transposeColumn + 3] * p->matrix[transposeRow][transposeColumn + 3];
                            sum += p->matrix[row][transposeColumn + 4] * p->matrix[transposeRow][transposeColumn + 4];
                            transposeColumn += 4;
                        } else if (transposeColumn + 3 < p->m) {
                            sum += p->matrix[row][transposeColumn] * p->matrix[transposeRow][transposeColumn];
                            sum += p->matrix[row][transposeColumn + 1] * p->matrix[transposeRow][transposeColumn + 1];
                            sum += p->matrix[row][transposeColumn + 2] * p->matrix[transposeRow][transposeColumn + 2];
                            sum += p->matrix[row][transposeColumn + 3] * p->matrix[transposeRow][transposeColumn + 3];
                            transposeColumn += 3;
                        } else if (transposeColumn + 2 < p->m) {
                            sum += p->matrix[row][transposeColumn] * p->matrix[transposeRow][transposeColumn];
                            sum += p->matrix[row][transposeColumn + 1] * p->matrix[transposeRow][transposeColumn + 1];
                            sum += p->matrix[row][transposeColumn + 2] * p->matrix[transposeRow][transposeColumn + 2];
                            transposeColumn += 2;
                        } else if (transposeColumn + 1 < p->m) {
                            sum += p->matrix[row][transposeColumn] * p->matrix[transposeRow][transposeColumn];
                            sum += p->matrix[row][transposeColumn + 1] * p->matrix[transposeRow][transposeColumn + 1];
                            transposeColumn += 1;
                        } else {
                            sum += p->matrix[row][transposeColumn] * p->matrix[transposeRow][transposeColumn];
                        }
                    }
                    int index = indexRow + transposeRow - row;
                    p -> blockingResult[index] = sum;
                }
            }
            if (++blockCounter == p->number){
                return NULL;
            }
            iterateA++;
        }
        if (iterateA == p -> blocks){
            iterateA = numA +1;
        }
    }
    return NULL;
}


int main(int argc, char *argv[]){
    struct timeval
            start,
            end;
    int N, M, T, B;
    if(argc == 5){
        N = atoi(argv[1]);
        M = atoi(argv[2]);
        T = atoi(argv[3]);
        B = atoi(argv[4]);
        printf("N = %d, M = %d, T = %d, B = %d\n\n", N, M, T, B);
    }
    else{
        printf("Usage: %s N M T B\n\n"
               " N: matrix row length\n"
               " M: matrix column length\n"
               " T: total threads\n"
               " B: block size\n\n",argv[0]);
        return 1;
    }
    struct thread_data thread_data_array_blocking[T];
    struct thread_data thread_data_array[T];

    double* matrix0 = (double*)malloc(N*M*sizeof(double));
    double** matrix = (double**)malloc(N*sizeof(double*));
    for (int i = 0; i < N; i++)
    {
        matrix[i] = matrix0 + i * M;
    }
    srand(time(0));
    for (int i = 0; i < N; i++){
        for (int j = 0; j < M; j++){
            matrix[i][j] = (double)rand()/RAND_MAX;
        }
    }

    double total = N*(N+1)/2;   //total number of multiplications
    double each = N*(N+1)/2/(double)T;   //number of multiplications for each thread
    int numberCeil = ceil(each);  //the ceiling number of multiplications for each thread
    int numberFloor = floor(each);   //the floor number of multiplications for each thread
    int number = 0;     //will be assigned with numberCeil/numberFloor
    int startN = 0, startSubN = 0;    //start index for the outer for loop and inner for loop
    int counter = 0;

    srand(time(0));
    gettimeofday(&start, 0);       //used to calculate the running time
    pthread_t tid[T];
    double *result = (double*)malloc(N*(N+1)/2*sizeof(double));

    for(int i=0; i<T; i++){
        if((total-counter)/(T-i) > each){
            number = numberCeil;
        }else{
            number = numberFloor;
        }
        counter += number;
        thread_data_array[i].matrix = matrix;
        thread_data_array[i].result = result;
        thread_data_array[i].n = N;
        thread_data_array[i].m = M;
        thread_data_array[i].b = B;
        thread_data_array[i].number = number;
        thread_data_array[i].startN = startN;
        thread_data_array[i].startSubN = startSubN;
        pthread_create(&tid[i], NULL, baseline, &thread_data_array[i]);
        startSubN += number;
        while(startSubN >= N && startN <= N-1){
            startN++;
            startSubN = startSubN - N + startN;
            if(startN >= N || (startN == N-1 && startSubN > N-1)){
                break;
            }
        }
    }
    for (int i = 0; i < T; i++){
        pthread_join(tid[i], NULL);
    }
    gettimeofday(&end, 0);
    long seconds = end.tv_sec - start.tv_sec;
    long microseconds = end.tv_usec - start.tv_usec;
    double period = seconds + microseconds/1000000.0;
    printf("time period for baseline:   %8f \n", period);

    int Ablocks = ceil((double)N/(double)B);
    double *blockingResult = (double*)malloc(N*(N+1)/2*sizeof(double));
    pthread_t thread[T];

    gettimeofday(&start, 0);
    for (int i = 0; i < T; i++){
        int startNumA, startIterateA;
        int blockNum = Ablocks*(Ablocks+1)/2;
        int each = blockNum / T;
        int start = blockNum % T == 0 ?  i * each : i * (each + 1);
        int number = blockNum % T == 0? each : each + 1;
        int numA = 0;
        while ( numA < Ablocks){
            int tempN = Ablocks;
            int index = 0;
            for(int i=0; i<numA; i++){
                index += tempN--;
            }
            int iterateA = numA;
            while (iterateA < Ablocks){
                if (index+iterateA-numA == start){
                    startNumA = numA;
                    startIterateA = iterateA;
                    iterateA = Ablocks;
                    numA = Ablocks;
                }
                iterateA++;
            }
            numA++;
        }
        thread_data_array_blocking[i].threadCounter = i;
        thread_data_array_blocking[i].b = B;
        thread_data_array_blocking[i].n = N;
        thread_data_array_blocking[i].m = M;
        thread_data_array_blocking[i].matrix = matrix;
        thread_data_array_blocking[i].blockingResult = blockingResult;
        thread_data_array_blocking[i].startNumA = startNumA;
        thread_data_array_blocking[i].startIterateA = startIterateA;
        thread_data_array_blocking[i].blocks = Ablocks;
        thread_data_array_blocking[i].number = number;
        pthread_create(&thread[i], NULL, blockingAndLoopUnrolling, &thread_data_array_blocking[i]);
    }
    for (int i = 0; i < T; i++){
        pthread_join(thread[i], NULL);
    }
    gettimeofday(&end, 0);

    seconds = end.tv_sec - start.tv_sec;
    microseconds = end.tv_usec - start.tv_usec;
    period = seconds + microseconds/1000000.0;
    printf("time period for blocking and loop unrolling:   %8f\n", period);

    double *testResult = (double*)malloc(N*(N+1)/2*sizeof(double));
    test(matrix, testResult, N, M);
    for(int i=0; i<N*(N+1)/2; i++){
        if(testResult[i] != blockingResult[i]){
            printf("\nblockingResult %d false", i);
        }else if(testResult[i] != result[i]){
            printf("\nresult %d false", i);
        }
    }
    return 0;
}



