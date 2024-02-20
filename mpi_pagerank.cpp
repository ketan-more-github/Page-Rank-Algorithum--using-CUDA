#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>
#include <ctime>

using namespace std;

// Function declarations
void get_adj_matrix(float** graph, int n, float d, FILE *inputFilePtr, int local_start, int local_end);
void manage_adj_matrix(float** graph, int n);
float norm(float *vect, int n);
void power_method(float **graph, float *r, int n, int max_iter = 1000, float eps = 0.000001);
void top_nodes(float *r, int n, int count = 3);

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc != 2) {
        if (rank == 0)
            cout << "Usage: mpirun -np <num_processes> ./executable <input_file>\n";
        MPI_Finalize();
        return 1;
    }

    FILE *inputFilePtr;
    char * inputfile = argv[1];
    inputFilePtr = fopen(inputfile, "r");
    if (!inputFilePtr) {
        if (rank == 0)
            cout << "Failed to open input file.\n";
        MPI_Finalize();
        return 1;
    }

    int n; 
    fscanf(inputFilePtr, "%d", &n);
    int local_n = n / size;
    int local_start = rank * local_n;
    int local_end = local_start + local_n;

    float d = 0.85; 

    float** graph = (float**)malloc(local_n * sizeof(float*));
    for (int i = 0; i < local_n; ++i) {
        graph[i] = (float*)malloc(n * sizeof(float));
    }

    float* r = (float*)malloc(local_n * sizeof(float));

    get_adj_matrix(graph, n, d, inputFilePtr, local_start, local_end);

    manage_adj_matrix(graph, local_n);

    double start_time = MPI_Wtime();

    power_method(graph, r, local_n);

    double end_time = MPI_Wtime();

    if (rank == 0) {
        float* global_r = (float*)malloc(n * sizeof(float));
        MPI_Gather(r, local_n, MPI_FLOAT, global_r, local_n, MPI_FLOAT, 0, MPI_COMM_WORLD);
        if (rank == 0) {
            top_nodes(global_r, n);
        }
        free(global_r);
    } else {
        MPI_Gather(r, local_n, MPI_FLOAT, NULL, 0, MPI_FLOAT, 0, MPI_COMM_WORLD);
    }

    free(r);
    for (int i = 0; i < local_n; ++i) {
        free(graph[i]);
    }
    free(graph);

    fclose(inputFilePtr);

    if (rank == 0) {
        cout << "Time taken: " << end_time - start_time << " seconds for parallel implementation with " << n << " nodes.\n";
    }

    MPI_Finalize();
    return 0;
}

void get_adj_matrix(float** graph, int n, float d, FILE *inputFilePtr, int local_start, int local_end) {
    // Read data from the input file
    int m, indexing;
    fscanf(inputFilePtr, "%d", &m);
    fscanf(inputFilePtr, "%d", &indexing);

    for(int i = local_start; i < local_end; ++i) {
        for(int j = 0; j < n; ++j) {
            graph[i - local_start][j] = (1 - d) / float(n);
        }
    }

    while(m--) {
        int source, destin;
        fscanf(inputFilePtr, "%d", &source);
        fscanf(inputFilePtr, "%d", &destin);
        if (indexing == 0) {
            if (source >= local_start && source < local_end) {
                graph[source - local_start][destin] += d * 1.0;
            }
        } else {
            if (source >= local_start && source < local_end) {
                graph[source - local_start][destin - 1] += d * 1.0;
            }
        }
    }
}

void manage_adj_matrix(float** graph, int n) {
    // Manage adjacency matrix
    for(int j = 0; j < n; ++j) {
        float sum = 0.0;

        for (int i = 0; i < n; ++i) {
            sum += graph[i][j];
        }

        for (int i = 0; i < n; ++i) {
            if (sum != 0.0) {
                graph[i][j] /= sum;
            } else {
                graph[i][j] = (1/(float)n);
            }
        }
    }
}

float norm(float *vect, int n) {
    float ans = 0.0;
    for (int i = 0; i < n; ++i) {
        ans += abs(vect[i]);
    }
    return ans;
}

void power_method(float **graph, float *r, int n, int max_iter, float eps) {
    float* r_last = (float*)malloc(n * sizeof(float));

    for(int i = 0; i < n; ++i) {
        r[i] = (1/(float)n);
    }

    while(max_iter--) {
        for(int i = 0; i < n; ++i) {
            r_last[i] = r[i];
        }
        for(int i = 0; i < n; ++i) {
            float sum = 0.0;

            for (int j = 0; j < n; ++j) {
                sum += r_last[j] * graph[i][j];
            }

            r[i] = sum;
        }

        for(int i = 0; i < n; ++i) {
            r_last[i] -= r[i];
        }

        if(norm(r_last, n) < eps) {
            free(r_last);
            return;
        }
    }
    free(r_last);
    return;
}

void top_nodes(float *r, int n, int count) {
    priority_queue<pair<float, int>> pq;

    for(int i = 0; i < n; ++i) {
        pq.push(make_pair(r[i], i + 1));
    }
    int rank = 1;
    while(rank <= count) {
        printf("Rank %d Node is %d\n", rank, pq.top().second);
        rank++;
        pq.pop();
    }
}
