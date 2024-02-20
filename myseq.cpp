#include <stdio.h>
#include <bits/stdc++.h>

using namespace std;

void get_adj_matrix(float **graph, int n, float d, FILE *inputFilePtr)
{

    if (inputFilePtr == NULL)
    {
        printf("input.txt file failed to open.");
        return;
    }

    int m, indexing;
    fscanf(inputFilePtr, "%d", &m);
    fscanf(inputFilePtr, "%d", &indexing);

    for (int i = 0; i < n; i++)
    {
        graph[i] = (float *)malloc(sizeof(float) * n);
        for (int j = 0; j < n; ++j)
        {
            graph[i][j] = (1 - d) / float(n);
        }
    }

    cout << "index val " << indexing << endl;
    while (m--)
    {
        int source, destin;
        fscanf(inputFilePtr, "%d", &source);
        fscanf(inputFilePtr, "%d", &destin);

        cout << "src " << source << " "
             << "dest " << destin << endl;
        if (indexing == 0)
        {
            graph[destin][source] += d * 1.0; // transpose
        }
        else
        {
            graph[destin - 1][source - 1] += d * 1.0;
        }
    }
}

void manage_adj_matrix(float **graph, int n)
{

    for (int j = 0; j < n; ++j)
    {
        float sum = 0.0;

        for (int i = 0; i < n; ++i)
        {
            cout << graph[i][j] << " ";
            sum += graph[i][j];
        }
        cout << endl;
        cout << sum << " " << endl;

        for (int i = 0; i < n; ++i)
        {
            if (sum != 0.0)
            {
                graph[i][j] /= sum;
            }
            else
            {
                graph[i][j] = (1 / (float)n);
            }
        }
    }
}

void print_graph(float **graph, int n)
{

    cout << endl;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; ++j)
        {
            cout << graph[i][j] << "  ";
        }
        cout << endl;
    }
}
void printvector(float *r, int n)
{

    for (int i = 0; i < n; i++)
    {
        cout << r[i] << " ";
    }
    cout << endl;
}

float norm(float *vect, int n)
{
    float ans = 0.0;
    for (int i = 0; i < n; ++i)
    {
        ans += abs(vect[i]);
    }
    return ans;
}

void power_method(float **graph, float *r, int n, int max_iter = 1000, float eps = 0.000001)
{

    float *r_last = (float *)malloc(n * sizeof(float));

    for (int i = 0; i < n; i++)
    {
        r[i] = (1 / (float)n);
    }

    while (max_iter--)
    {

        for (int i = 0; i < n; ++i)
        {
            r_last[i] = r[i]; // copy r int r_last
        }

        cout << max_iter << ")"
             << "printing vector---" << endl;
        printvector(r, n);
        cout << max_iter << ")"
             << "printing last vector" << endl;
        printvector(r_last, n);

        for (int i = 0; i < n; ++i)
        { // multiply graph * R_last and store result in r
            float temp = 0.0;

            for (int j = 0; j < n; ++j)
            {
                temp += r_last[j] * graph[i][j]; // after 60 iteration value of r and r last is same
            }

            r[i] = temp;
        }

        cout << "---------------------------------------------------------------" << endl;

        cout << max_iter << ")"
             << "printing vector" << endl;
        printvector(r, n);
        cout << max_iter << ")"
             << "printing last vector" << endl;
        printvector(r_last, n);

        for (int i = 0; i < n; ++i)
        { // calculate diffrence betw r_last and r and store in r_last
            r_last[i] -= r[i];
        }

        cout << "***********printing last vector***********" << endl;
        printvector(r_last, n);
        cout << endl;

        if (norm(r_last, n) < eps)
        { // call ans should less than eps
            return;
        }
        
    }
    return;
}

void top_nodes(float *r, int n, int count = 4){

    priority_queue<pair<float, int>> pq; 

    for(int i = 0; i< n; ++i){
        pq.push(make_pair(r[i], i+ 1));
    }
    int rank =1;
    while(rank <= count){
        printf("Rank %d Node is %d and %lf\n", rank, pq.top().second ,pq.top().first);
        rank++;
        pq.pop();
    }

}


int main(int argc, char **argv)
{
    FILE *Inputfileptr;
    char *inputfile = argv[1];
    Inputfileptr = fopen(inputfile, "r");
    int total_size;

    fscanf(Inputfileptr, "%d", &total_size);
    float d = 0.85;

    float **graph = (float **)malloc(total_size * sizeof(float *));

    float *r = (float *)malloc(total_size * sizeof(float));

    get_adj_matrix(graph, total_size, d, Inputfileptr);
    cout<<endl;
    print_graph(graph, total_size);
    manage_adj_matrix(graph, total_size);
    print_graph(graph, total_size);
    power_method(graph, r, total_size);
    printvector(r, total_size);
    top_nodes(r, total_size);

    return 0;
}
