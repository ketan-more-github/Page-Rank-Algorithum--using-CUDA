#include<iostream>
#include<fstream>
using namespace std;

int main()
{

    int arr[30];
    ifstream is("input.txt");

for (int i = 0; i < 36; i++)
{

    
    float a;
    is >> a ;
    cout<<a<< " ";    /* code */
}

    // int x;
    // int cnt =0;

    // while(cnt < arr[36]  && is >> x)
    // {
    //     arr[cnt++] = x;
    // }

    //  for(int i = 0 ; i <cnt ;i++)
    // {
    //     cout<<arr[i]<<" ";
    // }
    // cout<<endl;
    // is.close();
  
    return 0;
}
