#include<iostream>
using namespace std;
const int MAX=3;
main() {
int var[MAX]={10,100,200};
int *ptr[MAX];
for (int i=0;i<MAX;i++) {
    ptr[i]=&var[i];
    cout<< *ptr[i]<<endl;
}
return 0;
}
