#include<iostream>
using namespace std;
int main() {
int a;
double b=2.76;
double c;
for (int i=0;i<3;i++) {
    a=i>0;
    c=a*b+a;
    cout << c;
}
return 0;
}
