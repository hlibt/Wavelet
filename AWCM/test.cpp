#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <vector>
using namespace std;

class Point {
    public:
        int x;
};

int main() {
    int n=5;
    int m=10;
    Point** mypnt = new Point*[n];
    for (int i=0;i<n;i++) {
        mypnt[i] = new Point[m];
    }
    for (int i=0;i<n;i++) {
        for (int j=0;j<m;j++) {
            mypnt[i][j].x = i*10 + j;
        }
    }
    for (int i=0;i<n;i++) {
        for (int j=0;j<m;j++) {
            cout << mypnt[i][j].x << endl;
        }
    }
}
