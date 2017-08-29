#include <iostream>
#include <math.h>
#include <cmath>
using namespace std;

double* partial_pivot(double** A,int size);
double* complete_pivot(double** A,int size);

int main() {
	double** A=new double*[2];
	for (int i=0;i<2;i++) A[i]= new double[3];
	A[0][0]=2.;
	A[0][1]=4.;
	A[1][0]=6.1;
	A[1][1]=6.2;
	for (int i=0;i<2;i++) {
		for (int j=0;j<2;j++) {
			cout << *(*(A+i)+j) << endl;
		}
	}
	A[0][2]=1.5;
	A[1][2]=2.;
    cout << A[0][2] << endl;
    cout << A[1][2] << endl;
	double* x=new double[2];
	x=partial_pivot(A,2);
	for (int i=0;i<2;i++) {
		cout << x[i] << endl;
	}
}

double* complete_pivot(double** A,int size) {
    int i, j, p, k, q, tmp;
    int* NROW=new int[size];
    int* NCOL=new int[size];
    double** m=new double*[size];
    for (i=0;i<size;i++) m[i]=new double[size+1];
    double* x=new double[size];
    double sigma;
    double A_pq;
    double A_kj;
    bool swap_row;
    bool swap_col;
    
    for (i=0;i<size;i++) {                          //                                      //    
        NROW[i]=i;                                  // initialize row pointer               //
        NCOL[i]=i;                                  // initialize column pointer            //
    }
    for (i=0;i<size;i++) {
        p=i;
        q=i;
        swap_row=false;
        swap_col=false;
        A_pq=abs(A[NROW[p]][NCOL[q]]);
        for (k=i;k<size;k++) {
            for (j=k;j<size;j++) {
                A_kj=abs(A[NROW[k]][NCOL[j]]);
                if (A_kj>A_pq) {
                    if (k!=p) {
                        p=k;
                        swap_row=true;
                    }
                    if (j!=q) {
                        q=j;
                        swap_col=true;
                    }
                    A_pq=abs(A[NROW[p]][NCOL[q]]);
                }
            }
        }
        if (swap_row==true) {
            tmp=NROW[i];
            NROW[i]=NROW[p];
            NROW[p]=tmp;
        }
        if (swap_col==true) {
            tmp=NCOL[i];
            NCOL[i]=NCOL[q];
            NCOL[q]=tmp;
        }
        for (j=i+1;j<size;j++) {
            m[NROW[j]][NCOL[i]]=A[NROW[j]][NCOL[i]]/A[NROW[i]][NCOL[i]];
            for (tmp=0;tmp<size+1;tmp++) {
                A[NROW[j]][tmp]=A[NROW[j]][tmp]-m[NROW[j]][NCOL[i]]*A[NROW[i]][tmp];
            }
        }
    }
    x[NCOL[size-1]]=A[NROW[size-1]][size]/A[NROW[size-1]][NCOL[size-1]];
    for (i=size-2;i>=0;i--) {
        sigma=0.;
        for (j=i;j<size;j++) {
            sigma+=A[NROW[i]][NCOL[j]]*x[NCOL[j]];
        }
        x[NCOL[i]]=(A[NROW[i]][size]-sigma)/A[NROW[i]][NCOL[i]];
    }
    return x;
}

double* partial_pivot(double** A,int size) {
    int i, j, p, k, q, tmp;
    int* NROW=new int[size];
    double** m=new double*[size+1];
    for (i=0;i<size;i++) m[i]=new double[size+1];
    double* x=new double[size];
    double sigma;
    bool swap;
    
    for (i=0;i<size;i++) NROW[i]=i;
    for (j=0;j<size;j++) {
        p=j;
        swap=false;
        for (i=j;i<size;i++) {
            for (j=k;j<size;j++) {
                if (abs(A[i][j])>abs(A[p][j])) {
                    p=i;
                    swap=true;
                }
                if (swap==true) {
                    tmp=NROW[i];
                    NROW[i]=NROW[p];
                    NROW[p]=tmp;
                }  
            }
        }
        for (i=j+1;i<size;i++) {
            m[NROW[i]][j]=A[NROW[i]][j]/A[NROW[j]][j];
            for (tmp=0;tmp<=size;tmp++) {
                A[NROW[i]][tmp]=A[NROW[i]][tmp]-m[NROW[i]][j]*A[NROW[j]][tmp];
            }
        }
    }
    x[size-1]=A[NROW[size-1]][size]/A[NROW[size-1]][size-1];
    for (i=size-2;i>=0;i--) {
        sigma=0.;
        for (j=i+1;j<size;j++) sigma+=A[NROW[i]][j]*x[j];
        x[i]=(A[NROW[i]][size]-sigma)/A[NROW[i]][i];
    }
    return x;
}
