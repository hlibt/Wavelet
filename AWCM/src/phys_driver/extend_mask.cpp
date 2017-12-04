#include "../CollocationPoint.hpp"
using namespace std;

void extend_mask(CollocationPoint** collPnt) {
    for (j=1;j<=J;j++) {                                                // 
        int jstar=j-1;                                                  //
        int gridMultplr=pow(2,j-jstar);                                 //
        for (k=0;k<jPnts(jstar);k++) {                                  //
            if (mask[jstar][k]==true) {                                 //
                mask[j][k*gridMultplr]=true;                            //
            }                                                           //
        }                                                               //
    }                                                                   //
    for (j=J-1;j>0;j--) {                                               // 
        for (k=0;k<jPnts(j);k++) {                                      //
            if (k%2==1 && mask[j][k]==true) {                           //
                int leftPnt=-interpPnts+1;  
                int rightPnt=interpPnts;    
                while ( leftPnt < 0 ) {
                    leftPnt++;
                    rightPnt++;
                }
                while ( rightPnt > (jPnts(j-1)-1) ) {
                    leftPnt--;
                    rightPnt--;
                }
                for (int l=leftPnt;l<=rightPnt;l++) {                   //
                    mask[j-1][(k-1)/2+l]=true;                          //
                }                                                       //
            }                                                           //
        }                                                               //
    }                                                                   //
}
