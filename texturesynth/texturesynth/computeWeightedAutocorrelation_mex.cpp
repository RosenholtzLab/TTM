#include "mex.h" 
#include "math.h"
#include <assert.h>

#define min(X, Y)  ((X) < (Y) ? (X) : (Y))
#define max(X, Y)  ((X) < (Y) ? (Y) : (X))
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

bool validIndx(int x, int y, int wt, int ht) {
    if (x < 0)
        return false;
    if (x >= wt)
        return false;
    if (y < 0)
        return false;
    if (y >= ht)
        return false;
    return true;
}
int replicate(int d, int maxd) {
    if (d < 0)
        return 0;
    if (d >= maxd)
        return maxd-1;
    return d;
}
int reflect(int d, int maxd) {
    if (d < 0)
        return -d;
    if (d >= maxd)
        return maxd-(d-maxd+1);
    return d;
}
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
   
    
    if (nrhs != 4) {
        mexErrMsgTxt("Expected one input image, and Na, and mask (same size as image), and (0 or [default=1] for whether to normalize)!");
    }
    if ((nlhs > 2) || (nlhs < 1)){
        printf("nlhs = %d", nlhs);
        mexErrMsgTxt("Expected one or two output images!");
    }    
    // number of coefficients
    int Na =(int) mxGetScalar(prhs[1]);
    
    
    bool normalizeCorrelation = mxGetScalar(prhs[3]) > 0;
    
    
    /* Find the dimensions of the data */
    int ht = mxGetM(prhs[0]);
    int wt = mxGetN(prhs[0]);
    
    if ((ht != mxGetM(prhs[2])) || (wt != mxGetN(prhs[2]))) {
        mexErrMsgTxt("Image and mask need to have the same size");
    }
    
    
    
    /* Create an mxArray for the output data */
    plhs[0] = mxCreateDoubleMatrix(Na, Na, mxREAL);
    
    double* normFactMx = NULL;
    if (nlhs > 1) {
        plhs[1] = mxCreateDoubleMatrix(Na, Na, mxREAL);
        normFactMx = mxGetPr(plhs[1]);
    }
    
    
    // Data stuffs
    double* x = mxGetPr(prhs[0]);
    double* mask = mxGetPr(prhs[2]);
    double* corr = mxGetPr(plhs[0]);
    int K = (int)((Na*Na)/2) + 1; // index of the middle guy
    int half = (int) (Na/2);      // half of Na
    
    
    // make sure corr is all zero
    for (int ii = 0; ii < Na*Na; ii++) {
        corr[ii] = 0;
        if (nlhs > 1) {
            normFactMx[ii] = 0;
        }
    }
    

    //const double numPixels = ((ht-half) - (half+1) + 1) * ((wt-half) - (half+1) + 1);
    
    // compute the mean
    double mean = 0;
    double sumMask = 0;
    for (int ii = 0; ii < ht; ii++) {
    for (int jj = 0; jj < wt; jj++) {
        int indx_ij = jj * ht + ii; 
        mean += x[indx_ij]*mask[indx_ij];
        sumMask += mask[indx_ij];
    }
    }
    mean = mean/sumMask;

    
//     // subtract the mean
//     for (int ii = 0; ii < ht; ii++) {
//     for (int jj = 0; jj < wt; jj++) {
//         int indx_ij = jj * ht + ii; 
//         x[indx_ij] -= mean;
//     }
//     }
    double meansq = mean*mean;
    for (int direction = 0; direction < K; direction++) {
        int n = direction % Na - half;
        int m = direction / Na - half;
        
        int indx_nm = ( m+half) * Na + ( n+half);
        int indx_mn = (-m+half) * Na + (-n+half);
        
        //for (int ii = half; ii < ht-half; ii++) {
        //for (int jj = half; jj < wt-half; jj++) {
        double normfact = 0;
        for (int ii = 0; ii < ht; ii++) {
        for (int jj = 0; jj < wt; jj++) {
            int indx_ij = jj * ht + ii; 
            int iiCirc = (ii+n+ht)%ht;
            int jjCirc = (jj+m+wt)%wt;
            //int indx_ijdir = (jj+m) * ht + (ii+n);                                   
            int indx_ijdir = jjCirc * ht + iiCirc;            
            
            double mv1 = mask[indx_ij];
            double mv2 = mask[indx_ijdir];            
            
            // only add if its non zero
            //if ((mv1 != 0) && (mv2 != 0)) { 
            double nf = mv1*mv2;
            double v1 = x[indx_ij]-mean;
            double v2 = x[indx_ijdir]-mean;
            double val = v1*v2;
            corr[indx_nm] += val*nf;
            normfact += nf;
            //}
        }
        }
        
        
        //printf("NormFact for %d,%d = %f\n", n,m,normfact);
        // normalize
        if (normalizeCorrelation) {
            corr[indx_nm] = corr[indx_nm] / normfact;
        } else {
            //corr[indx_nm] = corr[indx_nm];
        }
        
        // be symmetric
        corr[indx_mn] = corr[indx_nm];
        
        if (nlhs > 1) {
            normFactMx[indx_nm] = normfact;
            normFactMx[indx_mn] = normfact;            
        }
        
    }   
    //printf("sum mask = %f\n", sumMask);
} 