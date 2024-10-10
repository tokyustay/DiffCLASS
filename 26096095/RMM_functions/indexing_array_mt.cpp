/*==========================================================
 * indexing_array.c
 * put all elements of the input array 'Tkk' into the output array 'Tdkk' 
 * according to the subscript positions referred by 'index' array.
 * 
 * usage Tdkk = indexing_array(Tkk, index, n_dk, n_kin)
 * input: 
 *      T_kk: input array or matrix
 *      index: 1D column vector which contains the subscripts of Tdkk
 * output:
 *      Tdkk
 *========================================================*/
#include <mex.h> 
#include "matrix.h"
#include <math.h>
#include <thread>
#include <vector>

using std::thread;
using std::vector;

void errorChk(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void worker(mwSize, mwSize, mxUint32*, mxComplexSingle*, mxComplexSingle*);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   errorChk(nlhs, plhs, nrhs, prhs);

    mxUint32 *index;
    mxComplexSingle *Tkk, *Tdkk;    
    mwSize num_of_index, n_dk, m_dk, len, end_tmp, i;

    mwSize num_of_workers = (mwSize) std::thread::hardware_concurrency()/2;
    vector<thread> workers;
    vector<mwSize> start_index;
    vector<mwSize> end_index;   

    Tkk = mxGetComplexSingles(prhs[0]);  /* pointer for input Tkk */
    
    index = mxGetUint32s(prhs[1]);                   /* index array */
	num_of_index = mxGetM(prhs[1])*mxGetN(prhs[1]);  /* length of index array */

    m_dk = (mwSize)mxGetScalar(prhs[2]);             /* number of rows in Tdkk */
    n_dk = (mwSize)mxGetScalar(prhs[3]);             /* number of columns in Tdkk */
    plhs[0] = mxCreateNumericMatrix(n_dk, m_dk, mxSINGLE_CLASS, mxCOMPLEX); /* Create an nr_dk-by-nc_k mxArray  */
    Tdkk = mxGetComplexSingles(plhs[0]);
    
    len = num_of_index/num_of_workers;
    for (i = 0; i < (num_of_workers-1); i++)
    {
        start_index.push_back(i*len);
        end_index.push_back((i+1)*len-1);
    }
    start_index.push_back((num_of_workers-1)*len);
    end_index.push_back(num_of_index-1);
        
    for (mwSize i = 0; i < num_of_workers; i++)
        workers.push_back(thread(worker, start_index[i], end_index[i], index, Tkk, Tdkk));

    for (mwSize i = 0; i < num_of_workers; i++)
        workers[i].join();
}

void errorChk(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mwSize n, m, n_index, m_index;
    
    /* check for the proper number of arguments */
    if(nrhs != 4)
        mexErrMsgIdAndTxt( "MATLAB:indexing_array:invalidNumInputs", "Four inputs required.");
    
    if(nlhs > 1)
        mexErrMsgIdAndTxt( "MATLAB:indexing_array:maxlhs", "Too many output arguments.");
    
    if( !mxIsSingle(prhs[0]) )
        mexErrMsgIdAndTxt( "MATLAB:indexing_array:inputNotFloatingPointNumbers", "Input must be a single precision array.");
    
    n = mxGetN(prhs[0]);                 /* number of columns of Tkk */
	m = mxGetM(prhs[0]);                 /* number of rows of Tkk */
    n_index = mxGetN(prhs[1]);           /* length of index array */
    m_index = mxGetM(prhs[1]);           /* length of index array */
    
    if( !mxIsUint32(prhs[1]) || n_index!=1 )
        mexErrMsgIdAndTxt( "MATLAB:indexing_array:inputNot1DInt64Array", "index input must be a 1D uint32 column vector.");
        
    if( (n*m) != m_index )
        mexErrMsgIdAndTxt( "MATLAB:indexing_array:inputIndexSizeDiff", "number of elements in index array must be the same as the number of Tkk elements.");    
}

void worker(mwSize start, mwSize end, mxUint32* index, mxComplexSingle* Tkk, mxComplexSingle* Tdkk) 
{
    for (mwSize i = start; i <= end; i++)
        Tdkk[index[i]-1] = Tkk[i];
}
