/*    hist3d_fast.c

      Mike Lawrence 2012-11-17

      To compile: (from Matlab prompt)
      >> cd /xchip/cga2/lawrence/cga/trunk/matlab/mike
      >> mex hist3d_fast.c

*/

#include <string.h>
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ]) 
{
  long xlen, ylen, zlen, len;
  long firstx, lastx, firsty, lasty, firstz, lastz, nxvals, nyvals, nzvals;
  int dims[3];
  double *x, *y, *z, *xvals, *yvals, *zvals, *h;
  long i, xi, yi, zi, idx;

  if (nrhs!=3 && nrhs!=9) mexErrMsgTxt("requires either three or nine input arguments: x, y, z, [firstx, lastx, firsty, lasty, firstz, lastz]");
  if (nlhs>1) mexErrMsgTxt("too many outputs requested");

  xlen = mxGetN(prhs[0]); if (xlen==1) xlen = mxGetM(prhs[0]);
  ylen = mxGetN(prhs[1]); if (ylen==1) ylen = mxGetM(prhs[1]);
  zlen = mxGetN(prhs[1]); if (zlen==1) zlen = mxGetM(prhs[2]);
  if (xlen!=ylen || xlen!=zlen) mexErrMsgTxt("error: x,y,z are of different lengths");
  len = xlen;
  x = mxGetPr(prhs[0]);
  y = mxGetPr(prhs[1]);
  z = mxGetPr(prhs[2]);

  if (nrhs==9) {
    /* user supplied ranges */
    firstx = (long)mxGetScalar(prhs[3]);
    lastx = (long)mxGetScalar(prhs[4]);
    firsty = (long)mxGetScalar(prhs[5]);
    lasty = (long)mxGetScalar(prhs[6]);
    firstz = (long)mxGetScalar(prhs[7]);
    lastz = (long)mxGetScalar(prhs[8]);
  } else {
    /* compute ranges as min and max of x,y,z */
    firstx = *x;
    lastx = *x;
    firsty = *y;
    lasty = *y;
    firstz = *z;
    lastz = *z;
    for (i=1;i<len;i++) {
      if ((*(x+i))<firstx) firstx = (*(x+i));
      if ((*(x+i))>lastx) lastx = (*(x+i));
      if ((*(y+i))<firsty) firsty = (*(y+i));
      if ((*(y+i))>lasty) lasty = (*(y+i));
      if ((*(z+i))<firstz) firstz = (*(z+i));
      if ((*(z+i))>lastz) lastz = (*(z+i));
    }
    printf("x range = %d to %d;  y range = %d to %d;  z range = %d to %d\n",firstx,lastx,firsty,lasty,firstz,lastz);
  }

  nxvals = lastx-firstx+1;
  nyvals = lasty-firsty+1;
  nzvals = lastz-firstz+1;

  if (nxvals<1 || nyvals<1 || nzvals<1) mexErrMsgTxt("cannot have firstx>lastx, firsty>lasty, or firstz>lastz");

  dims[0] = nxvals;
  dims[1] = nyvals;
  dims[2] = nzvals;  

  /*    x = rows        */
  /*    y = columns     */
  /*    z = pages       */

  /*  plhs[0] = mxCreateDoubleMatrix(nxvals,nyvals,mxREAL);  */
  plhs[0] = mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL);

  h = mxGetPr(plhs[0]);
  for(i=0;i<len;i++) {
    xi = (long)(*(x+i))-firstx;
    yi = (long)(*(y+i))-firsty;
    zi = (long)(*(z+i))-firstz;
    if (xi<0 || xi>=nxvals || yi<0 || yi>nyvals || zi<0 || zi>nzvals) continue;  /* ignore out-of-range values */
    /*    idx = (nxvals*yi+xi);  */
    idx = (nyvals*nxvals*zi)+(nxvals*yi)+xi;
    (*(h+idx))++;
  }

}


