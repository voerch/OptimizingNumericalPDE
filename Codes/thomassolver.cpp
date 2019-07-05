// TridiagonalSolver.cpp : Defines the exported functions for the DLL application.

//

 

#include <windows.h> // for performance counter

 

long __stdcall ASimpleFunction()

{

              return 234;

}

 

long __stdcall TridiagonalSolverC(double *a, double *b, double*c, double*d,

       double*results, long dimension)

{

//     __int64 ctr1(0);

//     QueryPerformanceCounter((LARGE_INTEGER*)&ctr1);

       int iCounter(0);

       //forward sweep

       double denominator=1/b[0];

       c[0]=c[0]*denominator;

       d[0]=d[0]*denominator;

       for(iCounter=1;iCounter<dimension;iCounter++)

       {

              denominator=1/(b[iCounter]-c[iCounter-1]*a[iCounter]);

              c[iCounter]=c[iCounter]*denominator;

              d[iCounter]=(d[iCounter]-d[iCounter-1]*a[iCounter])*denominator;

       }

       //backward sweep

       results[dimension-1]=d[dimension-1];

       for(iCounter=dimension-2;iCounter>=0;iCounter--)

       {

              results[iCounter]=d[iCounter]-results[iCounter+1]*c[iCounter];

       }

//     __int64 ctr2(0);

//     QueryPerformanceCounter((LARGE_INTEGER *)&ctr2);

//     return long(ctr2-ctr1);

       return 1;

}

 

long __stdcall QueryPerformanceCounterFrequency()

{

       __int64 freq(0);

       QueryPerformanceFrequency((LARGE_INTEGER *)&freq);

       return long(freq);

}

 

long __stdcall TridiagonalSolverParallelC(double *a, double *b, double*c, double*d,

       double*results, long dimension)

{

       __int64 ctr1(0);

       QueryPerformanceCounter((LARGE_INTEGER*)&ctr1);

       if(dimension%2!=0)

              return -1;//should be even dimension

//     double *newC = new double[dimension]; // no exstra storage really needed

//     double *newD = new double[dimension]; // bin after code is bug free

       int iCounter(0);

       //forward sweep

       double determinant                = b[0]*b[1]-c[0]*a[1];

       double oneoverdeterminant  = 1/determinant;

       double oldc0                      = c[0];

       c[0]                                     *= -c[1]*oneoverdeterminant;

       c[1]                                     *= +b[0]*oneoverdeterminant;

       double oldd0                      = d[0];

       d[0]                                     = (+b[1]*d[0]-oldc0*d[1])*oneoverdeterminant;

       d[1]                                     = (-a[1]*oldd0+b[0]*d[1])*oneoverdeterminant;

       double newB(0), tempD(0);

//     for(iCounter=1;iCounter<dimension/2-1;iCounter++)

       for(int iCount=2;iCount<=dimension/2-1;iCount++)

       {

//            iCounter=iCount-1; iCount=(iCounter+1)

//            2*iCount-2=2*(iCounter+1)-2=2*iCounter

//            2*iCount-1=2*(iCounter+1)-1=2*iCounter+1

//            2*iCount-3=2*(iCounter+1)-1=2*iCounter+1

              newB                       = b[2*iCount-2]-c[2*iCount-3]*a[2*iCount-2];

              determinant                = newB*b[2*iCount-1]-c[2*iCount-2]*a[2*iCount-1];

        oneoverdeterminant = 1/determinant;

              oldc0                      = c[2*iCount-2];

        c[2*iCount-2]             *= -c[2*iCount-1]*oneoverdeterminant;

        c[2*iCount-1]             *= newB*oneoverdeterminant;

              tempD                      = d[2*iCount-2]-a[2*iCount-2]*d[2*iCount-3];

        d[2*iCount-2]             = b[2*iCount-1]*tempD-oldc0*d[2*iCount-1];

        d[2*iCount-2]             *= oneoverdeterminant;

              d[2*iCount-1]        = -a[2*iCount-1]*tempD+newB*d[2*iCount-1];

        d[2*iCount-1]             *= oneoverdeterminant;

       }

    newB                          = b[dimension-2]-c[dimension-3]*a[dimension-2];

    determinant                   = newB*b[dimension-1]-c[dimension-2]*a[dimension-1];

    oneoverdeterminant     = 1/determinant;

    tempD                         = d[dimension-2]-a[dimension-2]*d[dimension-3];

       d[dimension-2]             = b[dimension-1]*tempD-c[dimension-2]*d[dimension-1];

    d[dimension-2]         *= oneoverdeterminant;

       d[dimension-1]             = -a[dimension-1]*tempD+newB*d[dimension-1];

    d[dimension-1]         *= oneoverdeterminant;

       //backward sweep

       results[dimension-1] = d[dimension-1];

    results[dimension-2]   = d[dimension-2];

       for(int iCount=dimension/2-1;iCount>0;iCount--)

       {

              results[2*iCount-2]  = d[2*iCount-2]-c[2*iCount-2]*results[2*iCount];

              results[2*iCount-1]  = d[2*iCount-1]-c[2*iCount-1]*results[2*iCount];

       }

//     delete []newC;

//     delete []newD;

       __int64 ctr2(0);

       QueryPerformanceCounter((LARGE_INTEGER *)&ctr2);

       return long(ctr2-ctr1);

}