#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <omp.h>
using namespace std;

//#define f(x,y)  (sin(x))
#define g(x,y) 	(cos(x))

int main()
{
    int N=10;
    double l = 4*M_PI; // длина интервала [0,l]
    double h = l/ (double) (N+1);
    int i,j;
    double f[N+2][N+2], u[N+2][N+2], un[N+2][N+2];
    double d, dm, dmax;
    double temp;
    double eps = 0.0001;
    int iter_count; // счетчик итераций алгоритма
    double start, stop; // время начала и конца параллельной области

    // матрица неоднородности ур. Лапласа
    for (i=1; i<N+1; i++)
        for (j=1; j<N+1; j++)
            f[i][j] = f (i*h,j*h);

    // краевые условия
    for (i=0; i<=N+1; i++)
    {
        u[i][0] = g (i*h,0*h);
        u[i][N+1] = g (i*h,(N+1)*h);
    }
    // краевые условия
    for (j=0; j<=N+1; j++)
    {
        u[0][j] = g (0*h,j*h);
        u[N+1][j] = g ((N+1)*h,j*h);
    }


    iter_count=0;
    start = omp_get_wtime();
    #pragma omp parallel shared(u,un) 
    do
    {
        dmax = 0; // максимальное изменение значений u
        #pragma omp parallel for shared(u,un) firstprivate(dmax) reduction(max: dmax)
        for ( i=1; i<N+1; i++ )
        {
            for ( j=1; j<N+1; j++ )
            {
                temp = u[i][j];
                un[i][j] = 0.25*(u[i-1][j] + u[i+1][j] + \
                                 u[i][j-1] + u[i][j+1] - h*h*f[i][j] );
                d = fabs(temp-un[i][j]);

                if ( dmax < d ) dmax = d;
            }
            printf("i = %d \t dmax = %f \n", i, dmax);
        } 
        printf("master: \t dmax = %f \n", dmax);
        
        for ( i=1; i<N+1; i++ ) // обновление данных
            for ( j=1; j<N+1; j++ )
                u[i][j] = un[i][j];

        iter_count++;
    } while ( dmax > eps ); // конец параллельной области
    stop = omp_get_wtime();


/*    for (i=0; i <=N+1; i++)
    {
        printf("x=%f ", i*h);
        for (j=0;j <=N+1; j++)
            printf("%f ", u[i][j]);
        printf("\n");
    }*/
    
    printf("eps = %f \t iter_count = %d \t time = %f \n", eps, iter_count, \
		stop-start);

    ofstream myfile;
    myfile.open("data.txt");    
    for (i=0; i <=N+1; i++)
    {
        for (j=0;j <=N+1; j++)
            myfile << u[i][j] << " ";
        myfile << endl;
    }
    myfile.close();
}

