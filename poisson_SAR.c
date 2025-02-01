// This code is largely derived from original C code from Sven Hyberts around 2011. 
// This version (3.0) is by Scott Anthony Robson Nov 2013. This one will shuffle.
// Arguments are as follows
// 0) name of program (poisson)
// 1) number of NUS dimensions
// 2) seed num (0 means seed based on execution time) 
// 3) sine portion (2, 1 or 0. 1 and 0 are not practical. Only here for completeness and experimental purposes)
// 4) number of sampled points
// 5) tolerance (1 = 100%)
// 6) total size of dimension 1 (the full range)
// 7) total size of dimension 2 (the full range) 
// 8) total size of dimension 3 (the full range) 
// 9) 0 = in order 1 = shuffled
//
// This version constructs outputs in slow -> slower -> slowest order 
// from left to right. Direct dimension is fast and linear so not 
// needed or produced by this code.
 
// Description of algorithm used is found in:
// Applications of non-uniform sampling and processing
// Sven G. Hyberts, Haribabu Arthanari, and Gerhard Wagner
// Top Curr Chem. 2012; 316: 125–148.


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include "poisson_SAR.h"


void shuffle(int *array, size_t n)
{

    struct timeval tv;
    gettimeofday(&tv, NULL);
    int usec = tv.tv_usec;
    srand48(usec);


    if (n > 1)
    {
        size_t i;
        for (i = n - 1; i > 0; i--)
        {
          size_t j = (unsigned int) (drand48()*(i+1));
	  if ( j != 0 ) { // don't shuffle the first point
          	int t = array[j];
          	array[j] = array[i];
          	array[i] = t;
		}
        }
    }
}



int	poisson ( double lmbd )
{
	double	L = exp( -lmbd );
	int	k = 0;
	double	p = 1;

	do {
		double	u = drand48();
		p *= u;
		k += 1;
	} while ( p >= L );

	return( k-1 );
}

int	poisson_gap( int  direction, int3 i_0, int3 i_n, int *v, float ld, float w, float sine_portion )
{
	int	i;
	int	k = 0;
	int	active = 0;

	int3	passive;

	switch	(direction ) {
		case 0:	{	active = i_0[0] ;
				passive[0] = i_0[1];
				passive[1] = i_0[2];
			} ; break;
		case 1:	{	active = i_0[1] ;
				passive[0] = i_0[0];
				passive[1] = i_0[2];
			} ; break;
		case 2:	{	active = i_0[2] ;
				passive[0] = i_0[0];
				passive[1] = i_0[1];
			} ; break;
	}

	while ( active < i_n[direction] ) {

		//  Now make a gap : //
		//printf("Active = %i\n", active);
		if (sine_portion == 0) { active += poisson( (ld-1.0)*w);}
		else {active += poisson( (ld-1.0)*w*sin((float)(active+passive[0]+passive[1])/(float)(i_n[0]+i_n[1]+i_n[2]-3)*M_PI/sine_portion) );}

		if ( active < i_n[direction] ) {

			v[k] = active;    //  store point to be acquired //
			k+= 1;            //  next index of acqured point //

		}
		active += 1;              //  put pointer at next point (gap's end)  //

	}

	return ( k );
}

int	poisson_01_gap( int3 s_0, int3 z, int **v2d, float ld, float w, float sine_portion )
{

	int	i;
	int	ii;
	int	n;
	int	z_min;
	int	z_max;
	int	d_min;

	static	int	*v = NULL;

	int3	s;
	int3	k;

	float3	fss;

//	fprintf( stderr, "poisson_01_gap\n" );

	z_max = z[0] > z[1] ? z[0] : z[1];

	//printf("z0: %i and z1: %i\n", z[0], z[1]);
	if ( !v ) v = (int*) malloc( z_max*sizeof( int ) );

	for( i = 0 ; i < 3 ; i++ ) s[i] = s_0[i];

	for ( k[0] = 0 ; k[0] < z[0] ; k[0]++ ) {
		for ( k[1] = 0 ; k[1] < z[1] ; k[1]++ ) v2d[k[0]][k[1]] = 0;
	}

	z_min = z[0] < z[1] ? z[0] : z[1];
	d_min = z[0] < z[1] ? 0 : 1;

	fss[0] = (float)z[0]/(float)z_min;
	fss[1] = (float)z[1]/(float)z_min;

	while ( (s[0] < z[0]) || (s[1] < z[1]) ) {

		int	times;

		times = ((int) ( (s[d_min]+1)*fss[1] ) - (int) ( (s[d_min]+0)*fss[1] ));
		for ( ii = 0 ; ii < times ; ii++ ) {
			//fprintf( stderr, "yeap\n" );
			if ( s[1] < z[1]) { // do thread along 0 //

				int3	origin;		// origin of the thread //

				for( i = 0 ; i < 3 ; i++ ) origin[i] = s[i];

				if ( origin[0] < z[0] ) {

					while( origin[0] > 0 && !v2d[origin[0]][origin[1]] ) origin[0] -=1;

					n = poisson_gap( 0, origin, z, v, ld, w, sine_portion );

					for ( i = origin[0] ; i < z[0] ; i++ ) v2d[i][origin[1]] = 0;
					for ( i = 0 ; i < n ; i++ ) v2d[v[i]][origin[1]] = 1;

				}

				s[1] += 1;  // increment orthogonal dimension //

			}

		}

		times = ((int) ( (s[d_min]+0)*fss[0] ) - (int) ( (s[d_min]-1)*fss[0] ));
		for ( ii = 0 ; ii < times ; ii++ ) {

			if ( s[0] < z[0]) { // do thread along 1 //
	
				int3	origin;		// origin of the thread //

				for( i = 0 ; i < 3 ; i++ ) origin[i] = s[i];

				if ( origin[1] < z[1] ) {

					while( origin[1] > 0 && !v2d[origin[0]][origin[1]] ) origin[1] -=1;

					n = poisson_gap( 1, origin, z, v, ld, w, sine_portion );

					for ( i = origin[1] ; i < z[1] ; i++ ) v2d[origin[0]][i] = 0;
					for ( i = 0 ; i < n ; i++ ) v2d[origin[0]][v[i]] = 1;

				}

				s[0] += 1;  // increment orthogonal dimension //

			}

		}

		times = ((int) ( (s[d_min]+1)*fss[0] ) - (int) ( (s[d_min]+0)*fss[0] ));
		for ( ii = 0 ; ii < times ; ii++ ) {

			if ( s[0] < z[0]) { // do thread along 1 //
	
				int3	origin;		// origin of the thread //

				for( i = 0 ; i < 3 ; i++ ) origin[i] = s[i];

				if ( origin[1] < z[1] ) {

					while( origin[1] > 0 && !v2d[origin[0]][origin[1]] ) origin[1] -=1;

					n = poisson_gap( 1, origin, z, v, ld, w, sine_portion );

					for ( i = origin[1] ; i < z[1] ; i++ ) v2d[origin[0]][i] = 0;
					for ( i = 0 ; i < n ; i++ ) v2d[origin[0]][v[i]] = 1;

				}

				s[0] += 1;  // increment orthogonal dimension //

			}

		}

		times = ((int) ( (s[d_min]+0)*fss[1] ) - (int) ( (s[d_min]-1)*fss[1] ));
		for ( ii = 0 ; ii < times ; ii++ ) {

			if ( s[1] < z[1]) { // do thread along 0 //

				int3	origin;		// origin of the thread //

				for( i = 0 ; i < 3 ; i++ ) origin[i] = s[i];

				if ( origin[0] < z[0] ) {

					while( origin[0] > 0 && !v2d[origin[0]][origin[1]] ) origin[0] -=1;

					n = poisson_gap( 0, origin, z, v, ld, w, sine_portion );

					for ( i = origin[0] ; i < z[0] ; i++ ) v2d[i][origin[1]] = 0;
					for ( i = 0 ; i < n ; i++ ) v2d[v[i]][origin[1]] = 1;

				}

				s[1] += 1;  // increment orthogonal dimension //

			}

		}
	}
	n = 0;

	for ( k[0] = 0 ; k[0] < z[0] ; k[0]++ ) {
		for ( k[1] = 0 ; k[1] < z[1] ; k[1]++ ) if (v2d[k[0]][k[1]] == 1) n += 1;
	}

	return( n );
}

int	poisson_12_gap( int3 s_0, int3 z, int **v2d, float ld, float w, float sine_portion )
{

	int	i;
	int	ii;
	int	n;
	int	z_min;
	int	z_max;
	int	d_min;

	static	int	*v = NULL;

	int3	s;
	int3	k;

	float3	fss;

	//fprintf( stderr, "poisson_12_gap\n" );

	z_max = z[1] > z[2] ? z[1] : z[2];

	if ( !v ) v = (int*) malloc( z_max*sizeof( int ) );

	for( i = 0 ; i < 3 ; i++ ) s[i] = s_0[i];

	for ( k[1] = 0 ; k[1] < z[1] ; k[1]++ ) {
		for ( k[2] = 0 ; k[2] < z[2] ; k[2]++ ) v2d[k[1]][k[2]] = 0;
	}

	z_min = z[1] < z[2] ? z[1] : z[2];
	d_min = z[1] < z[2] ? 1 : 2;

	fss[1] = (float)z[1]/(float)z_min;
	fss[2] = (float)z[2]/(float)z_min;

	//fprintf( stderr, "%d, %d, %f, %f\n", z[1], z[2], fss[1], fss[2] );

	//fprintf( stderr, "0\n" );

	while ( (s[1] < z[1]) || (s[2] < z[2]) ) {

		int	times;


		times = ((int) ( (s[d_min]+1)*fss[2] ) - (int) ( (s[d_min]+0)*fss[2] ));

	//fprintf( stderr, "1\n" );
	//fprintf( stderr, "times = %d\n", times );

		for ( ii = 0 ; ii < times ; ii++ ) {

			if ( s[2] < z[2]) { // do thread along 1 //

				int3	origin;		// origin of the thread //

				for( i = 0 ; i < 3 ; i++ ) origin[i] = s[i];

				if ( origin[1] < z[1] ) {

					while( origin[1] > 0 && !v2d[origin[1]][origin[2]] ) origin[1] -=1;

					n = poisson_gap( 1, origin, z, v, ld, w, sine_portion );
	
					for ( i = origin[1] ; i < z[1] ; i++ ) v2d[i][origin[2]] = 0;
	
					for ( i = 0 ; i < n ; i++ )  v2d[v[i]][origin[2]] = 1;


				}

				s[2] += 1;  // increment orthogonal dimension //

			}

		}

		times = ((int) ( (s[d_min]+0)*fss[1] ) - (int) ( (s[d_min]-1)*fss[1] ));

	//fprintf( stderr, "2\n" );
	//fprintf( stderr, "times = %d\n", times );

		for ( ii = 0 ; ii < times ; ii++ ) {

			if ( s[1] < z[1]) { // do thread along 2 //
	
				int3	origin;		// origin of the thread //

				for( i = 0 ; i < 3 ; i++ ) origin[i] = s[i];

				if ( origin[2] < z[2] ) {

					while( origin[2] > 0 && !v2d[origin[1]][origin[2]] ) origin[2] -=1;

					n = poisson_gap( 2, origin, z, v, ld, w, sine_portion );

					for ( i = origin[2] ; i < z[2] ; i++ ) v2d[origin[1]][i] = 0;
					for ( i = 0 ; i < n ; i++ ) v2d[origin[1]][v[i]] = 1;

				}

				s[1] += 1;  // increment orthogonal dimension //

			}

		}


		times = ((int) ( (s[d_min]+1)*fss[1] ) - (int) ( (s[d_min]+0)*fss[1] ));

	//fprintf( stderr, "3\n" );
	//fprintf( stderr, "times = %d\n", times );

		for ( ii = 0 ; ii < times ; ii++ ) {

			if ( s[1] < z[1]) { // do thread along 2 //
	
				int3	origin;		// origin of the thread //

				for( i = 0 ; i < 3 ; i++ ) origin[i] = s[i];

				if ( origin[2] < z[2] ) {

					while( origin[2] > 0 && !v2d[origin[1]][origin[2]] ) origin[2] -=1;

					n = poisson_gap( 2, origin, z, v, ld, w, sine_portion );

					for ( i = origin[2] ; i < z[2] ; i++ ) v2d[origin[1]][i] = 0;
					for ( i = 0 ; i < n ; i++ ) v2d[origin[1]][v[i]] = 1;

				}

				s[1] += 1;  // increment orthogonal dimension //

			}

		}

		times = ((int) ( (s[d_min]+0)*fss[2] ) - (int) ( (s[d_min]-1)*fss[2] ));

	//fprintf( stderr, "4\n" );
	//fprintf( stderr, "times = %d\n", times );

		for ( ii = 0 ; ii < times ; ii++ ) {

			if ( s[2] < z[2]) { // do thread along 1 //

				int3	origin;		// origin of the thread //

				for( i = 0 ; i < 3 ; i++ ) origin[i] = s[i];

				if ( origin[1] < z[1] ) {

					while( origin[1] > 0 && !v2d[origin[1]][origin[2]] ) origin[1] -=1;

					n = poisson_gap( 1, origin, z, v, ld, w, sine_portion );

					for ( i = origin[1] ; i < z[1] ; i++ ) v2d[i][origin[2]] = 0;
					for ( i = 0 ; i < n ; i++ ) v2d[v[i]][origin[2]] = 1;

				}

				s[2] += 1;  // increment orthogonal dimension //

			}

		}
	}

	//fprintf( stderr, "5\n" );

	n = 0;

	for ( k[1] = 0 ; k[1] < z[1] ; k[1]++ ) {
		for ( k[2] = 0 ; k[2] < z[2] ; k[2]++ ) if (v2d[k[1]][k[2]] == 1) n += 1;
	}

	return( n );
}

int	poisson_20_gap( int3 s_0, int3 z, int **v2d, float ld, float w, float sine_portion )
{

	int	i;
	int	ii;
	int	n;
	int	z_min;
	int	d_min;
	int	z_max;

	static	int	*v;

	int3	s;
	int3	k;

	float3	fss;

	//fprintf( stderr, "poisson_20_gap\n" );

	z_max = z[2] > z[0] ? z[2] : z[0];

	if ( !v ) v = (int*) malloc( z_max*sizeof( int ) );

	for( i = 0 ; i < 3 ; i++ ) s[i] = s_0[i];

	for ( k[2] = 0 ; k[2] < z[2] ; k[2]++ ) {
		for ( k[0] = 0 ; k[0] < z[0] ; k[0]++ ) v2d[k[2]][k[0]] = 0;
	}

	d_min = z[2] < z[0] ? 2 : 0;
	z_min = z[2] < z[0] ? z[2] : z[0];

	fss[2] = (float)z[2]/(float)z_min;
	fss[0] = (float)z[0]/(float)z_min;

	while ( (s[2] < z[2]) || (s[0] < z[0]) ) {

		int	times;

		times = ((int) ( (s[d_min]+1)*fss[0] ) - (int) ( (s[d_min]+0)*fss[0] ));

		for ( ii = 0 ; ii < times ; ii++ ) {

			if ( s[0] < z[0]) { // do thread along 2 //

				int3	origin;		// origin of the thread //

				for( i = 0 ; i < 3 ; i++ ) origin[i] = s[i];

				if ( origin[2] < z[2] ) {

					while( origin[2] > 0 && !v2d[origin[2]][origin[0]] ) origin[2] -=1;

					n = poisson_gap( 2, origin, z, v, ld, w, sine_portion );

					for ( i = origin[2] ; i < z[2] ; i++ ) v2d[i][origin[0]] = 0;
					for ( i = 0 ; i < n ; i++ ) v2d[v[i]][origin[0]] = 1;

				}

				s[0] += 1;  // increment orthogonal dimension //

			}

		}

		times = ((int) ( (s[d_min]+0)*fss[2] ) - (int) ( (s[d_min]-1)*fss[2] ));

		for ( ii = 0 ; ii < times ; ii++ ) {

			if ( s[2] < z[2]) { // do thread along 0 //
	
				int3	origin;		// origin of the thread //

				for( i = 0 ; i < 3 ; i++ ) origin[i] = s[i];

				if ( origin[0] < z[0] ) {

					while( origin[0] > 0 && !v2d[origin[2]][origin[0]] ) origin[0] -=1;

					n = poisson_gap( 0, origin, z, v, ld, w, sine_portion );

					for ( i = origin[0] ; i < z[0] ; i++ ) v2d[origin[2]][i] = 0;
					for ( i = 0 ; i < n ; i++ ) v2d[origin[2]][v[i]] = 1;

				}

				s[2] += 1;  // increment orthogonal dimension //

			}

		}

		times = ((int) ( (s[d_min]+1)*fss[2] ) - (int) ( (s[d_min]+0)*fss[2] ));

		for ( ii = 0 ; ii < times ; ii++ ) {

			if ( s[2] < z[2]) { // do thread along 0 //
	
				int3	origin;		// origin of the thread //

				for( i = 0 ; i < 3 ; i++ ) origin[i] = s[i];

				if ( origin[0] < z[0] ) {

					while( origin[0] > 0 && !v2d[origin[2]][origin[0]] ) origin[0] -=1;

					n = poisson_gap( 0, origin, z, v, ld, w, sine_portion );

					for ( i = origin[0] ; i < z[0] ; i++ ) v2d[origin[2]][i] = 0;
					for ( i = 0 ; i < n ; i++ ) v2d[origin[2]][v[i]] = 1;

				}

				s[2] += 1;  // increment orthogonal dimension //

			}

		}

		times = ((int) ( (s[d_min]+0)*fss[0] ) - (int) ( (s[d_min]-1)*fss[0] ));

		for ( ii = 0 ; ii < times ; ii++ ) {

			if ( s[0] < z[0]) { // do thread along 2 //

				int3	origin;		// origin of the thread //

				for( i = 0 ; i < 3 ; i++ ) origin[i] = s[i];

				if ( origin[2] < z[2] ) {

					while( origin[2] > 0 && !v2d[origin[2]][origin[0]] ) origin[2] -=1;

					n = poisson_gap( 2, origin, z, v, ld, w, sine_portion );

					for ( i = origin[2] ; i < z[2] ; i++ ) v2d[i][origin[0]] = 0;
					for ( i = 0 ; i < n ; i++ ) v2d[v[i]][origin[0]] = 1;

				}

				s[0] += 1;  // increment orthogonal dimension //

			}

		}
	}

	n = 0;

	for ( k[2] = 0 ; k[2] < z[2] ; k[2]++ ) {
		for ( k[0] = 0 ; k[0] < z[0] ; k[0]++ ) if (v2d[k[2]][k[0]] == 1) n += 1;
	}

	return( n );
}

int	poisson_012_gap( int3 s_0, int3 z, int ***v3d, float ld, float w, float sine_portion )
{
	int	i;
	int	ii;
	int	n;
	int	z_min;
	int	d_min;

	static	int	**v2d_01 = NULL;
	static	int	**v2d_12 = NULL;
	static	int	**v2d_20 = NULL;

	int3	s;
	int3	k;

	float3	fss;

	if ( !v2d_01 ) {
		v2d_01 = ( int** ) malloc( z[0]*sizeof(int*) );
		for ( i = 0 ; i < z[0] ; i++ ) v2d_01[i] = ( int* ) malloc ( z[1]*sizeof(int) );
	}

	if ( !v2d_12 ) {
		v2d_12 = ( int** ) malloc( z[1]*sizeof(int*) );
		for ( i = 0 ; i < z[1] ; i++ ) v2d_12[i] = ( int* ) malloc ( z[2]*sizeof(int) );
	}

	if ( !v2d_20 ) {
		v2d_20 = ( int** ) malloc( z[2]*sizeof(int*) );
		for ( i = 0 ; i < z[2] ; i++ ) v2d_20[i] = ( int* ) malloc ( z[0]*sizeof(int) );
	}


	for( i = 0 ; i < 3 ; i++ ) s[i] = s_0[i];

	for ( k[0] = 0 ; k[0] < z[0] ; k[0]++ ) {
		for ( k[1] = 0 ; k[1] < z[1] ; k[1]++ ) {
			for ( k[2] = 0 ; k[2] < z[2] ; k[2]++ )	v3d[k[0]][k[1]][k[2]] = 0;
		}	
	}

	d_min = z[0] < z[1] ? 0 : 1;
	d_min = z[d_min] < z[2] ? d_min: 2;

	z_min = z[0] < z[1] ? z[0] : z[1];
	z_min = z_min< z[2] ? z_min: z[2];

	fss[0] = (float)z[0]/(float)z_min;
	fss[1] = (float)z[1]/(float)z_min;
	fss[2] = (float)z[2]/(float)z_min;

	while ( (s[0] < z[0]) || (s[1] < z[1]) || (s[2] < z[2]) ) {

		int	times;

		times = ((int) ( (s[d_min]+1)*fss[2] ) - (int) ( (s[d_min]+0)*fss[2] ));

		for ( ii = 0 ; ii < times ; ii++ ) {

			if ( s[2] < z[2]) { // do plane along 01 //

				int3	origin;		// origin of the plane //

				for( i = 0 ; i < 3 ; i++ ) origin[i] = s[i];

				poisson_01_gap( origin, z, v2d_01, ld, w, sine_portion );

				for ( k[0] = origin[0] ; k[0] < z[0] ; k[0]++ ) {
					for ( k[1] = origin[1] ; k[1] < z[1] ; k[1]++ ) {
						v3d[k[0]][k[1]][origin[2]] = v2d_01[k[0]][k[1]];
					}
				}

				s[2] += 1;  // increment orthogonal dimension //

			}

		}

		switch (d_min) {
			case 0 : { times = ((int) ( (s[d_min]+1)*fss[0] ) - (int) ( (s[d_min]+0)*fss[0] ));
				 } break;
			case 1 : { times = ((int) ( (s[d_min]+1)*fss[0] ) - (int) ( (s[d_min]+0)*fss[0] ));
				 } break;
			case 2 : { times = ((int) ( (s[d_min]+0)*fss[0] ) - (int) ( (s[d_min]-1)*fss[0] )); 
				 } break;
		}


		for ( ii = 0 ; ii < times ; ii++ ) {

			if ( s[0] < z[0]) { // do plane along 12 //

				int3	origin;		// origin of the plane //

				for( i = 0 ; i < 3 ; i++ ) origin[i] = s[i];

				poisson_12_gap( origin, z, v2d_12, ld, w, sine_portion );

				for ( k[1] = origin[1] ; k[1] < z[1] ; k[1]++ ) {
					for ( k[2] = origin[2] ; k[2] < z[2] ; k[2]++ ) {
						v3d[origin[0]][k[1]][k[2]] = v2d_12[k[1]][k[2]];
					}
				}

				s[0] += 1;  // increment orthogonal dimension //

			}

		}

		times = ((int) ( (s[d_min]+0)*fss[1] ) - (int) ( (s[d_min]-1)*fss[1] ));

		for ( ii = 0 ; ii < times ; ii++ ) {

			if ( s[1] < z[1]) { // do plane along 20 //

				int3	origin;		// origin of the plane //

				for( i = 0 ; i < 3 ; i++ ) origin[i] = s[i];

				poisson_20_gap( origin, z, v2d_20, ld, w, sine_portion );

				for ( k[2] = origin[2] ; k[2] < z[2] ; k[2]++ ) {
					for ( k[0] = origin[0] ; k[0] < z[0] ; k[0]++ ) {
						v3d[k[0]][origin[1]][k[2]] = v2d_20[k[2]][k[0]];
					}
				}

				s[1] += 1;  // increment orthogonal dimension //

			}

		}


		times = ((int) ( (s[d_min]+1)*fss[1] ) - (int) ( (s[d_min]+0)*fss[1] ));

		for ( ii = 0 ; ii < times ; ii++ ) {

			if ( s[1] < z[1]) { // do plane along 20 //

				int3	origin;		// origin of the plane //

				for( i = 0 ; i < 3 ; i++ ) origin[i] = s[i];

				poisson_20_gap( origin, z, v2d_20, ld, w, sine_portion );

				for ( k[2] = origin[2] ; k[2] < z[2] ; k[2]++ ) {
					for ( k[0] = origin[0] ; k[0] < z[0] ; k[0]++ ) {
						v3d[k[0]][origin[1]][k[2]] = v2d_20[k[2]][k[0]];
					}
				}

				s[1] += 1;  // increment orthogonal dimension //

			}

		}

		switch (d_min) {
			case 0 : { times = ((int) ( (s[d_min]+1)*fss[0] ) - (int) ( (s[d_min]+0)*fss[0] )); 
				 } break;
			case 1 : { times = ((int) ( (s[d_min]+0)*fss[0] ) - (int) ( (s[d_min]-1)*fss[0] )); 
				 } break;
			case 2 : { times = ((int) ( (s[d_min]+1)*fss[0] ) - (int) ( (s[d_min]+0)*fss[0] )); 
				 } break;
		}

		for ( ii = 0 ; ii < times ; ii++ ) {

			if ( s[0] < z[0]) { // do plane along 12 //

				int3	origin;		// origin of the plane //

				for( i = 0 ; i < 3 ; i++ ) origin[i] = s[i];

				poisson_12_gap( origin, z, v2d_12, ld, w, sine_portion );

				for ( k[1] = origin[1] ; k[1] < z[1] ; k[1]++ ) {
					for ( k[2] = origin[2] ; k[2] < z[2] ; k[2]++ ) {
						v3d[origin[0]][k[1]][k[2]] = v2d_12[k[1]][k[2]];
					}
				}

				s[0] += 1;  // increment orthogonal dimension //

			}

		}

		times = ((int) ( (s[d_min]+0)*fss[2] ) - (int) ( (s[d_min]-1)*fss[2] ));

		for ( ii = 0 ; ii < times ; ii++ ) {

			if ( s[2] < z[2]) { // do plane along 01 //

				int3	origin;		// origin of the plane //

				for( i = 0 ; i < 3 ; i++ ) origin[i] = s[i];

				poisson_01_gap( origin, z, v2d_01, ld, w, sine_portion );

				for ( k[0] = origin[0] ; k[0] < z[0] ; k[0]++ ) {
					for ( k[1] = origin[1] ; k[1] < z[1] ; k[1]++ ) {
						v3d[k[0]][k[1]][origin[2]] = v2d_01[k[0]][k[1]];
					}
				}

				s[2] += 1;  // increment orthogonal dimension //

			}

		}


	}

	n = 0;

	for ( k[0] = 0 ; k[0] < z[0] ; k[0]++ ) {
		for ( k[1] = 0 ; k[1] < z[1] ; k[1]++ ) {
			for ( k[2] = 0 ; k[2] < z[2] ; k[2]++ )	if ( v3d[k[0]][k[1]][k[2]] == 1 ) n += 1;
		}	
	}


	return( n );

}

int	main_1d( int argc, char** argv )
{
	
	int3	i_0;
	int3	z;
	float	seed = atof( argv[2] );  //  initial seed     //
	int	p = atoi( argv[4] );  //  sampling points  //
	float	sine_portion = atof( argv[3] );  //  sine portion       //
	int	*v;                   //  vector of sampling points //
	float	ld;
	float	w = 2.0;              //  inital weight    //
	int	k = 0;                //  currently found sampling points //
	int     kk = 0; 	      //  number of sampled points final //
	float   tol = atof( argv[5] ); // tolerance 1 = 100%, 0.01 = 1%
	if (tol==0) {tol = 0.000001;}
	if (seed != 0 ) {
		srand48( seed );        	      //  initialize seed //
		}
	else {
		 srand48(time(NULL));
		}

	z[0] = atoi( argv[6] );       //  total size       //
	z[1] = 1;
	z[2] = 1;

	v = ( int* ) malloc( z[0]*sizeof( int ) );

	ld = (float)z[0]/(float)p;


	do {
		
		i_0[0] = 0;                //  N.B. first point always acqured //
		i_0[1] = 0;                //  N.B. first point always acqured //
		i_0[2] = 0;                //  N.B. first point always acqured //
                                      //  we use the nomenclature of first point == 0 //

		k = poisson_gap( 0, i_0, z, v, ld, w, sine_portion );

//  if more points T.B. acquired found than than wanted: try again with weight 2% larger //
// if fewer points T.B. acquired found than than wanted: try again with weight 2% smaller //

		if ( (k <= p*(1-tol)) || (k >= p*(1+tol)) ) w *= (1.0 + 0.5*(k-p)/p);

	} while ( (k <= p*(1-tol)) || (k >= p*(1+tol)) ); // try until correct number points to be acquired found //

	char a[k][32];
	char point[32];
	int index=0;
//  print the data on standard output //
/*	printf("# number of points = %d\n", k);i*/
	for ( kk = 0 ; kk < k ; kk++ ) {
		/*printf( "%d\n", v[kk] );*/
		sprintf(point, "%4d\n", v[kk]);
		strcpy(a[index], point);
		index++;
	}

	int pool[index];	
	int li = 0;
	for (li = 0 ; li < index ; li++ ) {
	        pool[li] = li;
        }

	if (atoi(argv[9]) == 1) {
		shuffle(pool, index);
	}

	int run = 0 ;
	for ( run = 0 ; run < index ; run++ ) {
		int temp = pool[run];
		printf("%s", a[pool[run]]);
//		printf("%i\n", pool[run]);
	}

	return( 0 );
}

int	main_2d( int argc, char** argv )
{

	int	i, k1, k2, n;
	int3	i_0;
	int3	z;
	float   seed = atof( argv[2] );		// input seed value
	int	z1 = atoi( argv[6] );		// input maximum coordinate along 1st dim
	int	z2 = atoi( argv[7] );		// input maximum coordinate along 2nd dim
	int	tn = atoi( argv[4] );		// input total number of data points in schedule  < z1 * z2
	float	sine_portion = atof( argv[3] ); //  sine portion       //
	float   w = 2.0; 	                //  inital weight    //
	int	**v2d;
	float	ld;
	float	tol = atof( argv[5] ); // tolerance 1 = 100%, 0.01 = 1%

	if (tol==0) {tol = 0.000001;}
        if (seed != 0 ) {
                srand48( seed );                      //  initialize seed //
                }
        else {
                 srand48(time(NULL));
                }

	z[0] = z1;
	z[1] = z2;
	z[2] = 0;

	v2d = ( int** ) malloc( z1*sizeof(int*) );
		for ( i = 0 ; i < z1 ; i++ ) v2d[i] = ( int* ) malloc ( z2*sizeof(int) );

	ld = ( (float)z[0]*z[1] / (float) tn );

	do {

		i_0[0] = 0;                //  N.B. first point always acqured //
		i_0[1] = 0;                //  N.B. first point always acqured //
		i_0[2] = 0;                //  N.B. first point always acqured //
		n = poisson_01_gap( i_0, z, v2d, ld, w, sine_portion );

		//printf("%i\n", n);


		if ( (n <= tn*(1-tol)) || (n >= tn*(1+tol)) ) w *= (1.0 + 0.5*(n-tn)/tn);

	} while ( (n <= tn*(1-tol)) || (n >= tn*(1+tol)) ); // try until correct number points to be acquired found // 


	char a[n][32];
	char point[32];
	int index = 0;
	for ( k2 = 0 ; k2 < z2 ; k2++ ) {
		for ( k1 = 0 ; k1 < z1 ; k1++ ) {
			if ( v2d[k1][k2] ) {
		//		printf( "%4d %4d\n", k1, k2); 
				sprintf(point,"%4d %4d\n", k1, k2);
				strcpy(a[index], point);
				index++;
			}
		}
	}


	int pool[index];	
	int li = 0;
	for (li = 0 ; li < index ; li++ ) {
	        pool[li] = li;
        }

	if (atoi(argv[9]) == 1) {
		shuffle(pool, index);
	}
	int run = 0 ;
	for ( run = 0 ; run < index ; run++ ) {
		int temp = pool[run];
		printf("%s", a[pool[run]]);
//		printf("%i\n", pool[run]);
	}

	return( 0 );
}


int	main_3d( int argc, char** argv )
{

	int	i, k, k1, k2, k3, n;
	int3	i_0;
	int3	z;
	float   seed = atof( argv[2] );		// input seed value
	float	w = 1.0;
	int	z1 = atoi( argv[6] );		// input maximum coordinate along 1st dim
	int	z2 = atoi( argv[7] );		// input maximum coordinate along 2nd dim
	int	z3 = atoi( argv[8] );		// input maximum coordinate along 2nd dim
	int	tn = atoi( argv[4] );		// input total number of data points in schedule  < z1 * z2
	float	sine_portion = atof( argv[3] );  //  sine portion       //
	int	***v3d;
	float	ld;
	FILE	*fp;

	float   tol = atof( argv[5] ); // tolerance 1 = 100%, 0.01 = 1%

        if (tol==0) {tol = 0.000001;}
        if (seed != 0 ) {
                srand48( seed );                      //  initialize seed //
                }
        else {
                 srand48(time(NULL));
                }


	z[0] = z1;
	z[1] = z2;
	z[2] = z3;

	v3d = ( int*** ) malloc( z1*sizeof(int**) );

	for ( i = 0 ; i < z1 ; i++ ) v3d[i] = ( int** ) malloc ( z2*sizeof(int*) );

	for ( i = 0 ; i < z1 ; i++ ) for ( k = 0 ; k < z2 ; k++ ) v3d[i][k] = ( int* ) malloc ( z3*sizeof(int) );

	ld = ( (float)z[0]*z[1]*z[2] / (float) tn );

	do {

		i_0[0] = 0;                //  N.B. first point always acqured //
		i_0[1] = 0;                //  N.B. first point always acqured //
		i_0[2] = 0;                //  N.B. first point always acqured //

		n = poisson_012_gap( i_0, z, v3d, ld, w, sine_portion );

		if ( (n <= tn*(1-tol)) || (n >= tn*(1+tol)) ) w *= (1.0 + 0.5*(n-tn)/tn);

	}  while ( (n <= tn*(1-tol)) || (n >= tn*(1+tol)) ); // try until correct number points to be acquired found //;



	char a[n][32];
	char point[32];
	int index = 0;
	for ( k3 = 0 ; k3 < z3 ; k3++ ) {
		for ( k2 = 0 ; k2 < z2 ; k2++ ) {
			for ( k1 = 0 ; k1 < z1 ; k1++ ) {
				if ( v3d[k1][k2][k3] ) {
				//	printf( "%4d %4d %4d\n", k1, k2, k3);
					sprintf(point, "%4d %4d %4d\n", k1, k2, k3);
					strcpy(a[index], point);
					index++;
				}
			}
		}
	}

	int pool[index];
	int li;
	for (li = 0; li < index; li++) {
		pool[li] = li ;
	}
	if (atoi(argv[9]) == 1) {
		shuffle(pool, index);
	}
	int run;
	for (run = 0; run < index; run++) {
		printf("%s", a[pool[run]]);
	}

	return( 0 );

}


int	main( int argc, char** argv )
{
	int	status;

	if ( argc != 10 ) {
		fprintf( stderr, "Wrong number of arguments (%d provided, 9 required).\n\n", argc-1);
		
		fprintf( stderr, "Expected arguments:\n");
		fprintf( stderr, "1) number of NUS dimensions (1, 2, or 3)\n");
		fprintf( stderr, "2) seed num (0 means seed based on execution time)\n");
		fprintf( stderr, "3) sine portion (2, 1 or 0. 1 and 0 are not practical)\n");
		fprintf( stderr, "4) number of sampled points\n");
		fprintf( stderr, "5) tolerance (1 = 100%%)\n");
		fprintf( stderr, "6) total size of dimension 1 (the full range)\n");
		fprintf( stderr, "7) total size of dimension 2 (the full range)\n");
		fprintf( stderr, "8) total size of dimension 3 (the full range)\n");
		fprintf( stderr, "9) 0 = in order, 1 = shuffled\n\n");
		
		fprintf( stderr, "Received arguments:\n");
		fprintf( stderr, "0) %s (program name)\n", argv[0]);
		for (int i = 1; i < argc; i++) {
			fprintf( stderr, "%d) %s\n", i, argv[i]);
		}
		exit( -1 );
	}
	else {
		int numdim = atoi(argv[1]);

		switch ( numdim ) {
			case 1 : { status = main_1d( argc, argv ); break; }
			case 2 : { status = main_2d( argc, argv ); break; }
			case 3 : { status = main_3d( argc, argv ); break; }
			default: { fprintf( stderr, "Must make 1, 2 or 3 poisson gap dimensions\n" );
					exit( -1 );
			 	}
		}
	}


	exit( 0 );
}
