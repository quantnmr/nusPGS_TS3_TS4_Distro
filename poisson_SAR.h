// Header file for poisson_v4_2d.c //

typedef	int	int3[3];

typedef	float	float3[3];

// input: lamda of the poission distribution; return: a poission random number

int	poisson( double );


// input: direction (dimension), i_0 (init coordinates), i_n (size of 3D matrix),
// v (1d vector of poisson gap sampling, *updated*), ld (lamda), w (weight),
// sine_portion (weight for sine function)
// return: number of sampled points

int	poisson_gap( int, int3, int3, int*, float, float, float );

// input: i_0 (init coordinates), i_n (size of 3D matrix),
// v2d (2d vector of poisson gap sampling, *updated*), ld (lamda), w (weight),
// sine_portion (weight for sine function)
// return: number of sampled points

int	poisson_01_gap( int3, int3, int**, float, float, float );

int	poisson_12_gap( int3, int3, int**, float, float, float );

int	poisson_20_gap( int3, int3, int**, float, float, float );

// input: i_0 (init coordinates), i_n (size of 3D matrix),
// v3d (3d vector of poisson gap sampling, *updated*), ld (lamda), w (weight),
// sine_portion (weight for sine function)
// return: number of sampled points

int	poisson_012_gap( int3, int3, int***, float, float, float );

// input: argc (number of argv variables), argv (array of strings)
// return: nothing

int	main_1d( int, char** );

int	main_2d( int, char** );

int	main_3d( int, char** );

