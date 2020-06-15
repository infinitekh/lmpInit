#include <stdio.h>
#include <inttypes.h>
#include <stdlib.h>
#include <math.h>
#include "mymath.h"
#include <snapshot.h>
typedef int smallint;
typedef int64_t imageint;
typedef int64_t tagint;
typedef int64_t bigint;

#define MAXSMALLINT INT_MAX
#define MAXTAGINT INT64_MAX
#define MAXBIGINT INT64_MAX

#define TAGINT_FORMAT "%" PRId64
#define BIGINT_FORMAT "%" PRId64

#define BIGINT_FORMAT "%" PRId64
#define VERSION  "version"
#define ATOM_STYLE "hybrid"
#define AllocMem(a, n, t)  a = (t *) malloc ((n) * sizeof (t))
#define IncreaseMem(a, n, t)  a = (t *) \
		realloc (a, ( N + n) * sizeof (t))

#include <common.h>
typedef double real;
typedef struct {real x,y;} VecR2;
typedef struct {real x,y;} d_complex;
typedef struct {int x,y,z;} VecI3;
typedef struct {real m1,m2; } Prop;

typedef struct {
	tagint *id;
	int *type;
	VecR3 *position ;

	// hybrid sphere
	real *diameter;
	real *mass;
	
	//hybrid dipole
	real *charge;
	VecR3 *dipole;

	// periodic cell number
	VecI3 *cell;
	// velocities
  VecR3 *velocity;
	VecR3 *angular_velocity;
} Particle;
static int flag_alloc              = 0 ;
static Particle ptls;
static bigint timestep             = 0;

static VecR3 boxlo, boxhi;
static VecR3 Space                 = {64.0,64.0,64.0};
static long int N                  = 32;
static int Ntypes                  = 2;
static long seed;
static FILE* fp;
static int flag_firstrandom        = 0;
static int flag_make_solvent_error = 0;

real  no_restriction = 0.0;
real  no_charge      = 0.0;
real** ntype_condition;
void write_header()
{
	
  fprintf(fp,"LAMMPS data file via write_data, version %s, "
          "timestep = " BIGINT_FORMAT "\n",
          VERSION,timestep);

	 fprintf(fp,"\n");

  fprintf(fp,BIGINT_FORMAT " atoms\n",N);
  fprintf(fp,"%d atom types\n",Ntypes);

	fprintf(fp,"\n");

  fprintf(fp,"%-1.16e %-1.16e xlo xhi\n",boxlo.x,boxhi.x);
  fprintf(fp,"%-1.16e %-1.16e ylo yhi\n",boxlo.y,boxhi.y);
  fprintf(fp,"%-1.16e %-1.16e zlo zhi\n",boxlo.z,boxhi.z);

}
void write_type_array()
{
    fprintf(fp,"\nMasses\n\n");
    for (int i = 1; i <= Ntypes; i++) fprintf(fp,"%d %g\n",i,ptls.mass[i]);
}

void write_atoms(){
	fprintf(fp,"\nAtoms # %s\n\n",ATOM_STYLE);
	tagint id=(tagint)0;
	for (int i = 0; i < N; i++) {
		//tagint id    = (tagint) ptls.id[i];
		int type     = ptls.type[i];
		VecR3* r     = & ptls.position[i];

		fprintf(fp,TAGINT_FORMAT " %d %-1.16e %-1.16e %-1.16e",
				(tagint) ptls.id[i],(int)  ptls.type[i],
				r->x,r->y,r->z);


		// Sphere
		real diameter           = ptls.diameter[i];
		real mass               = (real) ptls.mass[i];
		fprintf(fp," %-1.16e %-1.16e",diameter, mass);


		// Dipole
		real   q  = ptls.charge[i];
		VecR3* mu = & ptls.dipole[i];
		fprintf(fp," %-1.16e %-1.16e %-1.16e %-1.16e",q, mu->x, mu->y, mu->z);

		// Periodic cell
		VecI3* cell = &ptls.cell[i];
		fprintf(fp," %d %d %d\n",
				cell->x, cell->y, cell->z);
		id= id+1;

	}


}

void write_velocities(){
	fprintf(fp,"\nVelocities\n\n");
	
	for (int i = 0; i < N; i++) {
		tagint id    = (tagint) ptls.id[i];
		VecR3* v = & ptls.angular_velocity[i];
		VecR3* w = & ptls.angular_velocity[i];

		fprintf(fp,TAGINT_FORMAT
      " %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e\n",
			id, v->x, v->y, v->z,
			  w->x, w->y, w->z);
	}
}
void write_data( FILE* file) {
	fp = file;

	write_header();

	write_type_array();

	write_atoms();

	write_velocities();
}
void write_data_stdout() {
	write_data(stdout);
}

void delete_data();
void renew_data( int N2) {
	
	IncreaseMem(ptls.id,              N2, tagint);
	IncreaseMem(ptls.type,            N2, int);
	IncreaseMem(ptls.position,        N2, VecR3);
	IncreaseMem(ptls.diameter,        N2, real);
	IncreaseMem(ptls.mass,            N2, real);
	IncreaseMem(ptls.charge,          N2, real);
	IncreaseMem(ptls.dipole,          N2, VecR3);
	IncreaseMem(ptls.cell,            N2, VecI3);
	IncreaseMem(ptls.velocity,        N2, VecR3);
	IncreaseMem(ptls.angular_velocity, N2, VecR3);
	
}
void new_data( int _N) {
	N=_N;
	if (flag_alloc == 1)
		delete_data(); 
	
	AllocMem(ptls.id,N, tagint);
	AllocMem(ptls.type,N, int);
	AllocMem(ptls.position,N, VecR3);
	AllocMem(ptls.diameter,N, real);
	AllocMem(ptls.mass,N, real);
	AllocMem(ptls.charge,N, real);
	AllocMem(ptls.dipole,N, VecR3);
	AllocMem(ptls.cell,N, VecI3);
	AllocMem(ptls.velocity,N, VecR3);
	AllocMem(ptls.angular_velocity,N, VecR3);
	for (int i=0 ; i<N ; i++) {
		ptls.id[i] = i+1; ptls.type[i] =1;
		ptls.diameter[i] = 2;
		ptls.mass[i] =1;
		ptls.charge[i] = 0.0;

	}
	
	flag_alloc = 1;
}


void delete_data() {
	if (flag_alloc == 1){
		free(&ptls.id);
		free(&ptls.type);
		free(&ptls.position );

		// hybrid sphere
		free(&ptls.diameter);
		free(&ptls.mass);

		//hybrid dipole
		free(&ptls.charge);
		free(&ptls.dipole);

		// periodic cell number
		free(&ptls.cell);
		// velocities
		free(&ptls.velocity);
		// hybrid spherer
		free(&ptls.angular_velocity);
	}
}


void init_random_type_n_if(tagint t, bigint  N2)
{
	
	if ( flag_firstrandom ==0) {
		seed = good_seed();
		flag_firstrandom = 1;
	}
	if (flag_alloc != 1){
		fprintf(stderr, "No data\n" );
		exit(1);
	}
	bigint i,j;
	int retry=0;
	VecR3 r_ij;
	real r2;

	int N1 = N;
	int id_condition=0;
	int flag_retry=0;
	renew_data(N2);
	
	real initial_lowcut=  pow(ntype_condition[0][t],2);
	i = N1;
	while( i<N1+N2 ) {
		// initializing translation parameters.
		ptls.position[i].x = (boxhi.x-boxlo.x) * ran1(&seed) -  boxlo.x;
		ptls.position[i].y = (boxhi.y-boxlo.y) * ran1(&seed) -  boxlo.y;
		ptls.position[i].z = (boxhi.z-boxlo.z) * ran1(&seed) -  boxlo.z;
		/// hard sphere init position
		for(j=0; j<N1; j++)
		{
			VSub(r_ij, ptls.position[i], ptls.position[j] );
			VWrapAll(r_ij);

			r2 = SQ3D( r_ij);
			if( r2 < initial_lowcut ) {
				flag_retry = 1;
				break;
			}
			flag_retry =0;
		}

		if( flag_retry == 1) 
		{
			retry++;
			if( (real)retry> (real)pow(N1+N2,4)  ) {
				fprintf(stderr, "Random initialize errror retry\n ");
				exit(1);
			}
		 	continue;
		}

		ptls.id[i]         = i+1;
		ptls.type[i]       = t;
		ptls.mass[i]       = 1;
		ptls.diameter[i]   = 0;

		VZero(ptls.velocity[i]);
		VZero(ptls.cell[i]);
/* #ifdef DIPOLE
 *  initializing rotation parameters.
 * 		real theta = ran1(&seed) * atan(1)*4;
 * 		real phi = ran1(&seed) * atan(1)*8;
 * 		ptls.dipole[i].x = sin(theta)*cos(phi);
 * 		ptls.dipole[i].y = sin(theta)*sin(phi);
 * 		ptls.dipole[i].z = cos(theta);
 * 
 * 		VZero(ptls.angular_velocity[i]);
 * #endif
 */
		i++;
	} 

	N = N1+N2;

}

void init_pos_random (int ntype, real* ntype_condition)
{
	if ( flag_firstrandom ==0) {
		seed = good_seed();
		flag_firstrandom = 1;
	}
	if (flag_alloc != 1){
		fprintf(stderr, "No data\n" );
		exit(1);
	}
	bigint i,j;
	int retry=0;
	VecR3 r_ij;
	real r2;

	int id_condition=0;
	int flag_retry=0;
	
	real initial_lowcut=  ntype_condition[id_condition];
	for(i=0; i<N; i++) {
		// initializing translation parameters.
		ptls.position[i].x = (boxhi.x-boxlo.x) * ran1(&seed) -  boxlo.x;
		ptls.position[i].y = (boxhi.y-boxlo.y) * ran1(&seed) -  boxlo.y;
		ptls.position[i].z = (boxhi.z-boxlo.z) * ran1(&seed) -  boxlo.z;
		/// hard sphere init position
		for(j=0; j<i; j++)
		{
			VSub(r_ij, ptls.position[i], ptls.position[j] );
			VWrapAll(r_ij);

			r2 = SQ3D( r_ij);
			if( r2 < initial_lowcut ) {
				flag_retry = 1;
				break;
			}
			flag_retry =0;
		}

		if( flag_retry == 1) 
		{
			retry++;
			if( (real)retry> (real)pow(N,4)  ) {
				fprintf(stderr, "Random initialize errror retry\n ");
				exit(1);
			}
		}
		ptls.id[i]         = i+1;
		ptls.type[i]       = 1;
		ptls.mass[i]       = 1;
		ptls.diameter[i]   = 0;

		VZero(ptls.velocity[i]);
		VZero(ptls.cell[i]);
#ifdef DIPOLE
		// initializing rotation parameters.
		real theta = ran1(&seed) * atan(1)*4;
		real phi = ran1(&seed) * atan(1)*8;
		ptls.dipole[i].x = sin(theta)*cos(phi);
		ptls.dipole[i].y = sin(theta)*sin(phi);
		ptls.dipole[i].z = cos(theta);

		VZero(ptls.angular_velocity[i]);
#endif
	} 

	
	printf("#End call rand pos setting!!\n");
}

void init_pos_Pfcc (int _N,real initial_spacing)
{
#ifdef DEBUG
	puts("call init_pos_Pfcc");
#endif
	int i=0,j=0,k=0, l =0;
	real ri,rj,rk;
	int ni = Space.x / initial_spacing;
	real is = Space.x / ni;

	VecR3 fcc_basis_pos[4];
	int fcc_max = 4;
	real move_on = is/sqrt(2);
	real moveoff = 0;
	VSet(fcc_basis_pos[0], moveoff, moveoff, moveoff);
	VSet(fcc_basis_pos[1], move_on, move_on, moveoff);
	VSet(fcc_basis_pos[2], move_on, moveoff, move_on);
	VSet(fcc_basis_pos[3], moveoff, move_on, move_on);

	puts("#time\tpotential\tmean_sq_displacement\tmsd_noninteraction ");

	int pn=0;
	int leave_ptls = _N;
	bigint max_lattice_point = ni*ni*ni*fcc_max;
	bigint leave_lattice_point = max_lattice_point;
	real prob = 1.0*leave_ptls / leave_lattice_point;


	if( max_lattice_point < leave_ptls ){
		printf("Error Max lattice point is lower than leave ptls \n "
				" %7ld  < %7d\n",max_lattice_point,leave_ptls );
		exit(1);
	}

	for(i = 0 ; i<ni; i++) {
		ri = boxlo.x + i * is;
	for ( j = 0; j < ni; j += 1 ) {
		rj = boxlo.y + j * is;
	for ( k = 0; k < ni; k += 1 ) {
		rk = boxlo.z + k * is;
		for ( l = 0; l < fcc_max; l += 1 ) {
			if( leave_ptls==leave_lattice_point || 
					ran1(&seed) < prob ) {

				ptls.position[pn].x = ri +fcc_basis_pos[l].x;
				ptls.position[pn].y = rj +fcc_basis_pos[l].y;
				ptls.position[pn].z = rk +fcc_basis_pos[l].z;
				VZero(ptls.velocity[pn]);
				VZero(ptls.cell[pn]);
//#ifdef DIPOLE
		// initializing rotation parameters.
		real theta = ran1(&seed) * atan(1)*4;
		real phi = ran1(&seed) * atan(1)*8;
		ptls.dipole[pn].x = sin(theta)*cos(phi);
		ptls.dipole[pn].y = sin(theta)*sin(phi);
		ptls.dipole[pn].z = cos(theta);

		VZero(ptls.angular_velocity[pn]);
//#endif
				pn ++;  leave_ptls--;
			}
			leave_lattice_point -- ;
			if (leave_lattice_point == 0 ) break;
			prob = 1.0*leave_ptls / leave_lattice_point;

				}
			}
		}
	}


	if( N!=pn){
		printf("error N!=pn   %ld,%d\n ",N,pn );
		exit(1);
	}
}


void load_from_file (char* filename)
{

	int i;
	FILE *file = fopen(filename,"r");
	char dummy[100];
	if(fgets(dummy,100,fp) != NULL )
		printf("readline:%s\n", dummy);
	for(i=0;i<N;i++) 
		ptls.type[i]=1;
		if( 0<
				fscanf(fp,"%le %le %le %le %le %le", &(ptls.position[i].x), &(ptls.position[i].y), &(ptls.position[i].z),
					&(ptls.dipole[i].x),&(ptls.dipole[i].y),&(ptls.dipole[i].z)) )
		{  }
	fclose(fp);
}

void load_from_snapshot (char* filename)
{
	FILE* input = fopen(filename, "r");
	if (input == NULL)
	{
		printf("error : for filename %s\n",filename);
	}
	Snapshot* s = read_dump( input );
	if (s == NULL)
	{
		printf("error : for read_dump %s\n","snapshot");
	}

	
	Box3* b     = & s->box;
	boxlo.x     = b->xlow; boxhi.x = b->xhigh;
	boxlo.y     = b->ylow; boxhi.y = b->yhigh;
	boxlo.z     = b->zlow; boxhi.z = b->zhigh;
	Space.x     = ( boxhi.x - boxlo.x);
	Space.y     = ( boxhi.y - boxlo.y);
	Space.z     = ( boxhi.z - boxlo.z);
	
	new_data(s->n_atoms);
	int i;
	
	for ( i = 0; i < N; i += 1 ) {
		atom* a      = & s->atoms[i];
		ptls.id[i]   = a->id;
		ptls.type[i] = a->type;
		
		ptls.diameter[i] = 1;
		ptls.mass[i]     = 1;

		ptls.position[i].x = a->x;
		ptls.position[i].y = a->y;
		ptls.position[i].z = a->z;

		ptls.charge[i]   = no_charge;
		ptls.dipole[i].x = a->mux;
		ptls.dipole[i].y = a->muy;
		ptls.dipole[i].z = a->muz;

		VZero(ptls.cell[i]);	
	}
	fclose(input);
}



#include	<stdlib.h>

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  main
 *  Description:  
 * =====================================================================================
 */
	int
main ( int argc, char *argv[] )
{

/* file:
 * 	load_from_snapshot(argv[1]);
 */
/* fcc:
 */
	real diameter ;
	int n1,n2;
	real lo,hi;
	printf("diameter = ");
	scanf("%lf",&diameter);
	printf("\nN1 = ");
	scanf("%d",&n1);
	printf("\nN2 = ");
	scanf("%d",&n2);
	printf("\nBoxlo,Boxhi = ");
	scanf("%lf,%lf", &lo, &hi);
	boxlo.x = boxlo.y = boxlo.z = lo;
	boxhi.x = boxhi.y = boxhi.z = hi;
	Space.x = Space.y = Space.z = hi-lo;

	new_data(n1);
	init_pos_Pfcc(n1,diameter);
		

	AllocMem2D(ntype_condition,  5,5,real);
	ntype_condition[0][2] = diameter/2.;
	

	init_random_type_n_if(2, n2);

	FILE * f = fopen("default_output.out","w+");
	write_data(f);
	fclose(f);
	write_data_stdout();
	return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
