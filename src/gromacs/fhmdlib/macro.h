#ifndef MACRO_H_
#define MACRO_H_

#include "data_structures.h"

#define NX      fh->N[0]                /* Number of FH cells along X axis */
#define NY      fh->N[1]                /* Number of FH cells along Y axis */
#define NZ      fh->N[2]                /* Number of FH cells along Z axis */

#define C       I(ind, fh->N)           /* Point [i][j][k] */
#define L       I(ind, fh->N)           /* Point [i][j][k] */
#define R       I3d(i,j,k,+1,d,fh->N)   /* Point [i+1][j][k] */
#define CL      I3d(i,j,k,-1,d,fh->N)   /* Point [i-1][j][k] */
#define CR      I3d(i,j,k,+1,d,fh->N)   /* Point [i+1][j][k] */
#define LL      I3d(i,j,k,-1,d,fh->N)   /* Point [i-1][j][k] */
#define RR      I3d(i,j,k,+2,d,fh->N)   /* Point [i+2][j][k] */

#define L0      I3b(-1,j,k,d,fh->N)     /* Point [-1][j][k] */
#define L1      I3b( 0,j,k,d,fh->N)     /* Point [0][j][k] */

#define Cm      Im(ind, fh->N, fh->N_md, fh->N_shift)           /* Point [i][j][k] inside inner (MD) grid */
#define CLm     I3dm(i,j,k,-1,d,fh->N,fh->N_md,fh->N_shift)     /* Point [i-1][j][k] inside inner (MD) grid */
#define CRm     I3dm(i,j,k,+1,d,fh->N,fh->N_md,fh->N_shift)     /* Point [i+1][j][k] inside inner (MD) grid */

#define IR      I(indR, fh->N)
#define IL      I(indL, fh->N)

//#define nzbin   fh->nzbin

#define SUM(f)  (f[0] + f[1] + f[2])

#define ASSIGN_IND(ind, i, j, k) \
    ind[0] = i; \
    ind[1] = j; \
    ind[2] = k;

#define ASSIGN_DVEC(vec, x, y, z) \
    vec[0] = x; \
    vec[1] = y; \
    vec[2] = z;


static void PBC(dvec xn, const rvec x, const dvec box)
{
    for(int d = 0; d < DIM; d++)
    {
        xn[d] = x[d];

        if(fabs(xn[d]) > FHMD_MAX_LENGTH)
        {
            printf(MAKE_RED "\nFHMD: ERROR: Solution diverged. Atom's coordinates: (%g, %g, %g) nm\n" RESET_COLOR "\n", x[0], x[1], x[2]);
            exit(20);
        }

        while(xn[d] < 0)       xn[d] += box[d];
        while(xn[d] >= box[d]) xn[d] -= box[d];
    }
}


static int I(const ivec ind, const ivec N)
{
    ivec indn;

    copy_ivec(ind, indn);

    for(int d = 0; d < DIM; d++)
    {
        if(indn[d] < 0)     indn[d] += N[d];
        if(indn[d] >= N[d]) indn[d] -= N[d];
    }

    return indn[0] + indn[1]*N[0] + indn[2]*N[0]*N[1];
}


static int Im(const ivec ind, const ivec N, const ivec Nm, const ivec shift)
{
    ivec indn;

    copy_ivec(ind, indn);

    for(int d = 0; d < DIM; d++)
    {
        if(indn[d] <   shift[d])          indn[d] += Nm[d];
        if(indn[d] >= (shift[d] + Nm[d])) indn[d] -= Nm[d];
    }

    return indn[0] + indn[1]*N[0] + indn[2]*N[0]*N[1];
}


static int I3(const int i, const int j, const int k, const ivec N)
{
    ivec ind;

    ASSIGN_IND(ind, i, j, k);

    return I(ind, N);
}


static int I3d(const int i, const int j, const int k, const int dir, const int d, const ivec N)
{
    ivec ind;

    switch(d)
    {
    case 0:
        ASSIGN_IND(ind, i+dir, j, k);
        break;
    case 1:
        ASSIGN_IND(ind, i, j+dir, k);
        break;
    case 2:
        ASSIGN_IND(ind, i, j, k+dir);
        break;
    }

    return I(ind, N);
}


static int I3dm(const int i, const int j, const int k, const int dir, const int d, const ivec N, const ivec Nm, const ivec shift)
{
    ivec ind;

    switch(d)
    {
    case 0:
        ASSIGN_IND(ind, i+dir, j, k);
        break;
    case 1:
        ASSIGN_IND(ind, i, j+dir, k);
        break;
    case 2:
        ASSIGN_IND(ind, i, j, k+dir);
        break;
    }

    return Im(ind, N, Nm, shift);
}


static int I3b(const int index, const int j, const int k, const int d, const ivec N)
{
    ivec ind;

    switch(d)
    {
    case 0:
        ASSIGN_IND(ind, index, j, k);
        break;
    case 1:
        ASSIGN_IND(ind, k, index, j);
        break;
    case 2:
        ASSIGN_IND(ind, j, k, index);
        break;
    }

    return I(ind, N);
}

//for bins in the buffer zone
/*static int find_ind(int n, const rvec x[],FHMD *fh)
{
	int ibinz, d = 2;
	dvec xn;
	double prd,invdelta;

	PBC(xn, x[n], fh->box);


	if (xn[d] < fh->z1)
	{
		return 0;
	}
	else if (xn[d] > fh->z2)
	{
		return fh->nzbin -1;
	}

	else
	{
		prd = fh->z2 - fh->z1;
		invdelta = (fh->nzbin - 2)/prd;
		ibinz = static_cast<int> ((xn[d] - fh->z1)*invdelta) + 1;
		return ibinz;
	}

}*/

//for the bins in the buffer zone
/*static void find_bin_center(FHMD *fh)
{
	double prd,delta;
	int    d;

	d=2;

	prd = fh->z2 - fh->z1;
	delta = prd/(fh->nzbin - 2);

	fh->bin_c[0][d] = fh->z1/2;
	fh->bin_c[fh->nzbin -1][d] =fh->z2 + (fh->box[2] - fh->z2)/2;

	for (int i = 1; i< fh->nzbin - 1; i++)
	{
		fh->bin_c[i][d] = fh->z1 + delta/2 + (i-1)*delta;
	}

	for (int j=0; j<fh->nzbin; j++)
	{
		printf("bin number is %d and center is %g \n",j,fh->bin_c[j][d]);
	}
}*/


//for bins in the buffer zone and FH
/*static int find_ind(int n, const rvec x[],FHMD *fh)
{
	int ibinz, d = 2;
	dvec xn;
	double prd,invdelta;

	PBC(xn, x[n], fh->box);


	if (xn[d] < fh->z1)
	{
		return 0;
	}

	else
	{
		prd = fh->box[2] - fh->z1;
		invdelta = (fh->nzbin - 1)/prd;
		ibinz = static_cast<int> ((xn[d] - fh->z1)*invdelta) + 1;
		return ibinz;
	}

}*/

//for the bins in the buffer zone and FH
/*static void find_bin_center(FHMD *fh)
{
	double prd,delta;
	int    d;

	d=2;

	prd = fh->box[2] - fh->z1;
	delta = prd/(fh->nzbin - 1);

	fh->bin_c[0][d] = fh->z1/2;
	//fh->bin_c[fh->nzbin -1][d] =fh->z2 + (fh->box[2] - fh->z2)/2;

	for (int i = 1; i< fh->nzbin ; i++)
	{
		fh->bin_c[i][d] = fh->z1 + delta/2 + (i-1)*delta;
	}

	for (int j=0; j<fh->nzbin; j++)
	{
		printf("bin number is %d and center is %g \n",j,fh->bin_c[j][d]);
	}
}*/

/*static int find_ind(int n, const rvec x[],FHMD *fh)
{
	int ibinz, d = 2;
	dvec xn;
	double prd,invdelta;

	PBC(xn, x[n], fh->box);

		prd = fh->box[2];
		invdelta = fh->nzbin/prd;
		ibinz = static_cast<int> (xn[d] * invdelta);
		return ibinz;

}*/

//for the bins in the buffer zone and FH
static void find_bin_center(FHMD *fh)
{
	double prd,delta;
	int    d;

	d=2;

	prd = fh->box[2];
	delta = prd/fh->nzbin;


	for (int i = 0; i< fh->nzbin ; i++)
	{
		fh->bin_c[i][d] =  delta/2 + i*delta;
		//printf("bin number is %d and center is %g \n",i,fh->bin_c[i][d]);
	}

}





#endif /* MACRO_H_ */
