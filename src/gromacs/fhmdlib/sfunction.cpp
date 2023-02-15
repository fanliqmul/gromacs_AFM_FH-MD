#include "data_structures.h"
#include "macro.h"
//#include "sfunction.h"
#include "gromacs/topology/mtop_util.h"     /* This is for gmx_mtop_atominfo_global() */
#include "gromacs/utility/smalloc.h"

//#define FHMD_DEBUG_COM
//sphere function
double fhmd_Sxyz_r(const rvec x, const dvec c, FHMD *fh)
{
    dvec xd;

    for(int d = 0; d < DIM; d++)
    {
        xd[d] = fabs(x[d] - c[d]);
        if(xd[d] > fh->box05[d]) xd[d] -= fh->box[d];
    }

    double r2 = xd[0]*xd[0] + xd[1]*xd[1] + xd[2]*xd[2];

    if(r2 <= fh->R12) return fh->Smin;      // S = Smin inside the radius R1
    if(r2 >= fh->R22) return fh->Smax;      // S = Smax outside the radius R2

    double r = sqrt(r2);

    return (r - fh->R1)*fh->RS + fh->Smin;                      // Linear interpolation from Smin to Smax
//  return tanh(((r - hhmd->R1)*hhmd->RS - 0.5)*6.0)*0.5 + 0.5; // Smoother interpolation
}


/*
//half sphere function
double fhmd_Sxyz_r(const rvec x, const dvec c, FHMD *fh)
{
    dvec xd;
    dvec c0;

    //c0 = c;

    c0[0] = c[0];
    c0[1] = c[1];
    c0[2] = 4;

    for(int d = 0; d < DIM; d++)
    {
        xd[d] = fabs(x[d] - c0[d]);
        if(xd[d] > fh->box05[d]) xd[d] -= fh->box[d];
    }

    double r2 = xd[0]*xd[0] + xd[1]*xd[1] + xd[2]*xd[2];

    if(r2 <= fh->R12)
    {
    	if(x[2] < fh->box05[2])
    		return fh->Smin;      // S = Smin inside the radius R1
    	else
    		return fh->Smax;
    }
    if(r2 >= fh->R22) return fh->Smax;      // S = Smax outside the radius R2

    if(x[2] > fh->box05[2]) return fh->Smax;

    double r = sqrt(r2);

    return (r - fh->R1)*fh->RS + fh->Smin;                      // Linear interpolation from Smin to Smax
//  return tanh(((r - hhmd->R1)*hhmd->RS - 0.5)*6.0)*0.5 + 0.5; // Smoother interpolation
}
*/




//half sphere plus 1-d thin layer half sphere based on absolute and mixed S function
/*double fhmd_Sxyz_AFM(const rvec x, const dvec c, FHMD *fh)
{
    dvec xd;
    dvec c0;
    double p1,
		   p2,
		   r,
		   r2,
		   S1,//1-d s function
		   S2,//spherical s function
		   z1,//upper limit of MD zone for 1d function
		   z2,//upper limit of buffer zone for 1d function
		   //up_sphere_buffer_bd,
    	   //low_sphere_buffer_bd_down,
		   low_sphere_buffer_bd_up;

    double S1,//1-d s function
		   S2,
		   r,
		   r2,
		   low_sphere_buffer_bd_up;
    c0[0] = c[0];
    c0[1] = c[1];
    c0[2] = fh->z_c_ball;//centre of sphere

    //p1 = 0.27;
    //p2 = 0.37;
    //z1 = fh->box[2]*p1;
    //z2 = fh->box[2]*p2;
    S1 = (fabs(x[2]) - fh->z1)*fh->zs + fh->Smin;
    //up_sphere_buffer_bd = c0[2] + fh->box[2] - fh->R2;
    low_sphere_buffer_bd_up = c0[2] + fh->R2;
    //low_sphere_buffer_bd_down = c0[2] - fh->R2;

    if (fh->z1 > c0[2])
    {
    	printf(MAKE_RED "AFM error: The centre of the spherical function is lower than upper limit of MD zone for 1d function z1=%g c0[2]=%g\n", fh->z1,c0[2]);
    	exit (1);
    }

    if (fh->z2 <= fh->z1)
    {
    	printf(MAKE_RED "AFM error: No buffer zone for 1d function\n");
    	exit (2);
    }

    if (low_sphere_buffer_bd_up >= fh->box05[2])
    {
    	printf(MAKE_RED "AFM error: The buffer boundary of spherical function  is bigger than half of box in z direction\n");
    	exit (3);
    }

    if (z2 >= up_sphere_buffer_bd)
    {
    	printf(MAKE_RED "AFM error: The buffer zone of 1d function overlap the buffer zone of upper spherical function,either decrease p2 or increase c0[2]\n");
    	exit (4);
    }

    if (low_sphere_buffer_bd_down <= 0)
    {
    	printf(MAKE_RED "AFM error: The sphere zone goes to the upper part of the box,either reduce the sphere zone or increase the increase c0[2]\n");
    	exit (5);
    }
    for(int d = 0; d < DIM; d++)
    {
        xd[d] = fabs(x[d] - c0[d]);
        if(xd[d] > fh->box05[d]) xd[d] -= fh->box[d];
    }

    r2 = xd[0]*xd[0] + xd[1]*xd[1] + xd[2]*xd[2];
    r = sqrt(r2);
    S2 = (r - fh->R1)*fh->RS + fh->Smin;

    if(x[2] < fh->z1 || r2 < fh->R12)
    {
    	if (x[2] < fh->box05[2])
    		return fh->Smin;//MD zone
    	else
    		return fh->Smax;

    }

    else if((r2 >= fh->R12 && r2 <= fh->R22) && (x[2] >= fh->z1 && x[2] <= fh->z2))
    	return S1*S2/fh->Smax;//overlapped 1d and sphere buffer zone

    else if(x[2] >= fh->z1 && x[2] <= fh->z2)
    {
    	if (r2 < fh->R12 && x[2] < fh->box05[2])
    		return fh->Smin;//MD zone
    	else if (r2 >= fh->R12 && r2 <= fh->R22 )
    		return S1*S2/fh->Smax;//overlapped 1d and sphere buffer zone
    	else
    		return S1;//1-d buffer zone
    }

    else if (r2 >= fh->R12 && r2 <= fh->R22)
    {
    	if (x[2] >fh->z2 && x[2] < fh->box05[2])
    		return S2;//sphere buffer zone
    	else if (x[2] > fh->box05[2])
    		return fh->Smax;
    }

    else if ((r2 >= fh->R12 && r2 <= fh->R22) && x[2] > fh->z2)
        	return S2;//sphere buffer zone

    else
    	return fh->Smax;
}*/

double fhmd_Sxyz_AFM(const rvec x, const dvec c, FHMD *fh)
{
    dvec xd,xz;
    dvec c0;

    double S1,//1-d s function
		   S2,
		   r,
		   r2,
		   low_sphere_buffer_bd_up;
    c0[0] = c[0];
    c0[1] = c[1];
    c0[2] = fh->z_c_ball;//centre of sphere

    int m = 2;

    low_sphere_buffer_bd_up = c0[2] + fh->R2;


    //if (fh->z1 > c0[2])
    //{
    	//printf(MAKE_RED "AFM error: The centre of the spherical function is lower than upper limit of MD zone for 1d function z1=%g c0[2]=%g\n", fh->z1,c0[2]);
    	//exit (1);
    //}

    if (fh->z2 <= fh->z1)
    {
    	printf(MAKE_RED "AFM error: No buffer zone for 1d function\n");
    	exit (2);
    }

    if (low_sphere_buffer_bd_up >= fh->box05[2])
    {
    	printf(MAKE_RED "AFM error: The buffer boundary of spherical function  is bigger than half of box in z direction\n");
    	exit (3);
    }

   /* for(int d = 0; d < DIM; d++)
    {
        xd[d] = fabs(x[d] - c0[d]);
        if(xd[d] > fh->box05[d]) xd[d] -= fh->box[d];
    }

    r2 = xd[0]*xd[0] + xd[1]*xd[1] + xd[2]*xd[2];
    r = sqrt(r2);
    S2 = (r - fh->R1)*fh->RS + fh->Smin;

    xz = fabs(x[2] - fh->box[2]);
    S1 = (xz - fh->z1)*fh->zs + fh->Smin;*/

    for (int d = 0; d<DIM; d++)
    {
    	if (d == 2)
    	{
    		if (x[d] < c[d])
    			xz[d] = x[d];
    		else
    			xz[d] = fabs(x[d] - fh->box[d]);
    	}
    	else
    		xz[d] = x[d];
    }

    S1 = (xz[m] - fh->z1)*fh->zs + fh->Smin;

    for (int d = 0; d < DIM; d++)
    {
    	xd[d] = fabs(xz[d] - c0[d]);
    }

    r2 = xd[0]*xd[0] + xd[1]*xd[1] + xd[2]*xd[2];
    r = sqrt(r2);
    S2 = (r - fh->R1)*fh->RS + fh->Smin;

    if(xz[m] < fh->z1 || r2 < fh->R12)
    {
    		return fh->Smin;//MD zone
    }


    else if (xz[m] >= fh->z1 && xz[m] <= fh->z2)
    {
    	if (r2 < fh->R12)
    		return fh->Smin;//MD zone
    	else if (r2 >= fh->R12 && r2 <= fh->R22)
    		return S1*S2/fh->Smax;//overlapped 1d and sphere buffer zone
    	else
    		return S1;//1-d buffer zone;
    }

    else if (xz[m] > fh->z2)
    {
    	if (r2 < fh->R12)
    		return fh->Smin;
    	else if (r2 >= fh->R12 && r2 <= fh->R22)
    		return S2;
    	else
    	    	return fh->Smax;
    }

    else
    {
    	printf (MAKE_RED "AFM error: one atom does not have S function\n" RESET_COLOR "\n");
    	exit (4);
    }

}



//1-d function using absolute coordinate in z direction
/*double fhmd_Sxyz_AFM(const rvec x, const dvec c, FHMD *fh)
{

    int d = 2;

    if (x[d] < fh->z1) return fh->Smin;
    if (x[d] > fh->z2) return fh->Smax;


    return (fabs(x[d]) - fh->z1)*fh->zs + fh->Smin;

}*/


//1-d function symmetry respect to the Smax at the center
/*double fhmd_Sxyz_AFM(const rvec x, const dvec c, FHMD *fh)
{

	    double xd;
	    int d = 2;

	    if (x[d] < c[d])

	    {
	    	xd = x[d];
	    }
	    else
	    {
	    	xd = fabs(x[d] - fh->box[d]);
	    }
	    //if(xd > fh->box05[d]) xd -= fh->box[d];


	    if(xd <= fh->z1) return fh->Smin;      // S = Smin inside the radius R1
	    if(xd >= fh->z2) return fh->Smax;      // S = Smax outside the radius R2

	    return (xd - fh->z1)*fh->zs + fh->Smin;
	    //return fh->A*(xd - fh->z1)*(xd - fh->z1)*(xd - fh->z1) + fh->B*(xd - fh->z1)*(xd - fh->z1) + fh->D*(xd - fh->z1); //power 3

}*/





double fhmd_Sxyz_protein (const rvec x, const dvec *c, FHMD *fh)
{
	dvec xd[fh->prot_num];
	double r2[fh->prot_num],r[fh->prot_num];
	double sb_1, sb_2;


	for (int i = 0; i < fh->prot_num; i++)
	{
		for (int d = 0; d < DIM; d++)
		{
			xd[i][d] = fabs(x[d] - fh->prot_com[i][d]);
			if(xd[i][d] > fh->box05[d]) xd[i][d] -= fh->box[d];
		}

		//printf("coord %g %g %g \n", x[0],x[1],x[2]);

		r2[i] = xd[i][0]*xd[i][0] + xd[i][1]*xd[i][1] + xd[i][2]*xd[i][2];
		r[i]  = sqrt(r2[i]);
	}

	/*part to change to determing the number of protein*/

	if(fh->prot_num == 1)
	{
		if (r2[0] <= fh->R12)
			return fh->Smin;
		else if (r2[0] >= fh->R22)
			return fh->Smax;
		else
			return (r[0] - fh->R1)*fh->RS + fh->Smin;

	}

	else if (fh->prot_num == 2)
	{
		if(r2[0] <= fh->R22 && r2[0] >= fh->R12 && r2[1] >= fh->R22)
		{
			sb_1 = (r[0] - fh->R1)*fh->RS + fh->Smin;
			return sb_1;
		}

		else if(r2[1] <= fh->R22 && r2[1] >= fh->R12 && r2[0] >= fh->R22)
		{
			sb_2 = (r[1] - fh->R1)*fh->RS + fh->Smin;
			return sb_2;
		}

		else if(r2[0] <= fh->R22 && r2[0] >= fh->R12 && r2[1] <= fh->R22 && r2[1] >= fh->R12)
		{
			return sb_1*sb_2/fh->Smax;
		}

		else if(r2[0] <= fh->R12 && r2[1] >= fh->R22)
			return fh->Smin;

		else if(r2[1] <= fh->R12 && r2[0] >= fh->R22)
			return fh->Smin;

		else if (r2[0] <= fh->R12 && r2[1] >= fh->R12 && r2[1] <= fh->R22)
			return fh->Smin;

		else if (r2[1] <= fh->R12 && r2[0] >= fh->R12 && r2[0] <= fh->R22)
			return fh->Smin;

		else if (r2[0] < fh->R12 && r2[1] < fh->R12)
			return fh->Smin;
		else if (r2[0] > fh->R22 && r2[1] > fh->R22)
			return fh->Smax;

		else
		{
			printf(MAKE_RED "FHMD: one of the atoms does not have S function \n");
			exit(5);
		}
	}

	else
	{
		printf(MAKE_RED "FHMD: the code can only deal with 1 or 2 or 4 proteins" RESET_COLOR "\n");
	}
}


double fhmd_Sxyz_d(const dvec x, const dvec c, FHMD *fh)
{
    rvec xr;
    copy_dvec_to_rvec(x, xr);
    return fhmd_Sxyz_r(xr, c, fh);
}

double fhmd_Sxyz_d_AFM(const dvec x, const dvec c, FHMD *fh)
{
	rvec xr;
	copy_dvec_to_rvec(x, xr);
	return fhmd_Sxyz_AFM(xr, c, fh);
}

double fhmd_Sxyz_d_protein(const dvec x, const dvec *c, FHMD *fh)
{
	rvec xr;
	copy_dvec_to_rvec(x, xr);
	return fhmd_Sxyz_protein(xr, c, fh);
}

/*
 ******************** Estimate S in the cells and cell faces ********************
 */
void FH_S_weighted(FHMD *fh)
{
    for(int i = 0; i < fh->Ntot; i++)
    {
        if(fh->S_function != constant_S)
            fh->arr[i].S = 1 - fh->arr[i].ro_md_s/fh->arr[i].ro_md;
        else
            fh->arr[i].S = fh->S;
    }

    ivec ind;

    for(int k = 0; k < NZ; k++)
    {
        for(int j = 0; j < NY; j++)
        {
            for(int i = 0; i < NX; i++)
            {
                ASSIGN_IND(ind, i, j, k);

                for(int d = 0; d < DIM; d++)
                    fh->arr[L].Sf[d] = (fh->arr[CL].S + fh->arr[C].S)*0.5;
            }
        }
    }
}


void FH_S_precise(FHMD *fh)
{
    dvec coordl, coordr;
    ivec ind;

    for(int k = fh->N_shift[2]; k < (fh->N_md[2] + fh->N_shift[2]); k++)
    {
        for(int j = fh->N_shift[1]; j < (fh->N_md[1] + fh->N_shift[1]); j++)
        {
            for(int i = fh->N_shift[0]; i < (fh->N_md[0] + fh->N_shift[0]); i++)
            {
                ASSIGN_IND(ind, i, j, k);

                if(fh->S_function == moving_sphere)
                    //fh->arr[C].S = fhmd_Sxyz_d(fh->grid.c[C], fh->protein_com, fh);         // MD/FH sphere follows protein
                	fh->arr[C].S = fhmd_Sxyz_d_protein(fh->grid.c[C], fh->prot_com, fh);
                else if(fh->S_function == fixed_sphere)
                    fh->arr[C].S = fhmd_Sxyz_d(fh->grid.c[C], fh->box05, fh);               // Fixed MD/FH sphere
                else if(fh->S_function == AFM)                                              // AFM MD/FH
                	fh->arr[C].S = fhmd_Sxyz_d_AFM(fh->grid.c[C], fh->box05, fh);
                else
                    fh->arr[C].S = fh->S;                                                   // Constant S

                for(int d = 0; d < DIM; d++)
                {
                    if(fh->S_function != constant_S)
                    {
                        switch(d)
                        {
                        case 0:
                            ASSIGN_DVEC(coordl, fh->grid.n[C][0], fh->grid.c[C][1], fh->grid.c[C][2]);
                            ASSIGN_DVEC(coordr, fh->grid.n[CR][0], fh->grid.c[C][1], fh->grid.c[C][2]);
                            break;
                        case 1:
                            ASSIGN_DVEC(coordl, fh->grid.c[C][0], fh->grid.n[C][1], fh->grid.c[C][2]);
                            ASSIGN_DVEC(coordr, fh->grid.c[C][0], fh->grid.n[CR][1], fh->grid.c[C][2]);
                            break;
                        case 2:
                            ASSIGN_DVEC(coordl, fh->grid.c[C][0], fh->grid.c[C][1], fh->grid.n[C][2]);
                            ASSIGN_DVEC(coordr, fh->grid.c[C][0], fh->grid.c[C][1], fh->grid.n[CR][2]);
                            break;
                        }

                        if(fh->S_function == moving_sphere)
                        {
                            //fh->arr[L].Sf[d] = fhmd_Sxyz_d(coordl, fh->protein_com, fh);    // MD/FH sphere follows protein
                            fh->arr[L].Sf[d] = fhmd_Sxyz_d_protein(coordl, fh->prot_com,fh);
                            //fh->arr[R].Sf[d] = fhmd_Sxyz_d(coordr, fh->protein_com, fh);
                            fh->arr[R].Sf[d] = fhmd_Sxyz_d_protein(coordr, fh->prot_com,fh);
                        }
                        else if(fh->S_function == fixed_sphere)
                        {
                            fh->arr[L].Sf[d] = fhmd_Sxyz_d(coordl, fh->box05, fh);          // Fixed MD/FH sphere
                            fh->arr[R].Sf[d] = fhmd_Sxyz_d(coordr, fh->box05, fh);
                        }

                        else if(fh->S_function == AFM)
                        {
                        	fh->arr[L].Sf[d] = fhmd_Sxyz_d_AFM(coordl, fh->box05, fh);     // AFM MD/FH
                        	fh->arr[R].Sf[d] = fhmd_Sxyz_d_AFM(coordr, fh->box05, fh);
                        }

                        else
                        {
                        	printf(MAKE_RED "FHMD error: the input of S function is wrong,no suitable S function is found");
                        	exit(4);
                        }
                    }
                    else
                    {
                            fh->arr[L].Sf[d] = fh->S;                                       // Constant S
                            fh->arr[R].Sf[d] = fh->S;
                    }
                }

            }
        }
    }
}


void FH_S(FHMD *fh)
{
    for(int i = 0; i < fh->Ntot; i++)
    {
        if(fh->S_function == moving_sphere)
            //fh->arr[i].S = fhmd_Sxyz_d(fh->grid.c[i], fh->protein_com, fh);         // MD/FH sphere follows protein
        fh->arr[i].S = fhmd_Sxyz_d_protein(fh->grid.c[i], fh->prot_com,fh);
        else if(fh->S_function == fixed_sphere)
            fh->arr[i].S = fhmd_Sxyz_d(fh->grid.c[i], fh->box05, fh);               // Fixed MD/FH sphere
        else if(fh->S_function == AFM)
        	fh->arr[i].S = fhmd_Sxyz_d_AFM(fh->grid.c[i], fh->box05, fh);
        else
            fh->arr[i].S = fh->S;                                                   // Constant S

        if(fh->grid.md[i] == FH_zone)
            fh->arr[i].S = 1;                                                       // Pure FH zone
    }

    ivec ind;

    for(int k = 0; k < NZ; k++)
    {
        for(int j = 0; j < NY; j++)
        {
            for(int i = 0; i < NX; i++)
            {
                ASSIGN_IND(ind, i, j, k);

                for(int d = 0; d < DIM; d++)
                    fh->arr[L].Sf[d] = (fh->arr[CL].S + fh->arr[C].S)*0.5;
            }
        }
    }
/*
    for(int k = 0; k < NZ; k++)
    {
        for(int j = 0; j < NY; j++)
        {
            for(int i = 0; i < NX; i++)
            {
                ASSIGN_IND(ind, i, j, k);

                if(fh->grid.md[C] == FH_zone)
                {
                    for(int d = 0; d < DIM; d++)
                    {
                        fh->arr[L].Sf[d] = 1;
                        fh->arr[R].Sf[d] = 1;
                    }
                }
            }
        }
    }
*/
}


/*
 ******************** Find protein molecule in the box ********************
 */

/*original fhmd_find_protein function*/
/*void fhmd_find_protein(gmx_mtop_t *mtop, int N_atoms, real mass[], t_commrec *cr, FHMD *fh)
{
    int    ind, res_nr;
    char  *atomname, *resname;
    int    protein_n    = 0;
    double protein_mass = 0;

    for(int n = 0; n < N_atoms; n++)
    {
        if(PAR(cr) && DOMAINDECOMP(cr))
        {
            ind = cr->dd->gatindex[n];
        } else {
            ind = n;
        }

        gmx_mtop_atominfo_global(mtop, ind, &atomname, &res_nr, &resname);

        if(strcmp(resname, "SOL") && strcmp(resname, "CL") && strcmp(resname, "NA"))
        {
            protein_n++;
            protein_mass += mass[n];
        }
    }

    if(PAR(cr))
    {
        gmx_sumi(1, &protein_n, cr);
        gmx_sumd(1, &protein_mass, cr);
    }

    fh->protein_N    = protein_n;
    fh->protein_mass = protein_mass;
}*/

void fhmd_find_protein(gmx_mtop_t *mtop, int N_atoms, real mass[], t_commrec *cr, FHMD *fh)
{
    int     ind, res_nr;
    char    * atomname=NULL, *resname=NULL;
    int *protein_n ;
    double *protein_mass;

    snew(protein_n,fh->prot_num);
    snew(protein_mass,fh->prot_num);

  for (int j = 0; j<fh->prot_num; j++)
    {
    	protein_n[j] = 0;
    	protein_mass[j] = 0;


    	for(int n = 0; n < N_atoms; n++)
    	{
    		if(PAR(cr) && DOMAINDECOMP(cr))
    		{
    			ind = cr->dd->gatindex[n];
    		} else {
    			ind = n;
    		}

    		gmx_mtop_atominfo_global(mtop, ind, &atomname, &res_nr, &resname);


    		if(!strcmp(resname,fh->prot_name[j]))
    		{
    			protein_n[j]++;
    			protein_mass[j]  = protein_mass[j] + mass[n];
    			//fh->prot_N[j]++;
    			//fh->prot_mass[j] = fh->prot_mass[j] + mass[n];
    		}


    	}
    }

  if(PAR(cr))
  {
      gmx_sumi(fh->prot_num, &protein_n[0], cr);  // Sum fh->prot_num elements starting from element 0
      gmx_sumd(fh->prot_num, &protein_mass[0], cr);

	  //gmx_sumi(fh->prot_num,fh->prot_N,cr);
	  //gmx_sumd(fh->prot_num,fh->prot_mass,cr);
  }

    for (int k = 0; k <fh->prot_num; k++)
	{
    	fh->prot_N[k]    = protein_n[k];
    	fh->prot_mass[k] = protein_mass[k];
	}
}


/*void fhmd_find_protein(gmx_mtop_t *mtop, int N_atoms, real mass[], t_commrec *cr, FHMD *fh)
{
    int     ind, res_nr;
    char    * atomname=NULL, *resname=NULL;

  for (int j = 0; j<fh->prot_num; j++)
    {

    	for(int n = 0; n < N_atoms; n++)
    	{
    		if(PAR(cr) && DOMAINDECOMP(cr))
    		{
    			ind = cr->dd->gatindex[n];
    		} else {
    			ind = n;
    		}

    		gmx_mtop_atominfo_global(mtop, ind, &atomname, &res_nr, &resname);


    		if(!strcmp(resname,fh->prot_name[j]))
    		{
    			fh->prot_N[j]++;
    			fh->prot_mass[j] = fh->prot_mass[j] + mass[n];

    		}


    	}
    }

  if(PAR(cr))
  {
	printf("The MPI code stopted before gmx_sumi\n");
	gmx_sumi(fh->prot_num, fh->prot_N, cr);
	gmx_sumd(fh->prot_num, fh->prot_mass, cr);

  }

}*/


/*
 ******************** Find protein molecule centre of mass ********************
 */
/*original fhmd_find_protein_com function*/
/*void fhmd_find_protein_com(gmx_mtop_t *mtop, int N_atoms, rvec x[], real mass[], t_commrec *cr, FHMD *fh)
{
    int    ind, res_nr;
    char  *atomname , *resname;
    rvec   pcom;
    dvec   r, rm;
    double xd;

    //atomname = NULL;
    //resname  = NULL;

    clear_dvec(rm);

    for(int n = 0; n < N_atoms; n++)
    {
        if(PAR(cr) && DOMAINDECOMP(cr))
        {
            ind = cr->dd->gatindex[n];
        } else {
            ind = n;
        }

        gmx_mtop_atominfo_global(mtop, ind, &atomname, &res_nr, &resname);

        if(strcmp(resname, "SOL") && strcmp(resname, "CL") && strcmp(resname, "NA"))
        {
            for(int d = 0; d < DIM; d++)
            {
                pcom[d] = fh->protein_com[d];

                xd = x[n][d] - pcom[d];

                if(fabs(xd) <= fh->box05[d])
                    r[d] = x[n][d];
                else if(xd > fh->box05[d])
                    r[d] = x[n][d] - fh->box[d];
                else    // In case if(xd < -fh->box05[d])
                    r[d] = x[n][d] + fh->box[d];

                rm[d] += mass[n]*r[d];
            }
        }
    }

    if(PAR(cr))
    {
        gmx_sumd(3, rm, cr);
    }

    for(int d = 0; d < DIM; d++)
        pcom[d] = rm[d]/fh->protein_mass;

    PBC(fh->protein_com, pcom, fh->box);

#ifdef FHMD_DEBUG_COM
    if(MASTER(cr) && !(fh->step_MD % 10000))
        printf("FHMD DEBUG: Protein COM position: %g, %g, %g\n", fh->protein_com[0], fh->protein_com[1], fh->protein_com[2]);
#endif
}*/



void fhmd_find_protein_com(gmx_mtop_t *mtop, int N_atoms, rvec x[], real mass[], t_commrec *cr, FHMD *fh)
{
    int    ind, res_nr;
    char  *atomname=NULL , *resname=NULL;
    rvec   pcom[fh->prot_num];
    dvec   r[fh->prot_num], rm[fh->prot_num];
    double xd;


    for(int i = 0; i < fh->prot_num; i++)
    {
    	clear_dvec(rm[i]);
    	clear_dvec(r[i]);
    	clear_rvec(pcom[i]);
    }

    for(int n = 0; n < N_atoms; n++)
    {
    	if(PAR(cr) && DOMAINDECOMP(cr))
    	{
    		ind = cr->dd->gatindex[n];
    	} else {
    		ind = n;
    	}

    	gmx_mtop_atominfo_global(mtop, ind, &atomname, &res_nr, &resname);

    	for (int j = 0; j < fh->prot_num; j++)
    	{
    			if(!strcmp(resname,fh->prot_name[j]))
    			{
    				for(int d = 0; d < DIM; d++)
    				{
    					pcom[j][d] = fh->prot_com[j][d];

    					xd = x[n][d] - pcom[j][d];

    					if(fabs(xd) <= fh->box05[d])
    						r[j][d] = x[n][d];
    					else if(xd > fh->box05[d])
    						r[j][d] = x[n][d] - fh->box[d];
    					else    	// In case if(xd < -fh->box05[d])
    						r[j][d] = x[n][d] + fh->box[d];

    					rm[j][d] += mass[n]*r[j][d];
    				}
    			}
    	}
    }

    for (int k = 0; k <fh->prot_num; k++)
    {
    	if(PAR(cr))
    	{
    		gmx_sumd(3, rm[k], cr);
    	}

    	for(int d = 0; d < DIM; d++)
    	{
    		if (fh->prot_mass[k] > 0)
    			pcom[k][d] = rm[k][d]/fh->prot_mass[k];
    		else
    			pcom[k][d] = 0;
    	}

    	PBC(fh->prot_com[k], pcom[k], fh->box);

	#ifdef FHMD_DEBUG_COM
    	if(MASTER(cr) && !(fh->step_MD % 10000))
    		printf("FHMD DEBUG: Protein COM position: %g, %g, %g\n", fh->prot_com[k][0], fh->prot_com[k][1], fh->prot_com[k][2]);
	#endif
    }
}


/*
void copy_com(FHMD *fh)
{
	int ind ;

	for (int i = 0; i < fh->prot_num; i ++)
	{
		if(fh->prot_N[i] == 0)
	    {
			ind = i;
			if (ind == 0)
				copy_dvec(fh->prot_com[ind+1],fh->prot_com[ind]);
			else if (ind == 1)
				copy_dvec(fh->prot_com[ind-1],fh->prot_com[ind]);
			else
			{
				printf("Only one protein, but the center is not copied\n");
				exit(6);
			}

	    }

		else
			continue;

	}


}
*/
