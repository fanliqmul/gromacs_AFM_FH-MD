#include "data_structures.h"
#include "parser.h"


int parse_prm(char const *fname, FHMD *fh)
{
    const char c = ';';     // comment delimiter

    FILE *fprm;

    if((fprm = fopen(fname, "r")) == NULL) return -1;   // file not found

    char line[255];
    int  ok = 1;

    //char name_1[255];
    //char name_2[255];

    while(fscanf(fprm, "%s", line) != -1)
    {
        if(line[0] == c) skip_line(fprm);   // skip comment

        if(!strcmp(line, "Scheme") || !strcmp(line, "scheme"))
            ok = assign_int_value(&fh->scheme, line, fprm);
        else if(!strcmp(line, "S") || !strcmp(line, "s"))
            ok = assign_double_value(&fh->S, line, fprm);
        else if(!strcmp(line,"flow_type"))
        	ok = assign_int_value(&fh->flow, line, fprm);
        else if(!strcmp(line, "R1"))
            ok = assign_double_value(&fh->R1, line, fprm);
        else if(!strcmp(line, "R2"))
            ok = assign_double_value(&fh->R2, line, fprm);
        else if(!strcmp(line, "z1"))
        	ok = assign_double_value(&fh->z1, line, fprm);
        else if(!strcmp(line, "z2"))
        	ok = assign_double_value(&fh->z2, line, fprm);
        else if(!strcmp(line, "z_c_ball"))
        	ok = assign_double_value(&fh->z_c_ball, line, fprm);
        //else if (!strcmp(line, "z3"))
        	//ok = assign_double_value(&fh->z3, line, fprm);
        //else if (!strcmp(line, "z4"))
        	//ok = assign_double_value(&fh->z4, line, fprm);
        else if(!strcmp(line,"nzbin"))
        	ok = assign_int_value(&fh->nzbin, line, fprm);
        //else if(!strcmp(line,"md_c"))
        	//ok = assign_double_value(&fh->md_c, line, fprm);
        //else if(!strcmp(line,"nwater"))
        	//ok = assign_int_value(&fh->nwater, line, fprm);
        else if(!strcmp(line, "Smin"))
            ok = assign_double_value(&fh->Smin, line, fprm);
        else if(!strcmp(line, "Smax"))
            ok = assign_double_value(&fh->Smax, line, fprm);
        else if (!strcmp(line,"gamma"))
        	ok = assign_double_value(&fh->gamma, line, fprm);
        else if (!strcmp(line,"thermostat"))
        	ok = ok = assign_double_value(&fh->thermostat, line, fprm);
        else if (!strcmp(line,"eta"))
        	ok = assign_double_value(&fh->eta, line, fprm);
        else if(!strcmp(line, "alpha"))
            ok = assign_double_value(&fh->alpha, line, fprm);
        else if(!strcmp(line, "beta"))
            ok = assign_double_value(&fh->beta, line, fprm);
        else if(!strcmp(line, "gamma_x"))
            ok = assign_double_value(&fh->gamma_x, line, fprm);
        else if(!strcmp(line, "gamma_u"))
            ok = assign_double_value(&fh->gamma_u, line, fprm);
        else if(!strcmp(line, "eps_rho"))
            ok = assign_double_value(&fh->eps_rho, line, fprm);
        else if(!strcmp(line, "eps_mom"))
            ok = assign_double_value(&fh->eps_mom, line, fprm);
        else if(!strcmp(line, "S_berendsen") || !strcmp(line, "s_berendsen"))
            ok = assign_double_value(&fh->S_berendsen, line, fprm);
        else if(!strcmp(line, "protein_num"))
        	ok = assign_int_value(&fh->prot_num,line,fprm);
        else if(!strcmp(line,"protein_name"))
        	ok = assign_string_value(fh->tem,line,fprm);
        //else if(!strcmp(line,"protein_name_1"))
        	//ok = assign_string_value(fh->prot_name_1,line,fprm);
        //else if(!strcmp(line,"protein_name_2"))
            //ok = assign_string_value(fh->prot_name_2,line,fprm);
        //else if(!strcmp(line,"shear_velocity"))
        	//ok = assign_double_value(&fh->vm,line,fprm);
        //else if(!strcmp(line,"shear_boundary"))
        	//ok = assign_double_value(&fh->shear_boundary,line,fprm);
        else if(!strcmp(line, "Nx"))
            ok = assign_int_value(&fh->N[0], line, fprm);
        else if(!strcmp(line, "Ny"))
            ok = assign_int_value(&fh->N[1], line, fprm);
        else if(!strcmp(line, "Nz"))
            ok = assign_int_value(&fh->N[2], line, fprm);
        else if(!strcmp(line, "NxMD"))
            ok = assign_int_value(&fh->N_md[0], line, fprm);
        else if(!strcmp(line, "NyMD"))
            ok = assign_int_value(&fh->N_md[1], line, fprm);
        else if(!strcmp(line, "NzMD"))
            ok = assign_int_value(&fh->N_md[2], line, fprm);
        else if(!strcmp(line, "FH_EOS"))
            ok = assign_int_value(&fh->FH_EOS, line, fprm);
        else if(!strcmp(line, "FH_equil"))
            ok = assign_int_value(&fh->FH_equil, line, fprm);
        else if(!strcmp(line, "FH_step"))
            ok = assign_int_value(&fh->FH_step, line, fprm);
        else if(!strcmp(line, "FH_dens"))
            ok = assign_double_value(&fh->FH_dens, line, fprm);
        else if(!strcmp(line, "FH_temp"))
            ok = assign_double_value(&fh->FH_temp, line, fprm);
        else if(!strcmp(line, "FH_blend"))
            ok = assign_double_value(&fh->FH_blend, line, fprm);
        else if(!strcmp(line, "Noutput"))
            ok = assign_int_value(&fh->Noutput, line, fprm);

        if(!ok) return 0;   // error in prm-file
    }

    fclose(fprm);



    return 1;
}


void skip_line(FILE *fprm)
{
    int ok = 1;
    char c = '0';

    while((c != '\n') && (ok == 1)) ok = fscanf(fprm, "%c", &c);
}


int assign_int_value(int *v, char *line, FILE *fprm)
{
    int ok = fscanf(fprm, "%s", line);  // skip 'equal' delimiter

    ok = fscanf(fprm, "%s", line);

    if(sscanf(line, "%d", v) != 1) {
        return 0;
    } else {
        return 1;
    }
}


int assign_double_value(double *v, char *line, FILE *fprm)
{
    float f;

    int ok = fscanf(fprm, "%s", line);  // skip 'equal' delimiter

    ok = fscanf(fprm, "%s", line);

    if(sscanf(line, "%f", &f) != 1) {
        return 0;
    } else {
        *v = f;
        return 1;
    }
}

/*int assign_string_value(char  *v, char *line, FILE *fprm)
{

	int ok = fscanf(fprm,"%s",line); // skip 'equal' delimiter

	ok = fscanf(fprm,"%s",line);


	if(sscanf(line,"%s",v)!=1)
	{
		return 0;
	}
	else
	{
		return 1;
	}
}*/

int assign_string_value(char  *v, char *line, FILE *fprm)
{

	int ok = fscanf(fprm,"%s",line); // skip 'equal' delimiter
	char *find;

    fseek(fprm,1,SEEK_CUR);

	if(fgets(v,800,fprm)==NULL)
	{
		printf(MAKE_RED "error reading protein name" RESET_COLOR "\n");

		return 0;
	}

	else
	{
		find = strchr(v,'\n');

		if (find)
		*find = '\0';
		return 1;
	}
}

