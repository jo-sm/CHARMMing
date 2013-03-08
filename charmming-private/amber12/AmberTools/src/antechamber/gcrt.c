/* GCRT */
int rgcrt(char *filename, int *atomnum, ATOM * atom, CONTROLINFO cinfo,
		  MOLINFO minfo)
{
	FILE *fpin;
	int i, index;
	int nameindex;
	int numatom;
	int read_charge = 1;
	int overflow_flag = 0;
	char tmpchar[MAXCHAR];
	char line[MAXCHAR];


	if ((fpin = fopen(filename, "r")) == NULL) {
		fprintf(stdout, "Cannot open the gcrt file %s to read, exit\n", filename);
		exit(1);
	}
	initial(cinfo.maxatom, atom, minfo.resname);
	index = 0;
	numatom = 0;
	nameindex = -1;
	for (;;) {
		if (fgets(line, MAXCHAR, fpin) == NULL) {
			break;
		}
		if (spaceline(line) == 1 || strlen(line) <= 1) {
			index++;
			continue;
		}
		if (index == 0 || index == 1) continue;
		if (read_charge ==1 && index == 2) {
			sscanf(line, "%d%d", &(minfo.icharge), &(minfo.multiplicity));
			read_charge = 0;
			continue;
		}
		if (index == 3) {
			break;
		}
		sscanf(line, "%s", tmpchar);
		if (nameindex == -1) {
			if (tmpchar[0] >= '0' && tmpchar[0] <= '9')
				nameindex = 0;
			else
				nameindex = 1;
		}
		if (overflow_flag == 0) {
			if (nameindex == 0)
				sscanf(line, "%d%lf%lf%lf", &atom[numatom].atomicnum,
					   &atom[numatom].x, &atom[numatom].y,
					   &atom[numatom].z);
			else
				sscanf(line, "%s%lf%lf%lf", atom[numatom].name,
					   &atom[numatom].x, &atom[numatom].y,
					   &atom[numatom].z);
		}
		numatom++;
		if (numatom >= cinfo.maxatom && overflow_flag == 0) {
			printf
				("\nInfo: the atom number exceeds the MAXATOM, reallocate memory automatically");
			overflow_flag = 1;
		}
	}
	*atomnum = numatom;
	if (nameindex == 0) {
		element(*atomnum, atom);
		for (i = 0; i < *atomnum; i++)
			strcpy(atom[i].name, atom[i].element);
	}
	if (nameindex == 1)
		atomicnum(*atomnum, atom);
/* printf("\n atom number is  %5d", *atomnum); */
	fclose(fpin);
	return overflow_flag;
}

void wgcrt(char *filename, int atomnum, ATOM atom[], MOLINFO minfo)
{
	FILE *fpout;
	int i;
	int esp_flag; 
	int iodine_flag;
	char ckeyword[MAXCHAR];
	/* int index; */

	if ((fpout = fopen(filename, "w")) == NULL) {
		fprintf(stdout, "Cannot open a file %s to write in wgcrt(), exit\n", filename);
		exit(1);
	}
	fprintf(fpout, "%s\n", "--Link1--");
	if(strlen(minfo.gn) >= 5) 
		fprintf(fpout, "%s\n", minfo.gn);
	fprintf(fpout, "%s%s\n", "%chk=", minfo.chkfile);
	if(strlen(minfo.gm) >= 4) 
		fprintf(fpout, "%s\n", minfo.gm);
//	check ESP-related keyword
	esp_flag = 0;
        for(i=0;i<strlen(minfo.gkeyword);i++)
                ckeyword[i] = toupper(minfo.gkeyword[i]);
        if((strstr(ckeyword, "POP=") != 0 || strstr(ckeyword, "POP(")  != 0)  &&
           (strstr(ckeyword, "MK")   != 0 || strstr(ckeyword, "CHELP") != 0))
                esp_flag = 1;
	if(esp_flag == 1) {
//	check if the molecule contains Iodine 
		iodine_flag = 0;
		for(i=0;i<atomnum;i++) {
			if(strcmp(atom[i].element, "I") == 0 || atom[i].atomicnum == 53) {
				iodine_flag = 1;
				break;
			}
		}
		if(iodine_flag == 1) {
			if(strstr(minfo.gkeyword, "ReadRadii") == 0 &&
                           strstr(minfo.gkeyword, "READRADII") == 0 && 
                           strstr(minfo.gkeyword, "readradii") == 0) {
				strcat(minfo.gkeyword, " pop=ReadRadii");	
			}
		}
		if(minfo.gv == 1) {
			if(strstr(minfo.gkeyword, "6/50=1") == 0) {
				strcat(minfo.gkeyword, " iop(6/50=1)");
			}
		}
	}	
	fprintf(fpout, "%s\n\n", minfo.gkeyword);
	fprintf(fpout, "%s\n\n", "remark line goes here");
	fprintf(fpout, "%d%4d\n", minfo.icharge, minfo.multiplicity);
	element(atomnum, atom);

	for (i = 0; i < atomnum; i++)
//		fprintf(fpout, "%5s%12.4lf    %12.4lf    %12.4lf     \n",
//				atom[i].element, atom[i].x, atom[i].y, atom[i].z);
		fprintf(fpout, "%5s%16.10lf    %16.10lf    %16.10lf     \n",
				atom[i].element, atom[i].x, atom[i].y, atom[i].z);

	if(esp_flag == 1) {
		if(iodine_flag == 1 && minfo.gv == 1) {
			fprintf(fpout, "\nI     1.98\n");
			fprintf(fpout, "\n%s\n", minfo.gesp);
			fprintf(fpout, "\nI     1.98\n");
			fprintf(fpout, "\n%s\n", minfo.gesp);
		}
		if(iodine_flag == 1 && minfo.gv == 0) {
			fprintf(fpout, "\nI     1.98\n");
			fprintf(fpout, "\nI     1.98\n");
		}
		if(iodine_flag == 0 && minfo.gv == 1) {
			fprintf(fpout, "\n%s\n", minfo.gesp);
			fprintf(fpout, "\n%s\n", minfo.gesp);
		}
	}

	fprintf(fpout, "\n\n");
	fclose(fpout);
}
