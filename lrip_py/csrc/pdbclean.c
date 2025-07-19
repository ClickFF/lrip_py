# include <stdio.h>
# include <math.h>
# include <ctype.h>
# include <stdlib.h>
# include <string.h>
# include <Python.h>
# define MAXSS 512
# define MAXCYS 512
# define MAXNONRES 128
# define MAXRES 5120
# define MAXPARM 512
# define MAXCHAR 128
# define COLORTEXT "YES"
# define SS_BOND_CUTOFF 3.0
# define CN_BOND_CUTOFF 2.0
# define debug 0
char *affinityhome;
typedef struct {
	int id;
	char chain1;
	char chain2;
	int rid1;
	int rid2;
	int new_rid1;
	int new_rid2;
} SS;
typedef struct {
	int rid;
	char res[10];
} RES;

typedef struct {
	int rid;
	char chain;	
	double x;
	double y;
	double z;
} CYS;

typedef struct {
	char elem[10];
	char a_name[10]; 
	char r_name[10];
	char c_name;	
	int a_id;
	int r_id;
	int m_id;
	int a_no;
	int r_no;
	double x;
	double y;
	double z;
	double occupancy;
	double tempfactor;
	int ter;  /* =1 means there is a ter line printed before this atom*/
	int type; /*0 non-standard residue, 1: aa, 2:na, 3: non-standard residue specified in non_res*/
} ATOM;

typedef struct {
	char res[20];
	char atomnames[1280];
} PARM;

typedef struct {
        int resid;
	int resno;
        char res; /*according to info in SEQRES*/
        char chain;
} SEQ;

SEQ aaseq [MAXRES];
SEQ aaseq2 [MAXRES];
int nres = 0;
int nres2 = 0;

char ifilename[MAXCHAR];
char ofilename[MAXCHAR];
char fasta_filename[MAXCHAR];
char free_filename[MAXCHAR];
char pfilename[MAXCHAR];
char lfilename[MAXCHAR]="PDBCLEAN.LOG";
char leapfilename[MAXCHAR]="PDBCLEAN.LEAPIN";
char line[MAXCHAR];
char non_res[MAXCHAR];
FILE *fp, *fpin, *fpout, *fpparm, *fplog, *fpleap;
int iparmfile = 0;
int i_non_res = 0;
int iter = 1;
int irr  = 1;
int iss  = 1;
int icap = 1;
int i_fasta = 0;
int i_free = 0;
char conf_id[MAXCHAR]; /*conformation id from input*/
char cid;
int  i_conf_id = 0;
int  multiple_conf = 0;
char confids[128];

int  nss = 0;
SS ss[MAXSS];
int  nparm = 0;
PARM parm[MAXPARM];
int ncys=0;
CYS cys[MAXCYS];
ATOM at;
ATOM old_at;
ATOM at_C;
int nnonres = 0;
RES nonres[MAXNONRES];

/* for reading pdb */
int a_id;
int r_id;
int a_no = 1;
int r_no = 1;
int m_id = -999;
char a_name[MAXCHAR];
char a_name2[MAXCHAR];
char r_name[MAXCHAR];
char c_name;
char elem[MAXCHAR];
double x,y,z;
double tempfactor, occupancy;
int terindex =-1;

int a_count = 1;
int r_count = 1;
int i_found_gap=0;

/*
1 - 6 Record name "ATOM "
7 - 11 Integer serial Atom serial number.
13 - 16 Atom name Atom name.
17 Character altLoc Alternate location indicator.
18 - 20 Residue name resName Residue name.
22 Character chainID Chain identifier.
23 - 26 Integer resSeq Residue sequence number.
27 AChar iCode Code for insertion of residues.
31 - 38 Real(8.3) x Orthogonal coordinates for X in Angstroms
39 - 46 Real(8.3) y Orthogonal coordinates for Y in Angstroms
47 - 54 Real(8.3) z Orthogonal coordinates for Z in Angstroms
55 - 60 Real(6.2) occupancy Occupancy.
61 - 66 Real(6.2) tempFactor Temperature factor.
77 - 78 LString(2) element Element symbol, right-justified.
79 - 80 LString(2) charge Charge on the atom.
12345678901234567890123456789012345678901234567890123456789012345678901234567890
ATOM    145  N   VAL A  25      32.433  16.336  57.540  1.00 11.92           N
ATOM    146  CA  VAL A  25      31.132  16.439  58.160  1.00 11.85           C  
ATOM  15014  N   SER  1000      73.328  80.218  90.964  1.00  0.00           N
*/
char *extract(char *str, char *str2, int start, int end) {
int i;
int count = 0;
        strcpy(str2, "");
        str2[0] = '\0';
        for(i=start; i<=end;i++) {
                if(str[i] == ' ') continue;
                str2[count++] = str[i];
        }
        str2[count] = '\0';
        return str2;
}

void decipher (void) {
	int i;
	int count;
	count = 0;
	for(i=0;i<strlen(non_res); i++) {
		if(non_res[i]==',') non_res[i]=' ';
	}
	for(i=0;i<strlen(non_res); i++) {
		if(non_res[0] != ' ') {
			sscanf(non_res, "%s", nonres[nnonres].res); 
			nnonres++;
			continue;
		}
		if(non_res[i] == ' '&& non_res[i+1] != ' ' && non_res[i+1] != '\0' &&
		   non_res[i+1] != '\t' && non_res[i+1] != '\n') {
			sscanf(non_res, "%s", nonres[nnonres].res); 
			nnonres++;
			continue;
		}
	}
}
void scan() {
int i;
double x, y, z;
char tmpchar[MAXCHAR];
char tmpchar1[MAXCHAR];
char tmpchar2[MAXCHAR];
char tmpchar3[MAXCHAR];
char tmpchar4[MAXCHAR];
nss=0;
for (;;) {
	if (fgets(line, MAXCHAR, fpin) == NULL) 
        	break;
        if (strncmp("ANISOU", line, 6) == 0) continue;
        if (strncmp("CONECT", line, 6) == 0) continue;
	if (strncmp("MODEL", line, 5) == 0) {
		if(m_id > -990) 
			break; /* 	only the first model is read in*/
        	sscanf(&line[6], "%d", &m_id);
		a_no = 1;
		r_no = 1;
	}
/*      read in SSBOND */
        if (strncmp("SSBOND", line, 6) == 0) {
		sscanf(&line[6], "%d%s%s%d%s%s%d", &ss[nss].id, tmpchar1, tmpchar2, &ss[nss].rid1,
	                                                        tmpchar3, tmpchar4, &ss[nss].rid2);
		ss[nss].chain1=tmpchar2[0];
		ss[nss].chain2=tmpchar4[0];
		ss[nss].new_rid1 = -1;
		ss[nss].new_rid2 = -1;
		nss++;
	}
        if (strncmp("ATOM", line, 4) == 0) {
                a_name[0] = line[12];
                a_name[1] = line[13];
                a_name[2] = line[14];
                a_name[3] = line[15];
                a_name[4] = '\0';
                sscanf(a_name, "%s", a_name2);
                if(a_name2[0] >= '0' && a_name2[0] <= '9') {
                        for(i=1;i<strlen(a_name2);i++)
                                a_name[i-1] = a_name2[i];
                                a_name[i-1] = a_name2[0];
                                a_name[i] = '\0';
                }
                else
                        strcpy(a_name, a_name2);

                r_name[0] = line[17];
                r_name[1] = line[18];
                r_name[2] = line[19];
                if (r_name[1] == ' ')
                        r_name[1] = '\0';
                else if (r_name[2] == ' ')
                        r_name[2] = '\0';
                else
                        r_name[3] = '\0';
                c_name = line[21];
	        extract(line, tmpchar, 22, 25);
                r_id = atoi(tmpchar);

		if (strcmp(r_name, "CYS") ==0 || strcmp(r_name, "CYX") == 0)
			if(strcmp(a_name, "SG") ==0) {
				sscanf(&line[30], "%lf%lf%lf", &x, &y, &z);
				cys[ncys].x = x;	
				cys[ncys].y = y;	
				cys[ncys].z = z;	
				cys[ncys].rid = r_id;
				cys[ncys].chain =c_name;
				ncys ++;
			}
	}
}
rewind(fpin);
}

int judge() {
int i,j;
int suc;
char atname[MAXCHAR];
double dist;
	
suc = 1;
/* check residue */
at.type = 0;
if(strcmp(at.r_name, "ALA") == 0 || strcmp(at.r_name, "CYS") == 0 || strcmp(at.r_name, "CYX") == 0 ||
   strcmp(at.r_name, "ASP") == 0 || strcmp(at.r_name, "GLU") == 0 || strcmp(at.r_name, "PHE") == 0 ||
   strcmp(at.r_name, "GLY") == 0 || strcmp(at.r_name, "HIS") == 0 || strcmp(at.r_name, "HID") == 0 ||
   strcmp(at.r_name, "HIE") == 0 || strcmp(at.r_name, "ILE") == 0 || strcmp(at.r_name, "LYS") == 0 ||
   strcmp(at.r_name, "LEU") == 0 || strcmp(at.r_name, "MET") == 0 || strcmp(at.r_name, "ASN") == 0 ||
   strcmp(at.r_name, "PRO") == 0 || strcmp(at.r_name, "GLN") == 0 || strcmp(at.r_name, "ARG") == 0 ||
   strcmp(at.r_name, "SER") == 0 || strcmp(at.r_name, "THR") == 0 || strcmp(at.r_name, "VAL") == 0 ||
   strcmp(at.r_name, "TRP") == 0 || strcmp(at.r_name, "TYR") == 0 ) 
	at.type = 1;
else if(strcmp(at.r_name, "DA") == 0 || strcmp(at.r_name, "DA3") == 0 || strcmp(at.r_name, "DA5") == 0 ||
   strcmp(at.r_name, "DC") == 0 || strcmp(at.r_name, "DC3") == 0 || strcmp(at.r_name, "DC5") == 0 ||
   strcmp(at.r_name, "DG") == 0 || strcmp(at.r_name, "DG3") == 0 || strcmp(at.r_name, "DG5") == 0 ||
   strcmp(at.r_name, "DT") == 0 || strcmp(at.r_name, "DT3") == 0 || strcmp(at.r_name, "DT5") == 0 ||
   strcmp(at.r_name, "RA") == 0 || strcmp(at.r_name, "RA3") == 0 || strcmp(at.r_name, "RA5") == 0 ||
   strcmp(at.r_name, "RC") == 0 || strcmp(at.r_name, "RC3") == 0 || strcmp(at.r_name, "RC5") == 0 ||
   strcmp(at.r_name, "RG") == 0 || strcmp(at.r_name, "RG3") == 0 || strcmp(at.r_name, "RG5") == 0 ||
   strcmp(at.r_name, "RU") == 0 || strcmp(at.r_name, "RU3") == 0 || strcmp(at.r_name, "RU5") == 0) 
	at.type = 2;
if(at.type == 0) {
	for(i=0; i< nnonres; i++) 
		if(strcmp(at.r_name, nonres[i].res) == 0) {
			at.type = 3;
			break;	
		}
}	
if(at.type == 0) 
	return 0;

/* check atom name */
if (at.type != 3) {
	suc = 0;
	strcpy(atname, "_");
	strcat(atname, at.a_name);
	strcat(atname, "_");
	for(i=0; i< nparm; i++) 
		if(strcmp(at.r_name, parm[i].res) == 0) {
			if(strstr(parm[i].atomnames, atname) != NULL) {
				suc=1;
				break;
			}
		}	
	if(suc == 0 && icap == 1) {
		if(at.type == 1) {
			for(i=0; i< nparm; i++) 
				if(strcmp(parm[i].res, "AA") == 0) {
					if(strstr(parm[i].atomnames, atname) != NULL) 
						suc = 1;
					break;
				}
		}
		if(at.type == 2) {
			for(i=0; i< nparm; i++) 
				if(strcmp(parm[i].res, "NA") == 0) {
					if(strstr(parm[i].atomnames, atname) != NULL) 
						suc = 1;
					break;
				}
		}
	}
}
if(suc == 0) return 0;
if(old_at.r_id < -990) return 1; 
/* check if gap exists */
if(iter == 1 && at.type == 1 && strcmp(at.a_name, "N") == 0) {
	if((at.r_no - at_C.r_no) <= 1) 
		if(at_C.c_name == at.c_name) {
			dist =(at_C.x - at.x) * (at_C.x - at.x); 
			dist+=(at_C.y - at.y) * (at_C.y - at.y); 
			dist+=(at_C.z - at.z) * (at_C.z - at.z); 
			dist=sqrt(dist);
			if(dist >= CN_BOND_CUTOFF) {
				terindex = 1;
				i_found_gap=1;
				if(icap == 1) {
					fprintf(stdout, "C: %8.3lf %8.3lf %8.3lf resid=%5d N: %8.3lf %8.3lf %8.3lf resid=%5d\n", at_C.x, at_C.y, at_C.z, at_C.r_id, at.x, at.y, at.z,at.r_id);
					fprintf(stdout, "GAP found!, please re-run the command using '-cap 0' flag\n");
				}
			}
		}
}
/* check if duplicated with the immediate previous atom */
if(strcmp(old_at.a_name, at.a_name) == 0 && 
   strcmp(old_at.r_name, at.r_name) == 0 &&
                  old_at.r_id == at.r_id && 
                  old_at.c_name==at.c_name) 
	return 0;
/* check if SS bonds */
if(strcmp(at.r_name, "CYS") == 0 || strcmp(at.r_name, "CYX") == 0) { 
	for(i=0;i<nss;i++) {
		if(at.r_id == ss[i].rid1 && at.c_name == ss[i].chain1) {
			if(strcmp(at.a_name, "HG") == 0) 
				return 0;
			ss[i].new_rid1 = at.r_no;
			strcpy(at.r_name, "CYX"); 
		}
		if(at.r_id == ss[i].rid2 && at.c_name == ss[i].chain2) {
			if(strcmp(at.a_name, "HG") == 0) 
				return 0;
			ss[i].new_rid2 = at.r_no;
			strcpy(at.r_name, "CYX"); 
		}
	}
}
return suc;
}

char translate (char *res) {
char aa;
                        aa = '-';
	if(strcmp(res, "ALA") == 0) aa = 'A';
        if(strcmp(res, "CYS") == 0) aa = 'C';
        if(strcmp(res, "CYX") == 0) aa = 'C';
        if(strcmp(res, "ASP") == 0) aa = 'D';
        if(strcmp(res, "GLU") == 0) aa = 'E';
        if(strcmp(res, "PHE") == 0) aa = 'F';
        if(strcmp(res, "GLY") == 0) aa = 'G';
        if(strcmp(res, "HIS") == 0) aa = 'H';
        if(strcmp(res, "HID") == 0) aa = 'H';
        if(strcmp(res, "HIE") == 0) aa = 'H';
        if(strcmp(res, "HIP") == 0) aa = 'H';
        if(strcmp(res, "ILE") == 0) aa = 'I';
        if(strcmp(res, "LYS") == 0) aa = 'K';
        if(strcmp(res, "LEU") == 0) aa = 'L';
        if(strcmp(res, "MET") == 0) aa = 'M';
        if(strcmp(res, "ASN") == 0) aa = 'N';
        if(strcmp(res, "PRO") == 0) aa = 'P';
        if(strcmp(res, "GLN") == 0) aa = 'Q';
        if(strcmp(res, "ARG") == 0) aa = 'R';
        if(strcmp(res, "SER") == 0) aa = 'S';
        if(strcmp(res, "THR") == 0) aa = 'T';
        if(strcmp(res, "VAL") == 0) aa = 'V';
        if(strcmp(res, "TRP") == 0) aa = 'W';
        if(strcmp(res, "TYR") == 0) aa = 'Y';
return aa;
}

void read_sequence (void){
int i;
int tmpint1, tmpint2;
char chain[MAXCHAR];
char aa[13][MAXCHAR];
for (;;) {
        if (fgets(line, MAXCHAR, fpin) == NULL) break;
	if(strncmp(line, "SEQRES", 6) == 0) {
		for(i=0; i < 13; i++) {
			strcpy(aa[i], "");
			aa[i][0]='\0';	
		}
		sscanf(&line[7], "%d%s%d %s%s%s%s%s%s%s%s%s%s%s%s%s", 
			&tmpint1, chain, &tmpint2, aa[0], aa[1], aa[2], aa[3], aa[4], aa[5], aa[6],
			aa[7], aa[8], aa[9], aa[10], aa[11], aa[12]); 
		for(i=0; i < 13; i++) {
			if(strlen(aa[i]) >= 1) {	
				aaseq[nres].chain= chain[0];
				aaseq[nres].res = translate (aa[i]);	
				nres++;
			}
		}
	}
	if(strncmp(line, "ATOM", 4) == 0 || strncmp(line, "HETATM", 6) == 0) break;
}
	rewind(fpin);	
}

void read_parm() {
int i;

for (;;) {
        if (fgets(line, MAXCHAR, fpparm) == NULL) {
                break;
        }
	if(strncmp(line, "PARM", 4) == 0) {
		sscanf(&line[4], "%s%s", parm[nparm].res, parm[nparm].atomnames);
		nparm++;
	}
}
if(debug == 1) {
	for(i=0;i<nparm;i++)
		printf("PARM %5d %5s %s\n", i+1, parm[i].res, parm[i].atomnames); 
}
}

void read_pdb() {
int i;
int suc;
char tmpchar[MAXCHAR];
char tmpchar1[MAXCHAR];
char tmpchar2[MAXCHAR];
char satomno[10];
char sresno[10];

for (;;) {
	if (fgets(line, MAXCHAR, fpin) == NULL) 
        	break;
        if (strncmp("ANISOU", line, 6) == 0) continue;
        if (strncmp("CONECT", line, 6) == 0) continue;
        if (strncmp("TER", line, 3) == 0) {
		terindex = 1;
		continue;
	}
	if (strncmp("MODEL", line, 5) == 0) {
		if(m_id > -990) 
			break; /* 	only the first model is read in*/
        	sscanf(&line[6], "%d", &m_id);
		a_no = 1;
		r_no = 1;
	}
        if (strncmp("ATOM", line, 4) == 0
        	|| strncmp("HETATM", line, 6) == 0) {
/*          --- columns 13-16 have the Brookhaven-formatted name:    */
		cid       = line[16];
		if(cid != ' ') {
			suc = 0;
			for(i=0; i< multiple_conf; i++)
				if(cid == confids[i]) {
					suc = 1;
					break;
				}
			if(suc == 0) {
				confids[multiple_conf] = cid;
			 	multiple_conf ++; 
			}
		}

		if(cid != ' ') {
			if(i_conf_id == 1 && cid != conf_id[0]) continue;
			if(i_conf_id == 0 && cid != confids[0]) continue; 
/*if not specified, only the first conformation is selected*/
		}
		
                a_name[0] = line[12];
                a_name[1] = line[13];
                a_name[2] = line[14];
                a_name[3] = line[15];
                a_name[4] = '\0';
                sscanf(a_name, "%s", a_name2);
                if(a_name2[0] >= '0' && a_name2[0] <= '9') {
                	for(i=1;i<strlen(a_name2);i++)
                        	a_name[i-1] = a_name2[i];
                                a_name[i-1] = a_name2[0];
                                a_name[i] = '\0';
                }
                else
                	strcpy(a_name, a_name2);

                r_name[0] = line[17];
                r_name[1] = line[18];
                r_name[2] = line[19];
                if (r_name[1] == ' ')
                	r_name[1] = '\0';
                else if (r_name[2] == ' ')
                	r_name[2] = '\0';
                else
                	r_name[3] = '\0';

		c_name = line[21];
	        extract(line, tmpchar, 22, 25);
                r_id = atoi(tmpchar);
	        extract(line, tmpchar, 6, 10);
                a_id = atoi(tmpchar);
		
		occupancy = 1.0;
		tempfactor = 0.0;
		strcpy(elem, "");	
		sscanf(&line[30], "%lf%lf%lf%lf%lf%s", &x, &y, &z,&occupancy,&tempfactor,elem);
		
		at.x = x;	
		at.y = y;	
		at.z = z;	
		at.occupancy = occupancy;	
		at.tempfactor = tempfactor;	
		strcpy(at.a_name, a_name);
		strcpy(at.r_name, r_name);
		strcpy(at.elem, elem);
		at.c_name = c_name;
		at.a_id = a_id;	
		at.r_id = r_id;	
		at.m_id = m_id;	
		at.a_no = a_no;
		at.r_no = r_no;
		suc=judge();
		if(suc == 0) continue;
		if(old_at.r_id > -990) {
			if(old_at.r_id != at.r_id ||old_at.c_name != at.c_name) {
				r_no++;
				at.r_no ++;
				if(i_free == 1) {
					aaseq2[nres2].res = translate(r_name);
					aaseq2[nres2].resid=r_id;
					aaseq2[nres2].resno=r_no;
					aaseq2[nres2].chain=c_name;
					nres2++;
				}
			}
		}
		if(a_no == 1 && i_free == 1) {
			aaseq2[nres2].res = translate(r_name);
			aaseq2[nres2].resid=r_id;
			aaseq2[nres2].resno=r_no;
			aaseq2[nres2].chain=c_name;
			nres2++;
		}
		a_no++;
		at.ter = terindex;
		if(strcmp(at.a_name, "C") == 0 && at.type == 1) 
			at_C=at;

		if(at.ter == 1) 
			fprintf(fpout, "TER\n");
                if(irr == 0) 
			fprintf(fpout, "%s", line);
		else {
                	sprintf(satomno, "%d", at.a_no);
                        sprintf(sresno, "%d", at.r_no);
                        if(at.a_no > 9999999) strcpy(satomno, "*******");
                        if(at.r_no > 99999) strcpy(sresno, "*****");
                       	if(strlen(a_name) <= 3)
                        	fprintf(fpout, "ATOM%7s  %-3s %3s %5s", satomno, at.a_name, at.r_name, sresno);
                       	else
                            	fprintf(fpout, "ATOM%7s %4s %3s %5s", satomno, at.a_name, at.r_name, sresno);
                        fprintf(fpout, "%12.3lf%8.3lf%8.3lf%6.2lf%6.2lf%12s\n", at.x,at.y,at.z,at.occupancy,at.tempfactor,at.elem);
		}
		old_at = at;
		terindex = 0;
		continue;
        }
	fprintf(fpout, "%s", line);	
    }
}

static PyObject* pdbclean(PyObject *self, PyObject *args) {

	PyObject *argList;
	if (!PyArg_ParseTuple(args, "O!", &PyList_Type, &argList)) {
		return NULL;
	}

	int argc = (int)PyList_Size(argList);
	char **argv = (char **)malloc(argc * sizeof(char*));
	if (!argv) return PyErr_NoMemory();

	for (int i = 0; i < argc; i++) {
		PyObject *item = PyList_GetItem(argList, i);
		if (!PyUnicode_Check(item)) {
			free(argv);
			PyErr_SetString(PyExc_TypeError, "All list items must be strings");
			return NULL;
		}
		argv[i] = PyUnicode_AsUTF8(item);
	}


	int i,j,k;
	int flag;
	char tmpchar[MAXCHAR];
	char old_chain;
	double dist;
	if (strcmp(COLORTEXT, "YES") == 0 || strcmp(COLORTEXT, "yes") == 0) {
		if (argc == 2
			&& (strcmp(argv[1], "-h") == 0
				|| strcmp(argv[1], "-H") == 0)) {
			printf
				("[31mUsage: pdbclean -i   [0m input file in pdb format\n"
				 "[31m                -o   [0m output file name in pdb format\n"
				 "[31m                -ofa [0m output file name in fasta format, optional\n"
				 "[31m                -ofr [0m output file name in free format, optional\n"
				 "[31m                -p   [0m parameter file listed the atom names of standard residues, optional\n"
				 "[31m                     [0m the default is $AFFINITY/dat/BIO_ATOM_NAME.PARM\n"
				 "[31m                -nr  [0m non-standard residue list to be kept in \"\",default is none, optional\n"
				 "[31m                -r   [0m regenerate residue ids in output file, optional\n"
				 "[34m                      1: [0m yes, the default\n"
				 "[34m                      0: [0m no\n"
				 "[31m                -at  [0m add \"TER\" line when there is a gap, even the input pdb does not have one\n"
				 "[34m                      1: [0m yes, the default\n"
				 "[34m                      0: [0m no\n"
				 "[31m                -ss  [0m predict possible disulfide bonds when there is no SSBOND fields\n"
				 "[34m                      1: [0m yes, the default\n"
				 "[34m                      0: [0m no\n"
				 "[31m                -cap [0m keep cap atoms\n"
				 "[34m                      1: [0m yes, the default\n"
				 "[34m                      0: [0m no\n"
				 "[31m                -l   [0m log file, default is PDBCLEAN.LOG, optional\n"
				 "[31m                -leap[0m leap info file, default is PDBCLEAN.LEAPIN, optional\n"
				 "[31m                -cid [0m conformation id (indicated in Column 17), optional\n");
			exit(0);
		}
		if (argc != 5  && argc !=7 && argc !=9 && argc !=11 && argc != 13 && argc != 15 && argc !=17 && argc !=19 && argc != 21 && argc !=23 && argc != 25 && argc != 27) {
			printf
				("[31mUsage: pdbclean -i   [0m input file in pdb format\n"
				 "[31m                -o   [0m output file name in pdb format\n"
				 "[31m                -ofa [0m output file name in fasta format, optional\n"
				 "[31m                -ofr [0m output file name in free format, optional\n"
				 "[31m                -p   [0m parameter file listed the atom names of standard residues, optional\n"
				 "[31m                     [0m the default is $AFFINITY/dat/BIO_ATOM_NAME.PARM\n"
				 "[31m                -nr  [0m non-standard residue list to be kept in \"\",default is none, optional\n"
				 "[31m                -r   [0m regenerate residue ids in output file, optional\n"
				 "[34m                      1: [0m yes, the default\n"
				 "[34m                      0: [0m no\n"
				 "[31m                -at  [0m add \"TER\" line when there is a gap, even the input pdb does not have one\n"
				 "[34m                      1: [0m yes, the default\n"
				 "[34m                      0: [0m no\n"
				 "[31m                -ss  [0m predict possible disulfide bonds when there is no SSBOND fields\n"
				 "[34m                      1: [0m yes, the default\n"
				 "[34m                      0: [0m no\n"
				 "[31m                -cap [0m keep cap atoms\n"
				 "[34m                      1: [0m yes, the default\n"
				 "[34m                      0: [0m no\n"
				 "[31m                -l   [0m log file, default is PDBCLEAN.LOG, optional\n"
				 "[31m                -leap[0m leap info file, default is PDBCLEAN.LEAPIN, optional\n"
				 "[31m                -cid [0m conformation id (indicated in Column 17), optional\n");
			exit(0);
		}
	} else {
		if (argc == 2
			&& (strcmp(argv[1], "-h") == 0
				|| strcmp(argv[1], "-H") == 0)) {
			printf
				("Usage: pdbclean -i    input file in pdb format\n"
				 "                -o    output file name in pdb format\n"
				 "                -ofa  output file name in fasta format, optional\n"
				 "                -ofr  output file name in free format, optional\n"
				 "                -p    parameter file listed the atom names of standard residues, optional\n"
				 "                      the default is $AFFINITY/dat/BIO_ATOM_NAME.PARM\n"
				 "                -nr   non-standard residue list to be kept in \"\",default is none, optional\n"
				 "                -r    regenerate residue ids in output file, optional\n"
				 "                      1:  yes, the default\n"
				 "                      0:  no\n"
				 "                -at   add \"TER\" line when there is a gap, even the input pdb does not have one\n"
				 "                      1:  yes, the default\n"
				 "                      0:  no\n"
				 "                -ss   predict possible disulfide bonds when there is no SSBOND fields\n"
				 "                      1:  yes, the default\n"
				 "                      0:  no\n"
				 "                -cap  keep cap atoms\n"
				 "                 1:   yes, the default\n"
				 "                      0:  no\n"
				 "                -l    log file, default is PDBCLEAN.LOG, optional\n"
				 "                -leap leap info file, default is PDBCLEAN.LEAPIN, optional\n"
				 "                -cid  conformation id (indicated in Column 17), optional\n");
			exit(0);
		}
		if (argc != 5  && argc !=7 && argc !=9 && argc !=11 && argc != 13 && argc != 15 && argc !=17 && argc !=19 && argc != 21 && argc !=23 && argc != 25 && argc != 27) {
			printf
				("Usage: pdbclean -i    input file in pdb format\n"
				 "                -o    output file name in pdb format\n"
				 "                -ofa  output file name in fasta format, optional\n"
				 "                -ofr  output file name in free format, optional\n"
				 "                -p    parameter file listed the atom names of standard residues, optional\n"
				 "                      the default is $AFFINITY/dat/BIO_ATOM_NAME.PARM\n"
				 "                -nr   non-standard residue list to be kept in \"\",default is none, optional\n"
				 "                -r    regenerate residue ids in output file, optional\n"
				 "                      1:  yes, the default\n"
				 "                      0:  no\n"
				 "                -at   add \"TER\" line when there is a gap, even the input pdb does not have one\n"
				 "                      1:  yes, the default\n"
				 "                      0:  no\n"
				 "                -ss   predict possible disulfide bonds when there is no SSBOND fields\n"
				 "                      1:  yes, the default\n"
				 "                      0:  no\n"
				 "                -cap  keep cap atoms\n"
				 "                 1:   yes, the default\n"
				 "                      0:  no\n"
				 "                -l    log file, default is PDBCLEAN.LOG, optional\n"
				 "                -leap leap info file, default is PDBCLEAN.LEAPIN, optional\n"
				 "                -cid  conformation id (indicated in Column 17), optional\n");
			exit(0);
		}

	}

	for (i = 1; i < argc; i += 2) {
		if (strcmp(argv[i], "-i") == 0)  
			strcpy(ifilename, argv[i + 1]);
		if (strcmp(argv[i], "-o") == 0)  
			strcpy(ofilename, argv[i + 1]);
		if (strcmp(argv[i], "-p") == 0) {
			strcpy(pfilename, argv[i + 1]);
			iparmfile = 1;
		}
		if (strcmp(argv[i], "-r") == 0)  
			irr=atoi(argv[i+1]);
		if (strcmp(argv[i], "-l") == 0)  
			strcpy(lfilename, argv[i + 1]);
		if (strcmp(argv[i], "-leap") == 0)  
			strcpy(leapfilename, argv[i + 1]);
		if (strcmp(argv[i], "-nr") == 0) {
			strcpy(non_res, argv[i + 1]);
			i_non_res = 1;
		}
		if (strcmp(argv[i], "-at") == 0)  
			iter=atoi(argv[i+1]);
		if (strcmp(argv[i], "-ss") == 0)  
			iss=atoi(argv[i+1]);
		if (strcmp(argv[i], "-cap") == 0)  
			icap=atoi(argv[i+1]);
		if (strcmp(argv[i], "-cid") == 0) { 
			strcpy(conf_id, argv[i + 1]);
			i_conf_id = 1;
		}
		if (strcmp(argv[i], "-ofa") == 0) { 
			strcpy(fasta_filename, argv[i + 1]);
			i_fasta = 1;
		}
		if (strcmp(argv[i], "-ofr") == 0) { 
			strcpy(free_filename, argv[i + 1]);
			i_free = 1;
		}
	}
	if ((fpin = fopen(ifilename, "r")) == NULL) {
        	fprintf(stderr, "Cannot open pdb file %s to read, exit\n", ifilename);
        	exit(1);
	}
	if ((fpout = fopen(ofilename, "w")) == NULL) {
        	fprintf(stderr, "Cannot open file %s to write, exit\n", ofilename);
        	exit(1);
	}
	if ((fplog = fopen(lfilename, "w")) == NULL) {
        	fprintf(stderr, "Cannot open log file %s to write, exit\n", lfilename);
        	exit(1);
	}
	if ((fpleap = fopen(leapfilename, "w")) == NULL) {
        	fprintf(stderr, "Cannot open log file %s to write, exit\n", leapfilename);
        	exit(1);
	}
	if(iparmfile == 0) {
    		affinityhome = (char *) getenv("AFFINITY");
    		if( affinityhome == NULL ){
       			fprintf( stdout, "AFFINITY is not set!\n" );
       			exit(1);
    		}
		strcpy(pfilename, affinityhome);
		strcat(pfilename, "/dat/BIO_ATOM_NAME.PARM");
	}
	if ((fpparm = fopen(pfilename, "r")) == NULL) {
        	fprintf(stderr, "Cannot open file %s, exit\n", pfilename);
        	exit(1);
	}
	if(iter !=0 && iter != 1) iter = 1;
	if(irr  !=0 && irr  != 1) irr  = 1;
	if(iss  !=0 && iss  != 1) iss  = 1;
	if(icap !=0 && icap != 1) icap = 1;
	if(i_fasta == 1 || i_free == 1) {
		read_sequence();
		if(debug == 1) { 
			for(i=0; i< nres; i++)
			printf("%5d %5c %5c\n", i+1, aaseq[i].res, aaseq[i].chain);
		}
	}
	
	if(i_non_res == 1) 
		decipher();
	if(iss == 1) {
		scan();
		a_no = 1;
		r_no = 1;
		m_id = -999;
	}
	if(nss >= 1) 
		iss = 0;
	else {
		nss = 0;
		for(i=0; i<ncys-1; i++) 
			for(j=i+1; j<ncys; j++) {
				dist = (cys[i].x - cys[j].x) * (cys[i].x - cys[j].x);	
				dist+= (cys[i].y - cys[j].y) * (cys[i].y - cys[j].y);	
				dist+= (cys[i].z - cys[j].z) * (cys[i].z - cys[j].z);	
				dist = sqrt(dist);
				if(dist <= SS_BOND_CUTOFF) {
					ss[nss].id = nss;
					ss[nss].chain1=cys[i].chain;	
					ss[nss].chain2=cys[j].chain;	
					ss[nss].rid1 = cys[i].rid;
					ss[nss].rid2 = cys[j].rid;
					ss[nss].new_rid1 = -1;
					ss[nss].new_rid2 = -1;
					nss++;
				}
			}		
	}
	read_parm();
	old_at.r_id = -999;
	at_C.r_no = -999;
	
	read_pdb();
	if(i_free == 1 && debug == 1) { 
		for(i=0; i< nres2; i++)
		printf("%5d %5c %5c %5d %5d\n", i+1, aaseq2[i].res, aaseq2[i].chain, aaseq2[i].resno, aaseq2[i].resid);
	}
	if(multiple_conf > 0) {
		fprintf(stdout, "REMARK MULTIPLE CONFORMATION FOUND: "); 
		fprintf(fpout, "REMARK MULTIPLE CONFORMATION FOUND: "); 
		fprintf(fplog, "REMARK MULTIPLE CONFORMATION FOUND: "); 
		for(i = 0 ; i< multiple_conf; i++) {
			fprintf(stdout, " %c ", confids[i]);
			fprintf(fpout, " %c ", confids[i]);
			fprintf(fplog, " %c ", confids[i]);
		}
		fprintf(stdout, "\n");
		fprintf(fpout, "\n");
		fprintf(fplog, "\n");
	}
	if(i_found_gap == 1) {
		fprintf(fpleap, "clearpdbresmap\n"); 
		fprintf(fpout, "REMARK ADDITIONAL TER LINE ADDED AS GAP IS FOUND\n"); 
		fprintf(fplog, "REMARK Additional TER line added as gap is found, please use \"clearpdbresmap\" in leap\n"); 
	}
    	if(nss > 0) {
		fprintf(fplog, "There are %d S-S bonds and the leap commands of adding those S-S bonds are listed below:\n", nss);
		for(i=0;i<nss;i++) {
			fprintf(fpout, "REMARK LEAP_SSBOND %5d %5d %5d\tbond mol.%d.SG mol.%d.SG S\n", i+1, ss[i].new_rid1, ss[i].new_rid2, ss[i].new_rid1, ss[i].new_rid2);
			fprintf(fplog, "REMARK LEAP_SSBOND bond mol.%d.SG mol.%d.SG S\n", ss[i].new_rid1, ss[i].new_rid2);
			fprintf(fpleap, "bond mol.%d.SG mol.%d.SG S\n", ss[i].new_rid1, ss[i].new_rid2); 
		}
    	}
	if(i_fasta == 1) {
		if ((fp = fopen(fasta_filename, "w")) == NULL) {
        		fprintf(stderr, "Cannot open fasta file %s to write, exit\n", fasta_filename);
        		exit(1);
		}
		for(i=0; i<strlen(ifilename); i++) {
			tmpchar[i]=ifilename[i];
			if(ifilename[i] == '.') {
				tmpchar[i]='\0';
				break;
			}
		}
		old_chain=aaseq[0].chain;
		fprintf(fp, ">%s :%c\n", tmpchar, aaseq[i].chain);
		for(i=0;i<nres;i++) {
			if(old_chain != aaseq[i].chain) { 
				fprintf(fp, "\n\n");
				fprintf(fp, ">%s :%c\n", tmpchar, aaseq[i].chain);
			}
			fprintf(fp, "%c", aaseq[i].res);
			old_chain=aaseq[i].chain;			
		}
		fclose(fp);	
	}

        if(i_free == 1) {
                if ((fp = fopen(free_filename, "w")) == NULL) {
                        fprintf(stderr, "Cannot open free file %s to write, exit\n", free_filename);
                        exit(1);
                }
                for(i=0; i<strlen(ifilename); i++) {
                        tmpchar[i]=ifilename[i];
                        if(ifilename[i] == '.') {
                                tmpchar[i]='\0';
                                break;
                        }
                }
		if(nres == nres2) {  /*no gap at all*/
			fprintf(fp, "MOD\t");
			for(i=0; i< nres; i++)
				fprintf(fp, "%c", aaseq[i].res);
			fprintf(fp, "\n");
			fprintf(fp, "PDB\t");
			for(i=0; i< nres; i++)
				fprintf(fp, "%c", aaseq2[i].res);
			fprintf(fp, "\n");
		}
		if(nres != nres2) {  /*gaps exist*/
			fprintf(fp, "MOD\t");
			for(i=0; i< nres; i++)
				fprintf(fp, "%c", aaseq[i].res);
			fprintf(fp, "\n");
			fprintf(fp, "PDB\t");
			for(i=0; i< nres2; i++){
				fprintf(fp, "%c", aaseq2[i].res);
				j = i + 1;
				if(j<nres2 && aaseq2[i].chain == aaseq2[j].chain) {
					for(k=aaseq2[i].resid; k<aaseq2[j].resid-1; k++)						
						fprintf(fp, "%c", '-');
				}
			}	
			fprintf(fp, "\n");
		}
                fclose(fp);
        }

fclose(fpin);
fclose(fpparm);
fclose(fpout);
fclose(fplog);
free(argv);
Py_RETURN_NONE;
}

static PyMethodDef module_methods[] = {
    {"pdbclean", pdbclean, METH_VARARGS, "Run pdbclean"},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef pdbclean_module = {
    PyModuleDef_HEAD_INIT,
    "pdbclean",
    "Run pdbclean from Python",
    -1,
    module_methods
};

PyMODINIT_FUNC PyInit_pdbclean(void) {
    return PyModule_Create(&pdbclean_module);
}