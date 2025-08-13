# include <stdio.h>
# include <math.h>
# include <ctype.h>
# include <stdlib.h>
# include <string.h>
# include <Python.h>
# define MAXCHAR 9999
# define SURFTEN 0.0072
# define debug 0 

/* calculate means for interresidue energies*/
typedef struct {
	int id1;
	int id2;
	double internal;
	double vdw;
	double eel;
	double pol;
	double sas;
}IE;

char line[MAXCHAR];
char dirname[MAXCHAR] = {0};
char rootname[MAXCHAR] = {0};
char ofilename[MAXCHAR] = {0};
char ifilename[MAXCHAR] = {0};
int i,j;
FILE *fpin;

int start;
int end;
int step;
int nrg1 = 0;
int nrg2 = 0;
int maxres;
int nsnap;

IE **ie;

void i_write(void) {
int i,j;
FILE *fpout;
double vdweel;
double vdweelpol;
double eelpol;
double gbsa;
double gas;

if ((fpout = fopen(ofilename, "w")) == NULL) {
        printf("\n Cannot open file 1 to write, exit", ofilename);
        exit(1);
}       

fprintf(fpout, "#%4s %5s %9s %9s %9s %9s %9s %9s %9s %9s %9s %9s %9s\n", "ID1", "ID2", "INTERNAL", "VDW", "EEL", "POL", "SAS",
                "VDWEEL", "GAS", "EELPOL", "VDWEELPOL", "GBSA", "MMGBSA"); 
for(i=0;i<nrg1;i++)
	for(j=0;j<nrg2;j++) {
		ie[i][j].sas *= SURFTEN;
		vdweel        = ie[i][j].vdw + ie[i][j].eel;
		gas           = ie[i][j].internal + vdweel;
		vdweelpol     = vdweel + ie[i][j].pol;
		eelpol        = ie[i][j].eel + ie[i][j].pol;
		gbsa          = ie[i][j].pol + ie[i][j].sas;
		fprintf(fpout, "%5d %5d %9.4lf %9.4lf %9.4lf %9.4lf %9.4lf %9.4lf %9.4lf %9.4lf %9.4lf %9.4lf %9.4lf\n", 
			ie[i][j].id1, ie[i][j].id2, ie[i][j].internal, ie[i][j].vdw, ie[i][j].eel, ie[i][j].pol,
                        ie[i][j].sas, vdweel, gas, eelpol, vdweelpol, gbsa, gbsa+gas); 
	}
fclose(fpout);

}

void check_dimension(char *filename) {
FILE *fp;
int i;
int flag;
char tmpchar1[MAXCHAR];
char tmpchar2[MAXCHAR];
char tmpchar3[MAXCHAR];
char line2[MAXCHAR];

int id1; 
int id2;
int id1bak; 
int id2bak;
int count1; 
int count2;

if ((fp = fopen(filename, "r")) == NULL) {
        printf("\n Cannot open file: %s to read, exit", filename);
        exit(1);
}       
flag = 0;
id1bak = -999;
id2bak = -999;
count1 = -1;
count2 = -1;
for (;;) {
        if (fgets(line, MAXCHAR, fp) == NULL) break;
	sscanf(line, "%s%s%s", tmpchar1, tmpchar2, tmpchar3);
	if(strcmp(tmpchar1, "PRINT") == 0 && 
           strcmp(tmpchar2, "PAIR")  == 0 && 
           strcmp(tmpchar3, "DECOMP")== 0)
		flag = 1;
	if(flag == 1 && strncmp(line, "TDC", 3) == 0) {
		for(i=0;i<strlen(line);i++) {
			if(line[i]=='-' && line[i+1] == '>') { 
				line2[i]  =' ';
				line2[i+1]=' ';
				i++;
				continue;
			}
			line2[i]=line[i];	
		}
		sscanf(&line2[3], "%d %d",&id1, &id2);
		if(id1 != id1bak) {
			count1++;
			count2 = -1;
			id2bak = -999;
			id1bak = id1;
		}
		if(id2 != id2bak) {
			count2++;
			id2bak = id2;
		}
	}
}
nrg1 = count1+1;
nrg2 = count2+1;
if(nrg1 > nrg2) 
	maxres = nrg1;
else
	maxres = nrg2;	

fclose(fp);
}


void read_minout(char *filename) {
FILE *fp;
int i;
int flag;
char tmpchar1[MAXCHAR];
char tmpchar2[MAXCHAR];
char tmpchar3[MAXCHAR];
char line2[MAXCHAR];

int id1; 
int id2;
int id1bak; 
int id2bak;
int count1; 
int count2;
double internal, vdw, eel, pol, sas;

if ((fp = fopen(filename, "r")) == NULL) {
        printf("\n Cannot open file: %s to read, exit", filename);
        exit(1);
}       
flag = 0;
id1bak = -999;
id2bak = -999;
count1 = -1;
count2 = -1;
for (;;) {
        if (fgets(line, MAXCHAR, fp) == NULL) break;
        strcpy(tmpchar1, "");
        strcpy(tmpchar2, "");
        strcpy(tmpchar3, "");
        tmpchar1[0]='\0';
        tmpchar2[0]='\0';
        tmpchar3[0]='\0';
	sscanf(line, "%s%s%s", tmpchar1, tmpchar2, tmpchar3);
	if(flag == 0  &&
           strcmp(tmpchar1, "PRINT") == 0 && 
           strcmp(tmpchar2, "PAIR")  == 0 && 
           strcmp(tmpchar3, "DECOMP")== 0)
		flag = 1;
	if(flag == 1 && strncmp(line, "TDC", 3) == 0) 
		flag = 2;
        if(flag == 2 && strcmp(tmpchar1, "TDC") != 0)
		break;
	if(flag == 2) {
		for(i=0;i<strlen(line);i++) {
			if(line[i]=='-' && line[i+1] == '>') { 
				line2[i]  =' ';
				line2[i+1]=' ';
				i++;
				continue;
			}
			line2[i]=line[i];	
		}
		sscanf(&line2[3], "%d %d%lf%lf%lf%lf%lf",&id1, &id2, &internal, &vdw, &eel, &pol, &sas);
		if(id1 != id1bak) {
			count1++;
			count2 = -1;
			id2bak = -999;
			id1bak = id1;
		}
		if(id2 != id2bak) {
			count2++;
			id2bak = id2;
		}
		ie[count1][count2].id1       = id1;
		ie[count1][count2].id2       = id2;
		ie[count1][count2].internal += internal;
		ie[count1][count2].vdw      += vdw;
		ie[count1][count2].eel      += eel;
		ie[count1][count2].pol      += pol;
		ie[count1][count2].sas      += sas;
	}
}
nrg1 = count1+1;
nrg2 = count2+1;
fclose(fp);
}


static PyObject* ie_ave(PyObject *self, PyObject *args) {

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


int i,j;
char command[MAXCHAR];
char filename[MAXCHAR];
char sid[MAXCHAR];

if (argc == 1 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "-H") == 0) {
	printf("Usage ie_ave -d  direcotry for sander output files\n");
	printf("             -i  inputs, the root sander output files (file names take a format of root_id.out)\n");
	printf("             -o  output file\n");
	printf("             -s  starting id \n");
	printf("             -e  ending id \n");
	printf("             -st step \n");

	exit(0);
}
if (argc != 13) {
	printf("Usage ie_ave -d  direcotry for sander output files\n");
	printf("             -i  input, the root sander output files (file names take a format of root_id.out)\n");
	printf("             -o  output file\n");
	printf("             -s  starting id \n");
	printf("             -e  ending id \n");
	printf("             -st step \n");
	exit(0);
}

printf("argc = %d\n", argc);
for (int k = 0; k < argc; k++) {
    printf("argv[%d]: '%s'\n", k, argv[k] ? argv[k] : "NULL");
}

for (i = 1; i < argc; i += 2) {
	if (i+1 >= argc) {
        fprintf(stderr, "Missing argument value for %s\n", argv[i]);
        // clean up and return error
    }
    if (!argv[i+1]) {
        fprintf(stderr, "argv[%d+1] is NULL\n", i);
        // clean up and return error
    }
	if (strcmp(argv[i], "-d") == 0) {
		strcpy(dirname, argv[i + 1]);
	}
	if (strcmp(argv[i], "-i") == 0) { 
		strcpy(rootname, argv[i + 1]);
	}
	if (strcmp(argv[i], "-o") == 0) {
		strcpy(ofilename, argv[i + 1]);
	}
	if (strcmp(argv[i], "-s") == 0) {
		start=atoi(argv[i + 1]);
	}
	if (strcmp(argv[i], "-e") == 0) {
		end=atoi(argv[i + 1]);
	}
	if (strcmp(argv[i], "-st") == 0) {
		step=atoi(argv[i + 1]);
	}
}

strcpy(filename, "./");
strcat(filename, dirname);
strcat(filename, "/");
strcat(filename, rootname);
strcat(filename, "_");

/* check dimension*/
strcpy(ifilename, filename);
sprintf(sid, "%d", start);
strcat(ifilename, sid);
strcat(ifilename, ".out");
check_dimension(ifilename);
maxres += 10;

/* allocate memory*/
ie = (IE**)malloc(maxres*sizeof(IE*));
for(i = 0 ; i<maxres ; i++) {
        *(ie+i) = (IE*) malloc(maxres*sizeof(IE));
        if (*(ie+i) == NULL) {
                fprintf(stderr, "memory allocation error for *ie[i]\n");
                exit(0);
        }
}
if (ie == NULL) {
        fprintf(stderr, "memory allocation error for **ie\n");
        exit(0);
}

/* begin the real work*/
for(i=0; i<maxres; i++)
	for(j=0; j<maxres; j++) {
		ie[i][j].internal = 0.0;
		ie[i][j].vdw = 0.0;
		ie[i][j].eel = 0.0;
		ie[i][j].pol = 0.0;
		ie[i][j].sas = 0.0;
	}
nsnap = 0;
for(i=start; i<=end; i+=step) {
	strcpy(ifilename, filename);
	sprintf(sid, "%d", i);
	strcat(ifilename, sid);
	strcat(ifilename, ".out");
	read_minout(ifilename);
	printf("\n... Handing snapshot %d...", i);
	nsnap++;
}
for(i=0;i<nrg1;i++)
	for(j=i;j<nrg2;j++) {
		ie[i][j].internal=(ie[i][j].internal + ie[j][i].internal)/2.0/nsnap;
		ie[i][j].vdw     =(ie[i][j].vdw      + ie[j][i].vdw)     /2.0/nsnap;
		ie[i][j].eel     =(ie[i][j].eel      + ie[j][i].eel)     /2.0/nsnap;
		ie[i][j].pol     =(ie[i][j].pol      + ie[j][i].pol)     /2.0/nsnap;
		ie[i][j].sas     =(ie[i][j].sas      + ie[j][i].sas)     /2.0/nsnap;
		ie[j][i].internal =ie[i][j].internal;
		ie[j][i].vdw      =ie[i][j].vdw;
		ie[j][i].eel      =ie[i][j].eel;
		ie[j][i].pol      =ie[i][j].pol;
		ie[j][i].sas      =ie[i][j].sas;
	}
i_write();
free(argv);
Py_RETURN_NONE;
}

static PyMethodDef module_methods[] = {
    {"ie_ave", ie_ave, METH_VARARGS, "Run ie_ave"},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef ie_ave_module = {
    PyModuleDef_HEAD_INIT,
    "ie_ave",
    "Run ie_ave from Python",
    -1,
    module_methods
};

PyMODINIT_FUNC PyInit_ie_ave(void) {
    return PyModule_Create(&ie_ave_module);
}