/*
   File:  FNtest1.C
          This is a WR test in coarse blocks.  Large write followed
          by many cycles of open ->read->write->close.  

 */

#include <iostream.h>
#include <stdlib.h>
#include "FileNavigator.h"

#define NBLOCKS 2
#define DATALENGTH 1000
#define NSTRINGS   10
static Word  iArray[NBLOCKS*DATALENGTH];
static Word  oArray[NBLOCKS*DATALENGTH];
static char iString[NBLOCKS*DATALENGTH+1];
static char oString[NBLOCKS*DATALENGTH+1];
static char **iStrings;
static char **oStrings;

int main(int argc, char *argv[]) {

  char *pname=NULL;
  char *iFile=NULL, *oFile=NULL;
  FileNavigator * fnav = NULL;

  int err, err1, err2;
  int i, j, k, m;
  Word  iWord, indx[NBLOCKS*DATALENGTH];
  uWord ns,dl;

  pname = argv[0];

  if (argc < 2) {
    cerr << pname << ":  usage = -o <outputfile> " << endl;
    exit(1);
  }
  
  for (i=1; i < argc; i++) {
    if ( argv[i][0] == '-' ) {
      if (strcmp(argv[i],"-f") == 0 ) {
        i++;
        iFile = (char *) calloc(strlen(argv[i]) + 1, sizeof(char));
        strcpy(iFile,argv[i]);          
      } else if (strcmp(argv[i],"-o") == 0 ) {
        i++;
        oFile = (char *) calloc(strlen(argv[i]) + 1, sizeof(char));
        strcpy(oFile,argv[i]);          
      } else {
        cerr << pname << ":  usage = -o <outputfile> " << endl;
        exit(1);
      }
    } else {
      cerr << pname << ":  usage = -o <outputfile> " << endl;
      exit(1);
    }
  }

  j= 65;
  for (i=0; i < NBLOCKS*DATALENGTH; i++) {
    iArray[i] = i;  oArray[i] = i;
    iString[i] = j; oString[i] = j; j++;
    if (j >89) j = 65; 
  }
  iString[i]  = '\0';
  oString[i]  = '\0';

  oStrings = (char **) calloc(NSTRINGS, sizeof(char *));
  iStrings = (char **) calloc(NSTRINGS, sizeof(char *));


  for (k=0; k < NSTRINGS; k++) {
    oStrings[k] = (char *) calloc(NBLOCKS*DATALENGTH+1, sizeof(char));
    iStrings[k] = (char *) calloc(NBLOCKS*DATALENGTH+1, sizeof(char));
    j=65;
    for (i=0; i < NBLOCKS*DATALENGTH; i++) {
      iStrings[k][i] = j; oStrings[k][i] = j; j++;
      if (j >89) j = 65;
    }
    iStrings[k][i]  = '\0';
    oStrings[k][i]  = '\0';
  }
  

  fnav = new FileNavigator();
  if (((err1 = fnav->OpenFile(oFile, WRITE_MODE)) < 0) || ((err2 = fnav->ReadFileHeader()) < 0)) {
    cerr << "************** FILE OPEN ERROR *********************" << endl;
    cerr << "err1 = " << err1 << " err2 = " << err2 << endl;
    delete fnav;
      fnav = NULL;
      exit(1);
  }
    //
  m = 0;
  err = fnav->WriteString(oString, indx[m]); m++;
  if (err) fnav->PrintError(err);
  err = fnav->WriteStrings(oStrings, NSTRINGS, indx[m]); m++;
  if (err)   fnav->PrintError(err);
  
  err = fnav->WriteWords(oArray, DATALENGTH, indx[m]); m++;
  if (err)   fnav->PrintError(err);
  
  for (i = 0; i < DATALENGTH; i++) {
    err = fnav->WriteWord(i,indx[m]); m++;
    if (err)    fnav->PrintError(err);
    }
  fnav->WriteWords(oArray, DATALENGTH, indx[m]); m++;
  if (err)  fnav->PrintError(err);
  
  fnav->CloseFile();
  if (fnav) delete fnav;

  //
  //
  //  

  for (int kk = 0; kk < 5; kk++) {
#if 0
    cerr << " Begin cycle " << kk << endl;
#endif
    fnav = new FileNavigator();
    if (((err1 = fnav->OpenFile(oFile, UPDATE_MODE,1)) < 0) || ((err2 = fnav->ReadFileHeader()) < 0)) {
      cerr << "************** FILE OPEN ERROR *********************" << endl;
      cerr << "err1 = " << err1 << " err2 = " << err2 << endl;
      delete fnav;
      fnav = NULL;
      exit(1);
    }

    m = 0;
    char *tString = fnav->GetString(indx[m], err); m++;
    if (err) fnav->PrintError(err);
    if ( strcmp(tString,iString) ) {
      cerr << "DIED AT 1" << endl;
      cerr << "tString " << tString << endl;
      cerr << "iString " << iString << endl;
      exit(1);
    }
  
    ns = NSTRINGS;
    char **tStrings = fnav->GetStrings(indx[m], ns, err); m++;
    if (err)   fnav->PrintError(err);
    for (i=0; i < ns; i++) {
      if (strcmp(tStrings[i],iStrings[i])) {
	cerr << "DIED AT 2" << endl;
	exit(2); 
      }
    }

    dl = DATALENGTH;
    Word * tArray = fnav->GetWords(indx[m], dl, err); m++;
    if (err)   fnav->PrintError(err);
  
    for (i=0; i < dl; i++) { 
      if (tArray[i] != iArray[i]) {
	cerr << "tArray["<<i<<"] = "<<tArray[i]<<" iArray["<<i<<"] = "<<iArray[i]<<endl;
	cerr << "DIED AT 3" << endl;
	exit(3);
      }
    }

    for (i = 0; i < DATALENGTH; i++) {
      iWord = fnav->GetWord(indx[m],err); m++;
      if (err)    fnav->PrintError(err);
      if (iWord != i) {
	cerr << "DIED AT 4" << endl;
	exit(4);
      }
    }

    dl = DATALENGTH;
    Word * tArray1 = fnav->GetWords(indx[m], dl, err); m++;
    if (err)   fnav->PrintError(err);
    for (i=0; i < dl; i++) { 
      if (tArray1[i] != iArray[i]) {
	cerr << "DIED AT 5" << endl;
	exit(5);
      }
    }


//    if (tArray)  delete tArray;
//    if (tArray1) delete tArray1;
//    if (tString) delete tString;

    if (tArray)   free(tArray);
    if (tArray1)  free(tArray1);
    if (tString)  free(tString);
/*
    for (i=0; i < NSTRINGS; i++) {
      if (tStrings[i]) delete tStrings[i];
    }
    if (tStrings) delete tStrings;
*/
    for (i=0; i < NSTRINGS; i++) {
      if (tStrings[i]) free(tStrings[i]);
    }
    if (tStrings) free(tStrings);

    //
    m = 0;
    err = fnav->WriteString(oString, indx[m]); m++;
    if (err) fnav->PrintError(err);
    err = fnav->WriteStrings(oStrings, NSTRINGS, indx[m]); m++;
    if (err)   fnav->PrintError(err);

    err = fnav->WriteWords(oArray, DATALENGTH, indx[m]); m++;
    if (err)   fnav->PrintError(err);

    for (i = 0; i < DATALENGTH; i++) {
      err = fnav->WriteWord(i,indx[m]); m++;
      if (err)    fnav->PrintError(err);
    }
    fnav->WriteWords(oArray, DATALENGTH, indx[m]); m++;
    if (err)  fnav->PrintError(err);

    fnav->CloseFile();
    if (fnav) delete fnav;
#if 0
    cerr << " END CYCLE " << kk << endl;
#endif
  }
  exit(0);
}


