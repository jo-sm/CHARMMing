/*
   File:  FNtest3.C
          This is an alternating WR test using random order and pattern.
          by many cycles of open -> [read->write](n) ->  close.  

 */

#include <iostream.h>
#include <stdlib.h>
#include "FileNavigator.h"

static int ioTest(FileNavigator *fnav, const char *func, int pattern, int &ind);


int main(int argc, char *argv[]) {

  char *pname=NULL;
  char *iFile=NULL, *oFile=NULL;
  FileNavigator * fnav = NULL;

  int err,err1, err2;
  int i, j, k, m;
  int ind[100], iPat[100];
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

  fnav = new FileNavigator();
  if (((err1 = fnav->OpenFile(oFile, WRITE_MODE,1)) < 0) || ((err2 = fnav->ReadFileHeader()) < 0)) {
    cerr << "** File open error" << endl;
    cerr << "err1 = " << err1 << " err2 = " << err2 << endl;
    delete fnav;
    fnav = NULL;
    exit(1);
  }

  for (i=0; i < 7; i++) {
    iPat[i] = i % 7;
    ioTest(fnav,"W", iPat[i] , ind[i]);
  }

  fnav->CloseFile();
  if (fnav) delete fnav;

  for (k = 0; k < 5; k++) {

    fnav = new FileNavigator();
    if (((err1 = fnav->OpenFile(oFile, UPDATE_MODE,1)) < 0) || ((err2 = fnav->ReadFileHeader()) < 0)) {
      cerr << "** File open error" << endl;
      cerr << "err1 = " << err1 << " err2 = " << err2 << endl;
      delete fnav;
      fnav = NULL;
      exit(1);
    }
    fnav->PrintIndex();
    
    for (i=7; i < 100; i++) {
      iPat[i] = rand() % 7;
#if 0      
      cerr << endl << endl;
      cerr << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++";
      cerr << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
      cerr << "Cycle "<< i << " pattern " << iPat[i] << endl << endl;
#endif
      ioTest(fnav,"W", iPat[i] , ind[i]);
      fnav->PrintIndex();
      
      for (j = i; j >= 0; j--) {
	m = rand() % i;
	err = ioTest(fnav,"R", iPat[m] , ind[m]);
	if (err < 0) exit(1);
      }
    }
    
    fnav->CloseFile();
    if (fnav) delete fnav;
  }
  exit(0);
}


static int ioTest(FileNavigator *fnav, const char *func, int pattern, int &ind) {

  const int NBLOCKS=2;
  const int DATALENGTH=1000;
  const int NSTRINGS=10;

  static uWord iUWord;
  static Word  iWord;
  static Word *iArray;
  static uWord *iUArray;
  static char *iString;
  static char **iStrings;


  static int init = 0;
  int i, j, k, err;

  uWord ns,dl;
  char *tString = NULL;
  char **tStrings = NULL;
  Word tWord, index, *tArray = NULL, *tArray1 = NULL;
  uWord tUWord,  *tUArray1 = NULL;

  
  if (!init) {
    init = 1;
    iWord=655431;
    iUWord=55666;
    j= 65;
    iArray  = (Word *) calloc(NBLOCKS*DATALENGTH,sizeof(Word));
    iUArray = (uWord *) calloc(NBLOCKS*DATALENGTH,sizeof(uWord));
    iString = (char *) calloc(NBLOCKS*DATALENGTH+1,sizeof(char));

    for (i=0; i < NBLOCKS*DATALENGTH; i++) {
      iArray[i] = i; 
      iString[i] = j;  j++;
      if (j >89) j = 65; 
    }
    iString[i]  = '\0';


    iStrings = (char **) calloc(NSTRINGS, sizeof(char *));
    for (k=0; k < NSTRINGS; k++) {
      iStrings[k] = (char *) calloc(NBLOCKS*DATALENGTH+1, sizeof(char));
      j=65;
      for (i=0; i < NBLOCKS*DATALENGTH; i++) {
	iStrings[k][i] = j; j++;
	if (j >89) j = 65;
      }
      iStrings[k][i]  = '\0';
    }
  }

  if (!strcmp(func,"R")) {
    index = ind;

    switch (pattern) {
   
    case 0:
      tString = fnav->GetString(index, err);
      if (err) {
	cerr << "Error reading pattern " << pattern << " index " << index << endl;
	fnav->PrintError(err);
	return -200;
      }

      if (!tString) {
	  cerr << "Read  NULL string " << endl;
	  return -pattern;
      } else if ( tString && strcmp(tString,iString) ) {
	  cerr << "** Failed pattern " << pattern << " index = "<< index << endl;
	  cerr << "Write length " << strlen(iString) << " String = " << iString << endl;
	  cerr << "Read  length " << strlen(tString) << " String = " << tString << endl;
	  return -pattern;
      }
//      if (tString) delete tString;
      if (tString) free(tString);
      break;

    case 1:
      //  
      ns = NSTRINGS;
      tStrings = fnav->GetStrings(index, ns, err);
      if (err) {
	cerr << "Error reading pattern " << pattern << " index " << index << endl;
	fnav->PrintError(err);
	return -200;
      }
      
      for (i=0; i < ns; i++) {
	if (!tStrings[i]) {
	  cerr << "Read  NULL string " << endl;
	  return -pattern;
	} else if ( tStrings[i] && strcmp(tStrings[i],iStrings[i])) {
	  cerr << "** Failed pattern " << pattern << " index = "<< index << endl;
	  cerr << "Write length " << strlen(iStrings[i]) << " String = " << iStrings[i] << endl;
	  cerr << "Read  length " << strlen(tStrings[i]) << " String = " << tStrings[i] << endl;
	  return -pattern;
	}
      }
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

      break;

    case 2:
      dl = DATALENGTH;
      tArray = fnav->GetWords(index, dl, err);
      if (err) {
	cerr << "Error reading pattern " << pattern << " index " << index << endl;
	fnav->PrintError(err);
	return -200;
      }
      for (i=0; i < dl; i++) { 
	if (tArray[i] != iArray[i]) {
	  cerr << "** Failed pattern " << pattern << " index = "<< index << endl;
	  return -pattern;
	}
      }
//      if (tArray)  delete tArray;
      if (tArray)  free(tArray);
      break;

    case 3:

      tWord = fnav->GetWord(index,err);
      if (err) {
	cerr << "Error reading pattern " << pattern << " index " << index << endl;
	fnav->PrintError(err);
	return -200;
      }
      if (tWord != iWord) {
	cerr << "** Failed pattern " << pattern << " index = "<< index << endl;
	return -pattern;
      }
      break;

    case 4:
      dl = DATALENGTH;
      tArray1 = fnav->GetWords(index, dl, err);
      if (err) {
	cerr << "Error reading pattern " << pattern << " index " << index << endl;
	fnav->PrintError(err);
	return -200;
      }

      for (i=0; i < dl; i++) { 
	if (tArray1[i] != iArray[i]) {
	  cerr << "** Failed pattern " << pattern << " index = "<< index << endl;
	  return -pattern;
	}
      }
//      if (tArray1) delete tArray1;
      if (tArray1) free(tArray1);
      break;


    case 5:

      tUWord = fnav->GetUWord(index,err);
      if (err) {
	cerr << "Error reading pattern " << pattern << " index " << index << endl;
	fnav->PrintError(err);
	return -200;
      }
      if (tUWord != iUWord) {
	cerr << "** Failed pattern " << pattern << " index = "<< index << endl;
	return -pattern;
      }
      break;

    case 6:
      dl = DATALENGTH;
      tUArray1 = fnav->GetUWords(index, dl, err);
      if (err) {
	cerr << "Error reading pattern " << pattern << " index " << index << endl;
	fnav->PrintError(err);
	return -200;
      }

      for (i=0; i < dl; i++) { 
	if (tUArray1[i] != iUArray[i]) {
	  cerr << "** Failed pattern " << pattern << " index = "<< index << endl;
	  return -pattern;
	}
      }
//      if (tUArray1) delete tUArray1;
      if (tUArray1) free(tUArray1);
      break;


    default:
      dl = DATALENGTH;
      tArray1 = fnav->GetWords(index, dl, err);
      if (err) {
	cerr << "Error reading pattern " << pattern << " index " << index << endl;
	fnav->PrintError(err);
	return -200;
      }

      for (i=0; i < dl; i++) { 
	if (tArray1[i] != iArray[i]) {
	  cerr << "** Failed pattern " << pattern << " index = "<< index << endl;
	  return -pattern;
	}
      }
//      if (tArray1) delete tArray1;
      if (tArray1) free(tArray1);
      break;

    }

  } else { // Write functions ... 

    switch (pattern) {
    case 0:
      err = fnav->WriteString(iString, index);
      if (err) {
	cerr << "Error writing pattern " << pattern << " index " << index << endl;
	fnav->PrintError(err);
	return -101;
      }
      break;

    case 1:
      err = fnav->WriteStrings(iStrings, NSTRINGS, index);
      if (err) {
	cerr << "Error writing pattern " << pattern << " index " << index << endl;
	fnav->PrintError(err);
	return -101;
      }
      break;

    case 2:
      err = fnav->WriteWords(iArray, DATALENGTH, index);
      if (err) {
	cerr << "Error writing pattern " << pattern << " index " << index << endl;
	fnav->PrintError(err);
	return -101;
      }
      break;

    case 3:
      err = fnav->WriteWord(iWord,index);
      if (err) {
	cerr << "Error writing pattern " << pattern << " index " << index << endl;
	fnav->PrintError(err);
	return -101;
      }
      break;

    case 4:
      err=fnav->WriteWords(iArray, DATALENGTH, index);
      if (err) {
	cerr << "Error writing pattern " << pattern << " index " << index << endl;
	fnav->PrintError(err);
	return -101;
      }
      break;


    case 5:
      err = fnav->WriteUWord(iUWord,index);
      if (err) {
	cerr << "Error writing pattern " << pattern << " index " << index << endl;
	fnav->PrintError(err);
	return -101;
      }
      break;

    case 6:
      err=fnav->WriteUWords(iUArray, DATALENGTH, index);
      if (err) {
	cerr << "Error writing pattern " << pattern << " index " << index << endl;
	fnav->PrintError(err);
	return -101;
      }
      break;

    default:
      err=fnav->WriteWords(iArray, DATALENGTH, index);
      if (err) {
	cerr << "Error writing pattern " << pattern << " index " << index << endl;
	fnav->PrintError(err);
	return -101;
      }
      break;


    }
    ind = index;
  }
  return (1);
  
}
