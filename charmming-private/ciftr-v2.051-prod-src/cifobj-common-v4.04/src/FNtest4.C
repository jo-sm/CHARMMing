/* 
    File: FNtest4.C
*/

#include <stdlib.h>
#include <iostream.h>
#include "FileNavigator.h"

int main(int argc, char ** argv) {

  FileNavigator * fnav=NULL;
  char fname[100];
  char fname2[100];

  int err;
  Word ret=0, tWord=12345;

  //
  strcpy(fname,"./test/outfile4.db");
  fnav = new FileNavigator();
  err=fnav->OpenFile(fname, WRITE_MODE, 1);
  if (err) fnav->PrintError(err);
  err=fnav->ReadFileHeader();
  if (err) fnav->PrintError(err);
  fnav->PrintIndex();

  err = fnav->WriteString(fname, ret);
  if (err) fnav->PrintError(err);
  err = fnav->WriteString(fname, ret);
  if (err) fnav->PrintError(err);
  err = fnav->WriteString(fname, ret);
  if (err) fnav->PrintError(err);
  err = fnav->WriteString(fname, ret);
  if (err) fnav->PrintError(err);

  err = fnav->WriteWord(tWord, ret);
  if (err) fnav->PrintError(err);
  err = fnav->WriteWord(tWord, ret);
  if (err) fnav->PrintError(err);
  err = fnav->WriteWord(tWord, ret);
  if (err) fnav->PrintError(err);
  err = fnav->WriteWord(tWord, ret);
  if (err) fnav->PrintError(err);

  err=fnav->Delete(ret);
  if (err) fnav->PrintError(err);
  fnav->CloseFile();

  int err1,err2;
  fnav = new FileNavigator();
  strcpy(fname2,"./test/outfile4.db");
  if (((err1 = fnav->OpenFile(fname2, UPDATE_MODE,1)) < 0) || ((err2 = fnav->ReadFileHeader()) < 0)) {
    cerr << "************** FILE OPEN ERROR *********************" << endl;
    cerr << "err1 = " << err1 << " err2 = " << err2 << endl;
    delete fnav;
    fnav = NULL;
    exit(1);
  }
  fnav->GetWord(ret, err);
  if (err) fnav->PrintError(err);
  cerr << "err = " << err << endl;
  fnav->CloseFile();
  exit(0);
} 

