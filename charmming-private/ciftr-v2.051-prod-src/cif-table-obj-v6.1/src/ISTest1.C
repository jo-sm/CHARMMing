/* 
    File: ISTest1.C
*/

#include <stdlib.h>
#include <iostream.h>
#include "ISTable.h"


// prototypes

void FillTestTable(ISTable *s);

int main(int argc, char ** argv) {

  ISTable *ss;
  char fname[100];
  ReVarPCifArray<int> list;
  ReVarCifArray<CifString> list2;

  ss = new ISTable();
  FillTestTable(ss);
  strcpy(fname,"./test/outfile.db");

  FileNavigator *fnav;
  int err;
  fnav = new FileNavigator();
  err=fnav->OpenFile(fname, WRITE_MODE, 1);
  if (err) fnav->PrintError(err);
  err=fnav->ReadFileHeader();
  if (err) fnav->PrintError(err);
  err = ss->WriteObject(fnav);
  fnav->CloseFile();
  delete fnav;

  fnav = new FileNavigator();
  ISTable gnu;
  fnav->OpenFile("./test/outfile.db", READ_MODE,0);
  fnav->ReadFileHeader();
  gnu.GetObject(err, fnav);

  gnu.DeleteRow(30);
  gnu.PrintTable();

  delete ss;
  fnav->CloseFile();
  delete fnav;
  exit(0);
}


void FillTestTable(ISTable *s) {
 int i;
 ReVarCifArray<CifString>* ColStart;
 ReVarCifArray<CifString>* ColEnd;
 ReVarCifArray<CifString>* ColLen;
 ColStart = new ReVarCifArray<CifString>;
 ColEnd = new ReVarCifArray<CifString>;
 ColLen = new ReVarCifArray<CifString>;
 CifString length;
 char a[5];
 s->AddColumn("start_v");
 s->AddColumn("end_v");
 s->AddColumn("lp_lenght");

 for (i=0; i<10; i++) {
   sprintf(a,"%d",49-i);
   length.Copy(a);
	ColStart->Add(length);
 }
 for (i=0; i<12; i++) {
   sprintf(a,"%d",i);
   length.Copy(a);
	ColEnd->Add(length);
	ColLen->Add(length);
 }
 s->FillColumn(*ColStart,0);
 s->FillColumn(*ColEnd,1);
 s->FillColumn(*ColLen,2);

 delete ColStart;
 delete ColEnd;
 delete ColLen;
}
