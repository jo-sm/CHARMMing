/* 
    File: ISTest6.C
*/

#include <stdlib.h>
#include <iostream.h>
#include "ISTable.h"


// prototypes

void FillTestTable1(ISTable *s);
void FillTestTable2(ISTable *s);

int main(int argc, char ** argv) {
  ISTable *ss1, *ss2;
  ReVarPCifArray<int> colIndex;
  ReVarCifArray<CifString> colName;

  ss1 = new ISTable();
  FillTestTable1(ss1);

  ss2 = new ISTable();
  FillTestTable2(ss2);

  ss1->PrintTable();
  colName.Add("start_v");
  ss1->CreateKey(colName);
  ss2->PrintTable();
  ss2->CreateKey(colName);
  ss1->Merge(ss2,0);
  ss1->PrintTable();
  delete ss1;
  delete ss2;
  exit(0);
} 



void FillTestTable1(ISTable *s) {
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
// s->AddColumn("lp_lenght");

 for (i=0; i<10; i++) {
   sprintf(a,"%d",49-i);
   length.Copy(a);
	ColStart->Add(length);
   sprintf(a,"%d",i);
   length.Copy(a);
	ColEnd->Add(length);
//	ColLen->Add(length);
 }

 s->FillColumn(*ColStart,0);
 s->FillColumn(*ColEnd,1);
// s->FillColumn(*ColLen,2);

 delete ColStart;
 delete ColEnd;
 delete ColLen;
}


void FillTestTable2(ISTable *s) {
 int i;
 ReVarCifArray<CifString>* ColStart;
 ReVarCifArray<CifString>* ColEnd;
 ReVarCifArray<CifString>* ColLen;
 ColStart = new ReVarCifArray<CifString>;
 ColEnd = new ReVarCifArray<CifString>;
 ColLen = new ReVarCifArray<CifString>;
 CifString length;
 char a[5];
 s->AddColumn("lp_lenght");
 s->AddColumn("start_v");
 s->AddColumn("end_v");

 for (i=0; i<10; i++) {
   sprintf(a,"%d",42-i);
   length.Copy(a);
	ColStart->Add(length);
   sprintf(a,"%d",i);
   length.Copy(a);
	ColEnd->Add(length);
   sprintf(a,"%d",10-i);
   length.Copy(a);
	ColLen->Add(length);
 }

 s->FillColumn(*ColStart,1);
 s->FillColumn(*ColEnd,2);
 s->FillColumn(*ColLen,0);

 delete ColStart;
 delete ColEnd;
 delete ColLen;
}
