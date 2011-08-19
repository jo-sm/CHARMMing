/* 
    File: ISTest8.C
	 Example fore DeleteRow, RemoveRow, CompressTable
	 DeleteRow marks row as deleted, RemoveRow removes row phisicaly
	 CompressTable removes all deleted rows and recreates indices
*/

#include <stdlib.h>
#include <iostream.h>
#include "ISTable.h"


// prototypes

void FillTestTable(ISTable *s);
void PrintArray(ReVarPCifArray<int> * target);

int main(int argc, char ** argv) {

  ISTable *ss1;
  ISTable *ss2;
  ISTable *ss3;
  ss1 = new ISTable();
  ss2 = new ISTable();
  ss3 = new ISTable();
  FillTestTable(ss1);
  FillTestTable(ss2);
  FillTestTable(ss3);
  ss1->DeleteRow(3);
  ss1->DeleteRow(5);
  ss1->DeleteRow(7);
//  ss1->CompressTable();

  cout<<ss1->GetNumRows()<<endl;

  ss2->RemoveRow(3);
  ss2->RemoveRow(5);
  ss2->RemoveRow(7);
  cout<<ss2->GetNumRows()<<endl;

  ss1->PrintTable();
  ss2->PrintTable();

  delete ss1;
  delete ss2;
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



void PrintArray(ReVarPCifArray<int> * target) {
  if (target) {
    int len = target->Length();
    for (int i = 0; i < len - 1; i++) {
      cout << (*target)[i] << ", ";
    }
    if (len > 0)
      cout << (*target)[len - 1];
  }
  cout << endl;
}
