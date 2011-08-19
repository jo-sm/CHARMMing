/* 
    File: ISTest5.C
*/

#include <stdlib.h>
#include <iostream.h>
#include "ISTable.h"


// prototypes

void FillTestTable(ISTable *s);
void PrintArray(ReVarPCifArray<int> * target);

int main(int argc, char ** argv) {

  ISTable *ss;
  char fname[100];
  ReVarPCifArray<int> list;
  ReVarCifArray<CifString> list2;
  int errCode;
  ReVarPCifArray<int> list3;
  ReVarCifArray<CifString> list4;
  int recNo;
  ReVarPCifArray<int> *listOut;
  ReVarPCifArray<int> list5;
  ReVarCifArray<CifString> list6;
  Word ind;

  ss = new ISTable("pera");
  FillTestTable(ss);

  ss->SetFlags(ISTable::DT_INTEGER,0);
  ss->SetFlags(ISTable::DT_INTEGER,1);
  ss->SetFlags(ISTable::DT_INTEGER,2);

  list.Add(0);
  list.Add(1);

  list3.Add(1);
  list3.Add(2);

  ss->CreateIndex("index0",list);
  errCode=ss->CreateIndex("index1",list3);
  cout<<"CreateIndex = "<<errCode<<endl;

  ReVarCifArray<CifString>* newRow;
  newRow = new ReVarCifArray<CifString>;
  newRow->Add("49");
  newRow->Add("0");
  newRow->Add("new");

  ss->AddRow();
  ss->FillRow(*newRow,ss->GetLastRowIndex());

  ss->AddRow();
  ss->FillRow(*newRow,ss->GetLastRowIndex());
  list2.Add("48");
  list2.Add("1");
  list4.Add("0");
  list4.Add("new");

  ss->AddRow();
  cout<<"GetNumRows = "<<ss->GetNumRows()<<endl;
  ss->FillRow(*newRow,ss->GetLastRowIndex());

  list5.Add(1);
  list6.Add("0");


  recNo=ss->FindFirst(list6,list5,errCode);
  cout<<"recNo = "<<recNo<<"     errCode ="<<errCode<<endl;

  ss->DeleteRow(50);
  ss->CompressTable();


  listOut=ss->Search(list2,list,errCode);
  cout<<"     errCode ="<<errCode<<endl;
  PrintArray(listOut);

  delete listOut;
  listOut=ss->SearchLessThan("index1",list6,list5,errCode);
  cout<<"     errCode ="<<errCode<<endl;
  PrintArray(listOut);

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
  ISTable gnu,gnu2;
  fnav->OpenFile("./test/outfile.db", READ_MODE,0);
  fnav->ReadFileHeader();
  ind=0;
  gnu.GetObject(ind, fnav);
  recNo=gnu.FindFirst("index0",list2,list,errCode);
  cout<<"recNo = "<<recNo<<"     errCode ="<<errCode<<endl;
  recNo=gnu.FindFirst("index1",list4,list3,errCode);
  cout<<"recNo = "<<recNo<<"     errCode ="<<errCode<<endl;
//  gnu.PrintTable("index0");

  delete listOut;
  delete newRow;
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
 for (i=0; i<50; i++) {
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
