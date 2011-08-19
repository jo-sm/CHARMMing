#include "BTreeObj.h"
#include "CifString.h"

main() {
  int n;

  BTreeObj<CifString> btobj;
  n=0;

  btobj.Add("A");n++;
  btobj.Add("B");n++;
  btobj.Add("C");n++;
  btobj.Add("D");n++;
  btobj.Add("E");n++;
  btobj.Add("F");n++;
  btobj.Add("G");n++;
  btobj.Add("H");n++;
  btobj.Add("N");n++;
  btobj.Add("M");n++;
  btobj.Add("L");n++;
  btobj.Add("K");n++;
  btobj.Add("J");n++;
  btobj.Add("I");n++;

//  btobj.Add("");n++;

  btobj.PrintBTreeInorder();

  cout<<"btobj.Root()    = "<<btobj.Root()<<endl;
  cout<<"btobj.Current() = "<<btobj.Current()<<endl;
  btobj.GoToPrev();
  cout<<"********** Prev ***********"<<endl;
  cout<<"btobj.Root()    = "<<btobj.Root()<<endl;
  cout<<"btobj.Current() = "<<btobj.Current()<<endl;

  btobj.GoFirst();
  cout<<"********** First ***********"<<endl;
  cout<<"btobj.Root()    = "<<btobj.Root()<<endl;
  cout<<"btobj.Current() = "<<btobj.Current()<<endl;

  btobj.GoLast();
  cout<<"********** Last ***********"<<endl;
  cout<<"btobj.Root()    = "<<btobj.Root()<<endl;
  cout<<"btobj.Current() = "<<btobj.Current()<<endl;

  btobj.GoToNext();
  cout<<"********** Next ***********"<<endl;
  cout<<"btobj.Root()    = "<<btobj.Root()<<endl;
  cout<<"btobj.Current() = "<<btobj.Current()<<endl;
  btobj.Seek("A");
  btobj.ReplData("W");
  cout<<"btobj.Root()    = "<<btobj.Root()<<endl;
  cout<<"btobj.Current() = "<<btobj.Current()<<endl;
  btobj.Seek("H");
  btobj.Delete();
  btobj.PrintBTreeInorder();
  cout<<"btobj.Root()    = "<<btobj.Root()<<endl;
  cout<<"btobj.Current() = "<<btobj.Current()<<endl;
}
