#include "TblIndexObj.h"
#include "CifString.h"

main() {
  int n;
  int ret;


  TblIndexObj btobj;
  n=0;
  btobj.Add("WW");n++;
  btobj.Add("BB");n++;
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
  ret = btobj.Seek("WW");
  cout<<"ret = "<<ret<<endl;
  btobj.ReplData("W");
  cout<<"btobj.Root()    = "<<btobj.Root()<<endl;
  cout<<"btobj.Current() = "<<btobj.Current()<<endl;
  btobj.Seek("H");
  btobj.Delete();
  cout<<"btobj.Root()    = "<<btobj.Root()<<endl;
  cout<<"btobj.Current() = "<<btobj.Current()<<endl;

  btobj.PrintBTreeInorder();
}
