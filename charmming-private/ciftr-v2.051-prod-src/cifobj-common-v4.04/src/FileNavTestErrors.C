
#include <iostream.h>
#include <malloc.h>
#include <stdlib.h>
#include <stdio.h>
#include "BlockIO.h"
#include "FileNavigator.h"
#include "FileNavigatorError.h"

//const int NUMSTR = 512;

void PrintString(char * theString) {
  if (!theString) return;
  char * target = theString;
  while (*target != '\0')
    putchar(*target++);
  fflush(stdout);
}

int main(int argc, char ** argv) {
  FileNavigator fnav;
  char * s;
  int err;
  uWord num;
  Word index;
  cout << "Opening null file for read: ";
  err = fnav.OpenFile(NULL, READ_MODE);
  fnav.PrintError(err);
  cout << endl;
  cout << "GetString: ";
  s = fnav.GetString(0, err);
  fnav.PrintError(err);
  cout << endl;
  cout << "Opening file \"" << argv[1] << "\" for read... ";
  err = fnav.OpenFile(argv[1], READ_MODE);
  fnav.PrintError(err);
  cout << endl;
  cout << "Reading file header... ";
  err = fnav.ReadFileHeader();
  fnav.PrintError(err);
  cout << endl;
  cout << "GetString: ";
  s = fnav.GetString(0, err);
  fnav.PrintError(err);
  cout << endl;
  cout << "Closing file... ";
  err = fnav.CloseFile();
  fnav.PrintError(err);
  cout << endl;

  cout << "Opening file \"" << argv[1] << "\" for write... ";
  err =  fnav.OpenFile(argv[1], WRITE_MODE);
  fnav.PrintError(err);
  cout << endl;
  cout << "Reading file header... ";
  err = fnav.ReadFileHeader();
  fnav.PrintError(err);
  cout << endl;
  char s2[13] = "Hello World!";
  Word w[10] = {1,2,3,4,5,6,7,8,9,10};
  err = fnav.WriteWords(w, 10, index);
  cout << "WriteWords: ";
  fnav.PrintError(err);
  cout << " Index = " << index;
  cout << endl;
  Word myIndex = index;
  err = fnav.WriteString(s2, index);
  cout << "WriteString: ";
  fnav.PrintError(err);
  cout << " Index = " << index << endl;
  cout << "UpdateString at index = 0 ... ";
  fnav.PrintError(fnav.UpdateString("Hi", 0, index));
  cout << endl;
  cout << " New Index = " << index << endl;
  cout << "UpdateString at index = 0 ... ";
  fnav.PrintError(fnav.UpdateString("Hello", 0, index));
  cout << endl;
  cout << " New Index = " << index << endl;
  err = fnav.WriteStrings(argv, argc, index);
  cout << "WriteStrings: ";
  fnav.PrintError(err);
  cout << " Index: " << index << endl;
  cout << "Closing file... ";
  fnav.PrintError(fnav.CloseFile());
  cout << endl;

  cout << "Reopening file to read... ";
  fnav.PrintError(fnav.OpenFile(argv[1], READ_MODE));
  cout << endl;
  cout << "Reading file header... ";
  fnav.PrintError(fnav.ReadFileHeader());
  cout << endl;
  Word * w2 = fnav.GetWords(myIndex, num, err);
  cout << "GetWords: ";
  fnav.PrintError(err);
  cout << "   # " << num << endl;
  int i;
  if (err == NO_ERROR) 
    for (i = 0; i < num; i++) {
      cout << w2[i] << endl;
    }
  free(w2);
  s = fnav.GetString(0, err);
  cout << "GetString at Index = 0: ";
  fnav.PrintError(err);
  if (err == NO_ERROR) cout << " : " << s << endl;
  else cout << endl;
  cout << "GetString at Index = " << index - 1 << ": "; 
  s = fnav.GetString(index - 1, err);
  fnav.PrintError(err);
  if (err == NO_ERROR) cout << " : " << s << endl;
  else cout << endl;

  char ** st = fnav.GetStrings(index, num, err);

  fnav.PrintError(err);
  cout << " : " << endl;
  if (err >= NO_ERROR)
    for (i = 0; i < num; i++) {
      PrintString(st[i]);
      cout << endl;
      free(st[i]);
    }
  cout << "Closing file... ";
  fnav.PrintError(fnav.CloseFile());
  cout << endl;

  Word nwords[4] = {55, 122, 73, 23};

  cout << "Reopening file for write... ";
  fnav.PrintError(fnav.OpenFile(argv[1], WRITE_MODE));
  cout << endl;
  cout << "Reading file header... ";
  fnav.PrintError(fnav.ReadFileHeader());
  cout << endl;
  cout << "Update: ";
  fnav.PrintError(fnav.UpdateWords(nwords, 4, myIndex, myIndex));
  cout << endl;
  cout << "Closing file... ";
  fnav.PrintError(fnav.CloseFile());
  cout << endl;

  cout << "Reopening file for write... ";
  fnav.PrintError(fnav.OpenFile(argv[1], READ_MODE));
  cout << endl;
  cout << "Reading file header... ";
  fnav.PrintError(fnav.ReadFileHeader());
  cout << endl;
  Word * g = fnav.GetWords(myIndex, num, err);
  cout << "GetWord: ";
  fnav.PrintError(err);
  cout << " at index: " << myIndex << " num: " << num << endl;
  if (err == NO_ERROR) 
    for (i = 0; i < num; i++) {
      cout << g[i] << endl;
    }
  free(g);
  cout << "Closing file... ";
  fnav.PrintError(fnav.CloseFile());
  cout << endl;
  exit(0);
}
