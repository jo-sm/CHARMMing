#ifndef ManualHeadings
#define ManualHeadings

void HorizontalRule(FILE *outp, int n);

void PrintSplash(FILE *outp);

void PrintParagraph(char *vpar, int width, FILE *outp);

void PrintUsage();

void PrintCommandLineInputOptions();

void PrintInputFormat();

void PrintFilesNamelistVariables();

void PrintCntrlNamelistVariables();

void PrintEwaldNamelistVariables();

void PrintForceNamelistVariables();

void PrintFitNamelistVariables();

void PrintAttributions();

#endif
