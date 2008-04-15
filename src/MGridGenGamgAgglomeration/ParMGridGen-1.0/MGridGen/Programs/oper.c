/* To process the output files of mgridgen and onmetis */

#include <stdio.h>
#include <stdlib.h>

int MAXLINE = 100;
int main(int argc, char *argv[]) 
{

   char line1[MAXLINE+1], line2[MAXLINE+1], line[MAXLINE+MAXLINE+2];
   FILE *fpin1, *fpin2, *fpout;

   fpin1 = fopen(argv[1], "r");
   fpin2 = fopen(argv[2], "r");
   fpout = fopen(argv[3], "w");

   printf("%s %s %s\n", argv[1], argv[2], argv[3]);

   do {
     fgets(line1, MAXLINE, fpin1);
     fgets(line2, MAXLINE, fpin2);
     strcpy(line, line1);
     line[strlen(line1)-1]=' ';
     strcat(line, line2); 
     fputs(line, fpout);
   } while (line1 != NULL);

   fclose(fpin1);
   fclose(fpin2);
   fclose(fpout);

   return(0);
}

