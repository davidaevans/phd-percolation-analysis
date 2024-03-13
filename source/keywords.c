/*--------------------------------------------------------------------------------

  KEYWORDS LIBRARY
  v 1.0
  MM 7.vi.99

  A set of functions to facilitate keyword driven input files for passing
  instructions to a C program. The following functions are provided.

  int read_line(FILE *)
     Reads a line from the (already open) stream. Returns 1 if a line is read
     successfully, or 0 if not (e.g. if end-of-file is reached). Empty lines and
     lines containing only white space are ignored. The maximum length of line
     that can be read is MAXLENGTH-1. If the length of a line exceeds this limit,
     the remainder of the line is discarded.
     The remaining functions are used to process the contents of the most recently
     read line.

  int get_string(char [], int)
     Reads a string from the line. A string is any sequence of characters delimited
     by space, tab, or carriage return. The integer argument is the number of
     elements in the first argument. If the length of the string being read equals
     or exceeds this number, the remainder is discarded. The function returns the
     length of the string read in, or 0 if none could be read.

  void upper_case(char [])
     Returns the string converted to upper case.

  int get_int(long int *)
     Reads an integer from the line. Returns 1 if successful or 0 if not.

  int get_double(double *)
     Reads a double precision float from the line. Returns 1 if successful or 0 if
     not.

  The header file keywords.h contains prototypes for these functions.

  --------------------------------------------------------------------------------*/

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXLENGTH 500      /* Maximum length of line to read in. */
#define MAXITEM 100        /* Maximum length of a numerical iten to process. */

char input[MAXLENGTH];     /* Current line being read or processed. */
char *pos;                 /* Pointer to current position in input line. */


int read_line(FILE *infile)
{
   int again, c, i;

   again = 1;
   do {
      for(i=0; (c=fgetc(infile))!='\n' && c!=EOF && i<=MAXLENGTH-2; i++) {
         input[i] = c;
         if (!isspace(c)) {
            again = 0;
         }
      }
      input[i] = '\0';
   } while(again && c!=EOF);   /* Repeat if line was empty or contained only blanks. */

   for(; c!='\n' && c!=EOF; c=fgetc(infile))
      ;   /* Read to end of line even if input string is not long enough to hold it. */

   pos = input;

   return (again ? 0 : 1);
}


int get_string(char string[], int maxlen)
{
   int pstring;   /* Position in string. */

   for(; *pos==' ' || *pos=='\t'; pos++)   /* Read past empty space. */
      ;

   if(*pos=='\n' || *pos=='\0') {
      return 0;
   } else {
      pstring = 0;
      for(; !isspace(*pos) && *pos!='\0' && pstring<=maxlen-2; pos++, pstring++) {
         string[pstring] = *pos;
      }
      string[pstring] = '\0';
      /* Move to end of string if longer than maxlen-1. */
      if(pstring==maxlen-1) {
         for(; !isspace(*pos) && *pos!='\0'; pos++)
            ;
      }
      return pstring;
   }
}


int get_int(long int *value)
{
   char item[MAXITEM];

   if(get_string(item, MAXITEM)) {
      return sscanf(item, "%ld", value);
   } else {
      return 0;
   }
}


int get_double(double *value)
{
   char item[MAXITEM];

   if(get_string(item, MAXITEM)) {
      return sscanf(item, "%le", value);
   } else {
      return 0;
   }
}


void upper_case(char string[])
{
   while ((*string = toupper(*string)) != '\0') string++;
}
