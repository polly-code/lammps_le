/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://lammps.sandia.gov/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "utils.h"

#include "comm.h"
#include "compute.h"
#include "error.h"
#include "fix.h"
#include "memory.h"
#include "modify.h"
#include "text_file_reader.h"
#include "tokenizer.h"
#include "update.h"

#include <cctype>
#include <cerrno>
#include <cstring>

#if defined(__linux__)
#include <unistd.h>  // for readlink
#endif

/*! \file utils.cpp */

/*
 * Mini regex-module adapted from https://github.com/kokke/tiny-regex-c
 * which is in the public domain.
 *
 * Supports:
 * ---------
 *   '.'        Dot, matches any character
 *   '^'        Start anchor, matches beginning of string
 *   '$'        End anchor, matches end of string
 *   '*'        Asterisk, match zero or more (greedy)
 *   '+'        Plus, match one or more (greedy)
 *   '?'        Question, match zero or one (non-greedy)
 *   '[abc]'    Character class, match if one of {'a', 'b', 'c'}
 *   '[a-zA-Z]' Character ranges, the character set of the ranges { a-z | A-Z }
 *   '\s'       Whitespace, \t \f \r \n \v and spaces
 *   '\S'       Non-whitespace
 *   '\w'       Alphanumeric, [a-zA-Z0-9_]
 *   '\W'       Non-alphanumeric
 *   '\d'       Digits, [0-9]
 *   '\D'       Non-digits
 *   '\i'       Integer chars, [0-9], '+' and '-'
 *   '\I'       Non-integers
 *   '\f'       Floating point number chars, [0-9], '.', 'e', 'E', '+' and '-'
 *   '\F'       Non-floats
 *
 * *NOT* supported:
 *   '[^abc]'   Inverted class
 *   'a|b'      Branches
 *   '(abc)+'   Groups
 */

extern "C"
{
  /** Match text against a (simplified) regular expression
   * (regexp will be compiled automatically). */
  static int  re_match(const char *text, const char *pattern);
}

////////////////////////////////////////////////////////////////////////
// Merge sort support functions

static void do_merge(int *idx, int *buf, int llo, int lhi, int rlo, int rhi,
                     void *ptr, int (*comp)(int, int, void *));
static void insertion_sort(int *index, int num, void *ptr,
                           int (*comp)(int, int, void*));

////////////////////////////////////////////////////////////////////////

using namespace LAMMPS_NS;

/** More flexible and specific matching of a string against a pattern.
 *  This function is supposed to be a more safe, more specific and
 *  simple to use API to find pattern matches. The purpose is to replace
 *  uses of either strncmp() or strstr() in the code base to find
 *  sub-strings safely. With strncmp() finding prefixes, the number of
 *  characters to match must be counted, which can lead to errors,
 *  while using "^pattern" will do the same with less problems.
 *  Matching for suffixes using strstr() is not as specific as 'pattern$',
 *  and complex matches, e.g. "^rigid.*\/small.*", to match all small
 *  body optimized rigid fixes require only one test.
 *
 *  The use of std::string arguments allows for simple concatenation
 *  even with char * type variables.
 *  Example: utils::strmatch(text, std::string("^") + charptr)
 */
bool utils::strmatch(const std::string &text, const std::string &pattern)
{
  const int pos = re_match(text.c_str(),pattern.c_str());
  return (pos >= 0);
}

/** This function simplifies the repetitive task of outputting some
 * message to both the screen and/or the log file. In combination
 * with using fmt::format(), which returns the formatted text
 * in a std::string() instance, this can be used to reduce
 * operations previously requiring several lines of code to
 * a single statement. */

void utils::logmesg(LAMMPS *lmp, const std::string &mesg)
{
  if (lmp->screen)  fputs(mesg.c_str(), lmp->screen);
  if (lmp->logfile) fputs(mesg.c_str(), lmp->logfile);
}

/* define this here, so we won't have to include the headers
   everywhere and utils.h will more likely be included anyway. */

std::string utils::getsyserror()
{
  return std::string(strerror(errno));
}

/** On Linux the folder /proc/self/fd holds symbolic links to the actual
 * pathnames associated with each open file descriptor of the current process.
 *
 * This function is used to provide a filename with error messages in functions
 * where the filename is not passed as an argument, but the FILE * pointer.
 */
const char *utils::guesspath(char *buf, int len, FILE *fp)
{
  memset(buf,0,len);

#if defined(__linux__)
  int fd = fileno(fp);
  // get pathname from /proc or copy (unknown)
  if (readlink(fmt::format("/proc/self/fd/{}",fd).c_str(),buf,len-1) <= 0)
    strncpy(buf,"(unknown)",len-1);
#else
  strncpy(buf,"(unknown)",len-1);
#endif
  return buf;
}

#define MAXPATHLENBUF 1024
/* like fgets() but aborts with an error or EOF is encountered */
void utils::sfgets(const char *srcname, int srcline, char *s, int size,
                   FILE *fp, const char *filename, Error *error)
{
  char *rv = fgets(s,size,fp);
  if (rv == nullptr) { // something went wrong
    char buf[MAXPATHLENBUF];
    std::string errmsg;

    // try to figure out the file name from the file pointer
    if (!filename)
      filename = guesspath(buf,MAXPATHLENBUF,fp);

    if (feof(fp)) {
      errmsg = "Unexpected end of file while reading file '";
    } else if (ferror(fp)) {
      errmsg = "Unexpected error while reading file '";
    } else {
      errmsg = "Unexpected short read while reading file '";
    }
    errmsg += filename;
    errmsg += "'";

    if (error) error->one(srcname,srcline,errmsg);
    if (s) *s = '\0'; // truncate string to empty in case error is null pointer
  }
  return;
}

/* like fread() but aborts with an error or EOF is encountered */
void utils::sfread(const char *srcname, int srcline, void *s, size_t size,
                   size_t num, FILE *fp, const char *filename, Error *error)
{
  size_t rv = fread(s,size,num,fp);
  if (rv != num) { // something went wrong
    char buf[MAXPATHLENBUF];
    std::string errmsg;

    // try to figure out the file name from the file pointer
    if (!filename)
      filename = guesspath(buf,MAXPATHLENBUF,fp);

    if (feof(fp)) {
      errmsg = "Unexpected end of file while reading file '";
    } else if (ferror(fp)) {
      errmsg = "Unexpected error while reading file '";
    } else {
      errmsg = "Unexpected short read while reading file '";
    }
    errmsg += filename;
    errmsg += "'";

    if (error) error->one(srcname,srcline,errmsg);
  }
  return;
}

/* ------------------------------------------------------------------ */

std::string utils::check_packages_for_style(const std::string &style,
                                            const std::string &name,
                                            LAMMPS *lmp)
{
  std::string errmsg = "Unrecognized " + style + " style '" + name + "'";
  const char *pkg = lmp->match_style(style.c_str(),name.c_str());

  if (pkg) {
    errmsg += fmt::format(" is part of the {} package",pkg);
    if (lmp->is_installed_pkg(pkg))
      errmsg += ", but seems to be missing because of a dependency";
    else
      errmsg += " which is not enabled in this LAMMPS binary.";
  }
  return errmsg;
}


/* ----------------------------------------------------------------------
   read a floating point value from a string
   generate an error if not a legitimate floating point value
   called by various commands to check validity of their arguments
------------------------------------------------------------------------- */

double utils::numeric(const char *file, int line, const char *str,
                      bool do_abort, LAMMPS *lmp)
{
  int n = 0;

  if (str) n = strlen(str);
  if (n == 0) {
    if (do_abort)
      lmp->error->one(file,line,"Expected floating point parameter instead of"
                      " NULL or empty string in input script or data file");
    else
      lmp->error->all(file,line,"Expected floating point parameter instead of"
                      " NULL or empty string in input script or data file");
  }

  for (int i = 0; i < n; i++) {
    if (isdigit(str[i])) continue;
    if (str[i] == '-' || str[i] == '+' || str[i] == '.') continue;
    if (str[i] == 'e' || str[i] == 'E') continue;
    std::string msg("Expected floating point parameter instead of '");
    msg += str;
    msg += "' in input script or data file";
    if (do_abort)
      lmp->error->one(file,line,msg);
    else
      lmp->error->all(file,line,msg);
  }

  return atof(str);
}

/* ----------------------------------------------------------------------
   read an integer value from a string
   generate an error if not a legitimate integer value
   called by various commands to check validity of their arguments
------------------------------------------------------------------------- */

int utils::inumeric(const char *file, int line, const char *str,
                    bool do_abort, LAMMPS *lmp)
{
  int n = 0;

  if (str) n = strlen(str);
  if (n == 0) {
    if (do_abort)
      lmp->error->one(file,line,"Expected integer parameter instead of "
                      "NULL or empty string in input script or data file");
    else
      lmp->error->all(file,line,"Expected integer parameter instead of "
                      "NULL or empty string in input script or data file");
  }

  for (int i = 0; i < n; i++) {
    if (isdigit(str[i]) || str[i] == '-' || str[i] == '+') continue;
    std::string msg("Expected integer parameter instead of '");
    msg += str;
    msg += "' in input script or data file";
    if (do_abort)
      lmp->error->one(file,line,msg);
    else
      lmp->error->all(file,line,msg);
  }

  return atoi(str);
}

/* ----------------------------------------------------------------------
   read a big integer value from a string
   generate an error if not a legitimate integer value
   called by various commands to check validity of their arguments
------------------------------------------------------------------------- */

bigint utils::bnumeric(const char *file, int line, const char *str,
                       bool do_abort, LAMMPS *lmp)
{
  int n = 0;

  if (str) n = strlen(str);
  if (n == 0) {
    if (do_abort)
      lmp->error->one(file,line,"Expected integer parameter instead of "
                      "NULL or empty string in input script or data file");
    else
      lmp->error->all(file,line,"Expected integer parameter instead of "
                      "NULL or empty string in input script or data file");
  }

  for (int i = 0; i < n; i++) {
    if (isdigit(str[i]) || str[i] == '-' || str[i] == '+') continue;
    std::string msg("Expected integer parameter instead of '");
    msg += str;
    msg += "' in input script or data file";
    if (do_abort)
      lmp->error->one(file,line,msg);
    else
      lmp->error->all(file,line,msg);
  }

  return ATOBIGINT(str);
}

/* ----------------------------------------------------------------------
   read a tag integer value from a string
   generate an error if not a legitimate integer value
   called by various commands to check validity of their arguments
------------------------------------------------------------------------- */

tagint utils::tnumeric(const char *file, int line, const char *str,
                       bool do_abort, LAMMPS *lmp)
{
  int n = 0;

  if (str) n = strlen(str);
  if (n == 0) {
    if (do_abort)
      lmp->error->one(file,line,"Expected integer parameter instead of "
                      "NULL or empty string in input script or data file");
    else
      lmp->error->all(file,line,"Expected integer parameter instead of "
                      "NULL or empty string in input script or data file");
  }

  for (int i = 0; i < n; i++) {
    if (isdigit(str[i]) || str[i] == '-' || str[i] == '+') continue;
    std::string msg("Expected integer parameter instead of '");
    msg += str;
    msg += "' in input script or data file";
    if (do_abort)
      lmp->error->one(file,line,msg);
    else
      lmp->error->all(file,line,msg);
  }

  return ATOTAGINT(str);
}

/* ----------------------------------------------------------------------
   compute bounds implied by numeric str with a possible wildcard asterisk
------------------------------------------------------------------------- */
template<typename TYPE>
void utils::bounds(const char *file, int line, const std::string &str,
                   bigint nmin, bigint nmax, TYPE &nlo, TYPE &nhi, Error *error)
{
  size_t found = str.find_first_of("*");

  nlo = nhi = -1;
  if (found == std::string::npos) {    // contains no '*'
    nlo = nhi = strtol(str.c_str(),nullptr,10);
  } else if (str.size() == 1) {        // is only '*'
    nlo = nmin;
    nhi = nmax;
  } else if (found == 0) {             // is '*j'
    nlo = nmin;
    nhi = strtol(str.substr(1).c_str(),nullptr,10);
  } else if (str.size() == found+1) {  // is 'i*'
    nlo = strtol(str.c_str(),nullptr,10);
    nhi = nmax;
  } else {                             // is 'i*j'
    nlo = strtol(str.c_str(),nullptr,10);
    nhi = strtol(str.substr(found+1).c_str(),nullptr,10);
  }

  if (error) {
    if (nlo < nmin)
      error->all(file,line,fmt::format("Numeric index {} is out of bounds"
                                       "({}-{})",nlo,nmin,nmax));
    else if (nhi > nmax)
      error->all(file,line,fmt::format("Numeric index {} is out of bounds"
                                       "({}-{})",nhi,nmin,nmax));
    else if (nlo > nhi)
      error->all(file,line,fmt::format("Numeric index {} is out of bounds"
                                       "({}-{})",nlo,nmin,nhi));
  }
}

template void utils::bounds<>(const char *, int, const std::string &,
                              bigint, bigint, int &, int &, Error *);
template void utils::bounds<>(const char *, int, const std::string &,
                              bigint, bigint, long &, long &, Error *);
template void utils::bounds<>(const char *, int, const std::string &,
                              bigint, bigint, long long &, long long &, Error *);

/* -------------------------------------------------------------------------
   Expand list of arguments in arg to earg if arg contains wildcards
------------------------------------------------------------------------- */

int utils::expand_args(const char *file, int line, int narg, char **arg,
                       int mode, char **&earg, LAMMPS *lmp)
{
  int n,iarg,index,nlo,nhi,nmax,expandflag,icompute,ifix;
  char *ptr1,*ptr2,*str;

  ptr1 = nullptr;
  for (iarg = 0; iarg < narg; iarg++) {
    ptr1 = strchr(arg[iarg],'*');
    if (ptr1) break;
  }

  if (!ptr1) {
    earg = arg;
    return narg;
  }

  // maxarg should always end up equal to newarg, so caller can free earg

  int maxarg = narg-iarg;
  earg = (char **) lmp->memory->smalloc(maxarg*sizeof(char *),"input:earg");

  int newarg = 0;
  for (iarg = 0; iarg < narg; iarg++) {
    expandflag = 0;

    if (strncmp(arg[iarg],"c_",2) == 0 ||
        strncmp(arg[iarg],"f_",2) == 0) {

      ptr1 = strchr(&arg[iarg][2],'[');
      if (ptr1) {
        ptr2 = strchr(ptr1,']');
        if (ptr2) {
          *ptr2 = '\0';
          if (strchr(ptr1,'*')) {
            if (arg[iarg][0] == 'c') {
              *ptr1 = '\0';
              icompute = lmp->modify->find_compute(&arg[iarg][2]);
              *ptr1 = '[';

              // check for global vector/array, peratom array, local array

              if (icompute >= 0) {
                if (mode == 0 && lmp->modify->compute[icompute]->vector_flag) {
                  nmax = lmp->modify->compute[icompute]->size_vector;
                  expandflag = 1;
                } else if (mode == 1 && lmp->modify->compute[icompute]->array_flag) {
                  nmax = lmp->modify->compute[icompute]->size_array_cols;
                  expandflag = 1;
                } else if (lmp->modify->compute[icompute]->peratom_flag &&
                           lmp->modify->compute[icompute]->size_peratom_cols) {
                  nmax = lmp->modify->compute[icompute]->size_peratom_cols;
                  expandflag = 1;
                } else if (lmp->modify->compute[icompute]->local_flag &&
                           lmp->modify->compute[icompute]->size_local_cols) {
                  nmax = lmp->modify->compute[icompute]->size_local_cols;
                  expandflag = 1;
                }
              }
            } else if (arg[iarg][0] == 'f') {
              *ptr1 = '\0';
              ifix = lmp->modify->find_fix(&arg[iarg][2]);
              *ptr1 = '[';

              // check for global vector/array, peratom array, local array

              if (ifix >= 0) {
                if (mode == 0 && lmp->modify->fix[ifix]->vector_flag) {
                  nmax = lmp->modify->fix[ifix]->size_vector;
                  expandflag = 1;
                } else if (mode == 1 && lmp->modify->fix[ifix]->array_flag) {
                  nmax = lmp->modify->fix[ifix]->size_array_cols;
                  expandflag = 1;
                } else if (lmp->modify->fix[ifix]->peratom_flag &&
                           lmp->modify->fix[ifix]->size_peratom_cols) {
                  nmax = lmp->modify->fix[ifix]->size_peratom_cols;
                  expandflag = 1;
                } else if (lmp->modify->fix[ifix]->local_flag &&
                           lmp->modify->fix[ifix]->size_local_cols) {
                  nmax = lmp->modify->fix[ifix]->size_local_cols;
                  expandflag = 1;
                }
              }
            }
          }
          *ptr2 = ']';
        }
      }
    }

    if (expandflag) {
      *ptr2 = '\0';
      bounds(file,line,ptr1+1,1,nmax,nlo,nhi,lmp->error);
      *ptr2 = ']';
      if (newarg+nhi-nlo+1 > maxarg) {
        maxarg += nhi-nlo+1;
        earg = (char **)
          lmp->memory->srealloc(earg,maxarg*sizeof(char *),"input:earg");
      }
      for (index = nlo; index <= nhi; index++) {
        n = strlen(arg[iarg]) + 16;   // 16 = space for large inserted integer
        str = earg[newarg] = new char[n];
        strncpy(str,arg[iarg],ptr1+1-arg[iarg]);
        sprintf(&str[ptr1+1-arg[iarg]],"%d",index);
        strcat(str,ptr2);
        newarg++;
      }

    } else {
      if (newarg == maxarg) {
        maxarg++;
        earg = (char **)
          lmp->memory->srealloc(earg,maxarg*sizeof(char *),"input:earg");
      }
      n = strlen(arg[iarg]) + 1;
      earg[newarg] = new char[n];
      strcpy(earg[newarg],arg[iarg]);
      newarg++;
    }
  }

  //printf("NEWARG %d\n",newarg);
  //for (int i = 0; i < newarg; i++)
  //  printf("  arg %d: %s\n",i,earg[i]);

  return newarg;
}

/* ----------------------------------------------------------------------
   Return string without leading or trailing whitespace
------------------------------------------------------------------------- */

std::string utils::trim(const std::string &line) {
  int beg = re_match(line.c_str(),"\\S+");
  int end = re_match(line.c_str(),"\\s+$");
  if (beg < 0) beg = 0;
  if (end < 0) end = line.size();

  return line.substr(beg,end-beg);
}

/* ----------------------------------------------------------------------
   Return string without trailing # comment
------------------------------------------------------------------------- */

std::string utils::trim_comment(const std::string &line) {
  auto end = line.find_first_of("#");
  if (end != std::string::npos) {
    return line.substr(0, end);
  }
  return std::string(line);
}

/* ----------------------------------------------------------------------
   return number of words
------------------------------------------------------------------------- */

size_t utils::count_words(const char *text) {
  size_t count = 0;
  const char * buf = text;
  char c = *buf;

  while (c) {
    if (c == ' ' || c == '\t' || c == '\r' ||  c == '\n' || c == '\f') {
      c = *++buf;
      continue;
    };

    ++count;
    c = *++buf;

    while (c) {
      if (c == ' ' || c == '\t' || c == '\r' || c == '\n' || c == '\f') {
        break;
      }
      c = *++buf;
    }
  }

  return count;
}

/* ----------------------------------------------------------------------
   return number of words
------------------------------------------------------------------------- */

size_t utils::count_words(const std::string &text) {
  return utils::count_words(text.c_str());
}

/* ----------------------------------------------------------------------
   Return number of words
------------------------------------------------------------------------- */

size_t utils::count_words(const std::string &text, const std::string &separators) {
  size_t count = 0;
  size_t start = text.find_first_not_of(separators);

  while (start != std::string::npos) {
    size_t end = text.find_first_of(separators, start);
    ++count;

    if(end == std::string::npos) {
      return count;
    } else {
      start = text.find_first_not_of(separators, end + 1);
    }
  }
  return count;
}

/* ----------------------------------------------------------------------
   Trim comment from string and return number of words
------------------------------------------------------------------------- */

size_t utils::trim_and_count_words(const std::string &text, const std::string &separators) {
  return utils::count_words(utils::trim_comment(text), separators);
}

/* ----------------------------------------------------------------------
   Convert string into words on whitespace while handling single and
   double quotes.
------------------------------------------------------------------------- */
std::vector<std::string> utils::split_words(const std::string &text)
{
  std::vector<std::string> list;
  const char *buf = text.c_str();
  std::size_t beg = 0;
  std::size_t len = 0;
  std::size_t add = 0;
  char c = *buf;

  while (c) {
    // leading whitespace
    if (c == ' ' || c == '\t' || c == '\r' ||  c == '\n' || c == '\f') {
      c = *++buf;
      ++beg;
      continue;
    };
    len = 0;

    // handle escaped/quoted text.
    quoted:

    // handle single quote
    if (c == '\'') {
      ++beg;
      add = 1;
      c = *++buf;
      while (((c != '\'') && (c != '\0'))
             || ((c == '\\') && (buf[1] == '\''))) {
        if ((c == '\\') && (buf[1] == '\'')) {
          ++buf;
          ++len;
        }
        c = *++buf;
        ++len;
      }
      if (c != '\'') ++len;
      c = *++buf;

      // handle double quote
    } else if (c == '"') {
      ++beg;
      add = 1;
      c = *++buf;
      while (((c != '"') && (c != '\0'))
             || ((c == '\\') && (buf[1] == '"'))) {
        if ((c == '\\') && (buf[1] == '"')) {
          ++buf;
          ++len;
        }
        c = *++buf;
        ++len;
      }
      if (c != '"') ++len;
      c = *++buf;
    }

    // unquoted
    while (1) {
      if ((c == '\'') || (c == '"')) goto quoted;
      // skip escaped quote
      if ((c == '\\') && ((buf[1] == '\'') || (buf[1] == '"'))) {
        ++buf;
        ++len;
        c = *++buf;
        ++len;
      }
      if ((c == ' ') || (c == '\t') || (c == '\r') || (c == '\n')
          || (c == '\f') || (c == '\0')) {
          list.push_back(text.substr(beg,len));
          beg += len + add;
          break;
      }
      c = *++buf;
      ++len;
    }
  }
  return list;
}

/* ----------------------------------------------------------------------
   Return whether string is a valid integer number
------------------------------------------------------------------------- */

bool utils::is_integer(const std::string &str) {
  if (str.size() == 0) {
    return false;
  }

  for (auto c : str) {
    if (isdigit(c) || c == '-' || c == '+') continue;
    return false;
  }
  return true;
}

/* ----------------------------------------------------------------------
   Return whether string is a valid floating-point number
------------------------------------------------------------------------- */

bool utils::is_double(const std::string &str) {
  if (str.size() == 0) {
    return false;
  }

  for (auto c : str) {
    if (isdigit(c)) continue;
    if (c == '-' || c == '+' || c == '.') continue;
    if (c == 'e' || c == 'E') continue;
    return false;
  }
  return true;
}

/* ----------------------------------------------------------------------
   strip off leading part of path, return just the filename
------------------------------------------------------------------------- */

std::string utils::path_basename(const std::string &path) {
#if defined(_WIN32)
  size_t start = path.find_last_of("/\\");
#else
  size_t start = path.find_last_of("/");
#endif

  if (start == std::string::npos) {
    start = 0;
  } else {
    start += 1;
  }

  return path.substr(start);
}

/* ----------------------------------------------------------------------
   Return only the leading part of a path, return just the directory
------------------------------------------------------------------------- */

std::string utils::path_dirname(const std::string &path) {
#if defined(_WIN32)
  size_t start = path.find_last_of("/\\");
#else
  size_t start = path.find_last_of("/");
#endif

  if (start == std::string::npos) return ".";

  return path.substr(0,start);
}

/* ----------------------------------------------------------------------
   join two paths
------------------------------------------------------------------------- */

std::string utils::path_join(const std::string &a, const std::string &b) {
  #if defined(_WIN32)
    return fmt::format("{}\\{}", a, b);
  #else
    return fmt::format("{}/{}", a, b);
  #endif
}

/* ----------------------------------------------------------------------
   try to open file for reading
------------------------------------------------------------------------- */

bool utils::file_is_readable(const std::string &path) {
  FILE * fp = fopen(path.c_str(), "r");
  if(fp) {
    fclose(fp);
    return true;
  }
  return false;
}

/* ----------------------------------------------------------------------
   try to find potential file as specified by name
   search current directory and the LAMMPS_POTENTIALS directory if
   specified
------------------------------------------------------------------------- */
#if defined(_WIN32)
#define OS_PATH_VAR_SEP ";"
#else
#define OS_PATH_VAR_SEP ":"
#endif

std::string utils::get_potential_file_path(const std::string &path) {
  std::string filepath = path;
  std::string filename = utils::path_basename(path);

  if(utils::file_is_readable(filepath)) {
    return filepath;
  } else {
    // try the environment variable directory
    const char *var = getenv("LAMMPS_POTENTIALS");

    if (var != nullptr){
      Tokenizer dirs(var,OS_PATH_VAR_SEP);

      while (dirs.has_next()) {
        auto pot = utils::path_basename(filepath);
        auto path = dirs.next();
        filepath = utils::path_join(path, pot);

        if (utils::file_is_readable(filepath)) {
          return filepath;
        }
      }
    }
  }
  return "";
}
#undef OS_PATH_VAR_SEP

/* ----------------------------------------------------------------------
   read first line of potential file
   if it has a DATE field, return the following word
------------------------------------------------------------------------- */

std::string utils::get_potential_date(const std::string &path, const std::string &potential_name) {
  TextFileReader reader(path, potential_name);
  reader.ignore_comments = false;

  char *line = reader.next_line();
  Tokenizer words(line);
  while (words.has_next()) {
    if (words.next() == "DATE:") {
      if (words.has_next()) return words.next();
    }
  }
  return "";
}

/* ----------------------------------------------------------------------
   read first line of potential file
   if it has UNITS field, return following word
------------------------------------------------------------------------- */

std::string utils::get_potential_units(const std::string &path, const std::string &potential_name) {
  TextFileReader reader(path, potential_name);
  reader.ignore_comments = false;

  char *line = reader.next_line();
  Tokenizer words(line);
  while (words.has_next()) {
    if (words.next() == "UNITS:") {
      if (words.has_next()) return words.next();
    }
  }
  return "";
}

/* ----------------------------------------------------------------------
   return bitmask of supported conversions for a given property
------------------------------------------------------------------------- */
int utils::get_supported_conversions(const int property)
{
  if (property == ENERGY) {
    return METAL2REAL | REAL2METAL;
  }
  return NOCONVERT;
}

/* ----------------------------------------------------------------------
   return conversion factor for a given property and conversion setting
   return 0.0 if unknown.
------------------------------------------------------------------------- */

double utils::get_conversion_factor(const int property, const int conversion)
{
  if (property == ENERGY) {
    if (conversion == NOCONVERT) {
      return 1.0;
    } else if (conversion == METAL2REAL) {
      return 23.060549;
    } else if (conversion == REAL2METAL) {
      return 1.0/23.060549;
    }
  }
  return 0.0;
}

/* ----------------------------------------------------------------------
   open a potential file as specified by name
   if fails, search in dir specified by env variable LAMMPS_POTENTIALS
------------------------------------------------------------------------- */

FILE *utils::open_potential(const std::string &name, LAMMPS *lmp,
                            int *auto_convert)
{
  auto error = lmp->error;
  auto me = lmp->comm->me;

  std::string filepath = get_potential_file_path(name);

  if(!filepath.empty()) {
    std::string unit_style = lmp->update->unit_style;
    std::string date       = get_potential_date(filepath, "potential");
    std::string units      = get_potential_units(filepath, "potential");

    if(!date.empty() && (me == 0)) {
      logmesg(lmp, fmt::format("Reading potential file {} "
                               "with DATE: {}\n", name, date));
    }

    if (auto_convert == nullptr) {
      if (!units.empty() && (units != unit_style) && (me == 0)) {
        error->one(FLERR, fmt::format("Potential file {} requires {} units "
                                      "but {} units are in use", name, units,
                                      unit_style));
        return nullptr;
      }
    } else {
      if (units.empty() || units == unit_style) {
        *auto_convert = NOCONVERT;
      } else {
        if ((units == "metal") && (unit_style == "real")
            && (*auto_convert & METAL2REAL)) {
          *auto_convert = METAL2REAL;
        } else if ((units == "real") && (unit_style == "metal")
            && (*auto_convert & REAL2METAL)) {
          *auto_convert = REAL2METAL;
        } else {
          error->one(FLERR, fmt::format("Potential file {} requires {} units "
                                        "but {} units are in use", name,
                                        units, unit_style));
          return nullptr;
        }
      }
      if ((*auto_convert != NOCONVERT) && (me == 0))
        error->warning(FLERR, fmt::format("Converting potential file in "
                                          "{} units to {} units",
                                          units, unit_style));
    }
    return fopen(filepath.c_str(), "r");
  }
  return nullptr;
}

/* ----------------------------------------------------------------------
   convert a timespec ([[HH:]MM:]SS) to seconds
   the strings "off" and "unlimited" result in -1.0;
------------------------------------------------------------------------- */

double utils::timespec2seconds(const std::string &timespec)
{
  double vals[3];
  int i = 0;

  // first handle allowed textual inputs
  if (timespec == "off") return -1.0;
  if (timespec == "unlimited") return -1.0;

  vals[0] = vals[1] = vals[2] = 0;

  ValueTokenizer values(timespec, ":");

  try {
    for (i = 0; i < 3; i++) {
      if (!values.has_next()) break;
      vals[i] = values.next_int();
    }
  } catch (TokenizerException &e) {
    return -1.0;
  }

  if (i == 3) return (vals[0]*60 + vals[1])*60 + vals[2];
  else if (i == 2) return vals[0]*60 + vals[1];
  return vals[0];
}

/* ----------------------------------------------------------------------
   convert a LAMMPS version date (1Jan01) to a number
------------------------------------------------------------------------- */

int utils::date2num(const std::string &date)
{
  std::size_t found = date.find_first_not_of("0123456789 ");
  int num = strtol(date.substr(0,found).c_str(),nullptr,10);
  auto month = date.substr(found);
  found = month.find_first_of("0123456789 ");
  num += strtol(month.substr(found).c_str(),nullptr,10)*10000;
  if (num < 1000000) num += 20000000;

  if (strmatch(month,"^Jan")) num += 100;
  else if (strmatch(month,"^Feb")) num += 200;
  else if (strmatch(month,"^Mar")) num += 300;
  else if (strmatch(month,"^Apr")) num += 400;
  else if (strmatch(month,"^May")) num += 500;
  else if (strmatch(month,"^Jun")) num += 600;
  else if (strmatch(month,"^Jul")) num += 700;
  else if (strmatch(month,"^Aug")) num += 800;
  else if (strmatch(month,"^Sep")) num += 900;
  else if (strmatch(month,"^Oct")) num += 1000;
  else if (strmatch(month,"^Nov")) num += 1100;
  else if (strmatch(month,"^Dec")) num += 1200;
  return num;
}

/* ----------------------------------------------------------------------
 * Merge sort part 1: Loop over sublists doubling in size with each iteration.
 * Pre-sort small sublists with insertion sort for better overall performance.
------------------------------------------------------------------------- */

void utils::merge_sort(int *index, int num, void *ptr,
                       int (*comp)(int, int, void *))
{
  if (num < 2) return;

  int chunk,i,j;

  // do insertion sort on chunks of up to 64 elements

  chunk = 64;
  for (i=0; i < num; i += chunk) {
    j = (i+chunk > num) ? num-i : chunk;
    insertion_sort(index+i,j,ptr,comp);
  }

  // already done?

  if (chunk >= num) return;

  // continue with merge sort on the pre-sorted chunks.
  // we need an extra buffer for temporary storage and two
  // pointers to operate on, so we can swap the pointers
  // rather than copying to the hold buffer in each pass

  int *buf = new int[num];
  int *dest = index;
  int *hold = buf;

  while (chunk < num) {
    int m;

    // swap hold and destination buffer

    int *tmp = dest; dest = hold; hold = tmp;

    // merge from hold array to destination array

    for (i=0; i < num-1; i += 2*chunk) {
      j = i + 2*chunk;
      if (j > num) j=num;
      m = i+chunk;
      if (m > num) m=num;
      do_merge(dest,hold,i,m,m,j,ptr,comp);
    }

    // copy all indices not handled by the chunked merge sort loop

    for ( ; i < num ; i++ ) dest[i] = hold[i];
    chunk *= 2;
  }

  // if the final sorted data is in buf, copy back to index

  if (dest == buf) memcpy(index,buf,sizeof(int)*num);

  delete[] buf;
}

/* ------------------------------------------------------------------ */

/* ----------------------------------------------------------------------
 * Merge sort part 2: Insertion sort for pre-sorting of small chunks
------------------------------------------------------------------------- */

void insertion_sort(int *index, int num, void *ptr,
                           int (*comp)(int, int, void*))
{
  if (num < 2) return;
  for (int i=1; i < num; ++i) {
    int tmp = index[i];
    for (int j=i-1; j >= 0; --j) {
      if ((*comp)(index[j],tmp,ptr) > 0) {
        index[j+1] = index[j];
      } else {
        index[j+1] = tmp;
        break;
      }
      if (j == 0) index[0] = tmp;
    }
  }
}

/* ----------------------------------------------------------------------
 * Merge sort part 3: Merge two sublists
------------------------------------------------------------------------- */

static void do_merge(int *idx, int *buf, int llo, int lhi, int rlo, int rhi,
                     void *ptr, int (*comp)(int, int, void *))
{
  int i = llo;
  int l = llo;
  int r = rlo;
  while ((l < lhi) && (r < rhi)) {
    if ((*comp)(buf[l],buf[r],ptr) < 0)
      idx[i++] = buf[l++];
    else idx[i++] = buf[r++];
  }

  while (l < lhi) idx[i++] = buf[l++];
  while (r < rhi) idx[i++] = buf[r++];
}

/* ------------------------------------------------------------------ */

extern "C" {
  /* Typedef'd pointer to get abstract datatype. */
  typedef struct regex_t *re_t;

  /* Compile regex string pattern to a regex_t-array. */
  static re_t re_compile(const char *pattern);


  /* Find matches of the compiled pattern inside text. */
  static int  re_matchp(const char *text, re_t pattern);


/* Definitions: */

#define MAX_REGEXP_OBJECTS 30 /* Max number of regex symbols in expression. */
#define MAX_CHAR_CLASS_LEN 40 /* Max length of character-class buffer in.   */


  enum { UNUSED, DOT, BEGIN, END, QUESTIONMARK, STAR, PLUS,
         CHAR, CHAR_CLASS, INV_CHAR_CLASS, DIGIT, NOT_DIGIT,
         INTEGER, NOT_INTEGER, FLOAT, NOT_FLOAT,
         ALPHA, NOT_ALPHA, WHITESPACE, NOT_WHITESPACE /*, BRANCH */ };

  typedef struct regex_t {
    unsigned char  type;   /* CHAR, STAR, etc.                      */
    union {
      unsigned char  ch;   /*      the character itself             */
      unsigned char *ccl;  /*  OR  a pointer to characters in class */
    };
  } regex_t;

/* Private function declarations: */
  static int matchpattern(regex_t *pattern, const char *text);
  static int matchcharclass(char c, const char *str);
  static int matchstar(regex_t p, regex_t *pattern, const char *text);
  static int matchplus(regex_t p, regex_t *pattern, const char *text);
  static int matchone(regex_t p, char c);
  static int matchdigit(char c);
  static int matchint(char c);
  static int matchfloat(char c);
  static int matchalpha(char c);
  static int matchwhitespace(char c);
  static int matchmetachar(char c, const char *str);
  static int matchrange(char c, const char *str);
  static int ismetachar(char c);

/* Semi-public functions: */
  int re_match(const char *text, const char *pattern)
  {
    return re_matchp(text, re_compile(pattern));
  }

  int re_matchp(const char *text, re_t pattern)
  {
    if (pattern != 0) {
      if (pattern[0].type == BEGIN) {
        return ((matchpattern(&pattern[1], text)) ? 0 : -1);
      } else {
        int idx = -1;

        do {
          idx += 1;

          if (matchpattern(pattern, text)) {
            if (text[0] == '\0')
              return -1;

            return idx;
          }
        }
        while (*text++ != '\0');
      }
    }
    return -1;
  }

  re_t re_compile(const char *pattern)
  {
    /* The sizes of the two static arrays below substantiates the static RAM usage of this module.
       MAX_REGEXP_OBJECTS is the max number of symbols in the expression.
       MAX_CHAR_CLASS_LEN determines the size of buffer for chars in all char-classes in the expression. */
    static regex_t re_compiled[MAX_REGEXP_OBJECTS];
    static unsigned char ccl_buf[MAX_CHAR_CLASS_LEN];
    int ccl_bufidx = 1;

    char c;     /* current char in pattern   */
    int i = 0;  /* index into pattern        */
    int j = 0;  /* index into re_compiled    */

    while (pattern[i] != '\0' && (j+1 < MAX_REGEXP_OBJECTS)) {
      c = pattern[i];

      switch (c) {
        /* Meta-characters: */
      case '^': {    re_compiled[j].type = BEGIN;           } break;
      case '$': {    re_compiled[j].type = END;             } break;
      case '.': {    re_compiled[j].type = DOT;             } break;
      case '*': {    re_compiled[j].type = STAR;            } break;
      case '+': {    re_compiled[j].type = PLUS;            } break;
      case '?': {    re_compiled[j].type = QUESTIONMARK;    } break;

        /* Escaped character-classes (\s \w ...): */
      case '\\': {
        if (pattern[i+1] != '\0') {
          /* Skip the escape-char '\\' */
          i += 1;
          /* ... and check the next */
          switch (pattern[i]) {
            /* Meta-character: */
          case 'd': {    re_compiled[j].type = DIGIT;            } break;
          case 'D': {    re_compiled[j].type = NOT_DIGIT;        } break;
          case 'i': {    re_compiled[j].type = INTEGER;          } break;
          case 'I': {    re_compiled[j].type = NOT_INTEGER;      } break;
          case 'f': {    re_compiled[j].type = FLOAT;            } break;
          case 'F': {    re_compiled[j].type = NOT_FLOAT;        } break;
          case 'w': {    re_compiled[j].type = ALPHA;            } break;
          case 'W': {    re_compiled[j].type = NOT_ALPHA;        } break;
          case 's': {    re_compiled[j].type = WHITESPACE;       } break;
          case 'S': {    re_compiled[j].type = NOT_WHITESPACE;   } break;

            /* Escaped character, e.g. '.' or '$' */
          default: {
            re_compiled[j].type = CHAR;
            re_compiled[j].ch = pattern[i];
          } break;
          }
        }
        /* '\\' as last char in pattern -> invalid regular expression. */
      } break;

        /* Character class: */
      case '[': {
        /* Remember where the char-buffer starts. */
        int buf_begin = ccl_bufidx;

        /* Look-ahead to determine if negated */
        if (pattern[i+1] == '^') {
          re_compiled[j].type = INV_CHAR_CLASS;
          i += 1; /* Increment i to avoid including '^' in the char-buffer */
        } else {
          re_compiled[j].type = CHAR_CLASS;
        }

        /* Copy characters inside [..] to buffer */
        while ((pattern[++i] != ']') && (pattern[i] != '\0')) {
          /* Missing ] */
          if (pattern[i] == '\\') {
            if (ccl_bufidx >= MAX_CHAR_CLASS_LEN - 1) {
              return 0;
            }
            ccl_buf[ccl_bufidx++] = pattern[i++];
          } else if (ccl_bufidx >= MAX_CHAR_CLASS_LEN) {
            return 0;
          }
          ccl_buf[ccl_bufidx++] = pattern[i];
        }
        if (ccl_bufidx >= MAX_CHAR_CLASS_LEN) {
          /* Catches cases such as [00000000000000000000000000000000000000][ */
          return 0;
        }
        /* Null-terminate string end */
        ccl_buf[ccl_bufidx++] = 0;
        re_compiled[j].ccl = &ccl_buf[buf_begin];
      } break;

        /* Other characters: */
      default: {
        re_compiled[j].type = CHAR;
        re_compiled[j].ch = c;
      } break;
      }
      i += 1;
      j += 1;
    }
    /* 'UNUSED' is a sentinel used to indicate end-of-pattern */
    re_compiled[j].type = UNUSED;

    return (re_t) re_compiled;
  }


/* Private functions: */
  static int matchdigit(char c)
  {
    return ((c >= '0') && (c <= '9'));
  }

  static int matchint(char c)
  {
    return (matchdigit(c) || (c == '-') || (c == '+'));
  }

  static int matchfloat(char c)
  {
    return (matchint(c) || (c == '.') || (c == 'e') || (c == 'E'));
  }

  static int matchalpha(char c)
  {
    return ((c >= 'a') && (c <= 'z')) || ((c >= 'A') && (c <= 'Z'));
  }

  static int matchwhitespace(char c)
  {
    return ((c == ' ') || (c == '\t') || (c == '\n') || (c == '\r') || (c == '\f') || (c == '\v'));
  }

  static int matchalphanum(char c)
  {
    return ((c == '_') || matchalpha(c) || matchdigit(c));
  }

  static int matchrange(char c, const char *str)
  {
    return ((c != '-') && (str[0] != '\0')
            && (str[0] != '-') && (str[1] == '-')
            && (str[1] != '\0') && (str[2] != '\0')
            && ((c >= str[0]) && (c <= str[2])));
  }

  static int ismetachar(char c)
  {
    return ((c == 's') || (c == 'S')
            || (c == 'w') || (c == 'W')
            || (c == 'd') || (c == 'D'));
  }

  static int matchmetachar(char c, const char *str)
  {
    switch (str[0]) {
    case 'd': return  matchdigit(c);
    case 'D': return !matchdigit(c);
    case 'i': return  matchint(c);
    case 'I': return !matchint(c);
    case 'f': return  matchfloat(c);
    case 'F': return !matchfloat(c);
    case 'w': return  matchalphanum(c);
    case 'W': return !matchalphanum(c);
    case 's': return  matchwhitespace(c);
    case 'S': return !matchwhitespace(c);
    default:  return (c == str[0]);
    }
  }

  static int matchcharclass(char c, const char *str)
  {
    do {
      if (matchrange(c, str)) {
        return 1;
      } else if (str[0] == '\\') {
        /* Escape-char: increment str-ptr and match on next char */
        str += 1;
        if (matchmetachar(c, str)) {
          return 1;
        } else if ((c == str[0]) && !ismetachar(c)) {
          return 1;
        }
      } else if (c == str[0]) {
        if (c == '-') {
          return ((str[-1] == '\0') || (str[1] == '\0'));
        } else {
          return 1;
        }
      }
    }
    while (*str++ != '\0');

    return 0;
  }

  static int matchone(regex_t p, char c)
  {
    switch (p.type) {
    case DOT:            return 1;
    case CHAR_CLASS:     return  matchcharclass(c, (const char *)p.ccl);
    case INV_CHAR_CLASS: return !matchcharclass(c, (const char *)p.ccl);
    case DIGIT:          return  matchdigit(c);
    case NOT_DIGIT:      return !matchdigit(c);
    case INTEGER:        return  matchint(c);
    case NOT_INTEGER:    return !matchint(c);
    case FLOAT:          return  matchfloat(c);
    case NOT_FLOAT:      return !matchfloat(c);
    case ALPHA:          return  matchalphanum(c);
    case NOT_ALPHA:      return !matchalphanum(c);
    case WHITESPACE:     return  matchwhitespace(c);
    case NOT_WHITESPACE: return !matchwhitespace(c);
    default:             return  (p.ch == c);
    }
  }

  static int matchstar(regex_t p, regex_t *pattern, const char *text)
  {
    do {
      if (matchpattern(pattern, text))
        return 1;
    }
    while ((text[0] != '\0') && matchone(p, *text++));

    return 0;
  }

  static int matchplus(regex_t p, regex_t *pattern, const char *text)
  {
    while ((text[0] != '\0') && matchone(p, *text++)) {
      if (matchpattern(pattern, text))
        return 1;
    }
    return 0;
  }

  static int matchquestion(regex_t p, regex_t *pattern, const char *text)
  {
    if (p.type == UNUSED)
      return 1;
    if (matchpattern(pattern, text))
      return 1;
    if (*text && matchone(p, *text++))
      return matchpattern(pattern, text);
    return 0;
  }

/* Iterative matching */
  static int matchpattern(regex_t *pattern, const char *text)
  {
    do {
      if ((pattern[0].type == UNUSED) || (pattern[1].type == QUESTIONMARK)) {
        return matchquestion(pattern[0], &pattern[2], text);
      } else if (pattern[1].type == STAR) {
        return matchstar(pattern[0], &pattern[2], text);
      } else if (pattern[1].type == PLUS) {
        return matchplus(pattern[0], &pattern[2], text);
      } else if ((pattern[0].type == END) && pattern[1].type == UNUSED) {
        return (text[0] == '\0');
      }
    }
    while ((text[0] != '\0') && matchone(*pattern++, *text++));

    return 0;
  }

}
