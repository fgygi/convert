////////////////////////////////////////////////////////////////////////////////
//
//  convert: conversion of units
//
//  F.Gygi, Jun 1994, Feb 1995, revised Feb 2016
//  units are stored in a weighted tree
//  conversion is done using the depth first search algorithm
//
//  The definition of units and their relations must be given in a file
//  named convert.def which contains definitions like:
//
//  node Ha Hartree  // define a node
//  edge meV 0.001 eV FALSE  // define an edge
//
//  The current directory is first searched for a convert.def file
//  if none is found, the file HOME/bin/convert.def is searched
//
//  use: convert 25 meV K
//  converts from meV to Kelvin
//
//  constants in .def file from:
//  Physics Vade Mecum, ed. by H.L.Anderson, AIP (1981)
//
//  compilation: g++ -o cv convert.C
//
////////////////////////////////////////////////////////////////////////////////

#include<iostream>
#include<iomanip>
#include<cstdlib>
#include<cstdio>
#include<cstring>
#include<sys/stat.h>
using namespace std;

#define TRUE 1
#define FALSE 0

struct node { char *name; char *long_name; struct edge *adj_list;
              int visited; struct node *next; };
struct edge { struct node *to_node; double factor;
              int inverse; struct edge *next; };

FILE *defFile;
char *homedir,defFileName[64];
double convfac;
int    invflag;
char line[256],type[32],shortname[32],longname[32],
     from_name[32],to_name[32],invstr[32];

double value, result;
int    found;

struct node *unit_list = NULL;

void add_node( char *new_name, char *new_long_name );
void add_edge( char *name1, double fac12, char *name2, int inversion );
double convert( double value, char *from_unit, char *to_unit );
void connect ( struct node *n1, struct node *n2, double val );
struct node *find_node ( char *name, struct node *list );

int main( int argc, char **argv )
{
  struct node *t;
  char *from_unit, *to_unit;

  // locate definition file:
  // Look first in current directory
  strcpy(defFileName,"convert.def");
  struct stat statbuf;
  // stat() returns 0 if file is present
  int status = stat(defFileName,&statbuf);

  if ( status )
  {
    // File is not in current directory. Try in HOME/bin
#ifdef DEBUG
    cerr << " File " << defFileName
         << " not found " << endl;
    cerr << " Trying HOME/bin/convert.def " << endl;
#endif
    homedir = getenv("HOME");
    strcpy(defFileName,homedir);
#ifdef DEBUG
    cerr << " Home directory is " << defFileName << endl;
#endif
    strcat(defFileName,"/bin/convert.def");
#ifdef DEBUG
    cerr << " Target definition file is " << defFileName << endl;
#endif
  }

  // Read definitions from file convert.def
  defFile= fopen ( defFileName, "r" ) ;
  if ( !defFile )
  {
    cerr << " Cannot open definition file" << endl;
    exit(1);
  }

  while ( !feof(defFile) )
  {
    fgets( line, 256, defFile );
    if ( !feof(defFile) )
    {
      // check if comment line
      if ( line[0] == '#' )
      {
        // comment line: do nothing
#ifdef DEBUG
        cerr << " Comment: " << line << endl;
#endif
      }
      else
      {
        sscanf(line,"%s",type);
        // define node or edge
        if ( !strcmp(type,"node") )
        {
          sscanf(line,"%s %s %s",type,shortname,longname);
#ifdef DEBUG
          cerr << " defining node "
               << shortname << " "
               << longname << endl;
#endif
          add_node(shortname,longname);
        }
        else if ( !strcmp(type,"edge") )
        {
          sscanf(line,"%s %s %lf %s %s",type,from_name,&convfac,to_name,invstr);
#ifdef DEBUG
          cerr << " defining conversion from "
               << from_name << " to " << to_name
               << " factor: " << convfac
               << " inv: " << invstr << endl;
#endif
          if ( strcmp(invstr,"INVERT") && strcmp(invstr,"NOINVERT") )
          {
            cerr << " Error in definition file: inversion flag "
                 << "must be INVERT or NOINVERT" << endl;
            exit(1);
          }
          invflag = !strcmp(invstr,"INVERT");
          add_edge ( from_name, convfac, to_name, invflag );
        }
        else
        {
          cerr << " invalid type in definition file: " << type << endl;
          exit(1);
        }
      }
    }
  }

  if ( argc < 4 )
  {
    cerr << " cv: unit conversions: " << endl;
        cerr << " Current definition file is " << defFileName << endl;
        cerr << " use: cv value from_unit to_unit " << endl;
        cerr << " allowed units are: " << endl;
        t = unit_list;
        while ( t )
        {
          cerr << " "
          << setw(12)
          << setiosflags(ios::left)
          << t->name
          << setiosflags(ios::left)
          << t->long_name
          << endl;

          t = t->next;
        }
        exit ( EXIT_SUCCESS );
  }

  value = atof( argv[1] );
  from_unit = argv[2];
  to_unit = argv[3];

  result = convert( value, from_unit, to_unit );

  cout << " "
       << setprecision(8)
       << value << " "
       << from_unit << " = "
       << setprecision(8)
       << result << " " << to_unit << endl;

  return ( EXIT_SUCCESS );
}

void add_node( char *new_name, char *new_long_name )
{
  /* add unit named "new_name" to the unit list */
  struct node *t;
  t = find_node( new_name, unit_list );
  if ( t )
  {
    cerr << " warning: unit " << new_name << " is already defined" << endl;
  }
  else
  {
    t = ( struct node * ) malloc ( sizeof( *t ) );
    t->next = unit_list;
    t->adj_list = NULL;
    t->visited = FALSE;
    t->name = ( char * ) malloc ( (strlen(new_name)+1) * sizeof( char ) );
    strcpy ( t->name, new_name );
    t->long_name = ( char * )
      malloc ( (strlen(new_long_name)+1) * sizeof( char ) );
    strcpy ( t->long_name, new_long_name );
    unit_list = t;
  }
}

void add_edge( char *name1, double fac12, char *name2, int inversion )
{
  struct node *n1, *n2;
  struct edge *t;

  if ( fac12 == 0.0 )
  {
    cerr << " Conversion factor from " << name1 << " to "
             << name2 << " is zero" << endl;
    exit ( EXIT_FAILURE );
  }

  n1 = find_node( name1, unit_list );
  if ( !n1 )
  {
    cerr << " add_edge: unit " << name1 << " not found " << endl;
    exit ( EXIT_FAILURE );
  }
  n2 = find_node( name2, unit_list );
  if ( !n2 )
  {
    cerr << " add_edge: unit " << name2 << " not found " << endl;
    exit ( EXIT_FAILURE );
  }

  /* add edge to the adjacency lists of n1 and n2 */
  t = ( struct edge * ) malloc ( sizeof( *t ) );
  t->to_node = n2;
  t->factor = fac12;
  t->inverse = inversion;
  t->next = n1->adj_list;
  n1->adj_list = t;

  t = ( struct edge * ) malloc ( sizeof( *t ) );
  t->to_node = n1;
  if ( !inversion )
  {
    t->factor = 1.0 / fac12;
  }
  else
  {
    t->factor = fac12;
  }
  t->inverse = inversion;
  t->next = n2->adj_list;
  n2->adj_list = t;
}

double convert( double value, char *from_unit, char *to_unit )
{
  struct node *fu, *tu;
  fu = find_node( from_unit, unit_list );
  if ( !fu )
  {
    cerr << " convert: unit " << from_unit << " not found " << endl;
    exit ( EXIT_FAILURE );
  }
  tu = find_node( to_unit, unit_list );
  if ( !tu )
  {
    cerr << " convert: unit " << to_unit << " not found " << endl;
    exit ( EXIT_FAILURE );
  }

  /* connect from_unit to to_unit */

  found = FALSE;
  connect ( fu, tu, value );
  if ( !found )
  {
    cerr << " Cannot convert " << from_unit << " to "
             << to_unit << endl;
    exit ( EXIT_FAILURE );
  }

  return result;
}

void connect ( struct node *n1, struct node *n2, double val )
{
  struct edge *t;

  /* Check if destination is reached */
  if ( n1 == n2 )
  {
    result = val;
    found = TRUE;
  }

  n1->visited = TRUE;

  t = n1->adj_list;
  while ( t )
  {
    if ( !t->to_node->visited )
    {
      /* attempt connection from t->to_node */
      if ( t->inverse )
      {
        if ( val == 0 )
        {
          cerr << " Cannot convert zero value " << endl;
          exit ( EXIT_FAILURE );
        }
        connect ( t->to_node, n2, (t->factor) / val );
      }
      else
      {
        connect ( t->to_node, n2, (t->factor) * val );
      }
    }
    t = t->next;
  }
}

struct node *find_node ( char *name, struct node *list )
{
  /* find node named name in list */
  struct node *t = list;
  while ( t && strcmp( name, t->name ) )
    t = t->next;
  return t;
}
