/*----------------------------------------------------------------------------*
 * Discrete Algorithms Module - Access from Python to basic algorithms on     *
 * sequences and graphs                                                       *
 *                                                                            *
 * Copyright (C) 2004 Nils Weskamp                                            *
 * Portions Copyright (C) 1999-2004 Cambridge Crystallographic Data Centre    *
 * Portions Copyright (C) 2000-2004 Stefan Schmitt,Daniel Kuhn,Gerhard Klebe  *
 *                                                                            *
 * Permission to use, copy, modify, and distribute this software and its      * 
 * documentation for any purpose, without fee, and without a written          *
 * agreement is hereby granted, provided that the above copyright notice and  *
 * this paragraph and the following two paragraphs appear in all copies.      *
 *                                                                            *
 * IN NO EVENT SHALL THE AUTHOR BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,  *
 * SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,     *
 * ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE * 
 * AUTHOR HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                 *
 *                                                                            *
 * THE AUTHOR SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT       *
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A    *
 * PARTICULAR PURPOSE. THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS,*
 * AND THE AUTHOR HAS NO OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT,         *
 * UPDATES, ENHANCEMENTS, OR MODIFICATIONS.                                   *
 *                                                                            *
 * Contact information:                                                       *
 *                                                                            *
 *    Nils Weskamp                                                            *
 *    Philipps-University, Marburg                                            *
 *    Department of Mathematics and Computer Science                          *
 *    Intelligent Systems Group                                               *
 *    Hans-Meerwein Straﬂe                                                    *
 *    35032 Marburg, Germany                                                  *
 *                                                                            *
 *    weskamp@informatik.uni-marburg.de                                       *
 *----------------------------------------------------------------------------*/

#include <iostream>
#include <string>
/*----------------------------------------------------------------------------*/
#include "Python.h"
#include "clique/bk.h"
#include "clique/util.h"
#include "arrayobject.h"
//#include "smiles.h"
/*----------------------------------------------------------------------------*/
using namespace std;
/*----------------------------------------------------------------------------*/
#define MODULE_NAME               Algorithms     
#define MODULE_INIT               initAlgorithms 
#define MODULE_STR                "Algorithms"   
#define LASTCHANGE                "01.04.2004"
#define MAX_CLIQ_NUM              10000


/*----------------------------------------------------------------------------*/
// This class is used as a context object to store the cliques detected by the
// clique detection algorithms and is called from the respective callback-
// function. It can store up to MAX_CLIQ_NUM cliques. 
class CliqueCollector {
/*----------------------------------------------------------------------------*/
public:

  Set cliques[MAX_CLIQ_NUM];
  long cliques_counter;

/*----------------------------------------------------------------------------*/
  CliqueCollector () {
/*----------------------------------------------------------------------------*/
    cliques_counter = 0;
  }

};

/*----------------------------------------------------------------------------*/
static char find_max_clique__doc__[] = 
"This function takes a two-dimensional Numeric-array and \
interpretes it as the adjacency matrix of a graph. It then \
uses the algorithm of Bron and Kerbosch to detect maximal \
cliques in that graph and returns the nodes of a maximal clique. \
It is assumed that the nodes 0..n-1 are given in the order as they \
appear in the adjacency matrix.";
/*----------------------------------------------------------------------------*/
static PyObject *find_max_clique(PyObject *self, PyObject *args) {
/*----------------------------------------------------------------------------*/
  PyArrayObject *array;

  // Parameter Checking

  if (!PyArg_ParseTuple(args, "O!", &PyArray_Type, &array))
    return NULL;

  if (array->nd != 2 || array->descr->type_num != PyArray_LONG) {
    PyErr_SetString(PyExc_ValueError, 
		     "array must be two-dimensional and of type int");
    return NULL;
  }

  if (array->dimensions[0] != array->dimensions[1]) {
    PyErr_SetString(PyExc_ValueError, 
		     "array must be quadratic");
    return NULL;
  }

  // Copying the adjacency matrix
  
  int n = array->dimensions[0];
  char *connected[n];
  for(int i = 0; i < n; i++) {
    connected[i] = (char *) malloc(n * sizeof(char));
  }
  
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < n; j++) {
   
      connected[i][j] = *(long*)(array->data + i*array->strides[0] + 
				 j*array->strides[1]);

      if(i == j) {
	connected[i][j] = 1;
      }
    }
  }

  // Calling the Bron-Kerbosch-Algorithm
  
  Set cliq;
  init_Set(&cliq, n);

  bron_kerbosch(n, connected, NULL, &cliq, 
		(clique_callback)print_clique_callback, 
		NULL);

  // Building result list

  PyObject* erg = PyList_New(0);

  for (int i=0; i<cliq.size; i++) {
    cerr << cliq.vertex[i] << " ";
    PyList_Append(erg, PyInt_FromLong(cliq.vertex[i]));
  }

  return erg;
  
}

/*----------------------------------------------------------------------------*/
// This is the callback-function used to collect all the cliques found by the
// algorithm of bron-kerbosch in find_cliques(). It stores the cliques in the
// CliqueCollector-object that is passed as context.
int collect_cliques(const _set *st, void *context)
/*----------------------------------------------------------------------------*/
{

  CliqueCollector * collector = (CliqueCollector*) context;

  if(collector->cliques_counter < MAX_CLIQ_NUM) {
    init_Set(&(collector->cliques[collector->cliques_counter]), st->size);
    copy_Set(st, &(collector->cliques[collector->cliques_counter]));
    collector->cliques_counter++;
    return CLIQUE_CONTINUE;
  }

  return CLIQUE_FOUND;
}

/*----------------------------------------------------------------------------*/
static char find_cliques__doc__[] = 
"This function takes a two-dimensional Numeric-array and \
interpretes it as the adjacency matrix of a graph. It then \
uses the algorithm of Bron and Kerbosch to detect maximal \
cliques in that graph and returns a list of lists containing \
the nodes of all maximal cliques (up to 10000). It is assumed \
that the nodes 0..n-1 are given in the order as they \
appear in the adjacency matrix.";
/*----------------------------------------------------------------------------*/
static PyObject *find_cliques(PyObject *self, PyObject *args) {
/*----------------------------------------------------------------------------*/
  PyArrayObject *array;

  // Parameter Checking

  if (!PyArg_ParseTuple(args, "O!", &PyArray_Type, &array))
    return NULL;

  if (array->nd != 2 || array->descr->type_num != PyArray_LONG) {
    PyErr_SetString(PyExc_ValueError, 
		     "array must be two-dimensional and of type int");
    return NULL;
  }

  if (array->dimensions[0] != array->dimensions[1]) {
    PyErr_SetString(PyExc_ValueError, 
		     "array must be quadratic");
    return NULL;
  }

  // Copying the adjacency matrix
  
  int n = array->dimensions[0];
  char *connected[n];
  for(int i = 0; i < n; i++) {
    connected[i] = (char *) malloc(n * sizeof(char));
  }
  
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < n; j++) {
   
      connected[i][j] = *(long*)(array->data + i*array->strides[0] + 
				 j*array->strides[1]);

      if(i == j) {
	connected[i][j] = 1;
      }

    }
  }

  // Calling the Bron-Kerbosch-Algorithm
  
  Set cliq;
  init_Set(&cliq, n);

  CliqueCollector collector;

  bron_kerbosch(n, connected, NULL, &cliq, 
		(clique_callback)collect_cliques, 
		&collector);

  // Building result list

  PyObject* result = PyList_New(0);

  for(int i = 0; i < collector.cliques_counter; i++) {

    PyObject* cli = PyList_New(0);

    for (int j=0; j < collector.cliques[i].size; j++) {
      PyList_Append(cli, PyInt_FromLong(collector.cliques[i].vertex[j]));
    }

    PyList_Append(result, cli);
  }
  
  return result;
  
}

/*----------------------------------------------------------------------------*/
//static char smiles2graph__doc__[] = 
//"";
/*----------------------------------------------------------------------------*/
//static PyObject *smiles2graph(PyObject *self, PyObject *args) {
/*----------------------------------------------------------------------------*/
//
//  char* s; 
//  int n = -1;
//  int ok;
//  
//  ok = PyArg_ParseTuple(args, "s", &s);
//
//  THB_NCONNECT * t = sms_make_thb(s, &n);
//
//  cout << "n: " << n << endl;
//
//  for(int i = 0; i < n; i++) {
//    printf(" %d: ", t[i].id[0]);
//    for(int j = 0; j < t[i].n_con; j++) {
//      printf("%d, ", t[i].connect[j]);
//    }
//    printf("\n");
//  }
//
//  Py_INCREF(Py_None);
//  return Py_None;
//
//}

/*----------------------------------------------------------------------------*/
static PyMethodDef ModuleMethods[] =   {
/*----------------------------------------------------------------------------*/
  {"find_max_clique",find_max_clique,1,find_max_clique__doc__}, 
  {"find_cliques",   find_cliques,   1,find_cliques__doc__}, 
  //  {"smiles2graph",   smiles2graph,   1,smiles2graph__doc__}, 
  {NULL,NULL} // sentinel
};

/*----------------------------------------------------------------------------*/
extern "C" {
  void MODULE_INIT()
  {
    
    PyObject *ThisModule = Py_InitModule( MODULE_STR , ModuleMethods);
    PyObject *ModuleDict = PyModule_GetDict(ThisModule);  

    // for correct handling of Numpy
    import_array();           

    // display welcome message

/*----------------------------------------------------------------------------/
 *   cerr << "Discrete Algorithms Module " << VERSION
 *	 << " [Build "   << TIMESTAMP << "]" << endl 
 *	 << "Send complaints to: Nils Weskamp < weskamp@informatik.uni-marburg.de >" << endl;               *----------------------------------------------------------------------------*/
 
  };
}

/*----------------------------------------------------------------------------*/
