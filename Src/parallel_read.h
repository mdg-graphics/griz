#define PARALLEL_READER_SUPPORT
#ifdef PARALLEL_READER_SUPPORT

#ifndef PARALLEL_READ_H
#define PARALLEL_READ_H

#include "misc.h"
#include "viewer.h"

#define ZERO_PROCESSOR 0

/* Python program name */
const char* PROGRAM_NAME = "Griz-mili-python-reader";

const char* MILI_READER_MODULE_NAME = "mili.reader";

/* Mili Reader function names as constants */
const char* OPEN_DATABASE = "open_database";
const char* PARAMETER_GETTER = "get_parameter"; // NOTE: May not need. Definitely need to clean up
const char* QUERY_FUNCTION = "query";

/* Griz Struct function calls */
const char* GET_GEOM_CALL = "griz_get_geom_call";
const char* GET_ST_DESCRIPTORS_CALL = "griz_get_st_descriptors_call";
const char* GET_STATE_CALL = "griz_get_state_call";
const char* STATE_QTY = "state_qty";

const int QTY_CONNECTS[] = { 0, 0, 2, 3, 4, 4, 4, 5, 6, 8, 0, 0, 0, 1, 10 };

const int PROC_COUNT = 8; // FOR NOW...

/* Global reference to mili reader object */
PyObject *mili_db = NULL;
PyObject *mili_reader_module = NULL;

int bin_search_index(int target, int* array, int size){
    int low = 0;
    int high = size - 1;
    int mid;

    while(low <= high){
        mid = low + (high - low) / 2;
        if( target == array[mid])
            return mid;
        else if( target < array[mid])
            high = mid-1;
        else
            low = mid+1;
    }

    return -1;
}


/* TAG( call_mili_module_function )
 *
 * Given a specified function name and a set of arguments, check that the function
 * exists in the mili module and call it. Returns the PyObject pointer returned from
 * calling the python function.
 *
 */
PyObject*
call_mili_module_function(const char *function_name, PyObject *py_Arglist){
    PyObject *py_Dict,
             *py_Func;
    PyObject *py_ReturnValue = Py_None;

    /* Check that mili_reader_module is not NULL */
    if(mili_reader_module != NULL){
        /* Get __dict__ for mili reader module */
        py_Dict = PyModule_GetDict(mili_reader_module);
        if(py_Dict != NULL){
            /* Check for requested function and call it */
            py_Func = PyDict_GetItemString(py_Dict, function_name);
            if(py_Func != NULL && PyCallable_Check(py_Func)){
                py_ReturnValue = PyObject_CallObject(py_Func, py_Arglist);
                PyErr_Print();
            }
            else{
                PyErr_Print();
            }
        }
        else{
            PyErr_Print();
        }
    }
    else{
        fprintf(stderr, "mili_reader_module is NULL\n");
    }

    return py_ReturnValue;
}

/* TAG( string_to_pyobject )
 *
 * Converts a string to a PyObject representation of string.
 */
PyObject*
string_to_pyobject(const char *str){
    int length = strlen(str);
    PyObject * py_Str = PyUnicode_Decode(str, length, "ascii", "ignore");
    return py_Str;
}


/* TAG( call_mili_reader_function_noargs )
 *
 * TODO: 
 */
PyObject*
call_mili_reader_function_noargs(const char * function_name)
{
    PyObject * py_FuncName = string_to_pyobject( function_name );
    PyObject * py_Result = PyObject_CallMethodObjArgs( mili_db, py_FuncName, NULL );
    PyErr_Print();
    return py_Result;
}



/* TAG( pyobject_as_string )
 *
 * Converts a PyObject* to a char*
 */
char *
pyobject_as_string( PyObject *py_Str )
{
    char * str = NULL;
    PyObject * py_UniStr = PyUnicode_AsEncodedString(py_Str, "ascii", "ignore");
    if ( py_UniStr != NULL )
        str = PyBytes_AsString( py_UniStr );
    return str;
}


/* TAG( get_pyobject_attribute )
 *
 * Retrieve a specific attribute from a PyObject *.
 */
PyObject*
get_pyobject_attribute( PyObject* object, char* attr )
{
    PyObject * py_Attribute = Py_None;
    if( PyObject_HasAttrString(object, attr) )
        py_Attribute = PyObject_GetAttrString( object, attr );
    return py_Attribute;
}


/* TAG( get_pyobject_attribute_as_string )
 *
 * Retrieve a specific attribute from a PyObject * as a string.
 * Does no error checking for performance.
 */
char*
get_pyobject_attribute_as_string( PyObject* object, const char* attr )
{
    PyObject * py_Attribute = PyObject_GetAttrString( object, attr );
    PyObject * py_AttributeStr = PyUnicode_AsEncodedString(py_Attribute, "ascii", "ignore");
    char * str = PyBytes_AsString(py_AttributeStr);
    return str;
}

/* TAG( get_pyobject_attribute_as_int )
 *
 * Retrieve a specific attribute from a PyObject * as an int.
 * Does no error checking for performance.
 */
int
get_pyobject_attribute_as_int( PyObject* object, const char* attr )
{
    PyObject * py_Attribute = PyObject_GetAttrString( object, attr );
    PyErr_Print();
    int value = PyLong_AsLong(py_Attribute);
    return value;
}

/* TAG( get_pyobject_attribute_as_double )
 *
 * Retrieve a specific attribute from a PyObject * as a double.
 * Does no error checking for performance.
 */
double
get_pyobject_attribute_as_double( PyObject* object, const char* attr )
{
    PyObject * py_Attribute = PyObject_GetAttrString( object, attr );
    double value = PyFloat_AsDouble(py_Attribute);
    return value;
}


/* TAG( integer_pointer_from_pyobject )
 *
 * Convert pyobject list to int array
 */
void
integer_pointer_from_pyobject( PyObject * py_List, int size, int * values )
{
    int i;
    PyObject * py_LongValue;
    for( i = 0; i < size; i++ ){
        py_LongValue = PySequence_GetItem(py_List, i);
        *(values+i) = PyLong_AsLong(py_LongValue);
    }
}


/* TAG( float_pointer_from_pyobject )
 *
 * Convert pyobject list to float array
 */
void
float_pointer_from_pyobject( PyObject * py_List, int size, float * values )
{
    int i;
    PyObject * py_FloatValue;
    for( i = 0; i < size; i++ ){
        py_FloatValue = PySequence_GetItem(py_List, i);
        *(values+i) = PyFloat_AsDouble(py_FloatValue);
    }
}



/* TAG( mili_reader_get_string )
 *
 * Query a string parameter from the mili reader.
 */
Bool_type
mili_reader_get_string(char *parameter_name, char *str){
    PyObject * py_FuncName,
             * py_Arglist;
    PyObject * py_ParameterName,
             * py_ParameterStr,
             * py_Parameter;

    /* Build parameter name and processor arguments */
    py_ParameterName = Py_BuildValue("s", parameter_name);

    /* Build argument list */
    py_Arglist = PyTuple_New(2);
    PyTuple_SetItem( py_Arglist, 0, mili_db );
    PyTuple_SetItem( py_Arglist, 1, py_ParameterName );

    /* Lookup parameter in mili and convert to string */
    py_Parameter = call_mili_module_function( PARAMETER_GETTER, py_Arglist );

    /* Get string representation of parameter */
    if(py_Parameter != NULL){
        str = pyobject_as_string( py_Parameter );
        if ( str == NULL )
            return FALSE;
    }

    return FALSE;
}


/* TAG( mili_reader_get_double )
 *
 * Query a double parameter from the mili reader.
 */
Bool_type
mili_reader_get_double(char *parameter_name, double *result){
    PyObject * py_FuncName,
             * py_Arglist;
    PyObject * py_MiliParameter,
             * py_ParameterName;

    /* Build parameter name as a python string */
    py_ParameterName = Py_BuildValue("s", parameter_name);

    /* Build argument list */
    py_Arglist = PyTuple_New(2);
    PyTuple_SetItem( py_Arglist, 0, mili_db );
    PyTuple_SetItem( py_Arglist, 1, py_ParameterName );

    /* Lookup parameter in mili and convert to string */
    py_MiliParameter = call_mili_module_function( PARAMETER_GETTER, py_Arglist );

    /* Get double representation of parameter */
    if(py_MiliParameter != NULL){
        *result = PyFloat_AsDouble(py_MiliParameter);
        return OK;
    }

    return NOT_OK;
}

/* TAG( mili_reader_get_int )
 *
 * Query an integer parameter from the mili reader.
 */
Bool_type
mili_reader_get_int(char *parameter_name, int *result){
    PyObject * py_FuncName,
             * py_Arglist;
    PyObject * py_MiliParameter,
             * py_ParameterName;

    /* Build parameter name as a python string */
    py_ParameterName = Py_BuildValue("s", parameter_name);

    /* Build argument list */
    py_Arglist = PyTuple_New(2);
    PyTuple_SetItem( py_Arglist, 0, mili_db );
    PyTuple_SetItem( py_Arglist, 1, py_ParameterName );

    /* Lookup parameter in mili and convert to string */
    py_MiliParameter = call_mili_module_function( PARAMETER_GETTER, py_Arglist );

    /* Get double representation of parameter */
    if(py_MiliParameter != NULL){
        *result = PyLong_AsLong(py_MiliParameter);
        return OK;
    }

    return NOT_OK;
}


/* TAG( mili_reader_db_open )
 *
 * Open the specified plot files with the mili reader and set mili_db to reference the newly created
 * python reader object.
 */
int
mili_reader_db_open( char *path_root, int *p_dbid );


/* TAG( mili_reader_get_results )
 *
 * Query the Mili reader for some results
 */
Bool_type
mili_reader_get_results( int dbid, int state, int subrec_id, int qty,
                         char **results, void *data );


/* TAG( mili_reader_get_geom )
 *
 * Load in the geometry (mesh definition) data from the mili reader
 */
Bool_type
mili_reader_get_geom( Analysis * analy, int dbid, Mesh_data **p_mtable, int *p_mesh_qty );


/* TAG( compute_label_blocking_data )
 *
 * For a given set of labels, compute the blocking data and return block_qty.
 */
int *
compute_label_blocking_data(int * labels, int qty, int * num_blocks);


/* TAG( mili_reader_get_st_descriptors )
 *
 * Query mili reader and store information about the available
 * state record formats and referenced state variables.
 */
Bool_type
mili_reader_get_st_descriptors( Analysis *analy, int dbid );


/* TAG( mili_reader_get_state )
 *
 * Move to a particular state in the Mili database and update
 * nodal positions for the mesh.
 */
Bool_type
mili_reader_get_state( Analysis *analy, int state_no, State2 *p_st, State2 **pp_new_st, int *state_qty );


#endif

#endif