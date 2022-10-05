#ifndef PARALLEL_READ_H
#define PARALLEL_READ_H

#include "griz_config.h"

#ifdef HAVE_PARALLEL_READ

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include "misc.h"

#define OK 0
#define NOT_OK 1

/* Python program name */
const char* PROGRAM_NAME = "Griz-mili-python-reader";

/* Module to import for parallel read features */
const char* MILI_READER_MODULE_NAME = "mili.grizinterface";

/* Mili Reader GrizInterface function names as constants */
const char* OPEN_DATABASE = "open_griz_interface";
const char* PARAMETER_GETTER = "parameters";
const char* PARAM_WILDCARD_SEARCH = "parameter_wildcard_search";
const char* STATE_MAPS_GETTER = "state_maps";
const char* QUERY_FUNCTION = "query";
const char* RELOAD_FUNCTION = "reload";
const char* FREE_NODE_DATA = "free_node_data";

/* Griz Struct function calls */
const char* GET_GEOM_CALL = "griz_get_geom_call";
const char* GET_ST_DESCRIPTORS_CALL = "griz_get_st_descriptors_call";
const char* GET_STATE_CALL = "griz_get_state_call";

const int QTY_CONNECTS[] = { 0, 0, 2, 3, 4, 4, 4, 5, 6, 8, 0, 0, 0, 1, 10 };

/* Error messages */
const char * MR_FAILED_IMPORT = "\nFailed to import the Mili Reader module.\n";
const char * MR_DEFAULT_FAILED = "\n%s calling Mili Python (%s) Failed.\n";


/*****************************************************************
 * TAG( check_running_on_login_node )
 *
 * Check if Griz is running on a login node.
 */
Bool_type
check_running_on_login_node()
{
    const char* SLURM = "SLURM_JOBID";
    const char* LSF = "LSF_JOBID";
    Bool_type on_login_node = !getenv(SLURM) && !getenv(LSF);
    return on_login_node;
}


/*****************************************************************
 * TAG( bin_search_index )
 *
 * Binary search to find index of the target element in an array.
 */
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
 * Given a specified function name from the mili module and a set of arguments, check that the function
 * exists in the mili module and call it. Returns the PyObject pointer returned from
 * calling the python function. If something goes wrong, returns Py_None.
 *
 */
PyObject*
call_mili_module_function(PyObject * mili_reader_module, const char *function_name, PyObject *py_Arglist){
    PyObject *py_Dict,
             *py_Func;
    PyObject *py_ReturnValue = Py_None;
    Py_INCREF(Py_None);

    /* Check that mili_reader_module is not NULL */
    if(mili_reader_module != NULL){
        /* Get __dict__ for mili reader module */
        py_Dict = PyModule_GetDict(mili_reader_module);
        if(py_Dict != NULL){
            /* Check for requested function and call it */
            py_Func = PyDict_GetItemString(py_Dict, function_name);
            if(py_Func != NULL && PyCallable_Check(py_Func)){
                py_ReturnValue = PyObject_CallObject(py_Func, py_Arglist);
                Py_DECREF(Py_None);
            }
        }
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
 * Given a specified function name from the GrizInterface that takes no
 * arguments, call the function and return the result.
 */
PyObject*
call_mili_reader_function_noargs(PyObject * mili_db, const char * function_name)
{
    PyObject * py_FuncName = string_to_pyobject( function_name );
    PyObject * py_Result = PyObject_CallMethodObjArgs( mili_db, py_FuncName, NULL );
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
    Py_INCREF(Py_None);
    if( PyObject_HasAttrString(object, attr) ){
        py_Attribute = PyObject_GetAttrString( object, attr );
        Py_DECREF(Py_None);
    }
    return py_Attribute;
}


/* TAG( get_pyobject_attribute_as_string )
 *
 * Retrieve a specific attribute from a PyObject * as a string.
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
 */
int
get_pyobject_attribute_as_int( PyObject* object, const char* attr )
{
    PyObject * py_Attribute = PyObject_GetAttrString( object, attr );
    int value = PyLong_AsLong(py_Attribute);
    return value;
}

/* TAG( get_pyobject_attribute_as_double )
 *
 * Retrieve a specific attribute from a PyObject * as a double.
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
    int value;
    for( i = 0; i < size; i++ ){
        py_LongValue = PySequence_GetItem(py_List, i);
        value = PyLong_AsLong(py_LongValue);
        values[i] = value;
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
    float flt;
    for( i = 0; i < size; i++ ){
        py_FloatValue = PySequence_GetItem(py_List, i);
        flt = (float) PyFloat_AsDouble(py_FloatValue);
        values[i] = flt;
    }
}



/* TAG( mili_reader_get_string )
 *
 * Query a string parameter from the mili reader.
 */
int
mili_reader_get_string(PyObject * py_Parameters, char *parameter_name, char *str){
    char* res;
    PyObject * py_Parameter;
    PyObject * py_Name;

    /* Lookup parameter in mili and convert to string */
    py_Name = string_to_pyobject(parameter_name);
    if( PyDict_Contains( py_Parameters, py_Name) ){
        py_Parameter = PyDict_GetItem( py_Parameters, py_Name );
        /* Get string representation of parameter */
        if(py_Parameter != NULL){
            res = pyobject_as_string( py_Parameter );
            if ( res == NULL )
                return NOT_OK;
        }
    }
    // Copy out the result string
    strcpy(str, res);
    return OK;
}


/* TAG( mili_reader_get_double )
 *
 * Query a double parameter from the mili reader.
 */
int
mili_reader_get_double(PyObject * py_Parameters, char *parameter_name, double *result){
    PyObject * py_MiliParameter;
    PyObject * py_Name;
    /* Lookup parameter in mili and convert to double */
    py_Name = string_to_pyobject(parameter_name);
    if( PyDict_Contains( py_Parameters, py_Name) ){
        py_MiliParameter = PyDict_GetItem( py_Parameters, py_Name );
        /* Get double representation of parameter */
        if(py_MiliParameter != NULL){
            double value = PyFloat_AsDouble(py_MiliParameter);
            *result = value;
            return OK;
        }
    }
    return NOT_OK;
}

#endif // HAVE_PARALLEL_READ

#endif