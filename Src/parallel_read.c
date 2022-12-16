/*
 *  IO wrappers for parallel reader using Mili Python Reader
 *
 */
#include "griz_config.h"

#ifdef HAVE_PARALLEL_READ

#include <time.h>
#include "parallel_read.h"
#include "viewer.h"

#ifndef MILI_READER_TIMING
#define MILI_READER_TIMING
#endif


/* TAG( griz_python_setup )
 *
 * Setup Griz Parallel read using the Mili Reader by setting up path to mili reader source code
 * and correct python packages. Import the mili reader module.
 */
extern int
griz_python_setup(Analysis * analy)
{
    PyObject * py_SrcPath;
    PyObject * py_SitePackagesPath;

    /* Default to Null */
    analy->py_MiliDB = NULL;
    analy->py_MiliReaderModule = NULL;

    /* Not necessary, but recommended */
    wchar_t *program_name = Py_DecodeLocale(PROGRAM_NAME, NULL);
    Py_SetProgramName(program_name);

    /* Initialize Python */
    Py_Initialize();

    /* Get correct path to mili reader source. */
    py_SrcPath = string_to_pyobject( MR_SRC );
    
    /* Get correct path to mili reader's bin/site-packages */
    py_SitePackagesPath = string_to_pyobject( MR_BIN );

    /* Update Python Path */
    PyObject *py_SysPath = PySys_GetObject("path");
    PyList_Insert( py_SysPath, 0, py_SitePackagesPath );
    PyList_Insert( py_SysPath, 0, py_SrcPath );

    /* Import the Mili Reader */
    analy->py_MiliReaderModule = PyImport_ImportModule(MILI_READER_MODULE_NAME);
    /* Check that we successfully found the reader */
    if( analy->py_MiliReaderModule == NULL ){
        fprintf(stderr, MR_FAILED_IMPORT);
        return FALSE;
    }

    return TRUE;
}


/* TAG( mili_reader_db_open )
 *
 * Open the specified plot files with the mili reader and set mili_db to reference the newly created
 * python reader object.
 *
 * NOTE: The argument p_dbid is only necessary to match the function signature of mili_db_open/taurus_db_open.
         This allows the analysis struct to store a pointer to these function depending on the db type.
 */
extern int
mili_reader_db_open( char *path_root, int *p_dbid )
{
    PyObject *py_Arglist;
    PyObject *py_PlotFilePath;
    Analysis * p_analysis;

#ifdef MILI_READER_TIMING
    long start, end;
    start = prec_timer();
#endif

    /* Get the global analysis pointer */
    p_analysis= get_analy_ptr();

    /* Build argument list for function call */
    py_Arglist = PyTuple_New(4);

    /* Get plot file path */
    py_PlotFilePath = Py_BuildValue("s", path_root);

    /* Create tuple of arguments for mili reader */
    PyTuple_SetItem(py_Arglist, 0, py_PlotFilePath); // base_filename
    PyTuple_SetItem(py_Arglist, 1, PyList_New(0));   // procs = []
    PyTuple_SetItem(py_Arglist, 2, Py_False);         // suppress_parallel = True
    PyTuple_SetItem(py_Arglist, 3, Py_True);        // experimental = False

    /* Call open_database function */
    p_analysis->py_MiliDB = call_mili_module_function(p_analysis->py_MiliReaderModule, OPEN_DATABASE, py_Arglist);
    if(p_analysis->py_MiliDB == NULL){
        return NOT_OK;
    }

    /* Get the number of processors for the problem */
    p_analysis->proc_count = get_pyobject_attribute_as_int( p_analysis->py_MiliDB, "processor_count" );
    if( p_analysis->proc_count == -1 )
        return NOT_OK;

    /* Just set db id to 0*/
    p_dbid = 0;

#ifdef MILI_READER_TIMING
    end = prec_timer();
    printf("\n[mili_reader_db_open] elapsed = %ldms\n", (end-start));
#endif

    return OK;
}


/* TAG( build_query_args )
 *
 * Construct a PyObject Tuple containing all the arguments for a query to the mili reader.
 */
PyObject*
build_query_args( PyObject* py_Svar, PyObject* py_Class, PyObject* py_Labels, PyObject* py_States, PyObject* py_Ipts )
{
    PyObject * py_Args = PyTuple_New(7);
    PyTuple_SetItem( py_Args, 0, py_Svar );
    PyTuple_SetItem( py_Args, 1, py_Class );
    PyTuple_SetItem( py_Args, 2, Py_None );
    PyTuple_SetItem( py_Args, 3, py_Labels );
    PyTuple_SetItem( py_Args, 4, py_States );
    PyTuple_SetItem( py_Args, 5, py_Ipts );
    PyTuple_SetItem( py_Args, 6, Py_None );
    return py_Args;
}


/* TAG( build_query_kwargs )
 *
 * Construct a PyObject Dict containing all the keyword arguments for a query to the mili reader.
 */
PyObject*
build_query_kwargs( PyObject* py_Subrecord )
{
    PyObject * py_Kwargs = PyDict_New();
    PyDict_SetItemString( py_Kwargs, "output_object_labels", Py_False );
    PyDict_SetItemString( py_Kwargs, "subrec", py_Subrecord );
    return py_Kwargs;
}


/* TAG( get_query_function )
 *
 * Retrieve PyObject* for the Mili reader query function.
 */
PyObject*
get_query_function( Analysis * p_analysis )
{
    PyObject * py_Func, * py_FuncName;
    py_FuncName = string_to_pyobject( QUERY_FUNCTION );
    py_Func = PyObject_GetAttr( p_analysis->py_MiliDB, py_FuncName );
    Py_DECREF( py_FuncName );
    return py_Func;
}


/* TAG( mili_reader_get_svar_argument )
 *
 * Build PyObject* of list of svar names to query.
 */
PyObject*
mili_reader_get_svar_argument( int qty_results, char ** result_names )
{
    int i;
    PyObject * py_SvarNames;
    PyObject * py_TempStr;

    /* Build list of state variable names. */
    py_SvarNames = PyList_New( qty_results );
    for( i = 0 ; i < qty_results; i++ )
    {
        py_TempStr = string_to_pyobject( result_names[i] );
        PyList_SetItem( py_SvarNames, i, py_TempStr );
    }

    return py_SvarNames;
}


/* TAG( mili_reader_get_labels_argument )
 *
 * Build PyObject* of list of labels to query.
 */
PyObject*
mili_reader_get_labels_argument( int qty_labels, int * labels_list )
{
    PyObject * py_Labels;

    if( qty_labels == 0 || labels_list == NULL )
    {
        py_Labels = Py_None;
        Py_INCREF( Py_None );
    }
    else
    {
        int i;
        py_Labels = PyList_New( qty_labels );
        for( i = 0; i < qty_labels; i++ )
        {
            PyList_SetItem( py_Labels, i, PyLong_FromLong(labels_list[i]) );
        }
    }

    return py_Labels;
}


/* TAG( mili_reader_get_states_argument )
 *
 * Build PyObject* of list of states to query.
 */
PyObject*
mili_reader_get_states_argument( int analysis_total_states, int start_state, int end_state )
{
    int i, j;
    int qty_states;
    PyObject * py_States;
    
    qty_states = (end_state - start_state) + 1;
    if( qty_states == analysis_total_states )
    {
        py_States = Py_None;
        Py_INCREF( Py_None );
    }
    else
    {
        py_States = PyList_New( qty_states );
        j = 0;
        for( i = start_state; i <= end_state; i++ )
        {
            PyList_SetItem( py_States, j++, PyLong_FromLong(i) );
        }
    }

    return py_States;
}


/* TAG( mili_reader_get_ipts_argument )
 *
 * Build PyObject* for integration point to query.
 */
PyObject*
mili_reader_get_ipts_argument( int ipt )
{
    PyObject * py_Ipt;

    if( ipt == -1 )
    {
        py_Ipt = Py_None;
        Py_INCREF( Py_None );
    }
    else
    {
        py_Ipt = PyLong_FromLong( ipt );
    }

    return py_Ipt;
}


/* TAG( mili_reader_get_subrec_argument )
 *
 * Build PyObject* name of Subrecord to query.
 */
PyObject*
mili_reader_get_subrec_argument( char * subrecord_name )
{
    PyObject * py_Subrecord;

    if( subrecord_name == NULL )
    {
        py_Subrecord = Py_None;
        Py_INCREF( Py_None );
    }
    else
    {
        py_Subrecord = string_to_pyobject( subrecord_name );
    }

    return py_Subrecord;
}


/* TAG( is_previous_query )
 *
 * Check if the current query matches the previous query.
 */
Bool_type is_previous_query( Analysis * p_analysis, char* result_name, char* class_name, int ipt, int start_state, int end_state )
{
    Bool_type same_query;
    same_query = ( p_analysis->py_PrevQuery != NULL
                   && strcmp(p_analysis->prev_query_result_name, result_name) == 0
                   && strcmp(p_analysis->prev_query_class_name, class_name) == 0
                   && p_analysis->prev_query_ipt == ipt
                   && p_analysis->prev_query_start_state <= start_state
                   && p_analysis->prev_query_end_state >= end_state );
    return same_query;
}


/* TAG( store_query )
 *
 * Store the query for resuse.
 */
void store_query( Analysis * p_analysis, PyObject * py_Query, char* result_name, char* class_name, int ipt, int start_state, int end_state )
{
    /* Update the stored previous query */
    if( p_analysis->py_PrevQuery != NULL )
        Py_DECREF( p_analysis->py_PrevQuery );
    p_analysis->py_PrevQuery = py_Query;
    Py_INCREF( p_analysis->py_PrevQuery );

    /* Update previous query arguments */
    p_analysis->prev_query_result_name = (char*) strdup(result_name);
    p_analysis->prev_query_class_name = (char*) strdup(class_name);
    p_analysis->prev_query_start_state = start_state;
    p_analysis->prev_query_end_state = end_state;
    p_analysis->prev_query_ipt = ipt;
}


/* TAG( mili_reader_query )
 *
 * Query a result using the mili reader and return the PyObject pointer to
 * the result dictionary
 */
PyObject*
mili_reader_query( char* result, char* class_name, char* subrecord_name, int ipt, int start_state, int end_state )
{
    int i;
    /* Arguments we need for mili reader query function */
    PyObject * py_SvarNames,
             * py_ClassName,
             * py_Labels,
             * py_States,
             * py_Subrecord,
             * py_Ipt;
    PyObject * py_Func,
             * py_Args,
             * py_Kwargs;
    PyObject * py_QueryReturn;

    /* Get analysis pointer */
    Analysis * p_analysis = get_analy_ptr();

    /* Check if this is a preloaded result for plot/outth */
    if( p_analysis->py_PreloadedResult && strcmp( result, "nodpos" ) != 0 )
    {
        /* Just use preloaded result */
        py_QueryReturn = p_analysis->py_PreloadedResult;
        Py_INCREF( py_QueryReturn );
    }
    /* Check if this is a repeated query */
    else if( is_previous_query(p_analysis, result, class_name, ipt, start_state, end_state) )
    {
        /* Just use saved result */
        py_QueryReturn = p_analysis->py_PrevQuery;
        Py_INCREF( py_QueryReturn );
    }
    /* If not, perform query */
    else
    {
        /* Generate arguments for the Query */
        py_SvarNames = mili_reader_get_svar_argument( 1, &result );
        py_ClassName = string_to_pyobject( class_name );
        py_Labels = mili_reader_get_labels_argument( 0, NULL );
        py_States = mili_reader_get_states_argument( p_analysis->state_count, start_state, end_state );
        py_Subrecord = mili_reader_get_subrec_argument( subrecord_name );
        py_Ipt = mili_reader_get_ipts_argument( ipt );

        /* Get function name and call query */
        py_Func = get_query_function( p_analysis );
        py_Args = build_query_args( py_SvarNames, py_ClassName, py_Labels, py_States, py_Ipt );
        py_Kwargs = build_query_kwargs( py_Subrecord );
        py_QueryReturn = PyObject_Call( py_Func, py_Args, py_Kwargs );

        /* If we got back a result, save it for possible repeated use */
        if( py_QueryReturn != NULL )
            store_query( p_analysis, py_QueryReturn, result, class_name, ipt, start_state, end_state );

        /* Free python refereces */
        Py_DECREF( py_Args );
        Py_DECREF( py_Kwargs );
        Py_DECREF( py_Func );
    }

    return py_QueryReturn;
}


/* TAG( result_dictionary_unwrap_data )
 *
 * Unwrap mili reader result dictionary to get to list of data.
 */
PyObject*
result_dictionary_unwrap_data(PyObject * py_Dict, char * result_name, char * subrecord_name)
{
    const char * data_str = "data";
    PyObject * py_Result;

    py_Result = PyDict_GetItemString( py_Dict, result_name);
    py_Result = PyDict_GetItemString( py_Result, data_str);
    py_Result = PyDict_GetItemWithError( py_Result, string_to_pyobject(subrecord_name));

    return py_Result;
}


/* TAG( result_dictionary_unwrap_layout )
 *
 * Unwrap mili reader result dictionary to get to list of object ids.
 */
PyObject*
result_dictionary_unwrap_layout(PyObject * py_Dict, char * result_name, char * subrecord_name)
{
    const char * layout_str = "layout";
    PyObject * py_Result;

    py_Result = PyDict_GetItemString( py_Dict, result_name);
    py_Result = PyDict_GetItemString( py_Result, layout_str);
    py_Result = PyDict_GetItemWithError( py_Result, string_to_pyobject(subrecord_name));

    return py_Result;
}


/* TAG( result_dictionary_unwrap_states )
 *
 * Unwrap mili reader result dictionary to get to list of states.
 */
PyObject*
result_dictionary_unwrap_states(PyObject * py_Dict, char * result_name)
{
    const char * layout_str = "layout";
    const char * states_str = "states";
    PyObject * py_Result;

    py_Result = PyDict_GetItemString( py_Dict, result_name);
    py_Result = PyDict_GetItemString( py_Result, layout_str);
    py_Result = PyDict_GetItemWithError( py_Result, string_to_pyobject(states_str));

    return py_Result;
}


/* TAG( get_state_index )
 *
 * Get the index of the state's results from Dictionary
 */
int
get_state_index( PyObject * py_Dict, char * result_name, int state_no )
{
    int j;
    int st_qty;
    int state_idx = 0;
    PyObject * py_States;
    
    py_States = result_dictionary_unwrap_states( py_Dict, result_name );
    if( py_States != NULL)
    {
        /* Get quantity of states in the query */
        st_qty = PySequence_Length( py_States );
        if( st_qty == 1 )
        {
            state_idx = 0;
        }
        else
        {
            /* Find correct index of state */
            int * states = NEW_N( int, st_qty, "List of states in query" );
            integer_pointer_from_pyobject( py_States, st_qty, states );
            for( j = 0; j < st_qty; j++ )
            {
                if( states[j] == state_no )
                {
                    state_idx = j;
                    break;
                }
            }
        }
    }

    return state_idx;
}


/* TAG( extract_node_result )
 *
 * Extract the results from a PyObject pointer to Mili reader query result
 * that contains results for nodes.
 */
int
extract_node_result( PyObject * py_QueryReturn, char * result_name, int queried_state, char * subrecord_name, int comp_qty, float * result_buffer )
{
    int i, j, k;
    int state_idx;
    int obj, elem_qty;
    int out_idx, in_idx;
    PyObject * py_ProcResult;
    PyObject * py_Dict;
    int * index_map;

    PyObject * py_ByteArray;
    unsigned char * bytearray;
    float * floatValues;

    /* Get global analysis pointer */
    Analysis * p_analysis = get_analy_ptr();

    state_idx = -1;
    /* Loop over processors and get each result */
    for( i = 0; i < p_analysis->proc_count; i++ )
    {
        /* Get results from the ith processor */
        py_Dict = PySequence_ITEM( py_QueryReturn, i );
        /* Get results for current subrecord */
        py_ProcResult = result_dictionary_unwrap_data( py_Dict, result_name, subrecord_name );
        /* If results exist... */
        if( py_ProcResult != NULL )
        {
            /* Determine index of state data in query (Only done once) */
            if( state_idx == -1 )
            {
                state_idx = get_state_index( py_Dict, result_name, queried_state );
            }
            py_ProcResult = PySequence_ITEM( py_ProcResult, state_idx );

            /* Get object ids and element count */
            elem_qty = PySequence_Length( py_ProcResult );

            /* Get map from subrec node ids to global positions */
            index_map = MESH_P(p_analysis)->index_map->node_map[i];

            py_ByteArray = PyByteArray_FromObject( py_ProcResult ); 
            bytearray = PyByteArray_AS_STRING( py_ByteArray );
            floatValues = (float*) bytearray;

            /*  Extract and properly order results. */
            for( j = 0; j < elem_qty; j++)
            {
                in_idx = j * comp_qty;
                out_idx = (index_map[j] * comp_qty);
                for( k = 0; k < comp_qty; k++ )
                    result_buffer[out_idx+k] = floatValues[in_idx+k];
            }
            Py_DECREF( py_ByteArray );
        }
    }

    return OK;
}


/* TAG( extract_element_result )
 *
 * Extract the results from a PyObject pointer to Mili reader query result
 * that contains results for an element class (Not M_NODE).
 */
int
extract_element_result( PyObject * py_QueryReturn, char * result_name, int queried_state, char * subrecord_name, int comp_qty, float * result_buffer )
{
    int i, j, k;
    int state_idx;
    int out_idx;
    int elem_qty;
    int pos;
    int qty_values;
    PyObject * py_ProcResult;
    PyObject * py_Dict;

    PyObject * py_ByteArray;
    unsigned char * bytearray;
    float * floatValues;
    float * out_loc;

    /* Get global analysis pointer */
    Analysis * p_analysis = get_analy_ptr();

    state_idx = -1;
    out_idx = 0;
    pos = 0;
    /* Loop over processors and get each result */
    for( i = 0; i < p_analysis->proc_count; i++ )
    {
        py_Dict = PySequence_ITEM( py_QueryReturn, i );
        py_ProcResult = result_dictionary_unwrap_data( py_Dict, result_name, subrecord_name );
        if( py_ProcResult != NULL )
        {
            /* Determine index of state data in query (Only done once) */
            if( state_idx == -1 )
            {
                state_idx = get_state_index( py_Dict, result_name, queried_state );
            }
            py_ProcResult = PySequence_ITEM( py_ProcResult, state_idx );

            qty_values = PySequence_Length( py_ProcResult ) * comp_qty;

            py_ByteArray = PyByteArray_FromObject( py_ProcResult ); 
            bytearray = PyByteArray_AS_STRING( py_ByteArray );
            floatValues = (float*) bytearray;

            out_loc = result_buffer + pos;
            memcpy( out_loc, floatValues, qty_values * sizeof(float) );
            pos += qty_values;

            Py_DECREF( py_ByteArray );
        }
    }

    return OK;
}


/* TAG( extract_element_subset_result )
 *
 * Extract the results from a PyObject pointer to Mili reader query result
 * that contains results for an element class (Not M_NODE).
 */
int
extract_element_subset_result( PyObject * py_QueryReturn, char * result_name, int queried_state, char * subrecord_name,
                               int elem_class_index, int comp_qty, float * result_buffer )
{
    int i, j, k;
    int subrec_idx;
    int obj;
    int out_idx, in_idx;
    int elem_qty;
    PyObject * py_ProcResult;
    PyObject * py_Results;
    PyObject * py_Value;
    PyObject * py_Dict;
    PyObject * py_Layout;
    int state_idx;
    int ** index_map;
    int * object_ids;
    int cur_state, srec_id, subrec;
    int subrecord_obj_qty;
    int qty_series;
    State_rec_obj * p_state_rec;
    Subrec_obj * p_subrec;

    PyObject * py_ByteArray;
    unsigned char * bytearray;
    float * floatValues;

    /* Get global analysis pointer */
    Analysis * p_analysis = get_analy_ptr();

    cur_state = p_analysis->cur_state;
    srec_id = p_analysis->state_srec_fmt_ids[cur_state];
    p_state_rec = p_analysis->srec_tree + srec_id;
    qty_series = p_state_rec->series_qty;

    subrec = p_analysis->cur_result->subrecs[p_analysis->result_index];
    p_subrec = p_state_rec->subrecs + subrec;
    subrecord_obj_qty = p_subrec->subrec.qty_objects;

    /* Get map from object ids to global array indexes */
    index_map = MESH_P( p_analysis )->index_map->elem_map[elem_class_index];

    /* qty_series may be more than necessary, but always enough */
    object_ids = NEW_N( int, qty_series, "List of object ids" );

    out_idx = 0;
    state_idx = -1;

    /* Loop over processors and get each result */
    for( i = 0; i < p_analysis->proc_count; i++ )
    {
        py_Dict = PySequence_ITEM( py_QueryReturn, i );
        py_ProcResult = result_dictionary_unwrap_data( py_Dict, result_name, subrecord_name );
        if( py_ProcResult != NULL )
        {
            /* Determine index of state data in query (Only done once) */
            if( state_idx == -1 )
            {
                state_idx = get_state_index( py_Dict, result_name, queried_state );
            }
            py_ProcResult = PySequence_ITEM( py_ProcResult, state_idx );

            /* Get object ids and element count */
            py_Layout = result_dictionary_unwrap_layout( py_Dict, result_name, subrecord_name );
            elem_qty = PySequence_Length( py_Layout );
            integer_pointer_from_pytuple( py_Layout, elem_qty, object_ids );

            py_ByteArray = PyByteArray_FromObject( py_ProcResult ); 
            bytearray = PyByteArray_AS_STRING( py_ByteArray );
            floatValues = (float*) bytearray;

            /*  Extract and properly order results. */
            for( j = 0; j < elem_qty; j++)
            {
                obj = object_ids[j];
                for( subrec_idx = 0; subrec_idx < subrecord_obj_qty; subrec_idx++ )
                {
                    if( p_subrec->object_ids[subrec_idx] == index_map[i][obj] )
                        break;
                }

                in_idx = j * comp_qty;
                out_idx = subrec_idx * comp_qty;
                for( k = 0; k < comp_qty; k++ )
                    result_buffer[out_idx+k] = floatValues[in_idx+k];
            }
        }
    }

    if( object_ids )
        free( object_ids );

    return OK;
}


/* TAG( mili_reader_extract_query_data )
 *
 * Call proper helper function to extract the results from a PyObject pointer to Mili reader query result.
 */
int
mili_reader_extract_query_data( PyObject * py_QueryReturn, char * result_name, int queried_state,
                                char * subrecord_name, Bool_type node_result,
                                int elem_class_idx, int comp_qty, void* data )
{
    Analysis * p_analysis = get_analy_ptr();
    float * result_buf = (float*) data;

    if( node_result )
    {
        extract_node_result( py_QueryReturn, result_name, queried_state, subrecord_name, comp_qty, result_buf );
    }
    else if( !node_result && p_analysis->py_PreloadedResult)
    {
        /* Most likely only need to extract a subset of elements */
        extract_element_subset_result( py_QueryReturn, result_name, queried_state, subrecord_name,
                                       elem_class_idx, comp_qty, result_buf );
    }
    else
    {
        /* Need to extract all elements */
        extract_element_result( py_QueryReturn, result_name, queried_state, subrecord_name, comp_qty, result_buf );
    }

    return OK;
}



/* TAG( mili_reader_get_results )
 *
 * Query the Mili reader for some results
 */
extern int
mili_reader_get_results( int dbid, int state, int subrec_id, int qty, char **results, void *data )
{
    int i;
    int rval;
    int srec_id;
    int elem_class_idx;
    int ipt, comp_qty;
    Bool_type node_result;
    char* subrecord_name;
    char* class_name;
    State_rec_obj* p_sro;
    Subrec_obj* p_subrec;
    State_variable * p_sv;
    Htable_entry* p_hte;
    PyObject * py_QueryReturn;
    PyObject * py_Result;

    /* Get analysis pointer */
    Analysis * p_analysis = get_analy_ptr();

#ifdef MILI_READER_TIMING
    long start, end;
    long istart, iend;
    start = prec_timer();
#endif

    /* Get srec id and subrecord name */
    srec_id = p_analysis->state_srec_fmt_ids[state-1];
    p_sro = p_analysis->srec_tree + srec_id;
    p_subrec = p_sro->subrecs + subrec_id;
    subrecord_name = p_subrec->subrec.name;
    class_name = p_subrec->p_object_class->short_name;

    /* Look up state variable to get comp_qty */
    comp_qty = 1;
    rval = htable_search( p_analysis->st_var_table, results[0], FIND_ENTRY, &p_hte );
    if( rval == OK )
    {
        p_sv = (State_variable*) p_hte->data;
        comp_qty = p_sv->vec_size;
    }

    /* Get the integration point */
    ipt = get_subrecord_integration_point_num( p_subrec );

    /* This technically breaks if reading multiple results, but that currently never happens */
    for( i = 0; i < qty; i++ )
    {
        istart = prec_timer();
        py_QueryReturn = mili_reader_query( results[i], class_name, NULL, ipt, state, state );
        iend = prec_timer();
        printf("\tJust query = %ld\n", (iend-istart));
        if( py_QueryReturn == NULL )
            return NOT_OK;

        /* Extract/reorder the results from the Python dictionary structure */
        node_result = (p_subrec->p_object_class->superclass == G_NODE);
        elem_class_idx = p_subrec->p_object_class->elem_class_index;
        istart = prec_timer();
        rval = mili_reader_extract_query_data( py_QueryReturn, results[i], state, subrecord_name, node_result,
                                               elem_class_idx, comp_qty, data );
        iend = prec_timer();
        printf("\tExtracting results = %ld\n", (iend-istart));
        
        Py_DECREF( py_QueryReturn );
    }

#ifdef MILI_READER_TIMING
    end = prec_timer();
    //printf("[mili_reader_get_results] - [%s,%s,%s] - elapsed = %ldms\n",
    //        results[0], class_name, subrecord_name, (end-start));
#endif

    return rval;
}


extern int
mili_reader_preload_primal_th( Analysis* analy, Result_mo_list_obj* p_rmlo, Subrec_obj* p_subrec,
                               int obj_qty, int *labels, int first_st, int last_st )
{
    int ipt;
    char * result_name;
    PyObject * py_SvarNames,
             * py_ClassName,
             * py_Labels,
             * py_States,
             * py_Subrecord,
             * py_Ipt;
    PyObject * py_Func,
             * py_Args,
             * py_Kwargs;
    PyObject * py_QueryReturn;

    /* Get integration point */
    ipt = get_subrecord_integration_point_num( p_subrec );

    /* Get result name */
    result_name = (char*) p_rmlo->result->name;

    /* Generate arguments for the Query */
    py_SvarNames = mili_reader_get_svar_argument( 1, &result_name );
    py_ClassName = string_to_pyobject( p_subrec->p_object_class->short_name );
    py_Labels = mili_reader_get_labels_argument( obj_qty, labels );
    py_States = mili_reader_get_states_argument( analy->state_count, first_st, last_st );
    py_Subrecord = mili_reader_get_subrec_argument( p_subrec->subrec.name );
    py_Ipt = mili_reader_get_ipts_argument( ipt );

    /* Call the query function */
    py_Func = get_query_function( analy );
    py_Args = build_query_args( py_SvarNames, py_ClassName, py_Labels, py_States, py_Ipt );
    py_Kwargs = build_query_kwargs( py_Subrecord );
    analy->py_PreloadedResult = PyObject_Call( py_Func, py_Args, py_Kwargs );

    /* Free python refereces */
    Py_DECREF( py_Args );
    Py_DECREF( py_Kwargs );
    Py_DECREF( py_Func );

    return OK;
}


/* TAG( mili_reader_get_derived_svar_argument )
 *
 * Build PyObject* of list of svars to query for derived result.
 */
PyObject*
mili_reader_get_derived_svar_argument( Analysis * p_analysis, Result* p_result, Subrec_obj* p_subrec )
{
    int qty_results;
    char ** result_names;
    const char * stress = "stress";
    const char * strain = "strain";
    char* stress_components[6] = {"sx", "sy", "sz", "sxy", "syz", "szx"};
    char* strain_components[6] = {"ex", "ey", "ez", "exy", "eyz", "ezx"};
    PyObject * py_SvarNames;

    int result_index = p_analysis->result_index;
    char * primal = p_result->primals[result_index][0];

    if( strcmp(primal, stress) == 0 && p_subrec->element_set )
    {
        result_names = stress_components;
        qty_results = 6;
    }
    else if( strcmp(primal, strain) == 0 && p_subrec->element_set )
    {
        result_names = strain_components;
        qty_results = 6;
    }
    else
    {
        result_names = &primal;
        qty_results = 1;
    }

    py_SvarNames = mili_reader_get_svar_argument( qty_results, result_names );

    return py_SvarNames;
}


extern int
mili_reader_preload_derived_th( Analysis* analy, Result_mo_list_obj* p_rmlo, Subrec_obj* p_subrec, 
                                int obj_qty, int *labels, int first_st, int last_st )
{
    int ipt;
    PyObject * py_SvarNames,
             * py_ClassName,
             * py_Labels,
             * py_States,
             * py_Subrecord,
             * py_Ipt;
    PyObject * py_Func,
             * py_Args,
             * py_Kwargs;
    PyObject * py_QueryReturn;

    /* Get integration point */
    ipt = get_subrecord_integration_point_num( p_subrec );

    py_SvarNames = mili_reader_get_derived_svar_argument( analy, p_rmlo->result, p_subrec );
    py_ClassName = string_to_pyobject( p_subrec->p_object_class->short_name );
    py_Labels = mili_reader_get_labels_argument( obj_qty, labels );
    py_States = mili_reader_get_states_argument( analy->state_count, first_st, last_st );
    py_Subrecord = mili_reader_get_subrec_argument( p_subrec->subrec.name );
    py_Ipt = mili_reader_get_ipts_argument( ipt );

    /* Call the query function */
    py_Func = get_query_function( analy );
    py_Args = build_query_args( py_SvarNames, py_ClassName, py_Labels, py_States, py_Ipt );
    py_Kwargs = build_query_kwargs( py_Subrecord );
    analy->py_PreloadedResult = PyObject_Call( py_Func, py_Args, py_Kwargs );

    /* Free python refereces */
    Py_DECREF( py_Args );
    Py_DECREF( py_Kwargs );
    Py_DECREF( py_Func );

    return OK;
}


/* TAG( mili_reader_load_stress_strain_components )
 *
 * Query the Mili reader for the stress or strain components from a specified subrecord.
 * NOTE: This function makes a lot of assumptions about how the data should be loaded in
 *       and really only designed to be used for the compute_es_... function in stress.c and strain.c
 */
extern int
mili_reader_load_stress_strain_components( int state, int subrec_id, Bool_type is_stress, Bool_type is_pressure, float** data )
{
    int i;
    int rval;
    int srec_id;
    int elem_class_idx;
    int ipt;
    Bool_type node_result;
    char* subrecord_name;
    char* class_name;
    char* result;
    State_rec_obj* p_sro;
    Subrec_obj* p_subrec;
    PyObject * py_Result;
    int num_components;
    char* stress_components[6] = {"sx", "sy", "sz", "sxy", "syz", "szx"};
    char* strain_components[6] = {"ex", "ey", "ez", "exy", "eyz", "ezx"};
    char** components_to_query;

    /* Arguments we need for mili reader query function */
    PyObject * py_SvarNames,
             * py_ClassName,
             * py_Labels,
             * py_States,
             * py_Ipt,
             * py_Subrecord;
    PyObject * py_Func,
             * py_Args,
             * py_Kwargs;
    PyObject * py_QueryReturn;

    /* Get analysis pointer */
    Analysis * p_analysis = get_analy_ptr();

    /* Get srec id and subrecord name */
    srec_id = p_analysis->state_srec_fmt_ids[state-1];
    p_sro = p_analysis->srec_tree + srec_id;
    p_subrec = p_sro->subrecs + subrec_id;
    subrecord_name = p_subrec->subrec.name;
    class_name = p_subrec->p_object_class->short_name;

    /* Get number/name of components we need to query */
    num_components = (is_pressure) ? 3 : 6;
    components_to_query = (is_stress) ? stress_components : strain_components;

    if( p_analysis->py_PreloadedResult != NULL )
    {
        py_QueryReturn = p_analysis->py_PreloadedResult;
        Py_INCREF( py_QueryReturn );
    }
    else
    {
        /* Get integration point */
        ipt = get_subrecord_integration_point_num( p_subrec );

        /* Generate arguments for the Query */
        py_SvarNames = mili_reader_get_svar_argument( num_components, components_to_query );
        py_ClassName = string_to_pyobject( class_name );
        py_Labels = mili_reader_get_labels_argument( 0, NULL );
        py_States = mili_reader_get_states_argument( p_analysis->state_count, state, state );
        py_Subrecord = mili_reader_get_subrec_argument( p_subrec->subrec.name );
        py_Ipt = mili_reader_get_ipts_argument( ipt );

        /* Call the query function */
        py_Func = get_query_function( p_analysis );
        py_Args = build_query_args( py_SvarNames, py_ClassName, py_Labels, py_States, py_Ipt );
        py_Kwargs = build_query_kwargs( py_Subrecord );
        py_QueryReturn = PyObject_Call( py_Func, py_Args, py_Kwargs );

        /* Free python refereces */
        Py_DECREF( py_Args );
        Py_DECREF( py_Kwargs );
        Py_DECREF( py_Func );
    }

    if( py_QueryReturn == NULL ){
        return OK;
    }

    /* Parse out stress and strain values and store in correct locations
     * Assumption:
     *      data[0][0:elem_qty] = sx/ex
     *      data[1][0:elem_qty] = sy/ey
     *      data[2][0:elem_qty] = sz/ez
     *      data[3][0:elem_qty] = sxy/exy
     *      data[4][0:elem_qty] = syz/eyz
     *      data[5][0:elem_qty] = szx/ezx
     */
    /* Loop over processors and get each result */

    for( i = 0; i < num_components; i++ ){
        result = (is_stress) ? stress_components[i] : strain_components[i];
        elem_class_idx = p_subrec->p_object_class->elem_class_index;
        rval = mili_reader_extract_query_data( py_QueryReturn, result, state, subrecord_name, FALSE,
                                               elem_class_idx, 1, (void*) data[i] );
    }

    return OK;
}


/* TAG( compute_label_blocking_data )
 *
 * For a given set of labels, compute the blocking data and return block_qty.
 */
int *
compute_label_blocking_data(int * labels, int qty, int * num_blocks)
{
    int i;
    int block_qty;
    int block_range_index;
    int start_of_block, end_of_block;
    int * block_range_ptr = NULL;

    block_qty = 1;
    for( i = 0; i < qty - 1; i++ ){
        if( labels[i+1] != labels[i]+1 )
            block_qty++;
    }
    *num_blocks = block_qty;

    block_range_ptr = (int*) malloc(block_qty * 2 * sizeof(int));
    if( block_range_ptr == NULL ){
        return NULL;
    }

    start_of_block = labels[0];
    end_of_block = labels[0];
    block_range_index = 0;

    /* Get Label blocks */
    for( i = 0; i < qty-1; i++ ){
        if( end_of_block+1 != labels[i+1] ){
            block_range_ptr[block_range_index++] = start_of_block;
            block_range_ptr[block_range_index++] = end_of_block;
            start_of_block = labels[i+1];
            end_of_block = labels[i+1];
        }
        else{
            end_of_block++;
        }
    }

    /* Add last block */
    block_range_ptr[block_range_index++] = start_of_block;
    block_range_ptr[block_range_index++] = end_of_block;

    return block_range_ptr;
}



/* TAG( mili_reader_get_geom )
 *
 * Load in the geometry (mesh definition) data from the mili reader
 */
extern int
mili_reader_get_geom( Analysis * analy )
{
    int i, j, k, l, m;
    int dims;
    int rval, target;
    int mesh_qty;
    int obj_qty;
    int obj_id;
    int num_classes;
    int elem_class_qty;
    int elem_class_count;
    int nodes_on_proc, elems_on_proc, pos, temp;
    int node_idx, out_idx, in_idx;
    int conn;
    int superclass;
    int block_qty;
    int block_index;
    int start_ident, stop_ident;
    int class_qty;
    int * block_range_ptr;
    int *node_labels, *elem_labels;
    int * index_map;
    Bool_type have_nodal;
    char * short_name;
    char * long_name;
    Mesh_data *p_md;
    Hash_table *p_ht;
    Bool_type idents_exist;
    Htable_entry *p_hte;
    List_head *p_lh;
    MO_class_data *node_class, *elem_class;
    MO_class_data **mo_classes, **htable_mo_classes;
    Elem_data * p_ed;
    Material_data * p_matd;
    PyObject * py_ElementClasses,
             * py_ElementClass,
             * py_ClassDataDict,
             * py_ProcData,
             * py_MOData,
             * py_CoordsTuple;
    PyObject * py_ProcLabels, * py_ProcConns, * py_ProcNodes;
    PyObject * py_ProcMats, * py_ProcParts;
    PyObject * py_ElemClassMats, * py_ElemClassParts;
    PyObject * py_ByteArray;
    unsigned char * bytearray;
    float * floatValues;
    int * intValues;
    PyObject * py_MatBytes, * py_PartBytes;
    int * mats, *parts;
            
    if( analy->mesh_table != NULL)
    {
        popup_dialog( WARNING_POPUP, "Mesh table pointer not NULL as initialization." );
        return -1;
    }

    int proc_count = analy->proc_count;

#ifdef MILI_READER_TIMING
    long start, end;
    long istart, iend;
    start = prec_timer();
#endif

    /* Single call to Mili reader to get all the data needed by mili_reader_get_geom */
    PyObject * py_GetGeomCallData = call_mili_reader_function_noargs( analy->py_MiliDB, GET_GEOM_CALL );
    if ( py_GetGeomCallData == NULL )
    {
        fprintf( stderr, "mili_reader_get_geom calling GET_GEOM_CALL returned NULL\n");
        return GRIZ_FAIL;
    }

    /* Split out data by processor for easier access later on */
    PyObject *py_AllMeshObjectClasses = get_pyobject_attribute( py_GetGeomCallData, "mo_classes" );
    PyObject *py_AllNodes = get_pyobject_attribute( py_GetGeomCallData, "nodes" );
    PyObject *py_AllLabels = get_pyobject_attribute( py_GetGeomCallData, "labels" );
    PyObject *py_AllConns = get_pyobject_attribute( py_GetGeomCallData, "connectivity" );
    PyObject *py_AllMats = get_pyobject_attribute( py_GetGeomCallData, "materials" );
    PyObject *py_AllParts = get_pyobject_attribute( py_GetGeomCallData, "parts" );

    /* Split out data by processor for easier access later on */
    PyObject **py_MOclass_by_proc = (PyObject**) malloc(proc_count * sizeof(PyObject*));
    PyObject **py_nodes_by_proc = (PyObject**) malloc(proc_count * sizeof(PyObject*));
    PyObject **py_labels_by_proc = (PyObject**) malloc(proc_count * sizeof(PyObject*));
    PyObject **py_conns_by_proc = (PyObject**) malloc(proc_count * sizeof(PyObject*));

    for( i = 0; i < proc_count; i++ )
    {
        /* Get data */
        py_MOclass_by_proc[i] = PySequence_ITEM( py_AllMeshObjectClasses, i);
        py_nodes_by_proc[i]   = PySequence_ITEM( py_AllNodes, i);
        py_labels_by_proc[i] = PySequence_ITEM( py_AllLabels, i );
        py_conns_by_proc[i] = PySequence_ITEM( py_AllConns, i );
    }

    /* Get dimensions of the database */
    dims = get_pyobject_attribute_as_int( analy->py_MiliDB, "dimensions" );
    if ( dims != 2 && dims != 3 )
    {
        popup_dialog( WARNING_POPUP, "mili_reader_get_geom() db dimension should be 2 or 3.");
        return GRIZ_FAIL;
    }
    if( dims == 2 )
        analy->limit_rotations = TRUE;
    analy->dimension = dims;

    /* For now we'll just assume is one. */
    mesh_qty = 1;
    analy->mesh_qty = mesh_qty;

    /* Allocate array of pointer to mesh geom hash tables. */
    analy->mesh_table = NEW_N( Mesh_data, mesh_qty, "Mesh data array" );

    for( i = 0; i < mesh_qty; i++ )
    {
        p_md = analy->mesh_table + i;
        p_ht = htable_create( 151 );
        p_md->class_table = p_ht;

        /* Create initial MO_class_data structs for each element class and sum obj_qty across processors */
        for( j = 0; j < proc_count; j++ )
        {
            // Get Mesh Object class data for 'j' processor
            py_MOData = py_MOclass_by_proc[j];
            if ( py_MOData == NULL )
                continue;
            py_MOData = PyDict_Values( py_MOData );
            num_classes = PySequence_Length( py_MOData );

            for( k = 0; k < num_classes; k++ )
            {
                py_ElementClass = PySequence_ITEM( py_MOData, k );
                short_name = get_pyobject_attribute_as_string( py_ElementClass, "short_name" );
                rval = htable_search( p_ht, short_name, ENTER_MERGE, &p_hte );
                /* New MO class */
                if( rval == OK )
                {
                    long_name = get_pyobject_attribute_as_string( py_ElementClass, "long_name" );
                    superclass = get_pyobject_attribute_as_int( py_ElementClass, "sclass" );
                    obj_qty = get_pyobject_attribute_as_int( py_ElementClass, "elem_qty" );
                    /* Record whether this element class has labels associated with it in the mili reader. */
                    idents_exist = get_pyobject_attribute_as_int( py_ElementClass, "idents_exist" );

                    /* Generate MO_class_data object for element class */
                    elem_class = NEW( MO_class_data, "Elem geom table entry" );
                    elem_class->mesh_id = i;
                    elem_class->short_name = short_name;
                    elem_class->long_name = long_name;
                    elem_class->superclass = superclass;
                    elem_class->qty = obj_qty;
                    /* Not really what this was meant for, but we only use it for part
                     * of this function and then update it later */
                    elem_class->labels_found = idents_exist;

                    p_hte->data = (void*) elem_class;

                }
                /* Existing mo_class */
                else
                {
                    /* Sum up object quantity over processors. Can skip mat, glob, and unit */
                    elem_class = (MO_class_data*) p_hte->data;
                    if( elem_class->superclass != G_MESH && elem_class->superclass != G_MAT && elem_class->superclass != G_UNIT )
                    {
                        obj_qty = get_pyobject_attribute_as_int( py_ElementClass, "elem_qty" );
                        elem_class->qty += obj_qty;
                    }
                }

                Py_DECREF( py_ElementClass );
            }
        }

        Py_DECREF( py_MOData );

        /* ----------------------------------------------------------------------------------------- */
        istart = prec_timer();

        /* Look up node class */
        rval = htable_search( p_ht, "node", FIND_ENTRY, &p_hte );
        node_class = (MO_class_data*) p_hte->data;

        /* Initialize mapping from processor to global */
        p_md->index_map = NEW( ProcessorToGlobalMap, "Proc to Global index map");
        p_md->index_map->node_count = (int*) malloc( proc_count * sizeof(int));
        p_md->index_map->node_map = (int **) malloc( proc_count * sizeof(int*));

        /* Populate data for nodal class */
        if ( rval == OK )
        {
            obj_qty = node_class->qty;
            /* Get all labels into a single list */
            node_labels = (int*) malloc(obj_qty * sizeof(int));
            pos = 0;
            /* For each processor */
            for( j = 0; j < proc_count; j++ )
            {
                /* Get node labels on processor */
                py_ProcLabels = PyDict_GetItemString( py_labels_by_proc[j], node_class->short_name );
                nodes_on_proc = PySequence_Length( py_ProcLabels );
                if( nodes_on_proc > 0 )
                {
                    /* Convert node labels to integer array */
                    integer_pointer_from_pyobject( py_ProcLabels, nodes_on_proc, node_labels+pos );
                    pos += nodes_on_proc;

                    /* Update arrays */
                    p_md->index_map->node_count[j] = nodes_on_proc;
                    p_md->index_map->node_map[j] = (int*) malloc( nodes_on_proc * sizeof(int) );
                }
                else
                {
                    /* Update arrays */
                    p_md->index_map->node_count[j] = 0;
                    p_md->index_map->node_map[j] = NULL;
                }
            }

            /* Sort the labels */
            qsort( node_labels, obj_qty, sizeof(int), compare_int);

            /* Remove duplicates and get real number of nodes */
            k = 0;
            for( j = 0; j < obj_qty-1; j++ )
            {
                if ( node_labels[j] != node_labels[j+1])
                    node_labels[k++] = node_labels[j];
            }
            node_labels[k++] = node_labels[obj_qty-1];

            obj_qty = k;
            node_class->qty = obj_qty;

            /* Create per processor mapping for nodes to index in global node list */
            for( j = 0; j < proc_count; j++ )
            {
                py_ProcLabels = PyDict_GetItemString( py_labels_by_proc[j], node_class->short_name );
                nodes_on_proc = p_md->index_map->node_count[j];
                if( nodes_on_proc > 0 )
                {
                    /* Store label numbers in array temporarily, will overwrite next. */
                    integer_pointer_from_pyobject( py_ProcLabels, nodes_on_proc, p_md->index_map->node_map[j] );

                    for( k = 0; k < nodes_on_proc; k++ )
                    {
                        target = p_md->index_map->node_map[j][k];
                        /* ----- Binary Search ----- */
                        p_md->index_map->node_map[j][k] = bin_search_index(target, node_labels, obj_qty, target-1);
                    }
                }
            }

            if( dims == 3 )
                node_class->objects.nodes3d = NEW_N( GVec3D, obj_qty, "3D node coord array");
            else
                node_class->objects.nodes2d = NEW_N( GVec2D, obj_qty, "2D node coord array");

            node_class->objects.nodes = NEW_N(float, obj_qty * 3, "Node Positions");

            node_class->data_buffer = NEW_N( float, obj_qty, "Class data buffer" );
            if( node_class->data_buffer == NULL )
                popup_fatal( "Unable to allocate data buffer on class load" );

            float * node_coords = node_class->objects.nodes;
    
            /* For each processor */
            for( j = 0; j < proc_count; j++ )
            {
                nodes_on_proc = p_md->index_map->node_count[j];

                PyErr_Clear();
                py_ByteArray = PyByteArray_FromObject( py_nodes_by_proc[j] );
                bytearray = PyByteArray_AS_STRING( py_ByteArray );
                floatValues = (float*) bytearray;

                index_map = p_md->index_map->node_map[j];
                /* Loop over all nodes and coordinates */
                for( k = 0; k < nodes_on_proc; k++ )
                {
                    node_idx = index_map[k];
                    out_idx = node_idx * dims;
                    in_idx = k * dims;

                    for( l = 0; l < dims; l++ )
                    {
                        node_coords[out_idx+l] = floatValues[in_idx+l];
                    }
                }
                Py_DECREF( py_ByteArray );
            }

            node_class->labels_max = -1;
            node_class->labels_min = MAXINT;

            node_class->labels = NEW_N( MO_class_labels, obj_qty, "Class Labels" );
            if ( node_class->labels == NULL )
                popup_fatal( "Unable to allocate labels on class load" );

            node_class->labels_index = NEW_N( int, obj_qty, "Class Labels index" );
            if ( node_class->labels_index == NULL )
                popup_fatal( "Unable to allocate labels index on class load" );

            node_class->labels_found = TRUE;

            for( obj_id = 0; obj_id < obj_qty; obj_id++ )
            {
                node_class->labels[obj_id].local_id = obj_id;
                node_class->labels[obj_id].label_num = node_labels[obj_id];
                if ( node_labels[obj_id] > node_class->labels_max )
                    node_class->labels_max = node_labels[obj_id];
                if ( node_labels[obj_id] < node_class->labels_min )
                    node_class->labels_min = node_labels[obj_id];
            }

            /* Sort the labels */
            qsort( node_class->labels, obj_qty, sizeof(MO_class_labels), mili_compare_labels);

            /* Create a mapping for the 1-n label index */
            for ( obj_id = 0; obj_id < obj_qty; obj_id++ )
            {
                node_class->labels_index[node_class->labels[obj_id].local_id] = obj_id; 
            }

            /* Construct blocking data for node labels */
            block_range_ptr = compute_label_blocking_data( node_labels, obj_qty, &block_qty);
            if( block_range_ptr == NULL )
            {
                fprintf(stderr, "mili_reader_load_nodal_data call compute_label_blocking_data\n");
                return GRIZ_FAIL;
            }

            /* Construct the label blocking table of contents */
            if( node_class->labels_found && block_qty > 0 && block_range_ptr )
            {
                node_class->label_blocking.block_qty = block_qty;
                node_class->label_blocking.block_total_objects = 0;
                node_class->label_blocking.block_min = MAXINT;
                node_class->label_blocking.block_max = MININT;
                node_class->label_blocking.block_objects = NEW_N( Label_block_data, block_qty, "Node Class Label Blocking Objects");

                block_index = 0;
                for( k = 0; k < block_qty; k++ )
                {
                    // Update min and max labels for this block 
                    if( block_range_ptr[block_index] < node_class->label_blocking.block_min )
                        node_class->label_blocking.block_min = block_range_ptr[block_index];
                    if( block_range_ptr[block_index] > node_class->label_blocking.block_max )
                        node_class->label_blocking.block_max = block_range_ptr[block_index];

                    node_class->label_blocking.block_objects[k].label_start = block_range_ptr[block_index++];
                    node_class->label_blocking.block_objects[k].label_stop = block_range_ptr[block_index++];
                    node_class->label_blocking.block_total_objects += ( node_class->label_blocking.block_objects[k].label_stop -
                        node_class->label_blocking.block_objects[k].label_start) + 1;
                }
            }

            if ( block_range_ptr )
            {
                free( block_range_ptr );
                block_range_ptr = NULL;
            }
            if ( node_labels )
            {
                free( node_labels );
                node_labels = NULL;
            }

            /* Add node class to list of element classes */
            p_lh = p_md->classes_by_sclass + M_NODE;
            mo_classes = (MO_class_data**) p_lh->list;
            mo_classes = (void*) RENEW_N( MO_class_data*, mo_classes, p_lh->qty, 1, "Extend Node sclass array");
            mo_classes[p_lh->qty] = node_class;
            p_lh->qty++;
            p_lh->list = (void*) mo_classes; 

        }
        /* No nodal class, create one... */
        else
        {
            popup_dialog( WARNING_POPUP, "Node object class not found \"node\"." );

            // No nodal or element data, so create a fake nodal class
            node_class->mesh_id = i;
            griz_str_dup( &node_class->short_name, " " );
            griz_str_dup( &node_class->long_name, " " );
            node_class->superclass = G_UNIT;
            node_class->elem_class_index = -1;
            node_class->qty = 0;

            // Allocate the data buffer for I/O and result derivation
            node_class->data_buffer = NEW_N( float, 10000, "Class data buffer" );
        }

        /* Keep a reference to node geometry handy. */
        p_md->node_geom = node_class;

        iend = prec_timer();
        printf("Nodal class load = %ld\n", (iend-istart));

        /* ----------------------------------------------------------------------------------------- */

        istart = prec_timer();

        /* For each element class load labels, connectivity, parts, materials, etc. */
        htable_get_data( p_ht, (void***) &htable_mo_classes, &elem_class_qty);

        /* Initialize arrays for storing element counts per processor and processor to global index mapping */
        p_md->index_map->elem_count = (int**) malloc( elem_class_qty * sizeof(int*) );
        p_md->index_map->elem_map = (int***) malloc( elem_class_qty * sizeof(int**) );
        for( j = 0; j < elem_class_qty; j++ )
        {
            p_md->index_map->elem_count[j] = (int*) malloc( proc_count * sizeof(int) );
            p_md->index_map->elem_map[j] = (int**) malloc( proc_count * sizeof(int*) );
        }

        /* Populate element class data. */
        for( j = 0; j < elem_class_qty; j++ )
        {
            elem_class = htable_mo_classes[j];

            /* Node class was already handled above. */
            if ( elem_class->superclass == G_NODE )
                continue;

            /* Skip if no elements */
            if ( elem_class->qty == 0 ){
                continue;
            }
            
            elem_class->elem_class_index = j;
            obj_qty = elem_class->qty;

            /* Allocate data buffer for I/O and result derivation */
            elem_class->data_buffer = NEW_N( float, obj_qty, "Class data buffer" );
            if ( elem_class->data_buffer == NULL )
                popup_fatal( "Unable to alloc data buffer on class load" );

            /* Allocate space for all element labels */
            elem_labels = NEW_N( int, obj_qty, "Element labels list" );
            if( elem_labels == NULL )
                popup_fatal( "Unable to alloc label array on class load" );
            
            /* Check if labels exist for this element class */
            idents_exist = elem_class->labels_found;
            elem_class->labels_found = FALSE;
            
            elems_on_proc = 0;
            pos = 0;
            temp = 0;
            for( k = 0; k < proc_count; k++ )
            {
                /* Get Labels on processor */
                py_ProcLabels = PyDict_GetItemString( py_labels_by_proc[k], elem_class->short_name );
                if( py_ProcLabels == NULL )
                {
                    p_md->index_map->elem_count[elem_class->elem_class_index][k] = 0;
                    p_md->index_map->elem_map[elem_class->elem_class_index][k] = NULL;
                }
                else
                {
                    elems_on_proc = PySequence_Length( py_ProcLabels );
                    /* Update arrays */
                    p_md->index_map->elem_count[elem_class->elem_class_index][k] = elems_on_proc;
                    p_md->index_map->elem_map[elem_class->elem_class_index][k] = (int*) malloc( elems_on_proc * sizeof(int) );

                    /* If labels exist */
                    if( idents_exist )
                    {
                        /* Read in element labels */
                        integer_pointer_from_pyobject( py_ProcLabels, elems_on_proc, (elem_labels+pos) );
                    }
                    else
                    {
                        /* Need to create our own labels -- 1->N */
                        for(l = 0; l < elems_on_proc; l++)
                        {
                            elem_labels[temp] = temp+1;
                            temp++;
                        }
                    }

                    /* Map proc elems to global indices */
                    if( elem_class->superclass == M_MESH || elem_class->superclass == M_MAT || elem_class->superclass == M_UNIT )
                    {
                        /* This is an awkward solution to MAT, MESH, and M_UNIT (maybe) being the same on
                        * all processors, can definitely make this better later */
                        for( l = 0; l < elems_on_proc; l++ )
                        {
                            p_md->index_map->elem_map[j][k][l] = l;
                        }
                    }
                    else
                    {
                        for( l = 0; l < elems_on_proc; l++ )
                        {
                            p_md->index_map->elem_map[j][k][l] = pos++;
                        }
                    }
                }
            }

            /* Special case for particle classes -- Elements need to be put in sorted order. */
            if( idents_exist && is_particle_class( analy, elem_class->superclass, elem_class->short_name ) )
            {
                /* Sort element labels */
                qsort( elem_labels, obj_qty, sizeof(int), compare_int );
                for( k = 0; k < proc_count; k++ )
                {
                    /* Get Labels on processor */
                    py_ProcLabels = PyDict_GetItemString( py_labels_by_proc[k], elem_class->short_name );
                    elems_on_proc = p_md->index_map->elem_count[j][k];
                    if( elems_on_proc > 0 )
                    {
                        /* Store label numbers in array temporarily, will overwrite next. */
                        integer_pointer_from_pyobject( py_ProcLabels, elems_on_proc, p_md->index_map->elem_map[j][k] );
                        
                        for( l = 0; l < elems_on_proc; l++ )
                        {
                            target = p_md->index_map->elem_map[j][k][l];
                            p_md->index_map->elem_map[j][k][l] = bin_search_index( target, elem_labels, obj_qty, target-1 );
                        }
                    }
                }
            }

            if( elem_class->superclass == M_MESH || elem_class->superclass == M_MAT || elem_class->superclass == M_UNIT )
            {

                elem_class->simple_start = elem_labels[0];
                elem_class->simple_stop = elem_labels[elem_class->qty-1];

                p_lh = p_md->classes_by_sclass + elem_class->superclass;
                mo_classes = (MO_class_data **) p_lh->list;
                mo_classes = (void *) RENEW_N( MO_class_data *, mo_classes, p_lh->qty, 1, "Extend classes_by_sclass array");
                mo_classes[p_lh->qty] = elem_class;
                p_lh->qty++;
                p_lh->list = (void *) mo_classes;
            }
            else
            {
                /* Generate Element data */
                int qty_conns = QTY_CONNECTS[elem_class->superclass];
                p_ed = NEW( Elem_data, "Element Conn Struct" );
                elem_class->objects.elems = p_ed;
                p_ed->nodes = NEW_N( int, obj_qty * qty_conns, "Element Connectivites" );
                p_ed->mat = NEW_N( int, obj_qty, "Element Materials" );
                p_ed->part = NEW_N( int, obj_qty, "Element parts" );
                p_ed->volume = NULL;

                /* Load in nodal connectivity, material, part number, etc. for each element */
                PyObject * py_M, * py_P;
                int * class_index_map;
                int * node_index_map;

                py_ElemClassMats = PyDict_GetItemString( py_AllMats, elem_class->short_name );
                py_ElemClassParts = PyDict_GetItemString( py_AllParts, elem_class->short_name );

                for( k = 0; k < proc_count; k++ )
                {
                    if( p_md->index_map->elem_count[j][k] > 0 )
                    {
                        elems_on_proc = p_md->index_map->elem_count[j][k];
                        class_index_map = p_md->index_map->elem_map[j][k];
                        node_index_map = p_md->index_map->node_map[k];

                        /* Get elements nodal connectivity */
                        py_ProcConns = PyDict_GetItemString( py_conns_by_proc[k], elem_class->short_name );
                        PyErr_Clear();
                        py_ByteArray = PyByteArray_FromObject( py_ProcConns );
                        bytearray = PyByteArray_AS_STRING( py_ByteArray );
                        intValues = (int*) bytearray;
                        for( l = 0; l < elems_on_proc; l++ )
                        {
                            in_idx = l*qty_conns;
                            out_idx = class_index_map[l] * qty_conns;
                            for( m = 0; m < qty_conns; m++ )
                            {
                                conn = intValues[in_idx+m];
                                p_ed->nodes[out_idx+m] = node_index_map[conn];
                            }
                        }
                        Py_DECREF( py_ProcConns );
                        Py_DECREF( py_ByteArray );

                        /* Load in element materials and part numbers */
                        py_ProcMats  = PySequence_ITEM( py_ElemClassMats, k );
                        py_ProcParts = PySequence_ITEM( py_ElemClassParts, k );

                        PyErr_Clear();
                        py_MatBytes = PyByteArray_FromObject( py_ProcMats );
                        py_PartBytes = PyByteArray_FromObject( py_ProcParts );

                        bytearray = PyByteArray_AS_STRING( py_MatBytes );
                        mats = (int*) bytearray;

                        bytearray = PyByteArray_AS_STRING( py_PartBytes );
                        parts = (int*) bytearray;

                        for( l = 0; l < elems_on_proc; l++ )
                        {
                            out_idx = class_index_map[l];
                            p_ed->mat[out_idx] = mats[l] - 1;
                            p_ed->part[out_idx] = parts[l] - 1;
                        }
                        Py_DECREF( py_ProcMats );
                        Py_DECREF( py_MatBytes );
                        Py_DECREF( py_ProcParts );
                        Py_DECREF( py_PartBytes );
                    }
                }

                elem_class->labels_max = -1;
                elem_class->labels_min = MAXINT;
                elem_class->labels_found = TRUE;

                /* Construct element labels */
                elem_class->labels = NEW_N( MO_class_labels, obj_qty, "Class Labels" );
                if ( elem_class->labels == NULL )
                    popup_fatal( "Unable to allocate labels on class load" );

                elem_class->labels_index = NEW_N( int, obj_qty, "Class Labels index" );
                if ( elem_class->labels_index == NULL )
                    popup_fatal( "Unable to allocate labels index on class load" );

                for( obj_id = 0; obj_id < obj_qty; obj_id++ ){
                    elem_class->labels[obj_id].local_id = obj_id;
                    elem_class->labels[obj_id].label_num = elem_labels[obj_id];
                    if ( elem_labels[obj_id] > elem_class->labels_max )
                        elem_class->labels_max = elem_labels[obj_id];
                    if ( elem_labels[obj_id] < elem_class->labels_min )
                        elem_class->labels_min = elem_labels[obj_id];
                }

                /* Sort the labels */
                qsort( elem_class->labels, obj_qty, sizeof(MO_class_labels), mili_compare_labels);

                /* Create a mapping for the 1-n label index */
                for ( obj_id = 0; obj_id < obj_qty; obj_id++ ){
                    elem_class->labels_index[elem_class->labels[obj_id].local_id] = obj_id; 
                }

                /* Construct blocking data for elem labels */
                block_range_ptr = compute_label_blocking_data( elem_labels, obj_qty, &block_qty);
                if( block_range_ptr == NULL ){
                    fprintf(stderr, "mili_reader_get_geom call compute_label_blocking_data\n");
                    return GRIZ_FAIL;
                }

                /* Construct the label blocking table of contents */
                if( elem_class->labels_found && block_qty > 0 && block_range_ptr ){
                    elem_class->label_blocking.block_qty = block_qty;
                    elem_class->label_blocking.block_total_objects = 0;
                    elem_class->label_blocking.block_min = MAXINT;
                    elem_class->label_blocking.block_max = MININT;
                    elem_class->label_blocking.block_objects = NEW_N( Label_block_data, block_qty, "Elem Class Label Blocking Objects");

                    block_index = 0;
                    for( k = 0; k < block_qty; k++ ){
                        // Update min and max labels for this block 
                        if( block_range_ptr[block_index] < elem_class->label_blocking.block_min )
                            elem_class->label_blocking.block_min = block_range_ptr[block_index];
                        if( block_range_ptr[block_index] > elem_class->label_blocking.block_max )
                            elem_class->label_blocking.block_max = block_range_ptr[block_index];

                        elem_class->label_blocking.block_objects[k].label_start = block_range_ptr[block_index++];
                        elem_class->label_blocking.block_objects[k].label_stop = block_range_ptr[block_index++];
                        elem_class->label_blocking.block_total_objects += ( elem_class->label_blocking.block_objects[k].label_stop -
                            elem_class->label_blocking.block_objects[k].label_start) + 1;
                    }
                }

                if( block_range_ptr ){
                    free(block_range_ptr);
                    block_range_ptr = NULL;
                }

                /* Element superclass-specific actions. */
                switch ( elem_class->superclass )
                {
                case M_HEX:
                    check_degen_hexs( elem_class );
                    break;
                case M_WEDGE:
                    popup_dialog( INFO_POPUP, "%s\n%s (class \"%s\")",
                                    "Checking for degenerate wedge elements",
                                    "is not implemented.", short_name );
                    break;
                case M_PYRAMID:
                    popup_dialog( INFO_POPUP, "%s\n%s (class \"%s\")",
                                    "Checking for degenerate pyramid",
                                    "elements is not implemented.",
                                    short_name );
                    break;
                case M_TET:
                    popup_dialog( INFO_POPUP, "%s\n%s (class \"%s\")",
                                    "Checking for degenerate tet elements",
                                    "is not implemented.", short_name );
                    break;
                case M_QUAD:
                    check_degen_quads( elem_class );
                    break;
                case M_TRI:
                    check_degen_tris( elem_class );
                    if ( elem_class->objects.elems->has_degen )
                        popup_dialog( INFO_POPUP, "%s\n(class \"%s\").",
                                        "Degenerate tri element(s) detected",
                                        short_name );
                    break;
                default:
                    // do nothing
                    ;
                }

                /* Update Mesh_data classes_by_sclass list for elements superclass */
                p_lh = p_md->classes_by_sclass + elem_class->superclass;
                mo_classes = (MO_class_data **) p_lh->list;
                mo_classes = (void *) RENEW_N( MO_class_data *, mo_classes, p_lh->qty, 1,
                                                "Extend classes_by_sclass array" );
                mo_classes[p_lh->qty] = elem_class;
                p_lh->qty++;
                p_lh->list = (void *) mo_classes;

            }

            /* Free elem labels */
            if( elem_labels != NULL ){
                free( elem_labels );
                elem_labels = NULL;
            }
        }

        iend = prec_timer();
        printf("Nodal class load = %ld\n", (iend-istart));

        // Need to call gen_material_data on M_MAT classes after all other classes are processed.
        MO_class_data ** mat_classes = p_md->classes_by_sclass[M_MAT].list;
        for( j = 0; j < p_md->classes_by_sclass[M_MAT].qty; j++ ){
            elem_class = mat_classes[i];
            p_matd = NEW_N( Material_data, elem_class->qty, "Material data array" );
            elem_class->objects.materials = p_matd;
            gen_material_data( elem_class, p_md );
        }

        /* Update the Mesh_data struct with element class info */
        p_md->elem_class_qty = elem_class_qty;

    } // for ( i = 0 i < mesh_qty; i++ )

    get_hex_volumes( analy->db_ident, analy );

#ifdef MILI_READER_TIMING
    end = prec_timer();
    printf("[mili_reader_get_geom] elapsed = %ldms\n", (end-start));
#endif

    /* Freeing these for now, may need them later though. May have to move to analysis struct or global. */
    free( py_MOclass_by_proc );
    free( py_nodes_by_proc );
    free( py_labels_by_proc );
    free( py_conns_by_proc );

    return OK;
}


/* TAG( get_subrecord_from_pyobject )
 *
 * Convert pyobject * to a Subrecord
 */
Bool_type
get_subrecord_from_pyobject( PyObject *py_Subrecord, Subrecord *p_subrec )
{
    int i;
    PyObject * py_List,
             * py_SvarName;

    p_subrec->name = get_pyobject_attribute_as_string( py_Subrecord, "name" );
    p_subrec->class_name = get_pyobject_attribute_as_string( py_Subrecord, "class_name" );
    p_subrec->organization = get_pyobject_attribute_as_int( py_Subrecord, "organization" );
    p_subrec->qty_svars = get_pyobject_attribute_as_int( py_Subrecord, "qty_svars" );
    p_subrec->superclass = get_pyobject_attribute_as_int( py_Subrecord, "superclass" );
    p_subrec->qty_blocks = 0;
    p_subrec->qty_objects = 0;

    p_subrec->svar_names = NEW_N( char *, p_subrec->qty_svars, "Subrecord svar name ptrs" );
    py_List = get_pyobject_attribute( py_Subrecord, "svar_names" );
    for( i = 0; i < p_subrec->qty_svars; i++ ){
        py_SvarName = PyList_GET_ITEM( py_List, i );
        p_subrec->svar_names[i] = pyobject_as_string( py_SvarName );
        if ( p_subrec->svar_names[i] == NULL )
            return GRIZ_FAIL;
    }

    return OK;
}


/* TAG( get_state_variable_from_pyobject )
 *
 * Convert pyobject * to a State_variable
 */
Bool_type
get_state_variable_from_pyobject( PyObject *py_StateVariable, State_variable *p_sv )
{
    int i;
    PyObject * py_List,
             * py_SvarName;

    p_sv->short_name = get_pyobject_attribute_as_string( py_StateVariable, "name" );
    p_sv->long_name = get_pyobject_attribute_as_string( py_StateVariable, "title" );
    p_sv->num_type = get_pyobject_attribute_as_int( py_StateVariable, "data_type" );
    p_sv->agg_type = get_pyobject_attribute_as_int( py_StateVariable, "agg_type" );
    p_sv->rank = 1;
    p_sv->vec_size = 1;
    p_sv->dims = NULL;

    if( p_sv->agg_type != SCALAR )
    {
        if( p_sv->agg_type != VECTOR )
        {
            PyObject * py_Dims = get_pyobject_attribute( py_StateVariable, "dims" );
            p_sv->rank = get_pyobject_attribute_as_int( py_StateVariable, "order" );
            p_sv->dims = NEW_N( int, p_sv->rank, "State_variable arr dims" );
            if ( p_sv->rank > 0 && p_sv->dims == NULL ){
                return GRIZ_FAIL;
            }
            integer_pointer_from_pytuple( py_Dims, p_sv->rank, p_sv->dims );
        }

        if ( p_sv->agg_type == VEC_ARRAY || p_sv->agg_type == VECTOR )
        {
            /* Copy vector size */
            p_sv->vec_size = get_pyobject_attribute_as_int( py_StateVariable, "list_size" );

            PyObject * py_CompNames = get_pyobject_attribute( py_StateVariable, "comp_svar_names" );
            PyObject * py_CompTitles = get_pyobject_attribute( py_StateVariable, "comp_svar_titles" );
            if ( py_CompNames == NULL || py_CompTitles == NULL )
            {
                fprintf( stderr, "comp_svar_names or comp_svar_titles failed for %s\n", p_sv->short_name );
                return GRIZ_FAIL;
            }

            /* Copy component names and titles */
            p_sv->components = NEW_N( char*, p_sv->vec_size, "Vector Array component names" );
            p_sv->component_titles = NEW_N( char*, p_sv->vec_size, "Vector Array component titles" );

            if( p_sv->vec_size > 0 && (p_sv->components == NULL || p_sv->component_titles == NULL ))
            {
                return GRIZ_FAIL;
            }

            PyObject *py_Name, *py_Title;
            for( i = 0; i < p_sv->vec_size; i++ )
            {
                py_Name = PySequence_ITEM( py_CompNames, i );
                py_Title = PySequence_ITEM( py_CompNames, i );

                p_sv->components[i] = pyobject_as_string( py_Name );
                p_sv->component_titles[i] = pyobject_as_string( py_Title );
            }
        }
    }

    return OK;
}


/* TAG( create_st_variable_from_pyobject )
 *
 * Convert a PyObject * to mili reader StateVariable object into a Mili State_variable object
 * and add to the hash table by name.
 */
Bool_type
create_st_variable_from_pyobject( Hash_table *p_sv_ht, char *p_name,
                                  State_variable **pp_svar, PyObject *py_StateVariable, int* num_added )
{
    int rval;
    Hash_action op;
    Htable_entry *p_hte;
    State_variable *p_sv;
    int i;

    *num_added = 0;

    /* Only create svar if not already in table. Exception for sand variable. */
    if( strncmp( p_name, "sand", 4 ) == 0 )
        op = ENTER_ALWAYS;
    else
        op = ENTER_UNIQUE;

    rval = htable_search( p_sv_ht, p_name, op, &p_hte );

    /* If this is a new entry... */
    if ( rval == OK )
    {
        /* Create the State_variable and store it in the hash table */
        p_sv = NEW( State_variable, "New state var" );
        rval = get_state_variable_from_pyobject( py_StateVariable, p_sv );
        if ( rval != OK )
        {
            fprintf(stderr, "create_st_variable_from_pyobject call get_state_variable_from_pyobject(%s)\n", p_name);
            return GRIZ_FAIL;
        }

        p_hte->data = (void*) p_sv;

        if( pp_svar )
            *pp_svar = p_sv;
        
        *num_added = 1;
    }
    else if ( pp_svar )
    {
        rval = htable_search( p_sv_ht, p_name, FIND_ENTRY, &p_hte );
        if ( rval == OK )
            *pp_svar = (State_variable*) p_hte->data;
    }

    return OK;
}


/* TAG( mili_reader_load_element_sets )
 *
 * Load element set data from mili reader.
 */
int
mili_reader_load_element_sets(Analysis *analy, PyObject ** py_ESList)
{
    int i, j;
    int rval;
    int cnt;
    char * es_name;
    char **es_names;
    int ipt, ipt_count;
    ElementSet * element_set;
    Htable_entry *p_hte;
    int label_length = strlen("es_");

    PyObject * py_ESProcList;
    PyObject * py_ES;
    PyObject * py_DictItems;
    PyObject * py_Tuple,
             * py_Key,
             * py_Value,
             * py_Ipt;

    analy->Element_sets = NULL;
    analy->Element_set_names = NULL;
    analy->es_cnt = 0;

    if ( py_ESList == NULL )
    {
        fprintf( stderr, "mili_reader_load_element_sets: py_ESList is NULL\n" );
        return GRIZ_FAIL;
    }

    for( i = 0; i < analy->proc_count; i++ )
    {
        py_ESProcList = py_ESList[i];
        cnt = PyDict_Size( py_ESProcList );
        
        /* If there are element sets on this processor */
        if ( cnt > 0 )
        {
            /* Initialize hash table if not already done */
            if ( analy->Element_sets == NULL )
                analy->Element_sets = htable_create( 151 );

            /* Add Element_set entries to hash table */
            py_DictItems = PyDict_Items( py_ESProcList );

            for( j = 0 ; j < cnt; j++ )
            {
                /* Get element set name */
                py_Tuple = PySequence_ITEM( py_DictItems, j );
                py_Key = PyTuple_GetItem( py_Tuple, 0 );
                es_name = pyobject_as_string( py_Key );

                /* Check if we already loaded this element set on another processor */
                rval = htable_search( analy->Element_sets, es_name, ENTER_UNIQUE, &p_hte );
                
                /* If new element set... */
                if ( rval == OK )
                {
                    /* Create new element set entry */
                    element_set = NEW( ElementSet, "New Element Set entry" );

                    /* Load integration point data */
                    py_Value = PyTuple_GetItem( py_Tuple, 1 );
                    ipt_count = PyTuple_Size( py_Value );
                    element_set->integration_points = NEW_N( int, ipt_count, "Element set ipt list" );
                    integer_pointer_from_pytuple( py_Value, ipt_count, element_set->integration_points );

                    /* Initialize other element set data */
                    element_set->size = ipt_count;
                    element_set->material_number = atoi(es_name + label_length);
                    set_es_middle_index( element_set );
                    element_set->current_index = element_set->middle_index;
                    element_set->tempIndex = -1;

                    /* Update list of element set names */
                    if( analy->Element_set_names == NULL ){
                        analy->Element_set_names = NEW_N( char*, 1, "New element set names entry");
                    }
                    else{
                        es_names = analy->Element_set_names;
                        es_names = RENEW_N( char*, es_names, analy->es_cnt, 1, "Expand element set names list");
                        /* Reassign in case realloc moved memory location */
                        analy->Element_set_names = es_names;
                    }
                    analy->Element_set_names[analy->es_cnt] = es_name;

                    analy->es_cnt += 1;

                    p_hte->data = (void*) element_set;
                }
            }
        }
    }

    return OK;
}


/* TAG( mili_reader_get_st_descriptors )
 *
 * Query mili reader and store information about the available
 * state record formats and referenced state variables.
 */
extern int
mili_reader_get_st_descriptors( Analysis *analy, int dbid )
{
    int i, j, k, l;
    int rval, cnt;
    int total_subrec_qty;
    int total_svar_qty;
    int srec_qty, svar_qty, subrec_qty, mesh_node_qty;
    int qty_comp_svars;
    int agg_type;
    int srec_id;
    int mesh_id;
    int class_size;
    int obj_qty, block_qty;
    int gid_block_qty;
    int block_start, block_end;
    int gidx;
    int * elem_gids, * elem_blocks;
    int elem_class_index;
    int subrec_index = 0;
    int * block_range_ptr;
    int * mo_blocks;
    char * name;
    char * svar_name, * comp_svar_name;
    char * class_name;
    char **svar_names;
    Subrecord * p_subr;
    Subrecord **p_htable_subrecords;
    State_rec_obj* p_sro;
    Subrec_obj* p_subrecs;
    Htable_entry *p_hte, *p_hte2;
    MO_class_data *p_mocd;
    Bool_type nodal, particle;
    char *node_work_array=NULL;
    State_variable *p_svar;
    Mesh_data *p_mesh;

    Hash_table * p_subrecord_ht;
    Hash_table * p_subrecord_gids_ht;
    Hash_table * p_sv_ht, * p_primal_ht;
     
    PyObject * py_CallData;
    PyObject * py_ES, * py_SubrecordList;
    PyObject ** py_ESList;
    PyObject * py_Subrecords,
             * py_Subrecord,
             * py_SrecQty;
    PyObject * py_SubrecordStateVariables,
             * py_StateVariable,
             * py_ComponentSvars,
             * py_ComponentSvar;

#ifdef MILI_READER_TIMING
    long start, end;
    start = prec_timer();
#endif

    analy->num_bad_subrecs = 0;
    analy->bad_subrecs = NULL;
    analy->old_shell_stresses = FALSE;

    analy->py_PreloadedResult = NULL;

    /* Single call to mili reader to get all the data we need */
    py_CallData = call_mili_reader_function_noargs( analy->py_MiliDB, GET_ST_DESCRIPTORS_CALL );
    if( py_CallData == NULL ){
        fprintf(stderr, MR_DEFAULT_FAILED, "mili_reader_get_st_descriptors", GET_ST_DESCRIPTORS_CALL);
        return GRIZ_FAIL;
    }
    py_ES = get_pyobject_attribute( py_CallData, "element_sets" );
    py_SubrecordList = get_pyobject_attribute( py_CallData, "subrecords" );

    /* Get data for each processor */
    int proc_count = analy->proc_count;
    int * subrecord_count_per_proc = (int*) malloc( proc_count * sizeof(int));
    PyObject **py_subrecord_list_by_proc = (PyObject**) malloc(proc_count * sizeof(PyObject*));
    py_ESList = NEW_N( PyObject*, proc_count, "PyObject Element set List" );

    /* Store data for each processor */
    for( i = 0; i < proc_count; i++ ){
        py_ESList[i] = PySequence_ITEM( py_ES, i );
        py_subrecord_list_by_proc[i] = PySequence_ITEM( py_SubrecordList, i );
        subrecord_count_per_proc[i] = PySequence_Length( py_subrecord_list_by_proc[i] );
    }

    /* Get Mesh id and number of state record formats */
    mesh_id = 0;
    srec_qty = get_pyobject_attribute_as_int( analy->py_MiliDB, "srec_fmt_qty" );

    /* Load in element sets */
    rval = mili_reader_load_element_sets( analy, py_ESList );
    if( rval != OK ){
        fprintf( stderr, "mili_reader_get_st_descriptors call mili_reader_load_element_sets\n" );
        return GRIZ_FAIL;
    }

    /* Create hash table to store subrecord information while processing */
    p_subrecord_ht = htable_create( 151 );
    p_subrecord_gids_ht = htable_create( 151 );
    /* Create hash table to store state variables */
    p_sv_ht = htable_create ( 1009 );

    /* Get all subrecords from each processor and combine in hash table */
    total_svar_qty = 0;
    /* For each processor... */
    for( i = 0; i < proc_count; i++ ){
        py_Subrecords = py_subrecord_list_by_proc[i];
        subrec_qty = subrecord_count_per_proc[i];
        /* For each subrecord on the processor... */
        for( j = 0; j < subrec_qty; j++ ){
            py_Subrecord = PySequence_ITEM( py_Subrecords, j );
            name = get_pyobject_attribute_as_string( py_Subrecord, "name" );

            /* Check if subrecord already exists in hash table */
            rval = htable_search( p_subrecord_ht, name, ENTER_MERGE, &p_hte );

            /* New subrecord */
            if ( rval == OK ){
                /* Create hash table entry for subrecord */
                p_subr = NEW( Subrecord, "Subrecord obj");
                if( get_subrecord_from_pyobject( py_Subrecord, p_subr ) != OK ){    
                    fprintf( stderr, "[%d] mili_reader_get_st_descriptors call get_subrecord_from_pyobject(%s)\n", i, name);
                    return GRIZ_FAIL;
                }
                p_hte->data = (void *) p_subr;
                
                /* Get state variables in this subrecord. */
                py_SubrecordStateVariables = get_pyobject_attribute( py_Subrecord, "svars" );
                svar_qty = PySequence_Length( py_SubrecordStateVariables );

                /* Create State_variable objects and add entries to p_sv_ht hashtable */
                for( k = 0; k < svar_qty; k++ )
                {
                    py_StateVariable = PySequence_ITEM( py_SubrecordStateVariables, k );
                    svar_name = get_pyobject_attribute_as_string( py_StateVariable, "name" );

                    /* Create state variable from pyobject and add to p_sv_ht */
                    rval = create_st_variable_from_pyobject( p_sv_ht, svar_name, &p_svar, py_StateVariable, &cnt );
                    total_svar_qty += cnt;
                    if( rval != OK ){
                        fprintf( stderr, "[%d] mili_reader_get_st_descriptors call create_st_variable_from_pyobject(%s)\n", i, svar_name);
                        return GRIZ_FAIL;
                    }

                    if( p_svar->agg_type == VECTOR || p_svar->agg_type == VEC_ARRAY ){
                        py_ComponentSvars = get_pyobject_attribute( py_StateVariable, "svars" );
                        if ( py_ComponentSvars != NULL )
                        {
                            qty_comp_svars = PySequence_Length( py_ComponentSvars ); 
                            for ( l = 0; l < qty_comp_svars; l++ )
                            {
                                py_ComponentSvar = PySequence_ITEM( py_ComponentSvars, l );
                                rval = create_st_variable_from_pyobject( p_sv_ht, p_svar->components[l], NULL, py_ComponentSvar, &cnt );
                                total_svar_qty += cnt;
                                if( rval != OK ){
                                    fprintf( stderr, "[%d] mili_reader_get_st_descriptors call create_st_variable_from_pyobject(%s)\n", i, p_svar->components[l]);
                                    return GRIZ_FAIL;
                                }
                            }
                        }
                    }
                }
            }

            /* Handle subrecord MO_blocks */
            p_subr = (Subrecord*) p_hte->data;
            class_name = get_pyobject_attribute_as_string( py_Subrecord, "class_name" );
            rval = htable_search( analy->mesh_table[mesh_id].class_table, class_name, FIND_ENTRY, &p_hte );
            if( rval != OK ){
                fprintf( stderr, "mili_reader_get_st_descriptors call htable_search(%s)\n", class_name);
                return GRIZ_FAIL;
            }

            p_mocd = (MO_class_data*) p_hte->data;
            elem_class_index = p_mocd->elem_class_index;

            // Need to handle qty_objects, qty_blocks, mo_blocks
            PyObject* py_OrdinalBlocks = get_pyobject_attribute( py_Subrecord, "ordinal_blocks" );
            block_qty = PySequence_Length( py_OrdinalBlocks ) / 2;

            if( block_qty > 0 ){
                elem_blocks = NEW_N(int, block_qty * 2, "Subrecord ordinal blocks");
                integer_pointer_from_pyobject( py_OrdinalBlocks, block_qty * 2, elem_blocks );

                if( p_mocd->superclass == G_NODE ){
                    /* Special case for node subrecords because the same node can be on multiple processors */
                    /* For now assume all nodes are in the same subrecords, just use blocks from MO_class */
                    if( p_subr->mo_blocks == NULL ){
                        p_subr->mo_blocks = NEW_N( int, 2, "Subrecord MO blocks");
                        p_subr->mo_blocks[0] = 1;
                        p_subr->mo_blocks[1] = p_mocd->qty;
                        p_subr->qty_objects = p_mocd->qty;
                        p_subr->qty_blocks = 1;
                    }
                    else{
                        // FIXME
                        /* Skip for now */
                        continue;
                    }
                }
                else{
                    /* Convert local processor mo blocks to list of global ids for elem class */
                    for( k = 0; k < block_qty; k++ ){
                        block_start = elem_blocks[k*2];
                        block_end  = elem_blocks[k*2+1];
                        obj_qty = (block_end - block_start);
                        /* Update object quantity in subrecord object. */
                        p_subr->qty_objects += obj_qty;

                        elem_gids = NEW_N(int, obj_qty, "Subrecord Element Gids");
                        gidx = 0;
                        for( l = block_start; l < block_end; l++ ){
                            elem_gids[gidx++] = MESH_P(analy)->index_map->elem_map[elem_class_index][i][l];
                        }

                        block_range_ptr = compute_label_blocking_data( elem_gids, obj_qty, &gid_block_qty );
                        if( block_range_ptr == NULL ){
                            fprintf( stderr,
                                    "[%d] mili_reader_get_st_descriptors call compute_label_blocking for class '%s'\n",
                                    i, class_name );
                            return GRIZ_FAIL;
                        }

                        /* Create/Expand subrecord mo blocks and add new element blocks */
                        if( p_subr->mo_blocks == NULL ){
                            p_subr->mo_blocks = NEW_N( int, gid_block_qty*2, "Subrecord MO blocks");
                        }
                        else{
                            mo_blocks = p_subr->mo_blocks;
                            mo_blocks = RENEW_N( int, mo_blocks, p_subr->qty_blocks*2, gid_block_qty*2, "Expand Subrecord MO block");
                            p_subr->mo_blocks = mo_blocks;
                        }
                        gidx = p_subr->qty_blocks * 2;
                        for( l = 0; l < gid_block_qty; l++ ){
                            /* Plus 1 to convert from 0 based indexing to 1 based from mili reader. */
                            p_subr->mo_blocks[gidx++] = block_range_ptr[l*2] + 1;
                            p_subr->mo_blocks[gidx++] = block_range_ptr[l*2+1] + 1;
                        }
                        p_subr->qty_blocks += gid_block_qty;

                    }

                    free( elem_blocks );
                }
            }
        }
    }

    /* Initialize state variable hash table. */
    if ( total_svar_qty > 0 ){
        if ( total_svar_qty < 100 )
            p_primal_ht = htable_create( 151 );
        else
            p_primal_ht = htable_create( 1009 );
    }
    else
        return OK;

    /* Get merged subrecords */
    htable_get_data( p_subrecord_ht, (void***) &p_htable_subrecords, &total_subrec_qty);

    /* Allocate array of state record structs. */
    p_sro = NEW_N( State_rec_obj, srec_qty, "Srec tree" );

    srec_id = 0;
    /* Allocate array of Subrec_obj's. */
    p_subrecs = NEW_N( Subrec_obj, total_subrec_qty, "Srec subrec branch");
    p_sro[srec_id].subrecs = p_subrecs;
    p_sro[srec_id].qty = total_subrec_qty;

    /* Init nodal position and velocity subrec indices to invalid values. */
    p_sro[srec_id].node_pos_subrec = -1;
    p_sro[srec_id].node_vel_subrec = -1;
    p_sro[srec_id].particle_pos_subrec = -1;

    p_mesh = analy->mesh_table + mesh_id;
    p_sro[srec_id].mesh_id = mesh_id;

    /* Allocate a temporary working array for subrec node list creation. */
    mesh_node_qty = p_mesh->node_geom->qty;
    node_work_array = NEW_N( char, mesh_node_qty, "Temp node array" );

    /* Populate Griz subrecord and state variable data */
    for( i = 0; i < total_subrec_qty; i++ ){
        p_subr = p_htable_subrecords[i];
        p_subrecs[i].subrec = *p_subr;

        /* Look up the element class in the subrecord */
        htable_search( p_mesh->class_table, p_subr->class_name, FIND_ENTRY, &p_hte );
        p_mocd = (MO_class_data*) p_hte->data;
        p_subrecs[i].p_object_class = p_mocd;

        /* Create list of nodes referenced by objects bound to subrecord */
        create_subrec_node_list( node_work_array, mesh_node_qty, p_subrecs + i );

        /* Create ident array if indexing is required */
        class_size = p_mocd->qty;

        p_subrecs[i].object_ids = NULL;
        if( p_subr->qty_objects != class_size || (class_size > 1 && (p_subr->qty_blocks != 1 || p_subr->mo_blocks[0] != 1 )) )
        {
            p_subrecs[i].object_ids = NEW_N( int, p_subr->qty_objects, "Subrec ident map array" );
            blocks_to_list( p_subr->qty_blocks, p_subr->mo_blocks, p_subrecs[i].object_ids, TRUE );
        }
        else
        {
            p_subrecs[i].object_ids = NULL;
        }

        /* M_NODE class "node" is special - need it for node positions and velocities. */
        nodal = FALSE;
        if ( strcmp( p_subr->class_name, "node" ) == 0 )
            nodal = TRUE;

        /* Loop over svars and create state variable and primal result table entries. */
        svar_names = p_subr->svar_names;
        for( j = 0; j < p_subr->qty_svars; j++ )
        {
            /* Look up State_variable object so we can create the primal result */
            rval = htable_search( p_sv_ht, svar_names[j], FIND_ENTRY, &p_hte );
            if ( rval != OK ){
                fprintf( stderr, "mili_reader_get_st_descriptors call htable_search( p_sv_ht, %s )\n", svar_names[i]);
                return GRIZ_FAIL;
            }
            p_svar = (State_variable*) p_hte->data;
            
            // Create primal result
            create_primal_result( p_mesh, srec_id, i, p_subrecs+i, p_primal_ht, srec_qty,
                                  svar_names[j], p_sv_ht, analy );

            // Check for old shell stresses
            if(!strcmp(svar_names[j], "stress_in") || !strcmp(svar_names[j], "stress_mid") || !strcmp(svar_names[j], "stress_out"))
                analy->old_shell_stresses = TRUE;
            
            if ( nodal )
            {
                if ( strcmp( svar_names[j], "nodpos" ) == 0 )
                {
                    if ( p_sro[srec_id].node_pos_subrec != -1 )
                        popup_dialog( WARNING_POPUP, "Multiple \"node\" position subrecs." );

                    p_sro[srec_id].node_pos_subrec = i;
                    analy->stateDB = FALSE;
                    if( mesh_node_qty == p_subr->qty_objects )
                    {
                        /* This is a state database and not a time history database. */
                        analy->stateDB = TRUE;
                    }

                    /* Note if data is double precision. */
                    if ( p_svar->num_type == M_FLOAT8 )
                        p_mesh->double_precision_nodpos = TRUE;
                }
                else if ( strcmp( svar_names[j], "nodvel" ) == 0 )
                {
                    if ( p_sro[srec_id].node_vel_subrec != -1 )
                        popup_dialog( WARNING_POPUP, "Multiple \"node\" velocity subrecs." );
                    p_sro[srec_id].node_vel_subrec = i;
                }
            }
        }
    }

    free( node_work_array );

    /* Return subrecord tree, state variable hash table and primal result hash table. */
    analy->srec_tree = p_sro;
    analy->qty_srec_fmts = srec_qty;
    analy->st_var_table = p_sv_ht;
    analy->primal_results = p_primal_ht;


#ifdef MILI_READER_TIMING
    end = prec_timer();
    printf("[mili_reader_get_st_descriptors] elapsed = %ldms\n", (end-start));
#endif
    return OK;
}


/* TAG( combine_nodpos )
 *
 * Combine the nodal positions from multiple processors
 */
extern int 
combine_nodpos( Analysis *analy, int state_no, void * out_buffer ){
    int i, j;
    int out_idx;
    int dims;
    int node_qty;
    char *nodpos = "nodpos";
    char *data = "data";
    char *node = "node";
    PyObject * py_ProcNodes;
    PyObject * py_Nodes;
    PyObject * py_Coords;

    float* out_buf = (float*) out_buffer;

    dims = analy->dimension;

    for( i = 0; i < analy->proc_count; i++ ){
        py_ProcNodes = PySequence_ITEM( analy->py_Nodpos, i);

        py_Nodes = PyDict_GetItemString( py_ProcNodes, nodpos );
        py_Nodes = PyDict_GetItemString( py_Nodes, data );
        py_Nodes = PyDict_GetItemString( py_Nodes, node );

        /* Single state, so remove state array */
        py_Nodes = PySequence_ITEM( py_Nodes, state_no );
        node_qty = MESH_P(analy)->index_map->node_count[i];

        /* Get coordinates for each element */
        for( j = 0; j < node_qty; j++ ){
            py_Coords = PySequence_ITEM( py_Nodes, j );
            out_idx = MESH_P(analy)->index_map->node_map[i][j] * dims;
            float_pointer_from_pyobject( py_Coords, dims, out_buf+out_idx);
        }
    }
    return OK;
}


/* TAG( mili_reader_get_state )
 *
 * Move to a particular state in the Mili database and update
 * nodal positions for the mesh.
 */
extern int
mili_reader_get_state( Analysis *analy, int state_no, State2 *p_st, State2 **pp_new_st, int *state_qty )
{
    int i, j;
    int rval;
    int st_qty;
    int dims;
    int srec_id;
    int mesh_id;
    int dbid;
    float st_time;

    int subrec_size;
    int ec_index;
    double *p_double;
    float *p_single;
    void *input_buf, *ibuf;
    int *object_ids;
    char *primal, *svar;
    float *p_float;

    State_variable sv;
    State_rec_obj *p_sro;
    Subrec_obj *p_subrecs, *p_subrec;
    Mesh_data *p_md;

#ifdef MILI_READER_TIMING
    long tstart, tend;
    tstart = prec_timer();
#endif

    dbid = analy->db_ident;

    /* Query the number of states in the database */
    st_qty = analy->state_count;

    /* Get the number of dimensions */
    dims = analy->dimension;

    /* Pass back number of states */
    if ( state_qty != NULL ){
        *state_qty = st_qty;
    }
    
    if ( st_qty == 0 || analy->qty_srec_fmts == 0 )
    {
        if ( p_st == NULL )
            p_st = mk_state2( analy, NULL, dims, 0, st_qty, p_st );
        
        /* No states, so use node positions from geometry definition. */
        p_st->nodes = analy->mesh_table->node_geom->objects;
        p_st->position_constant = TRUE;

        *pp_new_st = p_st;
        return OK;
    }

    /* Bounds check requested state */
    if ( state_no < 0 || state_no >= st_qty )
    {
        popup_dialog( WARNING_POPUP, "Get-state request for nonexistent state." );
        *pp_new_st = p_st;
        return GRIZ_FAIL;
    }
    
    /* Should only be a single state record format */
    srec_id = analy->state_srec_fmt_ids[state_no];

    /* Get subrecord tree, subrecords, and mesh data */
    p_sro = analy->srec_tree + srec_id;
    p_subrecs = p_sro->subrecs;
    mesh_id = p_subrecs[0].p_object_class->mesh_id;
    p_md = analy->mesh_table + mesh_id;

    /* Update or create State2 struct. */
    p_st = mk_state2( analy, p_sro, dims, srec_id, st_qty, p_st );
    p_st->state_no = state_no;

    /* Store the state time */
    p_st->time = analy->state_times[state_no];

    /* Read node position array if it exists, re-ordering if necessary. */
    if ( analy->stateDB && analy->load_nodpos)
    {
        //rval = combine_nodpos( analy, state_no, p_st->nodes.nodes );
        primal = "nodpos";
        rval = mili_reader_get_results( analy->db_ident, state_no+1, p_sro->node_pos_subrec, 1, &primal, p_st->nodes.nodes );

        /* If unable to get first state, then no state files exist */
        if ( rval != OK )
        {
            st_qty = 0;
            /* Pass back the current quantity of states in the db. */
            if ( state_qty != NULL )
                *state_qty = st_qty;

            p_st = mk_state2( analy, NULL, dims, 0, st_qty, NULL );

            /* No states, so use node positions from geometry definition. */
            p_st->nodes = analy->mesh_table->node_geom->objects;

            p_st->position_constant = TRUE;

            *pp_new_st = p_st;

            return OK;
        }
    }

    /* Read in sand flags */
    if( analy->load_sand )
    {
        for ( i = 0; i < p_sro->qty; i++ )
        {
            p_subrec = p_subrecs + i;
            if( p_subrec->sand )
            {
                primal = "sand";
                ec_index = p_subrec->p_object_class->elem_class_index;
                object_ids = p_subrec->object_ids;

                if( object_ids == NULL )
                {
                    rval = mili_reader_get_results( dbid, state_no+1, i, 1, &primal, (void*) p_st->elem_class_sand[ec_index] );
                    if( rval != 0 )
                    {
                        return GRIZ_FAIL;
                    }
                }
                else
                {
                    subrec_size = p_subrec->subrec.qty_objects;
                    input_buf = get_st_input_buffer( analy, subrec_size, p_md->double_precision_sand, &ibuf );
                    rval = mili_reader_get_results( dbid, state_no+1, i, 1, &primal, input_buf);
                    if( rval != 0 )
                    {
                        return GRIZ_FAIL;
                    }
                    reorder_float_array( subrec_size, object_ids, 1, input_buf, p_st->elem_class_sand[ec_index]);

                    if( ibuf != NULL )
                        free( ibuf );
                }
            }
        }
    }

    /* Handle spheral */
    if( p_st->sph_present ){
        primal     = "sph_itype";
        p_subrec   = p_subrecs + p_st->sph_srec_id;
        object_ids = p_subrec->object_ids;
        subrec_size = p_subrec->subrec.qty_objects;
        input_buf = get_st_input_buffer( analy, subrec_size, TRUE, &ibuf );
        rval = mili_reader_get_results( analy->db_ident, state_no+1, p_st->sph_srec_id, 1, &primal, input_buf );

        p_float = (float*) input_buf;
        for( i = 0; i < subrec_size; i++ ){
            p_st->sph_class_itype[i] = (int) p_float[i];
        }

        if( ibuf != NULL )
            free( ibuf );
    }

    /* If nodal positions weren't part of state data, get from geometry. */
    if ( !analy->stateDB )
    {
        p_st->nodes = p_md->node_geom->objects;
        p_st->position_constant = TRUE;
    }

    /* Return new state */
    *pp_new_st = p_st;

#ifdef MILI_READER_TIMING
    tend = prec_timer();
    //printf("[mili_reader_get_state] elapsed = %ldms\n", (tend-tstart));
#endif
    return OK;
}


/* TAG( mili_reader_get_metadata )
 *
 * Get the metadata from the A file
 */
extern int
mili_reader_get_metadata( Analysis * analy )
{
    mili_reader_get_string( analy->py_MiliParameters, "lib version", analy->mili_version );
    mili_reader_get_string( analy->py_MiliParameters, "host name", analy->mili_host );
    mili_reader_get_string( analy->py_MiliParameters, "arch name", analy->mili_arch );
    mili_reader_get_string( analy->py_MiliParameters, "date", analy->mili_timestamp );
    mili_reader_get_string( analy->py_MiliParameters, "xmilics version", analy->xmilics_version );

    return OK;
}


/* TAG( mili_reader_load_parameters )
 *
 * Get the parameters dictionary of the 0th processor
 */
extern int
mili_reader_load_parameters( Analysis * analy )
{
    PyObject * py_Parameters = call_mili_reader_function_noargs( analy->py_MiliDB, PARAMETER_GETTER );
    if( py_Parameters == NULL )
        return NOT_OK;
    analy->py_MiliParameters = py_Parameters;
    return OK;
}


/* TAG( mili_reader_get_title )
 *
 * Get title from parameters dictionary.
 */
extern int
mili_reader_get_title( Analysis * analy, char * str ){
    mili_reader_get_string( analy->py_MiliParameters, "title", str );
    return OK;
}


/* TAG( mili_reader_reload_states )
 *
 * Reload the statemaps for the given data base.
 */
extern int
mili_reader_reload_states()
{
    Analysis * analy = get_analy_ptr();
    call_mili_reader_function_noargs( analy->py_MiliDB, RELOAD_FUNCTION );
    return OK;
}


/* TAG( mili_reader_db_close )
 *
 * Close the mili reader database and stop Python
 */
extern int
mili_reader_db_close( Analysis * analy )
{
#ifdef MILI_READER_TIMING
    long start, end;
    start = prec_timer();
#endif

    if( analy->py_MiliDB || analy->py_MiliReaderModule ){
        Py_FinalizeEx(); 
    }

#ifdef MILI_READER_TIMING
    end = prec_timer();
    printf("[mili_reader_db_close] elapsed = %ldms\n", (end-start));
#endif
    return OK;
}


/* TAG( mili_reader_load_state_data )
 *
 * Load the number of states and the times of each state from the Mili reader
 * Also Loads in nodal positions and sand flags for all elements/nodes.
 */
extern int
mili_reader_load_state_data( Analysis * analy )
{
#ifdef MILI_READER_TIMING
    long start, end;
    start = prec_timer();
#endif

    int i;
    int qty_states;
    float st_time;
    int states_added;
    int start_state;
    float * state_times;
    int * state_srecs;
    PyObject * py_SMap;
    PyObject * py_StateData;

    /* Get state maps from the mili reader */
    PyObject * py_StateMaps = call_mili_reader_function_noargs( analy->py_MiliDB, STATE_MAPS_GETTER );
    if( py_StateMaps == NULL ){
        fprintf( stderr, MR_DEFAULT_FAILED, "mili_reader_load_state_data", STATE_MAPS_GETTER );
        return NOT_OK;
    }

    qty_states = PySequence_Length( py_StateMaps );

    if( analy->state_times == NULL ){
        /* First time, so load up all states available */
        /* Allocate memory to store all states times in the Analysis struct */
        analy->state_times = NEW_N( float, qty_states, "State time array" );
        analy->state_srec_fmt_ids = NEW_N( int, qty_states, "State srec fmt id array");
        start_state = 0;
    }
    else if( analy->state_count != qty_states ){
        /* Called the reload command, and there were new states to we need to update arrays */
        start_state = analy->state_count;
        states_added = qty_states - analy->state_count;

        state_times = analy->state_times;
        state_times = RENEW_N( float, state_times, analy->state_count, states_added, "Extend state_times array" );
        analy->state_times = state_times;

        state_srecs = analy->state_srec_fmt_ids;
        state_srecs = RENEW_N( int, state_srecs, analy->state_count, states_added, "Extend state_srecs array" );
        analy->state_srec_fmt_ids = state_srecs;
    }
    else{
        /* Same number of states, no need to update */
        return OK;
    }

    analy->state_count = qty_states;

    for( i = start_state; i < qty_states; i++ ){
        py_SMap = PySequence_ITEM( py_StateMaps, i );
        st_time = get_pyobject_attribute_as_double( py_SMap, "time" );
        analy->state_times[i] = st_time;

        /* There should only be a single state record format */
        analy->state_srec_fmt_ids[i] = 0;
    }


#ifdef MILI_READER_TIMING
    end = prec_timer();
    printf("[mili_reader_load_state_times] elapsed = %ldms\n", (end-start));
#endif
    return OK;
}


/* TAG( mili_reader_read_string )
 *
 * Get a string parameter from the mili reader.
 */
int
mili_reader_read_string( int dbid, char* key, char* value)
{
    int status;
    Analysis * p_analysis = get_analy_ptr();
    status = mili_reader_get_string( p_analysis->py_MiliParameters, key, value);
    return status;
}


/* TAG( mili_reader_search_param_wildcard )
 *
 * Search for a parameter in the mili database
 * There are many arguments that are unused so that this matches the function
 * mc_ti_htable_wildcard_search
 */
int
mili_reader_search_param_wildcard( int fid, int list_len, Bool_type allow_duplicates,
                                   char * key1, char * key2, char * key3, char ** return_list )
{
    int i;
    int num_entries = 0;
    char * str;
    PyObject * py_Matches;
    PyObject * py_FuncName;
    PyObject * py_Key;
    PyObject * py_Param;
    Analysis * p_analysis = get_analy_ptr();

    if( strcmp(key1, "*") == 0 ){
        num_entries = PyDict_Size( p_analysis->py_MiliParameters );
        return num_entries;
    }

    /* Retrieve state data from mili reader */
    py_FuncName = string_to_pyobject( PARAM_WILDCARD_SEARCH );
    py_Key = string_to_pyobject( key1 );
    py_Matches = PyObject_CallMethodObjArgs( p_analysis->py_MiliDB, py_FuncName, py_Key, NULL );
    if( py_Matches == NULL ){
        num_entries = 0;
        return num_entries;
    }

    num_entries = PySequence_Length( py_Matches );
    if( return_list == NULL ){
        return num_entries;
    }

    if( num_entries > 0 ){
        for( i = 0; i < num_entries; i++ ){
            py_Param = PySequence_ITEM( py_Matches, i );
            str = pyobject_as_string( py_Param );
            griz_str_dup( &return_list[i], str );
        }
    }

    return num_entries;
}


/* TAG( mili_reader_read_ti_array )
 *
 * Get an array of floats from the mili reader.
 */
int
mili_reader_read_ti_array( int dbid, char* key, void** array, int * size)
{
    PyObject * py_Array;
    PyObject * py_Char;
    PyObject * py_Key;
    int status;
    int num_values;
    Analysis * p_analysis = get_analy_ptr();

    /* Lookup parameter in mili and convert to string */
    py_Key = string_to_pyobject(key);

    if( PyDict_Contains( p_analysis->py_MiliParameters, py_Key) ){
        status = 0;
        py_Array = PyDict_GetItem( p_analysis->py_MiliParameters, py_Key );
        num_values = PySequence_Length( py_Array );
        *size = num_values;
        float_pointer_from_pyobject( py_Array, num_values, (float*) array[0] );
    }
    else{
        status = 1;
    }
    return status;
}


/* TAG( mili_reader_read_param_array )
 *
 * Get an array of floats from the mili reader.
 */
int
mili_reader_read_param_array( int dbid, char* key, void** array )
{
    PyObject * py_Array;
    PyObject * py_Char;
    PyObject * py_Key;
    int status;
    int num_values;
    Analysis * p_analysis = get_analy_ptr();

    /* Lookup parameter in mili and convert to string */
    py_Key = string_to_pyobject(key);

    if( PyDict_Contains( p_analysis->py_MiliParameters, py_Key) ){
        status = 0;
        py_Array = PyDict_GetItem( p_analysis->py_MiliParameters, py_Key );
        num_values = PySequence_Length( py_Array );
        float_pointer_from_pyobject( py_Array, num_values, (float*) array[0] );
    }
    else{
        status = 1;
    }
    return status;
}


/* TAG( mili_reader_get_free_node_data )
 *
 * Load in the Free Node mass and vol (If it exists) from the mili reader.
 */
int
mili_reader_get_free_node_data( Analysis * analy, float ** free_node_mass, float ** free_node_volume )
{
    int i, j;
    int qty_nodes;
    int * node_counts;
    int ** node_map;
    Mesh_data * p_mesh;
    float * temp_free_node_mass;
    float * temp_free_node_volume;

    PyObject* py_FreeNodeData = get_pyobject_attribute( analy->py_MiliDB, FREE_NODE_DATA );
    PyObject* py_FNMass = get_pyobject_attribute( py_FreeNodeData, "mass" );
    PyObject* py_FNVol = get_pyobject_attribute( py_FreeNodeData, "vol" );
    PyObject * py_Mass,
             * py_Vol;
    PyObject * py_Float;

    p_mesh = MESH_P( analy );
    qty_nodes = p_mesh->node_geom->qty;
    node_counts = p_mesh->index_map->node_count;
    node_map = p_mesh->index_map->node_map;
    
    temp_free_node_mass = NEW_N( float, qty_nodes, "FN Mass Array" );
    temp_free_node_volume = NEW_N( float, qty_nodes, "FN Volume Array" );

    for( i = 0; i < analy->proc_count; i++ ){
        py_Mass = PySequence_ITEM( py_FNMass, i );
        py_Vol = PySequence_ITEM( py_FNVol, i );

        if( py_Mass == Py_None || py_Vol == Py_None ){
            free( temp_free_node_mass );
            free( temp_free_node_volume );
            return 0;
        }
        else{
            // Nodal Mass
            for( j = 0; j < node_counts[i]; j++ ){
                py_Float = PySequence_ITEM( py_Mass, j );
                temp_free_node_mass[node_map[i][j]] = (float) PyFloat_AsDouble( py_Float );
            }
            // Nodal Volume
            for( j = 0; j < node_counts[i]; j++ ){
                py_Float = PySequence_ITEM( py_Vol, j );
                temp_free_node_volume[node_map[i][j]] = (float) PyFloat_AsDouble( py_Float );
            }
        }
    }

    free_node_mass[0] = temp_free_node_mass;
    free_node_volume[0] = temp_free_node_volume;

    Py_DECREF( py_FreeNodeData );
    Py_DECREF( py_FNMass );
    Py_DECREF( py_FNVol );
    Py_DECREF( py_Mass );
    Py_DECREF( py_Vol );

    return 1;
}

#endif // HAVE_PARALLEL_READ