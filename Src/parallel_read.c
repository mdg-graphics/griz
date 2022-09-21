/*
 *  IO wrappers for parallel reader using Mili Python Reader
 *
 */
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

#ifdef MILI_READER_TIMING
    clock_t start, end;
    double elapsed;
    start = clock();
#endif

    /* Not necessary, but recommended */
    wchar_t *program_name = Py_DecodeLocale(PROGRAM_NAME, NULL);
    Py_SetProgramName(program_name);

    /* Initialize Python */
    Py_Initialize();

    /* Get correct path to mili reader source. */
    if( analy->mili_reader_src_path )
        py_SrcPath = string_to_pyobject( analy->mili_reader_src_path );
    else
        py_SrcPath = string_to_pyobject( DEFAULT_SRC_PATH );
    
    /* Get correct path to mili reader's bin/site-packages */
    if( analy->mili_reader_venv_bin )
        py_SitePackagesPath = string_to_pyobject( analy->mili_reader_venv_bin );
    else
        py_SitePackagesPath = string_to_pyobject( DEFAULT_VENV_BIN );

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

#ifdef MILI_READER_TIMING
    end = clock();
    elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("\n[griz_python_setup] elapsed = %f\n", elapsed);
#endif

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
mili_reader_db_open( char *path_root, int *p_dbid ){
    PyObject *py_Arglist;
    PyObject *py_PlotFilePath;
    Analysis * p_analysis;

#ifdef MILI_READER_TIMING
    clock_t start, end;
    double elapsed;
    start = clock();
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
    PyTuple_SetItem(py_Arglist, 2, Py_True);         // suppress_parallel = True
    PyTuple_SetItem(py_Arglist, 3, Py_False);        // experimental = False

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
    end = clock();
    elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("[mili_reader_db_open] elapsed = %f\n", elapsed);
#endif

    return OK;
}


/* TAG( mili_reader_query )
 *
 * Query a result using the mili reader and return the PyObject pointer to
 * the result dictionary
 */
PyObject*
mili_reader_query( int state, char* class_name, int ipt, char *result ){
    int i;
    /* Arguments we need for mili reader query function */
    PyObject * py_SvarNames,
             * py_ClassName,
             * py_Material = Py_None,
             * py_Labels = Py_None,
             * py_State,
             * py_Ipts,
             * py_WriteData = Py_None;
    PyObject * py_QueryReturn,
             * py_FuncName,
             * py_TempStr;

    /* Get analysis pointer */
    Analysis * p_analysis = get_analy_ptr();

    /* Check if this is a repeated query */
    if( p_analysis->py_PrevQuery != NULL
            && strcmp(p_analysis->prev_query_result_name, result) == 0
            && strcmp(p_analysis->prev_query_class_name, class_name) == 0
            && (p_analysis->prev_query_state == state || p_analysis->prev_query_state == -1)
            && p_analysis->prev_query_ipt == ipt){
        /* Just use saved result */
        py_QueryReturn = p_analysis->py_PrevQuery;
        Py_INCREF( py_QueryReturn );
    }
    /* If not, perform query */
    else{
        /* Increment reference counters for Py_None */
        Py_INCREF( py_Material );
        Py_INCREF( py_Labels );
        Py_INCREF( py_WriteData );

        /* Build list of state variable names. */
        py_SvarNames = PyList_New(1);
        py_TempStr = string_to_pyobject(result);
        PyList_SetItem(py_SvarNames, 0, py_TempStr);

        py_ClassName = string_to_pyobject(class_name);

        /* Convert state number to python integer */
        py_State = Py_BuildValue("i", state);
        /* Do Not Remove this: For some reason the above function call causes an error,
        * even though it correctly returns the expected pyobject*. We need to clear the error
        * Indicator so that the next call does not fail because of it. */
        PyErr_Clear();

        /* Handle integration point */
        py_Ipts = Py_None;
        Py_INCREF( py_Ipts );
        if( ipt != -1 ){
            Py_DECREF( py_Ipts );
            py_Ipts = Py_BuildValue("i", ipt);
            PyErr_Clear();
        }

        /* Get function name and call query */
        py_FuncName = string_to_pyobject( QUERY_FUNCTION );
        py_QueryReturn = PyObject_CallMethodObjArgs( p_analysis->py_MiliDB, py_FuncName, py_SvarNames,
                                                    py_ClassName, py_Material, py_Labels,
                                                    py_State, py_Ipts, py_WriteData, NULL );

        /* Free python refereces */
        Py_DECREF( py_TempStr );
        Py_DECREF( py_FuncName );
        Py_DECREF( py_SvarNames );
        Py_DECREF( py_ClassName );
        Py_DECREF( py_Material );
        Py_DECREF( py_Labels );
        Py_DECREF( py_State );
        Py_DECREF( py_Ipts );
        Py_DECREF( py_WriteData );
    }

    /* Record query information */
    p_analysis->prev_query_result_name = (char*) strdup(result);
    p_analysis->prev_query_class_name = (char*) strdup(class_name);
    p_analysis->prev_query_state = state;
    p_analysis->prev_query_ipt = ipt;

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
get_state_index( PyObject * py_Dict, char * result_name, int state_no ){
    int j;
    int st_qty;
    int state_idx = 0;
    PyObject * py_States;
    
    py_States = result_dictionary_unwrap_states( py_Dict, result_name );
    if( py_States != NULL){
        /* Get quantity of states in the query */
        st_qty = PySequence_Length( py_States );
        if( st_qty == 1 ){
            state_idx = 0;
        }
        else{
            /* Find correct index of state */
            int * states = NEW_N( int, st_qty, "List of states in query" );
            integer_pointer_from_pyobject( py_States, st_qty, states );
            for( j = 0; j < st_qty; j++ ){
                if( states[j] == state_no ){
                    state_idx = j;
                    break;
                }
            }
        }
    }

    return state_idx;
}


/* TAG( mili_reader_extract_query_data )
 *
 * Extract the results from a PyObject pointer to Mili reader query result.
 */
int
mili_reader_extract_query_data( PyObject * py_QueryReturn, char * result_name, int queried_state,
                                char * subrecord_name, Bool_type node_result,
                                int elem_class_idx, void* data )
{
    int i, j, k;
    int comp_qty;
    int out_idx;
    int elem_qty;
    int obj;
    int * object_ids;
    float * result_buffer;
    PyObject * py_ProcResult;
    PyObject * py_Results;
    PyObject * py_Value;
    PyObject * py_Dict;
    PyObject * py_Layout;
    PyObject * py_States;
    int st_qty, state_idx;
    int * index_map;

    /* Get global analysis pointer */
    Analysis * p_analysis = get_analy_ptr();

    out_idx = 0;
    state_idx = -1;
    result_buffer = (float*) data;

    if( node_result ){
        /* This array is always going to be bigger than the number of elements on a given processor.
        * Allocating a larger array may be more efficient than freeing and reallocating an exact sized
        * array for each processor. */
        int max_elem_qty = 0;
        for( i = 0; i < p_analysis->proc_count; i++ ){
            elem_qty = MESH_P(p_analysis)->index_map->node_count[i];
            if( elem_qty > 0 ){
                if( elem_qty > max_elem_qty )
                    max_elem_qty = elem_qty;
            }
        }
        object_ids = NEW_N( int, max_elem_qty, "Node result object ids");

        /* Loop over processors and get each result */
        for( i = 0; i < p_analysis->proc_count; i++ ){
            /* Get results from the ith processor */
            py_Dict = PySequence_GetItem( py_QueryReturn, i );
            /* Get results for current subrecord */
            py_ProcResult = result_dictionary_unwrap_data( py_Dict, result_name, subrecord_name );
            /* If results exist... */
            if( py_ProcResult != NULL ){
                /* Determine index of state data in query (Only done once) */
                if( state_idx == -1 ){
                    /* Get index of state in results dictionary */
                    state_idx = get_state_index( py_Dict, result_name, queried_state );
                }

                /* Get data for the requested state */
                py_ProcResult = PySequence_GetItem( py_ProcResult, state_idx );

                /* Get object ids and element count */
                py_Layout = result_dictionary_unwrap_layout( py_Dict, result_name, subrecord_name );
                elem_qty = PySequence_Length( py_ProcResult );
                integer_pointer_from_pyobject( py_Layout, elem_qty, object_ids );

                /* Get map from subrec node ids to global positions */
                index_map = MESH_P(p_analysis)->index_map->node_map[i];

                /*  Extract and properly order results. */
                for( j = 0; j < elem_qty; j++){
                    obj = object_ids[j];
                    py_Results = PySequence_GetItem( py_ProcResult, j );
                    comp_qty = PySequence_Length( py_Results );
                    for( k = 0; k < comp_qty; k++ ){
                        out_idx = (index_map[obj] * comp_qty) + k;
                        py_Value = PySequence_GetItem( py_Results, k );
                        result_buffer[out_idx] = (float) PyFloat_AsDouble( py_Value );
                    }
                }
            }
        }

        if( object_ids )
            free( object_ids );
    }
    else{
        /* Loop over processors and get each result */
        for( i = 0; i < p_analysis->proc_count; i++ ){
            py_Dict = PySequence_GetItem( py_QueryReturn, i );
            py_ProcResult = result_dictionary_unwrap_data( py_Dict, result_name, subrecord_name );
            if( py_ProcResult != NULL ){
                /* Determine index of state data in query (Only done once) */
                if( state_idx == -1 ){
                    /* Get index of state in results dictionary */
                    state_idx = get_state_index( py_Dict, result_name, queried_state );
                }

                py_ProcResult = PySequence_GetItem( py_ProcResult, state_idx );
                elem_qty = PySequence_Length( py_ProcResult );

                /*  Extract and properly order results. */
                for( j = 0; j < elem_qty; j++){
                    py_Results = PySequence_GetItem( py_ProcResult, j );
                    comp_qty = PySequence_Length( py_Results );
                    for( k = 0; k < comp_qty; k++ ){
                        py_Value = PySequence_GetItem( py_Results, k );
                        result_buffer[out_idx++] = (float) PyFloat_AsDouble( py_Value );
                    }
                }
            }
        }
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
    int ipt;
    Bool_type node_result;
    char* subrecord_name;
    char* class_name;
    State_rec_obj* p_sro;
    Subrec_obj* p_subrec;
    PyObject * py_QueryReturn;
    PyObject * py_NewResultFlag;
    PyObject * py_Result;
    Bool_type new_result_flag;

    /* Get analysis pointer */
    Analysis * p_analysis = get_analy_ptr();

#ifdef MILI_READER_TIMING
    clock_t start, end;
    double elapsed;
    start = clock();
#endif

    /* Get srec id and subrecord name */
    srec_id = p_analysis->state_srec_fmt_ids[state-1];
    p_sro = p_analysis->srec_tree + srec_id;
    p_subrec = p_sro->subrecs + subrec_id;
    subrecord_name = p_subrec->subrec.name;
    class_name = p_subrec->p_object_class->short_name;

    ipt = -1;
    if( p_subrec->element_set ){
        if(p_subrec->element_set->tempIndex < 0)
            ipt = p_subrec->element_set->integration_points[p_subrec->element_set->current_index];
        else
            ipt = p_subrec->element_set->integration_points[p_subrec->element_set->tempIndex];
    }

    /* This technically breaks if reading multiple results, but that currently never happens */
    for( i = 0; i < qty; i++ ){
        /* Query result using the Mili reader */
        py_QueryReturn = mili_reader_query( state, class_name, ipt, results[i] );
        if( py_QueryReturn == NULL ){
            return NOT_OK;
        }

        /* Update the stored previous query */
        if( p_analysis->py_PrevQuery != NULL )
            Py_DECREF( p_analysis->py_PrevQuery );
        p_analysis->py_PrevQuery = py_QueryReturn;
        Py_INCREF( p_analysis->py_PrevQuery );

        /* Extract/reorder the results from the Python dictionary structure */
        node_result = (p_subrec->p_object_class->superclass == G_NODE);
        elem_class_idx = p_subrec->p_object_class->elem_class_index;
        rval = mili_reader_extract_query_data( py_QueryReturn, results[i], state, subrecord_name, node_result,
                                               elem_class_idx, data );
        
        Py_DECREF( py_QueryReturn );
    }

#ifdef MILI_READER_TIMING
    end = clock();
    elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("[mili_reader_get_results] elapsed = %f\n", elapsed);
#endif


    return rval;
}


/* TAG( cmpfunc )
 *
 * Function to compare two integers
 */
int cmpfunc (const void * a, const void * b) {
   return ( *(int*)a - *(int*)b );
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
    int superclass;
    int block_qty;
    int block_index;
    int start_ident, stop_ident;
    int class_qty;
    int * block_range_ptr;
    int *node_labels, *elem_labels;
    int elem_nodes[20];
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
            
    if( analy->mesh_table != NULL)
    {
        popup_dialog( WARNING_POPUP, "Mesh table pointer not NULL as initialization." );
        return -1;
    }

    int proc_count = analy->proc_count;

#ifdef MILI_READER_TIMING
    clock_t start, end;
    double elapsed;
    start = clock();
#endif

    /* Single call to Mili reader to get all the data needed by mili_reader_get_geom */
    PyObject * py_GetGeomCallData = call_mili_reader_function_noargs( analy->py_MiliDB, GET_GEOM_CALL );
    if ( py_GetGeomCallData == NULL ){
        fprintf( stderr, "mili_reader_get_geom calling GET_GEOM_CALL returned NULL\n");
        return GRIZ_FAIL;
    }

    /* Split out data by processor for easier access later on */
    PyObject **py_MOclass_by_proc = (PyObject**) malloc(proc_count * sizeof(PyObject*));
    PyObject **py_nodes_by_proc = (PyObject**) malloc(proc_count * sizeof(PyObject*));
    PyObject **py_labels_by_proc = (PyObject**) malloc(proc_count * sizeof(PyObject*));
    PyObject **py_conns_by_proc = (PyObject**) malloc(proc_count * sizeof(PyObject*));
    PyObject **py_mats_by_proc = (PyObject**) malloc(proc_count * sizeof(PyObject*));
    PyObject **py_parts_by_proc = (PyObject**) malloc(proc_count * sizeof(PyObject*));

    for( i = 0; i < proc_count; i++ ){
        py_ProcData = PySequence_GetItem( py_GetGeomCallData, i );
        if( py_ProcData == NULL ){
            fprintf( stderr, "mili_reader_get_geom - Geometry data is Null for processor %d\n", i );
            return GRIZ_FAIL;
        }

        /* Get data */
        py_MOclass_by_proc[i] = get_pyobject_attribute( py_ProcData, "mo_classes");
        py_nodes_by_proc[i]   = get_pyobject_attribute( py_ProcData, "nodes");
        py_labels_by_proc[i]  = get_pyobject_attribute( py_ProcData, "labels");
        py_conns_by_proc[i]   = get_pyobject_attribute( py_ProcData, "connectivity");
        py_mats_by_proc[i]   = get_pyobject_attribute( py_ProcData, "materials");
        py_parts_by_proc[i]   = get_pyobject_attribute( py_ProcData, "parts");
    }

    /* Get dimensions of the database */
    dims = get_pyobject_attribute_as_int( analy->py_MiliDB, "dimensions" );
    if ( dims != 2 && dims != 3 ){
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
        for( j = 0; j < proc_count; j++ ){
            // Get Mesh Object class data for 'j' processor
            py_MOData = py_MOclass_by_proc[j];
            if ( py_MOData == NULL )
                continue;
            num_classes = PySequence_Length( py_MOData );

            for( k = 0; k < num_classes; k++ ){
                py_ElementClass = PySequence_GetItem( py_MOData, k );
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
                    if( elem_class->superclass != G_MESH && elem_class->superclass != G_MAT && elem_class->superclass != G_UNIT ){
                        obj_qty = get_pyobject_attribute_as_int( py_ElementClass, "elem_qty" );
                        elem_class->qty += obj_qty;
                    }
                }
            }
        }
        
        /* ----------------------------------------------------------------------------------------- */

        /* Look up node class */
        rval = htable_search( p_ht, "node", FIND_ENTRY, &p_hte );
        node_class = (MO_class_data*) p_hte->data;

        /* Initialize mapping from processor to global */
        p_md->index_map = NEW( ProcessorToGlobalMap, "Proc to Global index map");
        p_md->index_map->node_count = (int*) malloc( proc_count * sizeof(int));
        p_md->index_map->node_map = (int **) malloc( proc_count * sizeof(int*));

        /* Populate data for nodal class */
        if ( rval == OK ){
            obj_qty = node_class->qty;
            /* Get all labels into a single list */
            node_labels = (int*) malloc(obj_qty * sizeof(int));
            pos = 0;
            /* For each processor */
            for( j = 0; j < proc_count; j++ ){
                /* Get node labels on processor */
                py_ProcLabels = PyDict_GetItemString( py_labels_by_proc[j], node_class->short_name );
                nodes_on_proc = PySequence_Length( py_ProcLabels );
                if( nodes_on_proc > 0 ){
                    /* Convert node labels to integer array */
                    integer_pointer_from_pyobject( py_ProcLabels, nodes_on_proc, node_labels+pos);
                    pos += nodes_on_proc;

                    /* Update arrays */
                    p_md->index_map->node_count[j] = nodes_on_proc;
                    p_md->index_map->node_map[j] = (int*) malloc( nodes_on_proc * sizeof(int) );
                }
                else{
                    /* Update arrays */
                    p_md->index_map->node_count[j] = 0;
                    p_md->index_map->node_map[j] = NULL;
                }
            }

            /* Sort the labels */
            qsort( node_labels, obj_qty, sizeof(int), cmpfunc);

            /* Remove duplicates and get real number of nodes */
            k = 0;
            for( j = 0; j < obj_qty-1; j++ ){
                if ( node_labels[j] != node_labels[j+1])
                    node_labels[k++] = node_labels[j];
            }
            node_labels[k++] = node_labels[obj_qty-1];

            obj_qty = k;
            node_class->qty = obj_qty;

            /* Create per processor mapping for nodes to index in global node list */
            for( j = 0; j < proc_count; j++ ){
                py_ProcLabels = PyDict_GetItemString( py_labels_by_proc[j], node_class->short_name );
                nodes_on_proc = p_md->index_map->node_count[j];

                if( nodes_on_proc > 0 ){
                    /* Store label numbers in array temporarily, will overwrite next. */
                    integer_pointer_from_pyobject( py_ProcLabels, nodes_on_proc, p_md->index_map->node_map[j] );

                    for( k = 0; k < nodes_on_proc; k++ ){
                        target = p_md->index_map->node_map[j][k];
                        /* ----- Binary Search ----- */
                        p_md->index_map->node_map[j][k] = bin_search_index(target, node_labels, obj_qty);
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
            for( j = 0; j < proc_count; j++ ){
                nodes_on_proc = p_md->index_map->node_count[j];
                py_ProcNodes = py_nodes_by_proc[j];

                /* Loop over all nodes and coordinates */
                for( k = 0; k < nodes_on_proc; k++ ){
                    node_idx = p_md->index_map->node_map[j][k];
                    out_idx = node_idx * dims;
                    in_idx = k * dims;
                    py_CoordsTuple = PySequence_GetItem( py_ProcNodes, k );
                    float_pointer_from_pyobject( py_CoordsTuple, dims, node_coords+out_idx);
                }
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

            for( obj_id = 0; obj_id < obj_qty; obj_id++ ){
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
            for ( obj_id = 0; obj_id < obj_qty; obj_id++ ){
                node_class->labels_index[node_class->labels[obj_id].local_id] = obj_id; 
            }

            /* Construct blocking data for node labels */
            block_range_ptr = compute_label_blocking_data( node_labels, obj_qty, &block_qty);
            if( block_range_ptr == NULL ){
                fprintf(stderr, "mili_reader_load_nodal_data call compute_label_blocking_data\n");
                return GRIZ_FAIL;
            }

            /* Construct the label blocking table of contents */
            if( node_class->labels_found && block_qty > 0 && block_range_ptr ){
                node_class->label_blocking.block_qty = block_qty;
                node_class->label_blocking.block_total_objects = 0;
                node_class->label_blocking.block_min = MAXINT;
                node_class->label_blocking.block_max = MININT;
                node_class->label_blocking.block_objects = NEW_N( Label_block_data, block_qty, "Node Class Label Blocking Objects");

                block_index = 0;
                for( k = 0; k < block_qty; k++ ){
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

            if ( block_range_ptr ){
                free( block_range_ptr );
                block_range_ptr = NULL;
            }
            if ( node_labels ){
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
        else{
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

        /* ----------------------------------------------------------------------------------------- */

        /* For each element class load labels, connectivity, parts, materials, etc. */
        htable_get_data( p_ht, (void***) &htable_mo_classes, &elem_class_qty);

        /* Initialize arrays for storing element counts per processor and processor to global index mapping */
        p_md->index_map->elem_count = (int**) malloc( elem_class_qty * sizeof(int*) );
        p_md->index_map->elem_map = (int***) malloc( elem_class_qty * sizeof(int**) );
        for( j = 0; j < elem_class_qty; j++ ){
            p_md->index_map->elem_count[j] = (int*) malloc( proc_count * sizeof(int) );
            p_md->index_map->elem_map[j] = (int**) malloc( proc_count * sizeof(int*) );
        }

        /* Populate element class data. */
        for( j = 0; j < elem_class_qty; j++ ){
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
            for( k = 0; k < proc_count; k++ ){
                /* Get Labels on processor */
                py_ProcLabels = PyDict_GetItemString( py_labels_by_proc[k], elem_class->short_name );
                if( py_ProcLabels == NULL ){
                    p_md->index_map->elem_count[elem_class->elem_class_index][k] = 0;
                    p_md->index_map->elem_map[elem_class->elem_class_index][k] = NULL;
                }
                else{
                    elems_on_proc = PySequence_Length( py_ProcLabels );
                    /* Update arrays */
                    p_md->index_map->elem_count[elem_class->elem_class_index][k] = elems_on_proc;
                    p_md->index_map->elem_map[elem_class->elem_class_index][k] = (int*) malloc( elems_on_proc * sizeof(int) );

                    /* If labels exist */
                    if( idents_exist ){
                        /* Read in element labels */
                        integer_pointer_from_pyobject( py_ProcLabels, elems_on_proc, elem_labels+pos );
                    }
                    else{
                        /* Need to create our own labels -- 1->N */
                        for(l = 0; l < elems_on_proc; l++){
                            elem_labels[temp] = temp+1;
                            temp++;
                        }
                    }

                    /* Map proc elems to global indices */
                    if( elem_class->superclass == M_MESH || elem_class->superclass == M_MAT || elem_class->superclass == M_UNIT ){
                        /* This is an awkward solution to MAT, MESH, and M_UNIT (maybe) being the same on
                        * all processors, can definitely make this better later */
                        for( l = 0; l < elems_on_proc; l++ ){
                            p_md->index_map->elem_map[j][k][l] = l;
                        }
                    }
                    else{
                        for( l = 0; l < elems_on_proc; l++ ){
                            p_md->index_map->elem_map[j][k][l] = pos++;
                        }
                    }
                }
            }

            /* Special case for particle classes -- Elements need to be put in sorted order. */
            if( idents_exist && is_particle_class( analy, elem_class->superclass, elem_class->short_name ) ){
                /* Sort element labels */
                qsort( elem_labels, obj_qty, sizeof(int), cmpfunc );
                for( k = 0; k < proc_count; k++ ){
                    /* Get Labels on processor */
                    py_ProcLabels = PyDict_GetItemString( py_labels_by_proc[k], elem_class->short_name );
                    elems_on_proc = p_md->index_map->elem_count[j][k];
                    if( elems_on_proc > 0 ){
                        /* Store label numbers in array temporarily, will overwrite next. */
                        integer_pointer_from_pyobject( py_ProcLabels, elems_on_proc,
                                                        p_md->index_map->elem_map[j][k] );
                        
                        for( l = 0; l < elems_on_proc; l++ ){
                            target = p_md->index_map->elem_map[j][k][l];
                            p_md->index_map->elem_map[j][k][l] = bin_search_index( target, elem_labels, obj_qty );
                        }
                    }
                }
            }


            if( elem_class->superclass == M_MESH || elem_class->superclass == M_MAT || elem_class->superclass == M_UNIT ){

                elem_class->simple_start = elem_labels[0];
                elem_class->simple_stop = elem_labels[elem_class->qty-1];

                if ( elem_class->superclass == M_MAT ){
                    p_matd = NEW_N( Material_data, obj_qty, "Material data array" );
                    elem_class->objects.materials = p_matd;
                }

                p_lh = p_md->classes_by_sclass + elem_class->superclass;
                mo_classes = (MO_class_data **) p_lh->list;
                mo_classes = (void *) RENEW_N( MO_class_data *, mo_classes, p_lh->qty, 1, "Extend classes_by_sclass array");
                mo_classes[p_lh->qty] = elem_class;
                p_lh->qty++;
                p_lh->list = (void *) mo_classes;
            }
            else{
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
                for( k = 0; k < proc_count; k++ ){
                    if( p_md->index_map->elem_count[j][k] > 0 ){
                        elems_on_proc = p_md->index_map->elem_count[j][k];

                        /* Load connectivity */
                        py_ProcConns = PyDict_GetItemString( py_conns_by_proc[k], elem_class->short_name );
                        /* Load materials and parts */
                        py_ProcMats = PyDict_GetItemString( py_mats_by_proc[k], elem_class->short_name );
                        py_ProcParts = PyDict_GetItemString( py_parts_by_proc[k], elem_class->short_name );

                        int out_idx;
                        PyObject * py_ConnSlice;
                        for( l = 0; l < elems_on_proc; l++ ){
                            /* Get elements nodal connectivity */
                            py_ConnSlice = PySequence_GetItem( py_ProcConns, l );
                            integer_pointer_from_pyobject( py_ConnSlice, qty_conns, elem_nodes );

                            /* Elem nodes is local id's for nodes need to convert to global labels */
                            out_idx = p_md->index_map->elem_map[j][k][l] * qty_conns;
                            for( m = 0; m < qty_conns; m++ ){
                                p_ed->nodes[out_idx+m] = p_md->index_map->node_map[k][elem_nodes[m]];
                            } 

                            out_idx = p_md->index_map->elem_map[j][k][l];
                            /* Get elements material and number */
                            py_M = PySequence_GetItem( py_ProcMats, l );
                            py_P = PySequence_GetItem( py_ProcParts, l );
                            p_ed->mat[out_idx] = PyLong_AsLong(py_M) - 1;
                            p_ed->part[out_idx] = PyLong_AsLong(py_P) - 1;
                        }
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

        // Need to call gen_material_data on M_MAT classes after all other classes are processed.
        MO_class_data ** mat_classes = p_md->classes_by_sclass[M_MAT].list;
        for( j = 0; j < p_md->classes_by_sclass[M_MAT].qty; j++ ){
            elem_class = mat_classes[i];
            gen_material_data( elem_class, p_md );
        }

        /* Update the Mesh_data struct with element class info */
        p_md->elem_class_qty = elem_class_qty;

    } // for ( i = 0 i < mesh_qty; i++ )

    get_hex_volumes( analy->db_ident, analy );

#ifdef MILI_READER_TIMING
    end = clock();
    elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("[mili_reader_get_geom] elapsed = %f\n", elapsed);
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
            integer_pointer_from_pyobject( py_Dims, p_sv->rank, p_sv->dims );
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
                py_Name = PySequence_GetItem( py_CompNames, i );
                py_Title = PySequence_GetItem( py_CompNames, i );

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


/* TAG( set_es_middle_index )
 *
 * Given an ElementSet, compute and sets the middle index.
 */
Bool_type
set_es_middle_index( ElementSet * element_set )
{
    int middle, j;

    /* We need to determine the middle integration point for the default.*/
    /* Take the actual number of point used in the simulation to find the middle*/
    middle = element_set->integration_points[element_set->size]/2.0 + 
                element_set->integration_points[element_set->size]%2;
    
    for( j = 0; j < element_set->size; j++)
    {
        if( middle < element_set->integration_points[j] )
            break;
    }

    if( j > 0 && j < element_set->size )
    {
        if( j > 1)
            j--;

        if( abs(element_set->integration_points[j] - middle)< abs(element_set->integration_points[j-1] - middle) )
            element_set->middle_index = j;
        else
            element_set->middle_index = j-1;
    }
    else if( j == element_set->size)
        element_set->middle_index = j-1;
    else
        element_set->middle_index = j;

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

    if ( py_ESList == NULL ){
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
                py_Tuple = PySequence_GetItem( py_DictItems, j );
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
                    integer_pointer_from_pyobject( py_Value, ipt_count, element_set->integration_points );

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
    int *node_work_array=NULL;
    State_variable *p_svar;
    Mesh_data *p_mesh;

    Hash_table * p_subrecord_ht;
    Hash_table * p_subrecord_gids_ht;
    Hash_table * p_sv_ht, * p_primal_ht;
     
    PyObject * py_CallData, * py_Struct;
    PyObject ** py_ESList;
    PyObject * py_Subrecords,
             * py_Subrecord,
             * py_SrecQty;
    PyObject * py_SubrecordStateVariables,
             * py_StateVariable,
             * py_ComponentSvars,
             * py_ComponentSvar;

#ifdef MILI_READER_TIMING
    clock_t start, end;
    double elapsed;
    start = clock();
#endif

    analy->num_bad_subrecs = 0;
    analy->bad_subrecs = NULL;
    analy->old_shell_stresses = FALSE;

    /* Single call to mili reader to get all the data we need */
    py_CallData = call_mili_reader_function_noargs( analy->py_MiliDB, GET_ST_DESCRIPTORS_CALL );
    if( py_CallData == NULL ){
        fprintf(stderr, MR_DEFAULT_FAILED, "mili_reader_get_st_descriptors", GET_ST_DESCRIPTORS_CALL);
        return GRIZ_FAIL;
    }

    /* Get data for each processor */
    int proc_count = analy->proc_count;
    int * subrecord_count_per_proc = (int*) malloc( proc_count * sizeof(int));
    PyObject **py_subrecord_list_by_proc = (PyObject**) malloc(proc_count * sizeof(PyObject*));
    py_ESList = NEW_N( PyObject*, proc_count, "PyObject Element set List" );

    /* Store data for each processor */
    for( i = 0; i < proc_count; i++ ){
        py_Struct = PySequence_GetItem( py_CallData, i );
        if ( py_Struct != NULL && py_Struct != Py_None ){
            py_ESList[i] = get_pyobject_attribute( py_Struct, "element_sets");
            py_subrecord_list_by_proc[i] = get_pyobject_attribute( py_Struct, "subrecords" );
            subrecord_count_per_proc[i] = PySequence_Length( py_subrecord_list_by_proc[i] );
        }
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
            py_Subrecord = PySequence_GetItem( py_Subrecords, j );
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
                    py_StateVariable = PySequence_GetItem( py_SubrecordStateVariables, k );
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
                                py_ComponentSvar = PySequence_GetItem( py_ComponentSvars, l );
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
    node_work_array = NEW_N( int, mesh_node_qty, "Temp node array" );

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

                    p_sro[srec_id].node_pos_subrec = subrec_index;
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
    end = clock();
    elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("[mili_reader_get_st_descriptors] elapsed = %f\n", elapsed);
#endif
    return OK;
}


/* TAG( combine_nodpos )
 *
 * Combine the nodal positions from multiple processors
 */
extern Bool_type
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
        py_ProcNodes = PySequence_GetItem( analy->py_Nodpos, i);

        py_Nodes = PyDict_GetItemString( py_ProcNodes, nodpos );
        py_Nodes = PyDict_GetItemString( py_Nodes, data );
        py_Nodes = PyDict_GetItemString( py_Nodes, node );

        /* Single state, so remove state array */
        py_Nodes = PySequence_GetItem( py_Nodes, state_no );
        node_qty = MESH_P(analy)->index_map->node_count[i];

        /* Get coordinates for each element */
        for( j = 0; j < node_qty; j++ ){
            py_Coords = PySequence_GetItem( py_Nodes, j );
            out_idx = MESH_P(analy)->index_map->node_map[i][j] * dims;
            float_pointer_from_pyobject( py_Coords, dims, out_buf+out_idx);
        }
    }

    return OK;
}


/* TAG( combine_sand_result )
 *
 * Combine the sand flags from multiple processors
 */
extern Bool_type
combine_sand_flags( Analysis *analy, Subrec_obj * p_subrec, int state_no, void * out_buffer ){
    int i, j;
    int out_idx;
    int elem_qty;
    char * data = "data";
    char * sand = "sand";
    char * subrecord_name;
    char * class_name;
    int elem_class_index;
    int obj;
    PyObject * py_ClassName;
    PyObject * py_SubrecordName;
    PyObject * py_FloatValue;
    
    PyObject * py_ProcResult;
    PyObject * py_Result;
    PyObject * py_SandFlags;
    PyObject * py_Layout;
    PyObject * py_ClassSandData;

    float * out_buf = (float*) out_buffer;

    subrecord_name = p_subrec->subrec.name;
    class_name = p_subrec->p_object_class->short_name;
    elem_class_index = p_subrec->p_object_class->elem_class_index;

    py_ClassName = string_to_pyobject( class_name );
    py_SubrecordName = string_to_pyobject( subrecord_name );

    out_idx = 0;
    for( i = 0; i < analy->proc_count; i++ ){
        py_ProcResult = PySequence_GetItem( analy->py_Sand, i );

        if( PyDict_Contains( py_ProcResult, py_ClassName ) ){
            py_ClassSandData = PyDict_GetItemString( py_ProcResult, class_name );

            /* Get list of result values */
            py_Result = PyDict_GetItemString( py_ClassSandData, sand );
            py_SandFlags = PyDict_GetItemString( py_Result, data );

            /* Check to make sure this subrecord exists on the current processor */
            if( PyDict_Contains( py_SandFlags, py_SubrecordName) ){
                py_SandFlags = PyDict_GetItemString( py_SandFlags, subrecord_name );
                /* Single state, so remove state array */
                py_SandFlags = PySequence_GetItem( py_SandFlags, state_no );
            
                /* Get List of element ids */
                elem_qty = MESH_P( analy )->index_map->elem_count[elem_class_index][i];

                /* Get sand flag for each element */
                for( j = 0; j < elem_qty; j++ ){
                    py_FloatValue = PySequence_GetItem( py_SandFlags, j );
                    py_FloatValue = PySequence_GetItem( py_FloatValue, 0 ); // Stored as single value list
                    out_buf[out_idx++] = (float) PyFloat_AsDouble( py_FloatValue );
                }
            }
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
    clock_t start, end;
    double elapsed;
    start = clock();
#endif

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
    if ( analy->stateDB )
    {
        rval = combine_nodpos( analy, state_no, p_st->nodes.nodes );

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
    for( i = 0; i < p_sro->qty; i++ ){
        p_subrec = p_subrecs + i;
        ec_index = p_subrec->p_object_class->elem_class_index;
        if( p_subrec->sand ){
            rval = combine_sand_flags( analy, p_subrec, state_no, p_st->elem_class_sand[ec_index]);
            if( rval != OK ){
                fprintf( stderr, "mili_reader_get_state call combine_sand_flags\n");
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
        rval = mili_reader_get_results( analy->db_ident, state_no, p_st->sph_srec_id, 1, &primal, input_buf );

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
    end = clock();
    elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("[mili_reader_get_state] elapsed = %f\n", elapsed);
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
#ifdef MILI_READER_TIMING
    clock_t start, end;
    double elapsed;
    start = clock();
#endif

    mili_reader_get_string( analy->py_MiliParameters, "lib version", analy->mili_version );
    mili_reader_get_string( analy->py_MiliParameters, "host name", analy->mili_host );
    mili_reader_get_string( analy->py_MiliParameters, "arch name", analy->mili_arch );
    mili_reader_get_string( analy->py_MiliParameters, "date", analy->mili_timestamp );
    mili_reader_get_string( analy->py_MiliParameters, "xmilics version", analy->xmilics_version );

#ifdef MILI_READER_TIMING
    end = clock();
    elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("[mili_reader_get_metadata] elapsed = %f\n", elapsed);
#endif
    return OK;
}


/* TAG( mili_reader_load_parameters )
 *
 * Get the parameters dictionary of the 0th processor
 */
extern int
mili_reader_load_parameters( Analysis * analy )
{
#ifdef MILI_READER_TIMING
    clock_t start, end;
    double elapsed;
    start = clock();
#endif

    PyObject * py_Parameters = call_mili_reader_function_noargs( analy->py_MiliDB, PARAMETER_GETTER );
    if( py_Parameters == NULL )
        return NOT_OK;
    analy->py_MiliParameters = py_Parameters;

#ifdef MILI_READER_TIMING
    end = clock();
    elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("[mili_reader_load_parameters] elapsed = %f\n", elapsed);
#endif
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
#ifdef MILI_READER_TIMING
    clock_t start, end;
    double elapsed;
    start = clock();
#endif

    Analysis * analy = get_analy_ptr();
    call_mili_reader_function_noargs( analy->py_MiliDB, RELOAD_FUNCTION );

#ifdef MILI_READER_TIMING
    end = clock();
    elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("[mili_reader_reload_states] elapsed = %f\n", elapsed);
#endif
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
    clock_t start, end;
    double elapsed;
    start = clock();
#endif

    if( analy->py_MiliDB || analy->py_MiliReaderModule ){
        Py_FinalizeEx(); 
    }

#ifdef MILI_READER_TIMING
    end = clock();
    elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("[mili_reader_db_close] elapsed = %f\n", elapsed);
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
    clock_t start, end;
    double elapsed;
    start = clock();
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

    /* Get state data (Nodal positions and sand flags) from the mili reader. */
    py_StateData = call_mili_reader_function_noargs( analy->py_MiliDB, GET_STATE_CALL );
    if( py_StateData == NULL ){
        fprintf( stderr, MR_DEFAULT_FAILED, "mili_reader_load_state_data", GET_STATE_CALL );
        return NOT_OK;
    }

    /* Unpack nodal positions */
    analy->py_Nodpos = get_pyobject_attribute( py_StateData, "nodpos" );
    if( analy->py_Nodpos == NULL ){
        fprintf( stderr, "mili_reader_load_state_data unable to unpack nodpos\n" );
        return NOT_OK;
    }

    /* Unpack sand flags */
    analy->py_Sand = get_pyobject_attribute( py_StateData, "sand" );
    if( analy->py_Sand == NULL ){
        fprintf( stderr, "mili_reader_load_state_data unable to unpack sand flags\n" );
        return NOT_OK;
    }

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
        py_SMap = PySequence_GetItem( py_StateMaps, i );
        st_time = get_pyobject_attribute_as_double( py_SMap, "time" );
        analy->state_times[i] = st_time;

        /* There should only be a single state record format */
        analy->state_srec_fmt_ids[i] = 0;
    }


#ifdef MILI_READER_TIMING
    end = clock();
    elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("[mili_reader_load_state_times] elapsed = %f\n", elapsed);
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
            py_Param = PySequence_GetItem( py_Matches, i );
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
        py_Mass = PySequence_GetItem( py_FNMass, i );
        py_Vol = PySequence_GetItem( py_FNVol, i );

        if( py_Mass == Py_None || py_Vol == Py_None ){
            free( temp_free_node_mass );
            free( temp_free_node_volume );
            return 0;
        }
        else{
            // Nodal Mass
            for( j = 0; j < node_counts[i]; j++ ){
                py_Float = PySequence_GetItem( py_Mass, j );
                temp_free_node_mass[node_map[i][j]] = (float) PyFloat_AsDouble( py_Float );
            }
            // Nodal Volume
            for( j = 0; j < node_counts[i]; j++ ){
                py_Float = PySequence_GetItem( py_Vol, j );
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