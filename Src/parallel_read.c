#define PARALLEL_READER_SUPPORT
#ifdef PARALLEL_READER_SUPPORT

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <time.h>

#include "parallel_read.h"
#include "mili.h"
#include "gahl.h"


void griz_python_setup()
{
    /* Not necessary, but recommended */
    wchar_t *program_name = Py_DecodeLocale(PROGRAM_NAME, NULL);
    Py_SetProgramName(program_name);

    Py_Initialize();

    /* Update path so that we can find the mili reader lib */
    PyRun_SimpleString("import sys; sys.path.append('/g/g16/rhathaw/all_codes/mili-python/src/')");

    /* Put venv lib at front of path */
    PyRun_SimpleString("sys.path.insert(0, '/g/g16/rhathaw/all_codes/mili-python/embedded_python_c/lib/python3.7/site-packages/')");
}


//int
//main(int argc, char* argv[])
//{
//    char * plotfile_name = "./parallel/d3samp6.plt";
//    PyObject *pName, *pModule, *pDict, *pFunc;
//    int dbid;
//    int i;
//    int num_states;
//    Analysis *analy;
//
//    analy = (Analysis*) malloc(sizeof(Analysis));
//
//    griz_python_setup();
//
//    /* Import mili_reader_lib */
//    mili_reader_module = PyImport_ImportModule(MILI_READER_MODULE_NAME);
//    /* Check that we successfully found the reader */
//    if( mili_reader_module == NULL ){
//        PyErr_Print();
//        fprintf(stderr, "Failed to load mili reader module\n");
//        exit(1);
//    }
//
//    /* Attempt to open Mili Plot file */
//    if( mili_reader_db_open(plotfile_name, &dbid) != OK ){
//        fprintf(stderr, "Failed to open plotfile: %s\n", plotfile_name);
//        exit(1);
//    }
//
//    analy->dimension = 3; // Just for now...
//
//    int rval = mili_reader_get_geom( analy, dbid, &analy->mesh_table, &analy->mesh_qty );
//    printf( "[mili_reader_get_geom] status = %d\n", rval );
//
//    rval = mili_reader_get_st_descriptors( analy, dbid );
//    printf( "[mili_reader_get_st_descriptors] status = %d\n", rval );
//
//    rval = mili_reader_get_state( analy, 0, analy->state_p, &analy->state_p, &num_states );
//    printf( "[mili_reader_get_state] status = %d\n", rval );
//
//    if( Py_FinalizeEx() < 0 ){
//        exit(1);
//    }
//}


/* TAG( mili_reader_db_open )
 *
 * Open the specified plot files with the mili reader and set mili_db to reference the newly created
 * python reader object.
 */
int
mili_reader_db_open( char *path_root, int *p_dbid ){
    PyObject *py_Arglist;
    PyObject *py_Arg1,
             *py_Arg2,
             *py_Arg3,
             *py_Arg4;

#ifdef MILI_READER_TIMING
    clock_t start, end;
    double elapsed;
    start = clock();
#endif

    /* Build argument list for function call */
    py_Arglist = PyTuple_New(4);

    py_Arg1 = Py_BuildValue("s", path_root);
    py_Arg2 = PyList_New(0);
    py_Arg3 = Py_False;
    py_Arg4 = Py_False;

    PyTuple_SetItem(py_Arglist, 0, py_Arg1);
    PyTuple_SetItem(py_Arglist, 1, py_Arg2);
    PyTuple_SetItem(py_Arglist, 2, py_Arg3);
    PyTuple_SetItem(py_Arglist, 3, py_Arg4);

    /* Call open_database function */
    mili_db = call_mili_module_function(OPEN_DATABASE, py_Arglist);
    if(mili_db == NULL){
        return NOT_OK;
    }

    /* Just set db id to 0*/
    p_dbid = 0;

#ifdef MILI_READER_TIMING
    end = clock();
    elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("[mili_reader_db_open] elapsed = %f\n", elapsed);
#endif

    return OK;
}


/* TAG( mili_reader_get_results )
 *
 * Query the Mili reader for some results
 */
int
mili_reader_get_results( int dbid, int state, int subrec_id, int qty, char **results, void *data )
{
    int i, j;
    int pos, values_read;
    float * result_buffer;
    PyObject * py_FuncName,
             * py_QueryReturn,
             * py_ResultArray,
             * py_ResultName,
             * py_FloatValue,
             * py_TempStr;

    /* Arguments we need for mili reader query function */
    PyObject * py_SvarNames,
             * py_ClassName,
             * py_Material = Py_None,
             * py_Labels = Py_None,
             * py_State,
             * py_Ipts,
             * py_WriteData = Py_None;

#ifdef MILI_READER_TIMING
    clock_t start, end;
    double elapsed;
    start = clock();
#endif

    /* Get arguments as PyObjects */

    /* Build list of state variable names. */
    py_SvarNames = PyList_New(qty);
    for( i = 0; i < qty; i++ ){
        py_TempStr = string_to_pyobject(results[i]);
        PyList_SetItem(py_SvarNames, i, py_TempStr);
    }
    
    /* Get subrecord object so we can get class name and subrecord name */


    /* Convert state number to python integer */
    py_State = Py_BuildValue("i", state);

    /* Set Integration points to None for now...
     * I think I'll have to wait until moving this into griz to access the
     * integration points there. */
    py_Ipts = Py_None;

    /* Get function name and call query */
    py_FuncName = string_to_pyobject( QUERY_FUNCTION );
    py_QueryReturn = PyObject_CallMethodObjArgs( mili_db, py_FuncName, py_SvarNames,
                                                 py_ClassName, py_Material, py_Labels,
                                                 py_State, py_Ipts, py_WriteData, NULL );

    if( py_QueryReturn == NULL ){
        PyErr_Print();
        return NOT_OK;
    }

    /*          TODO  TODO  TODO
     * 1. Object ordered vs results ordered in GRIZ (Does this matter)
     * 2. Performance Evaluation of this
     *
     */

    /*  Extract and properly order results. */
    for( i = 0; i < PROC_COUNT; i++ ){

    }

#ifdef MILI_READER_TIMING
    end = clock();
    elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("[mili_reader_get_results] elapsed = %f\n", elapsed);
#endif

    return OK;
}

int cmpfunc (const void * a, const void * b) {
   return ( *(int*)a - *(int*)b );
}

/* TAG( mili_reader_get_geom )
 *
 * Load in the geometry (mesh definition) data from the mili reader
 */
int
mili_reader_get_geom( Analysis* analy, int dbid, Mesh_data **p_mtable, int *p_mesh_qty )
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
    int nodes_on_proc, elems_on_proc, pos;
    int node_idx, out_idx, in_idx;
    int superclass;
    int block_qty;
    int block_index;
    int start_ident, stop_ident;
    int * block_range_ptr;
    int *node_labels, *elem_labels;
    Bool_type have_nodal;
    char * short_name;
    char * long_name;
    int * node_count_per_proc;
    int ** proc_node_to_global_idx;
    int** elem_count_per_proc;
    int *** proc_elem_to_global_idx;
    Mesh_data *mesh_array;
    Mesh_data *p_md;
    Hash_table *p_ht;
    Htable_entry *p_hte;
    /*
    List_head *p_lh;
    */
    MO_class_data *node_class, *elem_class;
    MO_class_data **mo_classes, **htable_mo_classes;
    Elem_data * p_ed;
    Material_data * p_matd;
    /*
    Analysis *p_analysis;
    */
    PyObject * py_ElementClasses,
             * py_ElementClass,
             * py_ClassDataDict,
             * py_ProcMOData,
             * py_MOData;
    PyObject * py_FuncName;
    PyObject * py_Arg;
    PyObject * py_Mats, * py_Parts;
    PyObject * py_ProcLabels, * py_ProcConns, * py_ProcNodes;
    PyObject * py_ProcMats, * py_ProcParts;
            
    if( *p_mtable != NULL)
    {
        //popup_dialog( WARNING_POPUP, "Mesh table pointer not NULL as initialization." );
        return -1;
    }

#ifdef MILI_READER_TIMING
    clock_t start, end;
    double elapsed;
    start = clock();
#endif

    /* Single call to Mili reader to get all the data needed by mili_reader_get_geom */
    PyObject * py_MOClassData = call_mili_reader_function_noargs( GET_GEOM_CALL );
    if ( py_MOClassData == NULL ){
        fprintf( stderr, "py_MOClassData is NULL\n");
        return GRIZ_FAIL;
    }

    /* Split out data by processor for easier access later on */
    PyObject **py_MOclass_by_proc = (PyObject**) malloc(PROC_COUNT * sizeof(PyObject*));
    PyObject **py_nodes_by_proc = (PyObject**) malloc(PROC_COUNT * sizeof(PyObject*));
    PyObject **py_labels_by_proc = (PyObject**) malloc(PROC_COUNT * sizeof(PyObject*));
    PyObject **py_conns_by_proc = (PyObject**) malloc(PROC_COUNT * sizeof(PyObject*));
    PyObject **py_mats_by_proc = (PyObject**) malloc(PROC_COUNT * sizeof(PyObject*));
    PyObject **py_parts_by_proc = (PyObject**) malloc(PROC_COUNT * sizeof(PyObject*));

    for( i = 0; i < PROC_COUNT; i++ ){
        py_ProcMOData = PySequence_GetItem( py_MOClassData, i );
        if ( py_ProcMOData == NULL )
            continue;

        /* If processor 0 get mesh dimensions */
        if( i == 0 ){
            dims = get_pyobject_attribute_as_int( py_ProcMOData, "dims" );
            printf("dims = %d\n", dims);
        } 

        /* Get data */
        py_MOclass_by_proc[i] = get_pyobject_attribute( py_ProcMOData, "mo_classes" );
        py_nodes_by_proc[i] = get_pyobject_attribute( py_ProcMOData, "nodes" );
        py_labels_by_proc[i] = get_pyobject_attribute( py_ProcMOData, "labels" );
        py_conns_by_proc[i] = get_pyobject_attribute( py_ProcMOData, "connectivity" );
        py_mats_by_proc[i] = get_pyobject_attribute( py_ProcMOData, "materials" );
        py_parts_by_proc[i] = get_pyobject_attribute( py_ProcMOData, "parts" );
    }


    /* For now we'll just assume is one. */
    mesh_qty = 1;

    /* Allocate array of pointer to mesh geom hash tables. */
    mesh_array = NEW_N( Mesh_data, mesh_qty, "Mesh data array" );

    for( i = 0; i < mesh_qty; i++ )
    {
        p_md = mesh_array + i;
        p_ht = htable_create( 151 );
        p_md->class_table = p_ht;

        /* Create initial MO_class_data structs for each element class and sum obj_qty across processors */
        for( j = 0; j < PROC_COUNT; j++ ){
            // Get MO class data for 'j' processor
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

                    /* Generate MO_class_data object for element class */
                    elem_class = NEW( MO_class_data, "Elem geom table entry" );
                    elem_class->mesh_id = i;
                    elem_class->short_name = short_name;
                    elem_class->long_name = long_name;
                    elem_class->superclass = superclass;
                    elem_class->qty = obj_qty;

                    p_hte->data = (void*) elem_class;

                }
                /* Existing mo_class */
                else
                {
                    elem_class = (MO_class_data*) p_hte->data;
                    obj_qty = get_pyobject_attribute_as_int( py_ElementClass, "elem_qty" );
                    elem_class->qty += obj_qty;
                }
            }
        }

        /* ----------------------------------------------------------------------------------------- */

        /* Look up node class */
        rval = htable_search( p_ht, "node", FIND_ENTRY, &p_hte );
        node_class = (MO_class_data*) p_hte->data;

        /* Initialize arrays for storing node counts per processor and processor to global index mapping */
        node_count_per_proc = (int*) malloc( PROC_COUNT * sizeof(int));
        proc_node_to_global_idx = (int **) malloc( PROC_COUNT * sizeof(int*));

        /* Populate data for nodal class */
        if ( rval == OK ){
            obj_qty = node_class->qty;

            /* Get all labels into a single list */
            node_labels = (int*) malloc(obj_qty * sizeof(int));
            pos = 0;
            /* For each processor */
            for( j = 0; j < PROC_COUNT; j++ ){
                /* Get node labels on processor */
                py_ProcLabels = PyDict_GetItemString( py_labels_by_proc[j], node_class->short_name );
                nodes_on_proc = PySequence_Length( py_ProcLabels );

                /* Convert node labels to integer array */
                integer_pointer_from_pyobject( py_ProcLabels, nodes_on_proc, node_labels+pos);
                pos += nodes_on_proc;

                /* Update arrays */
                node_count_per_proc[j] = nodes_on_proc;
                proc_node_to_global_idx[j] = (int*) malloc( nodes_on_proc * sizeof(int) );
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
            for( j = 0; j < PROC_COUNT; j++ ){
                py_ProcLabels = PyDict_GetItemString( py_labels_by_proc[j], node_class->short_name );
                nodes_on_proc = node_count_per_proc[j];

                /* Store label numbers in array temporarily, will overwrite next. */
                integer_pointer_from_pyobject( py_ProcLabels, nodes_on_proc, proc_node_to_global_idx[j] );

                for( k = 0; k < nodes_on_proc; k++ ){
                    target = proc_node_to_global_idx[j][k];
                    /* ----- Binary Search ----- */
                    proc_node_to_global_idx[j][k] = bin_search_index(target, node_labels, obj_qty);
                }
            }

            if( dims == 3 )
                node_class->objects.nodes3d = NEW_N( GVec3D, obj_qty, "3D node coord array");
            else
                node_class->objects.nodes2d = NEW_N( GVec2D, obj_qty, "2D node coord array");

            node_class->objects.nodes = NEW_N(float, obj_qty * 3, "Node Positions");

            node_class->data_buffer = NEW_N( float, obj_qty, "Class data buffer" );
            //if( node_class->data_buffer == NULL )
            //    popup_fatal( "Unable to allocate data buffer on class load" );

            float * proc_coords;
            float * node_coords = node_class->objects.nodes;

            /* For each processor */
            for( j = 0; j < PROC_COUNT; j++ ){
                nodes_on_proc = node_count_per_proc[j];
                py_ProcNodes = py_nodes_by_proc[j];

                proc_coords = (float*) malloc(nodes_on_proc * dims * sizeof(float));
                float_pointer_from_pyobject( py_ProcNodes, dims*nodes_on_proc, proc_coords);

                /* Loop over all nodes and coordinates */
                for( k = 0; k < nodes_on_proc; k++ ){
                    node_idx = proc_node_to_global_idx[j][k];
                    out_idx = node_idx * dims;
                    in_idx = k * dims;
                    node_coords[out_idx] = proc_coords[in_idx];
                    node_coords[out_idx+1] = proc_coords[in_idx+1];
                    node_coords[out_idx+2] = proc_coords[in_idx+2];
                }

                free( proc_coords );
            }

            node_class->labels_max = -1;
            node_class->labels_min = MAXINT;

            node_class->labels = NEW_N( MO_class_labels, obj_qty, "Class Labels" );
            //if ( node_class->labels == NULL )
            //    popup_fatal( "Unable to allocate labels on class load" );

            node_class->labels_index = NEW_N( int, obj_qty, "Class Labels index" );
            //if ( node_class->labels_index == NULL )
            //    popup_fatal( "Unable to allocate labels index on class load" );

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

            if ( block_range_ptr )
                free( block_range_ptr );
            if ( node_labels )
                free( node_labels );

            /* Add node class to list of element classes */
            //p_lh = p_md->classes_by_sclass + M_NODE;
            //mo_classes = (MO_class_data**) p_lh->list;
            //mo_classes = (void*) RENEW_N( MO_class_data*, p_lh->qty+1, "Extend Node sclass array");
            //mo_classes[p_lh->qty] = node_class;
            //p_lh->qty++;
            //p_lh->list = (void*) mo_classes; 

        }
        /* No nodal class, create one... */
        else{
            //popup_dialog( WARNING_POPUP, "Node object class not found \"node\"." );

            //// No nodal or element data, so create a fake nodal class
            //node_class->mesh_id = i;
            //griz_str_dup( &node_class->short_name, " " );
            //griz_str_dup( &node_class->long_name, " " );
            //node_class->superclass = G_UNIT;
            //node_class->elem_class_index = -1;
            //node_class->qty = 0;

            //// Allocate the data buffer for I/O and result derivation
            //node_class->data_buffer = NEW_N( float, 10000, "Class data buffer" );
        }

        /* Keep a reference to node geometry handy. */
        p_md->node_geom = node_class;

        /* ----------------------------------------------------------------------------------------- */

        /* For each element class load labels, connectivity, parts, materials, etc. */
        htable_get_data( p_ht, (void***) &htable_mo_classes, &elem_class_qty);

        /* Initialize arrays for storing element counts per processor and processor to global index mapping */
        elem_count_per_proc = (int**) malloc( elem_class_qty * sizeof(int*) );
        proc_elem_to_global_idx = (int***) malloc( elem_class_qty * sizeof(int**) );
        for( j = 0; j < elem_class_qty; j++ ){
            elem_count_per_proc[j] = (int*) malloc( PROC_COUNT * sizeof(int) );
            proc_elem_to_global_idx[j] = (int**) malloc( PROC_COUNT * sizeof(int*) );
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

            /* Allocate data buffer for I/O and result derivation */
            elem_class->data_buffer = NEW_N( float, obj_qty, "Class data buffer" );
            //if ( elem_class->data_buffer == NULL )
            //    popup_fatal( "Unable to alloc data buffer on class load" );

            int * elem_labels = (int*) malloc( elem_class->qty * sizeof(int));
            elems_on_proc = 0;
            pos = 0;
            for( k = 0; k < PROC_COUNT; k++ ){
                /* Get Labels on processor */
                py_ProcLabels = PyDict_GetItemString( py_labels_by_proc[k], elem_class->short_name );
                if( py_ProcLabels == NULL ){
                    elem_count_per_proc[elem_class->elem_class_index][k] = 0;
                    proc_elem_to_global_idx[elem_class->elem_class_index][k] = NULL;
                }
                else{
                    elems_on_proc = PySequence_Length( py_ProcLabels );

                    /* Read in element labels */
                    integer_pointer_from_pyobject( py_ProcLabels, elems_on_proc, elem_labels+pos );

                    /* Update arrays */
                    elem_count_per_proc[elem_class->elem_class_index][k] = elems_on_proc;
                    proc_elem_to_global_idx[elem_class->elem_class_index][k] = (int*) malloc( elems_on_proc * sizeof(int) );

                    /* Map proc elems to global indices */
                    for( l = 0; l < elems_on_proc; l++ ){
                        proc_elem_to_global_idx[j][k][l] = pos++;
                    }
                }
            }

            /* Sort labels */
            //qsort( elem_labels, elem_class->qty, sizeof(int), cmpfunc );

            /* Create per processor mapping for elem labels to index in global element list */
            /* MAYBE
            for( k = 0; k < PROC_COUNT; k++ ){
                if( elem_count_per_proc[k] > 0 ){
                    py_ProcLabels = PyDict_GetItemString( py_labels_by_proc[k], elem_class->short_name );
                    elems_on_proc = elem_count_per_proc[j][k];

                    integer_pointer_from_pyobject( py_ProcLabels, elems_on_proc, proc_elem_to_global_idx[j][k] );

                    int target;
                    for( l = 0; l < elems_on_proc; l++ ){
                        target = proc_elem_to_global_idx[j][k][l];
                        proc_elem_to_global_idx[j][k][l] = bin_search_index( target, elem_labels, elem_class->qty );
                    }
                }
            }
            */

            if( elem_class->superclass == M_MESH || elem_class->superclass == M_MAT || elem_class->superclass == M_UNIT ){

                elem_class->simple_start = elem_labels[0];
                elem_class->simple_stop = elem_labels[elem_class->qty-1];
                /* GRIZ
                if ( superclass == M_MAT ){
                    p_matd = NEW_N( Material_data, obj_qty, "Material data array" );
                    elem_class->objects.material = p_matd;
                    gen_material_data( elem_class, p_matd );
                }

                p_lh = p_md->classes_by_sclass + superclass;
                mo_classes = (MO_class_data **) p_lh->list;
                mo_classes = (void *) RENEW_N( MO_class_data *, mo_classes, p_lh->qty, 1, "Extend classes_by_sclass array");
                mo_classes[p_lh->qty] = elem_class;
                p_lh->qty++;
                p_lh->list = (void *) mo_classes;
                */
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

                int * elem_nodes = (int*) malloc( qty_conns * sizeof(int) );


                /* Load in nodal connectivity, material, part number, etc. for each element */
                for( k = 0; k < PROC_COUNT; k++ ){
                    if( elem_count_per_proc[j][k] > 0 ){
                        elems_on_proc = elem_count_per_proc[j][k];

                        /* Load connectivity */
                        py_ProcConns = PyDict_GetItemString( py_conns_by_proc[k], elem_class->short_name );
                        /* Load materials and parts */
                        py_ProcMats = PyDict_GetItemString( py_mats_by_proc[k], elem_class->short_name );
                        py_ProcParts = PyDict_GetItemString( py_parts_by_proc[k], elem_class->short_name );

                        int * mats = (int*) malloc( elems_on_proc * sizeof(int) );
                        integer_pointer_from_pyobject( py_ProcMats, elems_on_proc, mats );
                        int * parts = (int*) malloc( elems_on_proc * sizeof(int) );
                        integer_pointer_from_pyobject( py_ProcParts, elems_on_proc, parts );


                        int in_idx, out_idx;
                        PyObject * py_ConnSlice;
                        for( l = 0; l < elems_on_proc; l++ ){
                            /* Get elements nodal connectivity */
                            py_ConnSlice = PySequence_GetItem( py_ProcConns, l );
                            integer_pointer_from_pyobject( py_ConnSlice, qty_conns, elem_nodes );
                            /* Elem nodes is local id's for nodes need to convert to global labels */
                            for( m = 0; m < qty_conns; m++ ){
                                elem_nodes[m] = proc_node_to_global_idx[k][elem_nodes[m]];
                            } 
                            out_idx = proc_elem_to_global_idx[j][k][l] * qty_conns;
                            for( m = 0; m < qty_conns; m++ ){
                                p_ed->nodes[out_idx+m] = elem_nodes[m];
                            } 

                            out_idx = proc_elem_to_global_idx[j][k][l];
                            /* Get elements material and number */
                            p_ed->mat[out_idx] = mats[l];
                            p_ed->part[out_idx] = parts[l];
                        }

                        free( mats );
                        free( parts );
                    }
                }

                free( elem_nodes );

                elem_class->labels_max = -1;
                elem_class->labels_min = MAXINT;
                elem_class->labels_found = TRUE;

                elem_class->labels = NEW_N( MO_class_labels, obj_qty, "Class Labels" );
                //if ( elem_class->labels == NULL )
                //    popup_fatal( "Unable to allocate labels on class load" );

                elem_class->labels_index = NEW_N( int, obj_qty, "Class Labels index" );
                //if ( elem_class->labels_index == NULL )
                //    popup_fatal( "Unable to allocate labels index on class load" );


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

                /* Construct blocking data for node labels */
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

                if( block_range_ptr )
                    free(block_range_ptr);

                /* Element superclass-specific actions. */
                /*
                switch ( superclass )
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
                    check_degen_tris( p_mocd );
                    if ( p_mocd->objects.elems->has_degen )
                        popup_dialog( INFO_POPUP, "%s\n(class \"%s\").",
                                        "Degenerate tri element(s) detected",
                                        short_name );
                    break;
                default:
                    // do nothing
                    ;
                }
                */

                /* Update Mesh_data classes_by_sclass list for elements superclass */
                /* GRIZ
                p_lh = p_md->classes_by_sclass + superclass;
                mo_classes = (MO_class_data **) p_lh->list;
                mo_classes = (void *) RENEW_N( MO_class_data *, mo_classes, p_lh->qty, 1,
                                                "Extend classes_by_sclass array" );
                mo_classes[p_lh->qty] = elem_class;
                p_lh->qty++;
                p_lh->list = (void *) mo_classes;
                */

                /* Free elem labels */
                if( elem_labels != NULL )
                    free( elem_labels );
            }
        }

        /* Update the Mesh_data struct with element class info */
        p_md->elem_class_qty = elem_class_qty;

        /* Initialize mapping from processor to global */
        analy->index_map = NEW( ProcessorToGlobalMap, "Proc to Global index map");
        analy->index_map->node_map = proc_node_to_global_idx;
        analy->index_map->node_count = node_count_per_proc;
        analy->index_map->elem_map = proc_elem_to_global_idx;
        analy->index_map->elem_count = elem_count_per_proc;
    } // for ( i = 0 i < mesh_qty; i++ )

    /* Pass back address of geometry hash table array. */
    *p_mtable = mesh_array;
    *p_mesh_qty = mesh_qty;

    /* GRIZ
    p_analysis = get_analy_ptr();
    status = get_hex_volumes( dbid, p_analysis );
    p_analysis->dimension = dims;
    */

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
    free( py_mats_by_proc );
    free( py_parts_by_proc );

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
Bool_type
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

    for( i = 0; i < PROC_COUNT; i++ )
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
Bool_type
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

    mesh_id = 0;

    /* Single call to mili reader to get all the data we need */
    py_CallData = call_mili_reader_function_noargs( GET_ST_DESCRIPTORS_CALL );

    /* Get out the individual parts of the data */
    int * subrecord_count_per_proc = (int*) malloc( PROC_COUNT * sizeof(int));
    PyObject **py_subrecord_list_by_proc = (PyObject**) malloc(PROC_COUNT * sizeof(PyObject*));
    py_ESList = NEW_N( PyObject*, PROC_COUNT, "PyObject Element set List" );

    for( i = 0; i < PROC_COUNT; i++ ){
        py_Struct = PySequence_GetItem( py_CallData, i );
        if( py_Struct == NULL ){
            fprintf(stderr, "GrizStDescriptorsCallData is NULL for processor %d\n", i);
            return GRIZ_FAIL;
        }

        /* If processor 0, get number of state record formats. */
        if( i == 0 ){
            srec_qty = get_pyobject_attribute_as_int( py_Struct, "srec_fmt_qty" );
        }

        /* Store data for each processor */
        py_ESList[i] = get_pyobject_attribute( py_Struct, "element_sets");
        py_subrecord_list_by_proc[i] = get_pyobject_attribute( py_Struct, "subrecords" );
        subrecord_count_per_proc[i] = PySequence_Length( py_subrecord_list_by_proc[i] );
    }

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

    /* Get all subrecord from each processor and combine in hash table */
    total_svar_qty = 0;
    for( i = 0; i < PROC_COUNT; i++ ){
        py_Subrecords = py_subrecord_list_by_proc[i];
        subrec_qty = subrecord_count_per_proc[i];
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
                    fprintf( stderr, "mili_reader_get_st_descriptors call get_subrecord_from_pyobject(%s)\n", name);
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
                        fprintf( stderr, "mili_reader_get_st_descriptors call create_st_variable_from_pyobject(%s)\n", svar_name);
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
                                    fprintf( stderr, "mili_reader_get_st_descriptors call create_st_variable_from_pyobject(%s)\n", p_svar->components[l]);
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

            elem_blocks = NEW_N(int, block_qty, "Subrecord ordinal blocks");
            integer_pointer_from_pyobject( py_OrdinalBlocks, block_qty * 2, elem_blocks );

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
                    elem_gids[gidx++] = analy->index_map->elem_map[elem_class_index][i][l];
                }

                block_range_ptr = compute_label_blocking_data( elem_gids, obj_qty, &gid_block_qty );
                if( block_range_ptr == NULL ){
                    fprintf( stderr,
                             "mili_reader_get_st_descriptors call compute_label_blocking processor %d, class '%s'",
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
                    p_subr->mo_blocks[gidx++] = block_range_ptr[l*2];
                    p_subr->mo_blocks[gidx++] = block_range_ptr[l*2+1];
                }
                p_subr->qty_blocks += gid_block_qty;

            }

            free( elem_blocks );
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
        //create_subrec_node_list( node_work_array, mesh_node_qty, p_subrecs + i );

        /* Create ident array if indexing is required */
        class_size = p_mocd->qty;

        if( p_subr->qty_objects != class_size || (class_size > 1 && (p_subr->qty_blocks != 1 || p_subr->mo_blocks[0] != 1 )) )
        {
            p_subrecs[i].object_ids = NEW_N( int, p_subr->qty_objects, "Subrec ident map array" );
            //blocks_to_list( p_subr->qty_blocks, p_subr->mo_blocks, p_subrecs[i].object_ids, TRUE );
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
            //create_primal_result( p_mesh, srec_id, i, p_subrecs+i, p_primal_ht, srec_qty,
            //                      svar_names[k], p_sv_ht, analy );

            // Check for old shell stresses
            if(!strcmp(svar_names[j], "stress_in") || !strcmp(svar_names[j], "stress_mid") || !strcmp(svar_names[j], "stress_out"))
                analy->old_shell_stresses = TRUE;
            
            if ( nodal )
            {
                if ( strcmp( svar_names[j], "nodpos" ) == 0 )
                {
                    //if ( p_sro[srec_id].node_pos_subrec != -1 )
                    //    popup_dialog( WARNING_POPUP, "Multiple \"node\" position subrecs." );

                    p_sro[srec_id].node_pos_subrec = subrec_index;
                    analy->stateDB = FALSE;
                    if( mesh_node_qty == p_subr->qty_objects )
                    {
                        /* This is a state database and not a time history database. */
                        analy->stateDB = TRUE;
                    }

                    /* Note if data is double precision. */
                    //if ( p_svar->num_type == M_FLOAT8 )
                    //    p_mesh->double_precision_nodpos = TRUE;
                }
                else if ( strcmp( svar_names[j], "nodvel" ) == 0 )
                {
                    //if ( p_sro[srec_id].node_vel_subrec != -1 )
                    //    popup_dialog( WARNING_POPUP, "Multiple \"node\" velocity subrecs." );
                    p_sro[srec_id].node_vel_subrec = i;
                }
            }
        }
    }

    free( node_work_array );

    /* Return subrecord tree, state variable hash table and primal result hash table. */
    analy->srec_tree = p_sro;
    analy->qty_srec_fmt = srec_qty;
    /* GRIZ
    analy->st_var_table = p_sv_ht;
    analy->primal_result = p_primal_ht;
    */


#ifdef MILI_READER_TIMING
    end = clock();
    elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("[mili_reader_get_st_descriptors] elapsed = %f\n", elapsed);
#endif
    return OK;
}


/* TAG( combined_nodpos )
 *
 * Combine the nodal positions from multiple processors
 */
Bool_type
combine_nodpos( Analysis *analy, PyObject **py_ProcNodePositions, float * out_buffer ){
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

    dims = analy->dimension;

    for( i = 0; i < PROC_COUNT; i++ ){
        py_ProcNodes = py_ProcNodePositions[i];

        py_Nodes = PyDict_GetItemString( py_ProcNodes, nodpos );
        py_Nodes = PyDict_GetItemString( py_Nodes, data );
        py_Nodes = PyDict_GetItemString( py_Nodes, node );

        /* Single state, so remove state array */
        py_Nodes = PySequence_GetItem( py_Nodes, 0 );

        node_qty = analy->index_map->node_count[i];
        
        /* Get coordinates for each element */
        for( j = 0; j < node_qty; j++ ){
            py_Coords = PySequence_GetItem( py_Nodes, j );
            out_idx = analy->index_map->node_map[i][j] * dims;

            float_pointer_from_pyobject( py_Coords, dims, out_buffer+out_idx);
        }
    }
    return OK;
}

/* TAG( mili_reader_get_state )
 *
 * Move to a particular state in the Mili database and update
 * nodal positions for the mesh.
 */
Bool_type
mili_reader_get_state( Analysis *analy, int state_no, State2 *p_st, State2 **pp_new_st, int *state_qty )
{
    int i, j;
    int rval;
    int st_qty, st;
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

    PyObject * py_AllStateData;
    PyObject * py_ProcStateData;
    PyObject * py_StateArg;
    PyObject * py_FuncName;
    PyObject ** py_ProcNodePositions;

#ifdef MILI_READER_TIMING
    clock_t start, end;
    double elapsed;
    start = clock();
#endif

    /* Single Query to Mili Reader to get:
     *
     *  1. Number of states
     *  2. The Time associated with the requested state
     *  3. The query results for Nodal positions for this state
     *  4. The SAND flags (This can wait for now)
     *  5. Spheral data (This can wait for now)
     */

    /* Query the number of states in the database */
    PyObject* py_StateQty = PySequence_GetItem(call_mili_reader_function_noargs( STATE_QTY ), 0);
    if( py_StateQty == NULL ){
        fprintf( stderr, "py_StateQty is NULL\n" );
        return GRIZ_FAIL;
    }
    st_qty = PyLong_AsLong( py_StateQty );

    /* Get the number of dimensions */
    dims = analy->dimension;

    /* Pass back number of states */
    if ( state_qty != NULL )
        *state_qty = st_qty;
    
    if ( st_qty == 0 || analy->qty_srec_fmt == 0 )
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
        //popup_dialog( WARNING_POPUP, "Get-state request for nonexistent state." );
        *pp_new_st = p_st;
        return GRIZ_FAIL;
    }
    else
        st = state_no + 1;
    
    /* Retrieve state data from mili reader */
    py_FuncName = string_to_pyobject( GET_STATE_CALL );
    py_StateArg = PyLong_FromLong( st );
    py_AllStateData = PyObject_CallMethodObjArgs( mili_db, py_FuncName, py_StateArg, NULL );
    if( py_AllStateData == NULL ){
        fprintf( stderr, "py_AllStateData is NULL\n" );
        return GRIZ_FAIL;
    }

    /* Get data for each processor */
    py_ProcNodePositions = NEW_N( PyObject*, PROC_COUNT, "Processor Node Positions" );
    for( i = 0; i < PROC_COUNT; i++ ){
        py_ProcStateData = PySequence_GetItem( py_AllStateData, i );

        /* If processor 0, get the states time. */
        if( i == 0 )
            st_time = get_pyobject_attribute_as_int( py_ProcStateData, "state_time");

        /* Get nodal positions for each processor */
        py_ProcNodePositions[i] = get_pyobject_attribute( py_ProcStateData, "nodpos" );
    }
    
    /* Should only be a single state record format */
    srec_id = 0;

    /* Get subrecord tree, subrecords, and mesh data */
    p_sro = analy->srec_tree + srec_id;
    p_subrecs = p_sro->subrecs;
    mesh_id = p_subrecs[0].p_object_class->mesh_id;
    p_md = analy->mesh_table + mesh_id;

    /* Update or create State2 struct. */
    p_st = mk_state2( analy, p_sro, dims, srec_id, st_qty, p_st );
    p_st->state_no = state_no;

    /* Store the state time */
    p_st->time = st_time;

    /* Read node position array if it exists, re-ordering if necessary. */
    if ( TRUE || analy->stateDB )
    {
        rval = combine_nodpos( analy, py_ProcNodePositions, p_st->nodes.nodes );
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

    /* Handle spheral */

    /* If nodal positions weren't part of state data, get from geometry. */
    if ( !analy->stateDB )
    {
        p_st->nodes = p_md->node_geom->objects;
        p_st->position_constant = TRUE;
    }

    /* mili_db_get_state has code to read "partpos" svar. Do I need to support that??? */

    /* Return new state */
    *pp_new_st = p_st;

#ifdef MILI_READER_TIMING
    end = clock();
    elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("[mili_reader_get_state] elapsed = %f\n", elapsed);
#endif
    return OK;
}

#endif