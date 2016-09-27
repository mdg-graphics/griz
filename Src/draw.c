/* $Id$ */
/*
 * draw.c - Routines for drawing meshes.
 *
 *      Donald J. Dovey
 *      Lawrence Livermore National Laboratory
 *      Oct 23 1991
 *
 ************************************************************************
 * Modifications:
 *
 *  I. R. Corey - Sept 27, 2004: Increased length of label from 30 to 64.
 *
 *  I. R. Corey - Dec  28, 2004: Added new function draw_free_nodes. This
 *                function is invoked using the on free_nodes command.
 *                See SRC#: 286
 *
 *  I. R. Corey - Feb 10, 2005: Added a new free_node option to enable
 *                volume disable scaling of the nodes.
 *
 *  I. R. Corey - Aug 24, 2005: Fixed error with free nodes displaying
 *                for hidden materials.
 *                See SRC#: 337
 *
 *  I. R. Corey - Aug 31, 2005: Added state to the time info message..
 *                See SRC#: 315
 *
 *  I. R. Corey - Feb 02, 2006: Add a new capability to display meshless,
 *                particle-based results.
 *                See SRC#367.
 *
 *  I. R. Corey - Mar 20, 2006: Fixed a problem with display of meshless
 *                results when other object types are also present (bricks
 *                and shells).
 *                See SRC#376.
 *
 *  I. R. Corey - Mar 23, 2006: Removed edge lines from reflected boundaries
 *                created with the sym command.
 *                See SRC#377.
 *
 *  I. R. Corey - Oct 17, 2006: Fixed problem with polygon edges not rendering
 *                correctly on some platforms (Xwin-32). Lines need to be offset
 *                from faces.
 *                See SCR#420.
 *
 *  I. R. Corey - Oct 26, 2006: Added new feature for enable/disable vis/invis
 *                that allow selection of an object type and material for that
 *                object.
 *                See SRC#421.
 *
 *  I. R. Corey - Nov 3, 2006: Modified particle scaling so that particles no
 *                longer scale with window.
 *                See SRC#422.
 *
 *  I. R. Corey - Dec 16, 2007: Added full path name to plot window and remove
 *                from titlebar.
 *                See SRC#431.
 *
 *  I. R. Corey - Feb 12, 2007: Added wireframe viewing capability and added
 *                logic to fix degenerate polygon faces.
 *                See SRC#437.
 *
 *  I. R. Corey - Aug 09, 2007:	Added date/time string to render window.
 *                See SRC#478.
 *
 *  I. R. Corey - Oct 18, 2007: Fixed problem with free for multiple
 *                reflection planes - in draw_line moved free outside
 *                of reflection plane loop.
 *                See SRC#497.
 *
 *  I. R. Corey - Nov 08, 2007:	Add Node/Element labeling.
 *                See SRC#418 (Mili) and 504 (Griz)
 *
 *  I. R. Corey - Dec 17, 2007:	Modified the edge bias control for win32
 *                using a new flag to identify running on win32.
 *                See SRC#509
 *
 *  I. R. Corey - Feb 21, 2008: Fixed a problem with reflecting beams
 *                in draw_line().
 *                See SRC#520
 *
 *  I. R. Corey - June 08, 2009: Fixed a problem with checking of snd_array
 *                that has not been allocated.
 *                See SRC#608
 *
 *  I. R. Corey - Dec 28, 2009: Fixed several problems releated to drawing
 *                ML particles.
 *                See SRC#648
 *
 *  I. R. Corey - March 24th, 2011: Added option to allow viewing of
 *                inactive elements.
 *                See TeamForge#14543
 *
 *  I. R. Corey - Sept 8th, 2011: Add particles to list of objects that
 *                can be identified with numclass.
 *
 *  I. R. Corey - July 26th, 2012: Added capability to plot a Modal
 *                database from Diablo.
 *                See TeamForge#18395 & 18396
 *
 *  I. R. Corey - September 20th, 2012: Fixed problem with rendering free
 *                nodes.
 *
 *  I. R. Corey - November 7th, 2012:  Fixed color scale problem with
 *                rendering damage result.
 *
 *  I. R. Corey - February 11th, 2013: For DBC particles, only render results
 *                for subrecords that are valid for the particle class.
 *                See TeamForge#19430
 *
 *  I. R. Corey - March 18th, 2013: Fixed problem with lighting being disabled
 *                when any particles are rendered.
 *                See TeamForge#19627
 *
 *  I. R. Corey - March 25th, 2013: Added check for particle nodes being
 *                enabled when drawing a result hilite. If pn enabled then
 *                use the pn result.
 *                See TeamForge#19627
 *
 *  I. R. Corey - March 27th, 2013: Corrected bbox usage comment for
 *                Warning: Near/far planes shouldn't go behind the eyepoint.
 *                See TeamForge#19778
 *
 *
 ************************************************************************
 */

#include "griz_config.h"
#include <stdlib.h>
#include <math.h>

#ifdef sgi
#define TAN tanf
#define ATAN2 atan2f
#endif
#ifdef hpux
#define TAN tanf
#define ATAN2 atan2f
#endif
#ifdef __alpha
#define TAN tanf
#define ATAN2 atan2f
#endif
#ifndef TAN
#define TAN tan
#define ATAN2 atan2
#endif

#include "viewer.h"
#include "draw.h"
#include "image.h"
#include "mdg.h"

#ifdef xSERIAL_BATCHx
#include <GL/gl_mangle.h>
#endif

/*
 * Data/time variables
 */
#include <time.h>
struct tm *curr_datetime;
time_t tm;

/*
* Added November 05, 2007: IRC
*
* Returns the index for for the specified label.
*
*/

int  label_compare( const int *key, const MO_class_labels *label );
int  get_class_label_index( MO_class_data *class, int label_num );
void dump_class_labels( MO_class_data *class );
void populate_result(int superclass, char map[], int size, MO_class_data * p_class, Analysis * analy);

int grayel = 0;
const float grayval = 0.7;
char particle_name[256];
/*
 * This if not defined statement takes care of the problem at WES
 * where the titles, time, etc are not displayed
 */
#ifndef O2000x
#include "hershey.h"
#endif

#ifdef PNG_SUPPORT
#include <png.h>        /* libpng header; includes zlib.h and setjmp.h */
#endif


#ifdef JPEG_SUPPORT
#include "jpeglib.h"
#include "jerror.h"
#endif

#define SCL_MAX ((float) CMAP_SIZE - 2.01)

typedef struct _Render_poly_obj
{
    struct _Render_poly_obj *next;
    struct _Render_poly_obj *prev;
    int cnt;
    float pts[4][3];
    float norm[4][3];
    float cols[4][4];
    float vals[4];
} Render_poly_obj;


/* Object which holds all the display information. */
View_win_obj *v_win = NULL;

/* Cache some constants that are needed for "good" interpolation. */
static Transf_mat cur_view_mat;
static float proj_param_x;
static float proj_param_y;

/* Give this file scope for efficiency during "good" interpolation. */
static Bool_type colorflag=FALSE;

/* Particle radius. */
static GLdouble particle_radius = 0.025;

/* Vector lengths in pixels for drawing vector plots. */
#define VEC_2D_LENGTH 20.0
#define VEC_3D_LENGTH 100.0

/* Vertical spacing coeefficient and offset for text. */
#define LINE_SPACING_COEFF (1.18)
#define TOP_OFFSET (1.0)

/* Some annotation strings. */
char *strain_label[] =
{
    "infinitesimal", "Green-Lagrange", "Almansi", "Rate"
};
char *ref_surf_label[] = { "middle", "inner", "outer" };
char *ref_frame_label[] = { "global", "local" };

GLfloat black[] = {0.0, 0.0, 0.0};
GLfloat white[] = {1.0, 1.0, 1.0};
GLfloat red[] = {1.0, 0.0, 0.0};
GLfloat green[] = {0.0, 1.0, 0.0};
GLfloat blue[] = {0.0, 0.0, 1.0};
GLfloat magenta[] = {1.0, 0.0, 1.0};
GLfloat cyan[] = {0.0, 1.0, 1.0};
GLfloat yellow[] = {1.0, 1.0, 0.0};

#define MATERIAL_COLOR_CNT 20
GLfloat material_colors[MATERIAL_COLOR_CNT][3] =
{
    {0.933, 0.867, 0.51},    /* Light goldenrod */
    {0.529, 0.808, 0.922},   /* Sky blue */
    {0.725, 0.333, 0.827},   /* Medium orchid */
    {0.804, 0.361, 0.361},   /* Indian red */
    {0.4, 0.804, 0.667},     /* Medium aquamarine */
    {0.416, 0.353, 0.804},   /* Slate blue */
    {1.0, 0.647, 0.0},       /* Orange */
    {0.545, 0.271, 0.075},   /* Brown */
    {0.118, 0.565, 1.0},     /* Dodger blue */
    {0.573, 0.545, 0.341},   /* SeaGreen */
    {0.961, 0.871, 0.702},   /* Wheat */
    {0.824, 0.706, 0.549},   /* Tan */
    {1.0, 0.412, 0.706},     /* Hot pink */
    {0.627, 0.125, 0.941},   /* Purple */
    {0.282, 0.82, 0.8},      /* Medium turqoise */
    {0.902, 0.157, 0.157},   /* Red */
    {0.157, 0.902, 0.157},   /* Green */
    {0.902, 0.902, 0.157},   /* Yellow */
    {0.157, 0.902, 0.902},   /* Cyan */
    {0.902, 0.157, 0.902},   /* Magenta */
};

#define SURFACE_COLOR_CNT 20
static GLfloat surface_colors[SURFACE_COLOR_CNT][3] =
{
    {0.902, 0.157, 0.902},   /* Magenta */
    {0.157, 0.902, 0.902},   /* Cyan */
    {0.902, 0.902, 0.157},   /* Yellow */
    {0.157, 0.902, 0.157},   /* Green */
    {0.902, 0.157, 0.157},   /* Red */
    {0.282, 0.82, 0.8},      /* Medium turqoise */
    {0.627, 0.125, 0.941},   /* Purple */
    {1.0, 0.412, 0.706},     /* Hot pink */
    {0.824, 0.706, 0.549},   /* Tan */
    {0.961, 0.871, 0.702},   /* Wheat */
    {0.573, 0.545, 0.341},   /* SeaGreen */
    {0.118, 0.565, 1.0},     /* Dodger blue */
    {0.545, 0.271, 0.075},   /* Brown */
    {1.0, 0.647, 0.0},       /* Orange */
    {0.416, 0.353, 0.804},   /* Slate blue */
    {0.4, 0.804, 0.667},     /* Medium aquamarine */
    {0.804, 0.361, 0.361},   /* Indian red */
    {0.725, 0.333, 0.827},   /* Medium orchid */
    {0.529, 0.808, 0.922},   /* Sky blue */
    {0.933, 0.867, 0.51},    /* Light goldenrod */
};

GLfloat material_colors_gs[MATERIAL_COLOR_CNT][3] =
{
    {0.88, 0.88, 0.88},    /* Very light */
    {0.85, 0.85, 0.85},
    {0.80, 0.80, 0.80},
    {0.75, 0.75, 0.75},
    {0.70, 0.70, 0.70},
    {0.68, 0.68, 0.68},
    {0.62, 0.62, 0.62},
    {0.60, 0.60, 0.60},
    {0.58, 0.58, 0.58},
    {0.55, 0.55, 0.55},
    {0.50, 0.50, 0.50},
    {0.45, 0.45, 0.45},
    {0.40, 0.40, 0.40},
    {0.35, 0.35, 0.35},
    {0.32, 0.32, 0.32},
    {0.30, 0.30, 0.30},
    {0.28, 0.28, 0.28},
    {0.25, 0.25, 0.25},
    {0.22, 0.22, 0.22},
    {0.18, 0.18, 0.18},   /* Very dark */
};


#define PARTICLE_COLOR_CNT 1
static GLfloat particle_colors[PARTICLE_COLOR_CNT][3] =
{
    {0.933, 0.867, 0.51}     /* Light goldenrod */
};

/* Default colors for various functions. */
GLfloat bg_default[] = { 1.0, 1.0, 1.0 };       /* White. */
GLfloat fg_default[] = { 0.0, 0.0, 0.0 };       /* Black. */
GLfloat text_default[] = { 0.0, 0.0, 0.0 };     /* Black. */
GLfloat mesh_default[] = { 0.0, 0.0, 0.0 };     /* Black. */
GLfloat edge_default[] = { 0.0, 0.0, 0.0 };     /* Black. */
GLfloat con_default[] = { 1.0, 0.0, 1.0 };      /* Magenta. */
GLfloat hilite_default[] = { 1.0, 0.0, 0.0 };   /* Red. */
GLfloat select_default[] = { 0.0, 1.0, 0.0 };   /* Green. */
GLfloat rmin_default[] = { 1.0, 0.0, 1.0 };     /* Magenta. */
GLfloat rmax_default[] = { 1.0, 0.0, 0.0 };     /* Red. */
GLfloat vec_default[] = { 0.0, 0.0, 0.0 };      /* Black. */
GLfloat vechd_default[] = { 0.0, 0.0, 0.0 };    /* Black. */

/* Local routines. */
static void look_rot_mat( float [3], float [3], float [3], Transf_mat * );
static void create_color_prop_arrays( Color_property *, int );
static void define_one_color_property( Color_property *, int, GLfloat [][3], int );
static void color_lookup( float [4], float, float, float, float, int, Bool_type, Bool_type );
static void surf_color_lookup( float [4], float, float, float, float, int, Bool_type );

static int check_for_tri_face( int, int, MO_class_data *, float [4][3],
                               int [4] );
static void get_min_max( Analysis *, Bool_type, float *, float * );

/*static void draw_grid( Analysis * );
static void draw_grid_2d( Analysis * ); */ 

static void draw_hexs( Bool_type, Bool_type, Bool_type, MO_class_data *,
                       Analysis * );

static void draw_pyramids( Bool_type, Bool_type, Bool_type, MO_class_data *,
                           Analysis * );

static void draw_wedges( Bool_type, Bool_type, Bool_type, MO_class_data *,
                         Analysis * );

static void draw_tets( Bool_type, Bool_type, Bool_type, MO_class_data *,
                       Analysis * );

static void draw_quads_2d( Bool_type, Bool_type, Bool_type, MO_class_data *,
                           Analysis * );
static void draw_quads_3d( Bool_type, Bool_type, Bool_type, MO_class_data *,
                           Analysis * );

static void draw_tris_2d( Bool_type, Bool_type, Bool_type, MO_class_data *,
                          Analysis * );
static void draw_tris_3d( Bool_type, Bool_type, Bool_type, MO_class_data *,
                          Analysis * );

static void draw_beams_2d( Bool_type, Bool_type, Bool_type, MO_class_data *, Analysis * );
static void draw_beams_3d( Bool_type, Bool_type, Bool_type, MO_class_data *, Analysis * );

static void draw_truss_2d( Bool_type, Bool_type, MO_class_data *, Analysis * );
static void draw_truss_3d( Bool_type, Bool_type, MO_class_data *, Analysis * );

static void draw_surfaces_2d( Color_property *, Bool_type, Bool_type, Bool_type,MO_class_data *, Analysis * );
static void draw_surfaces_3d( Color_property *, Bool_type, Bool_type, Bool_type,MO_class_data *, Analysis * );

static void draw_particles_2d( MO_class_data *, Analysis * );
static void draw_particles_3d( MO_class_data *, Analysis * );

static void draw_edges_2d( Analysis * );
static void draw_edges_3d( Analysis * );

static void draw_nodes_2d_3d( MO_class_data *, Analysis * );

static void draw_hilite( Bool_type, MO_class_data *, int, Analysis * );
static void draw_class_numbers( Analysis * );
static void draw_bbox( float [2][3] );
static void draw_extern_polys( Analysis * );
static void draw_ref_polys( Analysis * );
static void draw_vec_result_2d( Analysis * );
static void draw_vec_result_3d( Analysis * );
static void draw_vecs_2d( Vector_pt_obj *, float [3][3], float, float, float,
                          float, float, Mesh_data *,  Analysis * );
static void draw_vecs_3d( Vector_pt_obj *, float, float, float, float, float,
                          Mesh_data *,  Analysis * );
static void draw_node_vec_2d_3d( Analysis * );
static void load_vec_result( Analysis *, float *[3], float *, float * );
/*
static void find_front_faces();
static void draw_reg_carpet();
static void draw_vol_carpet();
static void draw_shell_carpet();
static void draw_ref_carpet();
*/
static void draw_traces( Analysis * );
static void draw_trace_list( Trace_pt_obj *, Analysis * );
static void draw_sphere( float [3], float, int);
static void draw_sphere_GL( float [3], float, int);
static void draw_poly( int, float [4][3], float [4][3], float [4][4],
                       float [4], int, Mesh_data *, Analysis *, Bool_type );
static void draw_edged_poly( int, float [4][3], float [4][3], float [4][4],
                             float [4], int, Mesh_data *, Analysis *, Bool_type );
static void draw_edged_wireframe_poly( int, float [4][3], float [4][3], float [4][4],
                                       float [4], int, Mesh_data *, Analysis * );
static void draw_plain_poly( int, float [4][3], float [4][3], float [4][4],
                             float [4], int, Mesh_data *, Analysis *, Bool_type );
static float length_2d( float [2], float [2] );
static void scan_poly( int, float [4][3], float [4][3], float [4], int,
                       Mesh_data *, Analysis * );
static void draw_poly_2d( int, float [4][3], float [4][4], float [4], int,
                          Mesh_data *, Analysis * );
static void draw_line( int, float *, int, Mesh_data *, Analysis *, Bool_type, Bool_type * );
static void draw_3d_text( float [3], char *, Bool_type );

static void draw_free_nodes( Analysis * );
void   check_for_free_nodes( Analysis * );

static void draw_foreground( Analysis * );
static void get_verts_of_bbox( float [2][3], float [8][3] );
static void hvec_copy( float [4], float [4] );
static void antialias_lines( Bool_type, Bool_type );
static void begin_draw_poly( Analysis * );
static void end_draw_poly( Analysis * );
static void linear_variable_scale( float, float, int, float *, float *,
                                   float *, int * );
extern void log_variable_scale( float, float, int, float *, float *,
                                float *, float *, int * );
extern void log_scale_data_shift( float val, float data_minimum, float data_maximum,
                                  float *new_val, float *new_data_minimum,
                                  float *new_data_maximum, float *data_shift,
                                  float *data_mult);
static void puthexpix( FILE *, unsigned char );
static void epilogue( FILE *, int );
static void prologue( FILE *, int, int, int, float, float, float, float );
static double griz_round( double, double );
static double tfloor( double, double );
static double machine_tolerance( void );
static void draw_locref( Analysis * );
static void draw_locref_hex( Analysis * );
static void memory_to_screen();

static void change_current_color_property( Color_property *, int );
void update_current_color_property( Color_property *, Material_property_type );

static float
get_free_node_result( Analysis *, MO_class_data *, int, Bool_type *, Bool_type * );

static float
get_ml_result( Analysis *, MO_class_data *, int, Bool_type * );

int
get_particle_node_num( Analysis *analy, MO_class_data *p_mo_class,
                       int elem_num );
int
get_particle_from_node( Analysis *analy, MO_class_data *p_mo_class,
                        int node_num );
void
select_meshless( Analysis *analy, Mesh_data *p_mesh,
                 MO_class_data *p_ml_class, int near_node,
                 int *p_near );

void DrawCone(float len, float base_diam, float cols[4][4], int res );

/*
 * SECTION_TAG( Grid window )
 */


/*****************************************************************
 * TAG( init_mesh_window )
 *
 * Initialize GL window for viewing the mesh.
 */
void
init_mesh_window( Analysis *analy )
{
    int i;
    int qty, mtl_qty;
    int sum,  rval;
    Htable_entry *p_hte;
    MO_class_data **p_classes;
    /**
        This creates a positional light as opposed to a directional light
        static GLfloat light1_position[] = { 10.0, 10.0, -25.0, 1.0 };

        This creates a directional light
        static GLfloat light1_position[] = { 10.0, 10.0, -25.0, 0.0 };
    **/
    static GLfloat lmodel_ambient[] = { 0.2, 0.2, 0.2, 1.0 };
    static GLfloat light_ambient[] = { 0.0, 0.0, 0.0, 1.0 };
    static GLfloat light_diffuse[] = { 1.0, 1.0, 1.0, 1.0 };
    static GLfloat light_specular[] = { 1.0, 1.0, 1.0, 1.0 };
    static GLfloat light0_position[] = { 0.0, 0.0, 10.0, 0.0 };
    static GLfloat light1_position[] = { 10.0, 10.0, -25.0, 1.0 };

    /* Initialize the view window data structures. */
    v_win = NEW( View_win_obj, "View window" );

    /* Lighting is on for 3D and off for 2D. */
    if ( analy->dimension == 3 )
        v_win->lighting = TRUE;
    else
        v_win->lighting = FALSE;

    /* Orthographic or perspective depends on dimension. */
    if ( analy->dimension == 2 )
        set_orthographic( TRUE );
    else
        set_orthographic( FALSE );

    /* Initialize view. */
    VEC_SET( v_win->look_from, 0.0, 0.0, 20.0 );
    VEC_SET( v_win->look_at, 0.0, 0.0, 0.0 );
    VEC_SET( v_win->look_up, 0.0, 1.0, 0.0 );
    v_win->cam_angle = RAD_TO_DEG( 2.0*ATAN2(1.0, v_win->look_from[2]) );
    v_win->near = v_win->look_from[2] - 2.0;
    v_win->far = v_win->look_from[2] + 2.0;
    v_win->aspect_correct = 1.0;

    mat_copy( &v_win->rot_mat, &ident_matrix );
    VEC_SET( v_win->trans, 0.0, 0.0, 0.0 );
    VEC_SET( v_win->scale, 1.0, 1.0, 1.0 );

    if ( MESH( analy ).node_geom->qty>0 )
    {
        bbox_nodes( analy, analy->bbox_source, TRUE, analy->bbox[0],
                    analy->bbox[1] );
        set_view_to_bbox( analy->bbox[0], analy->bbox[1], analy->dimension );
    }

    /* A little initial rotation to make things look nice. */
    if ( analy->dimension == 3 )
    {
        inc_mesh_rot( 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, DEG_TO_RAD( 30 ) );
        inc_mesh_rot( 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, DEG_TO_RAD( 15 ) );
    }

    /* Set the default colors and load the result colormap. */
    set_color( "bg", NULL );
    set_color( "fg", NULL );
    set_color( "con", NULL );
    set_color( "hilite", NULL );
    set_color( "select", NULL );
    set_color( "rmin", NULL );
    set_color( "rmax", NULL );
    set_color( "vec", NULL );
    set_color( "vechd", NULL );
    hot_cold_colormap();

    /* Select and load a Hershey font. */
    hfont( "futura.l" );

    /* Initialize video title array. */
    for ( i = 0; i < 4; i++ )
        v_win->vid_title[i][0] = '\0';

    /* Turn on Z-buffering for 3D but not for 2D. */
    if ( analy->dimension == 3 )
    {
        glEnable( GL_DEPTH_TEST );
        glDepthFunc( GL_LEQUAL );
    }
    else
        glDisable( GL_DEPTH_TEST );

    /* Normalize normal vectors after transformation.  This is crucial. */
    if ( analy->dimension == 3 )
        glEnable( GL_NORMALIZE );

    /* Set the background color for clears. */
    glClearDepth( 1.0 );
    glClearColor( v_win->backgrnd_color[0],
                  v_win->backgrnd_color[1],
                  v_win->backgrnd_color[2], 0.0 );
    glClearStencil( 0 );

    /* Set alignment for glReadPixels during raster output operations. */
    glPixelStorei( GL_PACK_ALIGNMENT, (GLint) 1 );

    /* Ask for best line anti-aliasing. */
    glHint( GL_LINE_SMOOTH_HINT, GL_NICEST );

    /* Set up for perspective or orthographic viewing. */
    set_mesh_view();

    /* Initialize the GL model view matrix stack. */
    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();

    /* Define material properties for each material. */
    /* If particle class present, add an extra material. */
    rval = htable_search( analy->mesh_table[0].class_table, particle_cname,
                          FIND_ENTRY, &p_hte );
    mtl_qty = ( rval == OK ) ? analy->max_mesh_mat_qty + 1
              : analy->max_mesh_mat_qty;
    create_color_prop_arrays( &v_win->mesh_materials, mtl_qty );
    define_color_properties( &v_win->mesh_materials, NULL, mtl_qty,
                             material_colors, MATERIAL_COLOR_CNT );

    if ( (qty = MESH_P( analy )->classes_by_sclass[G_SURFACE].qty) > 0 )
    {
        p_classes = (MO_class_data **) MESH_P( analy )->classes_by_sclass[G_SURFACE].list;

        sum = 0;
        for ( i = 0; i < qty; i++ )
            sum += p_classes[i]->qty;

        create_color_prop_arrays( &v_win->surfaces, sum );
        define_color_properties( &v_win->surfaces, NULL, sum, surface_colors, SURFACE_COLOR_CNT );
    }

    if ( (qty = MESH_P( analy )->classes_by_sclass[G_PARTICLE].qty) > 0 )
    {
        p_classes = (MO_class_data **) MESH_P( analy )->classes_by_sclass[G_PARTICLE].list;

        sum = 0;
        for ( i = 0; i < qty; i++ )
            sum += p_classes[i]->qty;

        create_color_prop_arrays( &v_win->particles, sum );
        define_color_properties( &v_win->particles, NULL, sum, particle_colors, PARTICLE_COLOR_CNT );
    }
    /* Init lighting even if 2D db so a subsequent 3D db load will be OK. */

    /* Default to one-sided lighting (this is the OpenGL default). */
    glLightModeli( GL_LIGHT_MODEL_TWO_SIDE, 0 );

    /* Enable lighting. */
    glEnable( GL_LIGHTING );

    /* Set the global ambient light level. */
    glLightModelfv( GL_LIGHT_MODEL_AMBIENT, lmodel_ambient );

    /* Enable first light. */
    glLightfv( GL_LIGHT0, GL_AMBIENT, light_ambient );
    glLightfv( GL_LIGHT0, GL_DIFFUSE, light_diffuse );
    glLightfv( GL_LIGHT0, GL_SPECULAR, light_specular );
    glLightfv( GL_LIGHT0, GL_POSITION, light0_position );
    glEnable( GL_LIGHT0 );
    v_win->light_active[0] = TRUE;
    hvec_copy( v_win->light_pos[0], light0_position );

    /* Enable second light. */
    glLightfv( GL_LIGHT1, GL_AMBIENT, light_ambient );
    glLightfv( GL_LIGHT1, GL_DIFFUSE, light_diffuse );
    glLightfv( GL_LIGHT1, GL_SPECULAR, light_specular );
    glLightfv( GL_LIGHT1, GL_POSITION, light1_position );
    glEnable( GL_LIGHT1 );
    v_win->light_active[1] = TRUE;
    hvec_copy( v_win->light_pos[1], light1_position );

    /* Initialize the current material (the number is arbitrary). */
    v_win->mesh_materials.current_index = 0;
    glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT, v_win->mesh_materials.ambient[0] );
    glMaterialfv( GL_FRONT_AND_BACK, GL_DIFFUSE, v_win->mesh_materials.diffuse[0] );
    glMaterialfv( GL_FRONT_AND_BACK, GL_SPECULAR, v_win->mesh_materials.specular[0] );
    glMaterialfv( GL_FRONT_AND_BACK, GL_EMISSION, v_win->mesh_materials.emission[0] );
    glMaterialf( GL_FRONT_AND_BACK, GL_SHININESS, v_win->mesh_materials.shininess[0] );

    /*** BEGIN HERE---

       1.  v_win->mtl_ambient (and friends) to:  v_win->mesh_materials.ambient[0], as required
       2.  update "change_current_material"
       3.  check existing particle rendering logic

    ***/

    /*** Is this correct?? LAS ***/
    v_win->current_color_property = &v_win->mesh_materials;


    /* By default, leave lighting disabled. */
    glDisable( GL_LIGHTING );

    /* Clear the mesh window. */
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT |
             GL_STENCIL_BUFFER_BIT );
}


/*****************************************************************
 * TAG( reset_mesh_window )
 *
 * Reset the GL window for viewing the mesh.  This routine is
 * called after a new analysis has been loaded.  The routine
 * resets the bbox scaling and translation.
 */
void
reset_mesh_window( Analysis *analy )
{
    int mtl_qty, rval;
    Htable_entry *p_hte;

    /* Lighting is on for 3D and off for 2D. */
    if ( analy->dimension == 3 )
        v_win->lighting = TRUE;
    else
        v_win->lighting = FALSE;

    /* Orthographic or perspective depends on dimension. */
    if ( analy->dimension == 2 )
        set_orthographic( TRUE );
    else
        set_orthographic( FALSE );

    bbox_nodes( analy, analy->bbox_source, TRUE, analy->bbox[0],
                analy->bbox[1] );
    set_view_to_bbox( analy->bbox[0], analy->bbox[1], analy->dimension );

    /* Turn on Z-buffering for 3D but not for 2D. */
    if ( analy->dimension == 3 )
    {
        glEnable( GL_DEPTH_TEST );
        glDepthFunc( GL_LEQUAL );
    }
    else
        glDisable( GL_DEPTH_TEST );

    /* Normalize normal vectors after transformation.  This is crucial. */
    /**/
    /* Strictly 3D requirement? */
    if ( analy->dimension == 3 )
        glEnable( GL_NORMALIZE );
    else
        glDisable( GL_NORMALIZE );

    /* In case data has changed from 3D to 2D or vice-versa. */
    set_mesh_view();

    /* Set default color properties for any new materials. */
    rval = htable_search( analy->mesh_table[0].class_table, particle_cname,
                          FIND_ENTRY, &p_hte );
    mtl_qty = ( rval == OK ) ? analy->max_mesh_mat_qty + 1
              : analy->max_mesh_mat_qty;
    extend_color_prop_arrays( &v_win->mesh_materials, mtl_qty,
                              material_colors, MATERIAL_COLOR_CNT );

    /* Set default color properties for any new surfaces. */
    extend_color_prop_arrays( &v_win->surfaces, analy->surface_qty,
                              surface_colors, SURFACE_COLOR_CNT );

}


/*****************************************************************
 * TAG( create_color_prop_arrays )
 *
 * Allocate arrays for color property data.
 */
static void
create_color_prop_arrays( Color_property *p_color_property, int size )
{
    float gslevel=.40, gslevel_increment=0.;
    int i;
    if ( p_color_property->ambient != NULL )
        free( p_color_property->ambient );
    if ( p_color_property->diffuse != NULL )
        free( p_color_property->diffuse );
    if ( p_color_property->specular != NULL )
        free( p_color_property->specular );
    if ( p_color_property->shininess != NULL )
        free( p_color_property->shininess );
    if ( p_color_property->emission != NULL )
        free( p_color_property->emission );
    if ( p_color_property->gslevel != NULL )
        free( p_color_property->gslevel );

    p_color_property->property_array_size = size;
    p_color_property->property_array_size = size;

    p_color_property->ambient   = NEW_N( GLVec4, size, "Color property ambient array" );
    p_color_property->diffuse   = NEW_N( GLVec4, size, "Color property diffuse array" );
    p_color_property->specular  = NEW_N( GLVec4, size, "Color property specular array" );
    p_color_property->emission  = NEW_N( GLVec4, size, "Color property emission array" );
    p_color_property->shininess = NEW_N( GLfloat, size, "Color property shininess array" );
    p_color_property->gslevel   = NEW_N( GLfloat, size, "Color property gslevel array" );

    /* Put gray level in range of (.40 - .90/size) */
    gslevel_increment = (0.90 - 0.40)/size;
    for ( i=0;
            i<size;
            i++ )
    {
        p_color_property->gslevel[i] = gslevel;
        gslevel+=gslevel_increment;
    }
}


/*****************************************************************
 * TAG( extend_color_prop_arrays )
 *
 * Extend existing arrays for color property data if new size
 * is larger than existing size of arrays.
 */
void
extend_color_prop_arrays( Color_property *p_cp, int new_size, GLfloat rgb_array[][ 3 ], int rgb_size )
{
    int i, add, old;
    int *tmp_color_property_indicies;

    if ( new_size <= p_cp->property_array_size )
        return;

    old = p_cp->property_array_size;
    add = new_size - old;

    /* Extend the arrays. */

    p_cp->ambient   = RENEW_N( GLVec4, p_cp->ambient, old, add,
                               "Material ambient array" );
    p_cp->diffuse   = RENEW_N( GLVec4, p_cp->diffuse, old, add,
                               "Material diffuse array" );
    p_cp->specular  = RENEW_N( GLVec4, p_cp->specular, old, add,
                               "Material specular array" );
    p_cp->emission  = RENEW_N( GLVec4, p_cp->emission, old, add,
                               "Material emission array" );
    p_cp->shininess = RENEW_N( GLfloat, p_cp->shininess, old, add,
                               "Material shininess array" );


    p_cp->property_array_size = new_size;

    /* Initialize the new materials. */
    tmp_color_property_indicies = NEW_N( int, add, "Temp matl nums" );
    for ( i = 0; i < add; i++ )
        tmp_color_property_indicies[i] = old + i;

    define_color_properties( p_cp, tmp_color_property_indicies, add, rgb_array, rgb_size );

    free( tmp_color_property_indicies );
}


/*****************************************************************
 * TAG( get_pick_ray )
 *
 * From pick coordinates on the screen, get a ray that goes through
 * the scene.  Returns the eyepoint as the point on the ray and
 * the ray direction is in the view direction.
 */
void
get_pick_ray( int sx, int sy, float pt[3], float dir[3] )
{
    Transf_mat view_mat;
    float ppt[3], eyept[3];
    float midx, midy, px, py, cx, cy;

    /*
     * Get the pick point and eye point in view coordinates.
     * Inverse transform the pick point and eye point into world
     * coordinates.  The ray is then computed from these two points.
     */
    midx = 0.5*v_win->vp_width + v_win->win_x;
    midy = 0.5*v_win->vp_height + v_win->win_y;
    px = (sx - midx)/(0.5 * v_win->vp_width);
    py = (sy - midy)/(0.5 * v_win->vp_height);
    get_foreground_window( &ppt[2], &cx, &cy );
    ppt[0] = px * cx;
    ppt[1] = py * cy;

    /* In Motif, Y direction of window is reverse that of GL. */
    ppt[1] = -ppt[1];

    /* Get the inverse viewing transformation matrix. */
    inv_view_transf_mat( &view_mat );

    if ( v_win->orthographic )
    {
        VEC_COPY( eyept, ppt );
        eyept[2] = 0.0;
    }
    else
    {
        VEC_SET( eyept, 0.0, 0.0, 0.0 );
    }

    point_transform( pt, ppt, &view_mat );
    VEC_COPY( ppt, pt );
    point_transform( pt, eyept, &view_mat );
    VEC_SUB( dir, ppt, pt );
}


/*****************************************************************
 * TAG( get_pick_pt )
 *
 * Transform a pick point in screen coordinates to world coordinates.
 */
void
get_pick_pt( int sx, int sy, float pick_pt[3] )
{
    Transf_mat view_mat;
    float ppt[3], pt[3];
    float midx, midy, px, py, cx, cy;

    /*
     * Get the pick point and eye point in view coordinates.
     * Inverse transform the pick point.
     */
    midx = 0.5*v_win->vp_width + v_win->win_x;
    midy = 0.5*v_win->vp_height + v_win->win_y;
    px = (sx - midx)/(0.5 * v_win->vp_width);
    py = (sy - midy)/(0.5 * v_win->vp_height);
    get_foreground_window( &ppt[2], &cx, &cy );
    ppt[0] = px * cx;
    ppt[1] = py * cy;

    /* In Motif, Y direction of window is reverse that of GL. */
    ppt[1] = -ppt[1];

    /* Get the inverse viewing transformation matrix. */
    inv_view_transf_mat( &view_mat );

    point_transform( pt, ppt, &view_mat );
    VEC_COPY( pick_pt, pt );
}


/*
 * SECTION_TAG( View control )
 */


/*****************************************************************
 * TAG( set_orthographic )
 *
 * Select a perspective or orthographic view.
 */
void
set_orthographic( Bool_type val )
{
    v_win->orthographic = val;
}


/*****************************************************************
 * TAG( set_mesh_view )
 *
 * Set up the perspective or orthographic view for the mesh window,
 * taking into account the dimensions of the current viewport and
 * the camera angle.
 */
void
set_mesh_view( void )
{
    GLint param[4];
    GLdouble xp, yp, cp;
    float aspect;
    int i;
    for(i = 0; i < 4; i++)
    {
        param[i] = 0;
    }

    /* Get the current viewport location and size. */
    glGetIntegerv( GL_VIEWPORT, param );
    if(param[2] == 0)
    {
        param[2] = 600;
    }
    if(param[3] == 0)
    {
        param[3] = 600;
    }
    v_win->win_x = param[0];
    v_win->win_y = param[1];
    v_win->vp_width = param[2];
    v_win->vp_height = param[3];

    /* Correct the aspect ratio for NTSC video, if requested. */
    aspect = v_win->aspect_correct * v_win->vp_width / (float)v_win->vp_height;

    /* Get the window dimensions at the near plane. */
    if ( v_win->orthographic )
        cp = 1.0;
    else
        cp = v_win->near * TAN( DEG_TO_RAD( 0.5*v_win->cam_angle ) );

    if ( aspect >= 1.0 )
    {
        xp = cp * aspect;
        yp = cp;
    }
    else
    {
        xp = cp;
        yp = cp / aspect;
    }

    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();

    /*
     * We could special case 2D data, but choose not to for now.
     */
    if ( v_win->orthographic )
        glOrtho( -xp, xp, -yp, yp, v_win->near, v_win->far );
    else
        glFrustum( -xp, xp, -yp, yp, v_win->near, v_win->far );

    glMatrixMode( GL_MODELVIEW );
}


/************************************************************
 * TAG( aspect_correct )
 *
 * If video is TRUE, correct the screen aspect ratio for NTSC.
 * We found NTSC to have an aspect ratio of about x/y = 0.92.
 * Note that this means the image on the workstation screen won't
 * look right.  If video is false, the aspect is set to 1.0.
 */
void
aspect_correct( Bool_type video )
{

    if ( video )
        v_win->aspect_correct = 0.92;
    else
        v_win->aspect_correct = 1.0;

    set_mesh_view();
}


/*****************************************************************
 * TAG( set_camera_angle )
 *
 * Set the view angle.
 */
void
set_camera_angle( float angle, Analysis *analy )
{
    v_win->cam_angle = angle;
    adjust_near_far( analy );
}


/*****************************************************************
 * TAG( set_view_to_bbox )
 *
 * This routine sets a scaling and translation of the mesh
 * so that it fits the viewport well.  Routine is passed the
 * bounding box of the mesh and the dimension of the mesh.
 * Routine can be called repeatedly.
 */
void
set_view_to_bbox( float bb_lo[3], float bb_hi[3], int view_dimen )
{
    int i;
    float dimx, dimy, dimz, dim;

    for ( i = 0; i < 3; i++ )
        v_win->bbox_trans[i] = - ( bb_lo[i] + 0.5*(bb_hi[i]-bb_lo[i]) );

    dimx = 0.5*(bb_hi[0] - bb_lo[0]);
    dimy = 0.5*(bb_hi[1] - bb_lo[1]);
    dimz = 0.5*(bb_hi[2] - bb_lo[2]);
    dim = MAX( MAX( dimx, dimy ), dimz );

    if ( view_dimen == 3 )
        v_win->bbox_scale = 0.75 / dim;
    else  /* view_dimen == 2 */
        v_win->bbox_scale = 0.9 / dim;
}


/*****************************************************************
 * TAG( inc_mesh_rot )
 *
 * Increment the current mesh rotation.  Rotates about an axis
 * specified by a point and vector.  Angle is in radians.
 */
void
inc_mesh_rot( float px, float py, float pz, float vx, float vy, float vz,
              float angle )
{
    float pt[3];
    float vec[3];

    pt[0] = px;
    pt[1] = py;
    pt[2] = pz;
    vec[0] = vx;
    vec[1] = vy;
    vec[2] = vz;

    mat_rotate( &v_win->rot_mat, pt, vec, angle );
}


/*****************************************************************
 * TAG( inc_mesh_trans )
 *
 * Increment the current mesh translation.
 */
void
inc_mesh_trans( float tx, float ty, float tz )
{
    v_win->trans[0] += tx;
    v_win->trans[1] += ty;
    v_win->trans[2] += tz;
}


/*****************************************************************
 * TAG( set_mesh_scale )
 *
 * Set the current mesh scaling.
 */
void
set_mesh_scale( float scale_x, float scale_y, float scale_z )
{
    v_win->scale[0] = scale_x;
    v_win->scale[1] = scale_y;
    v_win->scale[2] = scale_z;
}


/*****************************************************************
 * TAG( get_mesh_scale )
 *
 * Returns the current mesh scaling factor.
 */
void
get_mesh_scale( float scale[3] )
{
    VEC_COPY( scale, v_win->scale );
}


/*****************************************************************
 * TAG( reset_mesh_transform )
 *
 * Resets the mesh view to the default.
 */
void
reset_mesh_transform( void )
{
    int i;

    mat_copy( &v_win->rot_mat, &ident_matrix );
    for ( i = 0; i < 3; i++ )
        v_win->trans[i] = 0.0;
    VEC_SET( v_win->scale, 1.0, 1.0, 1.0 );
}


/*****************************************************************
 * TAG( look_from look_at look_up )
 *
 * Modify the mesh view so that the user is looking from one
 * point toward a second point, with the specified up vector in
 * the Y direction.
 */
void
look_from( float pt[3] )
{
    VEC_COPY( v_win->look_from, pt );
}

void
look_at( float pt[3] )
{
    VEC_COPY( v_win->look_at, pt );
}

void
look_up( float vec[3] )
{
    VEC_COPY( v_win->look_up, vec );
}


/************************************************************
 * TAG( move_look_from )
 *
 * Move the look_from point along the specified axis (X == 0, etc.).
 */
void
move_look_from( int axis, float incr )
{
    v_win->look_from[axis] += incr;
}


/************************************************************
 * TAG( move_look_at )
 *
 * Move the look_at point along the specified axis (X == 0, etc.).
 */
void
move_look_at( int axis, float incr )
{
    v_win->look_at[axis] += incr;
}


/*****************************************************************
 * TAG( look_rot_mat )
 *
 * Creates a rotation matrix for the look from/at transform.
 */
static void
look_rot_mat( float from_pt[3], float at_pt[3], float up_vec[3],
              Transf_mat *rot_mat )
{
    float v1[3], v2[3], v3[3];
    int i;

    /* Do the look at. */
    mat_copy( rot_mat, &ident_matrix );
    VEC_SUB( v3, from_pt, at_pt );
    vec_norm( v3 );
    VEC_CROSS( v1, up_vec, v3 );
    vec_norm( v1 );
    VEC_CROSS( v2, v3, v1 );
    vec_norm( v2 );

    for ( i = 0; i < 3; i++ )
    {
        rot_mat->mat[i][0] = v1[i];
        rot_mat->mat[i][1] = v2[i];
        rot_mat->mat[i][2] = v3[i];
    }
}


/*****************************************************************
 * TAG( adjust_near_far )
 *
 * Resets the near and far planes of the viewing fulcrum in order
 * to avoid cutting off any of the model.
 */
void
adjust_near_far( Analysis *analy )
{
    Transf_mat view_trans;
    float verts[8][3];
    float t_pt[3];
    float low_z, high_z;
    int i;
    float new_near, new_far;

    /*
     * Set up the viewing transformation matrix.
     */
    view_transf_mat( &view_trans );

    /*
     * Get the corners of the bounding box.
     */
    get_verts_of_bbox( analy->bbox, verts );

    /*
     * Transform corners and get near and far Z values.
     */
    point_transform( t_pt, verts[0], &view_trans );
    low_z = t_pt[2];
    high_z = t_pt[2];

    for ( i = 1; i < 8; i++ )
    {
        point_transform( t_pt, verts[i], &view_trans );
        if ( t_pt[2] < low_z )
            low_z = t_pt[2];
        else if ( t_pt[2] > high_z )
            high_z = t_pt[2];
    }

    new_near = -high_z - 1.0;
    new_far = -low_z + 1.0;

    /*
     * Near/far planes shouldn't go behind the eyepoint.
     */
    if ( !analy->adjust_NF_disable_warning && (new_near < 0.0 || new_far < new_near) )
        popup_dialog( WARNING_POPUP, "%s\n%s\n%s\n%s\n%s",
                      "Unable to reset near/far planes.",
                      "If using material invisibility, try \"bbox xmin ymin zmin xmax ymax zmax\"",
                      "to minimize the mesh bounding box then \"rnf\" again OR",
                      "use \"info\" to see previous near/far planes and",
                      "adjust manually using \"near\" and \"far\"." );

    if ( new_near < 0.0 || new_far < new_near )
    {
        set_mesh_view();

    }
    else
    {
        v_win->near = new_near;
        v_win->far = new_far;
        set_mesh_view();
    }
}


/*****************************************************************
 * TAG( set_near )
 *
 * Set the position of the near cutting plane.
 */
void
set_near( float near_dist )
{
    v_win->near = near_dist;
}


/*****************************************************************
 * TAG( set_far )
 *
 * Set the position of the far cutting plane.
 */
void
set_far( float far_dist )
{
    v_win->far = far_dist;
}


/************************************************************
 * TAG( center_view )
 *
 * Center the view about a point specified in one of three
 * ways.
 */
void
center_view( Analysis *analy )
{
    Transf_mat view_trans;
    float ctr_pt[3], opt[3];
    float *pt, *factors;
    int hi_num;
    int node;
    int i, j;
    int node_qty;
    int dim;
    float *coords, *ocoords;
    MO_class_data *p_mo_class;
    int *connects;
    Bool_type dscale=FALSE;

    pt = analy->view_center;
    dim = analy->dimension;
    coords = analy->state_p->nodes.nodes;

    dscale = analy->displace_exag;
    if ( analy->displace_exag )
    {
        dscale = TRUE;
        ocoords = analy->cur_ref_state_data;
        factors = analy->displace_scale;
    }
    else
        dscale = FALSE;

    /* Get the centering point. */
    switch ( analy->center_view )
    {
    case NO_CENTER:
        VEC_SET( pt, 0.0, 0.0, 0.0 );
        break;

    case HILITE:
        p_mo_class = analy->hilite_class;
        if ( p_mo_class == NULL )
            return;
        VEC_SET( pt, 0.0, 0.0, 0.0 );
        if ( dscale )
            VEC_SET( opt, 0.0, 0.0, 0.0 );
        hi_num = analy->hilite_num;

        /* Get the centering point. */
        switch ( p_mo_class->superclass )
        {
        case G_NODE:
            for ( i = 0; i < dim; i++ )
            {
                pt[i] = coords[hi_num * dim + i];
                if ( dscale )
                    opt[i] = ocoords[hi_num * dim + i];
            }
            break;

        case G_TRUSS:
        case G_BEAM:
        case G_TRI:
        case G_QUAD:
        case G_TET:
        case G_PYRAMID:
        case G_WEDGE:
        case G_HEX:
        case G_SURFACE:
        case G_PARTICLE:
            connects = p_mo_class->objects.elems->nodes;
            node_qty = qty_connects[p_mo_class->superclass];
            for ( i = 0; i < node_qty; i++ )
                for ( j = 0; j < dim; j++ )
                {
                    node = connects[hi_num * node_qty + i];
                    pt[j] += coords[node * dim + j] / (float) node_qty;

                    if ( dscale )
                        opt[j] += ocoords[node * dim + j]
                                  / (float) node_qty;
                }
            break;
        default:
            return;
        }
        break;

    case NODE:
        node = analy->center_node;
        for ( i = 0; i < dim; i++ )
        {
            pt[i] = coords[node * dim + i];
            if ( dscale )
                opt[i] = ocoords[node * dim + i];
        }
        break;

    case POINT:
        /* Do nothing, coords already in pt (i.e., analy->view_center). */
        if ( dscale )
            VEC_SET( opt, 0.0, 0.0, 0.0 );
        break;

    default:
        return;
    }

    if ( dscale )
    {
        /* Scale the point's displacements. */
        for ( i = 0; i < 3; i++ )
            pt[i] = opt[i] + factors[i] * (pt[i] - opt[i]);
    }

    view_transf_mat( &view_trans );
    point_transform( ctr_pt, pt, &view_trans );
    inc_mesh_trans( -ctr_pt[0], -ctr_pt[1],
                    -((v_win->near + v_win->far) * 0.5 + ctr_pt[2]) );
}


/*****************************************************************
 * TAG( view_transf_mat )
 *
 * Create a view transformation matrix for the current view.
 */
void
view_transf_mat( Transf_mat *view_trans )
{
    Transf_mat look_rot;
    float scal;

    /*
     * Set up the viewing transformation matrix.
     */
    mat_copy( view_trans, &ident_matrix );
    mat_translate( view_trans, v_win->bbox_trans[0], v_win->bbox_trans[1],
                   v_win->bbox_trans[2] );
    scal = v_win->bbox_scale;
    mat_scale( view_trans, scal*v_win->scale[0], scal*v_win->scale[1],
               scal*v_win->scale[2] );
    mat_mul( view_trans, view_trans, &v_win->rot_mat );
    mat_translate( view_trans, v_win->trans[0], v_win->trans[1],
                   v_win->trans[2] );
    mat_translate( view_trans, -v_win->look_from[0], -v_win->look_from[1],
                   -v_win->look_from[2] );
    look_rot_mat( v_win->look_from, v_win->look_at, v_win->look_up,
                  &look_rot );
    mat_mul( view_trans, view_trans, &look_rot );
}


/*****************************************************************
 * TAG( inv_view_transf_mat )
 *
 * Create the inverse of the current view transformation matrix.
 */
void
inv_view_transf_mat( Transf_mat *inv_trans )
{
    Transf_mat rot_mat, look_rot;
    float scal;

    /*
     * Set up the inverse viewing transformation matrix.
     */
    mat_copy( inv_trans, &ident_matrix );
    look_rot_mat( v_win->look_from, v_win->look_at, v_win->look_up,
                  &look_rot );
    invert_rot_mat( &rot_mat, &look_rot );
    mat_mul( inv_trans, inv_trans, &rot_mat );
    mat_translate( inv_trans, v_win->look_from[0], v_win->look_from[1],
                   v_win->look_from[2] );
    mat_translate( inv_trans, -v_win->trans[0], -v_win->trans[1],
                   -v_win->trans[2] );
    invert_rot_mat( &rot_mat, &v_win->rot_mat );
    mat_mul( inv_trans, inv_trans, &rot_mat );
    scal = 1.0 / v_win->bbox_scale;
    mat_scale( inv_trans, scal/v_win->scale[0], scal/v_win->scale[1],
               scal/v_win->scale[2] );
    mat_translate( inv_trans, -v_win->bbox_trans[0], -v_win->bbox_trans[1],
                   -v_win->bbox_trans[2] );
}


/************************************************************
 * TAG( print_view )
 *
 * Print information about the current view matrix.
 */
void
print_view( void )
{
    float vec[3];
    int i;

    mat_to_angles( &v_win->rot_mat, vec );
    wrt_text( "Current view parameters:\n" );
    wrt_text( "    Rotation, X: %f  Y: %f  Z: %f\n",
              RAD_TO_DEG(vec[0]), RAD_TO_DEG(vec[1]), RAD_TO_DEG(vec[2]) );
    wrt_text( "    Translation, X: %f  Y: %f  Z: %f\n",
              v_win->trans[0], v_win->trans[1], v_win->trans[2] );
    wrt_text( "    Scale: %f %f %f\n", v_win->scale[0], v_win->scale[1],
              v_win->scale[2] );
    wrt_text( "    Near plane: %f\n", v_win->near );
    wrt_text( "    Far plane: %f\n", v_win->far );
    wrt_text( "    Look from position: %f %f %f\n", v_win->look_from[0],
              v_win->look_from[1], v_win->look_from[2] );
    wrt_text( "    Look at position: %f %f %f\n", v_win->look_at[0],
              v_win->look_at[1], v_win->look_at[2] );
    wrt_text( "    Look up position: %f %f %f\n", v_win->look_up[0],
              v_win->look_up[1], v_win->look_up[2] );
    wrt_text( "    View angle: %f\n", v_win->cam_angle );

    for ( i = 0; i < 6; i++ )
        if ( v_win->light_active[i] )
            wrt_text( "    Light %d is at point: %f, %f, %f, %f\n", i + 1,
                      v_win->light_pos[i][0], v_win->light_pos[i][1],
                      v_win->light_pos[i][2], v_win->light_pos[i][3] );

    wrt_text( "\n" );
}


/*
 * SECTION_TAG( Materials and colors )
 */


/*****************************************************************
 * TAG( define_color_properties )
 *
 * FORMERLY:  define_materials
 *
 * Defines gl color properties for a Color_property.
 */
extern void
define_color_properties( Color_property *p_color_property, int *indices, int qty,
                         GLfloat rgb_array[][ 3 ], int rgb_array_size )
{
    int i;

    if ( indices == NULL )
        /* Define color properties for all materials. */
        for ( i = 0; i < qty; i++ )
            define_one_color_property( p_color_property, i, rgb_array, rgb_array_size );
    else
        /* Define color properties only for the listed materials. */
        for ( i = 0; i < qty; i++ )
            define_one_color_property( p_color_property, indices[i],
                                       rgb_array, rgb_array_size );
}


/*****************************************************************
 * TAG( define_one_color_property )
 *
 * Defines gl material properties for one material in the mesh.
 */
static void
define_one_color_property( Color_property *p_color_property, int index,
                           GLfloat rgb_array[][3], int rgb_array_size )
{
    int idx, j;

    /*
     * Here are the OpenGL material default values.
     *
     *     GL_AMBIENT  0.2, 0.2, 0.2, 1.0
     *     GL_DIFFUSE  0.8, 0.8, 0.8, 1.0
     *     GL_AMBIENT_AND_DIFFUSE
     *     GL_SPECULAR  0.0, 0.0, 0.0, 1.0
     *     GL_SHININESS  0.0
     *     GL_EMISSION  0.0, 0.0, 0.0, 1.0
     */

    for ( j = 0; j < 3; j++ )
    {
        idx = index % rgb_array_size;
        p_color_property->ambient[index][j] = rgb_array[idx][j];
        p_color_property->diffuse[index][j] = rgb_array[idx][j];
        p_color_property->specular[index][j] = 0.0;
        p_color_property->emission[index][j] = 0.0;
    }

    p_color_property->ambient[index][3] = 1.0;
    p_color_property->diffuse[index][3] = 1.0;
    p_color_property->specular[index][3] = 1.0;
    p_color_property->emission[index][3] = 1.0;
    p_color_property->shininess[index] = 0.0;

    /* If update affects current OpenGL material, update it. */
    if ( index == p_color_property->current_index )
    {
        update_current_color_property( p_color_property, AMBIENT );
        update_current_color_property( p_color_property, DIFFUSE );
        update_current_color_property( p_color_property, SPECULAR );
        update_current_color_property( p_color_property, EMISSIVE );
        update_current_color_property( p_color_property, SHININESS );
    }
}


/*****************************************************************
 * TAG( change_current_color_property )
 *
 * REVISED VERSION:
 *
 * Change the currently selected color properties.
 */
void
change_current_color_property( Color_property *p_cp, int num )
{
    Color_property *p_local_cp;

    int current_index;


    /* */

    p_local_cp = v_win->current_color_property;
    current_index = p_local_cp->current_index;

    /*
     * no changes in any color properties
     */

    if ( p_cp == p_local_cp  && current_index == num )
        return;

    /*
     * Update only the properties that are different.
     */

    /*
     * NOTE:  ambient and diffuse properties are overridden by the current color, anyway

    if ( p_cp->ambient[ num ][ 0 ] != p_local_cp->ambient[ current_index ][ 0 ] ||
         p_cp->ambient[ num ][ 1 ] != p_local_cp->ambient[ current_index ][ 1 ] ||
         p_cp->ambient[ num ][ 2 ] != p_local_cp->ambient[ current_index ][ 2 ] )
    {
       glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT, p_cp->ambient[ num ] );
    }

    if ( p_cp->diffuse[ num ][ 0 ] != p_local_cp->diffuse[ current_index ][ 0 ] ||
         p_cp->diffuse[ num ][ 1 ] != p_local_cp->diffuse[ current_index ][ 1 ] ||
         p_cp->diffuse[ num ][ 2 ] != p_local_cp->diffuse[ current_index ][ 2 ] )
    {
       glMaterialfv( GL_FRONT_AND_BACK, GL_DIFFUSE, p_cp->diffuse[ num ] );
    }
    */


    if ( p_cp->specular[ num ][ 0 ] != p_local_cp->specular[ current_index ][ 0 ] ||
            p_cp->specular[ num ][ 1 ] != p_local_cp->specular[ current_index ][ 1 ] ||
            p_cp->specular[ num ][ 2 ] != p_local_cp->specular[ current_index ][ 2 ] )
    {
        glMaterialfv( GL_FRONT_AND_BACK, GL_SPECULAR, p_cp->specular[ num ] );
    }

    if ( p_cp->emission[ num ][ 0 ] != p_local_cp->emission[ current_index ][ 0 ] ||
            p_cp->emission[ num ][ 1 ] != p_local_cp->emission[ current_index ][ 1 ] ||
            p_cp->emission[ num ][ 2 ] != p_local_cp->emission[ current_index ][ 2 ] )
    {
        glMaterialfv( GL_FRONT_AND_BACK, GL_EMISSION, p_cp->emission[ num ] );
    }

    if ( p_cp->shininess[ num ] != p_local_cp->shininess[ current_index ]  )
    {
        glMaterialf( GL_FRONT_AND_BACK, GL_SHININESS, p_cp->shininess[ num ] );
    }


    if ( p_cp != p_local_cp )
        v_win->current_color_property = p_cp;

    p_cp->current_index = num;
}


/*****************************************************************
 * TAG( update_current_color_property )
 *
 * Update a property for the currently selected material.
 */
void
update_current_color_property( Color_property *p_color_property, Material_property_type prop )
{
    int index;

    index = p_color_property->current_index;

    switch( prop )
    {
    case AMBIENT:
        glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT,
                      p_color_property->ambient[index] );
        break;
    case DIFFUSE:
        glMaterialfv( GL_FRONT_AND_BACK, GL_DIFFUSE,
                      p_color_property->diffuse[index] );
        break;
    case SPECULAR:
        glMaterialfv( GL_FRONT_AND_BACK, GL_SPECULAR,
                      p_color_property->specular[index] );
        break;
    case EMISSIVE:
        glMaterialfv( GL_FRONT_AND_BACK, GL_EMISSION,
                      p_color_property->emission[index] );
        break;
    case SHININESS:
        glMaterialf( GL_FRONT_AND_BACK, GL_SHININESS,
                     p_color_property->shininess[index] );
        break;
    }
}


/************************************************************
 * TAG( set_material )
 *
 * Define a material's rendering properties.
 */
Redraw_mode_type
set_material( int token_cnt, char tokens[MAXTOKENS][TOKENLENGTH], int max_qty )
{
    int matl_num, i;
    Bool_type rval;

    /* Get the material number. */
    sscanf( tokens[1], "%d", &matl_num );

    if ( matl_num > max_qty )
        return NO_VISUAL_CHANGE;
    else
        matl_num--;

    /*    matl_num = matl_num % MAX_MATERIALS; */

    rval = FALSE;

    for ( i = 2; i < token_cnt; i++ )
    {
        if ( strcmp( tokens[i], "amb" ) == 0 )
        {
            sscanf( tokens[i+1], "%f", &v_win->mesh_materials.ambient[matl_num][0] );
            sscanf( tokens[i+2], "%f", &v_win->mesh_materials.ambient[matl_num][1] );
            sscanf( tokens[i+3], "%f", &v_win->mesh_materials.ambient[matl_num][2] );
            i += 3;

            /* Need to update GL if this material is currently in use. */
            if ( v_win->current_color_property == &v_win->mesh_materials  &&
                    matl_num == v_win->mesh_materials.current_index )
            {
                glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT,
                              v_win->mesh_materials.ambient[matl_num] );
                rval = TRUE;
            }
        }
        else if ( strcmp( tokens[i], "diff" ) == 0 )
        {
            sscanf( tokens[i+1], "%f", &v_win->mesh_materials.diffuse[matl_num][0] );
            sscanf( tokens[i+2], "%f", &v_win->mesh_materials.diffuse[matl_num][1] );
            sscanf( tokens[i+3], "%f", &v_win->mesh_materials.diffuse[matl_num][2] );
            i += 3;

            if ( v_win->current_color_property == &v_win->mesh_materials  &&
                    matl_num == v_win->mesh_materials.current_index )
            {
                glMaterialfv( GL_FRONT_AND_BACK, GL_DIFFUSE,
                              v_win->mesh_materials.diffuse[matl_num] );
                rval = TRUE;
            }

        }
        else if ( strcmp( tokens[i], "spec" ) == 0 )
        {
            sscanf( tokens[i+1], "%f", &v_win->mesh_materials.specular[matl_num][0] );
            sscanf( tokens[i+2], "%f", &v_win->mesh_materials.specular[matl_num][1] );
            sscanf( tokens[i+3], "%f", &v_win->mesh_materials.specular[matl_num][2] );
            i += 3;

            if ( v_win->current_color_property == &v_win->mesh_materials  &&
                    matl_num == v_win->mesh_materials.current_index )
            {
                glMaterialfv( GL_FRONT_AND_BACK, GL_SPECULAR,
                              v_win->mesh_materials.specular[matl_num] );
                rval = TRUE;
            }
        }
        else if ( strcmp( tokens[i], "emis" ) == 0 )
        {
            sscanf( tokens[i+1], "%f", &v_win->mesh_materials.emission[matl_num][0] );
            sscanf( tokens[i+2], "%f", &v_win->mesh_materials.emission[matl_num][1] );
            sscanf( tokens[i+3], "%f", &v_win->mesh_materials.emission[matl_num][2] );
            i += 3;

            if ( v_win->current_color_property == &v_win->mesh_materials  &&
                    matl_num == v_win->mesh_materials.current_index )
            {
                glMaterialfv( GL_FRONT_AND_BACK, GL_EMISSION,
                              v_win->mesh_materials.emission[matl_num] );
                rval = TRUE;
            }
        }
        else if ( strcmp( tokens[i], "shine" ) == 0 )
        {
            sscanf( tokens[i+1], "%f", &v_win->mesh_materials.shininess[matl_num] );
            i++;

            if ( v_win->current_color_property == &v_win->mesh_materials  &&
                    matl_num == v_win->mesh_materials.current_index )
            {
                glMaterialf( GL_FRONT_AND_BACK, GL_SHININESS,
                             v_win->mesh_materials.shininess[matl_num] );
                rval = TRUE;
            }
        }
        else if ( strcmp( tokens[i], "alpha" ) == 0 )
        {
            /* The alpha value at a vertex is just the material's diffuse
             * alpha value -- the other alphas are ignored by OpenGL.
             */
            sscanf( tokens[i+1], "%f", &v_win->mesh_materials.diffuse[matl_num][3] );
            i++;

            if ( v_win->current_color_property == &v_win->mesh_materials  &&
                    matl_num == v_win->mesh_materials.current_index )
            {
                glMaterialfv( GL_FRONT_AND_BACK, GL_DIFFUSE,
                              v_win->mesh_materials.diffuse[matl_num] );
                rval = TRUE;
            }
        }
    }

    return rval ? BINDING_MESH_VISUAL : NO_VISUAL_CHANGE;
}


/*****************************************************************
 * TAG( set_color )
 *
 * Sets a color for the mesh window.  The options are:
 *
 *     "fg" -- Color of lines and text.
 *     "bg" -- Background color.
 *     "text" -- Text color.
 *     "mesh" -- Element face edge line color.
 *     "edges" -- Mesh edge color.
 *     "con" -- Color of contour lines.
 *     "hilite" -- Color of node & element pick highlights.
 *     "select" -- Color of node & element selection highlights.
 *     "vec" -- Color of 3D vectors for vector command.
 *     "vechd" -- Color for vector heads in vector carpets.
 *     "rmin" -- Color of result minimum threshold.
 *     "rmax" -- Color of result maximum threshold.
 */
void
set_color( char *select, float *rgb )
{
    float *rgb_src;

    /* Sanity check. */
    if ( rgb != NULL )
    {
        if ( rgb[0] > 1.0 || rgb[1] > 1.0 || rgb[2] > 1.0 )
        {
            /* Assume the user accidentally entered an integer in [0, 255]. */
            rgb[0] = rgb[0] / 255.0;
            rgb[1] = rgb[1] / 255.0;
            rgb[2] = rgb[2] / 255.0;
            popup_dialog( INFO_POPUP,
                          "Colors should be in the range (0.0, 1.0)" );
        }
    }

    if ( strcmp( select, "fg" ) == 0 )
    {
        rgb_src = ( rgb != NULL ) ? rgb : fg_default;
        VEC_COPY( v_win->foregrnd_color, rgb_src );
    }
    else if ( strcmp( select, "bg" ) == 0 )
    {
        rgb_src = ( rgb != NULL ) ? rgb : bg_default;
        VEC_COPY( v_win->backgrnd_color, rgb_src );

        /* Set the background color for clears. */
        glClearColor( v_win->backgrnd_color[0],
                      v_win->backgrnd_color[1],
                      v_win->backgrnd_color[2], 0.0 );
    }
    else if ( strcmp( select, "text" ) == 0 )
    {
        rgb_src = ( rgb != NULL ) ? rgb : text_default;
        VEC_COPY( v_win->text_color, rgb_src );
    }
    else if ( strcmp( select, "mesh" ) == 0 )
    {
        rgb_src = ( rgb != NULL ) ? rgb : mesh_default;
        VEC_COPY( v_win->mesh_color, rgb_src );
    }
    else if ( strcmp( select, "edges" ) == 0 )
    {
        rgb_src = ( rgb != NULL ) ? rgb : edge_default;
        VEC_COPY( v_win->edge_color, rgb_src );
    }
    else if ( strcmp( select, "con" ) == 0 )
    {
        rgb_src = ( rgb != NULL ) ? rgb : con_default;
        VEC_COPY( v_win->contour_color, rgb_src );
    }
    else if ( strcmp( select, "hilite" ) == 0 )
    {
        rgb_src = ( rgb != NULL ) ? rgb : hilite_default;
        VEC_COPY( v_win->hilite_color, rgb_src );
    }
    else if ( strcmp( select, "select" ) == 0 )
    {
        rgb_src = ( rgb != NULL ) ? rgb : select_default;
        VEC_COPY( v_win->select_color, rgb_src );
    }
    else if ( strcmp( select, "vec" ) == 0 )
    {
        rgb_src = ( rgb != NULL ) ? rgb : vec_default;
        VEC_COPY( v_win->vector_color, rgb_src );
    }
    else if ( strcmp( select, "vechd" ) == 0 )
    {
        /* Color for vector heads in vector carpets. */
        rgb_src = ( rgb != NULL ) ? rgb : vechd_default;
        VEC_COPY( v_win->vector_hd_color, rgb_src );
    }
    else if ( strcmp( select, "rmin" ) == 0 )
    {
        /* Color for result minimum cutoff. */
        rgb_src = ( rgb != NULL ) ? rgb : rmin_default;
        VEC_COPY( v_win->rmin_color, rgb_src );
    }
    else if ( strcmp( select, "rmax" ) == 0 )
    {
        /* Color for result maximum cutoff. */
        rgb_src = ( rgb != NULL ) ? rgb : rmax_default;
        VEC_COPY( v_win->rmax_color, rgb_src );
    }
    else
        popup_dialog( INFO_POPUP, "Color selection \"%s\" not recognized",
                      select );
}


/************************************************************
 * TAG( color_lookup )
 *
 * Lookup a color in the colormap, and store it in the first
 * argument.  If file scope variable colorflag is FALSE, the
 * routine just returns the color property color.
 */
static void
color_lookup( float col[4], float val, float result_min, float result_max,
              float threshold, int index, Bool_type logscale, Bool_type greyscale )
{
    int idx, idx_new;
    float scl_max;


    float new_val,     new_result_min,     new_result_max,     new_result_shift;
    float new_result_mult;
    float new_val_log, new_result_min_log, new_result_max_log;

    Color_property *color_map;
    Analysis *p_analy=NULL;
    p_analy = get_analy_ptr();

    color_map = &v_win->mesh_materials;

    if ( greyscale )
    {
        col[0] = (float ) color_map->gslevel[index];
        col[1] = (float ) color_map->gslevel[index];
        col[2] = (float ) color_map->gslevel[index];
        return;
    }

    /* Get alpha from the material. */

    col[3] = color_map->diffuse[index][3];

    /* If not doing colormapping or if result is near zero, use the
     * default material color.  Otherwise, do the color table lookup.
     */
    if ( !colorflag )
    {
        VEC_COPY( col, color_map->diffuse[index] );
    }
    else if ( val < threshold && val > -threshold )
    {
        VEC_COPY( col, color_map->diffuse[index] );
    }
    else if ( val > result_max )
    {
        VEC_COPY( col, v_win->colormap[255] );
    }
    else if ( val < result_min )
    {
        VEC_COPY( col, v_win->colormap[0] );
    }
    else if ( result_min == result_max )
    {
        VEC_COPY( col, v_win->colormap[1] );
        if ( p_analy->damage_result )
        {
            if ( val>0 )
                VEC_COPY( col, v_win->colormap[255]); /* Show damage as RED */
        }
    }
    else
    {
        scl_max = (float) CMAP_SIZE - 2.01;

        if (logscale)
        {
            idx = (int)( scl_max * ((val-result_min)/(result_max-result_min))) + 1;

            /* Shift data to positive range */
            log_scale_data_shift( val,      result_min, result_max,
                                  &new_val, &new_result_min,
                                  &new_result_max, &new_result_shift,
                                  &new_result_mult);

            new_val_log        = log((double)new_val);
            new_result_min_log = log((double)new_result_min);
            new_result_max_log = log((double)new_result_max);

            idx_new = (int)( scl_max * ((new_val_log-new_result_min_log)/
                                        (new_result_max_log-new_result_min_log)) ) + 1;
            idx = idx_new;
        }
        else
        {
            idx = (int)( scl_max * (val-result_min)/(result_max-result_min) ) + 1;
        }

        if ( idx>=0 )
            VEC_COPY( col, v_win->colormap[idx] );
    }
}


/************************************************************
 * TAG( surf_color_lookup )
 *
 * Lookup a color in the colormap, and store it in the first
 * argument.  If file scope variable colorflag is FALSE, the
 * routine just returns the color property color.
 */
static void
surf_color_lookup( float col[4], float val, float result_min, float result_max,
                   float threshold, int index, Bool_type logscale)
{
    int idx, idx_new;
    float scl_max;


    float new_val,     new_result_min,     new_result_max,     new_result_shift;
    float new_result_mult;
    float new_val_log, new_result_min_log, new_result_max_log;

    /* Get alpha. */

    col[3] = v_win->surfaces.diffuse[index][3];

    /* If not doing colormapping or if result is near zero, use the
     * default material color.  Otherwise, do the color table lookup.
     */
    if ( !colorflag )
    {
        VEC_COPY( col, v_win->surfaces.diffuse[index] );
    }
    else if ( val < threshold && val > -threshold )
    {
        VEC_COPY( col, v_win->surfaces.diffuse[index] );
    }
    else if ( val > result_max )
    {
        VEC_COPY( col, v_win->colormap[255] );
    }
    else if ( val < result_min )
    {
        VEC_COPY( col, v_win->colormap[0] );
    }
    else if ( result_min == result_max )
    {
        VEC_COPY( col, v_win->colormap[1] );
    }
    else
    {
        scl_max = (float) CMAP_SIZE - 2.01;

        if (logscale)
        {
            idx = (int)( scl_max * ((val-result_min)/(result_max-result_min))) + 1;

            /* Shift data to positive range */
            log_scale_data_shift( val,      result_min, result_max,
                                  &new_val, &new_result_min,
                                  &new_result_max, &new_result_shift,
                                  &new_result_mult);

            new_val_log        = log((double)new_val);
            new_result_min_log = log((double)new_result_min);
            new_result_max_log = log((double)new_result_max);

            idx_new = (int)( scl_max * ((new_val_log-new_result_min_log)/
                                        (new_result_max_log-new_result_min_log)) ) + 1;
            idx = idx_new;
        }
        else
        {
            idx = (int)( scl_max * (val-result_min)/(result_max-result_min) ) + 1;
        }

        VEC_COPY( col, v_win->colormap[idx] );
    }
}


/*
 * SECTION_TAG( Draw grid )
 */


/************************************************************
 * TAG( draw_mesh_view )
 *
 * Redraw the grid window, using the current viewing
 * transformation to draw the grid.  This routine assumes
 * that the calling routine has set the current GL window
 * to be the grid view window.
 */
void
draw_mesh_view( Analysis *analy )
{
    Transf_mat look_rot;
    float arr[16], scal;
    RGB_raster_obj *p_rro;

    if ( env.timing )
    {
        wrt_text( "Timing for rendering...\n" );
        check_timing( 0 );
    }

    /* Set cursor to let user know something is happening. */
#ifdef SERIAL_BATCH
#else
    set_alt_cursor( CURSOR_WATCH );
#endif

    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT |
             GL_STENCIL_BUFFER_BIT );

    if ( analy->show_background_image && analy->background_image != NULL )
    {
        p_rro = analy->background_image;

        glDisable( GL_DEPTH_TEST );
        memory_to_screen( FALSE, p_rro->img_width, p_rro->img_height,
                          p_rro->alpha, p_rro->raster );
        glEnable( GL_DEPTH_TEST );
    }

    glPushMatrix();

    /*
     * Set up all the viewing transforms.  Transformations are
     * specified in opposite order to the order in which they
     * are applied.
     */
    look_rot_mat( v_win->look_from, v_win->look_at,
                  v_win->look_up, &look_rot );
    mat_to_array( &look_rot, arr );
    glMultMatrixf( arr );
    glTranslatef( -v_win->look_from[0], -v_win->look_from[1],
                  -v_win->look_from[2] );
    glTranslatef( v_win->trans[0], v_win->trans[1], v_win->trans[2] );
    mat_to_array( &v_win->rot_mat, arr );
    glMultMatrixf( arr );
    scal = v_win->bbox_scale;
    glScalef( scal*v_win->scale[0], scal*v_win->scale[1],
              scal*v_win->scale[2] );
    glTranslatef( v_win->bbox_trans[0], v_win->bbox_trans[1],
                  v_win->bbox_trans[2] );

    /* Draw the grid. */
    if ( analy->dimension == 3 )
        draw_grid( analy );
    else
        draw_grid_2d( analy );

    glPopMatrix();

    /* Draw all the foreground stuff. */

#ifdef SERIAL_BATCH
    draw_foreground( analy );
#else
    draw_foreground( analy );

    gui_swap_buffers();

    unset_alt_cursor();
#endif

    if ( env.timing )
        check_timing( 1 );
}


/************************************************************
 * TAG( draw_grid )
 *
 * Draw the element grid.  This routine handles all procedures
 * associated with drawing the grid, except for the viewing
 * transformation.
 */
extern void
draw_grid( Analysis *analy )
{
    Triangle_poly *tri;
    Surface_data *p_sd;
    Bool_type show_node_result, show_mesh_result, showgs=FALSE;
    Bool_type show_mat_result, show_surf_result;
    Bool_type composite_show;
    float verts[4][3];
    float norms[4][3];
    float cols[4][4];
    float vals[4];
    float rdiff;
    float scl_max;
    float rmin, rmax;
    int i, j, k;
    Bool_type show_result;
    Mesh_data *p_mesh;
    unsigned char *disable_mtl;
    Specified_obj *p_so;
    Classed_list *p_cl;
    int qty_classes;
    MO_class_data **mo_classes;
    MO_class_data *p_mo_class;
    Htable_entry *p_hte;
    int rval;

    /* Cap for colorscale interpolation. */
    scl_max = SCL_MAX;

    /* Enable color to change AMBIENT & DIFFUSE property of the material. */
    glEnable( GL_COLOR_MATERIAL );

    /* Turn lighting on. */
    if ( v_win->lighting )
        glEnable( GL_LIGHTING );

    /* Bias the depth-buffer if edges are on so they'll render in front. */
    if ( analy->show_edges )
        glDepthRange( analy->edge_zbias, 1 );

    p_mesh = MESH_P( analy );

    show_node_result = result_has_superclass( analy->cur_result, G_NODE,
                       analy );
    show_mat_result = result_has_superclass( analy->cur_result, G_MAT, analy );
    show_mesh_result = result_has_superclass( analy->cur_result, G_MESH,
                       analy );
    show_surf_result = result_has_superclass( analy->cur_result, G_SURFACE, analy );

    /*
     * Draw iso-surfaces.
     */

    if ( analy->show_isosurfs /* && !analy->show_carpet */ )
    {
        /* Difference between result min and max. */
        rdiff = analy->result_mm[1] - analy->result_mm[0];
        if ( rdiff == 0.0 )
            rdiff = 1.0;

        /* Use two-sided lighting. */
        glLightModeli( GL_LIGHT_MODEL_TWO_SIDE, 1 );

        /* Set up for polygon drawing. */
        begin_draw_poly( analy );

        for ( tri = analy->isosurf_poly; tri != NULL; tri = tri->next )
        {
            /*
             * Only compute once per polygon, since polygons are constant
             * color.
             */
            k = (int)( (tri->result[0] - analy->result_mm[0]) *
                       scl_max / rdiff ) + 1;

            for ( i = 0; i < 3; i++ )
            {
                VEC_COPY( cols[i], v_win->colormap[k] );
                VEC_COPY( norms[i], tri->norm );
                VEC_COPY( verts[i], tri->vtx[i] );
                vals[i] = tri->result[0];
            }
            draw_poly( 3, verts, norms, cols, vals, -1, p_mesh, analy, FALSE );
        }

        /* Back to the defaults. */
        end_draw_poly( analy );
        glLightModeli( GL_LIGHT_MODEL_TWO_SIDE, 0 );
    }

    /*
     * Draw cut planes.
     */
    /**/

    /*  Set glMaterial to draw from the correct color property data base */
    change_current_color_property( &v_win->mesh_materials, v_win->mesh_materials.current_index );


    if ( analy->show_cut /* && !analy->show_carpet */ )
    {
        composite_show = show_node_result || show_mat_result
                         || show_mesh_result;

        get_min_max( analy, ( show_mat_result || show_mesh_result ), &rmin,
                     &rmax );

        /* Use two-sided lighting. */
        glLightModeli( GL_LIGHT_MODEL_TWO_SIDE, 1 );

        /* Set up for polygon drawing. */
        begin_draw_poly( analy );

        for ( p_cl = analy->cut_poly_lists; p_cl != NULL; NEXT( p_cl ) )
        {
            p_mesh = analy->mesh_table + p_cl->mo_class->mesh_id;
            disable_mtl = p_mesh->disable_material;

            show_result = composite_show
                          || result_has_class( analy->cur_result,
                                               p_cl->mo_class, analy );

            for ( tri = (Triangle_poly *) p_cl->list; tri != NULL;
                    NEXT( tri ) )
            {
                colorflag = show_result && !disable_mtl[tri->mat];
                if ( analy->material_greyscale && disable_mtl[tri->mat] )
                    showgs = TRUE;
                else
                    showgs = FALSE;

                for ( i = 0; i < 3; i++ )
                {
                    color_lookup( cols[i], tri->result[i], rmin, rmax,
                                  analy->zero_result, tri->mat,  analy->logscale,
                                  showgs );
                    VEC_COPY( norms[i], tri->norm );
                    VEC_COPY( verts[i], tri->vtx[i] );
                    vals[i] = tri->result[i];
                }
                draw_poly( 3, verts, norms, cols, vals, -1, p_mesh, analy, FALSE );
            }
        }

        /* Back to the defaults. */
        end_draw_poly( analy );
        glLightModeli( GL_LIGHT_MODEL_TWO_SIDE, 0 );
    }

    /* Back to the default. */
    glDisable( GL_COLOR_MATERIAL );


    /*
     * Draw external surface polygons.
     */


    /********************************************************************
     *IRC: Added Oct 17, 2006. Offset Polygon lines so that lines
     *      are rendered correctly on XWin-32 systems.
     */

    if ( analy->edge_zbias == (float) DFLT_ZBIAS && env.win32 )
    {
        glEnable(GL_POLYGON_OFFSET_FILL);
        glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
        glPolygonOffset( 1.0, 1 );
        glDepthRange( 0, 2 );

        /* analy->edge_zbias = 20;
           glDepthRange( 0, 1 - analy->edge_zbias ); */
    }

    /* IRC: Added Oct 17, 2006. Offset Polygon lines so that lines do not
     *      show through the faces.
     *********************************************************************
     */


    if ( analy->show_extern_polys )
        draw_extern_polys( analy );

    /*
     * Draw reference surface polygons.
     */
    if ( analy->show_ref_polys )
        draw_ref_polys( analy );

    /*
     * Draw free nodes.
     */
    if ( analy->free_nodes_enabled || analy->particle_nodes_enabled )
        draw_free_nodes( analy ); 
      /*if(analy->free_nodes_enabled)
      {
          draw_free_nodes(analy);
      } */

    /*
     * Draw the elements.
     */

    p_mesh = MESH_P( analy );

    /* Hex element classes. */
    qty_classes = p_mesh->classes_by_sclass[G_HEX].qty;
    mo_classes = (MO_class_data **) p_mesh->classes_by_sclass[G_HEX].list;
    for ( i = 0; i < qty_classes; i++ )
        draw_hexs( show_node_result, show_mat_result, show_mesh_result,
                   mo_classes[i], analy );

    if ( analy->ei_result && analy->result_active )
        show_mat_result = TRUE;

    /* Tet element classes. */
    qty_classes = p_mesh->classes_by_sclass[G_TET].qty;
    mo_classes = (MO_class_data **) p_mesh->classes_by_sclass[G_TET].list;
    for ( i = 0; i < qty_classes; i++ )
        draw_tets( show_node_result, show_mat_result, show_mesh_result,
                   mo_classes[i], analy );

    /* Quad element classes. */
    qty_classes = p_mesh->classes_by_sclass[G_QUAD].qty;
    mo_classes = (MO_class_data **) p_mesh->classes_by_sclass[G_QUAD].list;
    for ( i = 0; i < qty_classes; i++ )
        draw_quads_3d( show_node_result, show_mat_result, show_mesh_result,
                       mo_classes[i], analy );

    /* Tri element classes. */
    qty_classes = p_mesh->classes_by_sclass[G_TRI].qty;
    mo_classes = (MO_class_data **) p_mesh->classes_by_sclass[G_TRI].list;
    for ( i = 0; i < qty_classes; i++ )
        draw_tris_3d( show_node_result, show_mat_result, show_mesh_result,
                      mo_classes[i], analy );

    /* Pyramid element classes. */
    qty_classes = p_mesh->classes_by_sclass[G_PYRAMID].qty;
    mo_classes = (MO_class_data **) p_mesh->classes_by_sclass[G_PYRAMID].list;
    for ( i = 0; i < qty_classes; i++ )
        draw_pyramids( show_node_result, show_mat_result, show_mesh_result,
                       mo_classes[i], analy );

    /* Turn lighting off (back to the default). */
    if ( v_win->lighting )
        glDisable( GL_LIGHTING );

    /* Remove depth-buffer bias applied to polygons. */
    if ( analy->show_edges && !env.win32 )
        glDepthRange( 0, 1 );

    /* Beam element classes. */
    qty_classes = p_mesh->classes_by_sclass[G_BEAM].qty;
    mo_classes = (MO_class_data **) p_mesh->classes_by_sclass[G_BEAM].list;
    for ( i = 0; i < qty_classes; i++ )
        draw_beams_3d( show_node_result, show_mat_result, show_mesh_result, mo_classes[i],
                       analy );

    /* Discrete element classes. */
    qty_classes = p_mesh->classes_by_sclass[G_TRUSS].qty;
    mo_classes = (MO_class_data **) p_mesh->classes_by_sclass[G_TRUSS].list;
    for ( i = 0; i < qty_classes; i++ )
        draw_truss_3d( show_mat_result, show_mesh_result, mo_classes[i],
                       analy );


    /*  Set glMaterial to draw from the correct color property data base */
    if ( MESH_P( analy )->classes_by_sclass[G_SURFACE].qty > 0 )
    {
        change_current_color_property( &v_win->surfaces,
                                       v_win->surfaces.current_index );
    }


    /* Surface element classes. */
    qty_classes = p_mesh->classes_by_sclass[G_SURFACE].qty;
    mo_classes = (MO_class_data **) p_mesh->classes_by_sclass[G_SURFACE].list;
    for ( i = 0; i < qty_classes; i++ )
        draw_surfaces_3d( v_win->current_color_property,
                          show_node_result, show_surf_result, show_mesh_result,
                          mo_classes[ i ], analy );

    /*
     * Draw point cloud.  Colors are assigned from mesh materials.
     */

    /* Node element classes. */
    if ( analy->mesh_view_mode == RENDER_POINT_CLOUD )
    {
        qty_classes = p_mesh->classes_by_sclass[G_NODE].qty;
        mo_classes = (MO_class_data **) p_mesh->classes_by_sclass[G_NODE].list;
        for ( i = 0; i < qty_classes; i++ )
            draw_nodes_2d_3d( mo_classes[i], analy );
    }

    /*
     * Draw particle data.
     */

    qty_classes = p_mesh->classes_by_sclass[G_PARTICLE].qty;
    mo_classes = (MO_class_data **) p_mesh->classes_by_sclass[G_PARTICLE].list;
    if(qty_classes > 0)
    {
        /*analy->show_particle_class = FALSE; */
        analy->show_particle_class = TRUE;
        
    } 

    /*  Set glMaterial to draw from the correct color property data base */
    if ( analy->show_particle_class )
    {
        strcpy(particle_name, "DBC1");
        for(i=0; i < qty_classes; i++)
        {
            /*strcpy(particle_cname, mo_classes[i]->short_name);
         
            v_win->particles = v_win->current_color_property; 
            /*change_current_color_property( &v_win->particles,
                                           v_win->particles.current_index );

            rval = htable_search( p_mesh->class_table, particle_name, FIND_ENTRY,
                                  &p_hte );
            if ( rval == OK )
            {
                p_mo_class = (MO_class_data *) p_hte->data;

                rval = htable_search( analy->primal_results, "partpos", FIND_ENTRY,
                                      &p_hte );
                if ( rval == OK )
                {
                    draw_particles_3d( p_mo_class, analy );
                }
            } */
            draw_particles_3d( mo_classes[i], analy);
        } 
    }

    /*
     * Draw mesh edges.
     */
    if ( analy->show_edges || analy->show_edges_vec )
        draw_edges_3d( analy );

    /*
     * Draw contours.
     */
    if ( analy->show_contours )
    {
        Contour_obj *cont;

        /*        lsetdepth( v_win->z_front_near, v_win->z_front_far ); */

        antialias_lines( TRUE, analy->z_buffer_lines );
        glLineWidth( (GLfloat) analy->contour_width );

        glColor3fv( v_win->contour_color );
        for ( cont = analy->contours; cont != NULL; cont = cont->next )
            draw_line( cont->cnt, cont->pts, -1, p_mesh, analy, FALSE, NULL );

        glLineWidth( (GLfloat) 1.0 );
        antialias_lines( FALSE, 0 );

        /* Reset near/far planes to the defaults. */
        /*        lsetdepth( v_win->z_default_near, v_win->z_default_far ); */
    }

    /*
     * Draw vector grid.
     */
    if ( analy->show_vectors )
    {
        if ( analy->vectors_at_nodes == TRUE )
            draw_node_vec_2d_3d( analy );
        else if ( analy->have_grid_points )
            draw_vec_result_3d( analy );
        else
            popup_dialog( INFO_POPUP, "No grid defined for vectors." );
    }

    /*
     * Draw regular vector carpet.
     */
    /*
        if ( analy->show_carpet )
            draw_reg_carpet( analy );
    */

    /*
     * Draw grid vector carpet or iso or cut surface carpet.
     */
    /*
        if ( analy->show_carpet )
        {
            draw_vol_carpet( analy );
            draw_shell_carpet( analy );
        }
    */

    /*
     * Draw reference surface vector carpet.
     */
    /*
        if ( analy->show_carpet )
            draw_ref_carpet( analy );
    */

    /*
     * Draw particle traces.
     */
    if ( analy->show_traces )
        draw_traces( analy );

    /*
     * Highlight a node or element.
     */
    if ( analy->hilite_class != NULL )
        draw_hilite( TRUE, analy->hilite_class, analy->hilite_num, analy );

    /*
     * Highlight selected nodes and elements.
     */
    for ( p_so = analy->selected_objects; p_so != NULL; NEXT( p_so ) )
        draw_hilite( FALSE, p_so->mo_class, p_so->ident, analy );

    /*
     * Draw local coordinate frame for selected shell or hex elements.
     */
    if ( analy->loc_ref )
    {
        draw_locref( analy );
        draw_locref_hex( analy );
    }

    /*
     * Draw selected superclasses numbers.
     */
    if ( analy->show_num )
        draw_class_numbers( analy );

    /*
     * Draw wireframe bounding box.
     */
    if ( analy->show_bbox )
        draw_bbox( analy->bbox );

    glDisable(GL_POLYGON_OFFSET_FILL);

}


/************************************************************
 * TAG( draw_grid_2d )
 *
 * Draw the element grid.  This routine handles all procedures
 * associated with drawing the grid, except for the viewing
 * transformation.
 */
extern void
draw_grid_2d( Analysis *analy )
{
    int i;
    Specified_obj *p_so;
    MO_class_data **mo_classes;
    MO_class_data *p_mo_class;
    int qty_classes;
    Mesh_data *p_mesh;
    Bool_type show_node_result, show_mat_result, show_mesh_result;
    Bool_type show_surf_result;
    Htable_entry *p_hte;
    int rval;
    /*
     * Draw external surface polygons.
     */
    /*
        if ( analy->show_extern_polys )
            draw_extern_polys( analy );
    */

    /*
     * Draw the elements.
     */

    p_mesh = MESH_P( analy );
    show_node_result = result_has_superclass( analy->cur_result, G_NODE,
                       analy );
    show_mat_result = result_has_superclass( analy->cur_result, G_MAT, analy );
    show_mesh_result = result_has_superclass( analy->cur_result, G_MESH,
                       analy );

    /* Quad element classes. */
    qty_classes = p_mesh->classes_by_sclass[G_QUAD].qty;
    mo_classes = (MO_class_data **) p_mesh->classes_by_sclass[G_QUAD].list;
    for ( i = 0; i < qty_classes; i++ )
        draw_quads_2d( show_node_result, show_mat_result, show_mesh_result,
                       mo_classes[i], analy );

    /* Triangle element classes. */
    qty_classes = p_mesh->classes_by_sclass[G_TRI].qty;
    mo_classes = (MO_class_data **) p_mesh->classes_by_sclass[G_TRI].list;
    for ( i = 0; i < qty_classes; i++ )
        draw_tris_2d( show_node_result, show_mat_result, show_mesh_result,
                      mo_classes[i], analy );

    /* Beam element classes. */
    qty_classes = p_mesh->classes_by_sclass[G_BEAM].qty;
    mo_classes = (MO_class_data **) p_mesh->classes_by_sclass[G_BEAM].list;
    for ( i = 0; i < qty_classes; i++ )
        draw_beams_2d( show_node_result, show_mat_result, show_mesh_result, mo_classes[i],
                       analy );

    /* Discrete Element classes. */
    qty_classes = p_mesh->classes_by_sclass[G_BEAM].qty;
    mo_classes = (MO_class_data **) p_mesh->classes_by_sclass[G_BEAM].list;
    for ( i = 0; i < qty_classes; i++ )
        draw_truss_2d( show_mat_result, show_mesh_result, mo_classes[i],
                       analy );

    /* Surface element classes. */
    qty_classes = p_mesh->classes_by_sclass[G_SURFACE].qty;
    mo_classes = (MO_class_data **) p_mesh->classes_by_sclass[G_SURFACE].list;
    for ( i = 0; i < qty_classes; i++ )
        draw_surfaces_2d( v_win->current_color_property,
                          show_node_result, show_surf_result, show_mesh_result,
                          mo_classes[ i ], analy );

    /*
     * Draw point cloud.
     */
    if ( analy->mesh_view_mode == RENDER_POINT_CLOUD )
    {
        qty_classes = p_mesh->classes_by_sclass[G_NODE].qty;
        mo_classes = (MO_class_data **) p_mesh->classes_by_sclass[G_NODE].list;
        for ( i = 0; i < qty_classes; i++ )
            draw_nodes_2d_3d( mo_classes[i], analy );
    }

    /*
     * Draw particle data.
     */
    if ( analy->show_particle_class )
    {
        rval = htable_search( p_mesh->class_table, particle_cname, FIND_ENTRY,
                              &p_hte );
        if ( rval == OK )
        {
            p_mo_class = (MO_class_data *) p_hte->data;

            rval = htable_search( analy->primal_results, "partpos", FIND_ENTRY,
                                  &p_hte );
            if ( rval == OK )
                draw_particles_2d( p_mo_class, analy );
        }
    }

    /*
     * Draw mesh edges.
     */
    if ( analy->show_edges || analy->show_edges_vec )
        draw_edges_2d( analy );

    /*
     * Draw contours.
     */
    if ( analy->show_contours )
    {
        Contour_obj *cont;

        antialias_lines( TRUE, analy->z_buffer_lines );
        glLineWidth( (GLfloat) analy->contour_width );

        glColor3fv( v_win->contour_color );
        for ( cont = analy->contours; cont != NULL; cont = cont->next )
            draw_line( cont->cnt, cont->pts, -1, p_mesh, analy, FALSE, NULL );

        glLineWidth( (GLfloat) 1.0 );
        antialias_lines( FALSE, 0 );
    }

    /*
     * Draw vector result.
     */
    if ( analy->show_vectors )
    {
        if ( analy->vectors_at_nodes == TRUE )
            draw_node_vec_2d_3d( analy );
        else if ( analy->have_grid_points )
            draw_vec_result_2d( analy );
        else
            popup_dialog( INFO_POPUP, "No grid defined for vectors." );
    }

    /*
     * Draw grid vector carpet.
     */
    /*
        if ( analy->show_carpet )
            draw_shell_carpet( analy );
    */

    /*
     * Draw particle traces.
     */

    if ( analy->show_traces )
        draw_traces( analy );


    /*
     * Draw selected class numbers.
     */

    if ( analy->show_num )
        draw_class_numbers( analy );

    /*
     * Highlight a node or element.
     */
    if ( analy->hilite_class != NULL )
        draw_hilite( TRUE, analy->hilite_class, analy->hilite_num, analy );

    /*
     * Highlight selected nodes and elements.
     */
    for ( p_so = analy->selected_objects; p_so != NULL; NEXT( p_so ) )
        draw_hilite( FALSE, p_so->mo_class, p_so->ident, analy );

    /*
     * Draw wireframe bounding box.
     */
    if ( analy->show_bbox )
        draw_bbox( analy->bbox );
}


/************************************************************
 * TAG( check_for_tri_face )
 *
 * Check for a triangular brick face (i.e., one of the four
 * face nodes is repeated.)  if the face is triangular, reorder
 * the nodes so the repeated node is last.  The new ordering
 * is returned in "ord".  The routine returns the number of
 * nodes for the face (3 or 4).
 */
static int
check_for_tri_face( int fc, int el, MO_class_data *p_hex_class,
                    float verts[4][3], int ord[4] )
{
    float tverts[4][3];
    int shift, nodes[4];
    int i, j;
    int (*connects)[8];

    if ( p_hex_class->objects.elems->has_degen )
    {
        connects = (int (*)[8]) p_hex_class->objects.elems->nodes;

        for ( i = 0; i < 4; i++ )
            nodes[i] = connects[el][ fc_nd_nums[fc][i] ];

        for ( i = 0; i < 4; i++ )
            if ( nodes[i] == nodes[(i+1)%4] )
            {
                /* Triangle. */
                if ( i == 3 )
                    shift = 3;
                else
                    shift = 2 - i;

                for ( j = 0; j < 4; j++ )
                    ord[j] = ( j + 4 - shift ) % 4;

                /* Reorder the vertices. */
                for ( j = 0; j < 4; j++ )
                {
                    VEC_COPY( tverts[j], verts[ ord[j] ] );
                }
                for ( j = 0; j < 4; j++ )
                {
                    VEC_COPY( verts[j], tverts[j] );
                }

                return( 3 );
            }
    }

    /* Quad. */
    for ( j = 0; j < 4; j++ )
        ord[j] = j;

    return( 4 );
}


/************************************************************
 * TAG( get_min_max )
 *
 * Get the appropriate result min and max.
 *
 */
static void
get_min_max( Analysis *analy, Bool_type no_interp, float *p_min, float *p_max )
{
    /* Nodal results never get an element min/max. */
    if ( result_has_superclass( analy->cur_result, G_NODE, analy ) )
    {
        *p_min = analy->result_mm[0];
        *p_max = analy->result_mm[1];
    }
    else if ( analy->interp_mode == NO_INTERP || no_interp )
    {
        if ( analy->mm_result_set[0] )
            *p_min = analy->result_mm[0];
        else if ( analy->use_global_mm )
            *p_min = analy->elem_global_mm.object_minmax[0];
        else
            *p_min = analy->elem_state_mm.object_minmax[0];

        if ( analy->mm_result_set[1] )
            *p_max = analy->result_mm[1];
        else if ( analy->use_global_mm )
            *p_max = analy->elem_global_mm.object_minmax[1];
        else
            *p_max = analy->elem_state_mm.object_minmax[1];
    }
    else
    {
        *p_min = analy->result_mm[0];
        *p_max = analy->result_mm[1];
    }
}


/************************************************************
 * TAG( draw_hexs )
 *
 * Draw the external faces of hex volume elements in the model.
 */
static void
draw_hexs( Bool_type show_node_result, Bool_type show_mat_result,
           Bool_type show_mesh_result, MO_class_data *p_hex_class,
           Analysis *analy )
{
    Bool_type show_result;
    float verts[4][3];
    float norms[4][3];
    float cols[4][4];
    float res[4];  /* "val" was formerly declared but never referenced  LAS */
    float rmin, rmax;
    float tverts[4][3], v1[3], v2[3], nor[3], ray[3];
    float *activity=NULL;

    int matl, el, fc, nd, ord[4], cnt, face_qty, mesh_idx;
    int i, j, k;
    int (*connects)[8];
    int *face_el, *face_fc;
    int *mtl;
    int *p_index_source;
    float *data_array;
    unsigned char *disable_mtl;
    Mesh_data *p_mesh;
    Visibility_data *p_vd;
    Interp_mode_type save_interp_mode;
    unsigned char *hide_mtl;
    int *p_mats;
    Bool_type disable_gray = TRUE;
    Result *p_result;
    Subrec_obj *p_subrec;
    p_result = analy->cur_result;
/* Adding code to determine if a subrecord for these objects has the appropriate variable and if
 * not it sets that element to a gray scale prior to calling draw poly */
    int index, subrec, srec, qty_blocks, mesh_id;
    Bool_type same_as_last = FALSE;
    char **results_map;
    char * result_name = NULL;
    char * short_name = NULL;
    
    static int first = 0;

    mesh_id = p_hex_class->mesh_id;

    show_result = show_node_result
                  || result_has_class( analy->cur_result, p_hex_class, analy )
                  || show_mat_result
                  || show_mesh_result;

    if(result_has_class(analy->cur_result, p_hex_class, analy))
    {
        disable_gray = FALSE;
    }

    /*if(first == 0)
    {
        results_map = (char **)malloc(analy->mesh_qty*sizeof(char *));
        if(results_map == NULL)
        {
            popup_dialog(WARNING_POPUP, "Out of memory in function draw_hexs. Terminating!!-n");
            parse_command("quit", analy);
        }
        for(i = 0; i < analy->mesh_qty; i++)
        {
            results_map[i] = NULL;
        } 

        first = 1;
    } */
        results_map = (char **)malloc(analy->mesh_qty*sizeof(char *));
        if(results_map == NULL)
        {
            popup_dialog(WARNING_POPUP, "Out of memory in function draw_hexs. Terminating!!-n");
            parse_command("quit", analy);
        }
        for(i = 0; i < analy->mesh_qty; i++)
        {
            results_map[i] = NULL;
        } 
   
    if(show_result)
    { 
        if(result_name != NULL)
        {
            if(strcmp(result_name, p_result->name))
            {
                /* we are showing a different result than the last shown result */
                free(result_name);
                result_name = NULL;
                if(results_map[mesh_id] != NULL)
                {
                    free(results_map[mesh_id]);
                    results_map[mesh_id] = NULL;
                }
                same_as_last = FALSE;
            } else
            {
                same_as_last = TRUE;
            }
        }
        if(short_name != NULL )
        {
            if(!strcmp(short_name, p_hex_class->short_name))
            {
                free(short_name);
                short_name = NULL;
                if(results_map[mesh_id] != NULL)
                {
                    free(results_map[mesh_id]);
                    results_map[mesh_id] = NULL;
                }
            }
        }  
       
        if(results_map[mesh_id] == NULL)
        {
            results_map[mesh_id] = (char *) malloc(p_hex_class->qty+1);
            if( results_map[mesh_id] == NULL)
            {
                popup_dialog(WARNING_POPUP, "Out of memory in function draw_hexs. Terminating!!\n");
                parse_command("quit", analy);
            }
        }
        if(result_name == NULL)
        {
            result_name = (char *) malloc(strlen(p_result->name) + 1);
            if(result_name == NULL)
            {
                popup_dialog(WARNING_POPUP, "Out of memory in function draw_hexs. Terminating!!\n");
                parse_command("quit", analy);
            }
            strncpy(result_name, p_result->name, strlen(p_result->name) + 1);
        }
        if(short_name == NULL)
        {
            short_name = (char *) malloc(strlen(p_hex_class->short_name) + 1);
            if(short_name == NULL)
            {
                popup_dialog(WARNING_POPUP, "Out of memory in function draw_hexs. Terminating!!\n");
                parse_command("quit", analy);
            }
            strncpy(short_name, p_hex_class->short_name, strlen(p_hex_class->short_name) + 1);
        }
    }

    if(show_result && analy->cur_result != NULL) 
    { 
        populate_result(G_HEX, results_map[mesh_id], p_hex_class->qty + 1, p_hex_class, analy);
    }


 
    Bool_type hidden_poly=FALSE,
              hidden_poly_elem=FALSE,
              hidden_poly_elem_wft=FALSE,
              hidden_poly_mat=FALSE,
              disable_flag=FALSE;
    Bool_type showgs=TRUE;

    if ( analy->mesh_view_mode != RENDER_FILLED          && analy->mesh_view_mode != RENDER_WIREFRAME &&
            analy->mesh_view_mode != RENDER_WIREFRAMETRANS  && analy->mesh_view_mode != RENDER_HIDDEN )
        return;

    if ( is_particle_class(analy, p_hex_class->superclass, p_hex_class->short_name) )
        return; 

    p_mesh = MESH_P( analy );
    connects = (int (*)[8]) p_hex_class->objects.elems->nodes;
    mtl = p_hex_class->objects.elems->mat;
    p_vd = p_hex_class->p_vis_data;
    face_qty = p_vd->face_cnt;
    face_el = p_vd->face_el;
    face_fc = p_vd->face_fc;
    disable_mtl = p_mesh->disable_material;
    hide_mtl    = p_mesh->hide_material;
    p_mats      = p_mesh->particle_mats;

    /* Override interpolation mode for material and mesh results. */
    save_interp_mode = analy->interp_mode;
    if ( show_mat_result || show_mesh_result )
        analy->interp_mode = NO_INTERP;

    activity = analy->state_p->sand_present
               ? analy->state_p->elem_class_sand[p_hex_class->elem_class_index]
               : NULL;

    /* Abstract the source array for data values and its index. */
    if ( analy->interp_mode == GOOD_INTERP )
    {
        data_array = NODAL_RESULT_BUFFER( analy );
        p_index_source = &nd;
    }
    else if ( show_mat_result )
    {
        data_array = ((MO_class_data **)
                      p_mesh->classes_by_sclass[G_MAT].list)[0]->data_buffer;
        p_index_source = &matl;
    }
    else if ( show_mesh_result )
    {
        mesh_idx = 0;
        data_array = ((MO_class_data **)
                      p_mesh->classes_by_sclass[G_MESH].list)[0]->data_buffer;
        p_index_source = &mesh_idx;
    }
    else if ( analy->interp_mode == NO_INTERP && !show_node_result )
    {
        data_array = p_hex_class->data_buffer;
        p_index_source = &el;
    }
    else
    {
        data_array = NODAL_RESULT_BUFFER( analy );
        p_index_source = &nd;
    }

    /* If we are drawing particles then look for a particle result and map onto hexes
      *
      */
    /* if ( show_result && analy->free_particles )
           data_array = analy->pn_buffer_ptr[0]; */

    get_min_max( analy, FALSE, &rmin, &rmax );

    /* Enable color to change AMBIENT & DIFFUSE property of the material. */
    glEnable( GL_COLOR_MATERIAL );

    if(analy->mesh_view_mode !=RENDER_WIREFRAMETRANS)
    {
        /* Throw away backward-facing faces. */
        glEnable( GL_CULL_FACE );
    }

    /* Set up for polygon drawing. */
    begin_draw_poly( analy );

    for ( i = 0; i < face_qty; i++ )
    {
        el = face_el[i];
        fc = face_fc[i];

        /*
         * Remove faces that are shared with quad elements, so
         * the polygons are not drawn twice.
         */
        if ( analy->shared_faces
                && p_mesh->classes_by_sclass[G_QUAD].qty > 0
                && face_matches_quad( el, fc, p_hex_class, p_mesh, analy ) )
            continue;

        /* Check for inactive elements. */
        if ( activity )
        {
            if ( activity[el] != 0.0 && analy->show_only_deleted_elements )
                continue;
            if ( activity[el] == 0.0 && !analy->show_deleted_elements && !analy->show_only_deleted_elements )
                continue;
        }

        get_hex_face_verts( el, fc, p_hex_class, analy, verts );

        /* Cull backfacing polygons by hand to speed things up. */
        if ( !analy->reflect && analy->manual_backface_cull )
        {
            for ( j = 0; j < 4; j++ )
                point_transform( tverts[j], verts[j], &cur_view_mat );
            VEC_SUB( v1, tverts[2], tverts[0] );
            VEC_SUB( v2, tverts[3], tverts[1] );
            VEC_CROSS( nor, v1, v2 );
            if ( v_win->orthographic )
            {
                VEC_SET( ray, 0.0, 0.0, -1.0 );
            }
            else
            {
                for ( k = 0; k < 3; k++ )
                    ray[k] = 0.25 * ( tverts[0][k] + tverts[1][k] +
                                      tverts[2][k] + tverts[3][k] );
            }

            if ( VEC_DOT( ray, nor ) >= 0.0 )
                continue;
        }

        /*
         * Check for triangular (degenerate) face and reorder nodes
         * if needed.
         */
        cnt = check_for_tri_face( fc, el, p_hex_class, verts, ord );

        matl = mtl[el];

        if ( v_win->mesh_materials.current_index != matl )
            change_current_color_property( &v_win->mesh_materials, matl );


        /* Colorflag is TRUE if displaying a result, otherwise
         * polygons are drawn in the material color.
         */
        disable_flag = disable_by_object_type( p_hex_class, matl, el, analy, data_array );

        colorflag = show_result && !disable_mtl[matl] && !disable_by_object_type( p_hex_class, matl, el, analy, data_array );

        colorflag = show_result && !disable_mtl[matl];
        colorflag = show_result && !disable_mtl[matl] && !disable_by_object_type( p_hex_class, matl, el, analy, data_array );

        if ( analy->material_greyscale && disable_mtl[matl])
            showgs = TRUE;
        else
            showgs = FALSE;

        /*
        * Array "ord[4]" represents a map from the order of vertex data
        * specified by the per-face node ordering to the order needed for
        * rendering.  For quad's, there is no difference (i.e., "ord"
        * contains [0, 1, 2, 3].  For tri's, which are actually degenerate
        * quads with one node repeated, we want the repeated node to be last
        * as we will only render the first three, so there is potentially
        * a re-ordering.  This requires all the associated data to be re-
        * ordered, so the data arrays passed to draw_poly() all have their
        * data stored taking the original data at position set by "ord[j]" and
        * storing it in position "j".
        */
        if ( analy->interp_mode == GOOD_INTERP )
        {
            for ( j = 0; j < 4; j++ )
            {
                nd = connects[el][ fc_nd_nums[fc][ord[j]] ];

                for ( k = 0; k < 3; k++ )
                    norms[j][k] = p_vd->face_norm[ ord[j] ][k][i];

                res[j] = data_array[nd];
            }
        }
        else
        {
            for ( j = 0; j < 4; j++ )
            {
                nd = connects[el][ fc_nd_nums[fc][ord[j]] ];

                for ( k = 0; k < 3; k++ )
                    norms[j][k] = p_vd->face_norm[ ord[j] ][k][i];

                color_lookup( cols[j], data_array[*p_index_source],
                              rmin, rmax, analy->zero_result, matl,
                              analy->logscale, showgs );
            }
        }

        hidden_poly_mat = hide_mtl[matl];

        hidden_poly_elem = hide_by_object_type( p_hex_class, matl, el, analy, data_array );
        if ( analy->mesh_view_mode == RENDER_WIREFRAMETRANS )
        {
            hidden_poly_elem_wft = TRUE;  
           /*hidden_poly_elem_wft = disable_by_object_type( p_hex_class, matl, el, analy, data_array ); */
        }

        if ( hidden_poly_mat || hidden_poly_elem || hidden_poly_elem_wft )
            hidden_poly = TRUE;
        else
            hidden_poly = FALSE;

        if ( analy->particle_nodes_hide_background && p_mats )
            if ( p_mats[matl] )
                hidden_poly = TRUE;

        /* check to see if this element has the variable to be shown */
        if(show_result && !analy->material_greyscale && analy->auto_gray && !disable_gray)
        {
            grayel = 0;
            if((results_map[mesh_id][el+1] == '0') && !disable_mtl[matl] && !disable_flag)
            {
                grayel = 1;         
                for(j = 0; j < 4; j++)
                {
                    for(k = 0; k < 4; k++)
                    {
                        cols[j][k] = grayval;
                    }
                }
               
            }
        } else if(!show_result && !analy->showmat &&!disable_mtl[matl] && !disable_flag && analy->auto_gray && !disable_gray)
        {
            grayel = 1;
            for(j = 0; j < 4; j++)
            {
                for(k = 0; k < 4; k++)
                {
                    cols[j][k] = grayval;
                }
            }
             
        }
  
        if(analy->mesh_view_mode == RENDER_WIREFRAMETRANS)
        {
            draw_edges_3d(analy);
        } else
        {
            draw_poly( cnt, verts, norms, cols, res, matl, p_mesh, analy, hidden_poly );
        }
    }

    grayel = 0;

    /* Back to the defaults. */
    end_draw_poly( analy );
    glDisable( GL_CULL_FACE );
    glDisable( GL_COLOR_MATERIAL );
    analy->interp_mode = save_interp_mode;
    if(results_map != NULL)
    {
        for(i = 0; i < analy->mesh_qty; i++)
        {
            if(results_map[i] != NULL)
            {
                free(results_map[i]);
                results_map[i] = NULL;
            }
        }
        free(results_map);
        results_map = NULL;
    }
}


/************************************************************
 * TAG( draw_tets )
 *
 * Draw the external faces of tet volume elements in the model.
 */
static void
draw_tets( Bool_type show_node_result, Bool_type show_mat_result,
           Bool_type show_mesh_result, MO_class_data *p_tet_class,
           Analysis *analy )
{
    Bool_type show_result, showgs=FALSE;
    float verts[3][3];
    float norms[3][3];
    float cols[3][4];
    float res[3];  /* "val" formerly declared but never referenced  LAS */
    float rmin, rmax;
    float tverts[3][3], v1[3], v2[3], nor[3], ray[3];
    int matl, el, fc, nd, cnt, face_qty, mesh_idx;
    int i, j, k;
    int (*connects)[4];
    int *face_el, *face_fc;
    int *mtl;
    unsigned char *disable_mtl;
    Mesh_data *p_mesh;
    Visibility_data *p_vd;
    float *data_array;
    int *p_index_source;
    Interp_mode_type save_interp_mode;
    Bool_type hidden_poly;
    Bool_type disable_flag=FALSE;

    int mesh_id;
    static int first = 0;
    static char ** results_map;
    static char * result_name = NULL;
    static char * short_name = NULL;
    Result * p_result;
    Bool_type same_as_last = FALSE;
    Bool_type disable_gray = TRUE;

    mesh_id = p_tet_class->mesh_id;
    p_result = analy->cur_result;


    show_result = show_node_result
                  || result_has_class( analy->cur_result, p_tet_class, analy )
                  || show_mat_result
                  || show_mesh_result; 

    if(result_has_class( analy->cur_result, p_tet_class, analy))
    {
        disable_gray = FALSE;
    }

    if(first == 0)
    {
        results_map = (char **)malloc(analy->mesh_qty*sizeof(char *));
        if(results_map == NULL)
        {
            popup_dialog(WARNING_POPUP, "Out of memory in function draw_hexs. Terminating!!-n");
            parse_command("quit", analy);
        }
        for(i = 0; i < analy->mesh_qty; i++)
        {
            results_map[i] = NULL;
        } 

        first = 1;
    }
   
    if(show_result)
    { 
        if(result_name != NULL)
        {
            if(strcmp(result_name, p_result->name))
            {
                /* we are showing a different result than the last shown result */
                free(result_name);
                result_name = NULL;
                if(results_map[mesh_id] != NULL)
                {
                    free(results_map[mesh_id]);
                    results_map[mesh_id] = NULL;
                }
                same_as_last = FALSE;
            } else
            {
                same_as_last = TRUE;
            }
        }
        if(short_name != NULL )
        {
            if(strcmp(short_name, p_tet_class->short_name))
            {
                free(short_name);
                short_name = NULL;
                if(results_map[mesh_id] != NULL)
                {
                    free(results_map[mesh_id]);
                    results_map[mesh_id] = NULL;
                }
            }
        }  
       
        if(results_map[mesh_id] == NULL)
        {
            results_map[mesh_id] = (char *) malloc(p_tet_class->qty+1);
            if( results_map[mesh_id] == NULL)
            {
                popup_dialog(WARNING_POPUP, "Out of memory in function draw_tets. Terminating!!\n");
                parse_command("quit", analy);
            }
        }
        if(result_name == NULL)
        {
            result_name = (char *) malloc(strlen(p_result->name) + 1);
            if(result_name == NULL)
            {
                popup_dialog(WARNING_POPUP, "Out of memory in function draw_tets. Terminating!!\n");
                parse_command("quit", analy);
            }
            strncpy(result_name, p_result->name, strlen(p_result->name) + 1);
        }
        if(short_name == NULL)
        {
            short_name = (char *) malloc(strlen(p_tet_class->short_name) + 1);
            if(short_name == NULL)
            {
                popup_dialog(WARNING_POPUP, "Out of memory in function draw_tets Terminating!!\n");
                parse_command("quit", analy);
            }
            strncpy(short_name, p_tet_class->short_name, strlen(p_tet_class->short_name) + 1);
        }
    }

    if(show_result && analy->cur_result != NULL && !same_as_last) 
    { 
        populate_result(G_TET, results_map[mesh_id], p_tet_class->qty + 1, p_tet_class, analy);
    }

    /*show_result = show_node_result
                  || result_has_class( analy->cur_result, p_tet_class, analy )
                  || show_mat_result
                  || show_mesh_result; */

    if ( analy->mesh_view_mode != RENDER_FILLED          && analy->mesh_view_mode != RENDER_WIREFRAME &&
            analy->mesh_view_mode != RENDER_WIREFRAMETRANS  && analy->mesh_view_mode != RENDER_HIDDEN )
        return;

    p_mesh = MESH_P( analy );
    connects = (int (*)[4]) p_tet_class->objects.elems->nodes;
    mtl = p_tet_class->objects.elems->mat;
    p_vd = p_tet_class->p_vis_data;
    face_qty = p_vd->face_cnt;
    face_el = p_vd->face_el;
    face_fc = p_vd->face_fc;
    disable_mtl = p_mesh->disable_material;

    /* Override interpolation mode for material and mesh results.  */
    save_interp_mode = analy->interp_mode;
    if ( show_mat_result || show_mesh_result )
        analy->interp_mode = NO_INTERP;

    /* Abstract the source array for data values and its index. */
    if ( analy->interp_mode == GOOD_INTERP )
    {
        data_array = NODAL_RESULT_BUFFER( analy );
        p_index_source = &nd;
    }
    else if ( show_mat_result )
    {
        data_array = ((MO_class_data **)
                      p_mesh->classes_by_sclass[G_MAT].list)[0]->data_buffer;
        p_index_source = &matl;
    }
    else if ( show_mesh_result )
    {
        mesh_idx = 0;
        data_array = ((MO_class_data **)
                      p_mesh->classes_by_sclass[G_MESH].list)[0]->data_buffer;
        p_index_source = &mesh_idx;
    }
    else if ( analy->interp_mode == NO_INTERP && !show_node_result )
    {
        data_array = p_tet_class->data_buffer;
        p_index_source = &el;
    }
    else
    {
        data_array = NODAL_RESULT_BUFFER( analy );
        p_index_source = &nd;
    }

    show_result = show_node_result
                  || result_has_class( analy->cur_result, p_tet_class, analy )
                  || show_mat_result
                  || show_mesh_result; 

    get_min_max( analy, FALSE, &rmin, &rmax );

    /* Enable color to change AMBIENT & DIFFUSE property of the material. */
    glEnable( GL_COLOR_MATERIAL );

    /* Throw away backward-facing faces. */
    glEnable( GL_CULL_FACE );

    /* Set up for polygon drawing. */
    begin_draw_poly( analy );

    for ( i = 0; i < face_qty; i++ )
    {
        el = face_el[i];
        fc = face_fc[i];

        /*
         * Remove faces that are shared with triangle elements, so
         * the polygons are not drawn twice.
         */
        if ( analy->shared_faces
                && p_mesh->classes_by_sclass[G_TRI].qty > 0
                && face_matches_tri( el, fc, p_tet_class, p_mesh, analy ) )
            continue;

        get_tet_face_verts( el, fc, p_tet_class, analy, verts );

        /* Cull backfacing polygons by hand to speed things up. */
        if ( !analy->reflect && analy->manual_backface_cull )
        {
            for ( j = 0; j < 3; j++ )
                point_transform( tverts[j], verts[j], &cur_view_mat );
            VEC_SUB( v1, tverts[1], tverts[0] );
            VEC_SUB( v2, tverts[2], tverts[1] );
            VEC_CROSS( nor, v1, v2 );
            if ( v_win->orthographic )
            {
                VEC_SET( ray, 0.0, 0.0, -1.0 );
            }
            else
            {
                for ( k = 0; k < 3; k++ )
                    ray[k] = 0.333 * (tverts[0][k] + tverts[1][k] +
                                      tverts[2][k]);
            }

            if ( VEC_DOT( ray, nor ) >= 0.0 )
                continue;
        }

        /* No logic to deal with degenerate tet's. */
        cnt = 3;

        matl = mtl[el];


        if ( v_win->mesh_materials.current_index != matl )
            change_current_color_property( &v_win->mesh_materials, matl );


        /* Colorflag is TRUE if displaying a result, otherwise
         * polygons are drawn in the material color.
         */
        colorflag = show_result && !disable_mtl[matl];
        if ( analy->material_greyscale && disable_mtl[matl] )
            showgs = TRUE;
        else
            showgs = FALSE;

        if ( analy->interp_mode == GOOD_INTERP )
        {
            for ( j = 0; j < 3; j++ )
            {
                nd = connects[el][ tet_fc_nd_nums[fc][j] ];

                for ( k = 0; k < 3; k++ )
                    norms[j][k] = p_vd->face_norm[j][k][i];

                res[j] = data_array[nd];
            }
        }
        else
        {
            for ( j = 0; j < 3; j++ )
            {
                nd = connects[el][ fc_nd_nums[fc][j] ];

                for ( k = 0; k < 3; k++ )
                    norms[j][k] = p_vd->face_norm[j][k][i];

                color_lookup( cols[j], data_array[*p_index_source],
                              rmin, rmax, analy->zero_result, matl,
                              analy->logscale, showgs );
            }
        }

        disable_flag = disable_by_object_type( p_tet_class, matl, el, analy, data_array );

        /* check to see if this element has the variable to be shown */
        if(show_result && !analy->material_greyscale && !disable_gray)
        {
            grayel = 0;
            if((results_map[mesh_id][i+1] == '0') && !disable_mtl[matl] && !disable_flag && analy->auto_gray)
            {
                grayel = 1;         
                for(j = 0; j < 3; j++)
                {
                    for(k = 0; k < 4; k++)
                    {
                        cols[j][k] = grayval;
                    }
                }
               
            }
        } else if(!show_result && !analy->showmat &&!disable_mtl[matl] && !disable_flag && analy->auto_gray && !disable_gray)
        {
            grayel = 1;
            for(j = 0; j < 3; j++)
            {
                for(k = 0; k < 4; k++)
                {
                    cols[j][k] = grayval;
                }
            }
             
        }

        hidden_poly = hide_by_object_type( p_tet_class, matl, el, analy, data_array );
        draw_poly( cnt, verts, norms, cols, res, matl, p_mesh, analy, hidden_poly );
    }

    /* Back to the defaults. */
    end_draw_poly( analy );
    glDisable( GL_CULL_FACE );
    glDisable( GL_COLOR_MATERIAL );
    analy->interp_mode = save_interp_mode;
}


/************************************************************
 * TAG( draw_quads_3d )
 *
 * Draw the quad elements in the model.
 */
static void
draw_quads_3d( Bool_type show_node_result, Bool_type show_mat_result,
               Bool_type show_mesh_result, MO_class_data *p_quad_class,
               Analysis *analy )
{
    Bool_type show_result, showgs=FALSE;
    Bool_type has_degen;
    float *activity;
    float verts[4][3];
    float norms[4][3];
    float cols[4][4];
    float res[4];
    float rmin, rmax;
    int matl, cnt, nd, i, j, k, mesh_idx;
    int (*connects)[4];
    int *mtl;
    unsigned char *disable_mtl, *hide_mtl;
    Mesh_data *p_mesh;
    Visibility_data *p_vd;
    int quad_qty;
    float *data_array;
    int *p_index_source;
    Interp_mode_type save_interp_mode;

    Bool_type hidden_poly;
    Bool_type inactive_element=FALSE;
    Bool_type disable_flag=FALSE;
    Result * p_result;

    int mesh_id;
    static int first = 0;
    static char ** results_map;
    static char * result_name = NULL;
    static char * short_name = NULL;
    Bool_type same_as_last = FALSE;
    Bool_type disable_gray = TRUE;

    p_result = analy->cur_result;
    mesh_id = p_quad_class->mesh_id;

    show_result = show_node_result
                  || result_has_class( analy->cur_result, p_quad_class, analy )
                  || show_mat_result
                  || show_mesh_result;

    if(result_has_class( analy->cur_result, p_quad_class, analy))
    {
        disable_gray = FALSE;
    }

    if(first == 0)
    {
        results_map = (char **)malloc(analy->mesh_qty*sizeof(char *));
        if(results_map == NULL)
        {
            popup_dialog(WARNING_POPUP, "Out of memory in function draw_hexs. Terminating!!-n");
            parse_command("quit", analy);
        }
        for(i = 0; i < analy->mesh_qty; i++)
        {
            results_map[i] = NULL;
        } 

        first = 1;
    }
   
    if(show_result)
    { 
        if(result_name != NULL)
        {
            if(strcmp(result_name, p_result->name))
            {
                /* we are showing a different result than the last shown result */
                free(result_name);
                result_name = NULL;
                if(results_map[mesh_id] != NULL)
                {
                    free(results_map[mesh_id]);
                    results_map[mesh_id] = NULL;
                }
                same_as_last = FALSE;
            } else
            {
                same_as_last = TRUE;
            }
        }
        if(short_name != NULL )
        {
            if(strcmp(short_name, p_quad_class->short_name))
            {
                free(short_name);
                short_name = NULL;
                if(results_map[mesh_id] != NULL)
                {
                    free(results_map[mesh_id]);
                    results_map[mesh_id] = NULL;
                }
            }
        }  
       
        if(results_map[mesh_id] == NULL)
        {
            results_map[mesh_id] = (char *) malloc(p_quad_class->qty+1);
            if( results_map[mesh_id] == NULL)
            {
                popup_dialog(WARNING_POPUP, "Out of memory in function draw_hexs. Terminating!!\n");
                parse_command("quit", analy);
            }
        }
        if(result_name == NULL)
        {
            result_name = (char *) malloc(strlen(p_result->name) + 1);
            if(result_name == NULL)
            {
                popup_dialog(WARNING_POPUP, "Out of memory in function draw_hexs. Terminating!!\n");
                parse_command("quit", analy);
            }
            strncpy(result_name, p_result->name, strlen(p_result->name) + 1);
        }
        if(short_name == NULL)
        {
            short_name = (char *) malloc(strlen(p_quad_class->short_name) + 1);
            if(short_name == NULL)
            {
                popup_dialog(WARNING_POPUP, "Out of memory in function draw_hexs. Terminating!!\n");
                parse_command("quit", analy);
            }
            strncpy(short_name, p_quad_class->short_name, strlen(p_quad_class->short_name) + 1);
        }
    }

    if(show_result && analy->cur_result != NULL && !same_as_last) 
    { 
        populate_result(G_QUAD, results_map[mesh_id], p_quad_class->qty + 1, p_quad_class, analy);
    }

    if ( analy->mesh_view_mode != RENDER_FILLED          && analy->mesh_view_mode != RENDER_WIREFRAME &&
            analy->mesh_view_mode != RENDER_WIREFRAMETRANS  && analy->mesh_view_mode != RENDER_HIDDEN )
        return;

    activity = analy->state_p->sand_present
               ? analy->state_p->elem_class_sand[p_quad_class->elem_class_index]
               : NULL;

    p_mesh = MESH_P( analy );
    connects = (int (*)[4]) p_quad_class->objects.elems->nodes;
    mtl = p_quad_class->objects.elems->mat;
    p_vd = p_quad_class->p_vis_data;
    quad_qty = p_quad_class->qty;
    disable_mtl = p_mesh->disable_material;
    hide_mtl = p_mesh->hide_material;
    has_degen = p_quad_class->objects.elems->has_degen;

    /* Override interpolation mode for material and mesh results. */
    save_interp_mode = analy->interp_mode;
    if ( show_mat_result || show_mesh_result )
        analy->interp_mode = NO_INTERP;

    /* Abstract the source array for data values and its index. */
    if ( analy->interp_mode == GOOD_INTERP )
    {
        data_array = NODAL_RESULT_BUFFER( analy );
        p_index_source = &nd;
    }
    else if ( show_mat_result )
    {
        data_array = ((MO_class_data **)
                      p_mesh->classes_by_sclass[G_MAT].list)[0]->data_buffer;
        p_index_source = &matl;
    }
    else if ( show_mesh_result )
    {
        mesh_idx = 0;
        data_array = ((MO_class_data **)
                      p_mesh->classes_by_sclass[G_MESH].list)[0]->data_buffer;
        p_index_source = &mesh_idx;
    }
    else if ( analy->interp_mode == NO_INTERP && !show_node_result )
    {
        data_array = p_quad_class->data_buffer;
        p_index_source = &i;
    }
    else
    {
        data_array = NODAL_RESULT_BUFFER( analy );
        p_index_source = &nd;
    }

   /* show_result = show_node_result
                  || result_has_class( analy->cur_result, p_quad_class, analy )
                  || show_mat_result
                  || show_mesh_result; */

    get_min_max( analy, FALSE, &rmin, &rmax );

    /* Enable color to change AMBIENT & DIFFUSE property of the material. */
    glEnable( GL_COLOR_MATERIAL );

    /* Use two-sided lighting for shells. */
    glLightModeli( GL_LIGHT_MODEL_TWO_SIDE, 1 );

    /* Set up for polygon drawing. */
    begin_draw_poly( analy );

    for ( i = 0; i < quad_qty; i++ )
    {
        /* Check for inactive elements. */
        if ( activity )
        {
            if ( activity[i] != 0.0 && analy->show_only_deleted_elements )
                continue;
            if ( activity[i] == 0.0 && !analy->show_deleted_elements && !analy->show_only_deleted_elements )
                continue;
        }

        /* Skip if this material is invisible. */
        matl = mtl[i];
        if ( hide_mtl[matl] || hide_by_object_type( p_quad_class, matl, i, analy, data_array ))
            continue;

        if ( v_win->mesh_materials.current_index != matl )
            change_current_color_property( &v_win->mesh_materials, matl );


        /* Colorflag is TRUE if displaying a result, otherwise
         * polygons are drawn in the material color.
         */
        colorflag = show_result && !disable_mtl[matl] &&
                    !disable_by_object_type( p_quad_class, matl, i, analy, NULL );
        if ( analy->material_greyscale && disable_mtl[matl] )
            showgs = TRUE;
        else
            showgs = FALSE;

        get_quad_verts_3d( i, p_quad_class->objects.elems->nodes,
                           p_quad_class->mesh_id, analy, verts );

        if ( analy->interp_mode == GOOD_INTERP )
        {
            for ( j = 0; j < 4; j++ )
            {
                nd = connects[i][j];

                for ( k = 0; k < 3; k++ )
                    norms[j][k] = p_vd->face_norm[j][k][i];

                res[j] = data_array[nd];
            }
        }
        else
        {
            for ( j = 0; j < 4; j++ )
            {
                nd = connects[i][j];

                for ( k = 0; k < 3; k++ )
                    norms[j][k] = p_vd->face_norm[j][k][i];

                color_lookup( cols[j], data_array[*p_index_source],
                              rmin, rmax, analy->zero_result, matl,
                              analy->logscale, showgs );
            }
        }

        /* Check for triangular (degenerate) quad element. */
        cnt = 4;
        if ( has_degen &&
                connects[i][2] == connects[i][3] )
            cnt = 3;
      
        disable_flag = disable_by_object_type(p_quad_class, matl, i, analy, NULL);

        /* check to see if this element has the variable to be shown */
        if(show_result && !analy->material_greyscale && !disable_gray)
        {
            grayel = 0;
            if((results_map[mesh_id][i+1] == '0') && !disable_mtl[matl] && !disable_flag && analy->auto_gray)
            {
                grayel = 1;         
                for(j = 0; j < 4; j++)
                {
                    for(k = 0; k < 4; k++)
                    {
                        cols[j][k] = grayval;
                    }
                }
               
            }
        } else if(!show_result && !analy->showmat && !disable_mtl[matl] && !disable_flag && analy->auto_gray && !disable_gray)
        {
            grayel = 1;
            for(j = 0; j < 4; j++)
            {
                for(k = 0; k < 4; k++)
                {
                    cols[j][k] = grayval;
                }
            }
        }

        hidden_poly = hide_by_object_type( p_quad_class, matl, i, analy, data_array );
 
        if(analy->mesh_view_mode == RENDER_WIREFRAMETRANS)
        {
            draw_edges_3d(analy);

        } else
        {
            draw_poly( cnt, verts, norms, cols, res, matl, p_mesh, analy, hidden_poly );
        }
    }

    grayel = 0;

    /* Back to the defaults. */
    end_draw_poly( analy );
    glLightModeli( GL_LIGHT_MODEL_TWO_SIDE, 0 );
    glDisable( GL_COLOR_MATERIAL );
    analy->interp_mode = save_interp_mode;
}


/************************************************************
 * TAG( draw_tris_3d )
 *
 * Draw the triangular elements in the model.
 */
static void
draw_tris_3d( Bool_type show_node_result, Bool_type show_mat_result,
              Bool_type show_mesh_result, MO_class_data *p_tri_class,
              Analysis *analy )
{
    Bool_type show_result, showgs=FALSE;
    float *activity;
    float verts[3][3];
    float norms[3][3];
    float cols[4][4];
    float res[3];
    float rmin, rmax;
    int matl, cnt, nd, i, j, k, mesh_idx;
    int (*connects)[3];
    int *mtl;
    unsigned char *disable_mtl, *hide_mtl;
    Mesh_data *p_mesh;
    Visibility_data *p_vd;
    int tri_qty;
    float *data_array;
    int *p_index_source;
    Interp_mode_type save_interp_mode;

/* Adding code to determine if a subrecord for these objects has the appropriate variable and if
 * not it sets that element to a gray scale prior to calling draw poly */
    int index, subrec, srec, qty_blocks;
    int total_blocks = 0;
    Subrec_obj *p_subrec;
    Bool_type hidden_poly;
    Bool_type disable_flag=FALSE;
    int mesh_id;
    static char **results_map;
    static char *result_name = NULL;
    static char *short_name = NULL;
    static int first = 0;
    Result * p_result;
    Bool_type same_as_last = FALSE;
    Bool_type disable_gray = TRUE;

    mesh_id = p_tri_class->mesh_id;
    p_result = analy->cur_result;

    show_result = show_node_result
                  || result_has_class( analy->cur_result, p_tri_class, analy )
                  || show_mat_result
                  || show_mesh_result;

    if(result_has_class( analy->cur_result, p_tri_class, analy))
    {
        disable_gray = FALSE;
    }

    if(first == 0)
    {
        results_map = (char **)malloc(analy->mesh_qty*sizeof(char *));
        if(results_map == NULL)
        {
            popup_dialog(WARNING_POPUP, "Out of memory in function draw_hexs. Terminating!!-n");
            parse_command("quit", analy);
        }
        for(i = 0; i < analy->mesh_qty; i++)
        {
            results_map[i] = NULL;
        } 

        first = 1;
    }
   
    if(show_result)
    { 
        if(result_name != NULL)
        {
            if(strcmp(result_name, p_result->name))
            {
                /* we are showing a different result than the last shown result */
                free(result_name);
                result_name = NULL;
                if(results_map[mesh_id] != NULL)
                {
                    free(results_map[mesh_id]);
                    results_map[mesh_id] = NULL;
                }
                same_as_last = FALSE;
            } else
            {
                same_as_last = TRUE;
            }
        }
        if(short_name != NULL )
        {
            if(strcmp(short_name, p_tri_class->short_name))
            {
                free(short_name);
                short_name = NULL;
                if(results_map[mesh_id] != NULL)
                {
                    free(results_map[mesh_id]);
                    results_map[mesh_id] = NULL;
                }
            }
        }  
       
        if(results_map[mesh_id] == NULL)
        {
            results_map[mesh_id] = (char *) malloc(p_tri_class->qty+1);
            if( results_map[mesh_id] == NULL)
            {
                popup_dialog(WARNING_POPUP, "Out of memory in function draw_hexs. Terminating!!\n");
                parse_command("quit", analy);
            }
        }
        if(result_name == NULL)
        {
            result_name = (char *) malloc(strlen(p_result->name) + 1);
            if(result_name == NULL)
            {
                popup_dialog(WARNING_POPUP, "Out of memory in function draw_hexs. Terminating!!\n");
                parse_command("quit", analy);
            }
            strncpy(result_name, p_result->name, strlen(p_result->name) + 1);
        }
        if(short_name == NULL)
        {
            short_name = (char *) malloc(strlen(p_tri_class->short_name) + 1);
            if(short_name == NULL)
            {
                popup_dialog(WARNING_POPUP, "Out of memory in function draw_hexs. Terminating!!\n");
                parse_command("quit", analy);
            }
            strncpy(short_name, p_tri_class->short_name, strlen(p_tri_class->short_name) + 1);
        }
    }

    if(show_result && analy->cur_result != NULL && !same_as_last) 
    { 
        populate_result(G_TRI, results_map[mesh_id], p_tri_class->qty + 1, p_tri_class, analy);
    }


    /*show_result = show_node_result
                  || result_has_class( analy->cur_result, p_tri_class, analy )
                  || show_mat_result
                  || show_mesh_result; */

    if ( analy->mesh_view_mode != RENDER_FILLED          && analy->mesh_view_mode != RENDER_WIREFRAME &&
            analy->mesh_view_mode != RENDER_WIREFRAMETRANS  && analy->mesh_view_mode != RENDER_HIDDEN )
        return;

    activity = analy->state_p->sand_present
               ? analy->state_p->elem_class_sand[p_tri_class->elem_class_index]
               : NULL;

    p_mesh = MESH_P( analy );
    connects = (int (*)[3]) p_tri_class->objects.elems->nodes;
    mtl = p_tri_class->objects.elems->mat;
    p_vd = p_tri_class->p_vis_data;
    tri_qty = p_tri_class->qty;
    disable_mtl = p_mesh->disable_material;
    hide_mtl = p_mesh->hide_material;

    /* Override interpolation mode for material and mesh results.  */
    save_interp_mode = analy->interp_mode;
    if ( show_mat_result || show_mesh_result )
        analy->interp_mode = NO_INTERP;

    /* Abstract the source array for data values and its index. */
    if ( analy->interp_mode == GOOD_INTERP )
    {
        data_array = NODAL_RESULT_BUFFER( analy );
        p_index_source = &nd;
    }
    else if ( show_mat_result )
    {
        data_array = ((MO_class_data **)
                      p_mesh->classes_by_sclass[G_MAT].list)[0]->data_buffer;
        p_index_source = &matl;
    }
    else if ( show_mesh_result )
    {
        mesh_idx = 0;
        data_array = ((MO_class_data **)
                      p_mesh->classes_by_sclass[G_MESH].list)[0]->data_buffer;
        p_index_source = &mesh_idx;
    }
    else if ( analy->interp_mode == NO_INTERP && !show_node_result )
    {
        data_array = p_tri_class->data_buffer;
        p_index_source = &i;
    }
    else
    {
        data_array = NODAL_RESULT_BUFFER( analy );
        p_index_source = &nd;
    }

    show_result = show_node_result
                  || result_has_class( analy->cur_result, p_tri_class, analy )
                  || show_mat_result
                  || show_mesh_result;

    get_min_max( analy, FALSE, &rmin, &rmax );

    /* Enable color to change AMBIENT & DIFFUSE property of the material. */
    glEnable( GL_COLOR_MATERIAL );

    /* Use two-sided lighting for shells. */
    glLightModeli( GL_LIGHT_MODEL_TWO_SIDE, 1 );

    /* Set up for polygon drawing. */
    begin_draw_poly( analy );

    /* Set vertex qty; if ever allow for degenerates, move into loop. */
    cnt = 3;

    for ( i = 0; i < tri_qty; i++ )
    {
        /* Check for inactive elements. */
        if ( activity )
        {
            if ( activity[i] != 0.0 && analy->show_only_deleted_elements )
                continue;
            if ( activity[i] == 0.0 && !analy->show_deleted_elements && !analy->show_only_deleted_elements )
                continue;
        }

        /* Skip if this material is invisible. */
        matl = mtl[i];
        if ( hide_mtl[matl] )
            continue;


        if ( v_win->mesh_materials.current_index != matl )
            change_current_color_property( &v_win->mesh_materials, matl );


        /* Colorflag is TRUE if displaying a result, otherwise
         * polygons are drawn in the material color.
         */
        colorflag = show_result && !disable_mtl[matl];
        if ( analy->material_greyscale && disable_mtl[matl] )
            showgs = TRUE;
        else
            showgs = FALSE;

        get_tri_verts_3d( i, p_tri_class, analy, verts );

        if ( analy->interp_mode == GOOD_INTERP )
        {
            for ( j = 0; j < 3; j++ )
            {
                nd = connects[i][j];

                for ( k = 0; k < 3; k++ )
                    norms[j][k] = p_vd->face_norm[j][k][i];

                res[j] = data_array[nd];
            }
        }
        else
        {
            for ( j = 0; j < 3; j++ )
            {
                nd = connects[i][j];

                for ( k = 0; k < 3; k++ )
                    norms[j][k] = p_vd->face_norm[j][k][i];

                color_lookup( cols[j], data_array[*p_index_source],
                              rmin, rmax, analy->zero_result, matl,
                              analy->logscale, showgs );
            }
        }

        disable_flag = disable_by_object_type( p_tri_class, matl, i, analy, data_array );

        /* check to see if this element has the variable to be shown */
        if(show_result && !analy->material_greyscale && !disable_gray)
        {
            grayel = 0;
            if((results_map[mesh_id][i+1] == '0') && !disable_mtl[matl] && !disable_flag && analy->auto_gray)
            {
                grayel = 1;         
                for(j = 0; j < 4; j++)
                {
                    for(k = 0; k < 4; k++)
                    {
                        cols[j][k] = grayval;
                    }
                }
               
            }
        } else if(!show_result && !analy->showmat &&!disable_mtl[matl] && !disable_flag && analy->auto_gray && !disable_gray)
        {
            grayel = 1;
            for(j = 0; j < 4; j++)
            {
                for(k = 0; k < 4; k++)
                {
                    cols[j][k] = grayval;
                }
            }
             
        }

        hidden_poly = FALSE;
        draw_poly( cnt, verts, norms, cols, res, matl, p_mesh, analy, hidden_poly );
    }

    /* Back to the defaults. */
    end_draw_poly( analy );
    glLightModeli( GL_LIGHT_MODEL_TWO_SIDE, 0 );
    glDisable( GL_COLOR_MATERIAL );
    analy->interp_mode = save_interp_mode;
}

/*****************************************************************
 * TAG( bad_node_warn_once )
 *
 * Controls whether a warning about a bad beam node is given.
 * Some meshes will incur this warning over many elements and
 * the warning becomes very redundant and time-consuming.
 */
static Bool_type bad_node_warn_once=TRUE;


/************************************************************
 * TAG( draw_beams_3d )
 *
 * Draw the beam elements in the model.
 */
static void
draw_beams_3d( Bool_type show_node_result, Bool_type show_mat_result, Bool_type show_mesh_result,
               MO_class_data *p_beam_class, Analysis *analy )
{
    Bool_type verts_ok, show_result, showgs;
    float *data_array;
    float *activity;
    float col[4];
    float verts[2][3];
    float pts[6];
    float threshold, rmin, rmax;
    GLfloat nearfar[2];
    int matl, i, j, k, mesh_idx;
    int *mtl;
    int *p_index_source;
    int beam_qty;
    unsigned char *disable_mtl, *hide_mtl;
    Bool_type disable_flag=FALSE;
    Mesh_data *p_mesh;

    int mesh_id;
    static int first = 0;
    static char ** results_map;
    static char * result_name = NULL;
    static char * short_name = NULL;
    Result * p_result;
    Bool_type same_as_last = FALSE;
    Bool_type disable_gray = TRUE;
    mesh_id = p_beam_class->mesh_id;

    p_result = analy->cur_result;

    show_result = result_has_class( analy->cur_result, p_beam_class, analy )
                  || show_node_result
                  || show_mat_result
                  || show_mesh_result;
   
    if(result_has_class( analy->cur_result, p_beam_class, analy ))
    {
        disable_gray = FALSE;
    }

    if(first == 0)
    {
        results_map = (char **)malloc(analy->mesh_qty*sizeof(char *));
        if(results_map == NULL)
        {
            popup_dialog(WARNING_POPUP, "Out of memory in function draw_hexs. Terminating!!-n");
            parse_command("quit", analy);
        }
        for(i = 0; i < analy->mesh_qty; i++)
        {
            results_map[i] = NULL;
        } 

        first = 1;
    }
   
    if(show_result)
    { 
        if(result_name != NULL)
        {
            if(strcmp(result_name, p_result->name))
            {
                /* we are showing a different result than the last shown result */
                free(result_name);
                result_name = NULL;
                if(results_map[mesh_id] != NULL)
                {
                    free(results_map[mesh_id]);
                    results_map[mesh_id] = NULL;
                }
                same_as_last = FALSE;
            } else
            {
                same_as_last = TRUE;
            }
        }
        if(short_name != NULL )
        {
            if(strcmp(short_name, p_beam_class->short_name))
            {
                free(short_name);
                short_name = NULL;
                if(results_map[mesh_id] != NULL)
                {
                    free(results_map[mesh_id]);
                    results_map[mesh_id] = NULL;
                }
            }
        }  
       
        if(results_map[mesh_id] == NULL)
        {
            results_map[mesh_id] = (char *) malloc(p_beam_class->qty+1);
            if( results_map[mesh_id] == NULL)
            {
                popup_dialog(WARNING_POPUP, "Out of memory in function draw_beam_3d. Terminating!!\n");
                parse_command("quit", analy);
            }
        }
        if(result_name == NULL)
        {
            result_name = (char *) malloc(strlen(p_result->name) + 1);
            if(result_name == NULL)
            {
                popup_dialog(WARNING_POPUP, "Out of memory in function draw_beam_3d. Terminating!!\n");
                parse_command("quit", analy);
            }
            strncpy(result_name, p_result->name, strlen(p_result->name) + 1);
        }
        if(short_name == NULL)
        {
            short_name = (char *) malloc(strlen(p_beam_class->short_name) + 1);
            if(short_name == NULL)
            {
                popup_dialog(WARNING_POPUP, "Out of memory in function draw_beam_3d. Terminating!!\n");
                parse_command("quit", analy);
            }
            strncpy(short_name, p_beam_class->short_name, strlen(p_beam_class->short_name) + 1);
        }
    }

    if(show_result && analy->cur_result != NULL && !same_as_last) 
    { 
        populate_result(G_BEAM, results_map[mesh_id], p_beam_class->qty + 1, p_beam_class, analy);
    }


    if ( analy->mesh_view_mode != RENDER_FILLED          && analy->mesh_view_mode != RENDER_WIREFRAME &&
            analy->mesh_view_mode != RENDER_WIREFRAMETRANS  && analy->mesh_view_mode != RENDER_HIDDEN )
        return;

    p_mesh = MESH_P( analy );
    mtl = p_beam_class->objects.elems->mat;
    beam_qty = p_beam_class->qty;
    disable_mtl = p_mesh->disable_material;
    hide_mtl = p_mesh->hide_material;

    /* Abstract the source array for data values and its index. */
    if ( show_mat_result )
    {
        data_array = ((MO_class_data **)
                      p_mesh->classes_by_sclass[G_MAT].list)[0]->data_buffer;
        p_index_source = &matl;
    }
    else if ( show_mesh_result )
    {
        mesh_idx = 0;
        data_array = ((MO_class_data **)
                      p_mesh->classes_by_sclass[G_MESH].list)[0]->data_buffer;
        p_index_source = &mesh_idx;
    }
    else
    {
        data_array = p_beam_class->data_buffer;
        p_index_source = &i;
    }

    /*show_result = result_has_class( analy->cur_result, p_beam_class, analy )
                  || show_mat_result
                  || show_mesh_result; */

    activity = analy->state_p->sand_present
               ?  analy->state_p->elem_class_sand[p_beam_class->elem_class_index]
               : NULL;
    get_min_max( analy, TRUE, &rmin, &rmax );
    threshold = analy->zero_result;

    antialias_lines( TRUE, analy->z_buffer_lines );
    glLineWidth( 3 );

    /* Bias the depth-buffer so beams are in front. */
    if ( analy->zbias_beams )
    {
        glGetFloatv( GL_DEPTH_RANGE, nearfar );
        glDepthRange( 0, 1 - analy->beam_zbias );
    }

    for ( i = 0; i < beam_qty; i++ )
    {
        /* Check for inactive elements. */
        if ( activity )
            if ( activity[i] == 0.0 && !analy->show_deleted_elements )
                continue;
        if ( activity )
            if ( activity[i] == 0.0 && !analy->show_deleted_elements )
                continue;

        /* Skip if this material is invisible. */
        matl = mtl[i];
        if ( hide_mtl[matl] || hide_by_object_type( p_beam_class, matl, i, analy, data_array ))
            continue;

        verts_ok = get_beam_verts_2d_3d( i, p_beam_class, analy, verts );
        if ( bad_node_warn_once && !verts_ok )
        {
            popup_dialog( WARNING_POPUP,
                          "Warning - beam %d ignored; bad node(s)", i + 1 );
            bad_node_warn_once = FALSE;
            continue;
        }

        colorflag = show_result && !disable_mtl[matl] &&
                    !disable_by_object_type( p_beam_class, matl, i, analy, data_array );
        if ( analy->material_greyscale && disable_mtl[matl] )
            showgs = TRUE;
        else
            showgs = FALSE;

        color_lookup( col,  data_array[*p_index_source],
                      rmin, rmax, threshold,
                      matl, analy->logscale, showgs );

        glColor3fv(col);

        disable_flag = disable_by_object_type(p_beam_class, matl, i, analy, data_array);

        grayel = 0;

        if(show_result && !analy->material_greyscale && !disable_gray)
        {
            grayel = 1;

            if((results_map[mesh_id][i + 1] == '0') && !disable_mtl[matl] && !disable_flag && analy->auto_gray)
            {   
                for(j = 0; j < 4; j++)
                {
                    col[j] = 0.7;
                }
                
                glColor3fv(col);
            }
        } else if(!show_result && !analy->showmat && !disable_mtl[matl] && !disable_flag && analy->auto_gray && !disable_gray)
        {
            for(j = 0; j < 4; j++)
            {
                col[j] = 0.7;
            }
             
            glColor3fv(col);
        }


        for ( j = 0; j < 2; j++ )
            for ( k = 0; k < 3; k++ )
                pts[j*3+k] = verts[j][k];
        draw_line( 2, pts, matl, p_mesh, analy, FALSE, NULL );
    }
    
    grayel = 0;

    /* Remove depth bias. */
    if ( analy->zbias_beams )
        glDepthRange( nearfar[0], nearfar[1] );

    antialias_lines( FALSE, 0 );
    glLineWidth( 1.0 );
}

/************************************************************
 * TAG( draw_truss_3d )
 *
 * Draw the beam elements in the model.
 */
static void
draw_truss_3d( Bool_type show_mat_result, Bool_type show_mesh_result,
               MO_class_data *p_truss_class, Analysis *analy )
{
    Bool_type verts_ok, show_result, showgs=FALSE;
    float *data_array;
    float *activity;
    float col[4];
    float verts[2][3];
    float pts[6];
    float threshold, rmin, rmax;
    GLfloat nearfar[2];
    int matl, i, j, k, mesh_idx;
    int *mtl;
    int *p_index_source;
    int truss_qty;
    unsigned char *disable_mtl, *hide_mtl;
    Mesh_data *p_mesh;
    Bool_type disable_flag=FALSE;
   
    int grayel = 0; 
    int mesh_id;
    static int first = 0;
    static char **results_map;
    static char * result_name = NULL;
    static char * short_name = NULL;
    Bool_type same_as_last = FALSE;
    Bool_type disable_gray = TRUE;
    Result * p_result;

    mesh_id = p_truss_class->mesh_id;
    p_result = analy->cur_result; 

    show_result = result_has_class( analy->cur_result, p_truss_class, analy )
                  || show_mat_result
                  || show_mesh_result;

    if(result_has_class( analy->cur_result, p_truss_class, analy ))
    {
        disable_gray = FALSE;
    }

    if(first == 0)
    {
        results_map = (char **)malloc(analy->mesh_qty*sizeof(char *));
        if(results_map == NULL)
        {
            popup_dialog(WARNING_POPUP, "Out of memory in function draw_hexs. Terminating!!-n");
            parse_command("quit", analy);
        }
        for(i = 0; i < analy->mesh_qty; i++)
        {
            results_map[i] = NULL;
        } 

        first = 1;
    }
   
    if(show_result)
    { 
        if(result_name != NULL)
        {
            if(strcmp(result_name, p_result->name))
            {
                /* we are showing a different result than the last shown result */
                free(result_name);
                result_name = NULL;
                if(results_map[mesh_id] != NULL)
                {
                    free(results_map[mesh_id]);
                    results_map[mesh_id] = NULL;
                }
                same_as_last = FALSE;
            } else
            {
                same_as_last = TRUE;
            }
        }
        if(short_name != NULL )
        {
            if(strcmp(short_name, p_truss_class->short_name))
            {
                free(short_name);
                short_name = NULL;
                if(results_map[mesh_id] != NULL)
                {
                    free(results_map[mesh_id]);
                    results_map[mesh_id] = NULL;
                }
            }
        }  
       
        if(results_map[mesh_id] == NULL)
        {
            results_map[mesh_id] = (char *) malloc(p_truss_class->qty+1);
            if( results_map[mesh_id] == NULL)
            {
                popup_dialog(WARNING_POPUP, "Out of memory in function draw_hexs. Terminating!!\n");
                parse_command("quit", analy);
            }
        }
        if(result_name == NULL)
        {
            result_name = (char *) malloc(strlen(p_result->name) + 1);
            if(result_name == NULL)
            {
                popup_dialog(WARNING_POPUP, "Out of memory in function draw_hexs. Terminating!!\n");
                parse_command("quit", analy);
            }
            strncpy(result_name, p_result->name, strlen(p_result->name) + 1);
        }
        if(short_name == NULL)
        {
            short_name = (char *) malloc(strlen(p_truss_class->short_name) + 1);
            if(short_name == NULL)
            {
                popup_dialog(WARNING_POPUP, "Out of memory in function draw_hexs. Terminating!!\n");
                parse_command("quit", analy);
            }
            strncpy(short_name, p_truss_class->short_name, strlen(p_truss_class->short_name) + 1);
        }
    }

    if(show_result && analy->cur_result != NULL && !same_as_last) 
    { 
        populate_result(G_TRUSS, results_map[mesh_id], p_truss_class->qty + 1, p_truss_class, analy);
    }


    if ( analy->mesh_view_mode != RENDER_FILLED          && analy->mesh_view_mode != RENDER_WIREFRAME &&
            analy->mesh_view_mode != RENDER_WIREFRAMETRANS  && analy->mesh_view_mode != RENDER_HIDDEN )
        return;

    p_mesh = MESH_P( analy );
    mtl = p_truss_class->objects.elems->mat;
    truss_qty = p_truss_class->qty;
    disable_mtl = p_mesh->disable_material;
    hide_mtl = p_mesh->hide_material;

    /* Abstract the source array for data values and its index. */
    if ( show_mat_result )
    {
        data_array = ((MO_class_data **)
                      p_mesh->classes_by_sclass[G_MAT].list)[0]->data_buffer;
        p_index_source = &matl;
    }
    else if ( show_mesh_result )
    {
        mesh_idx = 0;
        data_array = ((MO_class_data **)
                      p_mesh->classes_by_sclass[G_MESH].list)[0]->data_buffer;
        p_index_source = &mesh_idx;
    }
    else
    {
        data_array = p_truss_class->data_buffer;
        p_index_source = &i;
    }

    /*show_result = result_has_class( analy->cur_result, p_truss_class, analy )
                  || show_mat_result
                  || show_mesh_result; */

    activity = analy->state_p->sand_present
               ?  analy->state_p->elem_class_sand[p_truss_class->elem_class_index]
               : NULL;
    get_min_max( analy, TRUE, &rmin, &rmax );
    threshold = analy->zero_result;

    antialias_lines( TRUE, analy->z_buffer_lines );
    glLineWidth( 3 );

    /* Bias the depth-buffer so trusses  are in front. */
    if ( analy->zbias_beams )
    {
        glGetFloatv( GL_DEPTH_RANGE, nearfar );
        glDepthRange( 0, 1 - analy->beam_zbias );
    }

    for ( i = 0; i < truss_qty; i++ )
    {
        /* Check for inactive elements. */
        if ( activity )
        {
            if ( activity[i] != 0.0 && analy->show_only_deleted_elements )
                continue;
            if ( activity[i] == 0.0 && !analy->show_deleted_elements && !analy->show_only_deleted_elements )
                continue;
        }

        /* Skip if this material is invisible. */
        matl = mtl[i];
        if ( hide_mtl[matl] || hide_by_object_type( p_truss_class, matl, i, analy, data_array ))
            continue;

        verts_ok = get_truss_verts_2d_3d( i, p_truss_class, analy, verts );
        if ( !verts_ok )
        {
            popup_dialog( WARNING_POPUP,
                          "Warning - Discrete Element %d ignored; bad node(s)",
                          i + 1 );
            continue;
        }

        colorflag = show_result && !disable_mtl[matl] &&
                    !disable_by_object_type( p_truss_class, matl, i, analy, data_array );
        if ( analy->material_greyscale && disable_mtl[matl] )
            showgs = TRUE;
        else
            showgs = FALSE;

        color_lookup( col,  data_array[*p_index_source],
                      rmin, rmax, threshold,
                      matl, analy->logscale, showgs );

        glColor3fv( col );

        for ( j = 0; j < 2; j++ )
            for ( k = 0; k < 3; k++ )
                pts[j*3+k] = verts[j][k];

        if(show_result && !analy->material_greyscale && !disable_gray)
        {
            grayel = 0;
            if((results_map[mesh_id][i+1] == '0') && !disable_mtl[matl] && !disable_flag && analy->auto_gray)
            {
                grayel = 1;         
                for(k = 0; k < 4; k++)
                {
                    col[k] = grayval;
                }
                glColor3fv( col );
            }
        } else if(!show_result && !analy->showmat &&!disable_mtl[matl] && !disable_flag && analy->auto_gray && !disable_gray)
        {
            grayel = 1;
            for(k = 0; k < 4; k++)
            {
                col[k] = grayval;
            }
            glColor3fv( col );
        }
        draw_line( 2, pts, matl, p_mesh, analy, FALSE, NULL );
    }

    /* Remove depth bias. */
    if ( analy->zbias_beams )
        glDepthRange( nearfar[0], nearfar[1] );

    antialias_lines( FALSE, 0 );
    glLineWidth( 1.0 );
}




/************************************************************
 * TAG( draw_pyramids )
 *
 * Draw the external faces of pyramid volume elements in the model.
 */
static void
draw_pyramids( Bool_type show_node_result, Bool_type show_mat_result,
               Bool_type show_mesh_result, MO_class_data *p_pyramid_class,
               Analysis *analy )
{
    Bool_type show_result;
    float verts[4][3];
    float norms[4][3];
    float cols[4][4];
    float res[4];  /* "val" was formerly declared but never referenced  LAS */
    float rmin, rmax;
    float tverts[4][3], v1[3], v2[3], nor[3], ray[3];
    float *activity=NULL;

    int matl, el, fc, nd, ord[4], cnt, face_qty, mesh_idx;
    int i, j, k;
    int (*connects)[5];
    int *face_el, *face_fc;
    int *mtl;
    int *p_index_source;
    float *data_array;
    unsigned char *disable_mtl;
    Mesh_data *p_mesh;
    Visibility_data *p_vd;
    Interp_mode_type save_interp_mode;
    unsigned char *hide_mtl;
    int *p_mats;

    Bool_type hidden_poly=FALSE,
              hidden_poly_elem=FALSE,
              hidden_poly_elem_wft=FALSE,
              hidden_poly_mat=FALSE,
              disable_flag=FALSE;
    Bool_type showgs=TRUE;

    if ( analy->mesh_view_mode != RENDER_FILLED          && analy->mesh_view_mode != RENDER_WIREFRAME &&
            analy->mesh_view_mode != RENDER_WIREFRAMETRANS  && analy->mesh_view_mode != RENDER_HIDDEN )
        return;

    if ( is_particle_class(analy, p_pyramid_class->superclass, p_pyramid_class->short_name) )
        return;

    p_mesh = MESH_P( analy );
    connects = (int (*)[5]) p_pyramid_class->objects.elems->nodes;
    mtl = p_pyramid_class->objects.elems->mat;
    p_vd = p_pyramid_class->p_vis_data;

    /* until the visibality data is there just return */
    /*return; */
    if(p_vd == NULL)
    {
        return;
    }
    face_qty = p_vd->face_cnt;
    face_el = p_vd->face_el;
    face_fc = p_vd->face_fc;
    disable_mtl = p_mesh->disable_material;
    hide_mtl    = p_mesh->hide_material;
    p_mats      = p_mesh->particle_mats;

    /* Override interpolation mode for material and mesh results. */
    save_interp_mode = analy->interp_mode;
    if ( show_mat_result || show_mesh_result )
        analy->interp_mode = NO_INTERP;

    activity = analy->state_p->sand_present
               ? analy->state_p->elem_class_sand[p_pyramid_class->elem_class_index]
               : NULL;

    /* Abstract the source array for data values and its index. */
    if ( analy->interp_mode == GOOD_INTERP )
    {
        data_array = NODAL_RESULT_BUFFER( analy );
        p_index_source = &nd;
    }
    else if ( show_mat_result )
    {
        data_array = ((MO_class_data **)
                      p_mesh->classes_by_sclass[G_MAT].list)[0]->data_buffer;
        p_index_source = &matl;
    }
    else if ( show_mesh_result )
    {
        mesh_idx = 0;
        data_array = ((MO_class_data **)
                      p_mesh->classes_by_sclass[G_MESH].list)[0]->data_buffer;
        p_index_source = &mesh_idx;
    }
    else if ( analy->interp_mode == NO_INTERP && !show_node_result )
    {
        data_array = p_pyramid_class->data_buffer;
        p_index_source = &el;
    }
    else
    {
        data_array = NODAL_RESULT_BUFFER( analy );
        p_index_source = &nd;
    }

    show_result = show_node_result
                  || result_has_class( analy->cur_result, p_pyramid_class, analy )
                  || show_mat_result
                  || show_mesh_result;

    /* If we are drawing particles then look for a particle result and map onto hexes
      *
      */
    /* if ( show_result && analy->free_particles )
           data_array = analy->pn_buffer_ptr[0]; */

    get_min_max( analy, FALSE, &rmin, &rmax );

    /* Enable color to change AMBIENT & DIFFUSE property of the material. */
    glEnable( GL_COLOR_MATERIAL );

    /* Throw away backward-facing faces. */
    glEnable( GL_CULL_FACE );

    /* Set up for polygon drawing. */
    begin_draw_poly( analy );

    for ( i = 0; i < face_qty; i++ )
    {
        el = face_el[i];
        fc = face_fc[i];

        /*
         * Remove faces that are shared with quad elements, so
         * the polygons are not drawn twice.
         */
        if ( analy->shared_faces
                && p_mesh->classes_by_sclass[G_QUAD].qty > 0
                && face_matches_quad( el, fc, p_pyramid_class, p_mesh, analy ) )
            continue;

        /* Check for inactive elements. */
        if ( activity )
        {
            if ( activity[el] != 0.0 && analy->show_only_deleted_elements )
                continue;
            if ( activity[el] == 0.0 && !analy->show_deleted_elements && !analy->show_only_deleted_elements )
                continue;
        }

        get_pyramid_face_verts( el, fc, p_pyramid_class, analy, verts );

        /* Cull backfacing polygons by hand to speed things up. */
        if ( !analy->reflect && analy->manual_backface_cull )
        {
            for ( j = 0; j < 4; j++ )
                point_transform( tverts[j], verts[j], &cur_view_mat );
            VEC_SUB( v1, tverts[2], tverts[0] );
            VEC_SUB( v2, tverts[3], tverts[1] );
            VEC_CROSS( nor, v1, v2 );
            if ( v_win->orthographic )
            {
                VEC_SET( ray, 0.0, 0.0, -1.0 );
            }
            else
            {
                for ( k = 0; k < 3; k++ )
                    ray[k] = 0.25 * ( tverts[0][k] + tverts[1][k] +
                                      tverts[2][k] + tverts[3][k] );
            }

            if ( VEC_DOT( ray, nor ) >= 0.0 )
                continue;
        }

        /*
         * Check for triangular (degenerate) face and reorder nodes
         * if needed.
         */
        cnt = check_for_tri_face( fc, el, p_pyramid_class, verts, ord );

        matl = mtl[el];

        if ( v_win->mesh_materials.current_index != matl )
            change_current_color_property( &v_win->mesh_materials, matl );


        /* Colorflag is TRUE if displaying a result, otherwise
         * polygons are drawn in the material color.
         */
        disable_flag = disable_by_object_type( p_pyramid_class, matl, el, analy, data_array );

        colorflag = show_result && !disable_mtl[matl] && !disable_by_object_type( p_pyramid_class, matl, el, analy, data_array );

        colorflag = show_result && !disable_mtl[matl];
        colorflag = show_result && !disable_mtl[matl] && !disable_by_object_type( p_pyramid_class, matl, el, analy, data_array );

        if ( analy->material_greyscale && disable_mtl[matl])
            showgs = TRUE;
        else
            showgs = FALSE;

        /*
        * Array "ord[4]" represents a map from the order of vertex data
        * specified by the per-face node ordering to the order needed for
        * rendering.  For quad's, there is no difference (i.e., "ord"
        * contains [0, 1, 2, 3].  For tri's, which are actually degenerate
        * quads with one node repeated, we want the repeated node to be last
        * as we will only render the first three, so there is potentially
        * a re-ordering.  This requires all the associated data to be re-
        * ordered, so the data arrays passed to draw_poly() all have their
        * data stored taking the original data at position set by "ord[j]" and
        * storing it in position "j".
        */
        if ( analy->interp_mode == GOOD_INTERP )
        {
            for ( j = 0; j < 4; j++ )
            {
                nd = connects[el][ fc_nd_nums[fc][ord[j]] ];

                for ( k = 0; k < 3; k++ )
                    norms[j][k] = p_vd->face_norm[ ord[j] ][k][i];

                res[j] = data_array[nd];
            }
        }
        else
        {
            for ( j = 0; j < 4; j++ )
            {
                nd = connects[el][ pyramid_fc_nd_nums[fc][ord[j]] ];

                for ( k = 0; k < 3; k++ )
                    norms[j][k] = p_vd->face_norm[ ord[j] ][k][i];

                color_lookup( cols[j], data_array[*p_index_source],
                              rmin, rmax, analy->zero_result, matl,
                              analy->logscale, showgs );
            }
        }

        hidden_poly_mat = hide_mtl[matl];

        hidden_poly_elem = hide_by_object_type( p_pyramid_class, matl, el, analy, data_array );
        if ( analy->mesh_view_mode == RENDER_WIREFRAMETRANS )
            hidden_poly_elem_wft = disable_by_object_type( p_pyramid_class, matl, el, analy, data_array );

        if ( hidden_poly_mat || hidden_poly_elem || hidden_poly_elem_wft )
            hidden_poly = TRUE;
        else
            hidden_poly = FALSE;

        if ( analy->particle_nodes_hide_background && p_mats )
            if ( p_mats[matl] )
                hidden_poly = TRUE;

        draw_poly( cnt, verts, norms, cols, res, matl, p_mesh, analy, hidden_poly );
    }

    /* Back to the defaults. */
    end_draw_poly( analy );
    glDisable( GL_CULL_FACE );
    glDisable( GL_COLOR_MATERIAL );
    analy->interp_mode = save_interp_mode;
}


/************************************************************
 * TAG( draw_quads_2d )
 *
 * Draw the quad elements for a class in the model.
 */
static void
draw_quads_2d( Bool_type show_node_result, Bool_type show_mat_result,
               Bool_type show_mesh_result, MO_class_data *p_quad_class,
               Analysis *analy )
{
    int i, j, k, nd, matl;
    Mesh_data *p_mesh;
    int (*connects)[4];
    GVec2D *nodes, *onodes;
    Bool_type show_result, showgs;
    float *activity;
    float orig;
    float rmin, rmax;
    float verts[4][3];
    float cols[4][4];
    float pts[12];
    float res[4];
    int quad_qty;
    unsigned char *hide_mtl, *disable_mtl;
    int *mtls;
    float *data_array;
    int *p_index_source;
    int mesh_idx;
    Interp_mode_type save_interp_mode;

    if ( analy->mesh_view_mode != RENDER_FILLED          && analy->mesh_view_mode != RENDER_WIREFRAME &&
            analy->mesh_view_mode != RENDER_WIREFRAMETRANS  && analy->mesh_view_mode != RENDER_HIDDEN )
        return;

    p_mesh = analy->mesh_table + p_quad_class->mesh_id;
    hide_mtl = p_mesh->hide_material;
    disable_mtl = p_mesh->disable_material;
    mtls = p_quad_class->objects.elems->mat;
    quad_qty = p_quad_class->qty;
    connects = (int (*)[4]) p_quad_class->objects.elems->nodes;
    nodes = analy->state_p->nodes.nodes2d;
    onodes = (GVec2D *) analy->cur_ref_state_data;
    activity = analy->state_p->sand_present
               ? analy->state_p->elem_class_sand[p_quad_class->elem_class_index]
               : NULL;

    /* 2D, so set Z coord of every vertex to zero. */
    verts[0][2] = verts[1][2] = verts[2][2] = verts[3][2] = 0.0;

    /* Override interpolation mode for material and mesh results.  */
    save_interp_mode = analy->interp_mode;
    if ( show_mat_result || show_mesh_result )
        analy->interp_mode = NO_INTERP;

    /* Abstract the source array for data values and its index. */
    if ( analy->interp_mode == GOOD_INTERP )
    {
        data_array = NODAL_RESULT_BUFFER( analy );
        p_index_source = &nd;
    }
    else if ( show_mat_result )
    {
        data_array = ((MO_class_data **)
                      p_mesh->classes_by_sclass[G_MAT].list)[0]->data_buffer;
        p_index_source = &matl;
    }
    else if ( show_mesh_result )
    {
        mesh_idx = 0;
        data_array = ((MO_class_data **)
                      p_mesh->classes_by_sclass[G_MESH].list)[0]->data_buffer;
        p_index_source = &mesh_idx;
    }
    else if ( analy->interp_mode == NO_INTERP && !show_node_result )
    {
        data_array = p_quad_class->data_buffer;
        p_index_source = &i;
    }
    else
    {
        data_array = NODAL_RESULT_BUFFER( analy );
        p_index_source = &nd;
    }

    show_result = show_node_result
                  || result_has_class( analy->cur_result, p_quad_class, analy )
                  || show_mat_result
                  || show_mesh_result;

    get_min_max( analy, FALSE, &rmin, &rmax );

    /* Set up for polygon drawing. */
    begin_draw_poly( analy );

    /*
     * Draw filled polygons for the tri elements.
     */
    for ( i = 0; i < quad_qty; i++ )
    {
        /* Check for inactive elements. */
        if ( activity )
        {
            if ( activity[i] != 0.0 && analy->show_only_deleted_elements )
                continue;
            if ( activity[i] == 0.0 && !analy->show_deleted_elements && !analy->show_only_deleted_elements )
                continue;
        }

        matl = mtls[i];

        /* Skip if this material is invisible. */
        if ( hide_mtl[matl] || hide_by_object_type( p_quad_class, matl, i, analy, data_array ))
            continue;

        /* If displaying a result, enable color to change
         * the DIFFUSE property of the material.
         */
        colorflag = show_result && !disable_mtl[matl] &&
                    !disable_by_object_type( p_quad_class, matl, i, analy, data_array );
        if ( analy->material_greyscale && disable_mtl[matl] )
            showgs = TRUE;
        else
            showgs = FALSE;

        if ( analy->interp_mode == GOOD_INTERP )
        {
            if ( analy->displace_exag )
            {
                for ( j = 0; j < 4; j++ )
                {
                    nd = connects[i][j];
                    for ( k = 0; k < 2; k++ )
                    {
                        orig = onodes[nd][k];
                        verts[j][k] = orig + analy->displace_scale[k]
                                      * (nodes[nd][k] - orig);
                    }

                    res[j] = data_array[nd];
                }
            }
            else
            {
                for ( j = 0; j < 4; j++ )
                {
                    nd = connects[i][j];
                    for ( k = 0; k < 2; k++ )
                        verts[j][k] = nodes[nd][k];

                    res[j] = data_array[nd];
                }
            }
        }
        else
        {
            if ( analy->displace_exag )
            {
                for ( j = 0; j < 4; j++ )
                {
                    nd = connects[i][j];
                    for ( k = 0; k < 2; k++ )
                    {
                        orig = onodes[nd][k];
                        verts[j][k] = orig + analy->displace_scale[k]
                                      * (nodes[nd][k] - orig);
                    }

                    color_lookup( cols[j], data_array[*p_index_source],
                                  rmin, rmax, analy->zero_result, matl,
                                  analy->logscale, showgs );
                }
            }
            else
            {
                for ( j = 0; j < 4; j++ )
                {
                    nd = connects[i][j];
                    for ( k = 0; k < 2; k++ )
                        verts[j][k] = nodes[nd][k];

                    color_lookup( cols[j], data_array[*p_index_source],
                                  rmin, rmax, analy->zero_result, matl,
                                  analy->logscale, showgs );
                }
            }
        }

        draw_poly_2d( 4, verts, cols, res, matl, p_mesh, analy );
    }

    end_draw_poly( analy );

    /*
     * Draw the outlines of the shell elements.
     */
    if ( analy->mesh_view_mode == RENDER_HIDDEN || analy->mesh_view_mode == RENDER_HIDDEN || analy->mesh_view_mode == RENDER_WIREFRAMETRANS ||
            analy->mesh_view_mode != RENDER_WIREFRAME )
    {
        for ( i = 2; i < 12; i += 3 )
            pts[i] = 0.0;

        glColor3fv( v_win->mesh_color );

        for ( i = 0; i < quad_qty; i++ )
        {
            /* Check for inactive elements. */
            if ( activity )
            {
                if ( activity[i] != 0.0 && analy->show_only_deleted_elements )
                    continue;
                if ( activity[i] == 0.0 && !analy->show_deleted_elements && !analy->show_only_deleted_elements )
                    continue;
            }

            /* Skip if this material is invisible. */
            if ( hide_mtl[mtls[i]] || hide_by_object_type( p_quad_class, matl, i, analy, data_array ))
                continue;

            if ( analy->displace_exag )
            {
                for ( j = 0; j < 4; j++ )
                {
                    nd = connects[i][j];
                    for ( k = 0; k < 2; k++ )
                    {
                        orig = onodes[nd][k];
                        pts[j*3+k] = orig + analy->displace_scale[k]*
                                     (nodes[nd][k] - orig);
                    }
                }
            }
            else
            {
                for ( j = 0; j < 4; j++ )
                {
                    nd = connects[i][j];
                    for ( k = 0; k < 2; k++ )
                        pts[j*3+k] = nodes[nd][k];
                }
            }

            draw_line( 4, pts, mtls[i], p_mesh, analy, TRUE, NULL );
        }
    }

    analy->interp_mode = save_interp_mode;
}


/************************************************************
 * TAG( draw_tris_2d )
 *
 * Draw the tri elements for a class in the model.
 */
static void
draw_tris_2d( Bool_type show_node_result, Bool_type show_mat_result,
              Bool_type show_mesh_result, MO_class_data *p_tri_class,
              Analysis *analy )
{
    int i, j, k, nd, matl;
    Mesh_data *p_mesh;
    int (*connects)[3];
    GVec2D *nodes, *onodes;
    Bool_type show_result, showgs=FALSE;
    float *activity;
    float orig;
    float rmin, rmax;
    float verts[3][3];
    float cols[4][4];
    float pts[9];
    float res[4];
    int tri_qty;
    unsigned char *hide_mtl, *disable_mtl;
    int *mtls;
    float *data_array;
    int *p_index_source;
    int mesh_idx;
    Interp_mode_type save_interp_mode;

    if ( analy->mesh_view_mode != RENDER_FILLED          && analy->mesh_view_mode != RENDER_WIREFRAME &&
            analy->mesh_view_mode != RENDER_WIREFRAMETRANS  && analy->mesh_view_mode != RENDER_HIDDEN )
        return;

    p_mesh = analy->mesh_table + p_tri_class->mesh_id;
    hide_mtl = p_mesh->hide_material;
    disable_mtl = p_mesh->disable_material;
    mtls = p_tri_class->objects.elems->mat;
    tri_qty = p_tri_class->qty;
    connects = (int (*)[3]) p_tri_class->objects.elems->nodes;
    nodes = analy->state_p->nodes.nodes2d;
    onodes = (GVec2D *) analy->cur_ref_state_data;
    activity = analy->state_p->sand_present
               ? analy->state_p->elem_class_sand[p_tri_class->elem_class_index]
               : NULL;

    /* 2D, so set Z coord of every vertex to zero. */
    verts[0][2] = verts[1][2] = verts[2][2] = 0.0;

    /* Override interpolation mode for material and mesh results.  */
    save_interp_mode = analy->interp_mode;
    if ( show_mat_result || show_mesh_result )
        analy->interp_mode = NO_INTERP;

    /* Abstract the source array for data values and its index. */
    if ( analy->interp_mode == GOOD_INTERP )
    {
        data_array = NODAL_RESULT_BUFFER( analy );
        p_index_source = &nd;
    }
    else if ( show_mat_result )
    {
        data_array = ((MO_class_data **)
                      p_mesh->classes_by_sclass[G_MAT].list)[0]->data_buffer;
        p_index_source = &matl;
    }
    else if ( show_mesh_result )
    {
        mesh_idx = 0;
        data_array = ((MO_class_data **)
                      p_mesh->classes_by_sclass[G_MESH].list)[0]->data_buffer;
        p_index_source = &mesh_idx;
    }
    else if ( analy->interp_mode == NO_INTERP && !show_node_result )
    {
        data_array = p_tri_class->data_buffer;
        p_index_source = &i;
    }
    else
    {
        data_array = NODAL_RESULT_BUFFER( analy );
        p_index_source = &nd;
    }

    show_result = show_node_result
                  || result_has_class( analy->cur_result, p_tri_class, analy )
                  || show_mat_result
                  || show_mesh_result;

    get_min_max( analy, FALSE, &rmin, &rmax );

    /* Set up for polygon drawing. */
    begin_draw_poly( analy );

    /*
     * Draw filled polygons for the tri elements.
     */
    for ( i = 0; i < tri_qty; i++ )
    {
        /* Check for inactive elements. */
        if ( activity )
        {
            if ( activity[i] != 0.0 && analy->show_only_deleted_elements )
                continue;
            if ( activity[i] == 0.0 && !analy->show_deleted_elements && !analy->show_only_deleted_elements )
                continue;
        }

        matl = mtls[i];

        /* Skip if this material is invisible. */
        if ( hide_mtl[matl] )
            continue;

        /* If displaying a result, enable color to change
         * the DIFFUSE property of the material.
         */
        colorflag = show_result && !disable_mtl[matl];
        if ( analy->material_greyscale && disable_mtl[matl] )
            showgs = TRUE;
        else
            showgs = FALSE;

        if ( analy->interp_mode == GOOD_INTERP )
        {
            if ( analy->displace_exag )
            {
                for ( j = 0; j < 3; j++ )
                {
                    nd = connects[i][j];
                    for ( k = 0; k < 2; k++ )
                    {
                        orig = onodes[nd][k];
                        verts[j][k] = orig + analy->displace_scale[k]
                                      * (nodes[nd][k] - orig);
                    }

                    res[j] = data_array[nd];
                }
            }
            else
            {
                for ( j = 0; j < 3; j++ )
                {
                    nd = connects[i][j];
                    for ( k = 0; k < 2; k++ )
                        verts[j][k] = nodes[nd][k];

                    res[j] = data_array[nd];
                }
            }
        }
        else
        {
            if ( analy->displace_exag )
            {
                for ( j = 0; j < 3; j++ )
                {
                    nd = connects[i][j];
                    for ( k = 0; k < 2; k++ )
                    {
                        orig = onodes[nd][k];
                        verts[j][k] = orig + analy->displace_scale[k]
                                      * (nodes[nd][k] - orig);
                    }

                    color_lookup( cols[j], data_array[*p_index_source],
                                  rmin, rmax, analy->zero_result, matl,
                                  analy->logscale, showgs );
                }
            }
            else
            {
                for ( j = 0; j < 3; j++ )
                {
                    nd = connects[i][j];
                    for ( k = 0; k < 2; k++ )
                        verts[j][k] = nodes[nd][k];

                    color_lookup( cols[j], data_array[*p_index_source],
                                  rmin, rmax, analy->zero_result, matl,
                                  analy->logscale, showgs );
                }
            }
        }

        draw_poly_2d( 3, verts, cols, res, matl, p_mesh, analy );
    }

    end_draw_poly( analy );

    /*
     * Draw the outlines of the tri elements.
     */
    if ( analy->mesh_view_mode == RENDER_HIDDEN || analy->mesh_view_mode == RENDER_WIREFRAME ||
            analy->mesh_view_mode != RENDER_WIREFRAME )
    {
        for ( i = 2; i < 9; i += 3 )
            pts[i] = 0.0;

        glColor3fv( v_win->mesh_color );

        for ( i = 0; i < tri_qty; i++ )
        {
            /* Check for inactive elements. */
            if ( activity )
            {
                if ( activity[i] != 0.0 && analy->show_only_deleted_elements )
                    continue;
                if ( activity[i] == 0.0 && !analy->show_deleted_elements && !analy->show_only_deleted_elements )
                    continue;
            }

            /* Skip if this material is invisible. */
            if ( hide_mtl[mtls[i]] )
                continue;

            if ( analy->displace_exag )
            {
                for ( j = 0; j < 3; j++ )
                {
                    nd = connects[i][j];
                    for ( k = 0; k < 2; k++ )
                    {
                        orig = onodes[nd][k];
                        pts[j*3+k] = orig + analy->displace_scale[k]*
                                     (nodes[nd][k] - orig);
                    }
                }
            }
            else
            {
                for ( j = 0; j < 3; j++ )
                {
                    nd = connects[i][j];
                    for ( k = 0; k < 2; k++ )
                        pts[j*3+k] = nodes[nd][k];
                }
            }

            draw_line( 3, pts, mtls[i], p_mesh, analy, TRUE, NULL );
        }
    }

    analy->interp_mode = save_interp_mode;
}


/************************************************************
 * TAG( draw_beams_2d )
 *
 * Draw the beam elements in the model.
 */
static void
draw_beams_2d( Bool_type show_node_result, Bool_type show_mat_result, Bool_type show_mesh_result,
               MO_class_data *p_beam_class, Analysis *analy )
{
    int i, j, k, nd, matl, mesh_idx;
    int (*connects)[3];
    int *p_index_source;
    Mesh_data *p_mesh;
    GVec2D *nodes, *onodes;
    Bool_type show_result, showgs=FALSE;
    float *data_array;
    float *activity;
    float col[4];
    float pts[6];
    float orig, threshold, rmin, rmax;
    unsigned char *hide_mtl, *disable_mtl;
    int *mtls;
    int beam_qty;

    if ( analy->mesh_view_mode != RENDER_FILLED          && analy->mesh_view_mode != RENDER_WIREFRAME &&
            analy->mesh_view_mode != RENDER_WIREFRAMETRANS  && analy->mesh_view_mode != RENDER_HIDDEN )
        return;

    p_mesh = analy->mesh_table + p_beam_class->mesh_id;
    connects = (int (*)[3]) p_beam_class->objects.elems->nodes;
    hide_mtl = p_mesh->hide_material;
    disable_mtl = p_mesh->disable_material;
    mtls = p_beam_class->objects.elems->mat;
    beam_qty = p_beam_class->qty;
    nodes = analy->state_p->nodes.nodes2d;
    onodes = (GVec2D *) analy->cur_ref_state_data;
    activity = analy->state_p->sand_present
               ? analy->state_p->elem_class_sand[p_beam_class->elem_class_index]
               : NULL;

    /* Abstract the source array for data values and its index. */
    if ( show_mat_result )
    {
        data_array = ((MO_class_data **)
                      p_mesh->classes_by_sclass[G_MAT].list)[0]->data_buffer;
        p_index_source = &matl;
    }
    else if ( show_mesh_result )
    {
        mesh_idx = 0;
        data_array = ((MO_class_data **)
                      p_mesh->classes_by_sclass[G_MESH].list)[0]->data_buffer;
        p_index_source = &mesh_idx;
    }
    else
    {
        data_array = p_beam_class->data_buffer;
        p_index_source = &i;
    }

    show_result = result_has_class( analy->cur_result, p_beam_class, analy )
                  || show_node_result
                  || show_mat_result
                  || show_mesh_result;

    get_min_max( analy, FALSE, &rmin, &rmax );
    threshold = analy->zero_result;

    antialias_lines( TRUE, analy->z_buffer_lines );
    glLineWidth( 3.0 );

    /* Set third dimension to zero. */
    pts[2] = pts[5] = 0.0;

    for ( i = 0; i < beam_qty; i++ )
    {
        /* Check for inactive elements. */
        if ( activity )
        {
            if ( activity[i] != 0.0 && analy->show_only_deleted_elements )
                continue;
            if ( activity[i] == 0.0 && !analy->show_deleted_elements && !analy->show_only_deleted_elements )
                continue;
        }

        matl = mtls[i];

        /* Skip if this material is invisible. */
        if ( hide_mtl[matl] || hide_by_object_type( p_beam_class, matl, i, analy, data_array ))
            continue;

        colorflag = show_result && !disable_mtl[matl] &&
                    !disable_by_object_type( p_beam_class, matl, i, analy, NULL );
        if ( analy->material_greyscale && disable_mtl[matl] )
            showgs = TRUE;
        else
            showgs = FALSE;

        color_lookup( col, data_array[*p_index_source],
                      rmin, rmax, threshold,
                      matl, analy->logscale, showgs );

        glColor3fv( col );

        for ( j = 0; j < 2; j++ )
        {
            nd = connects[i][j];

            if ( analy->displace_exag )
            {
                /* Scale the node displacements. */
                for ( k = 0; k < 2; k++ )
                {
                    orig = onodes[nd][k];
                    pts[j*3+k] = orig + analy->displace_scale[k]
                                 * (nodes[nd][k] - orig);
                }
            }
            else
            {
                for ( k = 0; k < 2; k++ )
                    pts[j*3+k] = nodes[nd][k];
            }
        }
        draw_line( 2, pts, matl, p_mesh, analy, FALSE, NULL );
    }

    antialias_lines( FALSE, 0 );
    glLineWidth( 1.0 );
}


/************************************************************
 * TAG( draw_truss_2d )
 *
 * Draw the beam elements in the model.
 */
static void
draw_truss_2d( Bool_type show_mat_result, Bool_type show_mesh_result,
               MO_class_data *p_truss_class, Analysis *analy )
{
    int i, j, k, nd, matl, mesh_idx;
    int (*connects)[3];
    int *p_index_source;
    Mesh_data *p_mesh;
    GVec2D *nodes, *onodes;
    Bool_type show_result, showgs=FALSE;
    float *data_array;
    float *activity;
    float col[4];
    float pts[6];
    float orig, threshold, rmin, rmax;
    unsigned char *hide_mtl, *disable_mtl;
    int *mtls;
    int truss_qty;

    if ( analy->mesh_view_mode != RENDER_FILLED          && analy->mesh_view_mode != RENDER_WIREFRAME &&
            analy->mesh_view_mode != RENDER_WIREFRAMETRANS  && analy->mesh_view_mode != RENDER_HIDDEN )
        return;

    p_mesh = analy->mesh_table + p_truss_class->mesh_id;
    connects = (int (*)[3]) p_truss_class->objects.elems->nodes;
    hide_mtl = p_mesh->hide_material;
    disable_mtl = p_mesh->disable_material;
    mtls = p_truss_class->objects.elems->mat;
    truss_qty = p_truss_class->qty;
    nodes = analy->state_p->nodes.nodes2d;
    onodes = (GVec2D *) analy->cur_ref_state_data;
    activity = analy->state_p->sand_present
               ? analy->state_p->elem_class_sand[p_truss_class->elem_class_index]
               : NULL;

    /* Abstract the source array for data values and its index. */
    if ( show_mat_result )
    {
        data_array = ((MO_class_data **)
                      p_mesh->classes_by_sclass[G_MAT].list)[0]->data_buffer;
        p_index_source = &matl;
    }
    else if ( show_mesh_result )
    {
        mesh_idx = 0;
        data_array = ((MO_class_data **)
                      p_mesh->classes_by_sclass[G_MESH].list)[0]->data_buffer;
        p_index_source = &mesh_idx;
    }
    else
    {
        data_array = p_truss_class->data_buffer;
        p_index_source = &i;
    }

    show_result = result_has_class( analy->cur_result, p_truss_class, analy )
                  || show_mat_result
                  || show_mesh_result;

    get_min_max( analy, FALSE, &rmin, &rmax );
    threshold = analy->zero_result;

    antialias_lines( TRUE, analy->z_buffer_lines );
    glLineWidth( 3.0 );

    /* Set third dimension to zero. */
    pts[2] = pts[5] = 0.0;

    for ( i = 0; i < truss_qty; i++ )
    {
        /* Check for inactive elements. */
        if ( activity )
        {
            if ( activity[i] != 0.0 && analy->show_only_deleted_elements )
                continue;
            if ( activity[i] == 0.0 && !analy->show_deleted_elements && !analy->show_only_deleted_elements )
                continue;
        }

        matl = mtls[i];

        /* Skip if this material is invisible. */
        if ( hide_mtl[matl] )
            continue;

        colorflag = show_result && !disable_mtl[matl] &&
                    !disable_by_object_type( p_truss_class, matl, i, analy, data_array );
        if ( analy->material_greyscale && disable_mtl[matl] )
            showgs = TRUE;
        else
            showgs = FALSE;

        color_lookup( col,  data_array[*p_index_source],
                      rmin, rmax, threshold,
                      matl, analy->logscale, showgs );

        glColor3fv( col );

        for ( j = 0; j < 2; j++ )
        {
            nd = connects[i][j];

            if ( analy->displace_exag )
            {
                /* Scale the node displacements. */
                for ( k = 0; k < 2; k++ )
                {
                    orig = onodes[nd][k];
                    pts[j*3+k] = orig + analy->displace_scale[k]
                                 * (nodes[nd][k] - orig);
                }
            }
            else
            {
                for ( k = 0; k < 2; k++ )
                    pts[j*3+k] = nodes[nd][k];
            }
        }
        draw_line( 2, pts, matl, p_mesh, analy, FALSE, NULL );
    }

    antialias_lines( FALSE, 0 );
    glLineWidth( 1.0 );
}


/************************************************************
 * TAG( draw_nodes_2d_3d )
 *
 * Draw the mesh nodes for "point cloud" rendering.
 */
static void
draw_nodes_2d_3d( MO_class_data *p_node_class, Analysis *analy )
{
    GVec3D *coords3;
    GVec2D *coords2;
    float pt[3];
    float col[4];
    float *nodal_data;
    int i, j;
    int node_qty;
    Bool_type no_result;

    node_qty = p_node_class->qty;

    /*
     * If "p_node_class" isn't the same class accessed by the
     * NODAL_RESULT_BUFFER macro (i.e., it isn't the Mesh_data "p_node_geom"
     * class), we have some work to do.  If true and data was elemental,
     * it was interpolated into NODAL_RESULT_BUFFER and doesn't exist for
     * p_node_class nodes.  In that case it's as if analy->cur_result were
     * NULL.  If data was nodal, it will only exist in
     * p_node_class->data_buffer, but that's OK for this call.
     */
    if ( analy->cur_result == NULL
            || ( p_node_class->data_buffer != NODAL_RESULT_BUFFER( analy )
                 && analy->cur_result->origin.is_elem_result ) )
        no_result = TRUE;
    else
        no_result = FALSE;

    nodal_data = p_node_class->data_buffer;

    colorflag = TRUE;

    if ( analy->dimension == 3 )
    {
        /* glEnable( GL_POINT_SMOOTH );
           glEnable( GL_BLEND );
           glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
        */
        coords3 = analy->state_p->nodes.nodes3d;
        if ( analy->point_diam>0.0 )
        {
            glEnable( GL_POINT_SMOOTH );
            glPointSize( analy->point_diam );
        }

        glBegin( GL_POINTS );

        for ( i = 0; i < node_qty; i++ )
        {
            for ( j = 0; j < 3; j++ )
                pt[j] = coords3[i][j];

            if ( no_result )
            {
                /* No result, color by background color */
                VEC_COPY( col, v_win->foregrnd_color );
            }
            else
                color_lookup( col, nodal_data[i], analy->result_mm[0],
                              analy->result_mm[1], analy->zero_result, -1,
                              analy->logscale, analy->material_greyscale );

            glColor3fv( col );
            glVertex3fv( pt );
        }

        glEnd();
    }
    else
    {
        coords2 = analy->state_p->nodes.nodes2d;

        glBegin( GL_POINTS );

        for ( i = 0; i < node_qty; i++ )
        {
            for ( j = 0; j < 2; j++ )
                pt[j] = coords2[i][j];

            if ( no_result )
            {
                VEC_COPY( col, v_win->foregrnd_color );
            }
            else
                color_lookup( col, nodal_data[i], analy->result_mm[0],
                              analy->result_mm[1], analy->zero_result, -1,
                              analy->logscale, analy->material_greyscale );


            glColor3fv( col );
            glVertex2fv( pt );
        }

        glEnd();
    }
}


/************************************************************
 * TAG( draw_surfaces_2d )
 *
 * Draw the surface elements for a class in the model.
 */
static void
draw_surfaces_2d( Color_property *p_cp,
                  Bool_type show_node_result, Bool_type show_surf_result,
                  Bool_type show_mesh_result, MO_class_data *p_surface_class,
                  Analysis *analy )
{
    Bool_type show_result, showgs=FALSE;
    float verts[4][3];
    float cols[4][4];
    float pts[12];
    float res[4];
    float rmin, rmax;
    int surf, matl, cnt, nd, i, j, k, mesh_idx;
    int (*connects)[4];
    unsigned char *hide_surf, *disable_surf;
    Mesh_data *p_mesh;
    Surface_data *p_surface;
    Visibility_data *p_vd;
    int facet_qty;
    int surface_qty;
    float *data_array;
    int *p_index_source;
    GVec2D *nodes, *onodes;
    Interp_mode_type save_interp_mode;

    if ( analy->mesh_view_mode != RENDER_FILLED &&
            analy->mesh_view_mode != RENDER_WIREFRAME &&
            analy->mesh_view_mode != RENDER_WIREFRAMETRANS &&
            analy->mesh_view_mode != RENDER_HIDDEN )
        return;

    p_mesh = MESH_P( analy );
    p_vd = p_surface_class->p_vis_data;
    p_surface = p_surface_class->objects.surfaces;
    surface_qty = p_surface_class->qty;
    disable_surf = p_mesh->disable_surface;
    hide_surf = p_mesh->hide_surface;

    nodes = analy->state_p->nodes.nodes2d;
    onodes = (GVec3D *) analy->cur_ref_state_data;
    /* 2D, so set Z coord of every vertex to zero. */
    verts[0][2] = verts[1][2] = verts[2][2] = verts[3][2] = 0.0;

    /* Override interpolation mode for material and mesh results.  */
    save_interp_mode = analy->interp_mode;
    if ( show_surf_result || show_mesh_result )
        analy->interp_mode = NO_INTERP;

    /* Abstract the source array for data values and its index. */
    if ( analy->interp_mode == GOOD_INTERP )
    {
        data_array = NODAL_RESULT_BUFFER( analy );
        p_index_source = &nd;
    }
    else if ( show_surf_result )
    {
        data_array = ((MO_class_data **)
                      p_mesh->classes_by_sclass[G_SURFACE].list)[0]->data_buffer;
        p_index_source = &matl;
    }
    else if ( show_mesh_result )
    {
        mesh_idx = 0;
        data_array = ((MO_class_data **)
                      p_mesh->classes_by_sclass[G_MESH].list)[0]->data_buffer;
        p_index_source = &mesh_idx;
    }
    else if ( analy->interp_mode == NO_INTERP && !show_node_result )
    {
        data_array = p_surface_class->data_buffer;
        p_index_source = &i;
    }
    else
    {
        data_array = NODAL_RESULT_BUFFER( analy );
        p_index_source = &nd;
    }

    show_result = show_node_result
                  || result_has_class( analy->cur_result, p_surface_class, analy )
                  || show_mesh_result;

    get_min_max( analy, FALSE, &rmin, &rmax );

    /* Set up for polygon drawing. */
    begin_draw_poly( analy );

    for ( surf = 0; surf < surface_qty; surf++ )
    {
        /* Skip if this surface is invisible. */
        if ( hide_surf[surf] )
            continue;

        connects = (int (*)[4]) p_surface[surf].nodes;
        p_vd = p_surface[surf].p_vis_data;
        if( v_win->surfaces.current_index != p_surface[surf].surface_id )
            change_current_color_property( &v_win->surfaces,
                                           p_surface[surf].surface_id );

        /* If displaying a result, enable color to change
         * the DIFFUSE property of the material.
         */
        colorflag = show_result && !disable_surf[surf];
        if ( analy->material_greyscale && disable_surf[surf] )
            showgs = TRUE;
        else
            showgs = FALSE;

        facet_qty = p_surface[surf].facet_qty;
        for ( i = 0; i < facet_qty; i++ )
        {
            get_surface_verts_2d( i, p_surface[surf].nodes,
                                  p_surface[surf].surface_id,
                                  analy, verts );

            if ( analy->interp_mode == GOOD_INTERP )
            {
                for ( j = 0; j < 4; j++ )
                {
                    nd = connects[i][j];
                    for ( k = 0; k < 2; k++ )
                        verts[j][k] = nodes[nd][k];

                    res[j] = data_array[nd];
                }
            }
            else
            {
                for ( j = 0; j < 4; j++ )
                {
                    nd = connects[i][j];
                    for ( k = 0; k < 2; k++ )
                        verts[j][k] = nodes[nd][k];

                    surf_color_lookup( cols[j], data_array[*p_index_source],
                                       rmin, rmax, analy->zero_result, matl,
                                       analy->logscale );
                }
            }
        }

        draw_poly_2d( 4, verts, cols, res, matl, p_mesh, analy );
    }

    end_draw_poly( analy );
    analy->interp_mode = save_interp_mode;
}


/************************************************************
 * TAG( draw_surfaces_3d )
 *
 * Draw the surface elements in the model.
 */
static void
draw_surfaces_3d( Color_property *p_cp,
                  Bool_type show_node_result, Bool_type show_surf_result,
                  Bool_type show_mesh_result, MO_class_data *p_surface_class,
                  Analysis *analy )
{
    Bool_type show_result, showgs=FALSE;
    Bool_type has_degen;
    float verts[4][3];
    float norms[4][3];
    float cols[4][4];
    float res[4];
    float rmin, rmax;
    int surf, matl, cnt, nd, i, j, k, mesh_idx;
    int (*connects)[4];
    unsigned char *disable_surf, *hide_surf;
    Mesh_data *p_mesh;
    Surface_data *p_surface;
    Visibility_data *p_vd;
    int facet_qty;
    int surface_qty;
    float *data_array;
    int *p_index_source;
    Interp_mode_type save_interp_mode;

    Bool_type hidden_poly;

    if ( analy->mesh_view_mode != RENDER_FILLED &&
            analy->mesh_view_mode != RENDER_WIREFRAME &&
            analy->mesh_view_mode != RENDER_WIREFRAMETRANS &&
            analy->mesh_view_mode != RENDER_HIDDEN )
        return;


    p_mesh = MESH_P( analy );
    p_vd = p_surface_class->p_vis_data;
    p_surface = p_surface_class->objects.surfaces;
    surface_qty = p_surface_class->qty;
    disable_surf = p_mesh->disable_surface;
    hide_surf = p_mesh->hide_surface;

    /*
     * This has not been implemented with respect to surfaces.
     *
    has_degen = p_surface->objects.elems->has_degen;
     */

    /*
     * Is a result constant over the surface or unique per facet?
     *
     * To be resolved by D. Speck
     */

    /* Override interpolation mode for mesh results. */
    save_interp_mode = analy->interp_mode;
    if ( show_mesh_result )
        analy->interp_mode = NO_INTERP;

    /* Abstract the source array for data values and its index. */
    if ( analy->interp_mode == GOOD_INTERP )
    {
        data_array = NODAL_RESULT_BUFFER( analy );
        p_index_source = &nd;
    }
    else if ( show_mesh_result )
    {
        mesh_idx = 0;
        data_array = ((MO_class_data **)
                      p_mesh->classes_by_sclass[G_MESH].list)[0]->data_buffer;
        p_index_source = &mesh_idx;
    }
    else if ( analy->interp_mode == NO_INTERP && show_surf_result )
    {
        data_array = ((MO_class_data **)
                      p_mesh->classes_by_sclass[G_SURFACE].list)[0]->data_buffer;
        p_index_source = &surf;
    }
    else if ( analy->interp_mode == NO_INTERP )
    {
        data_array = p_surface_class->data_buffer;
        p_index_source = &i;
    }
    else
    {
        data_array = NODAL_RESULT_BUFFER( analy );
        p_index_source = &nd;
    }

    show_result = show_node_result
                  || result_has_class( analy->cur_result, p_surface_class, analy )
                  || show_mesh_result;

    get_min_max( analy, FALSE, &rmin, &rmax );

    /* Use two-sided lighting for surfaces. */
    glLightModeli( GL_LIGHT_MODEL_TWO_SIDE, 1 );

    /* Set up for polygon drawing. */
    begin_draw_poly( analy );
    for( surf = 0; surf < surface_qty; surf++ )
    {
        /* Skip if surface is invisible. */
        if( hide_surf[surf] )
            continue;

        connects = (int (*)[4]) p_surface[surf].nodes;
        p_vd = p_surface[surf].p_vis_data;
        if( v_win->surfaces.current_index != p_surface[surf].surface_id )
            change_current_color_property( &v_win->surfaces,
                                           p_surface[surf].surface_id );

        /*
         * Colorflag is TRUE if displaying a result, otherwise
         * polygons are drawn in the surface color.
         */
        colorflag = show_result && !disable_surf[surf];
        if ( analy->material_greyscale && disable_surf[surf] )
            showgs = TRUE;
        else
            showgs = FALSE;

        facet_qty = p_surface[surf].facet_qty;
        for ( i = 0; i < facet_qty; i++ )
        {
            get_surface_verts_3d( i, p_surface[surf].nodes,
                                  p_surface[surf].surface_id,
                                  analy, verts );

            if ( analy->interp_mode == GOOD_INTERP )
            {
                for ( j = 0; j < 4; j++ )
                {
                    nd = connects[i][j];

                    for ( k = 0; k < 3; k++ )
                        norms[j][k] = p_vd->face_norm[j][k][i];

                    res[j] = data_array[nd];
                }
            }
            else
            {
                for ( j = 0; j < 4; j++ )
                {
                    nd = connects[i][j];

                    for ( k = 0; k < 3; k++ )
                        norms[j][k] = p_vd->face_norm[j][k][i];

                    surf_color_lookup( cols[j], data_array[*p_index_source],
                                       rmin, rmax, analy->zero_result,
                                       p_surface[surf].surface_id, analy->logscale );
                }
            }

            /* Check for triangular (degenerate) surface element. */
            cnt = 4;

            /*
             * See above comment regarding use of "has_degen"

            if ( has_degen  &&
                 connects[i][2] == connects[i][3] )
                 cnt = 3;
             */

            hidden_poly = FALSE;
            draw_poly( cnt, verts, norms, cols, res, p_surface[surf].surface_id,
                       p_mesh, analy, hidden_poly );
        }
    }

    /* Back to the defaults. */
    end_draw_poly( analy );
    glLightModeli( GL_LIGHT_MODEL_TWO_SIDE, 0 );
    analy->interp_mode = save_interp_mode;
}


/************************************************************
 * TAG( particle_error )
 *
 * GL utility library error callback for sphere rendering below.
 */
static void
particle_error( GLenum err_code )
{
    const GLubyte *err_string;

    err_string = gluErrorString( err_code );
    fprintf( stderr, "GLU error during particle quadric rendering: %s\n",
             (char *) err_string );
}


/************************************************************
 * TAG( draw_particles_3d )
 *
 * Render particles from a particle unit class as spheres.
 */
static void
draw_particles_3d( MO_class_data *p_particle_class, Analysis *analy )
{
    GVec3D *coords3;
    float *el_state_mm, *part_res;
    float col[4];
    int i;
    int part_qty;
    GLUquadricObj *sphere;
    GLuint display_list;
    int index, node, matl;
    float x, y, z, rmin, rmax;
    Bool_type show_result = FALSE;

    MO_class_data * p_node_class = NULL;

    p_node_class = MESH(analy).node_geom;

    part_qty = p_particle_class->qty;
    part_res = p_particle_class->data_buffer;
    el_state_mm = analy->elem_state_mm.object_minmax;

    index = 0;

    colorflag = analy->cur_result != NULL;

    /*coords3 = analy->state_p->particles.particles3d;*/
    coords3 = analy->state_p->nodes.particles3d;

    if(analy->cur_result != NULL)
    {
        show_result = result_has_class(analy->cur_result, p_particle_class, analy);
    }
  
    get_min_max(analy, FALSE, &rmin, &rmax);

    display_list = glGenLists( 1 );

    sphere = gluNewQuadric();

    if ( sphere == NULL )
    {
        popup_dialog( WARNING_POPUP,
                      "Unable to create sphere quadric for particle rendering."
                    );
        return;
    }

    gluQuadricCallback( sphere, (GLenum) GLU_ERROR, particle_error );

    gluQuadricDrawStyle( sphere, (GLenum) GLU_FILL );
    gluQuadricNormals( sphere, (GLenum) GLU_SMOOTH );
    gluQuadricTexture(sphere, (GLboolean) GLU_TRUE );
    glNewList( display_list, GL_COMPILE );
    /*gluSphere( sphere, particle_radius, 16, 8 ); */
    gluSphere( sphere, particle_radius, 32, 16 );
    glEndList();

    if ( v_win->lighting )
        glEnable( GL_LIGHTING );

    glEnable( GL_COLOR_MATERIAL );

    for ( i = 0; i < part_qty; i++ )
    {
        if(p_particle_class->p_vis_data->visib[i] == FALSE)
        {
            continue;
        }
 
        matl = p_particle_class->objects.elems->mat[i];

        if(v_win->mesh_materials.current_index != matl)
        {
            change_current_color_property(&v_win->mesh_materials, matl);
        }

        /*color_lookup( col, part_res[i], el_state_mm[0], el_state_mm[1],
                      analy->zero_result, index, analy->logscale,
                      analy->material_greyscale ); */

        color_lookup( col, part_res[i], rmin, rmax,
                      analy->zero_result, matl, analy->logscale,
                      analy->material_greyscale );

        glColor3fv( col );

        glPushMatrix();
        node = p_particle_class->objects.elems->nodes[i];
        x = p_node_class->objects.nodes3d[node][0]; 
        y = p_node_class->objects.nodes3d[node][1]; 
        z = p_node_class->objects.nodes3d[node][2];
  
        /*glTranslatef( coords3[i][0], coords3[i][1], coords3[i][2] ); */
        glTranslatef( x, y, z );

        glCallList( display_list );

        glPopMatrix();
    }

    glDisable( GL_COLOR_MATERIAL );

    if ( v_win->lighting )
        glDisable( GL_LIGHTING );

    gluDeleteQuadric( sphere );

    glDeleteLists( display_list, 1 );
}


/************************************************************
 * TAG( draw_particles_2d )
 *
 * Render particles from a particle unit class as circles.
 */
static void
draw_particles_2d( MO_class_data *p_particle_class, Analysis *analy )
{
    GVec2D *coords2;
    float *el_state_mm, *part_res;
    float col[4];
    int i;
    int part_qty;
    GLUquadricObj *circle;
    GLuint display_list;
    int mat_no;

    part_qty = p_particle_class->qty;
    part_res = p_particle_class->data_buffer;
    el_state_mm = analy->elem_state_mm.object_minmax;


    mat_no = (MESH_P( analy )->material_qty) % MATERIAL_COLOR_CNT;

    colorflag = analy->cur_result != NULL;

    coords2 = analy->state_p->particles.particles2d;

    display_list = glGenLists( 1 );

    circle = gluNewQuadric();

    if ( circle == NULL )
    {
        popup_dialog( WARNING_POPUP,
                      "Unable to create disk quadric for particle rendering." );
        return;
    }

    gluQuadricCallback( circle, (GLenum) GLU_ERROR, particle_error );

    gluQuadricDrawStyle( circle, (GLenum) GLU_FILL );
    gluQuadricNormals( circle, (GLenum) GLU_SMOOTH );

    glNewList( display_list, GL_COMPILE );
    gluDisk( circle, (GLdouble) 0.0, particle_radius, 16, 5 );
    glEndList();

    if ( v_win->lighting )
        glEnable( GL_LIGHTING );

    glEnable( GL_COLOR_MATERIAL );

    for ( i = 0; i < part_qty; i++ )
    {
        color_lookup( col, part_res[i], el_state_mm[0], el_state_mm[1],
                      analy->zero_result, mat_no, analy->logscale,
                      analy->material_greyscale );
        glColor3fv( col );

        glPushMatrix();
        glTranslatef( coords2[i][0], coords2[i][1], 0.0 );

        glCallList( display_list );

        glPopMatrix();
    }

    glDisable( GL_COLOR_MATERIAL );

    if ( v_win->lighting )
        glDisable( GL_LIGHTING );

    gluDeleteQuadric( circle );

    glDeleteLists( display_list, 1 );
}


/************************************************************
 * TAG( draw_edges_3d )
 *
 * Draw the mesh edges.
 */
static void
draw_edges_3d( Analysis *analy )
{
    GVec3D *nodes, *onodes;
    float pts[6];
    float orig;
    int nd, i, j, k, index;
    Mesh_data *p_mesh;
    int edge_qty;
    Edge_obj *eo_array;
    int enodes[2];

    /* Variables used for removing edges on reflection
     * planes.
     */
    float reflect_pts[10][6], num_refl_sets = 0, point_diff;
    Refl_plane_obj *plane;

    Bool_type draw_refl_edge[10] = {TRUE,TRUE,TRUE,TRUE,TRUE,
                                    TRUE,TRUE,TRUE,TRUE,TRUE
                                   };

    Bool_type point_diff_flag, draw_refl_flag;

    int  c_index, c1_index, c2_index, num_planes = 0, plane_index;

    p_mesh = MESH_P( analy );

    /* Sanity check. */
    if ( p_mesh->edge_list == NULL )
        return;

    eo_array = p_mesh->edge_list->list;
    edge_qty = p_mesh->edge_list->size;
    nodes = analy->state_p->nodes.nodes3d;
    onodes = (GVec3D *) analy->cur_ref_state_data;

    antialias_lines( TRUE, analy->z_buffer_lines );

    /* Bias the depth buffer so edges are in front of polygons. */
    if ( !env.win32 )
        glDepthRange( 0, 1 - analy->edge_zbias );

    glLineWidth( analy->edge_width );
    glColor3fv( v_win->edge_color  );

    if ( analy->reflect && analy->refl_planes )
        for ( plane = analy->refl_planes;
                plane != NULL;
                plane = plane->next )
            num_planes++;

    for ( i = 0; i < edge_qty; i++ )
    {
        enodes[0] = eo_array[i].node1;
        enodes[1] = eo_array[i].node2;

        for ( j = 0; j < 2; j++ )
        {
            nd = enodes[j];

            if ( analy->displace_exag )
            {
                /* Scale the node displacements. */
                for ( k = 0; k < 3; k++ )
                {
                    orig = onodes[nd][k];
                    pts[j*3+k] = orig + analy->displace_scale[k]*
                                 (nodes[nd][k] - orig);
                }
            }
            else
            {
                for ( k = 0; k < 3; k++ )
                    pts[j*3+k] = nodes[nd][k];
            }
        }

        /*
         * Draw edges for reflections.
         */

        draw_refl_flag = TRUE;

        if ( analy->reflect && analy->refl_planes )
        {
            if ( analy->refl_orig_only )
            {
                plane_index = 0;
                for ( plane = analy->refl_planes;
                        plane != NULL;
                        plane = plane->next )
                {
                    point_transform( &reflect_pts[0][0], &pts[0], &plane->pt_transf );
                    point_transform( &reflect_pts[0][3], &pts[3], &plane->pt_transf );

                    for (index=0;
                            index<6;
                            index++)
                    {
                        point_diff = fabs(reflect_pts[index][0]-pts[index]);
                        if ( point_diff > 1.0e-5 )
                            point_diff_flag = TRUE;
                    }
                    if (!point_diff_flag)
                        draw_refl_edge[plane_index] = FALSE;
                }
            }
            else
            {
                /* Cumulative reflections */

                plane_index = 0;
                plane       = analy->refl_planes;
                point_transform( &reflect_pts[0][0], &pts[0], &plane->pt_transf );
                point_transform( &reflect_pts[0][3], &pts[3], &plane->pt_transf );

                point_diff_flag = FALSE;
                draw_refl_flag  = TRUE;

                for (index=0;
                        index<6;
                        index++)
                {
                    point_diff = fabs(reflect_pts[0][index]-pts[index]);
                    if ( point_diff > 1.0e-5 )
                        point_diff_flag = TRUE;
                }

                if (!point_diff_flag)
                {
                    draw_refl_edge[plane_index] = FALSE;
                    draw_refl_flag = FALSE;
                }
                else
                {
                    draw_refl_edge[plane_index] = TRUE;
                }

                plane_index = 1;

                for ( plane = analy->refl_planes->next;
                        plane != NULL;
                        plane = plane->next )
                {
                    c1_index = num_planes*plane_index;
                    c2_index = 0;
                    for ( c_index=0;
                            c_index<plane_index;
                            c_index++ )
                    {
                        point_transform( &reflect_pts[c1_index][0], &reflect_pts[c2_index][0], &plane->pt_transf );
                        point_transform( &reflect_pts[c1_index][3], &reflect_pts[c2_index][3], &plane->pt_transf );

                        c1_index++;
                        c2_index++;
                    }

                    for (index=0;
                            index<6;
                            index++)
                    {
                        point_diff = fabs(reflect_pts[c1_index][index]-pts[index]);
                        if ( point_diff > 1.0e-5 )
                            point_diff_flag = TRUE;
                    }

                    if (!point_diff_flag)
                    {
                        draw_refl_edge[plane_index] = FALSE;
                        if ( num_planes>1)
                            draw_refl_flag = FALSE;
                    }
                    else
                    {
                        draw_refl_edge[plane_index] = TRUE;
                    }

                    point_diff_flag = FALSE;
                    plane_index++;
                }
            }
        }

        if ( draw_refl_flag )
            draw_line( 2, pts, eo_array[i].mtl, p_mesh, analy, FALSE, draw_refl_edge );
    }

    /* Remove depth bias. */
    glDepthRange( 0, 1 );

    glLineWidth( 1.0 );

    antialias_lines( FALSE, 0 );
}


/************************************************************
 * TAG( draw_edges_2d )
 *
 * Draw the mesh edges.
 */
static void
draw_edges_2d( Analysis *analy )
{
    GVec2D *nodes, *onodes;
    float pts[6];
    float orig;
    int nd, i, j, k;
    Mesh_data *p_mesh;
    int edge_qty;
    Edge_obj *eo_array;
    int enodes[2];

    p_mesh = MESH_P( analy );

    /* Sanity check. */
    if ( p_mesh->edge_list == NULL )
        return;

    eo_array = p_mesh->edge_list->list;
    edge_qty = p_mesh->edge_list->size;
    nodes = analy->state_p->nodes.nodes2d;
    onodes = (GVec2D *) analy->cur_ref_state_data;

    /* Set Z-coord of endpoints. */
    pts[2] = pts[5] = 0.0;

    antialias_lines( TRUE, analy->z_buffer_lines );

    /* Bias the depth buffer so edges are in front of polygons. */
    glDepthRange( 0, 1 - analy->edge_zbias );

    glLineWidth( analy->edge_width );

    glColor3fv( v_win->edge_color  );
    for ( i = 0; i < edge_qty; i++ )
    {
        enodes[0] = eo_array[i].node1;
        enodes[1] = eo_array[i].node2;

        for ( j = 0; j < 2; j++ )
        {
            nd = enodes[j];

            if ( analy->displace_exag )
            {
                /* Scale the node displacements. */
                for ( k = 0; k < 2; k++ )
                {
                    orig = onodes[nd][k];
                    pts[j*3+k] = orig + analy->displace_scale[k]*
                                 (nodes[nd][k] - orig);
                }
            }
            else
            {
                for ( k = 0; k < 2; k++ )
                    pts[j*3+k] = nodes[nd][k];
            }
        }

        draw_line( 2, pts, eo_array[i].mtl, p_mesh, analy, FALSE, NULL );
    }

    /* Remove depth bias. */
    glDepthRange( 0, 1 );

    glLineWidth( 1.0 );

    antialias_lines( FALSE, 0 );
}


/************************************************************
 * TAG( draw_hilite )
 *
 * If a node or element is highlighted, this routine draws
 * the highlight.
 */
static void
draw_hilite( Bool_type hilite, MO_class_data *p_mo_class, int hilite_num,
             Analysis *analy )
{
    float verts[8][3], leng[3], vec[3];
    float node_base_radius, rfac, radius;
    char label[64];
    int vert_cnt;
    int i, j;
    Bool_type verts_ok, show_label;
    float val;
    float *hilite_col;
    int fracsz;
    char *cname;
    int dim, obj_qty;
    Bool_type particle_hilite = FALSE;
    float *data_array;
    Surface_data *p_surface;
    int facet_qty, facet;

    Bool_type pn_hilite=FALSE, result_defined=FALSE, valid_free_node=FALSE;
    float particle_select_col[3] = {.4, .4, .4};
    int hilite_label;
    int p_near;
    Mesh_data *p_mesh;

    rfac = 0.1;                             /* Radius scale factor. */
    fracsz = analy->float_frac_size;
    cname = p_mo_class->long_name;
    dim = analy->dimension;
    obj_qty = p_mo_class->qty;
    data_array = p_mo_class->data_buffer;

    if ( p_mo_class->labels_found )
        hilite_label = get_class_label( p_mo_class, hilite_num );
    else
        hilite_label = hilite_num;

    if ( is_particle_class( analy, p_mo_class->superclass, p_mo_class->short_name ) &&
            analy->particle_nodes_enabled &&
            p_mo_class->labels_found )
    {
        p_mesh = MESH_P( analy );
        select_meshless_node( analy, p_mesh,
                              p_mo_class, hilite_num,
                              &p_near );
        analy->hilite_ml_node = p_near;
    }

    /**/
    /* Change superclass tested to G_PARTICLE when it exists. */
    if ( p_mo_class->short_name )
    {
        if ( p_mo_class->superclass == G_UNIT
                && strcmp( p_mo_class->short_name, particle_cname ) == 0 )
            particle_hilite = TRUE;
        else
            particle_hilite = FALSE;
    }

    /* Validity check. */
    if( p_mo_class->superclass !=  G_SURFACE)
    {
        if ( !is_particle_class( analy, p_mo_class->superclass, p_mo_class->short_name ) )
            if ( hilite_num < 0 || hilite_num > obj_qty )
            {
                if ( analy->mesh_view_mode != RENDER_POINT_CLOUD )
                    popup_dialog( INFO_POPUP, "%s ident out of range.", cname );
                return;
            }
    }
    else
    {
        p_surface = p_mo_class->objects.surfaces;
        facet_qty = p_surface[0].facet_qty;
        if ( hilite_num < 0 || hilite_num >= facet_qty )
        {
            popup_dialog( INFO_POPUP, "%s ident out of range.", cname );
            return;
        }
    }

    hilite_col = hilite ? v_win->hilite_color : v_win->select_color;

    show_label = hilite
                 || p_mo_class->superclass != G_QUAD
                 || !analy->loc_ref;

    /* Outline the hilighted mesh object. */
    switch ( p_mo_class->superclass )
    {
    case G_UNIT:
        if ( !particle_hilite )
            /* Don't hilite G_UNIT objects that aren't particles. */
            return;

        /* Highlight a particle (just label it). */
        if ( analy->cur_result != NULL )
        {
            val = analy->perform_unit_conversion
                  ? data_array[hilite_num]
                  * analy->conversion_scale + analy->conversion_offset
                  : data_array[hilite_num];
            sprintf( label, " %s %d (%.*e)", cname, hilite_label, fracsz,
                     val );
        }
        else
            sprintf( label, " %s %d", cname, hilite_label );
        break;

    case G_NODE:
        /* Highlight a node. */

        /* See discussion in draw_nodes_2d_3d(). */
        if ( analy->cur_result != NULL
                && ( !analy->cur_result->origin.is_elem_result
                     || data_array == NODAL_RESULT_BUFFER( analy ) ) )
        {
            if ( analy->particle_nodes_enabled &&
                    analy->free_nodes_list &&
                    analy->free_nodes_vals )
            {
                if ( analy->free_nodes_list[hilite_num]==TRUE )
                    val = analy->perform_unit_conversion
                          ? analy->free_nodes_vals[hilite_num] * analy->conversion_scale
                          + analy->conversion_offset
                          : analy->free_nodes_vals[hilite_num];
                if(data_array == NODAL_RESULT_BUFFER(analy))
                {
                    val = analy->perform_unit_conversion ? data_array[hilite_num] * analy->conversion_scale
                                                  + analy->conversion_offset
                                                  : data_array[hilite_num];
                    sprintf( label, " %s %d (%.*e)", cname, hilite_label, fracsz,
                             val );
                } else
                { 
                    sprintf( label, " %s %d", cname, hilite_label );
                }
            }
            else
            {
                val = analy->perform_unit_conversion
                      ? data_array[hilite_num] * analy->conversion_scale
                      + analy->conversion_offset
                      : data_array[hilite_num];

                sprintf( label, " %s %d (%.*e)", cname, hilite_label, fracsz,
                         val );
            }
        }
        else
            sprintf( label, " %s %d", cname, hilite_label );

        get_node_vert_2d_3d( hilite_num, p_mo_class, analy, verts[0] );
        vert_cnt = 1;

        for ( i = 0; i < dim; i++ )
            leng[i] = analy->bbox[1][i] - analy->bbox[0][i];
        radius = 0.01 * (leng[0] + leng[1] + leng[2]) / 3.0;
        radius *= 1.0 / v_win->scale[0];
        break;

    case G_TRUSS:
        /* Highlight a Truss element. */

        verts_ok = get_truss_verts_2d_3d( hilite_num, p_mo_class, analy,
                                          verts );
        if ( !verts_ok )
        {
            popup_dialog( WARNING_POPUP, "%s %d ignored; bad node(s).\n",
                          cname, hilite_label );
            return;
        }
        vert_cnt = 2;
        if ( dim == 2 )
            verts[0][2] = verts[1][2] = 0.0;

        if ( result_has_superclass( analy->cur_result, G_TRUSS, analy ) )
        {
            val = analy->perform_unit_conversion
                  ? data_array[hilite_num]
                  * analy->conversion_scale + analy->conversion_offset
                  : data_array[hilite_num];
            sprintf( label, " %s %d (%.*e)", cname, hilite_label, fracsz,
                     val );
        }
        else
            sprintf( label, " %s %d", cname, hilite_label );

        VEC_SUB( vec, verts[1], verts[0] );
        leng[0] = VEC_LENGTH( vec );
        radius = rfac * leng[0];

        /* Outline the element. */
        glDepthFunc( GL_ALWAYS );
        glColor3fv( hilite_col );
        glBegin( GL_LINES );
        glVertex3fv( verts[0] );
        glVertex3fv( verts[1] );
        glEnd();
        glDepthFunc( GL_LEQUAL );

        break;

    case G_BEAM:
        /* Highlight a beam element. */

        verts_ok = get_beam_verts_2d_3d( hilite_num, p_mo_class, analy,
                                         verts );
        if ( !verts_ok )
        {
            popup_dialog( WARNING_POPUP, "%s %d ignored; bad node(s).\n",
                          cname, hilite_label );
            return;
        }
        vert_cnt = 2;
        if ( dim == 2 )
            verts[0][2] = verts[1][2] = 0.0;

        if ( result_has_superclass( analy->cur_result, G_BEAM, analy ) )
        {
            val = analy->perform_unit_conversion
                  ? data_array[hilite_num]
                  * analy->conversion_scale + analy->conversion_offset
                  : data_array[hilite_num];
            sprintf( label, " %s %d (%.*e)", cname, hilite_label, fracsz,
                     val );
        }
        else
            sprintf( label, " %s %d", cname, hilite_label );

        VEC_SUB( vec, verts[1], verts[0] );
        leng[0] = VEC_LENGTH( vec );
        radius = rfac * leng[0];

        /* Outline the element. */
        glDepthFunc( GL_ALWAYS );
        glColor3fv( hilite_col );
        glBegin( GL_LINES );
        glVertex3fv( verts[0] );
        glVertex3fv( verts[1] );
        glEnd();
        glDepthFunc( GL_LEQUAL );

        break;

    case G_TRI:
        /* Highlight a tri element. */

        if ( result_has_superclass( analy->cur_result, G_TRI, analy ) )
        {
            val = analy->perform_unit_conversion
                  ? data_array[hilite_num]
                  * analy->conversion_scale + analy->conversion_offset
                  : data_array[hilite_num];
            sprintf( label, " %s %d (%.*e)", cname, hilite_label, fracsz,
                     val );
        }
        else
            sprintf( label, " %s %d", cname, hilite_label );

        if ( dim == 3 )
            get_tri_verts_3d( hilite_num, p_mo_class, analy, verts );
        else
        {
            get_tri_verts_2d( hilite_num, p_mo_class, analy, verts );
            for ( i = 0 ; i < 3; i++ )
                verts[i][2] = 0.0;
        }
        vert_cnt = 3;

        VEC_SUB( vec, verts[1], verts[0] );
        leng[0] = VEC_LENGTH( vec );
        VEC_SUB( vec, verts[2], verts[0] );
        leng[1] = VEC_LENGTH( vec );
        radius = rfac * (leng[0] + leng[1]) / 2.0;

        /* Outline the element. */
        glDepthFunc( GL_ALWAYS );
        glColor3fv( hilite_col );
        glBegin( GL_LINE_LOOP );
        for ( i = 0; i < 3; i++ )
            glVertex3fv( verts[i] );
        glEnd();
        glDepthFunc( GL_LEQUAL );
        break;

    case G_QUAD:
        /* Highlight a quad element. */

        if ( show_label )
        {
            if ( result_has_superclass( analy->cur_result, G_QUAD, analy ) )
            {
                val = analy->perform_unit_conversion
                      ? data_array[hilite_num]
                      * analy->conversion_scale + analy->conversion_offset
                      : data_array[hilite_num];
                sprintf( label, " %s %d (%.*e)", cname, hilite_label,
                         fracsz, val );
            }
            else
                sprintf( label, " %s %d", cname, hilite_label );
        }

        if ( dim == 3 )
            get_quad_verts_3d( hilite_num, p_mo_class->objects.elems->nodes,
                               p_mo_class->mesh_id, analy, verts );
        else
        {
            get_quad_verts_2d( hilite_num, p_mo_class, analy, verts );
            for ( i = 0 ; i < 4; i++ )
                verts[i][2] = 0.0;
        }
        vert_cnt = 4;

        VEC_SUB( vec, verts[1], verts[0] );
        leng[0] = VEC_LENGTH( vec );
        VEC_SUB( vec, verts[3], verts[0] );
        leng[1] = VEC_LENGTH( vec );
        radius = rfac * (leng[0] + leng[1]) / 2.0;

        /* Outline the element. */
        glDepthFunc( GL_ALWAYS );
        glColor3fv( hilite_col );
        glBegin( GL_LINE_LOOP );
        for ( i = 0; i < 4; i++ )
            glVertex3fv( verts[i] );
        glEnd();
        glDepthFunc( GL_LEQUAL );

        break;

    case G_TET:
        /* Highlight a tet element. */

        if ( result_has_superclass( analy->cur_result, G_TET, analy ) )
        {
#ifdef OLD_RES_BUFFERS
            val = analy->perform_unit_conversion
                  ? analy->tet_result[hilite_num] * analy->conversion_scale
                  + analy->conversion_offset
                  : analy->tet_result[hilite_num];
#endif
            val = analy->perform_unit_conversion
                  ? data_array[hilite_num] * analy->conversion_scale
                  + analy->conversion_offset
                  : data_array[hilite_num];
            sprintf( label, " %s %d (%.*e)", cname, hilite_label, fracsz,
                     val );
        }
        else
            sprintf( label, " %s %d", cname, hilite_label );

        get_tet_verts( hilite_num, p_mo_class, analy, verts );
        vert_cnt = 4;

        VEC_SUB( vec, verts[1], verts[0] );
        leng[0] = VEC_LENGTH( vec );
        VEC_SUB( vec, verts[2], verts[0] );
        leng[1] = VEC_LENGTH( vec );
        VEC_SUB( vec, verts[3], verts[0] );
        leng[2] = VEC_LENGTH( vec );
        radius = rfac * (leng[0] + leng[1] + leng[2]) / 3.0;

        /* Outline the element. */
        glDepthFunc( GL_ALWAYS );
        glColor3fv( hilite_col );
        glBegin( GL_LINES );
        for ( i = 0; i < 6; i++ )
        {
            glVertex3fv( verts[tet_edge_node_nums[i][0]] );
            glVertex3fv( verts[tet_edge_node_nums[i][1]] );
        }
        glEnd();
        glDepthFunc( GL_LEQUAL );
        break;

    case G_PARTICLE:
    case G_HEX:
        /* Highlight a hex element. */

        if ( analy->particle_nodes_enabled && is_particle_class(analy, p_mo_class->superclass, p_mo_class->short_name) )
        {
            pn_hilite = TRUE;

            if ( is_particle_class( analy, p_mo_class->superclass, p_mo_class->short_name ) && p_mo_class->superclass != G_HEX)
                val = get_free_node_result( analy, p_mo_class, hilite_num, &result_defined, &valid_free_node );
            else
                val = get_ml_result( analy, p_mo_class, hilite_num, &result_defined );

            if ( result_defined )
            {
                if ( analy->perform_unit_conversion )
                    val = val * analy->conversion_scale + analy->conversion_offset;
                sprintf( label, " %s %d (%.*e)", cname, hilite_label, fracsz, val );
            }
            else
                sprintf( label, " %s %d", cname, hilite_label );

            if ( analy->particle_nodes_enabled && is_particle_class( analy, p_mo_class->superclass, p_mo_class->short_name ) )
                get_node_vert_2d_3d( analy->hilite_ml_node, p_mo_class, analy, verts[0] );
            else
                get_node_vert_2d_3d( hilite_num, p_mo_class, analy, verts[0] );

            vert_cnt = 1;

            for ( i = 0; i < dim; i++ )
                leng[i] = analy->bbox[1][i] - analy->bbox[0][i];

            node_base_radius = 0.01 * (leng[0] + leng[1] + leng[2]) / 3.0;

            radius = node_base_radius*analy->free_nodes_scale_factor*1.2;
            break;
        }

        if ( result_has_superclass( analy->cur_result, G_HEX, analy ) )
        {
            val = analy->perform_unit_conversion
                  ? data_array[hilite_num] * analy->conversion_scale
                  + analy->conversion_offset
                  : data_array[hilite_num];
            sprintf( label, " %s %d (%.*e)", cname, hilite_label, fracsz,
                     val );
        }
        else
            sprintf( label, " %s %d", cname, hilite_label );

        get_hex_verts( hilite_num , p_mo_class, analy, verts );
        vert_cnt = 8;

        VEC_SUB( vec, verts[1], verts[0] );
        leng[0] = VEC_LENGTH( vec );
        VEC_SUB( vec, verts[3], verts[0] );
        leng[1] = VEC_LENGTH( vec );
        VEC_SUB( vec, verts[4], verts[0] );
        leng[2] = VEC_LENGTH( vec );
        radius = rfac * (leng[0] + leng[1] + leng[2]) / 3.0;

        /* Outline the element. */
        glDepthFunc( GL_ALWAYS );

        /* Set particle select color to a greyscale */
        glColor3fv( hilite_col );
        glBegin( GL_LINES );
        for ( i = 0; i < 12; i++ )
        {
            glVertex3fv( verts[edge_node_nums[i][0]] );
            glVertex3fv( verts[edge_node_nums[i][1]] );
        }
        glEnd();
        glDepthFunc( GL_LEQUAL );

        break;

    case G_PYRAMID:
        /* Highlight a pyramid element. */

        /*if ( analy->particle_nodes_enabled && is_particle_class(analy, p_mo_class->superclass, p_mo_class->short_name) )
        {
             pn_hilite = TRUE;

                 if ( is_particle_class( analy, p_mo_class->superclass, p_mo_class->short_name ) )
              val = get_free_node_result( analy, p_mo_class, hilite_num, &result_defined, &valid_free_node );
         else
              val = get_ml_result( analy, p_mo_class, hilite_num, &result_defined );

                 if ( result_defined )
         {
              if ( analy->perform_unit_conversion )
                   val = val * analy->conversion_scale + analy->conversion_offset;
              sprintf( label, " %s %d (%.*e)", cname, hilite_label, fracsz, val );
         }
         else
              sprintf( label, " %s %d", cname, hilite_label );

         if ( analy->particle_nodes_enabled && is_particle_class( analy, p_mo_class->superclass, p_mo_class->short_name ) )
              get_node_vert_2d_3d( analy->hilite_ml_node, p_mo_class, analy, verts[0] );
         else
              get_node_vert_2d_3d( hilite_num, p_mo_class, analy, verts[0] );

         vert_cnt = 1;

                 for ( i = 0; i < dim; i++ )
                       leng[i] = analy->bbox[1][i] - analy->bbox[0][i];

                 node_base_radius = 0.01 * (leng[0] + leng[1] + leng[2]) / 3.0;

                 radius = node_base_radius*analy->free_nodes_scale_factor*1.2;
         break;
        } */

        if ( result_has_superclass( analy->cur_result, G_PYRAMID, analy ) )
        {
            val = analy->perform_unit_conversion
                  ? data_array[hilite_num] * analy->conversion_scale
                  + analy->conversion_offset
                  : data_array[hilite_num];
            sprintf( label, " %s %d (%.*e)", cname, hilite_label, fracsz,
                     val );
        }
        else
            sprintf( label, " %s %d", cname, hilite_label );

        get_pyramid_verts( hilite_num, p_mo_class, analy, verts );
        vert_cnt = 5;

        VEC_SUB( vec, verts[1], verts[0] );
        leng[0] = VEC_LENGTH( vec );
        VEC_SUB( vec, verts[3], verts[0] );
        leng[1] = VEC_LENGTH( vec );
        VEC_SUB( vec, verts[4], verts[0] );
        leng[2] = VEC_LENGTH( vec );
        radius = rfac * (leng[0] + leng[1] + leng[2]) / 3.0;

        /* Outline the element. */
        glDepthFunc( GL_ALWAYS );

        /* Set particle select color to a greyscale */
        glColor3fv( hilite_col );
        glBegin( GL_LINES );
        for ( i = 0; i < 8; i++ )
        {
            glVertex3fv( verts[pyramid_edge_node_nums[i][0]] );
            glVertex3fv( verts[pyramid_edge_node_nums[i][1]] );
        }
        glEnd();
        glDepthFunc( GL_LEQUAL );

        break;
    case G_SURFACE:
        /* Highlight a surface facet. */
        p_surface = p_mo_class->objects.surfaces;

        if ( show_label )
        {
            if ( result_has_superclass( analy->cur_result, G_SURFACE, analy ) )
            {
                val = analy->perform_unit_conversion
                      ? data_array[hilite_num]
                      * analy->conversion_scale + analy->conversion_offset
                      : data_array[hilite_num];
                sprintf( label, " %s %d (%.*e)", cname, hilite_num,
                         fracsz, val );
            }
            else
                sprintf( label, " %s %d", cname, hilite_num );
        }

        if ( dim == 3 )
            get_surface_verts_3d( hilite_num, p_surface[0].nodes,
                                  p_mo_class->mesh_id, analy, verts );
        else
        {
            get_surface_verts_2d( hilite_num, p_mo_class, analy, verts );
            for ( i = 0 ; i < 4; i++ )
                verts[i][2] = 0.0;
        }
        vert_cnt = 4;

        VEC_SUB( vec, verts[1], verts[0] );
        leng[0] = VEC_LENGTH( vec );
        VEC_SUB( vec, verts[3], verts[0] );
        leng[1] = VEC_LENGTH( vec );
        radius = rfac * (leng[0] + leng[1]) / 2.0;

        /* Outline the facet. */
        glDepthFunc( GL_ALWAYS );
        glColor3fv( hilite_col );
        glBegin( GL_LINE_LOOP );
        for ( i = 0; i < 4; i++ )
            glVertex3fv( verts[i] );
        glEnd();
        glDepthFunc( GL_LEQUAL );

        break;

    default:
        /* Don't visually represent other superclasses. */
        return;
    }

    /* Draw spheres at the nodes of the element. */
    if ( analy->dimension == 3 && !particle_hilite )
    {
        if ( v_win->lighting )
            glEnable( GL_LIGHTING );
        glEnable( GL_COLOR_MATERIAL );

        if ( pn_hilite )
            glColor3fv( particle_select_col ); /* change to 'hilite_col' in a future release when
						 * we add a setcol command to modify particle hilite
						 * color.
						 */
        else
            glColor3fv( hilite_col );

        if ( p_mo_class->superclass != G_SURFACE )
        {
            for ( i = 0; i < vert_cnt; i++ )
                draw_sphere( verts[i], radius, 1);
        }
        else
        {
            p_surface = p_mo_class->objects.surfaces;
            if ( dim == 3 )
                get_surface_verts_3d( hilite_num, p_surface[0].nodes,
                                      p_mo_class->mesh_id, analy, verts );
            else
            {
                get_surface_verts_2d( hilite_num, p_mo_class, analy, verts );
                for ( i = 0 ; i < 4; i++ )
                    verts[i][2] = 0.0;
            }
            vert_cnt = 4;
            for ( i = 0; i < vert_cnt; i++ )
                draw_sphere( verts[i], radius, 1);
        }

        glDisable( GL_COLOR_MATERIAL );
        if ( v_win->lighting )
            glDisable( GL_LIGHTING );
    }

    /* Get position for the label. */
    if ( !particle_hilite )
    {
        for ( j = 0; j < 3; j++ )
            vec[j] = 0.0;
        for ( i = 0; i < vert_cnt; i++ )
            for ( j = 0; j < 3; j++ )
                vec[j] += verts[i][j] / vert_cnt;
    }
    else
    {
        if ( analy->dimension == 3 )
        {
            GVec3D *parts3d;
            parts3d = analy->state_p->particles.particles3d;
            for ( i = 0; i < 3; i++ )
                vec[i] = parts3d[hilite_num][i];
        }
        else
        {
            GVec2D *parts2d;
            parts2d = analy->state_p->particles.particles2d;
            for ( i = 0; i < 2; i++ )
                vec[i] = parts2d[hilite_num][i];
            vec[2] = 0.0;
        }
    }

    /* Draw the element label, it goes on top of everything. */
    if ( show_label )
    {
        glDepthFunc( GL_ALWAYS );
        glColor3fv( v_win->text_color );
        draw_3d_text( vec, label, FALSE );
        glDepthFunc( GL_LEQUAL );
    }
}



/************************************************************
 * TAG( draw_class_numbers )
 *
 * This routine numbers superclasses specified by the user.
 * For the volume elements, only the visible faces are labeled.
 */
static void
draw_class_numbers( Analysis *analy )
{
    MO_class_data *p_mo_class, **p_classes;
    List_head *p_lh;
    Mesh_data *p_mesh;
    int quant, dim, fcnt, *node_nums;
    Bool_type verts_ok;
    float pt[3], verts[8][3];
    float v1[3], v2[3], nor[3], ray[3], tverts[4][3];
    char label[64];
    int i, j, k, l, p, n, offset;
    Visibility_data *p_vd;
    int *face_el, *face_fc, *ndflag;
    unsigned char *hide_mat;
    int *mat;
    int rel_node_qty[QTY_SCLASS] = { 0, 1, 2, 2, 3, 4, 4, 0, 0, 8, 0, 0, 0, 1 };
    int elem_sclasses[] = { G_TRUSS, G_BEAM, G_TRI, G_QUAD, G_TET, G_HEX, G_PARTICLE };

    int class_label;

    dim = analy->dimension;
    p_mesh = MESH_P( analy );
    view_transf_mat( &cur_view_mat );

    glColor3fv( v_win->text_color );
    glDepthFunc( GL_ALWAYS );

    /*
     * Label selected or all superclasses, as requested by the user.
     */

    for ( l=0;
            l<analy->classqty;
            l++ )
    {
        p_mo_class = analy->classarray[l];
        quant = p_mo_class->qty;

        switch( p_mo_class->superclass )
        {
        case G_UNIT:
            popup_dialog( INFO_POPUP,
                          "G_UNIT:unsupported superclass.");
            break;
        case G_MAT:
            popup_dialog( INFO_POPUP,
                          "G_MAT:unsupported superclass.");
            break;
        case G_MESH:
            popup_dialog( INFO_POPUP,
                          "G_MESH:unsupported superclass.");
            break;
        case G_PYRAMID:
            popup_dialog( INFO_POPUP,
                          "G_PYRAMID:unsupported superclass.");
            break;
        case G_WEDGE:
            popup_dialog( INFO_POPUP,
                          "G_WEDGE:unsupported superclass.");
            break;
        case G_NODE:
        {
            /* Mark the nodes to be labeled.*/

            ndflag = NEW_N( int, quant, "Temp node label flags" );
            for ( n = 0; n < sizeof( elem_sclasses ) / sizeof( elem_sclasses[0] ); n++ )
            {
                p_lh = p_mesh->classes_by_sclass + elem_sclasses[n];
                if ( p_lh->qty != 0 )
                {
                    p_classes = (MO_class_data **) p_lh->list;
                    for ( j = 0; j < p_lh->qty; j++ )
                    {
                        p = rel_node_qty[p_classes[j]->superclass];
                        node_nums = p_classes[j]->objects.elems->nodes;
                        switch( p_classes[j]->superclass )
                        {
                        case G_TRUSS:
                        case G_BEAM:
                        case G_TRI:
                        case G_QUAD:
                        {
                            hide_mat = MESH_P(analy)->hide_material;
                            mat = p_classes[j]->objects.elems->mat;


                            for ( i = 0; i < p_classes[j]->qty; i++ )
                                if ( !(hide_mat[mat[i]]) )
                                    for ( k=0; k < p; k++ )
                                    {
                                        if ( p_classes[j]->superclass == G_BEAM )
                                            ndflag[ node_nums[i*(p+1)+k]]=1;
                                        else
                                            ndflag[ node_nums[i*p+k]]=1;
                                    }
                        }
                        break;

                        case G_TET:
                        {
                            hide_mat = MESH_P(analy)->hide_material;
                            mat = p_classes[j]->objects.elems->mat;

                            /*Only nodes on visible faces are labeled.*/

                            p_vd = p_classes[j]->p_vis_data;
                            face_el = p_vd->face_el;
                            face_fc = p_vd->face_fc;
                            fcnt = p_vd->face_cnt;

                            for ( i = 0; i < fcnt; i++ )
                            {
                                if ( !(hide_mat[mat[face_el[i]]]) && (dim == 3) )
                                    get_tet_face_verts( face_el[i], face_fc[i],
                                                        p_classes[j], analy, verts );

                                /* Cull back faces by hand. */

                                for ( k = 0; k < 3; k++ )
                                    point_transform( tverts[k], verts[k], &cur_view_mat );
                                VEC_SUB( v1, tverts[1], tverts[0] );
                                VEC_SUB( v2, tverts[2], tverts[0] );
                                VEC_CROSS( nor, v1, v2 );
                                if ( v_win->orthographic )
                                {
                                    VEC_SET( ray, 0.0, 0.0, 1.0 );
                                }
                                else
                                {
                                    for ( k = 0; k < 3; k++ )
                                        ray[k] = -1.0/3.0 * ( tverts[0][k] + tverts[1][k]
                                                              + tverts[2][k] );
                                }
                                if ( VEC_DOT( ray, nor ) >= 0.0 )
                                {
                                    for ( k = 0; k < 3; k++ )
                                    {
                                        offset = tet_fc_nd_nums[face_fc[i]][k];
                                        ndflag[node_nums[offset+face_el[i]*4]] = 1;
                                    }
                                }
                            }
                        }
                        break;

                        case G_HEX:
                        {
                            hide_mat = MESH_P(analy)->hide_material;
                            mat = p_classes[j]->objects.elems->mat;

                            /* Only nodes visible faces are labeled. */

                            p_vd = p_classes[j]->p_vis_data;
                            face_el = p_vd->face_el;
                            face_fc = p_vd->face_fc;
                            fcnt = p_vd->face_cnt;

                            for ( i = 0; i < fcnt; i++ )
                                if ( !(hide_mat[mat[face_el[i]]]) && (dim == 3) )
                                {
                                    get_hex_face_verts( face_el[i], face_fc[i],
                                                        p_classes[j], analy, verts );

                                    /*Cull back faces by hand.*/

                                    for ( k = 0; k < 4; k++ )
                                        point_transform( tverts[k], verts[k],
                                                         &cur_view_mat );
                                    VEC_SUB( v1, tverts[2], tverts[0] );
                                    VEC_SUB( v2, tverts[3], tverts[1] );
                                    VEC_CROSS( nor, v1, v2 );
                                    if ( v_win->orthographic )
                                    {
                                        VEC_SET( ray, 0.0, 0.0, 1.0 );
                                    }
                                    else
                                    {
                                        for ( k = 0; k < 3; k++ )
                                            ray[k] = -0.25 * ( tverts[0][k] + tverts[1][k] +
                                                               tverts[2][k] + tverts[3][k] );
                                    }
                                    if ( VEC_DOT( ray, nor ) >= 0.0 )
                                    {
                                        for ( k = 0; k < 4; k++ )
                                        {
                                            offset = fc_nd_nums[face_fc[i]][k];
                                            ndflag[node_nums[offset+face_el[i]*8]]=1;
                                        }
                                    }
                                }
                        }
                        break;

                        default:
                            return;
                        }
                    }
                }
            }

            /* Draw the node labels. */

            for ( i = 0; i < quant; i++ )
                if ( ndflag[i] )
                {
                    get_node_vert_2d_3d( i, p_mo_class, analy, pt );

                    class_label = get_class_label( p_mo_class, i );

                    sprintf( label, "%d", class_label );
                    draw_3d_text( pt, label, TRUE );
                }
            free( ndflag );
        }
        break;

        case G_TRUSS:
        case G_BEAM:
        {
            hide_mat = MESH_P(analy)->hide_material;
            mat = p_mo_class->objects.elems->mat;

            for ( i = 0; i < quant; i++ )
            {
                if ( !(hide_mat[mat[i]]) )
                {
                    verts_ok = get_beam_verts_2d_3d( i, p_mo_class,
                                                     analy, verts );
                    if ( !verts_ok )
                    {
                        popup_dialog( WARNING_POPUP,
                                      "%d ignored; bad node(s).\n", i+1 );
                        return;
                    }
                    if ( dim == 2 )
                        verts[0][2] = 0.0;
                    for ( k=0; k < 3; k++ )
                        pt[k] = .5*(verts[0][k]+verts[1][k]);

                    class_label = get_class_label( p_mo_class, i );

                    sprintf( label, "%d", class_label );
                    draw_3d_text( pt, label, TRUE );
                }
            }
        }
        break;

        case G_TRI:
        {
            hide_mat = MESH_P(analy)->hide_material;
            mat = p_mo_class->objects.elems->mat;

            for ( i = 0; i < quant; i++ )
            {
                if ( !(hide_mat[mat[i]]) )
                {
                    if ( dim == 3 )
                        get_tri_verts_3d( i, p_mo_class, analy,
                                          verts );
                    else
                    {
                        get_tri_verts_2d( i, p_mo_class, analy, verts );
                        for ( k=0; k<3; k++ )
                            verts[i][2] = 0.0;
                    }
                    for ( k=0; k < 3; k++ )
                    {
                        pt[k] = .5*(verts[0][k]+verts[2][k]);
                    }

                    class_label = get_class_label( p_mo_class, i );

                    sprintf( label, "%d", class_label );
                    draw_3d_text( pt, label, TRUE );
                }
            }
        }
        break;

        case G_QUAD:
        {
            hide_mat = MESH_P(analy)->hide_material;
            mat = p_mo_class->objects.elems->mat;

            for ( i = 0; i < quant; i++ )
            {
                if ( !(hide_mat[mat[i]]) )
                {
                    if ( dim == 3 )
                        get_quad_verts_3d( i, p_mo_class->objects.elems->nodes,
                                           p_mo_class->mesh_id, analy, verts );
                    else
                    {
                        get_quad_verts_2d( i, p_mo_class, analy,
                                           verts );
                        for ( k=0; k<4; k++ )
                            verts[i][2] = 0.0;
                    }
                    for ( k=0; k < 3; k++ )
                        pt[k] = .5*(verts[0][k]+verts[2][k]);

                    class_label = get_class_label( p_mo_class, i+1 );

                    sprintf( label, "%d", class_label );
                    draw_3d_text( pt, label, TRUE );
                }
            }
        }
        break;

        case G_TET:
        {
            hide_mat = MESH_P(analy)->hide_material;
            mat = p_mo_class->objects.elems->mat;

            /*Only visible faces are labeled.*/

            p_vd = p_mo_class->p_vis_data;
            face_el = p_vd->face_el;
            face_fc = p_vd->face_fc;
            fcnt = p_vd->face_cnt;

            for ( i = 0; i < fcnt; i++ )
            {
                if ( !(hide_mat[mat[face_el[i]]]) && (dim == 3) )
                {
                    get_tet_face_verts( face_el[i], face_fc[i],
                                        p_mo_class, analy, verts );

                    /*Cull back faces by hand.*/

                    for ( j = 0; j < 3; j++ )
                        point_transform( tverts[j], verts[j], &cur_view_mat );
                    VEC_SUB( v1, tverts[1], tverts[0] );
                    VEC_SUB( v2, tverts[2], tverts[0] );
                    VEC_CROSS( nor, v1, v2 );
                    if ( v_win->orthographic )
                    {
                        VEC_SET( ray, 0.0, 0.0, 1.0 );
                    }
                    else
                    {
                        for ( k = 0; k < 3; k++ )
                            ray[k] = -1.0/3.0 * ( tverts[0][k] + tverts[1][k] +
                                                  tverts[2][k] );
                    }
                    if ( VEC_DOT( ray, nor ) >= 0.0 )
                    {
                        for ( k=0; k < 3; k++ )
                            pt[k] = .5*(verts[0][k]+verts[3][k]);

                        class_label = get_class_label( p_mo_class, i);

                        sprintf( label, "%d", class_label );
                        draw_3d_text( pt, label, TRUE );
                    }
                }
            }
        }
        /* The old way - put the label in the middle of the element,
         * even if it is invisible.

           for ( i = 0; i < quant; i++ )
               if ( !(hide_mat[mat[i]]) && (dim == 3) )
               {
                   get_tet_verts( i, p_mo_class, analy, verts );
                   for ( k=0; k < 3; k++ )
                         pt[k] = .5*(verts[0][k]+verts[3][k]);

        class_label = get_class_label( p_mo_class, i+1 );

        sprintf( label, "%d", class_label );
        draw_3d_text( pt, label, TRUE );
               }
        */
        break;

        case G_HEX:
        {
            hide_mat = MESH_P(analy)->hide_material;
            mat = p_mo_class->objects.elems->mat;

            /* Display numbers for particles which have no faces.*/

            if (is_particle_class( analy, p_mo_class->superclass, p_mo_class->short_name ))
            {
                for ( i = 0;
                        i < quant;
                        i++ )
                    if ( !(hide_mat[mat[i]]) && (dim == 3) )
                    {
                        get_hex_verts( i, p_mo_class, analy, verts );
                        for ( k=0;
                                k < 3;
                                k++ )
                            pt[k] = .5*(verts[0][k]+verts[6][k]);

                        class_label = get_class_label( p_mo_class, i );

                        sprintf( label, "%d", class_label );

                        draw_3d_text( pt, label, TRUE );
                    }
            }
            else
            {
                /* Only visible faces are labeled.*/
                p_vd = p_mo_class->p_vis_data;
                face_el = p_vd->face_el;
                face_fc = p_vd->face_fc;
                fcnt = p_vd->face_cnt;

                for ( i = 0; i < fcnt; i++ )
                {
                    if ( !(hide_mat[mat[face_el[i]]]) && (dim == 3) )
                    {
                        get_hex_face_verts( face_el[i], face_fc[i],
                                            p_mo_class, analy, verts );

                        /* Cull back faces by hand.*/

                        for ( j = 0; j < 4; j++ )
                            point_transform( tverts[j], verts[j], &cur_view_mat );
                        VEC_SUB( v1, tverts[2], tverts[0] );
                        VEC_SUB( v2, tverts[3], tverts[1] );
                        VEC_CROSS( nor, v1, v2 );
                        if ( v_win->orthographic )
                        {
                            VEC_SET( ray, 0.0, 0.0, 1.0 );
                        }
                        else
                        {
                            for ( k = 0; k < 3; k++ )
                                ray[k] = -0.25 * ( tverts[0][k] + tverts[1][k] +
                                                   tverts[2][k] + tverts[3][k] );
                        }
                        if ( VEC_DOT( ray, nor ) >= 0.0 )
                        {
                            for ( k=0; k < 3; k++ )
                                pt[k] = .5*(verts[0][k]+verts[2][k]);

                            class_label = get_class_label( p_mo_class, face_el[i] );

                            sprintf( label, "%d", class_label);
                            draw_3d_text( pt, label, TRUE );
                        }
                    }
                }
            }
        }
        break;

        case G_PARTICLE:
        {
            hide_mat = MESH_P(analy)->hide_material;
            mat = p_mo_class->objects.elems->mat;

            /* Display numbers for particles which have no faces.*/

            for ( i = 0;
                    i < quant;
                    i++ )
                if ( !(hide_mat[mat[i]]) && (dim == 3) )
                {
                    get_particle_verts( i, p_mo_class, analy, verts );
                    for ( k=0;
                            k < 3;
                            k++ )
                        pt[k] = verts[0][k];

                    class_label = get_class_label( p_mo_class, i+1 );

                    sprintf( label, "%d", class_label );

                    draw_3d_text( pt, label, TRUE );
                }
        }
        break;

        default:
            return;
        }
    }

    glDepthFunc( GL_LEQUAL );
}


/**/
/* Not correct for 2D, but might work. */
/************************************************************
 * TAG( draw_bbox )
 *
 * Draw the bounding box of the mesh.
 */
static void
draw_bbox( float bbox[2][3] )
{
    float verts[8][3];
    int i;

    get_verts_of_bbox( bbox, verts );

    antialias_lines( TRUE, FALSE );

    glColor3fv( v_win->foregrnd_color );

    glBegin( GL_LINES );
    for ( i = 0; i < 12; i++ )
    {
        glVertex3fv( verts[edge_node_nums[i][0]] );
        glVertex3fv( verts[edge_node_nums[i][1]] );
    }
    glEnd();

    antialias_lines( FALSE, 0 );
}


/************************************************************
 * TAG( draw_extern_polys )
 *
 * Draw polygons that were read from an external file.
 */
static void
draw_extern_polys( Analysis *analy )
{
    Surf_poly *poly;
    float cols[4][4];
    float res[4];
    int matl, i;
    Mesh_data *p_mesh;

    p_mesh = MESH_P( analy );

    /* Dummy values. */
    for ( i = 0; i < 4; i++ )
        res[i] = 0.0;

    /* Enable color to change AMBIENT and DIFFUSE properties. */
    glEnable( GL_COLOR_MATERIAL );

    /* Set up for polygon drawing. */
    begin_draw_poly( analy );

    for ( poly = analy->extern_polys; poly != NULL; poly = poly->next )
    {
        matl = poly->mat;


        if ( v_win->mesh_materials.current_index != matl )
            change_current_color_property( &v_win->mesh_materials, matl );


        for ( i = 0; i < poly->cnt; i++ )
            hvec_copy( cols[i], v_win->mesh_materials.diffuse[matl] );

        draw_poly( poly->cnt, poly->vtx, poly->norm, cols, res, matl, p_mesh,
                   analy, FALSE );
    }

    /* Back to the default. */
    end_draw_poly( analy );

    glDisable( GL_COLOR_MATERIAL );
}


/************************************************************
 * TAG( draw_ref_polys )
 *
 * Draw reference surface polygons.
 */
static void
draw_ref_polys( Analysis *analy )
{
    extern float crease_threshold;
    Ref_poly *poly;
    Bool_type show_result, showgs=FALSE;
    float v1[3], v2[3], f_norm[3], n_norm[3], dot;
    float verts[4][3];
    float norms[4][3];
    float cols[4][4];
    float res[4];
    float *data_buffer;
    int matl;
    int i, j;
    GVec3D *nodes3d;
    unsigned char *disable_mtl;
    Mesh_data *p_mesh;

    nodes3d = analy->state_p->nodes.nodes3d;
    p_mesh = MESH_P( analy );
    disable_mtl = p_mesh->disable_material;
    show_result = analy->result_on_refs
                  && result_has_superclass( analy->cur_result, G_HEX, analy );
    data_buffer = NODAL_RESULT_BUFFER( analy );

    /* Enable color to change AMBIENT and DIFFUSE properties. */
    glEnable( GL_COLOR_MATERIAL );

    /* Set up for polygon drawing. */
    begin_draw_poly( analy );

    for ( poly = analy->ref_polys; poly != NULL; poly = poly->next )
    {
        /* Flip polygon orientation by reversing vertex order. */
        for ( i = 0; i < 4; i++ )
            for ( j = 0; j < 3; j++ )
                verts[3-i][j] = nodes3d[poly->nodes[i]][j];

        /* Get polygon average normal. */
        VEC_SUB( v1, verts[2], verts[0] );
        VEC_SUB( v2, verts[3], verts[1] );
        VEC_CROSS( f_norm, v1, v2 );
        vec_norm( f_norm );

        /* Get polygon vertex normals using edge detection. */
        for ( i = 0; i < 4; i++ )
        {
            for ( j = 0; j < 3; j++ )
                n_norm[j] = analy->node_norm[j][poly->nodes[i]];
            vec_norm( n_norm );

            dot = VEC_DOT( f_norm, n_norm );

            if ( dot < crease_threshold )
            {
                /* Edge detected, use face normal. */
                for ( j = 0; j < 3; j++ )
                    norms[i][j] = f_norm[j];
            }
            else
            {
                /* No edge, use node normal. */
                for ( j = 0; j < 3; j++ )
                    norms[i][j] = n_norm[j];
            }
        }

        matl = poly->mat;


        if ( v_win->mesh_materials.current_index != matl )
            change_current_color_property( &v_win->mesh_materials, matl );


        /* Colorflag is TRUE if displaying a result, otherwise
         * polygons are drawn in the material color.
         */
        colorflag = show_result && !disable_mtl[matl];
        if ( analy->material_greyscale && disable_mtl[matl] )
            showgs = TRUE;
        else
            showgs = FALSE;

        /*
         * analy->result is gone in favor of node class-specific data buffers,
         * but we have no way of knowing which is the right node class
         * data buffer to read data values from.  For now we will punt and
         * always use NODAL_RESULT_BUFFER.  When element mesh object class data
         * is extended to include a reference to the nodal class that the
         * elements are defined on, we should modify reference surface logic
         * so that a surface can only be defined on a single element class,
         * which will then make obvious the correct nodal data buffer to use
         * and eliminate the possibility of having to access multiple nodal
         * data buffers for a single surface.
         */
        for ( i = 0; i < 4; i++ )
        {
            if ( analy->interp_mode == GOOD_INTERP )
            {
                res[3-i] = data_buffer[poly->nodes[i]];
            }
            else
            {
                /* We don't know the element number, so there's no way
                 * to implement NO_INTERP on reference faces -- the
                 * vertex values are always interpolated.
                 */
                color_lookup( cols[3-i], data_buffer[poly->nodes[i]],
                              analy->result_mm[0], analy->result_mm[1],
                              analy->zero_result, matl, analy->logscale,
                              analy->material_greyscale );
            }
        }

        draw_poly( 4, verts, norms, cols, res, matl, p_mesh, analy, FALSE );
    }

    /* Back to the defaults. */
    end_draw_poly( analy );
    glDisable( GL_COLOR_MATERIAL );
}


/************************************************************
 * TAG( draw_vec_result_3d )
 *
 * Driver for vector drawing on 3D vgrids.
 */
static void
draw_vec_result_3d( Analysis *analy )
{
    float vec_leng, vmin, vmax, diff;
    float denom;
    int i;
    float scl_max;
    Mesh_data *p_mesh;
    int qty_classes;
    MO_class_data **mo_classes;

    /* Cap for colorscale interpolation. */
    scl_max = SCL_MAX;

    /* Make the max vector length about VEC_3D_LENGTH pixels long. */

    /*
     * This (unpublished) option added to scale and color vectors by
     * the current result.  This must be exercised very carefully and
     * in general should only be used when the result is the magnitude
     * function of the vector quantity.
     */
    if ( analy->scale_vec_by_result && analy->cur_result != NULL )
    {
        vmax = analy->result_mm[1];
        vmin = analy->result_mm[0];
    }
    else
    {
        /* Original assignments. */
        vmax = analy->vec_max_mag;
        vmin = analy->vec_min_mag;
    }

    diff = vmax - vmin;
    denom = v_win->vp_height * v_win->bbox_scale * vmax;

    vec_leng = ( denom != 0.0 ) ? analy->vec_scale * VEC_3D_LENGTH / denom
               : 0.0;

    p_mesh = MESH_P( analy );

    /* Hex element classes. */
    qty_classes = p_mesh->classes_by_sclass[G_HEX].qty;
    mo_classes = (MO_class_data **) p_mesh->classes_by_sclass[G_HEX].list;
    for ( i = 0; i < qty_classes; i++ )
        if ( mo_classes[i]->vector_pts != NULL )
            draw_vecs_3d( mo_classes[i]->vector_pts, scl_max, vmin, vmax, diff,
                          vec_leng, p_mesh, analy );

    /* Tet element classes. */
    qty_classes = p_mesh->classes_by_sclass[G_TET].qty;
    mo_classes = (MO_class_data **) p_mesh->classes_by_sclass[G_TET].list;
    for ( i = 0; i < qty_classes; i++ )
        if ( mo_classes[i]->vector_pts != NULL )
            draw_vecs_3d( mo_classes[i]->vector_pts, scl_max, vmin, vmax, diff,
                          vec_leng, p_mesh, analy );

}


/************************************************************
 * TAG( draw_vecs_3d )
 *
 * Draw a vector result with a bunch of little line segment
 * glyphs.
 */
static void
draw_vecs_3d( Vector_pt_obj *vec_list, float scl_max, float vmin, float vmax,
              float diff, float vec_leng, Mesh_data *p_mesh,  Analysis *analy )
{
    Vector_pt_obj *pt;
    float vmag;
    float tmp[3], pts[6];
    float leng[3], radius;
    int cnt, i;

    antialias_lines( TRUE, analy->z_buffer_lines );

    /* Draw the vectors. */
    for ( pt = vec_list, cnt = 0; pt != NULL; pt = pt->next, cnt++ )
    {
        if ( analy->vec_col_set )
            glColor3fv( v_win->vector_color );
        else
        {
            /* Set the color based on the magnitude.  Avoid divide-by-zero. */
            if ( APX_EQ( diff, 0.0 ) )
                vmag = 0.0;
            else
                vmag = (sqrt( (double) (VEC_DOT( pt->vec, pt->vec )) ) - vmin)
                       / diff;

            i = (int)(vmag * scl_max) + 1;
            glColor3fv( v_win->colormap[i] );
        }

        /* Scale the vector, using zero as the min and vmax as the max. */
        VEC_COPY( pts, pt->pt );
        VEC_ADDS( tmp, vec_leng, pt->vec, pt->pt );
        pts[3] = tmp[0];
        pts[4] = tmp[1];
        pts[5] = tmp[2];

        draw_line( 2, pts, -1, p_mesh, analy, FALSE, NULL );
    }

    antialias_lines( FALSE, 0 );

    /* Get sphere radius for base of vector. */
    for ( i = 0; i < 3; i++ )
        leng[i] = analy->bbox[1][i] - analy->bbox[0][i];
    radius = 0.01 * (leng[0] + leng[1] + leng[2]) / 3.0;
    radius *= 1.0 / v_win->scale[0];

    /* Draw spheres. */
    if ( analy->show_vector_spheres )
    {
        glEnable( GL_COLOR_MATERIAL );
        if ( v_win->lighting )
            glEnable( GL_LIGHTING );

        for ( pt = vec_list; pt != NULL; pt = pt->next )
        {
            if ( analy->vec_col_set )
                glColor3fv( v_win->vector_color );
            else
            {
                /* Set the color based on the magnitude. */
                if ( APX_EQ( diff, 0.0 ) )
                    vmag = 0.0;
                else
                    vmag = (sqrt((double)(VEC_DOT(pt->vec, pt->vec)))-vmin) /
                           diff;
                i = (int)(vmag * scl_max) + 1;
                glColor3fv( v_win->colormap[i] );
            }
            draw_sphere( pt->pt, radius, 1);
        }

        glDisable( GL_COLOR_MATERIAL );
        if ( v_win->lighting )
            glDisable( GL_LIGHTING );
    }
}


/************************************************************
 * TAG( draw_vec_result_2d )
 *
 * Draw a vector result with a bunch of little line segment
 * glyphs.  For 2D vectors, we can try to draw the vector
 * heads.
 */
static void
draw_vec_result_2d( Analysis *analy )
{
    float vpx, vpy, pixsize, vec_leng, vmin, vmax;
    float diff;
    float headpts[3][3];
    int i;
    float scl_max;
    Mesh_data *p_mesh;
    int qty_classes;
    MO_class_data **mo_classes;

    /* Cap for colorscale interpolation. */
    scl_max = SCL_MAX;

    /* See if there are some to draw. */
    if ( !analy->have_grid_points )
        return;

    vpx = (float) v_win->vp_width;
    vpy = (float) v_win->vp_height;
    if ( vpx >= vpy )
        pixsize = 2.0 / ( vpy * v_win->bbox_scale * v_win->scale[1] );
    else
        pixsize = 2.0 / ( vpx * v_win->bbox_scale * v_win->scale[0] );

    /* Make the max vector length about VEC_2D_LENGTH pixels long. */
    vmax = analy->vec_max_mag;
    vmin = analy->vec_min_mag;
    diff = vmax - vmin;
    vec_leng = ( vmax != 0.0 )
               ? analy->vec_scale * VEC_2D_LENGTH * pixsize / vmax
               : 0.0;

    /* Create the vector head points. */
    headpts[0][0] = 0.0;
    headpts[0][1] = -2.0*pixsize*analy->vec_head_scale;
    headpts[0][2] = 0.0;
    headpts[1][0] = 5.0*pixsize*analy->vec_head_scale;
    headpts[1][1] = 0.0;
    headpts[1][2] = 0.0;
    headpts[2][0] = 0.0;
    headpts[2][1] = 2.0*pixsize*analy->vec_head_scale;
    headpts[2][2] = 0.0;

    p_mesh = MESH_P( analy );

    /* Tri element classes. */
    qty_classes = p_mesh->classes_by_sclass[G_TRI].qty;
    mo_classes = (MO_class_data **) p_mesh->classes_by_sclass[G_TRI].list;
    for ( i = 0; i < qty_classes; i++ )
        if ( mo_classes[i]->vector_pts != NULL )
            draw_vecs_2d( mo_classes[i]->vector_pts, headpts, scl_max,
                          vmin, vmax, diff, vec_leng, p_mesh, analy );

    /* Quad element classes. */
    qty_classes = p_mesh->classes_by_sclass[G_QUAD].qty;
    mo_classes = (MO_class_data **) p_mesh->classes_by_sclass[G_QUAD].list;
    for ( i = 0; i < qty_classes; i++ )
        if ( mo_classes[i]->vector_pts != NULL )
            draw_vecs_2d( mo_classes[i]->vector_pts, headpts, scl_max,
                          vmin, vmax, diff, vec_leng, p_mesh, analy );

    /* Surface element classes. */
    qty_classes = p_mesh->classes_by_sclass[G_SURFACE].qty;
    mo_classes = (MO_class_data **) p_mesh->classes_by_sclass[G_SURFACE].list;
    for ( i = 0; i < qty_classes; i++ )
        if ( mo_classes[i]->vector_pts != NULL )
            draw_vecs_2d( mo_classes[i]->vector_pts, headpts, scl_max,
                          vmin, vmax, diff, vec_leng, p_mesh, analy );
}


/************************************************************
 * TAG( draw_vecs_2d )
 *
 * Draw a vector result with a bunch of little line segment
 * glyphs.  For 2D vectors, we can try to draw the vector
 * heads.
 */
static void
draw_vecs_2d( Vector_pt_obj *vec_list, float headpts[3][3], float scl_max,
              float vmin, float vmax, float diff, float vec_leng,
              Mesh_data *p_mesh,  Analysis *analy )
{
    Vector_pt_obj *pt;
    Transf_mat tmat;
    float vmag;
    float angle;
    float tmp[3], pts[6], res[3];
    float verts[4][3];
    float cols[4][4];
    int cnt, i;

    /* Dummy. */
    VEC_SET( res, 0.0, 0.0, 0.0 );

    antialias_lines( TRUE, analy->z_buffer_lines );

    /* Draw the vectors. */
    for ( pt = vec_list, cnt = 0; pt != NULL; pt = pt->next, cnt++ )
    {
        if ( analy->vec_col_set )
        {
            VEC_COPY( cols[0], v_win->vector_color );
            glColor3fv( cols[0] );
        }
        else
        {
            /* Set the color based on the magnitude.  Avoid divide-by-zero. */
            if ( APX_EQ( diff, 0.0 ) )
                vmag = 0.0;
            else
                vmag = (sqrt((double)(VEC_DOT(pt->vec, pt->vec))) - vmin)/diff;
            i = (int)(vmag * scl_max) + 1;
            VEC_COPY( cols[0], v_win->colormap[i] );
            glColor3fv( cols[0] );
        }

        /* Scale the vector, using zero as the min and vmax as the max. */
        VEC_COPY( pts, pt->pt );
        VEC_ADDS( tmp, vec_leng, pt->vec, pt->pt );
        pts[3] = tmp[0];
        pts[4] = tmp[1];
        pts[5] = tmp[2];

        draw_line( 2, pts, -1, p_mesh, analy, FALSE, NULL );

        /* Draw the vector head. */
        if ( DEF_GT( VEC_LENGTH( pt->vec ), 0.0 ) )
        {
            VEC_COPY( cols[1], cols[0] );
            VEC_COPY( cols[2], cols[0] );

            VEC_SCALE( tmp, vec_leng, pt->vec );
            angle = atan2( (double)tmp[1], (double)tmp[0] );

            mat_copy( &tmat, &ident_matrix );
            mat_translate( &tmat, VEC_LENGTH( tmp ), 0.0, 0.0 );
            mat_rotate_z( &tmat, angle );
            mat_translate( &tmat, pt->pt[0], pt->pt[1], pt->pt[2] );

            for ( i = 0; i < 3; i++ )
                point_transform( verts[i], headpts[i], &tmat );
            draw_poly_2d( 3, verts, cols, res, -1, p_mesh, analy );
        }
    }

    antialias_lines( FALSE, 0 );
}


/************************************************************
 * TAG( draw_node_vec_2d_3d )
 *
 * Draw a vector result at the grid nodes.
 */
static void
draw_node_vec_2d_3d( Analysis *analy )
{
    float *coords;
    Bool_type draw_heads;
    Transf_mat tmat;
    float *vec_result[3];
    float vec_leng, vmag, vmin, vmax, diff;
    float vpx, vpy, pixsize, angle;
    float vec[3], pts[6], tmp[3];
    float res[3], headpts[3][3], verts[4][3], cols[4][4];
    float leng[3];
    float radius;
    int i, j, dim;
    float scl_max;
    int node_qty, idx;
    Mesh_data *p_mesh;

    float *nodal_data, nodal_val=0., temp_val=0., val_range=0., range_percent=0.,
                       rmin=0., rmax=0.;

    /* Cap for colorscale interpolation. */
    scl_max = SCL_MAX;

    coords = analy->state_p->nodes.nodes;
    p_mesh = MESH_P( analy );
    node_qty = p_mesh->node_geom->qty;
    dim = analy->dimension;

    /* Load the three results into vector result array. */
    load_vec_result( analy, vec_result, &vmin, &vmax );

    nodal_data = NODAL_RESULT_BUFFER( analy );
    get_min_max( analy, FALSE, &rmin, &rmax );

    /* Make the max vector length about VEC_3D_LENGTH pixels long. */
    diff = vmax - vmin;
    vec_leng = ( vmax != 0.0 )
               ? analy->vec_scale * VEC_3D_LENGTH /
               (v_win->vp_height * v_win->bbox_scale * vmax)
               : 0.0;

    if ( dim == 2 )
        draw_heads = TRUE;
    else
        draw_heads = FALSE;

    /* Get vector head points. */
    if ( draw_heads )
    {
        /* Dummy. */
        VEC_SET( res, 0.0, 0.0, 0.0 );

        vpx = (float) v_win->vp_width;
        vpy = (float) v_win->vp_height;
        if ( vpx >= vpy )
            pixsize = 2.0 / ( vpy * v_win->bbox_scale * v_win->scale[1] );
        else
            pixsize = 2.0 / ( vpx * v_win->bbox_scale * v_win->scale[0] );

        /* Make the max vector length about VEC_2D_LENGTH pixels long. */
        vec_leng = ( vmax != 0.0 )
                   ? analy->vec_scale * VEC_2D_LENGTH * pixsize / vmax
                   : 0.0;

        /* Create the vector head points. */
        headpts[0][0] = 0.0;
        headpts[0][1] = -2.0*pixsize*analy->vec_head_scale;
        headpts[0][2] = 0.0;
        headpts[1][0] = 5.0*pixsize*analy->vec_head_scale;
        headpts[1][1] = 0.0;
        headpts[1][2] = 0.0;
        headpts[2][0] = 0.0;
        headpts[2][1] = 2.0*pixsize*analy->vec_head_scale;
        headpts[2][2] = 0.0;
    }
    else if ( analy->show_vector_spheres )
    {
        /* Get sphere radius for base of vector. */
        for ( i = 0; i < 3; i++ )
            leng[i] = analy->bbox[1][i] - analy->bbox[0][i];
        radius = 0.01 * (leng[0] + leng[1] + leng[2]) / 3.0;
        radius *= 1.0 / v_win->scale[0];

        /* Initialize lighting for sphere rendering. */
        glEnable( GL_COLOR_MATERIAL );
        if ( v_win->lighting )
            glEnable( GL_LIGHTING );
    }

    antialias_lines( TRUE, analy->z_buffer_lines );

    /* Draw the vectors. */
    for ( i = 0; i < node_qty; i++ )
    {
        /**/
        /*
                for ( j = 0; j < 3; j++ )
                {
                    pts[j] = coords[i * dim + j];
                    vec[j] = vec_result[j][i];
                }
        */
        idx = i * dim;
        pts[0] = coords[idx];
        vec[0] = vec_result[0][i];
        pts[1] = coords[idx + 1];
        vec[1] = vec_result[1][i];
        vec[2] = vec_result[2][i]; /* "vec" always defined with three comp's. */
        if ( dim == 3 )
            pts[2] = coords[idx + 2];
        else
            pts[2] = 0.0;

        if ( analy->vec_col_set )
        {
            VEC_COPY( cols[0], v_win->vector_color );
        }
        else
        {
            /* Set the color based on the magnitude.  Avoid divide-by-zero. */
            if ( APX_EQ( diff, 0.0 ) )
                vmag = 0.0;
            else
                vmag = (sqrt((double)(VEC_DOT(vec, vec))) - vmin)/diff;
            val_range = rmax - rmin;
            if ( nodal_data && val_range!= 0.0 )
            {
                nodal_val = nodal_data[i];
                temp_val  = nodal_val - rmin;
                range_percent = fabsf(nodal_val/val_range);
                j = (int) (range_percent*256.0);

                if ( j>255 )
                    j = (int)(vmag * scl_max) + 1;
            }
            else
                j = (int)(vmag * scl_max) + 1;

            VEC_COPY( cols[0], v_win->colormap[j] );
        }
        glColor3fv( cols[0] );

        /* Scale the vector, using zero as the min and vmax as the max. */
        VEC_ADDS( tmp, vec_leng, vec, pts );
        pts[3] = tmp[0];
        pts[4] = tmp[1];
        pts[5] = tmp[2];

        draw_line( 2, pts, -1, p_mesh, analy, FALSE, NULL );

        /* Draw the vector head. */
        if ( draw_heads )
        {
            if ( DEF_GT( VEC_LENGTH( vec ), 0.0 ) )
            {
                VEC_COPY( cols[1], cols[0] );
                VEC_COPY( cols[2], cols[0] );

                VEC_SCALE( tmp, vec_leng, vec );
                angle = atan2( (double)tmp[1], (double)tmp[0] );

                mat_copy( &tmat, &ident_matrix );
                mat_translate( &tmat, VEC_LENGTH( tmp ), 0.0, 0.0 );
                mat_rotate_z( &tmat, angle );
                mat_translate( &tmat, pts[0], pts[1], pts[2] );

                for ( j = 0; j < 3; j++ )
                    point_transform( verts[j], headpts[j], &tmat );
                draw_poly_2d( 3, verts, cols, res, -1, p_mesh, analy );
            }
        }
        else if ( analy->show_vector_spheres )
            draw_sphere( pts, radius, 1);
    }

    /*      len diam */
    VEC_COPY( cols[1], cols[0] );
    VEC_COPY( cols[2], cols[0] );

    VEC_SCALE( tmp, vec_leng, vec );
    angle = atan2( (double)tmp[1], (double)tmp[0] );

    mat_copy( &tmat, &ident_matrix );
    mat_translate( &tmat, VEC_LENGTH( tmp ), 0.0, 0.0 );
    mat_rotate_z( &tmat, angle );
    mat_translate( &tmat, pts[0], pts[1], pts[2] );
    DrawCone(.01*vec_leng, .01, cols, 8);

    antialias_lines( FALSE, 0 );

    if ( !draw_heads && analy->show_vector_spheres )
    {
        glDisable( GL_COLOR_MATERIAL );
        if ( v_win->lighting )
            glDisable( GL_LIGHTING );
    }

    /* Free temporary arrays. */
    for ( i = 0; i < 3; i++ )
        free( vec_result[i] );
}


/*
 * Carpetting routines
 *
 *     Carpeting reference surfaces
 *     Carpeting cut planes and isosurfaces
 *     Volume carpeting - regular sampling method
 *     Volume carpeting - density/sorting method (not implemented yet)
 */


/************************************************************
 * TAG( load_vec_result )
 *
 * Load in the current vector result.  Also, return the min and
 * max magnitude of the vector result.  The vec result array is
 * allocated by this routine and needs to be deallocated by the
 * caller.
 */
static void
load_vec_result( Analysis *analy, float *vec_result[3], float *vmin,
                 float *vmax )
{
    Result *tmp_res;
    float *tmp_data;
    Bool_type convert;
    float vec[3];
    float mag, min, max;
    int node_qty, i;
    float scale, offset;
    MO_class_data *p_node_class;

    p_node_class = MESH_P( analy )->node_geom;
    node_qty = p_node_class->qty;

    /*
     * Hack alert.  In replacing analy->result with node-class-specific
     * data buffers (and, eventually, permitting each element class to
     * specify a particular node class it's defined on), we've created the
     * _potential_ for impossible associations of data to form vectors.  We
     * have no more central array (analy->result) into which all data gets
     * funneled, and setting the temporary destination array here is wrong
     * because load_result() really needs to identify the correct destination
     * as it loops over subrecords and encounters the bound object classes.
     * This implementation will work as long as there is only a single node
     * class.  Hopefully there are some constraints that can be defined that
     * actually align with likely use cases that will limit the problem space
     * to make vector rendering tractable.  As is, it could be that each
     * vector component requires separate data arrays for each of the node
     * classes used by all of the element classes that support the vector.
     */

    for ( i = 0; i < 3; i++ )
        vec_result[i] = NEW_N( float, node_qty, "Vec result" );

    tmp_res = analy->cur_result;
    tmp_data = p_node_class->data_buffer;
    for ( i = 0; i < 3; i++ )
    {
        analy->cur_result = analy->vector_result[i];
        p_node_class->data_buffer = vec_result[i];
        load_result( analy, FALSE, TRUE, FALSE );
    }
    analy->cur_result = tmp_res;
    p_node_class->data_buffer = tmp_data;

    convert = analy->perform_unit_conversion;
    if ( convert )
    {
        scale = analy->conversion_scale;
        offset = analy->conversion_offset;
    }

    /* Get the min and max magnitude for the vector result. */
    vec[0] = vec_result[0][0];
    vec[1] = vec_result[1][0];
    vec[2] = vec_result[2][0];
    mag = VEC_DOT( vec, vec );
    min = max = mag;
    for ( i = 1; i < node_qty; i++ )
    {
        if ( convert )
        {
            vec_result[0][i] = vec_result[0][i] * scale + offset;
            vec_result[1][i] = vec_result[1][i] * scale + offset;
            vec_result[2][i] = vec_result[2][i] * scale + offset;
        }

        vec[0] = vec_result[0][i];
        vec[1] = vec_result[1][i];
        vec[2] = vec_result[2][i];
        mag = VEC_DOT( vec, vec );

        if ( mag > max )
            max = mag;
        else if ( mag < min )
            min = mag;
    }
    *vmin = sqrt( (double)min );
    *vmax = sqrt( (double)max );
}


#ifdef CARPET_STUFF

/************************************************************
 * TAG( find_front_faces )
 *
 * Determines which of the six faces of a 3D grid are front-
 * facing in the current view, and returns this information
 * in the front array.
 */
static void
find_front_faces( Bool_type front, Analysis *analy )
{
    static float face_norms[6][3] = { {-1.0,  0.0,  0.0},
        { 1.0,  0.0,  0.0},
        { 0.0, -1.0,  0.0},
        { 0.0,  1.0,  0.0},
        { 0.0,  0.0, -1.0},
        { 0.0,  0.0,  1.0}
    };
    Transf_mat inv_rot_mat;
    float eye[3], rot_eye[3];
    int i, j;

    /*
     * This doesn't take into account the look_at point or the
     * perspective, but we'll forget about that for now.
     */
    VEC_COPY( eye, v_win->look_from );

    /* Get the inverse of the view rotation matrix. */
    for ( i = 0; i < 4; i++ )
        for ( j = 0; j < 4; j++ )
            inv_rot_mat.mat[i][j] = v_win->rot_mat.mat[j][i];

    point_transform( rot_eye, eye, &inv_rot_mat );

    for ( i = 0; i < 6; i++ )
    {
        if ( VEC_DOT( face_norms[i], rot_eye ) > 0.0 )
            front[i] = TRUE;
        else
            front[i] = FALSE;
    }
}


/************************************************************
 * TAG( draw_carpet_vector )
 *
 * Draw a single carpet vector "hair".
 */
static void
draw_carpet_vector( float pt[3], float vec[3], float vec_leng, float vmin,
                    float vmax, Analysis *analy )
{
    float pts[2][3], cols[2][4];
    float vec2[3], vec3[3];
    float vmag, diff, scale;
    int i;
    float scl_max;

    /* Cap for colorscale interpolation. */
    scl_max = SCL_MAX;

    diff = vmax - vmin;

    if ( analy->vec_length_factor > 0.0 )
        vec_leng = vec_leng * vmax;

    /* Get vector color. */
    cols[0][3] = 1.0;
    cols[1][3] = 1.0;
    if ( analy->vec_col_set )
    {
        VEC_COPY( cols[0], v_win->vector_color );
    }
    else
    {
        /* Set the color based on the magnitude. */
        if ( APX_EQ( diff, 0.0 ) )
            vmag = 0.0;
        else
            vmag = (sqrt((double)VEC_DOT(vec, vec)) - vmin)/diff;
        i = (int)(vmag * scl_max) + 1;
        VEC_COPY( cols[0], v_win->colormap[i] );
    }
    if ( analy->vec_hd_col_set )
    {
        VEC_COPY( cols[1], v_win->vector_hd_color );
    }
    else
    {
        VEC_COPY( cols[1], cols[0] );
    }

    /* Jitter intensity of vector randomly. */
    if ( analy->vec_jitter_factor > 0.0 )
    {
        scale = analy->vec_jitter_factor * drand48() +
                1.0 - analy->vec_jitter_factor;
        VEC_SCALE( cols[0], scale, cols[0] );
        VEC_SCALE( cols[1], scale, cols[1] );
    }

    /* Modulate opacity based on importance. */
    if ( analy->vec_import_factor > 0.0 )
    {
        if ( APX_EQ( diff, 0.0 ) )
            vmag = 0.0;
        else
            vmag = (sqrt((double)VEC_DOT(vec, vec)) - vmin)/diff;

        /* The factor tells us the maximum opacity. */
        vmag = vmag * analy->vec_import_factor;

        /* Do a squared drop off. */
        cols[0][3] = vmag * vmag;
        cols[1][3] = cols[0][3];
    }

    /* Offset the vector scaling. */
    if ( analy->vec_length_factor > 0.0 )
    {
        vmag = sqrt((double)VEC_DOT(vec, vec));
        if ( vmag == 0.0 )
            return;
        VEC_SCALE( vec2, analy->vec_length_factor / vmag, vec );
        VEC_SCALE( vec3, (1.0 - analy->vec_length_factor) / vmax, vec );
        VEC_ADD( vec, vec2, vec3 );
    }

    /* Scale the vector. */
    VEC_ADDS( pts[0], -0.5 * vec_leng, vec, pt );
    VEC_ADDS( pts[1], vec_leng, vec, pts[0] );

    /* Draw a line segment -- no reflection is done. */
    glBegin( GL_LINES );
    glColor4fv( cols[0] );
    glVertex3fv( pts[0] );
    glColor4fv( cols[1] );
    glVertex3fv( pts[1] );
    glEnd();
}


/************************************************************
 * TAG( draw_reg_carpet )
 *
 * Draw a vector carpet on a jittered regular grid through the
 * volume of the model.
 */
static void
draw_reg_carpet( Analysis *analy )
{
    Hex_geom *bricks;
    Nodal *nodes;
    Bool_type front[6];
    float *vec_result[3];
    float vec_leng, vmin, vmax;
    float r, s, t, h[8];
    float vec[3], pt[3];
    int *reg_dim;
    int si, ei, incri, sj, ej, incrj, sk, ek, incrk;
    int ii, jj, kk, el, idx;
    int i, j, k;

    /* See if there are some to draw. */
    if ( analy->reg_dim[0] == 0 )
        return;

    bricks = analy->geom_p->bricks;
    nodes = analy->state_p->nodes;
    reg_dim = analy->reg_dim;

    /* Load the three results into vector result array. */
    load_vec_result( analy, vec_result, &vmin, &vmax );

    /* Make the max vector length about VEC_3D_LENGTH pixels long. */
    vec_leng = analy->vec_scale * VEC_3D_LENGTH /
               (v_win->vp_height * v_win->bbox_scale * vmax);

    /* Determine traversal directions. */
    find_front_faces( front, analy );

    if ( front[0] )
    {
        si = reg_dim[0] - 1;
        ei = -1;
        incri = -1;
    }
    else
    {
        si = 0;
        ei = reg_dim[0];
        incri = 1;
    }
    if ( front[2] )
    {
        sj = reg_dim[1] - 1;
        ej = -1;
        incrj = -1;
    }
    else
    {
        sj = 0;
        ej = reg_dim[1];
        incrj = 1;
    }
    if ( front[4] )
    {
        sk = reg_dim[2] - 1;
        ek = -1;
        incrk = -1;
    }
    else
    {
        sk = 0;
        ek = reg_dim[2];
        incrk = 1;
    }

    antialias_lines( TRUE, analy->z_buffer_lines );

    /* Loop through grid points from back to front. */
    for ( ii = si; ii != ei; ii += incri )
        for ( jj = sj; jj != ej; jj += incrj )
            for ( kk = sk; kk != ek; kk += incrk )
            {
                idx = ii*reg_dim[1]*reg_dim[2] + jj*reg_dim[2] + kk;

                /* Calculate the shape functions. */
                el = analy->reg_carpet_elem[idx];

                r = analy->reg_carpet_coords[0][idx];
                s = analy->reg_carpet_coords[1][idx];
                t = analy->reg_carpet_coords[2][idx];
                shape_fns_hex( r, s, t, h );

                /* Get the physical coordinates of the vector point. */
                VEC_SET( pt, 0.0, 0.0, 0.0 );
                for ( j = 0; j < 8; j++ )
                    for ( k = 0; k < 3; k++ )
                        pt[k] += h[j]*nodes->xyz[k][ bricks->nodes[j][el] ];

                /* Interpolate the vector quantity to the vector point. */
                VEC_SET( vec, 0.0, 0.0, 0.0 );
                for ( j = 0; j < 8; j++ )
                    for ( k = 0; k < 3; k++ )
                        vec[k] += h[j]*vec_result[k][ bricks->nodes[j][el] ];

                /* Draw the vector glyph. */
                draw_carpet_vector( pt, vec, vec_leng, vmin, vmax, analy );
            }

    antialias_lines( FALSE, 0 );

    /* Free temporary arrays. */
    for ( i = 0; i < 3; i++ )
        free( vec_result[i] );
}


/************************************************************
 * TAG( draw_vol_carpet )
 *
 * Draw carpet points contained in volume elements.
 */
static void
draw_vol_carpet( Analysis *analy )
{
    Hex_geom *bricks;
    Nodal *nodes;
    float *vec_result[3];
    float vec_leng, vmin, vmax;
    float r, s, t, h[8];
    float vec[3], pt[3];
    int *index;
    int el, i, j, k;

    /* See if there are some to draw. */
    if ( analy->vol_carpet_cnt == 0 )
        return;

    bricks = analy->geom_p->bricks;
    nodes = analy->state_p->nodes;

    /* Load the three results into vector result array. */
    load_vec_result( analy, vec_result, &vmin, &vmax );

    /* Avoid divide-by-zero, nothing to show in this case. */
    if ( vmax == 0.0 )
        return;

    /* Make the max vector length about VEC_3D_LENGTH pixels long. */
    vec_leng = analy->vec_scale * VEC_3D_LENGTH /
               (v_win->vp_height * v_win->bbox_scale * vmax);

    antialias_lines( TRUE, analy->z_buffer_lines );

    /* Sort the vector points. */
    index = NEW_N( int, analy->vol_carpet_cnt, "Carpet tmp" );
    if ( TRUE )
    {
        sort_carpet_points( analy, analy->vol_carpet_cnt,
                            analy->vol_carpet_coords,
                            analy->vol_carpet_elem, index, TRUE );
    }
    else
    {
        for ( i = 0; i < analy->vol_carpet_cnt; i++ )
            index[i] = i;
    }

    for ( i = analy->vol_carpet_cnt - 1; i >= 0; i-- )
    {
        /* Calculate the shape functions. */
        el = analy->vol_carpet_elem[index[i]];

        r = analy->vol_carpet_coords[0][index[i]];
        s = analy->vol_carpet_coords[1][index[i]];
        t = analy->vol_carpet_coords[2][index[i]];
        shape_fns_hex( r, s, t, h );

        /* Get the physical coordinates of the vector point. */
        VEC_SET( pt, 0.0, 0.0, 0.0 );
        for ( j = 0; j < 8; j++ )
            for ( k = 0; k < 3; k++ )
                pt[k] += h[j]*nodes->xyz[k][ bricks->nodes[j][el] ];

        /* Interpolate the vector quantity to the vector point. */
        VEC_SET( vec, 0.0, 0.0, 0.0 );
        for ( j = 0; j < 8; j++ )
            for ( k = 0; k < 3; k++ )
                vec[k] += h[j]*vec_result[k][ bricks->nodes[j][el] ];

        /* Draw the vector glyph. */
        draw_carpet_vector( pt, vec, vec_leng, vmin, vmax, analy );
    }

    antialias_lines( FALSE, 0 );

    /* Free temporary arrays. */
    for ( i = 0; i < 3; i++ )
        free( vec_result[i] );
    free( index );
}


/************************************************************
 * TAG( draw_shell_carpet )
 *
 * Draw carpet points contained in shell elements.
 */
static void
draw_shell_carpet( Analysis *analy )
{
    Shell_geom *shells;
    Nodal *nodes;
    float *vec_result[3];
    float vec_leng, vmin, vmax;
    float r, s, h[4];
    float vec[3], pt[3];
    int *index;
    int el, i, j, k;

    /* See if there are some to draw. */
    if ( analy->shell_carpet_cnt == 0 )
        return;

    shells = analy->geom_p->shells;
    nodes = analy->state_p->nodes;

    /* Load the three results into vector result array. */
    load_vec_result( analy, vec_result, &vmin, &vmax );

    /* Avoid divide-by-zero, nothing to show in this case. */
    if ( vmax == 0.0 )
        return;

    /* Make the max vector length about VEC_3D_LENGTH pixels long. */
    vec_leng = analy->vec_scale * VEC_3D_LENGTH /
               (v_win->vp_height * v_win->bbox_scale * vmax);

    antialias_lines( TRUE, analy->z_buffer_lines );

    /* Sort the vector points. */
    index = NEW_N( int, analy->shell_carpet_cnt, "Carpet tmp" );
    if ( analy->dimension == 3 )
    {
        sort_carpet_points( analy, analy->shell_carpet_cnt,
                            analy->shell_carpet_coords,
                            analy->shell_carpet_elem, index, FALSE );
    }
    else
    {
        for ( i = 0; i < analy->shell_carpet_cnt; i++ )
            index[i] = i;
    }

    for ( i = analy->shell_carpet_cnt - 1; i >= 0; i-- )
    {
        /* Calculate the shape functions. */
        el = analy->shell_carpet_elem[index[i]];
        r = analy->shell_carpet_coords[0][index[i]];
        s = analy->shell_carpet_coords[1][index[i]];
        shape_fns_quad( r, s, h );

        /* Get the physical coordinates of the vector point. */
        VEC_SET( pt, 0.0, 0.0, 0.0 );
        for ( j = 0; j < 4; j++ )
            for ( k = 0; k < 3; k++ )
                pt[k] += h[j]*nodes->xyz[k][ shells->nodes[j][el] ];

        /* Interpolate the vector quantity to the vector point. */
        VEC_SET( vec, 0.0, 0.0, 0.0 );
        for ( j = 0; j < 4; j++ )
            for ( k = 0; k < 3; k++ )
                vec[k] += h[j]*vec_result[k][ shells->nodes[j][el] ];

        /* Draw the vector glyph. */
        draw_carpet_vector( pt, vec, vec_leng, vmin, vmax, analy );
    }

    antialias_lines( FALSE, 0 );

    /* Free temporary arrays. */
    for ( i = 0; i < 3; i++ )
        free( vec_result[i] );
    free( index );
}


/************************************************************
 * TAG( draw_ref_carpet )
 *
 * Draw a vector carpet on the reference surface.
 */
static void
draw_ref_carpet( Analysis *analy )
{
    Ref_poly *poly;
    Nodal *nodes;
    float *vec_result[3];
    float vec_leng, vmin, vmax;
    float pr, ps, h[4];
    float verts[4][3], vec[3], pt[3];
    float vecs_per_area, vec_num, remain;
    float scale;
    int vec_cnt, seed;
    int i, j, k;

    /* See if there are some to draw. */
    if ( analy->ref_polys == NULL )
        return;

    nodes = analy->state_p->nodes;

    /* Load the three results into vector result array. */
    load_vec_result( analy, vec_result, &vmin, &vmax );

    /* Make the max vector length about VEC_3D_LENGTH pixels long. */
    vec_leng = analy->vec_scale * VEC_3D_LENGTH /
               (v_win->vp_height * v_win->bbox_scale * vmax);
    vecs_per_area = 1.0 / ( analy->vec_cell_size * analy->vec_cell_size );

    antialias_lines( TRUE, analy->z_buffer_lines );

    /* Seed the random function so jitter is consistent across passes. */
    seed = 50731;
    srand48( seed );

    for ( poly = analy->ref_polys; poly != NULL; poly = poly->next )
    {
        for ( i = 0; i < 4; i++ )
            for ( j = 0; j < 3; j++ )
                verts[i][j] = nodes->xyz[j][poly->nodes[i]];

        vec_num = area_of_quad( verts ) * vecs_per_area;
        vec_cnt = (int) vec_num;
        remain = vec_num - vec_cnt;

        if ( drand48() < remain )
            ++vec_cnt;

        /* Draw the vectors for this face. */
        for ( i = 0; i < vec_cnt; i++ )
        {
            pr = 2.0*drand48() - 1.0;
            ps = 2.0*drand48() - 1.0;

            /* Map these to physical coordinates, using shape functions. */
            shape_fns_quad( pr, ps, h );

            for ( k = 0; k < 3; k++ )
                pt[k] = h[0]*verts[0][k] + h[1]*verts[1][k] +
                        h[2]*verts[2][k] + h[3]*verts[3][k];

            /* Interpolate the vector quantities to the display point. */
            VEC_SET( vec, 0.0, 0.0, 0.0 );
            for ( j = 0; j < 4; j++ )
                for ( k = 0; k < 3; k++ )
                    vec[k] += h[j]*vec_result[k][poly->nodes[j]];

            /* Draw the vector glyph. */
            draw_carpet_vector( pt, vec, vec_leng, vmin, vmax, analy );
        }
    }

    antialias_lines( FALSE, 0 );

    /* Free temporary array. */
    for ( i = 0; i < 3; i++ )
        free( vec_result[i] );
}


#endif


/************************************************************
 * TAG( draw_traces )
 *
 * Draw all particle trace paths.
 */
static void
draw_traces( Analysis *analy )
{

    /* Draw existing traces. */
    if ( analy->trace_pts != NULL )
        draw_trace_list( analy->trace_pts, analy );

    /* Draw traces currently under construction. */
    if ( analy->new_trace_pts != NULL )
        draw_trace_list( analy->new_trace_pts, analy );
}


/************************************************************
 * TAG( draw_trace_list )
 *
 * Draw particle trace paths in a list.
 */
static void
draw_trace_list( Trace_pt_obj *ptlist, Analysis *analy )
{
    Trace_pt_obj *pt;
    float time;
    float rgb[3];
    int i;
    int skip, limit, new_skip, new_limit;
    Bool_type init_subtrace;
    Mesh_data *p_mesh;

    p_mesh = MESH_P( analy );

    limit = analy->ptrace_limit;

    antialias_lines( TRUE, analy->z_buffer_lines );
    glLineWidth( (GLfloat) analy->trace_width );

    time = analy->state_p->time;

    for ( pt = ptlist; pt != NULL; pt = pt->next )
    {
        if ( pt->color[0] < 0 )
            glColor3fv( v_win->foregrnd_color );
        else
        {
            VEC_COPY( rgb, pt->color );
            glColor3fv( rgb );
        }

        /* For static field traces, only draw if current time is trace time. */
        if ( pt->time[0] == pt->time[1] )
        {
            if ( pt->time[0] != time )
                continue;
            else
                i = pt->cnt;
        }
        else
            /* Only draw the trace up to the current time. */
            for ( i = 0; i < pt->cnt; i++ )
                if ( pt->time[i] > time )
                    break;

        skip = 0;
        if ( limit > 0 )
            if ( i > limit )
            {
                skip = i - limit;
                i = limit;
            }

        if ( analy->trace_disable != NULL )
        {
            init_subtrace = TRUE;
            while ( find_next_subtrace( init_subtrace, skip, i, pt, analy,
                                        &new_skip, &new_limit ) )
            {
                init_subtrace = FALSE;
                draw_line( new_limit, pt->pts + 3 * new_skip, -1, p_mesh, analy,
                           FALSE, NULL );
            }
        }
        else
            draw_line( i, pt->pts + 3 * skip, -1, p_mesh, analy, FALSE, NULL );
    }

    glLineWidth( (GLfloat) 1.0 );
    antialias_lines( FALSE, 0 );
}


/************************************************************
 * TAG( draw_sphere )
 *
 * Draw a polygonized sphere.  The center of the sphere and
 * the radius are specified by parameters.  There are better
 * ways to do this -- e.g. call the GL utility code.
 */
static void
draw_sphere( float ctr[3], float radius, int res_factor)
{
    float latincr, longincr, latangle, longangle, length;
    float vert[3], norm[3];
    int latres, longres, i, j;

    latres  = 8*res_factor;
    longres = 5*res_factor;

    latincr = 2.0*PI / latres;
    longincr = PI / longres;

    /* Draw a unit sphere. */
    glBegin( GL_QUADS );
    for ( i = 0, latangle = 0.0; i < latres; i++, latangle += latincr )
    {
        for ( j = 0, longangle = -PI/2.0;
                j < longres;
                j++, longangle += longincr )
        {
            length = cos( (double)longangle );
            norm[0] = cos( (double)latangle ) * length;
            norm[1] = sin( (double)latangle ) * length;
            norm[2] = sin( (double)longangle );
            VEC_ADDS( vert, radius, norm, ctr );
            glNormal3fv( norm );
            glVertex3fv( vert );
            norm[0] = cos( (double)(latangle+latincr) ) * length;
            norm[1] = sin( (double)(latangle+latincr) ) * length;
            VEC_ADDS( vert, radius, norm, ctr );
            glNormal3fv( norm );
            glVertex3fv( vert );
            length = cos( (double)(longangle+longincr) );
            norm[0] = cos( (double)(latangle+latincr) ) * length;
            norm[1] = sin( (double)(latangle+latincr) ) * length;
            norm[2] = sin( (double)(longangle+longincr) );
            VEC_ADDS( vert, radius, norm, ctr );
            glNormal3fv( norm );
            glVertex3fv( vert );
            norm[0] = cos( (double)latangle ) * length;
            norm[1] = sin( (double)latangle ) * length;
            VEC_ADDS( vert, radius, norm, ctr );
            glNormal3fv( norm );
            glVertex3fv( vert );
        }
    }
    glEnd();
}

/************************************************************
 * TAG( draw_sphere_GL )
 *
 * This function will draw a solid sphere using the
 * GL function.
 */
static void
draw_sphere_GL( float ctr[3], float radius, int res_factor)
{

    int num_slices=1, num_stacks=1;

    /* Variables related to drawing spheres */
    GLUquadricObj *sphereObj;
    sphereObj = gluNewQuadric();

    num_slices*=res_factor;
    num_stacks*=res_factor;

    glPushMatrix();
    glTranslatef( ctr[0], ctr[1], ctr[2] );
    gluSphere(sphereObj, radius, num_slices, num_stacks);
    glPopMatrix();
    glEnd();
}


/************************************************************
 * TAG( draw_poly )
 *
 * Draw a polygon, reflecting it as needed.  The routine expects
 * either a quadrilateral or a triangle.  The matl flag is the
 * material of the polygon.  A value of -1 indicates
 * that no material is associated with the polygon.
 *
 * NOTE: the calling routine needs to set up and shut down the
 * stencil buffer for polygon edge drawing (in hidden line mode.)
 */
static void
draw_poly( int cnt, float pts[4][3], float norm[4][3], float cols[4][4],
           float vals[4], int matl, Mesh_data *p_mesh, Analysis *analy,
           Bool_type hidden_poly )
{
    Refl_plane_obj *plane;
    Render_poly_obj *poly, *rpoly, *origpoly, *reflpoly;
    float rpts[4][3];
    float rnorm[4][3];
    float rcols[4][4];
    float rvals[4];
    int i, ii, j;

    if ( p_mesh->translate_material )
    {
        if ( matl >= 0 )
            for ( i = 0; i < cnt; i++ )
                for ( j = 0; j < 3; j++ )
                    pts[i][j] += p_mesh->mtl_trans[j][matl];
    }

    /* Draw the initial polygon. */
    if ( analy->mesh_view_mode == RENDER_HIDDEN )
        draw_edged_poly( cnt, pts, norm, cols, vals, matl, p_mesh, analy, hidden_poly );
    else if ( analy->mesh_view_mode == RENDER_WIREFRAME )
        draw_edged_wireframe_poly( cnt, pts, norm, cols, vals, matl, p_mesh, analy );
    else if ( analy->mesh_view_mode == RENDER_WIREFRAMETRANS )
        draw_edged_poly( cnt, pts, norm, cols, vals, matl, p_mesh, analy, hidden_poly );
    else
        draw_plain_poly( cnt, pts, norm, cols, vals, matl, p_mesh, analy, hidden_poly );

    if ( !analy->reflect || analy->refl_planes == NULL )
        return;

    if ( analy->refl_planes->next == NULL )
    {
        /* Just one reflection plane; blast it out.
         * The polygon's node order is reversed by reflecting
         * the individual points, so we have to correct the order.
         */
        plane = analy->refl_planes;
        for ( i = 0; i < cnt; i++ )
        {
            ii = cnt - 1 - i;
            point_transform( rpts[i], pts[ii], &plane->pt_transf );
            point_transform( rnorm[i], norm[ii], &plane->norm_transf );
            hvec_copy( rcols[i], cols[ii] );
            rvals[i] = vals[ii];
        }

        if ( analy->mesh_view_mode == RENDER_HIDDEN )
            draw_edged_poly( cnt, rpts, rnorm, rcols, rvals, matl, p_mesh,
                             analy, hidden_poly );
        else if ( analy->mesh_view_mode == RENDER_WIREFRAME )
            draw_edged_wireframe_poly( cnt, rpts, rnorm, rcols, rvals, matl, p_mesh,
                                       analy );
        else if ( analy->mesh_view_mode == RENDER_WIREFRAMETRANS )
            draw_edged_poly( cnt, pts, norm, cols, vals, matl, p_mesh, analy, hidden_poly );
        else
            draw_plain_poly( cnt, rpts, rnorm, rcols, rvals, matl, p_mesh,
                             analy, hidden_poly );
    }
    else
    {
        /* More than one reflection plane, do it the hard way. */

        origpoly = NEW( Render_poly_obj, "Reflection polygon" );
        origpoly->cnt = cnt;
        for ( i = 0; i < cnt; i++ )
        {
            VEC_COPY( origpoly->pts[i], pts[i] );
            VEC_COPY( origpoly->norm[i], norm[i] );
            hvec_copy( origpoly->cols[i], cols[i] );
            origpoly->vals[i] = vals[i];
        }

        reflpoly = NULL;

        for ( plane = analy->refl_planes;
                plane != NULL;
                plane = plane->next )
        {
            for ( poly = origpoly; poly != NULL; poly = poly->next )
            {
                rpoly = NEW( Render_poly_obj, "Reflection polygon" );

                for ( i = 0; i < cnt; i++ )
                {
                    ii = cnt - 1 - i;
                    point_transform( rpoly->pts[i], poly->pts[ii],
                                     &plane->pt_transf );
                    point_transform( rpoly->norm[i], poly->norm[ii],
                                     &plane->norm_transf );
                    hvec_copy( rpoly->cols[i], poly->cols[ii] );
                    rpoly->vals[i] = poly->vals[ii];
                }

                if ( analy->mesh_view_mode == RENDER_HIDDEN )
                    draw_edged_poly( cnt, rpoly->pts, rpoly->norm, rpoly->cols,
                                     rpoly->vals, matl, p_mesh, analy, FALSE );
                else if ( analy->mesh_view_mode == RENDER_WIREFRAME )
                    draw_edged_poly( cnt, rpoly->pts, rpoly->norm, rpoly->cols,
                                     rpoly->vals, matl, p_mesh, analy, FALSE );
                else
                    draw_plain_poly( cnt, rpoly->pts, rpoly->norm, rpoly->cols,
                                     rpoly->vals, matl, p_mesh, analy, NULL );

                INSERT( rpoly, reflpoly );
            }

            if ( analy->refl_orig_only )
                free( reflpoly );
            else
                APPEND( reflpoly, origpoly );
            reflpoly = NULL;
        }

        DELETE_LIST( origpoly );
    }
}


/************************************************************
 * TAG( draw_edged_poly )
 *
 * Draw a polygon with borders, for hidden line drawings.
 * The routine expects either a quadrilateral or a triangle.
 *
 * NOTE: stenciling must be set up by the calling routine.
 */
static void
draw_edged_poly( int cnt, float pts[4][3], float norm[4][3], float cols[4][4],
                 float vals[4], int matl, Mesh_data *p_mesh, Analysis *analy,
                 Bool_type hidden_poly )
{
    int i, j;

    /* The hidden poly feature is not used for now - when activated we will be able to
     * render wireframe elements solid elements together.
     */

    if ( hidden_poly && analy->mesh_view_mode != RENDER_WIREFRAMETRANS )
        return;

    /* Fix degenerate faces */

    for (i=0; i<4; i++)
        for (j=0; j<4; j++)
            if (j != i)
                if (pts[i][0]==pts[j][0] && pts[i][1]==pts[j][1] && pts[i][2]==pts[j][2])
                {
                    pts[i][0]+= 0.000001;
                    pts[i][1]+= 0.000001;
                    pts[i][2]+= 0.000001;
                }

    /* Outline the polygon and write into the stencil buffer. */
    if ( hidden_poly )
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    glBegin( GL_LINE_LOOP );
    for ( i = 0; i < cnt; i++ )
        glVertex3fv( pts[i] );
    glEnd();

    glStencilFunc( GL_EQUAL, 0, 1 );
    glStencilOp( GL_KEEP, GL_KEEP, GL_KEEP );

    /* Fill the polygon. */
    if ( !hidden_poly )
    {
        if ( analy->interp_mode == GOOD_INTERP )
        {
            if(!grayel)
            {
                /* Scan convert the polygon by hand. */
                scan_poly( cnt, pts, norm, vals, matl, p_mesh, analy );
            } else
            {
                glBegin( GL_POLYGON );
                for ( i = 0; i < cnt; i++ )
                {
                    glColor3fv( cols[i] );
                    glNormal3fv( norm[i] );
                    glVertex3fv( pts[i] );
                }
                glEnd();
            }
        }
        else
        {
            glBegin( GL_POLYGON );
            for ( i = 0; i < cnt; i++ )
            {
                glColor3fv( cols[i] );
                glNormal3fv( norm[i] );
                glVertex3fv( pts[i] );
            }
            glEnd();
        }

        /* Clear the stencil buffer. */
        glStencilFunc( GL_ALWAYS, 0, 1 );
        glStencilOp( GL_INVERT, GL_INVERT, GL_INVERT );
    }

    /* Draw Outline #2 */

    glColor3fv( v_win->mesh_color  );
    glStencilFunc( GL_ALWAYS, 0, 1 );
    glStencilOp( GL_INVERT, GL_INVERT, GL_INVERT );

    glColor3fv( v_win->mesh_color  );
    glBegin( GL_LINE_LOOP );
    for ( i = 0; i < cnt; i++ )
        glVertex3fv( pts[i] );
    glEnd();


    if (hidden_poly )
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}


/************************************************************
 * TAG( draw_edged_wireframe_poly )
 *
 * Draw a polygon with borders, for hidden line drawings.
 * The routine expects either a quadrilateral or a triangle.
 *
 * NOTE: stenciling must be set up by the calling routine.
 */
static void
draw_edged_wireframe_poly( int cnt, float pts[4][3], float norm[4][3], float cols[4][4],
                           float vals[4], int matl, Mesh_data *p_mesh, Analysis *analy )
{
    int i, j;

    for (i=0; i<4; i++)
        for (j=0; j<4; j++)
            if (j != i)
                if (pts[i][0]==pts[j][0] && pts[i][1]==pts[j][1] && pts[i][2]==pts[j][2])
                {
                    pts[i][0]+= 0.000000001;
                    pts[i][1]+= 0.000000001;
                    pts[i][2]+= 0.000000001;
                }

    /* Outline the polygon and write into the stencil buffer. */

    glColor3fv( v_win->mesh_color );

    /* Outline #1 */
    glBegin( GL_LINE_LOOP );
    for ( i = 0; i < cnt; i++ )
        glVertex3fv( pts[i] );
    glEnd();


    /* Draw Polygons */

    /* Fill the polygon. */
    glStencilFunc( GL_EQUAL, 0, 1 );
    glStencilOp( GL_KEEP, GL_KEEP, GL_KEEP );

    if ( analy->interp_mode == GOOD_INTERP )
    {
        if(!grayel)
        { 
            /* Scan convert the polygon by hand. */
            scan_poly( cnt, pts, norm, vals, matl, p_mesh, analy );
        } else
        {
            glBegin( GL_POLYGON );
            for ( i = 0; i < cnt; i++ )
            {
                glColor3fv( cols[i] );
                glColor3fv( v_win->backgrnd_color  );
                glNormal3fv( norm[i] );
                glVertex3fv( pts[i] );
            }
            glEnd();
       }
    }
    else
    {
        glBegin( GL_POLYGON );
        for ( i = 0; i < cnt; i++ )
        {
            glColor3fv( cols[i] );
            glColor3fv( v_win->backgrnd_color  );
            glNormal3fv( norm[i] );
            glVertex3fv( pts[i] );
        }
        glEnd();
    }

    /* Draw Outline #2 */
    glColor3fv( v_win->mesh_color  );
    glStencilFunc( GL_ALWAYS, 0, 1 );
    glStencilOp( GL_INVERT, GL_INVERT, GL_INVERT );

    glBegin( GL_LINE_LOOP );
    for ( i = 0; i < cnt; i++ )
        glVertex3fv( pts[i] );
    glEnd();

    /* Clear the stencil buffer. */
    glStencilFunc( GL_ALWAYS, 0, 1 );
    glStencilOp( GL_INVERT, GL_INVERT, GL_INVERT );
    glColor3fv( v_win->mesh_color  );
}


/************************************************************
 * TAG( draw_plain_poly )
 *
 * Draw a polygon with no borders.  The routine expects
 * either a quadrilateral or a triangle.
 */
static void
draw_plain_poly( int cnt, float pts[4][3], float norm[4][3], float cols[4][4],
                 float vals[4], int matl, Mesh_data *p_mesh, Analysis *analy,
                 Bool_type hidden_poly )
{
    int i;

    if ( hidden_poly )
        return;

    if ( analy->interp_mode == GOOD_INTERP )
    {
        if(!grayel)
        {
            /* Scan convert the polygon by hand. */
            scan_poly( cnt, pts, norm, vals, matl, p_mesh, analy );
        } else
        {    
            glBegin( GL_POLYGON );
            for ( i = 0; i < cnt; i++ )
            {
                glColor3fv( cols[i] );
                glNormal3fv( norm[i] );
                glVertex3fv( pts[i] );
            }
            glEnd();
        
        }
    }
    else
    {
        glBegin( GL_POLYGON );
        for ( i = 0; i < cnt; i++ )
        {
            glColor3fv( cols[i] );
            glNormal3fv( norm[i] );
            glVertex3fv( pts[i] );
        }
        glEnd();
    }
}


/************************************************************
 * TAG( length_2d )
 *
 * Return the length between two 2D points.  Helper function
 * for scan_poly.
 */
static float
length_2d( float pt1[2], float pt2[2] )
{
    double dx, dy;

    dx = (double) pt2[0] - pt1[0];
    dy = (double) pt2[1] - pt1[1];
    return sqrt( (double)( dx*dx + dy*dy ) );
}


/************************************************************
 * TAG( scan_poly )
 *
 * Draw a polygon by subdividing it into pixel-sized fragments.
 * This routine allows us to implement correct (but slow) linear
 * interpolation in the polygon interior.  Handles 3D and 2D
 * polygons.
 */
static void
scan_poly( int cnt, float pts[4][3], float norm[4][3], float vals[4], int matl,
           Mesh_data *p_mesh, Analysis *analy )
{
    float pt[3], proj_pts[4][2];
    float d1, d2;
    float r, s, h[4], val;
    float *dpts, *dcol, *dnorm;
    float rmin, rmax, threshold;
    int dist, ni, nj;
    int i, ii, j, k;
 
    /*
     * For triangular faces, duplicate the last vertex value
     * and then treat as a quad.  Is interpolation correct in
     * this case?
     */
    if ( cnt == 3 )
    {
        VEC_COPY( pts[3], pts[2] );
        VEC_COPY( norm[3], norm[2] );
        vals[3] = vals[2];
    }

    if ( vals[0] == vals[1] && vals[0] == vals[2] && vals[0] == vals[3] )
    {
        /*
         * If vertex values are all the same, no need to subdivide polygon.
         */
        ni = 2;
        nj = 2;
    }
    else
    {
        /* Transform the vertices and project to the screen. */
        for ( i = 0; i < 4; i++ )
        {
            /* Run the vertex through the view matrix. */
            point_transform( pt, pts[i], &cur_view_mat );

            /*
             * See glFrustrum() and glOrtho() in the OpenGL Reference Manual
             * to see where these factors come from.
             */
            if ( v_win->orthographic )
            {
                proj_pts[i][0] = proj_param_x * pt[0];
                proj_pts[i][1] = proj_param_y * pt[1];
            }
            else
            {
                proj_pts[i][0] = (double) proj_param_x * pt[0] / -pt[2];
                proj_pts[i][1] = (double) proj_param_y * pt[1] / -pt[2];
            }
        }

        /* Estimate the number of subdivisions needed in each direction. */
        d1 = length_2d( proj_pts[0], proj_pts[1] );
        d2 = length_2d( proj_pts[2], proj_pts[3] );
        /*
         * Force subdivision by incrementing distances.  This is really only
         * important for polygons with dimensions of approx. 2 pixels, but
         * will have little computational impact for larger polygons where
         * it's unnecessary.
         */
        /* dist = MAX( d1, d2 ); */
        dist = MAX( d1 + 1.5, d2 + 1.5 );
        ni = MAX( dist, 2 );
        d1 = length_2d( proj_pts[1], proj_pts[2] );
        d2 = length_2d( proj_pts[0], proj_pts[3] );
        dist = MAX( d1 + 1.5, d2 + 1.5 );
        nj = MAX( dist, 2 );
    }

    rmin = analy->result_mm[0];
    rmax = analy->result_mm[1];
    threshold = analy->zero_result;


    /* Dice the polygon and render the fragments. */
    dpts = NEW_N( float, ni*nj*3, "Tmp scan poly" );
    dcol = NEW_N( float, ni*nj*4, "Tmp scan poly" );
    dnorm = NEW_N( float, ni*nj*3, "Tmp scan poly" );

    for ( i = 0; i < ni; i++ )
        for ( j = 0; j < nj; j++ )
        {
            r = 2.0 * i / (ni - 1.0) - 1.0;
            s = 2.0 * j / (nj - 1.0) - 1.0;
            shape_fns_quad( r, s, h );

            for ( k = 0; k < 3; k++ )
                dpts[i*nj*3 + j*3 + k] = (double) h[0]*pts[0][k]
                                         + (double) h[1]*pts[1][k]
                                         + (double) h[2]*pts[2][k]
                                         + (double) h[3]*pts[3][k];
            for ( k = 0; k < 3; k++ )
                dnorm[i*nj*3 + j*3 + k] = (double) h[0]*norm[0][k]
                                          + (double) h[1]*norm[1][k]
                                          + (double) h[2]*norm[2][k]
                                          + (double) h[3]*norm[3][k];

            val = (double) h[0]*vals[0] + (double) h[1]*vals[1]
                  + (double) h[2]*vals[2] + (double) h[3]*vals[3];

            color_lookup( &dcol[i*nj*4 + j*4], val, rmin, rmax,
                          threshold, matl, analy->logscale,
                          analy->material_greyscale );
        }

    if ( analy->dimension == 3 )
    {
        for ( i = 0, ii = 1; i < ni - 1; i++, ii++ )
        {
            glBegin( GL_QUAD_STRIP );
            for ( j = 0; j < nj; j++ )
            {
                glColor3fv( &dcol[i*nj*4 + j*4] );
                glNormal3fv( &dnorm[i*nj*3 + j*3] );
                glVertex3fv( &dpts[i*nj*3 + j*3] );
                if(!grayel)
                {
                    glColor3fv( &dcol[ii*nj*4 + j*4] );
                }
                glNormal3fv( &dnorm[ii*nj*3 + j*3] );
                glVertex3fv( &dpts[ii*nj*3 + j*3] );
            }
            glEnd();
        }
    }
    else
    {
        /* We make this a separate loop only for efficiency. */
        for ( i = 0, ii = 1; i < ni - 1; i++, ii++ )
        {
            glBegin( GL_QUAD_STRIP );
            for ( j = 0; j < nj; j++ )
            {
                glColor3fv( &dcol[i*nj*4 + j*4] );
                glVertex2fv( &dpts[i*nj*3 + j*3] );
                glColor3fv( &dcol[ii*nj*4 + j*4] );
                glVertex2fv( &dpts[ii*nj*3 + j*3] );
            }
            glEnd();
        }
    }

    free( dpts );
    free( dcol );
    free( dnorm );
}


/************************************************************
 * TAG( draw_poly_2d )
 *
 * Draw a 2D polygon, reflecting it as needed.  The routine expects
 * either a quadrilateral or a triangle.  The matl flag is the
 * material of the polygon.  A value of -1 indicates that no
 * material is associated with the polygon.
 */
static void
draw_poly_2d( int cnt, float pts[4][3], float cols[4][4], float vals[4],
              int matl, Mesh_data *p_mesh, Analysis *analy )
{
    Bool_type good_interp;
    Refl_plane_obj *plane;
    Render_poly_obj *poly, *rpoly, *origpoly, *reflpoly;
    float rpts[4][3];
    float rcols[4][4];
    float rvals[4];
    float norm[4][3];
    int i, ii, j;

    good_interp = ( analy->interp_mode == GOOD_INTERP && matl >= 0 );

    if ( p_mesh->translate_material )
    {
        if ( matl >= 0 )
            for ( i = 0; i < cnt; i++ )
                for ( j = 0; j < 2; j++ )
                    pts[i][j] += p_mesh->mtl_trans[j][matl];
    }

    /* Draw the initial polygon. */
    if ( good_interp )
    {
        /* This is just a dummy argument in this case. */
        for ( i = 0; i < cnt; i++ )
        {
            VEC_SET( norm[i], 0.0, 0.0, 1.0 );
        }
        if(!grayel)
        {
            scan_poly( cnt, pts, norm, vals, matl, p_mesh, analy );
        } else
        {
            glBegin( GL_POLYGON );
            for ( i = 0; i < cnt; i++ )
            {
                glColor3fv( cols[i] );
                glVertex2fv( pts[i] );
            }
            glEnd();
        }
    }
    else
    {
        glBegin( GL_POLYGON );
        for ( i = 0; i < cnt; i++ )
        {
            glColor3fv( cols[i] );
            glVertex2fv( pts[i] );
        }
        glEnd();
    }

    if ( !analy->reflect || analy->refl_planes == NULL )
        return;

    if ( analy->refl_planes->next == NULL )
    {
        /* Just one reflection plane; blast it out.
         * The polygon's node order is reversed by reflecting
         * the individual points, so we have to correct the order.
         */
        plane = analy->refl_planes;
        for ( i = 0; i < cnt; i++ )
        {
            ii = cnt - 1 - i;
            point_transform( rpts[i], pts[ii], &plane->pt_transf );
            hvec_copy( rcols[i], cols[ii] );
            rvals[i] = vals[ii];
        }

        if ( good_interp )
            scan_poly( cnt, rpts, norm, rvals, matl, p_mesh, analy );
        else
        {
            glBegin( GL_POLYGON );
            for ( i = 0; i < cnt; i++ )
            {
                glColor3fv( rcols[i] );
                glVertex2fv( rpts[i] );
            }
            glEnd();
        }
    }
    else
    {
        /* More than one reflection plane, do it the hard way. */

        origpoly = NEW( Render_poly_obj, "Reflection polygon" );
        origpoly->cnt = cnt;
        for ( i = 0; i < cnt; i++ )
        {
            VEC_COPY( origpoly->pts[i], pts[i] );
            hvec_copy( origpoly->cols[i], cols[i] );
            origpoly->vals[i] = vals[i];
        }

        reflpoly = NULL;

        for ( plane = analy->refl_planes;
                plane != NULL;
                plane = plane->next )
        {
            for ( poly = origpoly; poly != NULL; poly = poly->next )
            {
                rpoly = NEW( Render_poly_obj, "Reflection polygon" );

                for ( i = 0; i < cnt; i++ )
                {
                    ii = cnt - 1 - i;
                    point_transform( rpoly->pts[i], poly->pts[ii],
                                     &plane->pt_transf );
                    hvec_copy( rpoly->cols[i], poly->cols[ii] );
                    rpoly->vals[i] = poly->vals[ii];
                }

                if ( good_interp )
                    scan_poly( cnt, rpoly->pts, norm, rpoly->vals,
                               matl, p_mesh, analy );
                else
                {
                    glBegin( GL_POLYGON );
                    for ( i = 0; i < cnt; i++ )
                    {
                        glColor3fv( rpoly->cols[i] );
                        glVertex2fv( rpoly->pts[i] );
                    }
                    glEnd();
                }

                INSERT( rpoly, reflpoly );
            }

            if ( analy->refl_orig_only )
                free( reflpoly );
            else
                APPEND( reflpoly, origpoly );
            reflpoly = NULL;
        }

        DELETE_LIST( origpoly );
    }
}


/************************************************************
 * TAG( draw_line )
 *
 * Draw a polyline, reflecting it as needed.  The "close" flag
 * should be TRUE for a closed polyline and FALSE for an open
 * polyline.  The matl flag is the material of the line.
 * A value of -1 indicates that no material is associated with
 * the line.
 */
static void
draw_line( int cnt, float *pts, int matl, Mesh_data *p_mesh, Analysis *analy,
           Bool_type close, Bool_type *draw_reflection_flag )
{
    GLenum line_mode;
    Refl_plane_obj *plane;
    float rpt[3];
    float *origpts[128], *reflpts[128];
    int origcnt, reflcnt;
    int i, j;
    float col[4];

    int plane_index=0;
    Bool_type draw_reflection_flag_temp[10], draw_line=TRUE;

    if(analy->mesh_view_mode == RENDER_WIREFRAMETRANS)
    {
        for(i = 0; i < 4; i++)
        {
           col[i] = 0.0;
        }
        glColor3fv(col);
    }

    /* Handles up to 6 cumulative reflection planes. */

    if ( p_mesh->translate_material )
    {
        if ( matl >= 0 )
            for ( i = 0; i < cnt; i++ )
                for ( j = 0; j < 3; j++ )
                    pts[i*3+j] += p_mesh->mtl_trans[j][matl];
    }

    if ( close )
        line_mode = GL_LINE_LOOP;
    else
        line_mode = GL_LINE_STRIP;

    if ( draw_reflection_flag )
        if ( draw_reflection_flag[0] )
            draw_line = TRUE;
        else
            draw_line = FALSE;

    /* Draw the original line. */
    if ( draw_line )
    {
        glBegin( line_mode );
        for ( i = 0; i < cnt; i++ )
            glVertex3fv( &pts[i*3] );
        glEnd();
    }

    if ( !analy->reflect || analy->refl_planes == NULL )
        return;

    if ( draw_reflection_flag == NULL)
        for ( i=0;
                i<10;
                i++ )
            draw_reflection_flag_temp[i] = TRUE;
    else
        for ( i=0;
                i<10;
                i++ )
            draw_reflection_flag_temp[i] = draw_reflection_flag[i];

    if ( analy->refl_planes->next == NULL )
    {
        /* Just one reflection plane; blast it out. */

        if (draw_reflection_flag_temp[0])
        {
            plane = analy->refl_planes;
            glBegin( line_mode );
            for ( i = 0; i < cnt; i++ )
            {
                point_transform( rpt, &pts[i*3], &plane->pt_transf );
                glVertex3fv( rpt );
            }
            glEnd();
        }
    }
    else if ( analy->refl_orig_only )
    {
        for ( plane = analy->refl_planes;
                plane != NULL;
                plane = plane->next )
        {
            if (draw_reflection_flag_temp[plane_index])
            {
                glBegin( line_mode );
                for ( i = 0; i < cnt; i++ )
                {
                    point_transform( rpt, &pts[i*3], &plane->pt_transf );
                    glVertex3fv( rpt );
                }
                glEnd();
            }
            plane_index++;
        }
    }
    else
    {
        origpts[0] = pts;
        origcnt = 1;
        reflcnt = 0;

        for ( plane = analy->refl_planes;
                plane != NULL;
                plane = plane->next )
        {
            if (draw_reflection_flag_temp[plane_index])
            {
                for ( j = 0; j < origcnt; j++ )
                {
                    reflpts[reflcnt] = NEW_N( float, cnt*3, "Reflect line" );
                    glBegin( line_mode );
                    for ( i = 0; i < cnt; i++ )
                    {
                        point_transform( &reflpts[reflcnt][i*3], &origpts[j][i*3],
                                         &plane->pt_transf );
                        glVertex3fv( &reflpts[reflcnt][i*3] );
                    }
                    glEnd();
                    reflcnt++;
                }

                for ( j = 0; j < reflcnt; j++ )
                {
                    if ( origcnt >= 128 )
                        popup_fatal( "Too many reflection planes!\n" );

                    origpts[origcnt] = reflpts[j];
                    origcnt++;
                }
                reflcnt = 0;
            }

            plane_index++;
        }

        for ( j = 1; j < origcnt; j++ )
        {
            if ( origpts[j]!=NULL )
            {
                free( origpts[j] );
                origpts[j] = NULL;
            }
        }
    }
}


/************************************************************
 * TAG( draw_3d_text )
 *
 * Draw 3D text at the position specified by pt.  This routine
 * orients the text toward the viewer and draws it with the
 * proper scaling and rotation.
 */
static void
draw_3d_text( float pt[3], char *text, Bool_type center_text )
{
    Transf_mat mat;
    float arr[16], tpt[3], spt[3];
    float zpos, cx, cy, text_height;
    int i, j;

    /* Run the point through the model transform matrix. */
    glGetFloatv( GL_MODELVIEW_MATRIX, arr );
    for ( i = 0; i < 4; i++ )
        for ( j = 0; j < 4; j++ )
            mat.mat[i][j] = arr[i*4 + j];
    point_transform( tpt, pt, &mat );

    /* Get drawing window and position. */
    get_foreground_window( &zpos, &cx, &cy );

    /* Project the point to the drawing window. */
    if ( v_win->orthographic )
    {
        spt[0] = tpt[0];
        spt[1] = tpt[1];
        spt[2] = zpos;
    }
    else
    {
        spt[0] = tpt[0] * zpos / tpt[2];
        spt[1] = tpt[1] * zpos / tpt[2];
        spt[2] = zpos;
    }

    /* Get the text size. */
    text_height = 14.0 * 2.0*cy / v_win->vp_height;

    /* Clear the model view matrix and draw the text. */
    glPushMatrix();
    glLoadIdentity();

    hmove( spt[0], spt[1], spt[2] );
    htextsize( text_height, text_height );
    if ( center_text )
        hcentertext( TRUE );
    antialias_lines( TRUE, TRUE );
    glLineWidth( 1.25 );
    hcharstr( text );
    glLineWidth( 1.0 );
    antialias_lines( FALSE, 0 );

    /* Restore the model view matrix. */
    glPopMatrix();
}


/************************************************************
 * TAG( get_foreground_window )
 *
 * Calculate a Z position for the foreground that won't get
 * clipped by the near/far planes.  Also returns the view
 * coordinates of the upper right corner of the window at
 * the Z position.
 */
void
get_foreground_window( float *zpos, float *cx, float *cy )
{
    float aspect, cp;

    /* Get a Z position that won't get clipped by near/far planes.
     * Translate all the text to that position.  Have to scale in
     * X and Y if perspective viewing is used.
     */
    *zpos = - (v_win->far + v_win->near)/2.0;

    aspect = v_win->aspect_correct * v_win->vp_width / (float)v_win->vp_height;

    /* Get the viewport dimensions at the drawing (zpos) plane. */
    if ( v_win->orthographic )
        cp = 1.0;
    else
        cp = - *zpos * TAN( DEG_TO_RAD( 0.5*v_win->cam_angle ) );

    if ( aspect >= 1.0 )
    {
        *cx = cp * aspect;
        *cy = cp;
    }
    else
    {
        *cx = cp;
        *cy = cp / aspect;
    }
}


/************************************************************
 * TAG( draw_foreground )
 *
 * Draw the foreground stuff in the mesh window (colormap,
 * information strings, world coordinate system, etc.)
 * stuff in the mesh window (colormap, information strings, etc.)
 */
static void
draw_foreground( Analysis *analy )
{
    Transf_mat mat, tmat;
    float pt[3], pto[3], pti[3], ptj[3], ptk[3], ptv[3];
    float (*ttmat)[3];
    float world_to_vp[2], vp_to_world[2];
    float zpos, cx, cy;
    float xpos=0.0, ypos=0.0, xp=0.0, yp=0.0, xsize, ysize;
    float scalax, scalay, scalaz;
    float text_height, b_width, b_height, cw, ch;
    float low, high;
    float leng, sub_leng;
    Bool_type extend_colormap, raw_minmax, show_dirvec, found_data;
    char str[90];
    static char *st_max = "State Maximum:";
    static char *st_min = "State Minimum:";
    static char *glob_max = "Global Maximum:";
    static char *glob_min = "Global Minimum:";
    char *maximum_label;
    char *minimum_label;
    int nstripes, i, frst, lst;
    static char *strain_label[] =
    {
        "infinitesimal", "Green-Lagrange", "Almansi", "Rate"
    };
    static char *ref_surf_label[] = { "middle", "inner", "outer" };
    static char *ref_frame_label[] = { "global", "local" };
    char *vec_x, *vec_y, *vec_z;
    char *zero = "0";
    float *el_mm;
    char **classes;
    int *sclasses;
    int *el_id;
    Result_spec *p_rs;
    int fracsz;
    float scale_y;
    double value;
    float scale_minimum, scale_maximum, distance;
    int qty_of_intervals, scale_error_status;
    float rmin_offset, rmax_offset;
    int ntrips;
    float low_text_bound, high_text_bound;
    float comparison_tolerance;
    float scl_max;
    Mesh_data *p_mesh;
    float map_text_offset;
    int dim;

    int ylines=2;

    /* Error Indicator (EI) variables */
    Bool_type ei_labels = FALSE;

    /* Class labeling variables */
    int class_label;
    MO_class_data *class_ptr;

    if ( analy->ei_result && analy->result_active )
        ei_labels = TRUE;

    /* Cap for colorscale interpolation. */
    scl_max = SCL_MAX;

    /* Foreground overwrites everything else in the scene. */
    glDepthFunc( GL_ALWAYS );

    /* Get drawing window and position. */
    get_foreground_window( &zpos, &cx, &cy );

    world_to_vp[0] = v_win->vp_width / (2.0*cx);
    world_to_vp[1] = v_win->vp_height / (2.0*cy);
    vp_to_world[0] = 2.0*cx / v_win->vp_width;
    vp_to_world[1] = 2.0*cy / v_win->vp_height;

    /* Translate everything to the drawing plane. */
    glPushMatrix();
    glTranslatef( 0.0, 0.0, zpos );

    /* For the textual information. */
    text_height = 14.0 * vp_to_world[1];
    htextsize( text_height, text_height );

    if ( analy->cur_result != NULL
            && strcmp( analy->cur_result->name, "pvmag" ) == 0 )
    {
        show_dirvec = TRUE;
        map_text_offset = text_height * LINE_SPACING_COEFF;
    }
    else
    {
        show_dirvec = FALSE;
        map_text_offset = 0.0;
    }

    dim = analy->dimension;

    /*
     * Pull out min/max info now since from this we can tell if there was
     * _any_ result evaluated from an enabled material.
     */
    if ( analy->cur_result != NULL )
    {
        if ( !analy->use_global_mm )
        {
            /* Take min/max from current state's data. */

            if ( analy->cur_result->origin.is_node_result )
            {
                el_mm    = analy->state_mm;
                classes  = analy->state_mm_class_long_name;
                sclasses = analy->state_mm_sclass;
                el_id    = analy->state_mm_nodes;
            }
            else
            {
                el_mm    = analy->elem_state_mm.object_minmax;
                classes  = analy->elem_state_mm.class_long_name;
                sclasses = analy->elem_state_mm.sclass;
                el_id    = analy->elem_state_mm.object_id;
            }

            minimum_label = st_min;
            maximum_label = st_max;
        }
        else
        {
            /* Take min/max from global store. */
            if ( analy->cur_result->origin.is_node_result )
            {
                el_mm    = analy->global_mm;
                classes  = analy->global_mm_class_long_name;
                sclasses = analy->global_mm_sclass;
                el_id    = analy->global_mm_nodes;
            }
            else
            {
                el_mm    = analy->elem_global_mm.object_minmax;
                classes  = analy->elem_global_mm.class_long_name;
                sclasses = analy->elem_global_mm.sclass;
                el_id    = analy->elem_global_mm.object_id;
            }

            minimum_label = glob_min;
            maximum_label = glob_max;
        }

        /*
         * Conceptually, we could set found_data by OR'ing in
         * ( analy->cur_result != NULL ), but by setting found_data purely
         * on the basis of whether or not actual data was evaluated we
         * can distinguish the case where the result is valid at the
         * current state but for whatever reason (mainly, all supporting
         * materials were disabled) no data was touched.  In rendering, we
         * will use this distinction to correctly avoid rendering the
         * min/max display (which requires a valid class name pointer)
         * while leaving the colormap as an indicator that the user has
         * created a meaningless state in which there's a current result but
         * nothing available to render it on.
         */
        if ( classes[0] != NULL )
            found_data = TRUE;
        else
            found_data = FALSE;
    }

    /* Compute the fraction size for numbers in the display. */
    if ( analy->cur_result != NULL )
    {
        /* First get the result extremes. */
        raw_minmax = result_has_superclass( analy->cur_result, G_BEAM,
                                            analy )
                     || result_has_superclass( analy->cur_result, G_MAT,
                                               analy )
                     || result_has_superclass( analy->cur_result, G_MESH,
                                               analy )
                     || result_has_superclass( analy->cur_result, G_UNIT,
                                               analy );
        get_min_max( analy, raw_minmax, &low, &high );

        if ( analy->perform_unit_conversion )
        {
            low = (analy->conversion_scale * low)
                  + analy->conversion_offset;
            high = (analy->conversion_scale * high)
                   + analy->conversion_offset;
        }

        /* Now compute the fraction size. */
        fracsz = 0;
        if ( low != high )
        {
            /* Compute legend scale intervals */
            qty_of_intervals = 5;
            linear_variable_scale( low, high, qty_of_intervals,
                                   &scale_minimum, &scale_maximum, &distance,
                                   &scale_error_status );

            /* If scaling succeeded, estimate new quantity of intervals. */
            if ( FALSE == scale_error_status )
            {
                comparison_tolerance = machine_tolerance();
                ntrips = MAX( griz_round( (double)((scale_maximum
                                                    - scale_minimum + distance)
                                                   / distance),
                                          comparison_tolerance ), 0.0 );
            }
            else
            {
                /*
                 * Punt; add 1 because we'll actually want 1 less than
                 * the "ideal" value calculated above.
                 */
                ntrips = qty_of_intervals + 1;
            }

            fracsz = calc_fracsz( low, high, ntrips - 1 );
        }

        if ( fracsz < analy->float_frac_size )
            fracsz = analy->float_frac_size;
    }

    /* Colormap. */
    if ( analy->show_colormap && analy->cur_result != NULL )
    {
        /* Extend the cutoff colors if cutoff is being used. */
        if ( analy->mm_result_set[0] || analy->mm_result_set[1] )
            extend_colormap = TRUE;
        else
            extend_colormap = FALSE;

        /* Width of black border. */
        b_width = 3.5 * vp_to_world[0];
        b_height = 3.5 * vp_to_world[1];

        if ( analy->use_colormap_pos )
        {
            xpos = cx * analy->colormap_pos[0];
            ypos = cy * analy->colormap_pos[1];
            xsize = cx * analy->colormap_pos[2];
            ysize = cy * analy->colormap_pos[3];
        }
        else
        {
            xpos = cx - vp_to_world[0]*40 - b_width;

            /* "255" is somewhat arbitrary and should become a parameter of the
               colorscale size.
                        ypos = cy - vp_to_world[1]*255 - b_height;
            */
            ypos = cy - vp_to_world[1]*240 - b_height;
            xsize = vp_to_world[0]*25;
            ysize = vp_to_world[1]*200;
        }

        /* Put a black border on the colormap. */
        glColor3fv( black );
        glRectf( xpos - b_width, ypos - b_height,
                 xpos + xsize + b_width, ypos + ysize + b_height );

        /* Draw the colormap. */
        nstripes = (int)(ysize*world_to_vp[1] + 0.5);
        yp = ypos;
        if ( extend_colormap )
        {
            /*
             * When using cutoff, may want to extend the cutoff
             * colors a little.
             */
            if ( analy->mm_result_set[0] )
                frst = nstripes/20;
            else
                frst = 0;
            if ( analy->mm_result_set[1] )
                lst = nstripes-nstripes/20;
            else
                lst = nstripes;
            for ( i = 0; i < frst; i++ )
            {
                glColor3fv( v_win->colormap[0] );
                glRectf( xpos, yp, xpos + xsize, yp + vp_to_world[1] );
                yp += vp_to_world[1];
            }
            for ( i = frst; i < lst; i++ )
            {
                glColor3fv( v_win->colormap[ (int)((i-frst)*scl_max /
                                                   (float)(lst-frst-1)) + 1 ] );
                glRectf( xpos, yp, xpos + xsize, yp + vp_to_world[1] );
                yp += vp_to_world[1];
            }
            for ( i = lst; i < nstripes; i++ )
            {
                glColor3fv( v_win->colormap[CMAP_SIZE - 1] );
                glRectf( xpos, yp, xpos + xsize, yp + vp_to_world[1] );
                yp += vp_to_world[1];
            }
        }
        else
        {
            scl_max = (float) CMAP_SIZE - 0.01;

            for ( i = 0; i < nstripes; i++ )
            {
                glColor3fv( v_win->colormap[ (int)(i * scl_max /
                                                   (float)(nstripes - 1)) ] );
                glRectf( xpos, yp, xpos + xsize, yp + vp_to_world[1] );
                yp += vp_to_world[1];
            }
        }

        glColor3fv( v_win->text_color );

        /* Draw the writing (scale) next to the colormap. */
        if ( analy->show_colorscale )
        {
            antialias_lines( TRUE, TRUE );
            glLineWidth( 1.25 );
            hrightjustify( TRUE );

            /* Result title. */
            /* Temporary adjustment.  Should become an explicit function of colorscale
               position.
                        hmove2( xpos + xsize, ypos+ysize+b_height+2.0*text_height );
            */

            hmove2( xpos + xsize,
                    ypos + ysize + b_height + 1.0 * text_height
                    + map_text_offset );
            hcharstr( analy->cur_result->title );

            if ( show_dirvec )
            {
                vec_x = ( analy->vector_result[0] != NULL )
                        ? analy->vector_result[0]->name : zero;
                vec_y = ( analy->vector_result[1] != NULL )
                        ? analy->vector_result[1]->name : zero;
                vec_z = ( analy->vector_result[2] != NULL )
                        ? analy->vector_result[2]->name : zero;

                if ( dim == 3 )
                    sprintf( str, "(%s, %s, %s)", vec_x, vec_y, vec_z );
                else
                    sprintf( str, "(%s, %s)", vec_x, vec_y );

                hmove2( xpos + xsize,
                        ypos + ysize + b_height + 1.0 * text_height );
                hcharstr( str );
            }

            /* Account for min/max result thresholds. */
            if ( TRUE == analy->mm_result_set[0] )
                rmin_offset = 0.05 * ysize;
            else
                rmin_offset = 0.0;

            if ( TRUE == analy->mm_result_set[1] )
                rmax_offset = 0.05 * ysize;
            else
                rmax_offset = 0.0;

            /* Always label the low value. */
            if ( low==high && analy->damage_result )
            {
                scale_y = ((ypos + ysize - rmax_offset) - (ypos + rmin_offset));

                low_text_bound = yp - (0.6 * text_height) + text_height;
                if ( low>0 )
                {
                    xp = xpos - (2.0 * b_width) - (0.6 * text_height);
                    yp = ypos + rmin_offset + (scale_y * high);

                    high_text_bound = yp - (0.6 * text_height);

                    hmove2( xp, yp - (0.6 * text_height) );
                    sprintf( str, "%.*e", fracsz, high );
                    hcharstr( str );
                }
                else
                {
                    xp = xpos - (2.0 * b_width) - (0.6 * text_height);
                    yp = ypos + rmin_offset;

                    hmove2( xp, yp - (0.6 * text_height) );
                    sprintf( str, "%.*e", fracsz, low );
                    hcharstr( str );
                }
            }
            else
            {
                xp = xpos - (2.0 * b_width) - (0.6 * text_height);
                yp = ypos + rmin_offset;

                /* hmove2( xp, yp - (0.6 * text_height) ); */
               if((low > high) || (low == 0.0 && high == 0.0))
                {
                    glColor3fv(v_win->backgrnd_color);
                } else
                {
                    hmove2( xp, yp - (0.6 * text_height) );
                    sprintf( str, "%.*e", fracsz, low );
                    hcharstr( str );
                } 
            }

            xp = xpos - (2.0 * b_width);

            glBegin( GL_TRIANGLES );
            glVertex2f( xp, yp );
            glVertex2f( xp - (0.4 * text_height), yp + (0.2 * text_height) );
            glVertex2f( xp - (0.4 * text_height), yp - (0.2 * text_height) );
            glEnd();

            /* Conditionally render the high value and scale. */
            if ( low != high )
            {
                scale_y = ((ypos + ysize - rmax_offset) - (ypos + rmin_offset))
                          / (high - low);

                low_text_bound = yp - (0.6 * text_height) + text_height;

                /* Label the high value. */
                xp = xpos - (2.0 * b_width) - (0.6 * text_height);
                yp = ypos + rmin_offset + (scale_y * (high - low));

                high_text_bound = yp - (0.6 * text_height);

                hmove2( xp, yp - (0.6 * text_height) );
                sprintf( str, "%.*e", fracsz, high );
                hcharstr( str );

                xp = xpos - (2.0 * b_width);

                glBegin( GL_TRIANGLES );
                glVertex2f( xp, yp );
                glVertex2f( xp - (0.4 * text_height),
                            yp + (0.2 * text_height) );
                glVertex2f( xp - (0.4 * text_height),
                            yp - (0.2 * text_height) );
                glEnd();

                /* If scaling succeeded, render the scale. */
                if ( FALSE == scale_error_status )
                {
                    /* Label scaled values at computed intervals */
                    for ( i = 0; i < ntrips; i++ )
                    {
                        value = scale_minimum + (i * (double) distance);

                        /*
                         * NOTE:  scaled values MAY exceed bounds of tightly
                         * restricted legend limits.
                         */

                        yp = ypos + rmin_offset + (scale_y * (value - low));

                        if ( (low < value)
                                && (value < high)
                                && (low_text_bound <= yp - (0.6 * text_height))
                                && (yp - (0.6 * text_height)
                                    + text_height <= high_text_bound) )
                        {
                            xp = xpos - (2.0 * b_width) - (0.6 * text_height);

                            hmove2( xp, yp - (0.6 * text_height) );
                            sprintf( str, "%.*e", fracsz, value );
                            hcharstr( str );

                            xp = xpos - (2.0 * b_width);

                            glBegin( GL_TRIANGLES );
                            glVertex2f( xp, yp );
                            glVertex2f( xp - (0.4 * text_height),
                                        yp + (0.2 * text_height) );
                            glVertex2f( xp - (0.4 * text_height),
                                        yp - (0.2 * text_height) );
                            glEnd();

                            low_text_bound = yp - (0.6 * text_height)
                                             + text_height;
                        }
                    }
                }
            }

            /* Set up for left side text. */
            hleftjustify( TRUE );
            xp = -cx + text_height;
            yp = cy - (LINE_SPACING_COEFF + TOP_OFFSET) * text_height;

            /* Allow for presence of minmax. */
            if ( analy->show_minmax && found_data )
                yp -= 2.0 * LINE_SPACING_COEFF * text_height;

            if ( ei_labels )
                yp -= LINE_SPACING_COEFF * 2*text_height;

            /* Allow for presence of Displacement Scale */
            if ( analy->show_scale && found_data )
                yp -= LINE_SPACING_COEFF * text_height;
            /*
             * Strain type, reference surface, and reference frame,
             * if applicable.
             */

            if ( analy->cur_result != NULL && found_data )
            {
                p_rs = &analy->cur_result->modifiers;

                if ( p_rs->use_flags.use_strain_variety )
                {
                    hmove2( xp, yp );
                    sprintf( str, "Strain type: %s",
                             strain_label[p_rs->strain_variety] );
                    hcharstr( str );
                    yp -= LINE_SPACING_COEFF * text_height;
                }

                if ( p_rs->use_flags.use_ref_surface )
                {
                    hmove2( xp, yp );
                    sprintf( str, "Surface: %s",
                             ref_surf_label[p_rs->ref_surf] );
                    hcharstr( str );
                    yp -= LINE_SPACING_COEFF * text_height;
                }

                if ( p_rs->use_flags.use_ref_frame )
                {
                    hmove2( xp, yp );
                    sprintf( str, "Ref frame: %s",
                             ref_frame_label[p_rs->ref_frame] );
                    hcharstr( str );
                    yp -= LINE_SPACING_COEFF * text_height;
                    ylines++;
                }

                if ( p_rs->use_flags.coord_transform )
                {
                    hmove2( xp, yp );
                    sprintf( str, "Coord transform: yes" );
                    hcharstr( str );
                    yp -= LINE_SPACING_COEFF * text_height;
                    ylines++;
                }

                if ( p_rs->use_flags.time_derivative )
                {
                    hmove2( xp, yp );
                    sprintf( str, "Result type: Time derivative" );
                    hcharstr( str );
                    yp -= LINE_SPACING_COEFF * text_height;
                    ylines++;
                }

                if ( p_rs->use_flags.use_ref_state )
                {
                    hmove2( xp, yp );
                    if ( p_rs->ref_state == 0 )
                        sprintf( str, "Ref state: initial geom" );
                    else
                        sprintf( str, "Ref state: %d", p_rs->ref_state );
                    hcharstr( str );
                    yp -= LINE_SPACING_COEFF * text_height;
                    ylines++;
                }
            }


            /* Result scale and offset, if applicable. */
            if ( analy->perform_unit_conversion )
            {
                hmove2( xp, yp );
                sprintf( str, "Scale/offset: %.2f/%.2f",
                         analy->conversion_scale, analy->conversion_offset );
                hcharstr( str );
                yp -= LINE_SPACING_COEFF * text_height;
                ylines++;
            }

            antialias_lines( FALSE, 0 );
            glLineWidth( 1.0 );
        }
    }
    else if ( analy->show_colormap )
    {
        antialias_lines( TRUE, TRUE );
        glLineWidth( 1.25 );

        /* Draw the material colors and label them. */
        if ( analy->use_colormap_pos )
        {
            xpos = analy->colormap_pos[0];
            ypos = analy->colormap_pos[1];
            xsize = analy->colormap_pos[2];
            ysize = analy->colormap_pos[3];
        }
        else
        {
            xpos = cx - vp_to_world[0]*40;
            ypos = cy - vp_to_world[1]*250;
            xsize = vp_to_world[0]*25;
            ysize = vp_to_world[1]*200;
        }
        xp = xpos + xsize - text_height;
        yp = ypos + ysize;

        hrightjustify( TRUE );

        p_mesh = MESH_P( analy );

        for ( i = 0; i < p_mesh->material_qty; i++ )
        {
            /* Don't show material in legend if it's invisible. */
            if ( p_mesh->hide_material[i] )
                continue;

            glColor3fv( v_win->mesh_materials.diffuse[i] );
            glRectf( xp, yp, xp + text_height, yp + text_height );

            glColor3fv( v_win->text_color );
            hmove2( xp - 0.5*text_height, yp );
            sprintf( str, "%d", i+1 );
            hcharstr( str );

            yp -= 1.5*text_height;
        }
        hmove2( xp + text_height, ypos + ysize + 1.5*text_height );
        hcharstr( "Materials" );

        antialias_lines( FALSE, 0 );
        glLineWidth( 1.0 );
    }

    antialias_lines( TRUE, TRUE );
    glLineWidth( 1.25 );

    /* File title. */
    if ( analy->show_title )
    {
        glColor3fv( v_win->text_color );
        hcentertext( TRUE );
        hmove2( 0.0, -cy + 2.28 * text_height );
        hcharstr( analy->title );
        hcentertext( FALSE );
    }

    /* File path. */
    if ( analy->show_title_path )
    {
        /* Added current path to plot window */
        hgetfontsize(&cw, &ch);
        htextsize(cw*0.75, ch*0.75);

        hcentertext( TRUE );
        hmove2( 0.0, cy - 0.7 * text_height );
        hcharstr( analy->path_name );
        htextsize(cw, ch);
        hcentertext( FALSE );
    }

    /* Time. */
    if ( analy->show_time )
    {
        char time_name[256], mode_name[32];
        strcpy( mode_name, "State" );
        strcpy( time_name, "t" );
        if ( analy->time_name!=NULL )
            strcpy( time_name, analy->time_name );
        if ( analy->analysis_type == MODAL )
            strcpy( mode_name, "Mode" );

        glColor3fv( v_win->text_color );
        hcentertext( TRUE );
        hmove2( 0.0, -cy + 1.0 * text_height );
        sprintf( str, "%s = %.5e [%s = %d/%d]", time_name, analy->state_p->time, mode_name,
                 analy->cur_state+1, get_max_state(analy)+1);
        hcharstr( str );
        hcentertext( FALSE );

        /* Date and time */

    }

    if ( analy->show_datetime )
    {
        hleftjustify( TRUE );
        hgetfontsize(&cw, &ch);
        htextsize(cw*0.75, ch*0.75);
        hmove2( -cx + .002, -cy + .005);
        tm = time(NULL);
        curr_datetime = localtime(&tm);
        sprintf( str, "%s", asctime(curr_datetime));
        hcharstr( str );
    }


    /* Print an Error Indicator (EI) message if error indicator
     * result is enabled.
     */
    if ( ei_labels )
    {
        glColor3fv( material_colors[15] ); /* Red */
        hcentertext( TRUE );
        hgetfontsize(&cw, &ch);
        htextsize(cw*1.7, ch*1.7);
        hmove2( cx-.5, cy-.038);
        sprintf( str, "Error Indicator for Result:");
        hcharstr( str );
        glColor3fv( v_win->text_color );
    }


    /* Global coordinate system. */
    if ( analy->show_coord )
    {
        leng = 35*vp_to_world[0];
        sub_leng = leng / 10.0;

        /* Rotate the coord system properly, then translate it down
         * to the lower right corner of the view window.
         */
        look_rot_mat( v_win->look_from, v_win->look_at, v_win->look_up, &mat );
        mat_mul( &tmat, &v_win->rot_mat, &mat );
        mat_translate( &mat, cx - 60*vp_to_world[0],
                       -cy + 60*vp_to_world[0], 0.0 );
        mat_mul( &tmat, &tmat, &mat );

        /* Draw the axes. */
        VEC_SET( pt, 0.0, 0.0, 0.0 );
        point_transform( pto, pt, &tmat );
        VEC_SET( pt, leng, 0.0, 0.0 );
        point_transform( pti, pt, &tmat );
        VEC_SET( pt, 0.0, leng, 0.0 );
        point_transform( ptj, pt, &tmat );
        VEC_SET( pt, 0.0, 0.0, leng );
        point_transform( ptk, pt, &tmat );

        glColor3fv( v_win->foregrnd_color );

        glBegin( GL_LINES );
        glVertex3fv( pto );
        glVertex3fv( pti );
        glVertex3fv( pto );
        glVertex3fv( ptj );
        if ( analy->dimension == 3 )
        {
            glVertex3fv( pto );
            glVertex3fv( ptk );
        }
        glEnd();

        if ( show_dirvec )
        {
            VEC_SET( pt, analy->dir_vec[0] * leng * 0.75,
                     analy->dir_vec[1] * leng * 0.75,
                     ( dim == 3 ) ? analy->dir_vec[2] * leng * 0.75
                     : 0.0 );
            point_transform( ptv, pt, &tmat );

            glColor3fv( material_colors[15] ); /* Red */
            glLineWidth( 2.25 );

            glBegin( GL_LINES );
            glVertex3fv( pto );
            glVertex3fv( ptv );
            glEnd();

            glLineWidth( 1.25 );
        }

        /* Label the axes. */
        VEC_SET( pt, leng + sub_leng, sub_leng, sub_leng );
        point_transform( pti, pt, &tmat );
        VEC_SET( pt, sub_leng, leng + sub_leng, sub_leng );
        point_transform( ptj, pt, &tmat );
        if ( analy->dimension == 3 )
        {
            VEC_SET( pt, sub_leng, sub_leng, leng + sub_leng );
            point_transform( ptk, pt, &tmat );
        }

        antialias_lines( FALSE, 0 );
        glLineWidth( 1.0 );

        glColor3fv( v_win->text_color );

        draw_3d_text( pti, "X", FALSE );
        draw_3d_text( ptj, "Y", FALSE );
        if ( analy->dimension == 3 )
            draw_3d_text( ptk, "Z", FALSE );

        antialias_lines( TRUE, TRUE );
        glLineWidth( 1.25 );

        /* Show tensor transformation coordinate system if on. */
        if ( analy->do_tensor_transform
                && analy->tensor_transform_matrix != NULL )
        {
            ttmat = analy->tensor_transform_matrix;

            /* Draw the axis lines. */
            VEC_SET( pt, 0.0, 0.0, 0.0 );
            point_transform( pto, pt, &tmat );
            VEC_SET( pt, ttmat[0][0] * leng, ttmat[1][0] * leng,
                     ttmat[2][0] * leng );
            point_transform( pti, pt, &tmat );
            VEC_SET( pt, ttmat[0][1] * leng, ttmat[1][1] * leng,
                     ttmat[2][1] * leng );
            point_transform( ptj, pt, &tmat );
            VEC_SET( pt, ttmat[0][2] * leng, ttmat[1][2] * leng,
                     ttmat[2][2] * leng );
            point_transform( ptk, pt, &tmat );
            /**/
            /* Hard-code brown/red for now. */
            glColor3f( 0.6, 0.2, 0.0 );

            glBegin( GL_LINES );
            glVertex3fv( pto );
            glVertex3fv( pti );
            glVertex3fv( pto );
            glVertex3fv( ptj );
            if ( dim == 3 )
            {
                glVertex3fv( pto );
                glVertex3fv( ptk );
            }
            glEnd();

            /* Draw the axis labels. */
            VEC_SET( pt, ttmat[0][0] * leng + sub_leng,
                     ttmat[1][0] * leng - sub_leng,
                     ttmat[2][0] * leng - sub_leng );
            point_transform( pti, pt, &tmat );
            VEC_SET( pt, ttmat[0][1] * leng - sub_leng,
                     ttmat[1][1] * leng + sub_leng,
                     ttmat[2][1] * leng - sub_leng );
            point_transform( ptj, pt, &tmat );
            VEC_SET( pt, ttmat[0][2] * leng - sub_leng,
                     ttmat[1][2] * leng - sub_leng,
                     ttmat[2][2] * leng + sub_leng );
            point_transform( ptk, pt, &tmat );

            antialias_lines( FALSE, 0 );
            glLineWidth( 1.0 );

            draw_3d_text( pti, "x", FALSE );
            draw_3d_text( ptj, "y", FALSE );
            if ( dim == 3 )
                draw_3d_text( ptk, "z", FALSE );

            antialias_lines( TRUE, TRUE );
            glLineWidth( 1.25 );
            glColor3fv( v_win->text_color );
        }
    }

    /* Result value min/max. */
    if ( analy->show_minmax  &&  analy->cur_result != NULL && found_data )
    {
        if ( analy->perform_unit_conversion )
        {
            low = (analy->conversion_scale * el_mm[0])
                  + analy->conversion_offset;
            high = (analy->conversion_scale * el_mm[1])
                   + analy->conversion_offset;
        }
        else
        {
            low = el_mm[0];
            high = el_mm[1];
        }

        xp = -cx + text_height;
        yp = cy - (LINE_SPACING_COEFF + TOP_OFFSET) * text_height;

        glColor3fv( v_win->text_color );
        hleftjustify( TRUE );

        hmove2( xp, yp );
        if ( !ei_labels )
        {
            class_ptr = mili_get_class_ptr( analy, sclasses[1], classes[1] );

            class_label = get_class_label( class_ptr, el_id[1] );

            sprintf( str, "%s %.*e, %s %d", maximum_label, fracsz, high,
                     classes[1], class_label );
            hcharstr( str );

            class_ptr = mili_get_class_ptr( analy, sclasses[0], classes[0] );

            class_label = get_class_label( class_ptr, el_id[0] );
            hmove2( xp, yp - LINE_SPACING_COEFF * text_height );
            sprintf( str, "%s %.*e, %s %d", minimum_label, fracsz, low,
                     classes[0], class_label );
            hcharstr( str );
            ylines+=2;
        }
        else
        {


            /* Do not display EI global values until algorithm
             * has been finalized
             */

#ifdef EI_TOTALS

            yp -= 0.025;
            hmove2( xp, yp );
            glColor3fv( material_colors[15] ); /* Red */
            sprintf( str, "%s %4.3e", "EI Global Error Norm:    ", analy->ei_error_norm );
            hcharstr( str );

            hmove2( xp, yp - LINE_SPACING_COEFF * text_height );
            sprintf( str, "%s %4.3e", "EI Global Norm Qty:      ", analy->ei_norm_qty );
            hcharstr( str );

            yp = yp - LINE_SPACING_COEFF * text_height;
            hmove2( xp, yp - LINE_SPACING_COEFF * text_height );
            sprintf( str, "%s %4.3e", "EI Global Error Indicator: ", analy->ei_global_indicator );
            hcharstr( str );

            glColor3fv( v_win->text_color );

            yp = yp - 2*(LINE_SPACING_COEFF * text_height);
#endif
        }
    }

    /* Displacement Scale */
    if ( analy->show_scale  &&  analy->cur_result != NULL && found_data &&
            !ei_labels )
    {
        ylines++;
        xp = -cx + text_height;
        yp = cy - (LINE_SPACING_COEFF + TOP_OFFSET) * text_height;

        /* Allow for presence of minmax. */
        if ( analy->show_minmax && found_data )
            yp -= 2.0 * LINE_SPACING_COEFF * text_height;

        glColor3fv( v_win->text_color );
        hleftjustify( TRUE );

        scalax = analy->displace_scale[0];
        scalay = analy->displace_scale[1];
        scalaz = analy->displace_scale[2];

        hmove2( xp, yp );

        sprintf( str, "Displacement Scale: %.1f/%.1f",
                 scalax, scalay );
        if ( analy->dimension == 3 )
            sprintf( str + strlen( str ), "/%.1f", scalaz );

        hcharstr( str );
    }

    /* Extreme Min/Max */
    if ( analy->extreme_result && analy->cur_result != NULL )
    {
        if ( yp==0.0 )
            yp = cy - (2.0*LINE_SPACING_COEFF * text_height);
        else
        {
            yp =  yp - (LINE_SPACING_COEFF + TOP_OFFSET) * text_height/2;
        }


        glColor3fv( material_colors[15] ); /* Red */

        hcentertext( TRUE );
        hmove2( 0.0, yp );
        hgetfontsize(&cw, &ch);

        htextsize(cw*1.25, ch*1.25);

        if ( analy->extreme_min )
            sprintf( str, "** Extreme Min Result ** " );
        else
            sprintf( str, "** Extreme Max Result ** " );
        hcharstr( str );

        hcentertext( FALSE );
        glColor3fv( v_win->text_color );
    }

    if ( analy->show_tinfo && analy->infoMsgCnt>0 )
    {
        for ( i=0;
                i<analy->infoMsgCnt;
                i++ )
        {
            xp = -cx + text_height;
            yp = cy-((ylines+1) * text_height);

            glColor3fv( v_win->text_color );
            hleftjustify( TRUE );

            hmove2( xp, yp );

            sprintf( str, "Info [%d]: %s", i+1, analy->infoMsg[i] );
            hcharstr( str );
        }
    }

    antialias_lines( FALSE, 0 );
    glLineWidth( 1.0 );


    /* Draw a "safe action area" box for the conversion to NTSC video.
     * The region outside the box may get lost.  This box was defined
     * using trial and error for a 486 x 720 screen size.
     */
    if ( analy->show_safe )
    {
        glColor3fv( red );
        glBegin( GL_LINE_LOOP );
        glVertex2f( -0.91*cx, -0.86*cy );
        glVertex2f( 0.91*cx, -0.86*cy );
        glVertex2f( 0.91*cx, 0.93*cy );
        glVertex2f( -0.91*cx, 0.93*cy );
        glEnd();
    }

    glPopMatrix();

    /* End of foreground overwriting. */
    glDepthFunc( GL_LEQUAL );
}


/************************************************************
 * TAG( quick_draw_mesh_view )
 *
 * Fast display refresh for interactive view control.  Only the
 * edges of the mesh and the foreground are drawn.
 */
void
quick_draw_mesh_view( Analysis *analy )
{
    Transf_mat look_rot;
    float arr[16], scal;

    /* Center the view on a point before rendering. */
    if ( analy->center_view )
        center_view( analy );

    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT |
             GL_STENCIL_BUFFER_BIT );

    glPushMatrix();

    /*
     * Set up all the viewing transforms.  Transformations are
     * specified in opposite order to the order in which they
     * are applied.
     */
    look_rot_mat( v_win->look_from, v_win->look_at, v_win->look_up, &look_rot );
    mat_to_array( &look_rot, arr );
    glMultMatrixf( arr );
    glTranslatef( -v_win->look_from[0], -v_win->look_from[1],
                  -v_win->look_from[2] );
    glTranslatef( v_win->trans[0], v_win->trans[1], v_win->trans[2] );
    mat_to_array( &v_win->rot_mat, arr );
    glMultMatrixf( arr );
    scal = v_win->bbox_scale;
    glScalef( scal*v_win->scale[0], scal*v_win->scale[1],
              scal*v_win->scale[2] );
    glTranslatef( v_win->bbox_trans[0], v_win->bbox_trans[1],
                  v_win->bbox_trans[2] );

    /* Draw the grid edges. */
    if ( MESH_P( analy )->edge_list != NULL
            && MESH_P( analy )->edge_list->size > 0 &&
            !analy->particle_nodes_enabled)
    {
        if ( analy->dimension == 3 )
            draw_edges_3d( analy );
        else
            draw_edges_2d( analy );
    }
    else
        draw_bbox( analy->bbox );

    glPopMatrix();

    /* Draw the foreground. */
    if ( analy->show_bbox || analy->show_coord || analy->show_time ||
            analy->show_title || analy->show_safe || analy->show_colormap )
    {
        draw_foreground( analy );
    }

    gui_swap_buffers();
}


/*
 * SECTION_TAG( Image saving )
 */


/************************************************************
 * TAG( rgb_to_screen )
 *
 * Load an rgb file image into memory.
 */
void
rgb_to_screen( char *fname, Bool_type background, Analysis *analy )
{
    IMAGE *img;
    int img_width, img_height, channels;
    int i, j;
    unsigned char *raster, *prgb_byte;
    unsigned short *rbuf, *gbuf, *bbuf, *abuf;
    unsigned short *pr_short, *pg_short, *pb_short, *pa_short;
    Bool_type alpha;
    RGB_raster_obj *p_rro;

    img = iopen( fname, "r", 0, 0, 0, 0, 0 );
    if ( img == NULL )
    {
        popup_dialog( INFO_POPUP,
                      "Unable to open rgb image file \"%s\".", fname );
        return;
    }

    channels = img->zsize;
    if ( channels != 3 && channels != 4 )
    {
        popup_dialog( INFO_POPUP, "File \"%s\" %s%s%s", fname,
                      "does not contain a 3 channel (RGB) or\n",
                      "4 channel (RGBA) image.  Other image types are not",
                      " supported." );
        iclose( img );
        return;
    }

    alpha = ( channels == 4 );
    img_width = img->xsize;
    img_height = img->ysize;

    raster = NEW_N( unsigned char, img_width * img_height * channels,
                    "Image memory" );

    rbuf = NEW_N( unsigned short, img_width, "RGB red row buf" );
    gbuf = NEW_N( unsigned short, img_width, "RGB green row buf" );
    bbuf = NEW_N( unsigned short, img_width, "RGB blue row buf" );
    if ( alpha )
        abuf = NEW_N( unsigned short, img_width, "RGB alpha row buf" );

    /*
     * Multiplex r, g, b, and alpha values into the raster buffer
     * on a row-by-row basis.
     */
    prgb_byte = raster;
    for( j = 0; j < img_height; j++ )
    {
        getrow( img, rbuf, j, 0 );
        getrow( img, gbuf, j, 1 );
        getrow( img, bbuf, j, 2 );
        if ( alpha )
            getrow( img, abuf, j, 3 );

        pr_short = rbuf;
        pg_short = gbuf;
        pb_short = bbuf;
        if(alpha)
        {
            pa_short = abuf;
        }

        for( i = 0; i < img_width; i++ )
        {
            *prgb_byte++ = (unsigned char) *pr_short++;
            *prgb_byte++ = (unsigned char) *pg_short++;
            *prgb_byte++ = (unsigned char) *pb_short++;
            if ( alpha )
                *prgb_byte++ = (unsigned char) *pa_short++;
        }
    }

    iclose( img );

    free( rbuf );
    free( gbuf );
    free( bbuf );
    if ( alpha )
        free( abuf );

    if ( background )
    {
        /*
         * Prepare to use during background initialization, but don't
         * load the image now.
         */
        if ( analy->background_image != NULL )
        {
            p_rro = analy->background_image;
            free( p_rro->raster );
        }
        else
            p_rro = NEW( RGB_raster_obj, "Background raster" );

        p_rro->raster = raster;
        p_rro->img_width = img_width;
        p_rro->img_height = img_height;
        p_rro->alpha = alpha;

        analy->background_image = p_rro;
    }
    else
    {
        /* Not for a background, so load it now and release it. */
        memory_to_screen( TRUE, img_width, img_height, alpha, raster );
        free( raster );
    }
}


/************************************************************
 * TAG( memory_to_screen )
 *
 * Load the framebuffer with pixels from an in-memory image
 * raster.
 */
static void
memory_to_screen( clear, width, height, alpha, raster )
Bool_type clear;
int width;
int height;
Bool_type alpha;
unsigned char *raster;
{
    int upa;

    /* Save the current value of the unpack alignment. */
    glGetIntegerv( GL_UNPACK_ALIGNMENT, &upa );

    glPixelStorei( GL_UNPACK_ALIGNMENT, 1 );

    /* glDrawBuffer( GL_BACK ); */

    /*
     * Clear everything to the default background color in case the
     * image has a smaller dimension than the current window.
     */
    if ( clear )
        glClear( GL_COLOR_BUFFER_BIT );

    /* Now copy the pixels... */
    glDrawPixels( width, height, (alpha ? GL_RGBA : GL_RGB), GL_UNSIGNED_BYTE,
                  raster );

    if ( clear )
        gui_swap_buffers();

    /* Reset the unpack alignment if necessary. */
    if ( upa != 1 )
        glPixelStorei( GL_UNPACK_ALIGNMENT, upa );
}


/************************************************************
 * TAG( screen_to_memory )
 *
 * Copy the screen image to an array in memory.  The save_alpha
 * parameter specifies whether to save the alpha channel.  The
 * return array is given by the parameter ipdat, and should be of
 * size xsize*ysize*4 if the alpha is requested and xsize*ysize*3
 * otherwise.
 */
void
screen_to_memory( Bool_type save_alpha, int xsize, int ysize,
                  unsigned char *ipdat )
{
    if ( v_win->vp_width != xsize || v_win->vp_height != ysize )
    {
        popup_dialog( WARNING_POPUP,
                      "Viewport width and height don't match memory dimensions." );
        return;
    }

    glReadBuffer( GL_FRONT );

    if ( save_alpha )
        glReadPixels( 0, 0, xsize, ysize, GL_RGBA, GL_UNSIGNED_BYTE, ipdat );
    else
        glReadPixels( 0, 0, xsize, ysize, GL_RGB, GL_UNSIGNED_BYTE, ipdat );
}


/************************************************************
 * TAG( screen_to_rgb )
 *
 * Simple snap to rle-encoded rgb/rgba file.
 */
void
screen_to_rgb( char fname[], Bool_type alpha )
{
    unsigned char *ipdat, *prgb_byte;
    unsigned short *rbuf, *gbuf, *bbuf, *abuf;
    unsigned short *pr_short, *pg_short, *pb_short, *pa_short;
    int i, j;
    int vp_width, vp_height;
    int bufsiz;
    int components;
    IMAGE *img;

    components = alpha ? 4 : 3;

    vp_width = v_win->vp_width;
    vp_height = v_win->vp_height;

    img = iopen( fname, "w", RLE( 1 ), 3, vp_width, vp_height, components );
    if ( img == 0 )
    {
        popup_dialog( WARNING_POPUP, "Unable to open rgb image file \"%s\".\n",
                      fname );
        return;
    }

    bufsiz = vp_width * vp_height * components;
    ipdat = NEW_N( unsigned char, bufsiz, "Screen to rgb" );

    /* Read the pixel data into ipdat. */
    screen_to_memory( alpha, vp_width, vp_height, ipdat );

    wrt_text( "Writing data to image file %s\n\n", fname );

    rbuf = NEW_N( unsigned short, vp_width, "RGB red row buf" );
    gbuf = NEW_N( unsigned short, vp_width, "RGB green row buf" );
    bbuf = NEW_N( unsigned short, vp_width, "RGB blue row buf" );
    if ( alpha )
        abuf = NEW_N( unsigned short, vp_width, "RGB alpha row buf" );

    /*
     * Demultiplex r, g, b, and alpha values into separate buffers
     * for output to rgb file on a raster row-by-row basis.
     */
    prgb_byte = ipdat;
    for( j = 0; j < vp_height; j++ )
    {
        pr_short = rbuf;
        pg_short = gbuf;
        pb_short = bbuf;
        pa_short = abuf;

        for( i = 0; i < vp_width; i++ )
        {
            *pr_short++ = (unsigned short) *prgb_byte++;
            *pg_short++ = (unsigned short) *prgb_byte++;
            *pb_short++ = (unsigned short) *prgb_byte++;
            if ( alpha )
                *pa_short++ = (unsigned short) *prgb_byte++;
        }

        putrow( img, rbuf, j, 0 );
        putrow( img, gbuf, j, 1 );
        putrow( img, bbuf, j, 2 );
        if ( alpha )
            putrow( img, abuf, j, 3 );
    }

    iclose( img );

    free( ipdat );
    free( rbuf );
    free( gbuf );
    free( bbuf );
    if ( alpha )
        free( abuf );
}


/************************************************************
 * TAG( screen_to_ps )
 *
 * Simple snap to postscript raster file.  By Doug Speck.
 *
 * Utility funcs prologue, epilogue, and puthexpix swiped from Utah Raster
 * Toolkit program "rletops", necessitating copyright info below.
 */

/*
 * rletops.c - Convert RLE to postscript file.
 *
 * Author:      Rod Bogart & John W. Peterson
 *              Computer Science Dept.
 *              University of Utah
 * Date:        Tue Nov  4 1986
 * Copyright (c) 1986 Rod Bogart
 *
 * Modified by: Gregg Townsend
 *              Department of Computer Science
 *              University of Arizona
 * Date:        June 22, 1990
 * Changes:     50% speedup using putc(), add -C option, translate to page ctr
 *              Fix comments to conform to EPS Version 2.0  (Al Clark)
 *
 * Based on "tobw.c" by Spencer Thomas, and
 * "rps" by Marc Majka (UBC Vision lab)
 *
 * EPSF code from Mike Zyda (Naval Postgraduate School)
 */

#define DFLT_PAGE_HEIGHT        (11.0)
#define DFLT_PAGE_WIDTH         (8.5)
#define DFLT_MARGIN             (1.0)

static int gencps = 1;          /* Only generate color PostScript. */

void
screen_to_ps( char fname[] )
{
    unsigned char *ipdat;
    int vp_width, vp_height;
    int nrow, ncol;
    int i, j, rotate, add_line;
    FILE *outfile;
    float page_height = DFLT_PAGE_HEIGHT,
          page_width = DFLT_PAGE_WIDTH,
          margin = DFLT_MARGIN,
          max_prn_height, max_prn_width,
          prn_aspect, img_aspect,
          prn_img_height, prn_img_width,
          page_ymin, page_ymax, page_xmin, page_xmax;

    if ( (outfile = fopen( fname, "w")) == 0 )
    {
        popup_dialog( WARNING_POPUP, "Unable to open output file %s", fname );
        return;
    }

    /* calc aspect ratios for print page and screen image */
    max_prn_height = page_height - 2.0 * margin;
    max_prn_width = page_width - 2.0 * margin;
    prn_aspect = max_prn_height / max_prn_width;

    vp_width = v_win->vp_width;
    vp_height = v_win->vp_height;

    img_aspect = vp_height / (float) vp_width;

    /* Compare aspect ratios to see if rotation necessary. */
    if ( prn_aspect >= 1.0 && img_aspect < 1.0 ||
            prn_aspect < 1.0 && img_aspect >= 1.0 )
    {
        rotate = TRUE;
        img_aspect = 1.0 / img_aspect;
        nrow = vp_width;
        ncol = vp_height;
        add_line = vp_width & 0x1;
    }
    else
    {
        rotate = FALSE;
        nrow = vp_height;
        ncol = vp_width;
        add_line = vp_height & 0x1;
    }

    /* scale image as large as possible */
    if (img_aspect <= prn_aspect)
    {
        prn_img_width = max_prn_width;
        prn_img_height = img_aspect * max_prn_width;
    }
    else
    {
        prn_img_height = max_prn_height;
        prn_img_width = max_prn_height / img_aspect;
    }

    ipdat = NEW_N( unsigned char, vp_width*vp_height*3, "Screen to ps" );

    /* Read the pixel data into ipdat. */
    screen_to_memory( FALSE, vp_width, vp_height, ipdat );

    /* Calc bounding box dimensions (inches). */
    page_xmin = (page_width - prn_img_width) / 2.0;
    page_xmax = page_xmin + prn_img_width;
    page_ymin = (page_height - prn_img_height) / 2.0;
    page_ymax = page_ymin + prn_img_height;

    /* Write to postscript file. */
    prologue( outfile, 0, nrow + add_line, ncol, page_xmin, page_ymin,
              page_xmax, page_ymax);

    /* Write image raster. */
    if (rotate)
    {
        for( j = 0; j < nrow; j++)
            for( i = ncol - 1; i >= 0; i--)
            {
                puthexpix( outfile, *(ipdat + i*nrow*3 + j*3) );
                puthexpix( outfile, *(ipdat + i*nrow*3 + j*3 + 1) );
                puthexpix( outfile, *(ipdat + i*nrow*3 + j*3 + 2) );
            }
    }
    else
    {
        for( i = 0; i < nrow; i++)
            for( j = 0; j < ncol; j++)
            {
                puthexpix( outfile, *(ipdat + i*ncol*3 + j*3) );
                puthexpix( outfile, *(ipdat + i*ncol*3 + j*3 + 1) );
                puthexpix( outfile, *(ipdat + i*ncol*3 + j*3 + 2) );
            }
    }

    /* Laserwriters need an even number of scan lines (from rletops). */
    if (add_line)
        for( j = 0; j < ncol; j++)
        {
            puthexpix( outfile, 0xFF);
            puthexpix( outfile, 0xFF);
            puthexpix( outfile, 0xFF);
        }

    epilogue( outfile, 0 );

    fclose( outfile );
    free( ipdat );
}

static void
prologue( FILE *outfile, int scribe_mode, int nr, int nc, float x1, float y1,
          float x2, float y2)
{
    fprintf(outfile,"%%!\n");
    fprintf(outfile, "%%%%BoundingBox: %d %d %d %d\n",
            (int)(x1 * 72), (int)(y1 * 72),
            (int)(x2 * 72), (int)(y2 * 72));
    fprintf(outfile, "%%%%EndComments\n");
    fprintf(outfile,"gsave\n");
    if ( !scribe_mode )
        fprintf(outfile,"initgraphics\n");
    fprintf(outfile,"72 72 scale\n");
    fprintf(outfile,"/imline %d string def\n",nc*2*(gencps?3:1));
    fprintf(outfile,"/drawimage {\n");
    fprintf(outfile,"    %d %d 8\n",nc,nr);
    fprintf(outfile,"    [%d 0 0 %d 0 %d]\n",nc,-1*nr,nr);
    fprintf(outfile,"    { currentfile imline readhexstring pop } ");
    if ( gencps )
        fprintf(outfile,"false 3 colorimage\n");
    else
        fprintf(outfile,"image\n");
    fprintf(outfile,"} def\n");
    fprintf(outfile,"%f %f translate\n",x1,y2);
    fprintf(outfile,"%f %f scale\n",x2-x1,y1-y2);
    fprintf(outfile,"drawimage\n");
}

static void
epilogue( FILE *outfile, int scribemode )
{
    fprintf(outfile,"\n");
    if (!scribemode)
        fprintf(outfile,"showpage\n");
    fprintf(outfile,"grestore\n");
}

static void
puthexpix( FILE *outfile, unsigned char p )
{
    static int npixo = 0;
    static char tohex[] = "0123456789ABCDEF";

    putc(tohex[(p>>4)&0xF],outfile);
    putc(tohex[p&0xF],outfile);
    npixo += 1;
    if ( npixo >= 32 )
    {
        putc('\n',outfile);
        npixo = 0;
    }
}


/*
 * SECTION_TAG( Colormap )
 */


/************************************************************
 * TAG( read_text_colormap )
 *
 * Read a colormap from a text file.  The file should contain
 * CMAP_SIZE float triples for r,g,b, in the range [0, 1].
 */
void
read_text_colormap( char *fname )
{
    FILE *infile;
    float r, g, b;
    int i;

    if ( ( infile = fopen(fname, "r") ) == NULL )
    {
        popup_dialog( WARNING_POPUP, "Unable to open file %s", fname );
        return;
    }

    for ( i = 0; i < CMAP_SIZE; i++ )
    {
        fscanf( infile, "%f%f%f", &r, &g, &b );
        v_win->colormap[i][0] = r;
        v_win->colormap[i][1] = g;
        v_win->colormap[i][2] = b;
        if ( feof( infile ) )
        {
            fclose( infile );
            popup_dialog( WARNING_POPUP,
                          "Palette file end reached too soon." );
            hot_cold_colormap();
            return;
        }
    }

    /* Save the colormap's native min and max colors. */
    VEC_COPY( v_win->cmap_min_color, v_win->colormap[0] );
    VEC_COPY( v_win->cmap_max_color, v_win->colormap[CMAP_SIZE - 1] );

    fclose( infile );
}


/************************************************************
 * TAG( hot_cold_colormap )
 *
 * Fills the colormap array with a hot-to-cold colormap.
 * By Jeffery W. Long.
 */
void
hot_cold_colormap( void )
{
    float delta, rv, gv, bv;
    int group_sz, i, j;

    group_sz = CMAP_SIZE / 4;   /* Number of indices per group. */
    i = 0;                      /* Next color index. */
    delta = 1.0 / (group_sz-1);

    /* Blue to cyan group. */
    rv = 0.0;
    gv = 0.0;
    bv = 1.0;
    for ( j = 0; j < group_sz; j++ )
    {
        v_win->colormap[i][0] = rv;
        v_win->colormap[i][1] = gv;
        v_win->colormap[i][2] = bv;
        gv = MIN( 1.0, gv + delta );
        i++;
    }

    /* Cyan to Green group. */
    gv = 1.0;
    bv = 1.0;
    for ( j = 0; j < group_sz; j++ )
    {
        v_win->colormap[i][0] = rv;
        v_win->colormap[i][1] = gv;
        v_win->colormap[i][2] = bv;
        bv = MAX( 0.0, bv - delta );
        i++;
    }

    /* Green to yellow group. */
    bv = 0.0;
    for ( j = 0; j < group_sz; j++ )
    {
        v_win->colormap[i][0] = rv;
        v_win->colormap[i][1] = gv;
        v_win->colormap[i][2] = bv;
        rv = MIN( 1.0, rv + delta );
        i++;
    }

    /* Yellow to red group. */
    rv = 1.0;
    for ( j = 0; j < group_sz; j++ )
    {
        v_win->colormap[i][0] = rv;
        v_win->colormap[i][1] = gv;
        v_win->colormap[i][2] = bv;
        gv = MAX( 0.0, gv - delta );
        i++;
    }

    /* Take care of any remainder. */
    while ( i < CMAP_SIZE )
    {
        v_win->colormap[i][0] = rv;
        v_win->colormap[i][1] = gv;
        v_win->colormap[i][2] = bv;
        i++;
    }
    
    /* set the initial colormap to go back from grayscale to the initial colormap */
    for(i = 0; i < CMAP_SIZE; i++)
    {
       v_win->initial_colormap[i][0] = v_win->colormap[i][0];
       v_win->initial_colormap[i][1] = v_win->colormap[i][1];
       v_win->initial_colormap[i][2] = v_win->colormap[i][2];
    }

    /* Save the colormap's native min and max colors. */
    VEC_COPY( v_win->cmap_min_color, v_win->colormap[0] );
    VEC_COPY( v_win->cmap_max_color, v_win->colormap[CMAP_SIZE - 1] );
}


/************************************************************
 * TAG( invert_colormap )
 *
 * Invert the colormap.  Just flip it upside down.
 */
void
invert_colormap( void )
{
    float tmp;
    int i, j;
    int loop_max;
    int idx_max;

    loop_max = CMAP_SIZE / 2 - 1;
    idx_max = CMAP_SIZE - 1;

    for ( i = 0; i <= loop_max; i++ )
        for ( j = 0; j < 3; j++ )
        {
            tmp = v_win->colormap[i][j];
            v_win->colormap[i][j] = v_win->colormap[idx_max-i][j];
            v_win->colormap[idx_max-i][j] = tmp;
        }

    /* Update the cached native min and max colors. */
    VEC_COPY( v_win->cmap_min_color, v_win->colormap[0] );
    VEC_COPY( v_win->cmap_max_color, v_win->colormap[idx_max] );
}


/************************************************************
 * TAG( set_cutoff_colors )
 *
 * Load cutoff threshold colors into the colormap or reset
 * those colormap entries to their native values.
 */
void
set_cutoff_colors( Bool_type set, Bool_type cut_low, Bool_type cut_high )
{
    if ( set )
    {
        /* Set the low cutoff color. */
        if ( cut_low )
            VEC_COPY( v_win->colormap[0], v_win->rmin_color );

        /* Set the high cutoff color. */
        if ( cut_high )
            VEC_COPY( v_win->colormap[CMAP_SIZE - 1], v_win->rmax_color );
    }
    else
    {
        /* Reset the low cutoff color. */
        if ( cut_low )
            VEC_COPY( v_win->colormap[0], v_win->cmap_min_color );

        /* Reset the high cutoff color. */
        if ( cut_high )
            VEC_COPY( v_win->colormap[CMAP_SIZE - 1], v_win->cmap_max_color );
    }
}


/************************************************************
 * TAG( cutoff_colormap )
 *
 * Fills the given array with a cutoff colormap.  The parameters
 * tell whether the low or high or both cutoffs are active.
 * The cutoff colormap is a hot-to-cold map with out-of-bounds
 * values at the ends.
 */
void
cutoff_colormap( Bool_type cut_low, Bool_type cut_high )
{
    float delta, rv, gv, bv;
    int group_sz, i, j;

    group_sz = CMAP_SIZE / 4;   /* Number of indices per group. */
    i = 0;                      /* Next color index. */
    delta = 1.0 / (group_sz-1);

    /* Blue to cyan group. */
    rv = 0.0;
    gv = 0.0;
    bv = 1.0;
    for ( j = 0; j < group_sz; j++ )
    {
        v_win->colormap[i][0] = rv;
        v_win->colormap[i][1] = gv;
        v_win->colormap[i][2] = bv;
        gv = MIN( 1.0, gv + delta );
        i++;
    }

    /* Cyan to Green group. */
    gv = 1.0;
    bv = 1.0;
    for ( j = 0; j < group_sz; j++ )
    {
        v_win->colormap[i][0] = rv;
        v_win->colormap[i][1] = gv;
        v_win->colormap[i][2] = bv;
        bv = MAX( 0.0, bv - delta );
        i++;
    }

    /* Green to yellow group. */
    bv = 0.0;
    for ( j = 0; j < group_sz; j++ )
    {
        v_win->colormap[i][0] = rv;
        v_win->colormap[i][1] = gv;
        v_win->colormap[i][2] = bv;
        rv = MIN( 1.0, rv + delta );
        i++;
    }

    /* Yellow to red group. */
    rv = 1.0;
    for ( j = 0; j < group_sz; j++ )
    {
        v_win->colormap[i][0] = rv;
        v_win->colormap[i][1] = gv;
        v_win->colormap[i][2] = bv;
        gv = MAX( 0.0, gv - delta );
        i++;
    }

    /* Take care of any remainder. */
    while ( i < CMAP_SIZE )
    {
        v_win->colormap[i][0] = rv;
        v_win->colormap[i][1] = gv;
        v_win->colormap[i][2] = bv;
        i++;
    }

    /* Set the low cutoff color. */
    if ( cut_low )
        VEC_COPY( v_win->colormap[0], v_win->rmin_color );

    /* Set the high cutoff color. */
    if ( cut_high )
        VEC_COPY( v_win->colormap[CMAP_SIZE - 1], v_win->rmax_color );
}

/************************************************************
 * TAG( restore_colormap )
 *
 * restores the colormap to that which was the colormap on startup 
 * By William Oliver 
 */
void
restore_colormap()
{
    int sz, i;
    
    sz = CMAP_SIZE;
    for(i = 0; i < sz; i++)
    {
        v_win->colormap[i][0] = v_win->initial_colormap[i][0];      
        v_win->colormap[i][1] = v_win->initial_colormap[i][1];      
        v_win->colormap[i][2] = v_win->initial_colormap[i][2];      
    }

    VEC_COPY( v_win->cmap_min_color, v_win->colormap[0] );
    VEC_COPY( v_win->cmap_max_color, v_win->colormap[CMAP_SIZE - 1] );

}

/************************************************************
 * TAG( gray_colormap )
 *
 * Fills the colormap array with a grayscale or inverted
 * grayscale colormap.
 * By Jeffery W. Long.
 */
void
gray_colormap( Bool_type inverse )
{
    float start, step;
    int sz, i;

    sz = CMAP_SIZE;

    if ( inverse )
    {
        step = -1.0 / (float)(sz-1);
        start = 1.0;
    }
    else
    {
        step = 1.0 / (float)(sz-1);
        start = 0.0;
    }

    for ( i = 0; i < sz; i++ )
    {
        v_win->colormap[i][0]          = start + i*step;
        v_win->colormap[i][1]          = start + i*step;
        v_win->colormap[i][2]          = start + i*step;
    }

    /* Save the colormap's native min and max colors. */
    VEC_COPY( v_win->cmap_min_color, v_win->colormap[0] );
    VEC_COPY( v_win->cmap_max_color, v_win->colormap[CMAP_SIZE - 1] );
}


/************************************************************
 * TAG( contour_colormap )
 *
 * Create N contours in the colormap by sampling the current
 * colormap and creating color bands.
 */
void
contour_colormap( int ncont )
{
    float band_sz, col[3], thresh;
    int idx, lo, hi, max;
    int i, j;

    max = CMAP_SIZE - 2;
    thresh = (float) max - 0.001;
    band_sz = thresh / ncont;

    for ( i = 0; i < ncont; i++ )
    {
        idx = (int)(thresh * i / (ncont - 1.0)) + 1;
        VEC_COPY( col, v_win->colormap[idx] );

        lo = i * band_sz + 1;
        hi = (i+1) * band_sz + 1;
        if ( hi > max )
            hi = max;

        for ( j = lo; j <= hi; j++ )
        {
            VEC_COPY( v_win->colormap[j], col );
        }
    }
}


/*
 * SECTION_TAG( Lights )
 */


/************************************************************
 * TAG( delete_lights )
 *
 * Delete any currently bound lights.
 */
void
delete_lights( void )
{
    int i;
    static GLenum lights[] =
    {
        GL_LIGHT0, GL_LIGHT1, GL_LIGHT2, GL_LIGHT3, GL_LIGHT4, GL_LIGHT5
    };

    for ( i = 0; i < 6; i++ )
    {
        glDisable( lights[i] );
        v_win->light_active[i] = FALSE;
    }
}


/************************************************************
 * TAG( set_light )
 *
 * Define and enable a light.
 */
void
set_light( int token_cnt, char tokens[MAXTOKENS][TOKENLENGTH] )
{
    float vec[4], val;
    int light_num, i;
    GLenum light_id;

    /* Get the light number. */
    sscanf( tokens[1], "%d", &light_num );
    --light_num;
    v_win->light_active[light_num] = TRUE;

    switch ( light_num )
    {
    case 0:
        light_id = GL_LIGHT0;
        break;
    case 1:
        light_id = GL_LIGHT1;
        break;
    case 2:
        light_id = GL_LIGHT2;
        break;
    case 3:
        light_id = GL_LIGHT3;
        break;
    case 4:
        light_id = GL_LIGHT4;
        break;
    case 5:
        light_id = GL_LIGHT5;
        break;
    default:
        popup_dialog( INFO_POPUP, "Light number %d out of range",
                      light_num + 1 );
        return;
    }

    /* Process the light position. */
    sscanf( tokens[2], "%f", &vec[0] );
    sscanf( tokens[3], "%f", &vec[1] );
    sscanf( tokens[4], "%f", &vec[2] );
    sscanf( tokens[5], "%f", &vec[3] );
    glLightfv( light_id, GL_POSITION, vec );
    v_win->light_pos[light_num][0] = vec[0];
    v_win->light_pos[light_num][1] = vec[1];
    v_win->light_pos[light_num][2] = vec[2];
    v_win->light_pos[light_num][3] = vec[3];

    for ( i = 6; i < token_cnt; i++ )
    {
        if ( strcmp( tokens[i], "amb" ) == 0 )
        {
            sscanf( tokens[i+1], "%f", &vec[0] );
            sscanf( tokens[i+2], "%f", &vec[1] );
            sscanf( tokens[i+3], "%f", &vec[2] );
            vec[3] = 1.0;
            glLightfv( light_id, GL_AMBIENT, vec );
            i += 3;
        }
        else if ( strcmp( tokens[i], "diff" ) == 0 )
        {
            sscanf( tokens[i+1], "%f", &vec[0] );
            sscanf( tokens[i+2], "%f", &vec[1] );
            sscanf( tokens[i+3], "%f", &vec[2] );
            vec[3] = 1.0;
            glLightfv( light_id, GL_DIFFUSE, vec );
            i += 3;
        }
        else if ( strcmp( tokens[i], "spec" ) == 0 )
        {
            sscanf( tokens[i+1], "%f", &vec[0] );
            sscanf( tokens[i+2], "%f", &vec[1] );
            sscanf( tokens[i+3], "%f", &vec[2] );
            vec[3] = 1.0;
            glLightfv( light_id, GL_SPECULAR, vec );
            i += 3;
        }

        else if ( strcmp( tokens[i], "spotdir" ) == 0 )
        {
            sscanf( tokens[i+1], "%f", &vec[0] );
            sscanf( tokens[i+2], "%f", &vec[1] );
            sscanf( tokens[i+3], "%f", &vec[2] );
            glLightfv( light_id, GL_SPOT_DIRECTION, vec );
            i += 3;
        }
        else if ( strcmp( tokens[i], "spot" ) == 0 )
        {
            sscanf( tokens[i+1], "%f", &val );
            glLightf( light_id, GL_SPOT_EXPONENT, val );
            i++;
            sscanf( tokens[i+1], "%f", &val );
            glLightf( light_id, GL_SPOT_CUTOFF, val );
            i++;
        }
        else
            popup_dialog( INFO_POPUP, "Invalid light parameter \"%s\" ignored",
                          tokens[i] );
    }
    glEnable( light_id );
}


/************************************************************
 * TAG( move_light )
 *
 * Move a light along the specified axis (X == 0, etc.).
 */
void
move_light( int light_num, int axis, float incr )
{
    GLenum light_id;

    if ( light_num < 0 || light_num > 5 )
    {
        popup_dialog( INFO_POPUP, "Light number %d out of range",
                      light_num + 1 );
        return;
    }

    v_win->light_pos[light_num][axis] += incr;

    switch ( light_num )
    {
    case 0:
        light_id = GL_LIGHT0;
        break;
    case 1:
        light_id = GL_LIGHT1;
        break;
    case 2:
        light_id = GL_LIGHT2;
        break;
    case 3:
        light_id = GL_LIGHT3;
        break;
    case 4:
        light_id = GL_LIGHT4;
        break;
    case 5:
        light_id = GL_LIGHT5;
        break;
    }

    glLightfv( light_id, GL_POSITION, v_win->light_pos[light_num] );
}


/*
 * SECTION_TAG( Title screens )
 */


/************************************************************
 * TAG( set_vid_title )
 *
 * Set a video title line.  The argument nline is the number
 * of the line to be set (0-3).
 */
void
set_vid_title( int nline, char *title )
{
    if ( nline >= 0 && nline < 4 )
        strcpy( v_win->vid_title[nline], title );
    else
        popup_dialog( INFO_POPUP, "Illegal line number" );
}


/************************************************************
 * TAG( draw_vid_title )
 *
 * Put up a video title screen.
 */
void
draw_vid_title( void )
{
    float zpos, cx, cy, text_height;
    int i;

    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glColor3fv( v_win->text_color );

    antialias_lines( TRUE, TRUE );
    glLineWidth( 1.5 );

    /* Get drawing window and position. */
    get_foreground_window( &zpos, &cx, &cy );

    text_height = 0.1*cy;
    hcentertext( TRUE );
    htextsize( text_height, text_height );

    for ( i = 0; i < 4; i++ )
    {
        hmove( 0.0, 0.45*cy - i*0.2*cy, zpos );
        hcharstr( v_win->vid_title[i] );
    }

    hmove( 0.0, -0.45*cy, zpos );
    hcharstr( env.date );

    hcentertext( FALSE );

    antialias_lines( FALSE, 0 );
    glLineWidth( 1.0 );

    gui_swap_buffers();
}


/************************************************************
 * TAG( copyright )
 *
 * Put up the program name and copyright screen.
 */
void
copyright( void )
{
    float pos[3], cx, cy, vp_to_world_y, text_height;

    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glColor3fv( v_win->text_color  );

    antialias_lines( TRUE, TRUE );
    glLineWidth( 1.5 );

    /* Get drawing window and position. */
    get_foreground_window( &pos[2], &cx, &cy );

    pos[0] = 0.0;

    vp_to_world_y = 2.0*cy / v_win->vp_height;
    text_height = 18.0 * vp_to_world_y;
    hcentertext( TRUE );

    htextsize( 3.0*text_height, 3.0*text_height );
    pos[1] = cy - 5.0*text_height;
    hmove( pos[0], pos[1], pos[2] );
    hcharstr( "GRIZ" );

    pos[1] -= 3.0*text_height;
    hmove( pos[0], pos[1], pos[2] );
    htextsize( text_height, text_height );
    hcharstr( "Visualization of Finite Element Analysis Results" );
    pos[1] -= 1.5*text_height;
    hmove( pos[0], pos[1], pos[2] );
    hcharstr( "" );

    pos[1] -= 3.0*text_height;
    hmove( pos[0], pos[1], pos[2] );
    hcharstr( "Methods Development Group" );
    pos[1] -= 1.5*text_height;
    hmove( pos[0], pos[1], pos[2] );
    hcharstr( "Lawrence Livermore National Laboratory" );

    pos[1] -= 2.0*text_height;
    hmove( pos[0], pos[1], pos[2] );
    hcharstr( "Authors:" );
    pos[1] -= 1.5*text_height;
    hmove( pos[0], pos[1], pos[2] );
    hcharstr( "Don Dovey" );
    pos[1] -= 1.5*text_height;
    hmove( pos[0], pos[1], pos[2] );
    hcharstr( "Tom Spelce" );
    pos[1] -= 1.5*text_height;
    hmove( pos[0], pos[1], pos[2] );
    hcharstr( "Doug Speck" );

    pos[1] -= 3.0*text_height;
    hmove( pos[0], pos[1], pos[2] );
    text_height = 14.0 * vp_to_world_y;
    htextsize( text_height, text_height );
    hcharstr(
        "Copyright (c) 1992-2009 Lawrence Livermore National Laboratory");
    hcentertext( FALSE );

    antialias_lines( FALSE, 0 );
    glLineWidth( 1.0 );

    gui_swap_buffers();
}


/*
 * SECTION_TAG( Swatch )
 */


/*****************************************************************
 * TAG( init_swatch )
 *
 * Init the swatch rendering context.
 */
void
init_swatch( void )
{
    glClearColor( 0.0, 0.0, 0.0, 0.0 );
    glClear( GL_COLOR_BUFFER_BIT );
    glFlush();
    gui_swap_buffers();
}


/*****************************************************************
 * TAG( draw_mtl_swatch )
 *
 * Draw the color swatch with a black border.
 */
void
draw_mtl_swatch( float mtl_color[3] )
{
    glClear( GL_COLOR_BUFFER_BIT );
    glColor3f( (GLfloat) mtl_color[0], (GLfloat) mtl_color[1],
               (GLfloat) mtl_color[2] );
    gluOrtho2D( -1.0, 1.0, -1.0, 1.0 );
    glBegin( GL_POLYGON );
    glVertex2f( -0.75, -0.75 );
    glVertex2f( -0.75, 0.75 );
    glVertex2f( 0.75, 0.75 );
    glVertex2f( 0.75, -0.75 );
    glEnd();
    glFlush();

    gui_swap_buffers();
}


/*
 * SECTION_TAG( Miscellaneous )
 */


/************************************************************
 * TAG( get_verts_of_bbox )
 *
 * Return the eight vertices of a bounding box from the
 * min and max vertices.
 */
static void
get_verts_of_bbox( float bbox[2][3], float verts[8][3] )
{
    int i;

    for ( i = 0; i < 8; i++ )
    {
        VEC_COPY( verts[i], bbox[0] );
    }
    verts[1][0] = bbox[1][0];
    verts[2][0] = bbox[1][0];
    verts[2][1] = bbox[1][1];
    verts[3][1] = bbox[1][1];
    verts[4][2] = bbox[1][2];
    verts[5][0] = bbox[1][0];
    verts[5][2] = bbox[1][2];
    VEC_COPY( verts[6], bbox[1] );
    verts[7][1] = bbox[1][1];
    verts[7][2] = bbox[1][2];
}


/************************************************************
 * TAG( hvec_copy )
 *
 * Copy a homogenous vector.
 */
static void
hvec_copy( float b[4], float a[4] )
{
    b[0] = a[0];
    b[1] = a[1];
    b[2] = a[2];
    b[3] = a[3];
}


/************************************************************
 * TAG( antialias_lines )
 *
 * Switch on or off line antialiasing.
 */
static void
antialias_lines( Bool_type select, Bool_type depth_buffer_lines )
{
    if ( select == TRUE )
    {
        /* Antialias the lines. */
        glEnable( GL_LINE_SMOOTH );
        glEnable( GL_BLEND );
        glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
        /*
         * For correct antialiasing, this should be uncommented
         * AND line segments should be depth-sorted and rendered
         * back-to-front (but there's no code for this yet).
         * Disabling depth masking (i.e., making the Z-buffer
         * read only instead of read/write as GRIZ did previously)
         * without depth-sorting will allow incorrect overlaps
         * between lines.
         */
        if ( !depth_buffer_lines )
            glDepthMask( GL_FALSE );
    }
    else
    {
        /* Antialiasing off. */
        glDisable( GL_LINE_SMOOTH );
        glDisable( GL_BLEND );
        glDepthMask( GL_TRUE );
    }
}


/************************************************************
 * TAG( begin_draw_poly )
 *
 * Set up the stencil buffer for hidden line drawing, and
 * initialize some variables needed for polygon projection
 * when doing "good" interpolation.
 */
static void
begin_draw_poly( Analysis *analy )
{
    float xp, yp, cp;
    float aspect;

    /* Cache the view matrix for later use. */
    view_transf_mat( &cur_view_mat );

    /* Stencil buffer setup -- not used for 2D grids. */
    if ( (analy->mesh_view_mode == RENDER_HIDDEN || analy->mesh_view_mode == RENDER_WIREFRAME ||
            analy->mesh_view_mode == RENDER_WIREFRAMETRANS )
            && analy->dimension == 3 )
    {
        glEnable( GL_STENCIL_TEST );
        glStencilFunc( GL_ALWAYS, 0, 1 );
        glStencilOp( GL_INVERT, GL_INVERT, GL_INVERT );
        glColor3fv( v_win->mesh_color  );
        glLineWidth( (GLfloat) analy->hidden_line_width );
    }

    /* Calculate constants needed for projecting polygons to screen. */
    if ( analy->interp_mode == GOOD_INTERP )
    {
        /* Correct the aspect ratio for NTSC video, if requested. */
        aspect = v_win->aspect_correct * v_win->vp_width /
                 (float)v_win->vp_height;

        /* Get the window dimensions at the near plane. */
        if ( v_win->orthographic )
            cp = 1.0;
        else
            cp = v_win->near * TAN( DEG_TO_RAD( 0.5*v_win->cam_angle ) );

        if ( aspect >= 1.0 )
        {
            xp = cp * aspect;
            yp = cp;
        }
        else
        {
            xp = cp;
            yp = cp / aspect;
        }

        /* Projection parameters. */
        if ( v_win->orthographic )
        {
            /*
             * See glOrtho() and glFrustrum() in the OpenGL Reference
             * Manual to see where these factors come from.
             */
            proj_param_x = v_win->vp_width / ( xp * 2.0 );
            proj_param_y = v_win->vp_height / ( yp * 2.0 );
        }
        else
        {
            proj_param_x = v_win->vp_width * v_win->near / ( xp * 2.0 );
            proj_param_y = v_win->vp_height * v_win->near / ( yp * 2.0 );
        }
    }
}


/************************************************************
 * TAG( end_draw_poly )
 *
 * Turn off the stencil buffer to end hidden line drawing.
 */
static void
end_draw_poly( Analysis *analy )
{
    if ( analy->mesh_view_mode == RENDER_HIDDEN || analy->mesh_view_mode == RENDER_WIREFRAME ||
            analy->mesh_view_mode == RENDER_WIREFRAMETRANS )
    {
        glDisable( GL_STENCIL_TEST );
        glLineWidth( (GLfloat) 1.0 );
    }
}


/************************************************************
 * TAG( set_particle_radius )
 *
 * Set the radius to use for rendering particle class data.
 */
extern void
set_particle_radius( float radius )
{
    if ( radius > 0.0 )
        particle_radius = radius;
}


/************************************************************
 * TAG( linear_variable_scale )
 *
 * Establish a linear scale, i.e., a scale having an interval size
 * which is a product of an integer power of 10 and 1, 2, or 5,
 * and scale values which are integer multiples of the interval
 * size.
 *
 * Given "data_minimum", "data_maximum", and "qty_of_intervals"
 * the routine computes a new scale range from "scale_minimum" to
 * "scale_maximum" divisible into approximately "qty_of_intervals"
 * linear intervals of size "distance".
 *
 */
void
linear_variable_scale( float data_minimum, float data_maximum,
                       int qty_of_intervals, float *scale_minimum,
                       float *scale_maximum, float *distance,
                       int *error_status )
{
    float
    approximate_interval
    ,float_magnitude
    ,geometric_mean [3]
    ,scaled_interval;

    int
    evaluate
    ,i
    ,log_approximate_interval
    ,magnitude;

    /* Establish compensation for computer round-off */
    static float epsilon = 0.00002;

    /*
     * Establish acceptable scale interval values (times an integer power
     * of 10
     */
    static float interval_size[4] = { 1.0, 2.0, 5.0, 10.0 };

    if ( (data_maximum > data_minimum) && (qty_of_intervals > 0) )
    {
        *error_status = FALSE;

        /* Establish geometric means between adjacent scale interval values */
        geometric_mean[0] = sqrt( (double) interval_size[0] *
                                  (double) interval_size[1] );
        geometric_mean[1] = sqrt( (double) interval_size[1] *
                                  (double) interval_size[2] );
        geometric_mean[2] = sqrt( (double) interval_size[2] *
                                  (double) interval_size[3] );

        /* Compute approximate interval size of data */
        approximate_interval = ( data_maximum - data_minimum )
                               / qty_of_intervals;
        log_approximate_interval = (int) log10( approximate_interval );
        if ( approximate_interval < 1.0 )
            log_approximate_interval--;

        /* Scale approximate interval to fall within range:  1 -- 10 */
        scaled_interval = approximate_interval /
                          pow( 10.0, (double) log_approximate_interval );

        /* Establish closest permissible value for scaled interval */
        evaluate = TRUE;
        i = 0;
        while ( (TRUE == evaluate) && (i < 3) )
        {
            if ( scaled_interval < geometric_mean[i] )
                evaluate = FALSE;
            else
                i++;
        }

        *distance = interval_size[i]
                    * pow( 10.0, (double) log_approximate_interval );

        /* Compute new scale minimum */
        float_magnitude = data_minimum / *distance;
        magnitude = (int) float_magnitude;

        if ( float_magnitude < 0.0 )
            magnitude--;

        if ( fabs( (float) magnitude + 1.0 - float_magnitude ) < epsilon )
            magnitude++;

        *scale_minimum = *distance * (float) magnitude;


        /* Compute new scale maximum */
        float_magnitude = data_maximum / *distance;
        magnitude       = (int) ( float_magnitude + 1.0 );

        if ( float_magnitude < -1.0 )
            magnitude--;

        if ( fabs( float_magnitude + 1.0 - (float) magnitude) < epsilon )
            magnitude--;

        *scale_maximum = *distance * magnitude;

        /* Adjust, as required, scale limits to account for round-off */
        *scale_minimum = MIN( *scale_minimum, data_minimum );
        *scale_maximum = MAX( *scale_maximum, data_maximum );
    }
    else
        *error_status = TRUE;
}

/************************************************************
 * TAG( log_scale_data_shift )
 *
 * This function will compute a log scale of n equal distance
 * points.
 *
 *
 */
void
log_scale_data_shift( float val, float data_minimum, float data_maximum,
                      float *new_val,  float *new_data_minimum,
                      float *new_data_maximum, float *data_shift,
                      float *data_mult)

{
    double data_shift_dbl, mult_factor = 1.0;

    /* Shift data to be in positive range */
    *data_shift = 0.;

    if (data_minimum==0.)
    {
        data_shift_dbl = 1.0;
        data_minimum   = 1.0;
    }

    if (data_minimum<0.)
    {
        data_shift_dbl = ((double)-data_minimum)+1;
        data_minimum   = 1.0;
    }

    if (data_minimum<1.0)
    {
        data_minimum = 1.0;
        mult_factor  = 1/data_minimum;
    }

    data_maximum = (data_maximum+data_shift_dbl)*mult_factor;

    val = (val+data_shift_dbl)*mult_factor;
    if (val<=0.)
    {
        val = 1.0;
    }

    *data_shift       = data_shift_dbl;
    *data_mult        = mult_factor;
    *new_data_minimum = data_minimum;
    *new_data_maximum = data_maximum;
    *new_val          = val;
}


/************************************************************
 * TAG( log_variable_scale )
 *
 * This function will compute a log scale of n equal distance
 * points.
 *
 * Given "data_minimum", "data_maximum", and "qty_of_intervals"
 * the routine computes a new scale range from "scale_minimum" to
 * "scale_maximum" divisible into approximately "qty_of_intervals"
 * log intervals of size "distance".
 *
 */
void
log_variable_scale( float data_minimum, float data_maximum,
                    int qty_of_intervals, float *scale_minimum,
                    float *scale_maximum, float *distance,
                    float *data_shift, int *error_status)
{

    Bool_type logFound = FALSE;

    float data_range;

    double a,
           al,
           b,
           dist, distl,
           fm1,  fm2,
           xmin, xminl, xminp,
           xmax, xmaxl, xmaxp;

    float vint[11] = { 10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.5 };

    int   fn,
          i,
          index,
          m1, m2,
          nal,
          magnitude,
          n,
          np, nx;


    float val, new_val, new_data_minimum, new_data_maximum,
          new_data_shift, new_data_mult;

    /* Establish compensation for computer round-off */
    static float epsilon = 0.00002;
    /*
     * Establish acceptable scale interval values (times an integer power
     * of 10
     */
    static float interval_size[4] = { 1.0, 2.0, 5.0, 10.0 };

    /* Shift data to be in a positive range */
    val = 0.0;
    log_scale_data_shift( val,      data_minimum, data_maximum,
                          &new_val, &new_data_minimum,
                          &new_data_maximum, &new_data_shift,
                          &new_data_mult);

    *data_shift = new_data_shift;

    if ( (data_minimum>0.) &&
            (data_maximum > data_minimum) && (qty_of_intervals > 0) )
    {
        xmin  = new_data_minimum;
        xmax  = new_data_maximum;

        n     = qty_of_intervals;
        xminl = log10(xmin);
        xmaxl = log10(xmax);
        fn    = n;

        /* Find an interval size */

        a   = (xmaxl - xminl)/fn;
        al  = log10(a);
        nal = al;

        if (a<1.0) nal--;
        b = a/pow(10.0, (double)nal);

        index = -1;
        for (i=0;
                i<9;
                i++)
        {
            if((b==(10.0/vint[i])) ||
                    b<((10.0/vint[i])+epsilon))
            {
                index = i;
                break;
            }
        }

        if (index<0)
            index = 9;

        while (!logFound)
        {
            distl = pow(10.0,(double)(nal+1.0))/vint[index];
            fm1   = xminl/distl;
            m1    = fm1;
            if (fm1<0.) m1--;
            if (abs((float)m1+1.0-fm1)<epsilon) m1++;

            xminp = distl*(float)m1;
            fm2   = xmaxl/distl;
            m2    = fm2 + 1.0;

            if (fm2<(-1.0)) m2--;
            if (abs(fm2+1.-(float)m2)<epsilon) m2--;
            xmaxp = distl*(float)m2;

            np = m2 - m1;
            if (np<=n)
            {
                logFound = TRUE;
                break;
            }
            else
                index++;
        }

        nx = (n-np)/2.0;

        xminp = xminp - (float)nx*distl;
        xmaxp = xminp + (float)n*distl;

        *distance = pow(10.0, distl);
        xminp     = pow(10.0, xminp);
        xmaxp     = pow(10.0, xmaxp);

        if (xminp>xmin) xminp = xmin;
        if (xmaxp<xmax) xmaxp = xmax;

        *scale_minimum = xminp;
        *scale_maximum = xmaxp;

        *error_status = FALSE;
    }
    else
    {
        *distance      = 0.0;
        *scale_minimum = 0.0;
        *scale_maximum = 0.0;

        *error_status = TRUE;
    }
}


/************************************************************
 * TAG( round )
 *
 * Perform "fuzzy" rounding
 *
 */
double
griz_round( double x, double comparison_tolerance )
{
    return( tfloor( (x + 0.5), comparison_tolerance ) );
}


/************************************************************
 * TAG( tfloor )
 *
 * Perform "fuzzy" floor operation
 *
 */
#define DFLOOR(a)  (DINT( (a) ) - (double)fmod( (2.0 + DSIGN( 1.0, (a) )), 3.0 ))
#define DINT(a)    ((a) - (double)fmod( (a), 1.0 ))
#define DSIGN(a,b) ((b) >= 0.0 ? (double)fabs( (a) ) : -(double)fabs( (a) ))

static double
tfloor( double x, double comparison_tolerance )
{
    double q, rmax, eps5;
    double temp_tfloor;


    if ( x >= 0.0 )
        q = 1.0;
    else
        q = 1.0 - comparison_tolerance;

    rmax = q / (2.0 - comparison_tolerance);

    eps5 = comparison_tolerance / q;

    temp_tfloor = DFLOOR( x + MAX( comparison_tolerance,
                                   MIN( rmax, eps5 * (double)fabs( 1.0 + DFLOOR( x ) ))));

    if ( (x <= 0.0) || ((temp_tfloor - x) < rmax) )
        return( temp_tfloor );
    else
        return( temp_tfloor - 1.0 );
}


/************************************************************
 * TAG( machine_tolerance )
 *
 * Compute machine tolerance
 *
 */
double
machine_tolerance( void )
{
    double a;
    double b;
    double c;
    double tolerance;


    a = 4.0 / 3.0;

label:

    b = a - 1.0;

    c = b + b + b;

    tolerance = (double)fabs( c - 1.0);

    if ( 0.0 == tolerance )
        goto label;

    return( tolerance );
}


/************************************************************
 * TAG( check_interp_mode )
 *
 * To be called when isocontours or isosurfaces are turned
 * on to notify users of inconsistency with non-interpolated
 * results.
 *
 */
void
check_interp_mode( Analysis *analy )
{
    if ( analy->interp_mode == NO_INTERP
            && analy->cur_result->origin.is_elem_result )
        popup_dialog( INFO_POPUP, "%s\n%s\n%s\n%s",
                      "You are currently visualizing non-interpolated results.",
                      "Isocontours/Isosurfaces will be generated from results",
                      "which have been interpolated to the nodes to provide",
                      "variation over the elements." );
}


/************************************************************
 * TAG( draw_locref )
 *
 * This routine draws the local coordinate frame for
 * selected shell elements.
 *
 */
static void
draw_locref( Analysis *analy )
{
    int j;
    float pti[3], ptj[3], ptk[3], cvec[3];
    MO_class_data *p_mo_class;
    Specified_obj *p_so;
    float leng1, leng2, cx, cy, zpos;
    float x[4], y[4], z[4];
    float xn[3], yn[3], zn[3];
    GVec3D *coords;
    int (*connects)[4];
    float vp_to_world[2];

    /* Hard-code nice bright blue for now. */
    glColor3f( 0.2, 0.2, 0.9 );

    /* Calculate the size of the local coordinate system. */
    get_foreground_window( &zpos, &cx, &cy );
    vp_to_world[0] = 2.0*cx / v_win->vp_width;
    vp_to_world[1] = 2.0*cy / v_win->vp_height;
    leng1 = 50*vp_to_world[0];
    leng2 = 20*vp_to_world[0];

    for ( p_so = analy->selected_objects; p_so != NULL; NEXT( p_so ) )
    {
        p_mo_class = p_so->mo_class;

        if ( p_mo_class->superclass == G_QUAD )
        {
            connects = (int (*)[4]) p_mo_class->objects.elems->nodes;
            coords = analy->state_p->nodes.nodes3d;

            /* Get the quad element geometry. */
            for ( j = 0; j < 4; j++ )
            {
                x[j] = coords[ connects[p_so->ident][j] ][0];
                y[j] = coords[ connects[p_so->ident][j] ][1];
                z[j] = coords[ connects[p_so->ident][j] ][2];
            }

            /* Compute the center point of the element */
            cvec[0] = .5 * (x[0] + x[2]);
            cvec[1] = .5 * (y[0] + y[2]);
            cvec[2] = .5 * (z[0] + z[2]);

            /*Compute the local coordinate matrix.*/

            /* Get average vectors for the sides of the element. */

            xn[0] = -x[0] + x[1] + x[2] - x[3];
            xn[1] = -y[0] + y[1] + y[2] - y[3];
            xn[2] = -z[0] + z[1] + z[2] - z[3];
            yn[0] = -x[0] - x[1] + x[2] + x[3];
            yn[1] = -y[0] - y[1] + y[2] + y[3];
            yn[2] = -z[0] - z[1] + z[2] + z[3];

            /* Get normal to element (assumes element is planar). */
            VEC_CROSS( zn, xn, yn );

            /* An orthonormal basis is what we need. */
            vec_norm( zn );
            vec_norm( yn );
            VEC_CROSS( xn, yn, zn);

            /* Scale the axes and translate to element center. */
            for ( j = 0; j < 3; j++ )
            {
                pti[j] = xn[j] * leng1 + cvec[j];
                ptj[j] = yn[j] * leng1 + cvec[j];
                ptk[j] = zn[j] * leng2 + cvec[j];
            }

            /*Draw the local coordinates and label the x-axis.*/

            glDepthFunc( GL_ALWAYS );

            antialias_lines( TRUE, FALSE );

            glBegin( GL_LINES );
            glVertex3fv( cvec );
            glVertex3fv( pti );
            glVertex3fv( cvec );
            glVertex3fv( ptj );
            glVertex3fv( cvec );
            glVertex3fv( ptk );
            glEnd();

            draw_3d_text( pti, "x", FALSE);

            antialias_lines( FALSE, FALSE );

            glDepthFunc( GL_LEQUAL );
        }
    }
}

/************************************************************
 * TAG( draw_locref_hex )
 *
 * This routine draws the local coordinate frame for
 * selected shell elements.
 *
 */
static void
draw_locref_hex( Analysis *analy )
{
    int j;
    float vi[3], vj[3], vk[3], vd[3];
    float len;
    Transf_mat tmat;
    float pt[3], pto[3], pti[3], ptj[3], ptk[3];
    MO_class_data *p_mo_class;
    Specified_obj *p_so;
    float leng1, cx, cy, zpos;
    float x[8], y[8], z[8];
    GVec3D *coords;
    int (*connects)[8];
    float world_to_vp[2], vp_to_world[2];
    int n1, n2, n3;

    n1=0;
    n2=1;
    n3=3;


    /* Hard-code nice dodger blue for now. */
    glColor3f( 0.25, 0.25, 0.9 );

    /* Calculate the size of the local coordinate system. */
    get_foreground_window( &zpos, &cx, &cy );
    world_to_vp[0] = v_win->vp_width / (2.0*cx);
    world_to_vp[1] = v_win->vp_height / (2.0*cy);
    vp_to_world[0] = 2.0*cx / v_win->vp_width;
    vp_to_world[1] = 2.0*cy / v_win->vp_height;
    leng1 = 35*vp_to_world[0];

    for ( p_so = analy->selected_objects; p_so != NULL; NEXT( p_so ) )
    {
        p_mo_class = p_so->mo_class;

        if ( p_mo_class->superclass == G_HEX )
        {
            connects = (int (*)[8]) p_mo_class->objects.elems->nodes;
            coords = analy->state_p->nodes.nodes3d;

            /* Get the hex element geometry. */
            for ( j = 0; j < 8; j++ )
            {
                x[j] = coords[ connects[p_so->ident][j] ][0];
                y[j] = coords[ connects[p_so->ident][j] ][1];
                z[j] = coords[ connects[p_so->ident][j] ][2];
            }

            /* Compute the center point of the element */
            pto[0] = 0.0;
            pto[1] = 0.0;
            pto[2] = 0.0;

            for ( j = 0; j < 8; j++ )
            {
                pto[0] += .125*x[j];
                pto[1] += .125*y[j];
                pto[2] += .125*z[j];
            }

            /*Compute the local coordinate matrix.*/
            /*
             *  "i" axis is the direction from node 1 to node 2
             */
            vi[0] = x[n2] - x[n1];
            vi[1] = y[n2] - y[n1];
            vi[2] = z[n2] - z[n1];
            len = VEC_LENGTH( vi );

            /*
             *  Vector "D" is the direction from node 1 to node 4
             */
            vd[0] = x[n3] - x[n1];
            vd[1] = y[n3] - y[n1];
            vd[2] = z[n3] - z[n1];

            /*
             * Normalize the reference vectors
             */
            vec_norm(vi);
            vec_norm(vd);

            /*
             *  "k" axis is cross-product of "i" and vector "D"
             */
            VEC_CROSS( vk, vi, vd );

            /*
             *  "j" axis is cross-product of "k" and "i"
             */
            VEC_CROSS( vj, vk, vi );

            for ( j = 0; j < 3; j++ )
            {
                tmat.mat[0][j] = vi[j];
                tmat.mat[1][j] = vj[j];
                tmat.mat[2][j] = vk[j];
                tmat.mat[3][j] = pto[j];
            }

            /* Draw the axes. */
            VEC_SET( pt, 0.0, 0.0, 0.0 );
            point_transform( pto, pt, &tmat );

            VEC_SET( pt, leng1, 0.0, 0.0 );
            point_transform( pti, pt, &tmat );

            VEC_SET( pt, 0.0, leng1, 0.0 );
            point_transform( ptj, pt, &tmat );

            VEC_SET( pt, 0.0, 0.0, leng1 );
            point_transform( ptk, pt, &tmat );


            /*Draw the local coordinates and label the x-axis.*/

            glDepthFunc( GL_ALWAYS );

            antialias_lines( TRUE, FALSE);

            glColor3fv( v_win->foregrnd_color );

            glBegin( GL_LINES );
            glVertex3fv( pto );
            glVertex3fv( pti );
            glVertex3fv( pto );
            glVertex3fv( ptj );
            glVertex3fv( pto );
            glVertex3fv( ptk );
            glEnd();

            draw_3d_text( pti, "x", FALSE);
            draw_3d_text( ptj, "y", FALSE);
            draw_3d_text( ptk, "z", FALSE);

            antialias_lines( FALSE, FALSE );

            glDepthFunc( GL_LEQUAL );
        }
    }
}


/************************************************************
 * TAG( get_free_node_result )
 *
 * This function will return the value of a ml for
 * a specified node number.
 */
float
get_free_node_result( Analysis  *analy, MO_class_data *p_mo_class, int particle_num, Bool_type *result_defined,
                      Bool_type *valid_free_node )
{
    int   node_num=0, (*particle_nodes)[8];
    float val=0.0, *nodal_data;

    Bool_type nodal_result=TRUE;

    *valid_free_node = FALSE;

    particle_nodes = p_mo_class->objects.elems->nodes;
    nodal_data     = NODAL_RESULT_BUFFER( analy );

    if ( !nodal_data )
    {
        nodal_data = p_mo_class->data_buffer; /* Do not get the particle data from the nodal buffer */
        nodal_result = FALSE;
    }

    if ( nodal_data!=NULL )
        *result_defined = TRUE;

    node_num = particle_nodes[particle_num][0];
    *valid_free_node = TRUE;
    if ( is_particle_class( analy, p_mo_class->superclass, p_mo_class->short_name ) )
        val = p_mo_class->data_buffer[particle_num];
    else if ( nodal_data )
        val = nodal_data[node_num];

    return ( val );
}


/**************************************************************
 * TAG( get_ml_result )
 *
 * This function will return the value of a ml for a specified
 * a specified node number.
 */
float
get_ml_result( Analysis  *analy, MO_class_data *p_mo_class, int elem_num, Bool_type *result_defined )
{
    float val=0., num_nodes, *nodal_data, *activity;
    int   node_num;
    MO_class_data *p_node_class;
    Bool_type nodal_result=TRUE;
    float ** sand_arrays;

    int i,j;
    int (*connects_hex)[8], (*connects_particle)[1];


    *result_defined  = FALSE;

    if(analy->cur_result == NULL)
    {
        return(0.0);
    }
  
    if(analy->cur_result != NULL)
    {
        *result_defined = TRUE;
         return p_mo_class->data_buffer[elem_num];
    }

    p_node_class = MESH_P( analy )->node_geom;

    num_nodes = p_node_class->qty;
    if ( elem_num<0 || elem_num>p_mo_class->qty )
        return( 0.0 );

    nodal_data = NODAL_RESULT_BUFFER( analy );

    if ( !nodal_data )
    {
        nodal_data = p_mo_class->data_buffer; /* Do not get the particle data from the nodal buffer */
        nodal_result = FALSE;
    }

    if ( nodal_data!=NULL )
        *result_defined = TRUE;

    if ( p_mo_class->superclass==G_HEX )
    {
        connects_hex = p_mo_class->objects.elems->nodes;
        node_num     = connects_hex[elem_num][0];
    }
    else
    {
        connects_particle = p_mo_class->objects.elems->nodes;
        node_num          = connects_particle[elem_num][0];
    }

    if ( !nodal_result )
    {
        val = p_mo_class->data_buffer[elem_num];
        *result_defined = TRUE;
    }
    else if ( nodal_data )
    {
        val = nodal_data[node_num];
        *result_defined = TRUE;
    }
    /*if(analy->show_deleted_elements || analy->show_only_deleted_elements)
    {
        val = p_mo_class->data_buffer[elem_num];
        *result_defined = TRUE;
    }*/
    /*if(analy->state_p->sand_present)
    {
        sand_arrays = analy->state_p->elem_class_sand;
        if(sand_arrays[p_mo_class->elem_class_index] != NULL)
        {
            activity = sand_arrays[p_mo_class->elem_class_index];
        } else
        {
            activity = NULL;
            *result_defined = TRUE;
            return val;
        }

        if(analy->show_deleted_elements)
        {
            val = activity[elem_num];
            *result_defined = TRUE;
        } else if(analy->show_only_deleted_elements)
        {
            if(activity[elem_num] == 0.0)
            {
                val = 0.0;
                *result_defined = TRUE;
            } else
            {
                *result_defined = TRUE;
            }
             
        }
    } */ 
    return ( val );
}


/**************************************************************
 * TAG( get_particle_node_num )
 *
 * This function will return the first node number for a particle
 * element.
 */
int
get_particle_node_num( Analysis *analy, MO_class_data *p_mo_class,
                       int elem_num )
{
    int  node_num=0, (*connects)[8];

    connects =  (int(*)[8]) p_mo_class->objects.elems->nodes;

    if ( elem_num<0 || elem_num>p_mo_class->qty )
        return( -1 );

    node_num = connects[elem_num][0];

    return ( node_num );
}


/************************************************************
 * TAG( check_for_free_nodes )
 *
 * This function will check if this database has either free-
 * nodes or free-particles. The appropriate status flag will
 * set in the Analysis data structure if either are found.
 */
void
check_for_free_nodes( Analysis *analy )
{
    MO_class_data *p_element_class,
                  *p_mo_class,
                  **mo_classes;

    Mesh_data     *p_mesh;

    List_head     *p_lh;

    int i, j;

    int status_mass=OK, status_vol=OK;

    float *free_nodes_mass=NULL, *free_nodes_vol=NULL;

    analy->free_nodes_found     = FALSE;
    analy->particle_nodes_found = FALSE;

    p_mesh       = MESH_P( analy );

    /* First check to see if we have free nodes */

    status_mass = mili_db_get_param_array(analy->db_ident,     "Nodal Mass", (void *) &free_nodes_mass);
    status_vol  = mili_db_get_param_array(analy->db_ident,   "Nodal Volume", (void *) &free_nodes_vol);

    if ( status_mass )
    {
        if ( free_nodes_mass )
            free ( free_nodes_mass );
    }
    if ( status_vol )
    {
        if ( free_nodes_vol )
            free ( free_nodes_vol );
    }

    if ( free_nodes_mass && free_nodes_vol )
        analy->free_nodes_found = TRUE;

    /* Loop over each element superclass. */
    for ( i = 0; i < QTY_SCLASS; i++ )
    {
        if ( !IS_ELEMENT_SCLASS( i ) )
            continue;

        p_lh = p_mesh->classes_by_sclass + i;

        /* If classes exist from current superclass... */
        if ( p_lh->qty > 0 )
        {
            mo_classes = (MO_class_data **) p_mesh->classes_by_sclass[i].list;

            /* Loop over each class. */
            for ( j = 0;
                    j < p_lh->qty;
                    j++ )
            {
                p_mo_class = mo_classes[j];

                if ( is_particle_class(analy, p_mo_class->superclass, p_mo_class->short_name) )
                {
                    analy->particle_nodes_found = TRUE;
                    break;
                }
            }
        }
    }

    /* Until we fully test both free-nodes and particle-nodes in
     * same problem do not allow both concurrently.
     */
    if ( analy->particle_nodes_found )
        analy->free_nodes_found = FALSE;
}

/************************************************************
 * TAG( draw_free_nodes )
 *
 * This function will display all free nodes as spheres with
 * a specified scale factor.
 */
static void
draw_free_nodes( Analysis *analy )
{
    MO_class_data   *p_element_class=NULL,
                     *p_node_class=NULL,
                      *p_mo_class=NULL,
                       *p_ml_class=NULL,
                        **mo_classes;

    Mesh_data     *p_mesh;

    List_head     *p_lh;

    Visibility_data *p_vd;
    unsigned char *part_visib;

    char *cname;

    int   *free_nodes_list=NULL, *free_nodes_elem_list=NULL;
    int   *part_nodes_list=NULL;
    int   *part_nodes_result=NULL;

    int   nd;

    float *activity, activity_flag=0.0, temp_activity_flag=0.0;
    float *data_array;
    float col[4], *hilite_col;
    float verts[3], leng[3];
    float rfac, node_base_radius, node_radius;
    float **sand_arrays;
    int   *sph_type;
    float *free_nodes_mass=NULL, *free_nodes_vol=NULL;
    float  mass_scale_factor, mass_scale_factor_max = 0, mass_scale_factor_min = MAXFLOAT;

    float rmin, rmax;
    float poly_pts[4][3],
    poly_norms[6][3] = {{-1., 0., 0},  {0., 1., 0.}, {1., 0., 0.},
        {0., -1., 0.}, {0., 0., 1.}, {0., 0., -1}
    };

    int  elem_qty;
    int  fn_count = 0;
    int  node_index, node_num, num_nodes, node_qty;
    int  class_index, class_qty;
    int  conn_qty;
    int  vert_cnt;

    int i, j, k, l;
    int dim, obj_qty;

    int  *connects;

    int  free_nodes_found = FALSE;
    int  part_nodes_found = FALSE;
    int  dbc_nodes_found  = FALSE;
    int  mass_scaling     = FALSE;
    int  vol_scaling      = FALSE;
    int  free_node_data_index=0;
    Bool_type skip_node=False;

    int  status;

    float *nodal_data;

    Bool particle_class       = FALSE, particle_class_found = FALSE,
         showgs=FALSE,
         particle_node_found  = FALSE, free_node_found=FALSE;
    int  particle_count;
    int  *particle_nodes;
    int  particle_hide=0;

    /* Variables related to materials */
    unsigned char *disable_mtl,  *hide_mtl, hide_one_mat;
    unsigned char *disable_part, *hide_part;

    int  *mat, mat_num, elem_id, ode_cnt = 0;

    Bool_type result_defined=FALSE;

    /* Used for fast sphere drawing */
    GLUquadricObj *sphere;
    GLuint display_list;
    Bool_type fastSpheres=TRUE;
    GLuint texMaps[1];
    static char * texNames[1] = { "maps/earth-seas-small.rgb"};
    IMAGE *img;
    GLint gluerr;

    Refl_plane_obj *plane;
    float refl_verts[3];

    float val=0.;
    Bool_type show_result=TRUE;

    p_mesh      = MESH_P( analy );

    disable_part = p_mesh->disable_particle;
    hide_mtl     = p_mesh->hide_material;
    hide_part    = p_mesh->hide_particle;

    rfac       = analy->free_nodes_scale_factor; /* Radius scale factor. */
    dim        = analy->dimension;
    data_array = NODAL_RESULT_BUFFER( analy );

    get_min_max( analy, FALSE, &rmin, &rmax );

    p_node_class = MESH_P( analy )->node_geom;
    num_nodes    = p_node_class->qty;

    /* Free node list is an array of node length used to identify
     * the free nodes - free_nodes_list[i] is set to mat# if node is free.
     */
    free_nodes_list      = NEW_N( int,   num_nodes, "free_nodes_list" );
    free_nodes_elem_list = NEW_N( int,   num_nodes, "free_nodes_elem_list" );
    part_nodes_list      = NEW_N( int,   num_nodes, "part_nodes_list" );
    part_nodes_result    = NEW_N( Bool_type, num_nodes, "part_nodes_result" );

    if ( analy->free_nodes_list==NULL && analy->free_nodes_vals==NULL )
    {
        analy->free_nodes_list = NEW_N( Bool_type, num_nodes, "free_nodes_list_saved" );
        analy->free_nodes_vals = NEW_N( float,     num_nodes, "free_nodes_vals" );
    }

    for (i=0;
            i<num_nodes;
            i++)
    {
        free_nodes_list[i]      = 0;
        free_nodes_elem_list[i] = 0;
        part_nodes_list[i]      = -1;
        part_nodes_result[i]    = TRUE;   /* Set to FALSE if current result is not
					    * valid for this node.
					    */
        analy->free_nodes_list[i] = FALSE;
        analy->free_nodes_vals[i] = 0.0;
    }

    /* If mass scaling option is enabled, then look for the nodal masses */

    if (analy->free_nodes_mass_scaling)
    {
        status = mili_db_get_param_array(analy->db_ident,     "Nodal Mass",   (void *) &free_nodes_mass);
        if (status==0)
        {
            mass_scaling = TRUE;
            vol_scaling  = FALSE;
            status = mili_db_get_param_array(analy->db_ident,  "Nodal Volume", (void *) &free_nodes_vol);
        }
    }

    /* Read the nodal volumes if volume scaling is enabled */

    if (analy->free_nodes_vol_scaling)
    {
        status = mili_db_get_param_array(analy->db_ident,  "Nodal Volume", (void *) &free_nodes_vol);
        if (status==0)
        {
            vol_scaling  = TRUE;
            mass_scaling = FALSE; /* Vol scaling will override mass scaling */
        }
    }

    /* Set up for polygon drawing. */

    if ( analy->free_nodes_sphere_res_factor<=2 && ! fastSpheres )
    {
        begin_draw_poly( analy );
    }

    /* Loop over each element superclass. */
    for ( i = 0; i < QTY_SCLASS; i++ )
    {
        if ( !IS_ELEMENT_SCLASS( i ) )
            continue;

        p_lh = p_mesh->classes_by_sclass + i;

        /* If classes exist from current superclass... */
        if ( p_lh->qty > 0 )
        {
            mo_classes = (MO_class_data **) p_mesh->classes_by_sclass[i].list;

            sand_arrays = analy->state_p->elem_class_sand;
            sph_type    = analy->state_p->sph_class_itype;

            node_qty    = qty_connects[i];
            conn_qty    = ( i == G_BEAM ) ? node_qty - 1 : node_qty;

            /* Loop over each class. */
            for ( j = 0;
                    j < p_lh->qty;
                    j++ )
            {
                p_mo_class = mo_classes[j];
                if ( is_particle_class( analy, p_mo_class->superclass, p_mo_class->short_name ) )
                {
                    p_ml_class = p_mo_class;
                    particle_class       = TRUE;
                    particle_class_found = TRUE;
                    data_array     = p_mo_class->data_buffer; /* Do not get the particle data from the nodal buffer */
                    part_visib     = p_mo_class->p_vis_data->visib;

                    if ( analy->pn_nodal_result && analy->pn_ref_nodes[0]==NULL )
                    {
                        if ( p_mo_class->referenced_nodes==NULL)
                        {
                            create_elem_class_node_list( analy, p_mo_class,
                                                         &p_mo_class->referenced_node_qty,
                                                         &p_mo_class->referenced_nodes );
                        }
                        analy->pn_ref_node_count[0] = p_mo_class->referenced_node_qty;
                        analy->pn_ref_nodes[0]      = p_mo_class->referenced_nodes;
                    }
                }
                else
                    particle_class = FALSE;

                connects   = p_mo_class->objects.elems->nodes;
                mat        = p_mo_class->objects.elems->mat;

                if ( analy->state_p->sand_present && sand_arrays[p_mo_class->elem_class_index]!=NULL )
                    activity = sand_arrays[p_mo_class->elem_class_index];
                else
                    activity = NULL;

                if ( particle_class && is_dbc_class( analy, p_mo_class->superclass, p_mo_class->short_name ) )
                {
                    dbc_nodes_found = TRUE;
                    show_result = result_has_class( analy->cur_result, p_mo_class, analy );
                }
                else show_result = TRUE;

                /* Loop over each element */
                for ( k = 0;
                        k < p_mo_class->qty;
                        k++ )
                {
                    if( activity )
                    {
                        activity_flag = activity[k];
                    } else
                    {
                        activity_flag = 1.0;
                    }
                    
                    mat_num = mat[k];

                     if ( !particle_class)
                    {
                        if ( analy->show_deleted_elements )
                            activity_flag = 1.0;

                        if ( analy->show_only_deleted_elements )
                        {
                            if  ( activity_flag == 0.0 )
                            {
                                activity_flag = 1.0;
                            }
                        }
                    } else
                    {
                        if ( analy->show_deleted_elements )
                            activity_flag = 1.0;

                        if ( analy->show_only_deleted_elements )
                        {
                            if  ( activity != NULL && activity[k] == 0.0 )
                            {
                                activity_flag = 1.0;
                                ode_cnt++;
                            } else
                            {
                                activity_flag = 0.0;
                            }
                        }
                    }

                    if ( sph_type && particle_class )
                    {
                        if ( analy->show_sph_ghost==FALSE && sph_type[k]==1 )
                            activity_flag = 0.0;
                    }

                    free_node_found = !hide_by_object_type( p_mo_class, mat_num, k, analy, data_array ) && activity_flag == 0.0 && analy->free_nodes_enabled;
                    if ( particle_class )
                        particle_node_found = !hide_by_object_type( p_mo_class, mat_num, k, analy, data_array ) && activity_flag > 0.0 && analy->particle_nodes_enabled &&
                                              part_visib[k];
                    else particle_node_found = FALSE;

                    if ( !particle_class )
                        if (free_node_found)
                        {
                            for ( l = 0;
                                    l < conn_qty;
                                    l++ )
                            {
                                nd = connects[k * node_qty + l];
                                if (free_nodes_list[nd]<0 )
                                    continue;

                                free_nodes_found         = TRUE;
                                free_nodes_list[nd]      = mat_num;
                                free_nodes_elem_list[nd] = k;;

                                if (mass_scaling)
                                {
                                    mass_scale_factor = free_nodes_mass[nd];
                                    if (mass_scale_factor_max < mass_scale_factor)
                                        mass_scale_factor_max = mass_scale_factor;
                                    if (mass_scale_factor_min > mass_scale_factor)
                                        mass_scale_factor_min = mass_scale_factor;
                                }

                                if (vol_scaling)
                                {
                                    mass_scale_factor = free_nodes_vol[nd];
                                    if (mass_scale_factor < 0.0)
                                        mass_scale_factor*=(-1.0);

                                    if (mass_scale_factor_max < mass_scale_factor)
                                        mass_scale_factor_max = mass_scale_factor;
                                    if (mass_scale_factor_min > mass_scale_factor)
                                        mass_scale_factor_min = mass_scale_factor;
                                }
                            } /* End For on l */
                        }
                        else
                        {
                            for ( l = 0;
                                    l < conn_qty;
                                    l++ )
                            {
                                nd = connects[k * node_qty + l];
                                free_nodes_list[nd] = -1;
                            }
                        }

                    if (particle_node_found)
                        for ( l = 0;
                                l < conn_qty;
                                l++ )
                        {
                            nd = connects[k * node_qty + l];

                            part_nodes_found         = TRUE;
                            part_nodes_list[nd]      = mat_num;
                            free_nodes_elem_list[nd] = k;
                            part_nodes_result[nd]    = show_result;
                        } /* End For on l */
                }
            }
        }
    }

    /* Make a quick exit if no free nodes are to
     * be displayed.
     */

    if (!free_nodes_found && !part_nodes_found)
    {
        return;
    }

    /* If we are viewing a nodal result then map the nodal
      * result onto the particles.
      *
      */
    if ( analy->pn_nodal_result )
    {
        particle_count = analy->pn_ref_node_count[0];
        particle_nodes = analy->pn_ref_nodes[0];
        nodal_data     = analy->pn_node_ptr[0];
    }

    /* Determine a scaling range factor for the Mass/Vol data */

    mass_scale_factor = 0.5;

    if (mass_scaling || vol_scaling)
    {
        mass_scale_factor = 0.5/mass_scale_factor_max;
    }

    for ( i = 0; i < dim; i++ )
        leng[i] = analy->bbox[1][i] - analy->bbox[0][i];

    node_base_radius  = 0.01 * (leng[0] + leng[1] + leng[2]) / 3.0;
    if ( node_base_radius==0 )
        node_base_radius=1.;

    /* Modified: Nov 3, 2006: IRC - Do not scale particles with window */
    /* node_base_radius *= 1.0 / v_win->scale[0]; */

    /*    if ( v_win->lighting )
      glEnable( GL_LIGHTING ); */

    glEnable( GL_COLOR_MATERIAL );

    if ( fastSpheres)
    {
        display_list = glGenLists( 1 );
        sphere       = gluNewQuadric();
        gluQuadricCallback( sphere, (GLenum) GLU_ERROR, particle_error );

        gluQuadricDrawStyle( sphere, (GLenum) GLU_FILL );
        gluQuadricNormals( sphere, (GLenum) GLU_SMOOTH );

        /*******************************/
        /* Future texture mapping code */
        /* glGenTextures(1,&(texMaps[0]));
        gluQuadricTexture(sphere, GL_TRUE);
               glBindTexture(GL_TEXTURE_2D,texMaps[0]);

         img=ImageLoad(texNames[0]);
               glPixelStorei(GL_UNPACK_ALIGNMENT,4);
               gluerr=gluBuild2DMipmaps(GL_TEXTURE_2D, 4, 100, 100,
        		 	  GL_RGBA, GL_UNSIGNED_BYTE,
        		          (GLvoid *)(img->data)); */

        /*******************************/

        glNewList( display_list, GL_COMPILE );
        gluSphere( sphere, node_base_radius*rfac, analy->free_nodes_sphere_res_factor*2, analy->free_nodes_sphere_res_factor );
        glEndList();
    }

    for (node_index=0;
            node_index<num_nodes;
            node_index++)
    {

        skip_node = TRUE;
        if (free_nodes_found && free_nodes_list[node_index] >= 0)
            skip_node = FALSE;
        if (part_nodes_found && part_nodes_list[node_index] >= 0)
            skip_node = FALSE;
        if ( skip_node )
            continue;

        analy->free_nodes_list[node_index] = TRUE;
        if ( free_nodes_list[node_index] < 0 )
            mat_num = free_nodes_list[node_index];
        mat_num = free_nodes_list[node_index];
        if ( mat_num <0 )
            mat_num = part_nodes_list[node_index];
        elem_id = free_nodes_elem_list[node_index];

        /* Colorflag is TRUE if displaying a result, otherwise
         * polygons are drawn in the material color.
         */
        colorflag = TRUE;
        if ( analy->cur_result==NULL || disable_part[mat_num] )
            colorflag = FALSE;
        if ( dbc_nodes_found && part_nodes_result[node_index] == FALSE )
            colorflag = FALSE;
        if ( analy->material_greyscale && disable_part[mat_num] )
            showgs = TRUE;
        else
            showgs = FALSE;

        get_node_vert_2d_3d( node_index, NULL, analy, verts );

        node_radius = node_base_radius;

        /* Scale node to mass and volume */
        if (mass_scaling)
        {
            node_radius = node_radius * free_nodes_mass[node_index] *
                          mass_scale_factor;
        }

        /* Scale node to volume only */
        if (vol_scaling)
        {
            node_radius = node_radius * fabs((double)free_nodes_vol[node_index]) *
                          mass_scale_factor;
        }

        /* Scale node radius by user selected scale factor */
        node_radius *= rfac;

        /* Scale node radius by user selected scale factor */
        node_radius *= rfac;

        /* If a result is available, then color node by the current
         * result.
         */

        if (data_array!=NULL && colorflag)
        {
            if ( p_ml_class )
            {
                int j;
                val = get_ml_result( analy, p_ml_class, elem_id, &result_defined );
                for(i = 0; i < MESH(analy).qty_class_selections; i++)
                {
                    if(!strcmp(p_ml_class->short_name, MESH(analy).by_class_select[i].p_class->short_name))
                    {
                        break;
                    }
                }
                /* if any of the elements are disabled show the material color */
                if(MESH(analy).by_class_select[i].disable_class_elem[elem_id] == TRUE)
                {
                    col[3] = v_win->mesh_materials.diffuse[mat_num][3];

                    VEC_COPY( col, v_win->mesh_materials.diffuse[mat_num] );
                    for (j=0;
                            j<4;
                            j++)
                        col[j] = col[j]- col[j]*.5;
                } else
                {
                    color_lookup( col, val,
                                  rmin, rmax, analy->zero_result, mat_num,
                                  analy->logscale, analy->material_greyscale );
                }
            }
            else
            {
                val = data_array[free_node_data_index];
                color_lookup( col, val,
                              rmin, rmax, analy->zero_result, mat_num,
                              analy->logscale, analy->material_greyscale );
            }
        }
        else
        {
            /* No result, color by material */
            col[3] = v_win->mesh_materials.diffuse[mat_num][3];

            VEC_COPY( col, v_win->mesh_materials.diffuse[mat_num] );
            for (i=0;
                    i<4;
                    i++)
                col[i] = col[i]- col[i]*.5;
        }

        analy->free_nodes_vals[node_index] = val;

        /* Draw spheres at the nodes of the element. If the res for the sphere is
         * <=2, then draw a box instead.
         */
        if ( analy->dimension == 3 )
            if (analy->free_nodes_sphere_res_factor>2)
            {
                glColor3fv( col );
                if ( !fastSpheres )
                    draw_sphere_GL( verts, node_radius, analy->free_nodes_sphere_res_factor);
                else
                {
                    /* gluSphere( sphere, node_base_radius*rfac, analy->free_nodes_sphere_res_factor*2, analy->free_nodes_sphere_res_factor ); */

                    glPushMatrix();
                    glTranslatef( verts[0], verts[1], verts[2] );
                    glCallList( display_list );
                    glPopMatrix();
                }

                /*
                 * If reflection planes are active, then reflect free nodes.
                 */
                if ( analy->reflect )
                    for ( plane = analy->refl_planes;
                            plane != NULL;
                            plane = plane->next )
                    {
                        point_transform( refl_verts, verts, &plane->pt_transf );
                        draw_sphere_GL( refl_verts, node_radius, analy->free_nodes_sphere_res_factor);
                    }
            }

        /* Draw a hex for the free node */

            else
            {
                glBegin( GL_QUADS );
                glColor3fv( col );

                for (j=0; j<6; j++)
                {

                    switch ( j )
                    {
                    case 0:  /* +z */
                        poly_pts[0][0] = verts[0]+node_radius;
                        poly_pts[0][1] = verts[1]+node_radius;
                        poly_pts[0][2] = verts[2]+node_radius;

                        poly_pts[1][0] = verts[0]-node_radius;
                        poly_pts[1][1] = verts[1]+node_radius;
                        poly_pts[1][2] = verts[2]+node_radius;

                        poly_pts[2][0] = verts[0]+node_radius;
                        poly_pts[2][1] = verts[1]-node_radius;
                        poly_pts[2][2] = verts[2]+node_radius;

                        poly_pts[3][0] = verts[0]-node_radius;
                        poly_pts[3][1] = verts[1]-node_radius;
                        poly_pts[3][2] = verts[2]+node_radius;
                        break;

                    case 1:  /* -z */
                        poly_pts[0][0] = verts[0]+node_radius;
                        poly_pts[0][1] = verts[1]+node_radius;
                        poly_pts[0][2] = verts[2]-node_radius;

                        poly_pts[1][0] = verts[0]-node_radius;
                        poly_pts[1][1] = verts[1]+node_radius;
                        poly_pts[1][2] = verts[2]-node_radius;

                        poly_pts[2][0] = verts[0]+node_radius;
                        poly_pts[2][1] = verts[1]-node_radius;
                        poly_pts[2][2] = verts[2]-node_radius;

                        poly_pts[3][0] = verts[0]-node_radius;
                        poly_pts[3][1] = verts[1]-node_radius;
                        poly_pts[3][2] = verts[2]-node_radius;
                        break;

                    case 2:  /* +x */
                        poly_pts[0][0] = verts[0]+node_radius;
                        poly_pts[0][1] = verts[1]+node_radius;
                        poly_pts[0][2] = verts[2]+node_radius;

                        poly_pts[1][0] = verts[0]+node_radius;
                        poly_pts[1][1] = verts[1]+node_radius;
                        poly_pts[1][2] = verts[2]-node_radius;

                        poly_pts[2][0] = verts[0]+node_radius;
                        poly_pts[2][1] = verts[1]-node_radius;
                        poly_pts[2][2] = verts[2]-node_radius;

                        poly_pts[3][0] = verts[0]+node_radius;
                        poly_pts[3][1] = verts[1]-node_radius;
                        poly_pts[3][2] = verts[2]+node_radius;
                        break;

                    case 3:  /* -x */
                        poly_pts[0][0] = verts[0]-node_radius;
                        poly_pts[0][1] = verts[1]+node_radius;
                        poly_pts[0][2] = verts[2]+node_radius;

                        poly_pts[1][0] = verts[0]-node_radius;
                        poly_pts[1][1] = verts[1]+node_radius;
                        poly_pts[1][2] = verts[2]-node_radius;

                        poly_pts[2][0] = verts[0]-node_radius;
                        poly_pts[2][1] = verts[1]-node_radius;
                        poly_pts[2][2] = verts[2]-node_radius;

                        poly_pts[3][0] = verts[0]-node_radius;
                        poly_pts[3][1] = verts[1]-node_radius;
                        poly_pts[3][2] = verts[2]+node_radius;
                        break;

                    case 4:  /* +y */
                        poly_pts[0][0] = verts[0]+node_radius;
                        poly_pts[0][1] = verts[1]+node_radius;
                        poly_pts[0][2] = verts[2]+node_radius;

                        poly_pts[1][0] = verts[0]+node_radius;
                        poly_pts[1][1] = verts[1]+node_radius;
                        poly_pts[1][2] = verts[2]-node_radius;

                        poly_pts[2][0] = verts[0]-node_radius;
                        poly_pts[2][1] = verts[1]+node_radius;
                        poly_pts[2][2] = verts[2]-node_radius;

                        poly_pts[3][0] = verts[0]-node_radius;
                        poly_pts[3][1] = verts[1]+node_radius;
                        poly_pts[3][2] = verts[2]+node_radius;
                        break;

                    case 5:  /* -y */
                        poly_pts[0][0] = verts[0]+node_radius;
                        poly_pts[0][1] = verts[1]-node_radius;
                        poly_pts[0][2] = verts[2]+node_radius;

                        poly_pts[1][0] = verts[0]+node_radius;
                        poly_pts[1][1] = verts[1]-node_radius;
                        poly_pts[1][2] = verts[2]-node_radius;

                        poly_pts[2][0] = verts[0]-node_radius;
                        poly_pts[2][1] = verts[1]-node_radius;
                        poly_pts[2][2] = verts[2]-node_radius;

                        poly_pts[3][0] = verts[0]-node_radius;
                        poly_pts[3][1] = verts[1]-node_radius;
                        poly_pts[3][2] = verts[2]+node_radius;
                        break;
                    }
                    glNormal3fv( poly_norms[j] );
                    for (k=0; k<4; k++)
                        glVertex3fv( poly_pts[k] );
                }
                glEnd();
            }
        free_node_data_index++;
    }


    /*if ( v_win->lighting )
      glDisable( GL_LIGHTING ); */

    glDisable( GL_COLOR_MATERIAL );

    if ( fastSpheres )
    {
        gluDeleteQuadric( sphere );
        glDeleteLists( display_list, 1 );
    }

    free(free_nodes_list);
    free(free_nodes_elem_list);
    free(free_nodes_mass);
    free(free_nodes_vol);
    free(part_nodes_list);
    free(part_nodes_result);
}

#ifdef JPEG_SUPPORT
/*****************************************************************
 * TAG( write_jpeg_file )
 *
 * Write a jpeg file.
 *
 * NOTE:  This procedure references "include files", data structures
 *        and procedure calls obtained from:
 *
 *        The Independent JPEG Group:  JPEG v.6-b
 */
/*
 * IMAGE DATA FORMATS:
 *
 * The standard input image format is a rectangular array of pixels with
 * each pixel having the same number of "component" values (color channels).
 * Each pixel row is an array of:
 *
 *      unsigned char
 *
 *      (or "char" if "unsigned char" is unavailable )
 *
 * When working with color data the color values for each pixel must be
 * adjacent in the row, e.g.:
 *
 *     R, G, B, R, G, B, R, G, B, ...
 *
 * for 24-bit RGB color
 *
 * This procedure assumes such a data structure matches the way GRIZ has
 * stored the image in memory.  As such, a pointer can be passed to the
 * image buffer.  In particular, the image is RGB color and described by:
 * RGB color and is described by:
 *
 *      Pointer to array of RGB-order data:
 *
 *             image_buffer = NEW_N( char, buffer_size, "screen to rgb" );
 *
 *     Number of rows and columns in image:
 *
 *             cinfo.image_width      = v_win->vp_width;
 *             cinfo.image_height     = v_win->vp_height;
 *
 *
 * Scanlines MUST be supplied in top-to-bottom order if the JPEG files are to be
 * compatible with those written by others.
 */


int
write_JPEG_file( char *filename, Bool_type alpha_flag, int quality )
{
    /*
     * This struct contains the JPEG compression parameters and pointers to
     * working space (which is allocated as needed by the JPEG library).
     */

    struct jpeg_compress_struct
            cinfo;


    /*
     * This struct represents a JPEG error handler.  It is declared separately
     * because applications often want to supply a specialized error handler.
     * This procedure will use the standard error handler to print a message
     * on stderr and call exit() if compression fails.
     *
     * NOTE:  This struct must live as long as the main JPEG parameter
     *        struct to avoid dangling-pointer problems.
     */

    struct jpeg_error_mgr
            jerr;


    FILE
    *outfile;

    JSAMPLE
    *image_buffer;     /* pointer to buffer array of R,G,B-order data */

    JSAMPROW
    row_pointer[1];   /* pointer to JSAMPLE row[s] */

    int
    buffer_size            /* size of image_buffer */
    ,index
    ,row_stride;            /* physical row width in image buffer */


    /* */

    /*
     * Step 0:  establishment of error handler
     *
     *          This step must be done first in case the initialization step fails,
     *          This routine fills contents of "struct jerr" and returns address of "jerr"
     *          which is placed into the link field of cinfo.
     */

    cinfo.err = jpeg_std_error( &jerr );


    /*
     * Step 1:  allocate and initialize JPEG compression object
     */

    jpeg_create_compress( &cinfo );


    /*
     * Step 2:  specify data destination file
     *
     *          JPEG library code write compressed data to a stdio stream.
     *
     *          NOTE:  VERY IMPORTANT!
     *
     *                 Use "b" option to fopen() on all platforms requiring such
     *                 a parameter setting in order to write BINARY files.
     */

    if ( (outfile = fopen( filename, "wb" )) == NULL )
{
        popup_dialog( WARNING_POPUP
                      ,"Unable to open JPEG output file:  %s\n"
                      ,filename );

        return( EXIT_FAILURE );
    }

    jpeg_stdio_dest( &cinfo, outfile );


    /*
     * Step 3:  set default parameters for compression
     *
     *          NOTE:  "cinfo.in_color_space" MUST be set PRIOR to calling
     *                 "jpeg_set_defaults" as defaults depend upon the source
     *                 color space.
     *
     *          description of input image:
     *
     *          image width (pixels);
     *          image height (pixels);
     *          number of color components per pixel;
     *          colorspace of imput image
     *
     *          These four fields of the "struct cinfo" MUST be set.
     */

    cinfo.image_width      = v_win->vp_width;
    cinfo.image_height     = v_win->vp_height;
    cinfo.input_components = alpha_flag ? 4 : 3;
    cinfo.in_color_space   = JCS_RGB;

    jpeg_set_defaults( &cinfo );

    /*
     * Step 3a:  incorporate user-specified modifications to default settings
     *
     *           quality:  quantization table scaling
     */

    jpeg_set_quality( &cinfo, quality, TRUE );

    /*
     * Step 4:  start compression
     *
     *          NOTE:  Parameter "TRUE" ensures a complete interchange-JPEG file
     *                 will be written.
     */

    jpeg_start_compress( &cinfo, TRUE );

    /*
     * Step 5:  establish image buffer containing screen pixel data
     */

    buffer_size = cinfo.image_width * cinfo.image_height * cinfo.input_components;

    image_buffer = NEW_N( JSAMPLE, buffer_size, "screen to rgb" );

    screen_to_memory( alpha_flag, cinfo.image_width, cinfo.image_height, image_buffer );

    wrt_text( "Writing data to image file %s\n\n", filename );

    /*
     * Step 6:  convert image scanlines to JPEG format
     *
     *          "component" number of JSAMPLE's per row in image_buffer
     */

    row_stride = cinfo.image_width * cinfo.input_components;

    /*
     * "jpeg_write_scanlines" expects an array of pointers to scanlines
     */

    index = cinfo.image_height - 1;

    while ( index >= 0 )
    {
        row_pointer[0] = &image_buffer[index * row_stride];

        (void) jpeg_write_scanlines( &cinfo, row_pointer, 1 );

        index--;
    }

    /*
     * Step 7:  finish compression
     */

    jpeg_finish_compress( &cinfo );

    fclose( outfile );

    /*
     * Step 8:  release JPEG compression object
     */

    jpeg_destroy_compress( &cinfo );

    return( EXIT_SUCCESS );
}
#else

/* JPEG not enabled */

write_JPEG_file( char *filename, Bool_type alpha_flag, int quality )
{
    printf("\nJPEG is not enabled.");
    return( EXIT_SUCCESS );
}
#endif

#ifdef PNG_SUPPORT
/*****************************************************************
 * TAG( write_png_file )
 *
 * Write a PNG (The Portable Network Graphics) file.
 *
 *
 * NOTE:  This procedure references "include files", data structures
 *        and procedure calls obtained from:
 *
 *        PNG Development Group
 */
/*
 * IMAGE DATA FORMATS:
 *
 *       Conceptually, a PNG image is a rectangular pixel array, with
 *       pixels appearing left-to-right within each scanline, and scanlines
 *       appearing top-to-bottom.  (For progressive display purposes, the
 *       data may actually be transmitted in a different order.  The size
 *       of each pixel is determined by the bit depth, which is the number
 *       of bits per sample in the image data.
 *
 *       Three types of pixel are supported:
 *
 *       * An indexed-color pixel is represented by a single sample
 *         that is an index into a supplied palette.  The image bit
 *         depth determines the maximum number of palette entries, but
 *         not the color precision within the palette.
 *
 *       * A grayscale pixel is represented by a single sample that is
 *         a grayscale level, where zero is black and the largest value
 *         for the bit depth is white.
 *
 *       * A truecolor pixel is represented by three samples: red (zero
 *         = black, max = red) appears first, then green (zero = black,
 *         max = green), then blue (zero = black, max = blue).  The bit
 *         depth specifies the size of each sample, not the total pixel
 *         size.
 *
 *       Optionally, grayscale and truecolor pixels can also include an
 *       alpha sample.
 *
 *       Pixels are always packed into scanlines with no wasted bits
 *       between pixels.  Pixels smaller than a byte never cross byte
 *       boundaries; they are packed into bytes with the leftmost pixel in
 *       the high-order bits of a byte, the rightmost in the low-order
 *       bits.  Permitted bit depths and pixel types are restricted so that
 *       in all cases the packing is simple and efficient.
 *
 *       Scanlines always begin on byte boundaries.  When pixels have fewer
 *       than 8 bits and the scanline width is not evenly divisible by the
 *       number of pixels per byte, the low-order bits in the last byte of
 *       each scanline are wasted.  The contents of these wasted bits are
 *       unspecified.
 *
 *       A PNG image can be stored in interlaced order to allow progressive
 *       display.  The purpose of this feature is to allow images to "fade
 *       in" when they are being displayed on-the-fly.  Interlacing
 *       slightly expands the file size on average, but it gives the user a
 *       meaningful display much more rapidly.  Note that decoders are
 *       required to be able to read interlaced images, whether or not they
 *       actually perform progressive display.
 *
 *       With interlace method 0, pixels are stored sequentially from left
 *       to right, and scanlines sequentially from top to bottom (no
 *       interlacing).
 *
 *       Interlace method 1, known as Adam7 after its author, Adam M.
 *       Costello, consists of seven distinct passes over the image.  Each
 *       pass transmits a subset of the pixels in the image.  The pass in
 *       which each pixel is transmitted is defined by replicating the
 *       following 8-by-8 pattern over the entire image, starting at the
 *       upper left corner:
 *
 *       1 6 4 6 2 6 4 6
 *       7 7 7 7 7 7 7 7
 *       5 6 5 6 5 6 5 6
 *       7 7 7 7 7 7 7 7
 *       3 6 4 6 3 6 4 6
 *       7 7 7 7 7 7 7 7
 *       5 6 5 6 5 6 5 6
 *       7 7 7 7 7 7 7 7
 *
 *       A PNG file consists of a PNG signature followed by a series of
 *       chunks.  This chapter defines the signature and the basic properties
 *       of chunks.  Individual chunk types are discussed in the next chapter.
 *
 *       The first eight bytes of a PNG file always contain the following
 *       (decimal) values:
 *
 *       137 80 78 71 13 10 26 10
 *
 *       This signature indicates that the remainder of the file contains a
 *       single PNG image, consisting of a series of chunks beginning with
 *       an IHDR chunk and ending with an IEND chunk.
 *
 *       Each chunk consists of four parts:
 *
 *       Length
 *          A 4-byte unsigned integer giving the number of bytes in the
 *          chunk's data field.  The length counts only the data field, not
 *          itself, the chunk type code, or the CRC.  Zero is a valid
 *          length.  Although encoders and decoders should treat the length
 *          as unsigned, its value must not exceed (2^31)-1 bytes.
 *
 *       Chunk Type
 *          A 4-byte chunk type code.  For convenience in description and
 *          in examining PNG files, type codes are restricted to consist of
 *          uppercase and lowercase ASCII letters (A-Z and a-z, or 65-90
 *          and 97-122 decimal).  However, encoders and decoders must treat
 *          the codes as fixed binary values, not character strings.  For
 *          example, it would not be correct to represent the type code
 *          IDAT by the EBCDIC equivalents of those letters.
 *
 *       Chunk Data
 *          The data bytes appropriate to the chunk type, if any.  This
 *          field can be of zero length.
 *
 *       CRC
 *          A 4-byte CRC (Cyclic Redundancy Check) calculated on the
 *          preceding bytes in the chunk, including the chunk type code and
 *          chunk data fields, but not including the length field.  The CRC
 *          is always present, even for chunks containing no data.
 *
 *
 *   All implementations must understand and successfully render the
 *   standard critical chunks.  A valid PNG image must contain an IHDR
 *   chunk, one or more IDAT chunks, and an IEND chunk.
 *
 *      The IHDR chunk must appear FIRST.  It contains:

 *         Width:              4 bytes
 *         Height:             4 bytes
 *         Bit depth:          1 byte
 *         Color type:         1 byte
 *         Compression method: 1 byte
 *         Filter method:      1 byte
 *         Interlace method:   1 byte
 *
 *      Width and height give the image dimensions in pixels.  They are
 *      4-byte integers.  Zero is an invalid value.  The maximum for
 *      each is (2^31)-1 in order to accommodate languages that have
 *      difficulty with unsigned 4-byte values.
 *
 *      Bit depth is a single-byte integer giving the number of bits
 *      per sample or per palette index (not per pixel).  Valid values
 *      are 1, 2, 4, 8, and 16, although not all values are allowed for
 *      all color types.
 *
 *      Color type is a single-byte integer that describes the
 *      interpretation of the image data.  Color type codes represent
 *      sums of the following values: 1 (palette used), 2 (color used),
 *      and 4 (alpha channel used).  Valid values are 0, 2, 3, 4, and
 *      6.
 *
 *      Bit depth restrictions for each color type are imposed to
 *      simplify implementations and to prohibit combinations that do
 *      not compress well.  Decoders must support all valid
 *      combinations of bit depth and color type.  The allowed
 *      combinations are:
 *
 *         Color    Allowed    Interpretation
 *         Type    Bit Depths
 *
 *         0       1,2,4,8,16  Each pixel is a grayscale sample.
 *
 *         2       8,16        Each pixel is an R,G,B triple.
 *
 *         3       1,2,4,8     Each pixel is a palette index;
 *                             a PLTE chunk must appear.
 *
 *         4       8,16        Each pixel is a grayscale sample,
 *                             followed by an alpha sample.
 *
 *         6       8,16        Each pixel is an R,G,B triple,
 *                             followed by an alpha sample.
 *
 *      The sample depth is the same as the bit depth except in the
 *      case of color type 3, in which the sample depth is always 8
 *      bits.
 *
 *      Compression method is a single-byte integer that indicates the
 *      method used to compress the image data.  At present, only
 *      compression method 0 (deflate/inflate compression with a
 *      sliding window of at most 32768 bytes) is defined.  All
 *      standard PNG images must be compressed with this scheme.  The
 *      compression method field is provided for possible future
 *      expansion or proprietary variants.  Decoders must check this
 *      byte and report an error if it holds an unrecognized code.
 *
 *      Filter method is a single-byte integer that indicates the
 *      preprocessing method applied to the image data before
 *      compression.  At present, only filter method 0 (adaptive
 *      filtering with five basic filter types) is defined.  As with
 *      the compression method field, decoders must check this byte and
 *      report an error if it holds an unrecognized code.
 *
 *      Interlace method is a single-byte integer that indicates the
 *      transmission order of the image data.  Two values are currently
 *      defined: 0 (no interlace) or 1 (Adam7 interlace).
 *
 * When working with color data the color values for each pixel must be
 * adjacent in the row, e.g.:
 *
 *     R, G, B, R, G, B, R, G, B, ...
 *
 * for 24-bit RGB color
 *
 * This procedure assumes such a data structure matches the way GRIZ has
 * stored the image in memory.  As such, a pointer can be passed to the
 * image buffer.  In particular, the image is RGB color and described by:
 *
 *      Pointer to array of RGB-order data:
 *
 *             pix_ptr = NEW_N( char, buffer_size, "screen to rgb" );
 *
 *     Number of rows and columns in image:
 *
 *             image_width      = v_win->vp_width;
 *             image_height     = v_win->vp_height;
 *
 *
 * Scanlines MUST be supplied in top-to-bottom order
 */

int
write_PNG_file( char *filename, Bool_type alpha )
{
    FILE          *png_file;
    png_structp   png_ptr = NULL;
    png_infop     info_ptr = NULL;
    png_byte      *png_pixels = NULL;
    png_byte      **row_pointers = NULL;
    png_uint_32   row_bytes;
    png_uint_32   width;
    png_uint_32   height;
    int           buffer_size;
    int           color_type;
    int           interlace_type;
    int           bit_depth = 0;
    int           channels;
    int           i;


    width    = v_win->vp_width;
    height   = v_win->vp_height;
    channels = alpha ? 4 : 3;

    if ( (png_file = fopen( filename, "wb" )) == NULL )
    {
        popup_dialog( WARNING_POPUP
                      ,"Unable to open PNG output file:  %s\n"
                      ,filename );

        return( EXIT_FAILURE );
    }

    /* prepare the standard PNG structures */
    png_ptr = png_create_write_struct (PNG_LIBPNG_VER_STRING, NULL, NULL, NULL );
    if (!png_ptr)
    {
        fclose( png_file );
        return( EXIT_FAILURE );
    }
    info_ptr = png_create_info_struct (png_ptr);
    if (!info_ptr)
    {
        png_destroy_write_struct (&png_ptr, (png_infopp) NULL);
        fclose( png_file );
        return( EXIT_FAILURE );
    }

    /* setjmp() must be called in every function that calls a PNG-reading libpng function */
    if (setjmp (png_jmpbuf(png_ptr)))
    {
        png_destroy_write_struct (&png_ptr, (png_infopp) NULL);
        fclose( png_file );
        return( EXIT_FAILURE );
    }

    /* initialize the png structure */
    png_init_io (png_ptr, png_file);


    /* set the zlib compression level */
    png_set_compression_level(png_ptr, Z_DEFAULT_COMPRESSION);
    png_set_compression_strategy(png_ptr, Z_DEFAULT_STRATEGY);

    bit_depth = 8;
    color_type = alpha ? PNG_COLOR_TYPE_RGB_ALPHA : PNG_COLOR_TYPE_RGB;
    interlace_type = PNG_INTERLACE_NONE;

    png_set_IHDR (png_ptr, info_ptr, width, height, bit_depth, color_type,
                  interlace_type, PNG_COMPRESSION_TYPE_DEFAULT,
                  PNG_FILTER_TYPE_DEFAULT);

    /* write the file header information */
    png_write_info (png_ptr, info_ptr);

    png_set_packing( png_ptr );

    /*
     * Image buffer containing screen pixel data
     */

    buffer_size = width * height * channels;
    png_pixels = NEW_N( png_byte, buffer_size, "screen to rgb" );

    screen_to_memory( alpha, width, height, png_pixels );

    wrt_text( "Writing data to image file %s\n\n", filename );

    /* if needed we will allocate memory for an new array of row-pointers */
    if (row_pointers == (unsigned char**) NULL)
    {
        if ((row_pointers = (png_byte **) malloc (height * sizeof (png_bytep))) == NULL)
        {
            png_destroy_write_struct (&png_ptr, (png_infopp) NULL);
            fclose( png_file );
            return( EXIT_FAILURE );
        }
    }

    /* set the individual row_pointers to point at the correct offsets */
    /* row_bytes is the width x number of channels x (bit-depth / 8) */
    row_bytes = width * channels;
    for (i = 0; i < height; i++)
        row_pointers[i] = png_pixels + ((height-1) - i) * row_bytes;

    /* write out the entire image data in one call */
    png_write_image (png_ptr, row_pointers);

    /* write the additional chuncks to the PNG file (not really needed) */
    png_write_end (png_ptr, info_ptr);

    /* clean up after the write, and free any memory allocated */
    png_destroy_write_struct (&png_ptr, (png_infopp) NULL);

    if (row_pointers != (unsigned char**) NULL)
        free (row_pointers);
    if (png_pixels != (unsigned char*) NULL)
        free (png_pixels);

    /* close output file */
    fclose (png_file);

    return( EXIT_SUCCESS );
}
#endif


/************************************************************
 * TAG( hide_by_object_type )
 *
 * Added October 26, 2006: IRC
 *
 * Determines if an element is to be hidden based on the elements
 * object type.
 *
 */
Bool_type
hide_by_object_type( MO_class_data *p_class, int mat_num, int elm, Analysis *analy, float *data_array )
{
    Bool_type hide_elem   = FALSE;
    Bool_type elem_result = TRUE;
    Mesh_data *p_mesh;
    int class_select_index=0;
    int obj_type;
    int i;

    obj_type = p_class->superclass;

    if (is_particle_class( analy, p_class->superclass, p_class->short_name ) )
        obj_type = M_PARTICLE;

    p_mesh = MESH_P( analy );

    if ( !analy->interp_mode == NO_INTERP )
        elem_result = FALSE;

    switch ( obj_type )
    {
    case M_HEX:

        /* Hide by result range */
        if ( data_array != NULL )
            if ( p_mesh->hide_brick_by_result && elem_result )
                if ( data_array[elm] >= p_mesh->brick_result_min &&
                        data_array[elm] <= p_mesh->brick_result_max )
                    hide_elem = TRUE;

        /* Hide by material number */
        if ( p_mesh->hide_brick != NULL && !hide_elem )
            if ( p_mesh->hide_brick[mat_num] )
                hide_elem = TRUE;

        /* Hide by element number */
        if ( !p_mesh->hide_brick_by_mat && !hide_elem )
        {
            if ( p_mesh->hide_brick_elem[elm] )
                hide_elem = TRUE;
            else
                hide_elem = FALSE;
        }
        break;

    case M_QUAD:

        /* Hide by result range */
        if ( data_array != NULL )
            if ( p_mesh->hide_brick_by_result && elem_result )
                if ( data_array[elm] >= p_mesh->shell_result_min &&
                        data_array[elm] <= p_mesh->shell_result_max )
                    hide_elem = TRUE;

        /* Hide by material number */
        if ( p_mesh->hide_shell != NULL && !hide_elem )
            if ( p_mesh->hide_shell[mat_num] )
                hide_elem = TRUE;

        /* Hide by element number */
        if ( !p_mesh->hide_shell_by_mat && !hide_elem )
        {
            if ( p_mesh->hide_shell_elem[elm] )
                hide_elem = TRUE;
            else
                hide_elem = FALSE;
        }
        break;

    case M_TRUSS:

        /* Hide by result range */
        if ( data_array != NULL )
            if ( p_mesh->hide_brick_by_result && elem_result )
                if ( data_array[elm] >= p_mesh->truss_result_min &&
                        data_array[elm] <= p_mesh->truss_result_max )
                    hide_elem = TRUE;

        /* Hide by material number */
        if ( p_mesh->hide_truss != NULL)
            if ( p_mesh->hide_truss[mat_num] )
                hide_elem = TRUE;

        /* Hide by element number */
        if ( !p_mesh->hide_truss_by_mat && !hide_elem )
        {
            if ( p_mesh->hide_truss_elem[elm] )
                hide_elem = TRUE;
            else
                hide_elem = FALSE;
        }
        break;

    case M_BEAM:

        /* Hide by result range */
        if ( data_array != NULL )
            if ( p_mesh->hide_beam_by_result && elem_result )
                if ( data_array[elm] >= p_mesh->beam_result_min &&
                        data_array[elm] <= p_mesh->beam_result_max )
                    hide_elem = TRUE;

        /* Hide by material number */
        if ( p_mesh->hide_beam != NULL)
            if ( p_mesh->hide_beam[mat_num] )
                hide_elem = TRUE;

        /* Hide by element number */
        if ( !p_mesh->hide_beam_by_mat && !hide_elem )
        {
            if ( p_mesh->hide_beam_elem[elm] )
                hide_elem = TRUE;
            else
                hide_elem = FALSE;
        }
        break;

    case M_PARTICLE:
        /* Hide by result range */
        if ( data_array != NULL )
            if ( p_mesh->hide_particle_by_result && elem_result )
                if ( data_array[elm] >= p_mesh->particle_result_min &&
                        data_array[elm] <= p_mesh->particle_result_max )
                    hide_elem = TRUE;

        if ( p_mesh->hide_particle[mat_num] )
            hide_elem = TRUE;
        else
            hide_elem = FALSE;
        break;
    }


    /* Check for class selections */
    class_select_index = get_class_select_index( analy, p_class->short_name );

    if ( !hide_elem && class_select_index>=0 && elm<p_class->qty)
    {
        if ( p_mesh->by_class_select[class_select_index].hide_class_elem[elm] || 
             p_mesh->by_class_select[class_select_index].exclude_class_elem[elm])
            hide_elem = TRUE;
    }
    return( hide_elem );
}


/************************************************************
 * TAG( disable_by_object_type )
 *
 * Added October 26, 2006: IRC
 *
 * Determines if an element is to be disabled based on the elements
 * object type.
 *
 */

Bool_type
disable_by_object_type( MO_class_data *p_class, int mat_num, int elm, Analysis *analy, float *data_array )
{
    Bool_type elem_result  = TRUE;
    Bool_type disable_elem = FALSE;
    Mesh_data *p_mesh;
    int obj_type;
    int i;
    int class_select_index = 0;

    obj_type = p_class->superclass;
    if (is_particle_class( analy, p_class->superclass, p_class->short_name ) )
        obj_type = M_PARTICLE;

    p_mesh = MESH_P( analy );
    if ( !analy->interp_mode == NO_INTERP )
        elem_result = FALSE;

    switch ( obj_type )
    {
    case M_HEX:
        /* Hide by material number */
        if ( p_mesh->disable_brick != NULL && p_mesh->disable_brick )
            if ( p_mesh->disable_brick[mat_num] )
                disable_elem = TRUE;

        /* Disable by element number */
        if ( !p_mesh->disable_brick_by_mat && !disable_elem && p_mesh->disable_brick_elem )
        {
            if ( p_mesh->disable_brick_elem[elm] )
                disable_elem = TRUE;
        }
        break;

    case M_QUAD:
        /* Hide by material number */
        if ( p_mesh->disable_shell != NULL)
            if ( p_mesh->disable_shell[mat_num] )
                disable_elem = TRUE;

        /* Disable by element number */
        if ( !p_mesh->disable_shell_by_mat && p_mesh->disable_shell_elem )
        {
            if ( p_mesh->disable_shell_elem[elm] && !disable_elem )
                disable_elem = TRUE;
        }
        break;

    case M_TRUSS:
        /* Hide by material number */
        if ( p_mesh->disable_truss != NULL)
            if ( p_mesh->disable_truss[mat_num] )
                disable_elem = TRUE;

        /* Disable by element number */
        if ( !p_mesh->disable_truss_by_mat && !disable_elem && p_mesh->disable_truss_elem )
        {
            if ( p_mesh->disable_truss_elem[elm] )
                disable_elem = TRUE;
        }
        break;

    case M_BEAM:
        /* Hide by material number */
        if ( p_mesh->disable_beam != NULL)
            if ( p_mesh->disable_beam[mat_num] )
                disable_elem = TRUE;

        /* Disable by element number */
        if ( !p_mesh->disable_beam_by_mat && !disable_elem && p_mesh->disable_beam_elem )
        {
            if ( p_mesh->disable_beam_elem[elm] )
                disable_elem = TRUE;
        }
        break;

    case M_PARTICLE:
        if ( p_mesh->disable_particle )
            if ( p_mesh->disable_particle[mat_num] && !disable_elem )
                disable_elem = TRUE;

        break;
    }


    class_select_index = get_class_select_index( analy, p_class->short_name );

    if ( !disable_elem && class_select_index>=0 && elm<p_class->qty)
    {
        if ( p_mesh->by_class_select[class_select_index].disable_class_elem[elm] ||
             p_mesh->by_class_select[class_select_index].exclude_class_elem[elm])
            disable_elem = TRUE;
    }
    return disable_elem;
}


/************************************************************
 * TAG( get_class_label )
 *
 * Added November 05, 2007: IRC
 *
 * Returns the label for the specified object index.
 *
 */

int
get_class_label( MO_class_data *class, int object_index )
{
    int label_index, i;

    /* If we have no labels then the index is still the label, so
     * just return!
    */

    if ( !class )
        return( object_index );

    if ( !class->labels_found )
        return( object_index );

    if ( object_index>class->qty )
        return (M_INVALID_LABEL);

    label_index = class->labels_index[object_index];

    /* We have a valid label - return it */
    return ( class->labels[label_index].label_num );
}

/************************************************************
 * TAG( dump_class_labels )
 *
 * Added August 19, 2011: IRC
 *
 * Returns the label for the specified object index.
 *
 */

void
dump_class_labels( MO_class_data *class )
{
    FILE *fp;
    int label_index, i;
    char labelsFilename[64]="Labels-";
    strcat(labelsFilename, class->short_name);

    fp = fopen( labelsFilename, "w" );
    for ( i=0;
            i<class->qty;
            i++ )
    {
        label_index = class->labels_index[i];
        fprintf( fp, "\nId=%d \tLabel=%8d \tIndex=%8d",
                 class->labels[i].local_id, class->labels[i].label_num-1, label_index );
    }
    fclose(fp);
}


/************************************************************
 * TAG( get_class_label_index )
 *
 * Added November 05, 2007: IRC
 *
 * Returns the index for for the specified label.
 *
 */

int
label_compare( const int *key, const MO_class_labels *label )
{
    int test;
    test = label->label_num;
    return( *key - label->label_num );
}



int
get_class_label_index( MO_class_data *class, int label_num )
{
    int label_index;
    MO_class_labels *result_ptr;

    if ( !class->labels_found )
        return( label_num -1);

    /* Perform a binary search on the label array to detrermine
     * the index number for this label.
     */
    result_ptr = bsearch( &label_num, &class->labels[0], class->qty,
                          sizeof(class->labels[0]), label_compare );
    if ( result_ptr )
        return( result_ptr->local_id );
    else
        return (M_INVALID_LABEL);

}


/************************************************************
 * TAG( populate_result )
 *
 * Added Jan 17 2014: William Oliver 
 * Takes a result array and loops thru the subrecords and the map array and puts a 1
 * in the indexes representing the element number if that element has that subrecord variable.
 *
 * Assumptions:
 * By the time this function is called the current result has already been processed
 *
 */

void populate_result(int superclass, char map[], int size, MO_class_data * p_class, Analysis * analy)
{
    int qty_blocks = 0;
    int total_blocks = 0;
    int i, j, k, mtl;
    int start, end;
    int qty_classes;
    int *range = NULL;
    Result *p_result;
    Subrec_obj *p_subrec;

    /* Add code to return FALSE if the element is associated with a material that has been disabled */

    MO_class_data ** mo_classes;
    Mesh_data *p_mesh;

    p_result = analy->cur_result;

    /* Just a precautionary step */
    if(p_result == NULL)
    {
        return;
    }
    /*if(p_result->origin.is_derived)
    {j
       return TRUE;
    } */
    if(map == NULL)
    {
        return;
    }
    p_mesh = MESH_P(analy);
    qty_classes = p_mesh->classes_by_sclass[superclass].qty;
    mo_classes = (MO_class_data **) p_mesh->classes_by_sclass[superclass].list;

    /* initialize the map array */
    for(i = 0; i < size; i++)
    {
        map[i] = '0';
    }
    
    /* loop thru the subrecs associated with this result to get the total
 *     number of subrec blocks*/ 
    for(i = 0; i < p_result->qty; i++)
    {
        /*p_subrec = analy->srec_tree[0].subrecs + p_result->subrecs[p_result->subrecs[i]]; */
        p_subrec = analy->srec_tree[0].subrecs + p_result->subrecs[i];
        /* make sure we have the correct superclass */
        qty_blocks += p_subrec->subrec.qty_blocks;
    }


    total_blocks = qty_blocks*2;

    /*now that we know the total quantity of subrecord blocks allocate the range array*/ 
    range = (int *) malloc(2*qty_blocks*sizeof(int));
    /* check to make sure we are not out of memory */ 
    if(range == NULL)
    {
        abort();
    }

    j = 0;
    qty_blocks = 0;
    /* now loop thru the subrecs again and populate the range array */
    for(i = 0; i < p_result->qty; i++)
    {
        /*p_subrec = analy->srec_tree[0].subrecs + p_result->subrecs[p_result->subrecs[i]]; */
        p_subrec = analy->srec_tree[0].subrecs + p_result->subrecs[i];

        qty_blocks = p_subrec->subrec.qty_blocks;
         /* make sure we have the correct superclass*/ 
        for(k = 0; k < qty_blocks*2; k+=2)
        {
            range[j] = p_subrec->subrec.mo_blocks[k];
            j++;
            range[j] = p_subrec->subrec.mo_blocks[k+1];
            j++;
        }
    }

    for(i = 0; i < total_blocks; i+= 2)
    {
        start = range[i];
        end = range[i+1];
        if(end >= size)
        {
            end = size - 1;
        }
        for(j = start; j <= end; j++)
        {
            map[j] = '1';
        }
    }

    free(range);
   
}


/************************************************************
 * TAG( is_particle_class )
 *
 * Added November 05, 2009: IRC
 *
 * Returns TRUE if this class is a particle class.
 *
 */

Bool_type
is_particle_class( Analysis *analy, int superclass, char *class_name )
{
    char short_name[M_MAX_NAME_LEN],
         short_name_upper[M_MAX_NAME_LEN],
         class_name_upper[M_MAX_NAME_LEN];
    int i;

    if ( superclass==G_PARTICLE )
        return( TRUE );

    strcpy( short_name, class_name );

    /* Convert to uppercase */

    string_to_upper( short_name, short_name_upper );

    if ( analy->mesh_table->num_particle_classes==0 )
    {
        if ( !strcmp( short_name_upper, "PARTICLE" )       ||
                !strcmp( short_name_upper, "PARTICLE_ELEM" )  ||
                !strcmp( short_name_upper, "ML" )             ||
                !strncmp( short_name_upper, "SPH", 3 )        ||
                !strncmp( short_name_upper, "DBC", 3 )  )
            return ( TRUE );
        return ( FALSE );
    }

    for ( i=0;
            i<analy->mesh_table->num_particle_classes;
            i++ )
    {
        string_to_upper( analy->mesh_table->particle_class_names[i], class_name_upper );
        if ( !strcmp( short_name_upper, class_name_upper ) )
            return( TRUE );
    }

    if ( !strcmp( short_name_upper, "PARTICLE" )       ||
            !strcmp( short_name_upper, "PARTICLE_ELEM" )  ||
            !strcmp( short_name_upper, "NTET" ))
        return ( TRUE );

    return ( FALSE );
}

Bool_type
is_dbc_class( Analysis *analy, int superclass, char *class_name )
{
    char short_name[M_MAX_NAME_LEN],
         short_name_upper[M_MAX_NAME_LEN];

    strcpy( short_name, class_name );

    /* Convert to uppercase */

    string_to_upper( short_name, short_name_upper );

    if ( !strncmp( short_name_upper, "DBC", 3 ) )
        return ( TRUE );

    return ( FALSE );
}

/************************************************************
 * TAG( is_brick_class )
 *
 * Added November 05, 2009: IRC
 *
 * Returns TRUE if this class is a brick class.
 *
 */

Bool_type
is_brick_class( Analysis *analy, char *short_name )
{
    char short_name_upper[M_MAX_NAME_LEN];
    int i;

    /* Convert to uppercase */
    string_to_upper( short_name, short_name_upper );

    if ( strstr(short_name_upper, "BRICK_") )
        return( TRUE );

    return( FALSE );
}

/************************************************************
 * TAG( is_elem_class )
 *
 * Added November 05, 2012: IRC
 *
 * Returns TRUE if this class is an element type class.
 *


Bool_type
is_elem_class( Analysis *analy, char *short_name )
{
  char short_name_upper[M_MAX_NAME_LEN], class_name_upper[M_MAX_NAME_LEN];
  MO_class_data **class_array, *p_mocd;
  Mesh_data *p_md;
  int class_qty=0;
  int i;
  int status=OK;
  string_to_upper( short_name, short_name_upper ); /* Make case insensitive

  p_md = MESH_P( analy );
  status = htable_get_data( p_md->class_table,
			    (void ***) &class_array,
			    &class_qty);
  for ( i = 0;
	i < class_qty;
	i++ )
  {
        p_mocd = class_array[i];
        if ( !strcmp( p_mocd->short_name, short_name ) ||
	     !strcmp( p_mocd->short_name, short_name_upper ) ) {
	     return( TRUE );
	}
  }

  return( FALSE );
} */


Bool_type
is_elem_class(int superclass)
{
    if((superclass > M_NODE) && (superclass <= M_HEX))
    {
        return TRUE;
    }
    else if(superclass == M_PARTICLE)
    {
        return TRUE;
    }

    return FALSE;

}

/************************************************************
 * TAG( calc_particle_radius )
 *
 * Returns the rendered radius for a particle or free-node.
 *
 */
float calc_particle_radius( Analysis *analy, float scale_factor_input )
{
    int i;
    float leng[3], radius, node_base_radius, scale_factor ;

    if ( scale_factor_input<=0 )
        scale_factor = analy->free_nodes_scale_factor;
    else
        scale_factor = scale_factor_input;

    for ( i = 0;
            i < analy->dimension;
            i++ )
        leng[i] = analy->bbox[1][i] - analy->bbox[0][i];

    node_base_radius = 0.01 * (leng[0] + leng[1] + leng[2]) / 3.0;

    radius = node_base_radius*scale_factor*1.2;
    return( radius );
}


/************************************************************
 * TAG( hide_particles )
 *
 * Hides particles my material or for all materials.
 *
 */
void hide_particles( Analysis *analy,
                     Bool_type all_flag, Bool_type invis,
                     int mat )
{
    Mesh_data     *p_mesh;
    MO_class_data *p_class;
    MO_class_data *p_node_geom;
    MO_class_data **mo_classes;

    int i, num_mats,
        mat_index;
    int qty_classes=0;

    p_mesh      = MESH_P( analy );
    num_mats    = p_mesh->material_qty;
    if ( mat<0 || mat>num_mats )
        return;
    p_node_geom = p_mesh->node_geom;
    qty_classes = p_mesh->classes_by_sclass[G_HEX].qty;
    mo_classes  = (MO_class_data **) p_mesh->classes_by_sclass[G_HEX].list;

    for ( i = 0;
            i < qty_classes;
            i++ )
    {
        p_class = mo_classes[i];
        if ( is_particle_class( analy, p_class->superclass, p_class->short_name ) )
        {
            if ( all_flag )
            {
                for ( mat_index=0;
                        mat_index<num_mats;
                        mat_index++ )
                    p_mesh->hide_particle[mat_index] = invis;
            }
            else
            {
                p_mesh->hide_particle[mat] = invis;
            }
        }
    }

    qty_classes = p_mesh->classes_by_sclass[G_PARTICLE].qty;
    mo_classes  = (MO_class_data **) p_mesh->classes_by_sclass[G_PARTICLE].list;

    for ( i = 0;
            i < qty_classes;
            i++ )
    {
        p_class = mo_classes[i];
        if ( is_particle_class( analy, p_class->superclass, p_class->short_name ) )
        {
            if ( all_flag )
            {
                for ( mat_index=0;
                        mat_index<num_mats;
                        mat_index++ )
                    p_mesh->hide_particle[mat_index] = invis;
            }
            else
            {
                p_mesh->hide_particle[mat] = invis;
            }
        }
    }
}

/************************************************************
 * TAG( disable_particles )
 *
 * Disables particles my material or for all materials.
 *
 */
void disable_particles( Analysis *analy,
                        Bool_type all_flag, Bool_type invis,
                        int mat )
{
    Mesh_data     *p_mesh;
    MO_class_data *p_class;
    MO_class_data *p_node_geom;
    MO_class_data **mo_classes;

    int i, num_mats,
        mat_index;
    int qty_classes=0;

    p_mesh      = MESH_P( analy );
    num_mats    = p_mesh->material_qty;
    if ( mat<0 || mat>num_mats )
        return;
    p_node_geom = p_mesh->node_geom;
    qty_classes = p_mesh->classes_by_sclass[G_HEX].qty;
    mo_classes = (MO_class_data **) p_mesh->classes_by_sclass[G_HEX].list;

    for ( i = 0;
            i < qty_classes;
            i++ )
    {
        p_class = mo_classes[i];
        if ( is_particle_class( analy, p_class->superclass, p_class->short_name ) )
        {
            if ( all_flag )
            {
                for ( mat_index=0;
                        mat_index<num_mats;
                        mat_index++ )
                    p_mesh->disable_particle[mat_index] = invis;
            }
            else
            {
                p_mesh->disable_particle[mat] = invis;
            }
        }
    }

    qty_classes = p_mesh->classes_by_sclass[G_PARTICLE].qty;
    mo_classes = (MO_class_data **) p_mesh->classes_by_sclass[G_PARTICLE].list;

    for ( i = 0;
            i < qty_classes;
            i++ )
    {
        p_class = mo_classes[i];
        if ( is_particle_class( analy, p_class->superclass, p_class->short_name ) )
        {
            if ( all_flag )
            {
                for ( mat_index=0;
                        mat_index<num_mats;
                        mat_index++ )
                    p_mesh->disable_particle[mat_index] = invis;
            }
            else
            {
                p_mesh->disable_particle[mat] = invis;
            }
        }
    }
}

/************************************************************
 * TAG( DrawCone )
 *
 * Draws a cone at specified location.
 *
 */
void DrawCone(float len, float base_diam, float cols[4][4], int res )
{
    int i;
    GLdouble base = .15, height=10., top_diam=0.;
    GLint    slices=2, stacks=2;

    float temp1=0.0;

    temp1=slices*res;
    slices = (int) temp1;
    temp1=stacks*res;
    stacks = (int) temp1;

    top_diam = base_diam*.05;

    glBegin(GL_LINE_LOOP);
    GLUquadricObj* quadric = gluNewQuadric();
    glColor3fv( cols[0] );
    gluQuadricDrawStyle(quadric, GLU_FILL);
    gluQuadricOrientation(quadric, GLU_OUTSIDE);
    gluCylinder(quadric, base_diam, top_diam, len, slices, stacks);
    gluDisk (quadric, 0., base_diam, slices, 10);
    gluDeleteQuadric(quadric);
    glEnd();
}
