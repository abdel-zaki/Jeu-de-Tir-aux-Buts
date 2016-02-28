/*****************************************************************************
File: KP-anim-modele.cpp

Virtual Humans
Master in Computer Science
Christian Jacquemin, University Paris 11

Copyright (C) 2008 University Paris 11 
This file is provided without support, instruction, or implied
warranty of any kind.  University Paris 11 makes no guarantee of its
fitness for a particular purpose and is not liable under any
circumstances for any damages or loss whatsoever arising from the use
or inability to use this file or items derived from it.
******************************************************************************/
//////////////////////////////////////////////////////////////////
// EXTERNAL LIBRARY urotate.cpp
//////////////////////////////////////////////////////////////////

#include "glutil.h"

//////////////////////////////////////////////////////////////////
// CONSTANTS FOR BONES
//////////////////////////////////////////////////////////////////
#define    NBMAXBONES     200
#define    MAXBONEWEIGHTS   6

//////////////////////////////////////////////////////////////////
// MATHEMATICAL LIBRARY
//////////////////////////////////////////////////////////////////

#include <GL/gl.h>           
#include <GL/glu.h>         
#include <GL/glut.h>    

#include <stdio.h>      
#include <stdlib.h>     
#include <string.h>     
#include <math.h>
#include <values.h>
#include <ctype.h>
#include <limits>
#include <iostream>
#include <sstream>

#include "texpng.h"
#include <png.h>
#define    MAX_TEXTURES 10
#include "pos_vert.cpp"
#define    KEY_ESC  27
#define    PI       3.1415926535898

#define    WINDOW_WIDTH   700
#define    WINDOW_HEIGHT  700

#define    NBMAXVERTICES  300000
#define    NBMAXFACES     100000
#define    NBMAXMESH      30
#define    NBMAXKP        100

#define    MAXKPWEIGHTS   8

#define    STRINGSIZE 80

#include <tiffio.h>     /* Sam Leffler's libtiff library. */

#include <fcntl.h>
#include <errno.h>
#include <map>
#ifndef _WIN32
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <netdb.h>
#include <unistd.h>
#else
//#define socklen_t int
#include <winsock2.h>
#include <Ws2tcpip.h>
//#include <wspiapi.h>
#endif

using namespace std;
enum DisplayType{ SURFACE = 0 , MESH , FULL , EmptyDisplayType };
enum WeightingType{ NoWeighting = 0 , Weighting , WeightSubstitution };

DisplayType TypeOfSurfaceDisplay = MESH;

// light properties of infinite light sources
static GLfloat ambient_light0[4] = { 0.35 , 0.35 , 0.35 , 1.0 };
static GLfloat diffuse_light0[4] = { 0.7 , 0.7 , 0.7 , 1.0 };
static GLfloat specular_light0[4] = { 0.1 , 0.1 , 0.1 , 1.0 };
static GLfloat lightDir[] = { 3.0f, 6.0f, 5.0f };
static GLfloat ambient_light1[4] = { 0.35 , 0.35 , 0.35 , 1.0 };
static GLfloat diffuse_light1[4] = { 1.0 , 1.0 , 1.0 , 1.0 };
static GLfloat specular_light1[4] = { 0.1 , 0.1 , 0.1 , 1.0 };

static GLfloat ambient_light2[4] = { 0.35 , 0.35 , 0.35 , 1.0 };
static GLfloat diffuse_light2[4] = { 0.50 , 0.50 , 0.50 , 1.0 };
static GLfloat specular_light2[4] = { 0.1 , 0.1 , 0.1 , 1.0 };

// infinite light direction i+k
GLfloat light_position0[4] = { 4.0 , 0.0 , 6.0 , 0.0 };
// infinite light direction -i+k
GLfloat light_position1[4] = { -1.0 , 0.0 , 1.0 , 0.0 };
// infinite light direction j-0.5k
GLfloat light_position2[4] = { 0.0 , 1.0 , -0.5 , 0.0 };

// material properties: silver
float Ambient_silver2[4] = { .2 , .05  , 0.03 , 1.0};
float Ambient_silver[4] = {.0, .0, .0, 0};
float Diffuse_silver[4] = {.50754, .50754, .50754, 1.0};
float Specular_silver[4] = {.0, .0, .0, 1.0};
float Emission_silver[4] = {.0, .0, .0, .0};
float Shininess_silver = 51.2;

float angle_x = -90, angle_y = 0;
int mouse_pos_x = 0, mouse_pos_y = 0;
int MeshID = -1;

// png screenshots
int ScreenShot = false;
int ShotRank = 0;
int MaxRank = 300;

int Trace = false;
float Zoom = 0.6;
char MeshFileName[STRINGSIZE];
char KPFileName[STRINGSIZE];
char MaterialFileName[STRINGSIZE];

int CurrentActiveKeyPoint = 0;

int main(int argc, char **argv);

void parse_mesh_obj( char *fileName );

void locate_KP_in_mesh( void );
float linear_weight( float distance , float radius , int exponent );
WeightingType weight_one_vertex( int indVertex , int indKP ,
				 float radius , int exponent , 
				 float (*pt2Function)(float,float,int) );
void weight_vertices_on_KP_in_mesh( float radius , int exponent  , 
				    float (*pt2Function)(float,float,int) );

void animate_vertices_in_mesh( void );

void displayFace( int indFace );
void make_mesh( void );

void init_scene( void );
void initGL( void );

void window_display( void );
void window_reshape(GLsizei width, GLsizei height); 
void window_key(unsigned char key, int x, int y); 
void window_mouseFunc(int button, int state, int x, int y);
void window_motionFunc(int x, int y);
void window_idle( void ); 
void window_special_key(int key, int x, int y);

void render_scene( void );
void animate_vertices_in_mesh_poses(float a );


void animate_vertices_in_mesh( void );


void compute_bone_transformations( void );
void parse_KP_obj( FILE *file );

void render_bones( void );

void parse_one_Bone_obj( FILE *file , int level , char *nextTag );
void parse_Bone_obj( FILE *file );

// local server address
  unsigned int Local_server_port = 1979;

  // local server socket
  int SocketToLocalServer = -1;
      char    message[1024];
  ///////////////////////////////
  // local server creation

  struct sockaddr_in localServAddr;

//////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////

unsigned int g_Texture[MAX_TEXTURES];						       
 //----------------------------------------
void dessinerGoal(){
int longuerbarre=7;
glPushMatrix();
{

	glDisable(GL_COLOR_MATERIAL) ;
	glEnable(GL_LIGHTING) ;
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, Diffuse_silver) ;
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, Specular_silver) ;
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, Ambient_silver) ;
	glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, &Shininess_silver) ;
    	//glColor3ub(255,255,255);
	glRotated(180,0,1,0);
	glTranslated(5,5 + 30,-3);

	//glRotated(90,0,1,0);

  //glColor3f( 0.0, 0.0, 1.0 );
  GLUquadricObj *quadratic1;
  quadratic1 = gluNewQuadric();
  gluCylinder(quadratic1,0.1,0.1,0.5*longuerbarre+2,32,32);

  glTranslated(-1.375*longuerbarre,0,0);

  //glColor3f( 0.0, 0.0, 1.0 );
  GLUquadricObj *quadratic;
  quadratic = gluNewQuadric();
  gluCylinder(quadratic,0.1,0.1,0.5*longuerbarre+2,32,32);

   glTranslated(-0.1,0,0.0);
   glRotated(90,0,1,0);

  //glColor3f( 0.0, 0.0, 1.0 );
  GLUquadricObj *quadratic2;
  quadratic2 = gluNewQuadric();
  gluCylinder(quadratic2,0.1,0.1,1.4*longuerbarre,32,32);
}
    glColor3ub(255,255,255);
    glPopMatrix();
}

int angle = 0;
#define BALLON 2
int i = 0;
float r = .5;

void drawBitmapText(char *string,float x,float y,float z) 
{
  char *c;
  glRasterPos3f(x, y,z);
  for (c=string; *c != '\0'; c++) 
  {
    glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24, *c);
  }
}

void dessinerCaisse(){
	glPushMatrix();{  
		glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);  
		glEnable(GL_DEPTH_TEST);
		glEnable(GL_TEXTURE_2D);
		int w=10;
		int h=10;
		glRotated(90,0,0,1);
		glTranslated(12,0,0);
	    		glPushMatrix();{
				glBindTexture(GL_TEXTURE_2D,g_Texture[0]);
				glRotated(90,0,0,1);
				glBegin(GL_QUADS);
					glTexCoord2d(0,0);
					glVertex3d(-w,w,-h);
					glTexCoord2d(1,0);
					glVertex3d(w,w,-h);
					glTexCoord2d(1,1);
					glVertex3d(w,w,h);
					glTexCoord2d(0,1);
					glVertex3d(-w,w,h);
				glEnd();

			}
			glPopMatrix();
		glDisable(GL_DEPTH_TEST);
		glDisable(GL_TEXTURE_2D);
	}
	glPopMatrix();
}
//-------------------------------------------------------------------

//////////////////////////////////////////////////////////////////
// MATHEMATICAL LIBRARY
//////////////////////////////////////////////////////////////////

// 3D VECTOR
class Vector {
public:
  float x, y, z;
  Vector( void ) {
    init();
  };
  Vector( float x0 , float y0 , float z0 ) {
    init( x0 , y0 , z0 );
  };
  ~Vector( void ) {
  };
  void init( void ) {
    x = 0;
    y = 0;
    z = 0;
  };
  void init( float x0 , float y0 , float z0 ) {
    x = x0;
    y = y0;
    z = z0;
  };
  void normalize( void ) {
    if( x == 0 && y == 0 && z == 0 ) {
      return;
    }
    float norm = 1.0 / sqrt( x*x + y*y + z*z );
    x *= norm;
    y *= norm;
    z *= norm;
  }
  // 1 vector
  float prodScal( Vector &v2 ) {
    return x * v2.x + y * v2.y + z * v2.z;
  }
  // average
  void averageVectors( Vector *vectors , int nbVectors ) {
    x = 0; y = 0; z = 0;
    if( nbVectors <= 0 ) {
      return;
    }
    for( int ind = 0 ; ind < nbVectors ; ind++ ) {
      x += vectors[ ind ].x;
      y += vectors[ ind ].y;
      z += vectors[ ind ].z;
    }
    float inv = 1.0 / (float)nbVectors;
    x *= inv;
    y *= inv;
    z *= inv;
  }
  void operator*=(double d) {
    x *= d;
    y *= d;
    z *= d;
  }
  void operator/=(double d) {
    x /= d;
    y /= d;
    z /= d;
  }
  void operator+=(Vector& v) {
    x += v.x;
    y += v.y;
    z += v.z;
  }
  int operator==(Vector& v) {
    return((x == v.x) && (y == v.y) && (z == v.z));
  }
  float norm(void) {
    return sqrt(x*x + y*y + z*z);
  }
  float norm2(void) {
    return (x*x + y*y + z*z);
  }
};

// 3D POINT
class Point {
public:
  float      x, y, z;
  Point( void ) {
    init();
  };
  Point( float x0 , float y0  , float z0  ) {
    x = x0;
    y = y0;
    z = z0;
  };
  ~Point( void ) {
  };
  void init( void ) {
    x = 0;
    y = 0;
    z = 0;
  };
  void operator=(Point& v) {
    x = v.x;
    y = v.y;
    z = v.z;
  }
  void operator*=(double d) {
    x *= d;
    y *= d;
    z *= d;
  }
  void operator/=(double d) {
    x /= d;
    y /= d;
    z /= d;
  }
  void operator+=(Vector& v) {
    x += v.x;
    y += v.y;
    z += v.z;
  }
  void operator+=(Point& v) {
    x += v.x;
    y += v.y;
    z += v.z;
  }
  void operator*(double f) {
    x = f * x;
    y = f * y;
    z = f * z;
  }
  int operator==(Point& v) {
    return((x == v.x) && (y == v.y) && (z == v.z));
  }
  float distance(Point& p) {
    float dx, dy, dz;
    dx = p.x - x;
    dy = p.y - y;
    dz = p.z - z;
    return sqrt(dx*dx + dy*dy + dz*dz);
  }
  void product( Point &p , float * matrixValue ) {
    x = p.x * matrixValue[0] + p.y * matrixValue[4] 
      + p.z * matrixValue[8] + matrixValue[12];
    y = p.x * matrixValue[1] + p.y * matrixValue[5]
      + p.z * matrixValue[9] + matrixValue[13];
    z = p.x * matrixValue[2] + p.y * matrixValue[6] 
      + p.z * matrixValue[10] + matrixValue[14];
  }
};

//////////////////////////////////////////////////////////////////
// VERTICES WITH 2 WEIGHTING SCHEMES
//////////////////////////////////////////////////////////////////

class Vertex {
public:
  // initial coordinates
  Point location;
  // coordinates after KeyPoint transformations
  Point curLocation;
  // 4 weights on keypoints (other vertices in the mesh)
  float      wKP[MAXKPWEIGHTS];
  // 4 indices of the weighted keypoints (other vertices in the mesh)
  int        indKP[MAXKPWEIGHTS];
  // boolean working variable to memorize weighting
  bool       weighted;

  //******************** NEW **********************
  // 4 weights on bones
  float      wBones[MAXBONEWEIGHTS];
  // 4 indices of the weighted bones
  int        indBones[MAXBONEWEIGHTS];
  // boolean working variable to memorize animation
  bool       updated;
  //******************** NEW **********************

  Vertex( void ) {
    location.init();
    curLocation.init();
    for( int ind = 0 ; ind < MAXKPWEIGHTS ; ind++ ) {
      wKP[ind] = 0.0;
      indKP[ind] = -1;
    }
    for( int ind = 0 ; ind < MAXBONEWEIGHTS ; ind++ ) {
      wBones[ind] = 0.0;
      indBones[ind] = -1;
    }
    weighted = false;
  };
  ~Vertex( void ) {
  };
};

void animate_one_vertex( int indMesh , int indVertex , Vertex *ptVertex );
class Face {
public:
  int indVertex1; 
  int indVertex2; 
  int indVertex3;
  int indNormal1; 
  int indNormal2; 
  int indNormal3;
  Face( void ) {
    indNormal1 = -1;
    indNormal2 = -1;
    indNormal3 = -1;
    indVertex1 = -1;
    indVertex2 = -1;
    indVertex3 = -1;
  }
  ~Face( void ) {
  };
};

class KP {
public:
  // keypoint ID (also reported in the vertex)
  char *id;
  // initial coordinates
  Point location;
  //  index of the mesh
  int indMesh;
  //  index of vertex
  int indVertex;
  // current translation
  Vector translation;

  KP( void ) {
    id = new char[STRINGSIZE];
    location.init();
    indMesh = -1;
    indVertex = -1;
    translation.init();
  }
  ~KP( void ) {
    delete [] id;
  }
};

class Mesh {
public:
  char *id;
  char *matId;
  int indFaceIni;
  int indFaceEnd;
  int nbKPs;
  Mesh( void ) {
    id = new char[STRINGSIZE];
    id[0] = 0;
    matId = new char[STRINGSIZE];
    matId[0] = 0;
    indFaceIni = 0;
    indFaceEnd = 0;
    nbKPs = 0;
  }
  ~Mesh( void ) {
    delete [] id;
    delete [] matId;
  }
};

//******************** NEW **********************
//
//                3D QUATERNION
//
//******************** NEW **********************
class Quaternion {
public:
  float      w;
  Vector     im;
  Vector     axe;
  float      angle;
  Quaternion( void ) {
    init();
  };
  ~Quaternion( void ) {
  };
  void init( void ) {
    im.init();
    w = 1;
    axe.init();
    angle = 0;
  };
  void operator=(Quaternion& v) {
    im = v.im;
    w = v.w;
  }
  void operator*=(double d) {
    im *= d;
    w *= d;
  }
  void operator+=(Quaternion& v) {
    im += v.im;
    w += v.w;
  }
  int operator==(Quaternion& v) {
    return((im == v.im) && (w == v.w));
  }
  ////////////////////////////////////////////////////////////
  // quaternion to oepnGL rotation transformation
  ////////////////////////////////////////////////////////////
  void QuaternionToAngleAxis( void ) {
    //            angle= 2 * acos(q_w)
    //            x= q_x / scale
    //            y= q_y / scale
    //            z= q_z / scale
    // where scale = sqrt (x2 + y2 + z2)    
    float length = im.norm2();
    if( length > 0.0000001 ) {
      angle = acos( w ) * 360.0 / M_PI;
      float scale = 1.0 / sqrt( length );
      axe = im;
      axe *= scale;
    }
    else {
      // null angle axis undetermined (here Ox)
      angle = 0;
      axe.x = 1.0;
      axe.y = 0.0;
      axe.z = 0.0;
    }
  }
  
  // spherical linear interpolation
  void Slerp( Quaternion *q1 , 
		    Quaternion *q2 ,
		    float coef , bool slerp_shortest_path ) {
    // Slerp(q1, q2, t) = q1*(q1^-1*q2)^t
    // if theta is the angle between q1 and q2
    //    = [q1 sin((1-t) * theta) + q2 sin (t * theta)] / sin(theta)
    
    // scalar product of q1 and q2
    float cosq1q2 = q1->w * q2->w 
      + q1->im.prodScal( q2->im );
    if( cosq1q2 > 1.0 ) {
      cosq1q2 = 1.0;
    }
    else if( cosq1q2 < -1.0 ) {
      cosq1q2 = -1.0;
    }
    
    // angle between q1 and q2
    float angle = acos( cosq1q2 );
    
    if ( fabs(angle) < 0.0001 ) {
      im = q1->im;
      w = q1->w;
      return;
    }
    // printf("cosq1q2 %f angle = %.2f \n" , cosq1q2 , angle );
    
    float sinq1q2 = sin(angle);
    float inverseSinq1q2 = 1.0 / sinq1q2;
    float coef0 = sin(coef * angle) * inverseSinq1q2;
    float coef1 = sin((1.0 - coef) * angle) * inverseSinq1q2;
    
    // http://www.kopete.org/spherical-linear-interpolation.html
    // because the covering is double (q and -q map to the same
    // rotation), the rotation path may turn either the "short way"
    // (less than 180°) or the "long way" (more than 180°). Long paths
    // can be prevented by negating one end if the dot product,
    // cosq1q2, is negative, thus ensuring that -90°<angle<90°.
    
    // Do we need to invert rotation?
    if( cosq1q2 < 0.0 && slerp_shortest_path ) {
      coef0 = -coef0;
      w = coef0 * q1->w + coef1 * q2->w;
      // same calculus for im
      im = q1->im;
      im *= coef0;
      Vector v2 = q2->im;
      v2 *= coef1;
      im += v2;
      
      // taking the opposite coefficient requires normalization
      normalize();
    }
    else {
      w = coef0 * q1->w + coef1 * q2->w;
      // same calculus for im
      im = q1->im;
      im *= coef0;
      Vector v2 = q2->im;
      v2 *= coef1;
      im += v2;
    }
    // printf( "Slerp out  x=%.3f y=%.3f z=%.3f w=%.3f \n",
    // im.x, im.y, im.z, w);
  }
  
  void normalize( void ) {
    float length = w * w + im.norm2();
    if( length > 0.0 ) {
      w /= length;
      im /= length;
    }    
  }
};

//******************** NEW **********************
//
//                MATRIX ALGEBRA
//
//******************** NEW **********************

void product_matrix(float *out , float *G , float*D) {
  out[0]  = G[0]*D[0] + G[4]*D[1] + G[8]*D[2] + G[12]*D[3];
  out[1]  = G[1]*D[0] + G[5]*D[1] + G[9]*D[2] + G[13]*D[3];
  out[2]  = G[2]*D[0] + G[6]*D[1] + G[10]*D[2] + G[14]*D[3];
  out[3]  = G[3]*D[0] + G[7]*D[1] + G[11]*D[2] + G[15]*D[3];
  
  out[4]  = G[0]*D[4] + G[4]*D[5] + G[8]*D[6] + G[12]*D[7];
  out[5]  = G[1]*D[4] + G[5]*D[5] + G[9]*D[6] + G[13]*D[7];
  out[6]  = G[2]*D[4] + G[6]*D[5] + G[10]*D[6] + G[14]*D[7];
  out[7]  = G[3]*D[4] + G[7]*D[5] + G[11]*D[6] + G[15]*D[7];
  
  out[8]  = G[0]*D[8] + G[4]*D[9] + G[8]*D[10] + G[12]*D[11];
  out[9]  = G[1]*D[8] + G[5]*D[9] + G[9]*D[10] + G[13]*D[11];
  out[10] = G[2]*D[8] + G[6]*D[9] + G[10]*D[10] + G[14]*D[11];
  out[11] = G[3]*D[8] + G[7]*D[9] + G[11]*D[10] + G[15]*D[11];
  
  out[12] = G[0]*D[12] + G[4]*D[13] + G[8]*D[14] + G[12]*D[15];
  out[13] = G[1]*D[12] + G[5]*D[13] + G[9]*D[14] + G[13]*D[15];
  out[14] = G[2]*D[12] + G[6]*D[13] + G[10]*D[14] + G[14]*D[15];
  out[15] = G[3]*D[12] + G[7]*D[13] + G[11]*D[14] + G[15]*D[15];
}

GLboolean invert_matrix(const GLfloat * m, GLfloat (&out)[16]) {
/* NB. OpenGL Matrices are COLUMN major. */
#define SWAP_ROWS(a, b) { GLfloat *_tmp = a; (a)=(b); (b)=_tmp; }
#define MAT(m,r,c) (m)[(c)*4+(r)]

   GLfloat wtmp[4][8];
   GLfloat m0, m1, m2, m3, s;
   GLfloat *r0, *r1, *r2, *r3;

   r0 = wtmp[0], r1 = wtmp[1], r2 = wtmp[2], r3 = wtmp[3];

   r0[0] = MAT(m, 0, 0), r0[1] = MAT(m, 0, 1),
      r0[2] = MAT(m, 0, 2), r0[3] = MAT(m, 0, 3),
      r0[4] = 1.0, r0[5] = r0[6] = r0[7] = 0.0,
      r1[0] = MAT(m, 1, 0), r1[1] = MAT(m, 1, 1),
      r1[2] = MAT(m, 1, 2), r1[3] = MAT(m, 1, 3),
      r1[5] = 1.0, r1[4] = r1[6] = r1[7] = 0.0,
      r2[0] = MAT(m, 2, 0), r2[1] = MAT(m, 2, 1),
      r2[2] = MAT(m, 2, 2), r2[3] = MAT(m, 2, 3),
      r2[6] = 1.0, r2[4] = r2[5] = r2[7] = 0.0,
      r3[0] = MAT(m, 3, 0), r3[1] = MAT(m, 3, 1),
      r3[2] = MAT(m, 3, 2), r3[3] = MAT(m, 3, 3),
      r3[7] = 1.0, r3[4] = r3[5] = r3[6] = 0.0;

   /* choose pivot - or die */
   if (fabs(r3[0]) > fabs(r2[0]))
      SWAP_ROWS(r3, r2);
   if (fabs(r2[0]) > fabs(r1[0]))
      SWAP_ROWS(r2, r1);
   if (fabs(r1[0]) > fabs(r0[0]))
      SWAP_ROWS(r1, r0);
   if (0.0 == r0[0])
      return GL_FALSE;

   /* eliminate first variable     */
   m1 = r1[0] / r0[0];
   m2 = r2[0] / r0[0];
   m3 = r3[0] / r0[0];
   s = r0[1];
   r1[1] -= m1 * s;
   r2[1] -= m2 * s;
   r3[1] -= m3 * s;
   s = r0[2];
   r1[2] -= m1 * s;
   r2[2] -= m2 * s;
   r3[2] -= m3 * s;
   s = r0[3];
   r1[3] -= m1 * s;
   r2[3] -= m2 * s;
   r3[3] -= m3 * s;
   s = r0[4];
   if (s != 0.0) {
      r1[4] -= m1 * s;
      r2[4] -= m2 * s;
      r3[4] -= m3 * s;
   }
   s = r0[5];
   if (s != 0.0) {
      r1[5] -= m1 * s;
      r2[5] -= m2 * s;
      r3[5] -= m3 * s;
   }
   s = r0[6];
   if (s != 0.0) {
      r1[6] -= m1 * s;
      r2[6] -= m2 * s;
      r3[6] -= m3 * s;
   }
   s = r0[7];
   if (s != 0.0) {
      r1[7] -= m1 * s;
      r2[7] -= m2 * s;
      r3[7] -= m3 * s;
   }

   /* choose pivot - or die */
   if (fabs(r3[1]) > fabs(r2[1]))
      SWAP_ROWS(r3, r2);
   if (fabs(r2[1]) > fabs(r1[1]))
      SWAP_ROWS(r2, r1);
   if (0.0 == r1[1])
      return GL_FALSE;

   /* eliminate second variable */
   m2 = r2[1] / r1[1];
   m3 = r3[1] / r1[1];
   r2[2] -= m2 * r1[2];
   r3[2] -= m3 * r1[2];
   r2[3] -= m2 * r1[3];
   r3[3] -= m3 * r1[3];
   s = r1[4];
   if (0.0 != s) {
      r2[4] -= m2 * s;
      r3[4] -= m3 * s;
   }
   s = r1[5];
   if (0.0 != s) {
      r2[5] -= m2 * s;
      r3[5] -= m3 * s;
   }
   s = r1[6];
   if (0.0 != s) {
      r2[6] -= m2 * s;
      r3[6] -= m3 * s;
   }
   s = r1[7];
   if (0.0 != s) {
      r2[7] -= m2 * s;
      r3[7] -= m3 * s;
   }

   /* choose pivot - or die */
   if (fabs(r3[2]) > fabs(r2[2]))
      SWAP_ROWS(r3, r2);
   if (0.0 == r2[2])
      return GL_FALSE;

   /* eliminate third variable */
   m3 = r3[2] / r2[2];
   r3[3] -= m3 * r2[3], r3[4] -= m3 * r2[4],
      r3[5] -= m3 * r2[5], r3[6] -= m3 * r2[6], r3[7] -= m3 * r2[7];

   /* last check */
   if (0.0 == r3[3])
      return GL_FALSE;

   s = 1.0 / r3[3];		/* now back substitute row 3 */
   r3[4] *= s;
   r3[5] *= s;
   r3[6] *= s;
   r3[7] *= s;

   m2 = r2[3];			/* now back substitute row 2 */
   s = 1.0 / r2[2];
   r2[4] = s * (r2[4] - r3[4] * m2), r2[5] = s * (r2[5] - r3[5] * m2),
      r2[6] = s * (r2[6] - r3[6] * m2), r2[7] = s * (r2[7] - r3[7] * m2);
   m1 = r1[3];
   r1[4] -= r3[4] * m1, r1[5] -= r3[5] * m1,
      r1[6] -= r3[6] * m1, r1[7] -= r3[7] * m1;
   m0 = r0[3];
   r0[4] -= r3[4] * m0, r0[5] -= r3[5] * m0,
      r0[6] -= r3[6] * m0, r0[7] -= r3[7] * m0;

   m1 = r1[2];			/* now back substitute row 1 */
   s = 1.0 / r1[1];
   r1[4] = s * (r1[4] - r2[4] * m1), r1[5] = s * (r1[5] - r2[5] * m1),
      r1[6] = s * (r1[6] - r2[6] * m1), r1[7] = s * (r1[7] - r2[7] * m1);
   m0 = r0[2];
   r0[4] -= r2[4] * m0, r0[5] -= r2[5] * m0,
      r0[6] -= r2[6] * m0, r0[7] -= r2[7] * m0;

   m0 = r0[1];			/* now back substitute row 0 */
   s = 1.0 / r0[0];
   r0[4] = s * (r0[4] - r1[4] * m0), r0[5] = s * (r0[5] - r1[5] * m0),
      r0[6] = s * (r0[6] - r1[6] * m0), r0[7] = s * (r0[7] - r1[7] * m0);

   MAT(out, 0, 0) = r0[4];
   MAT(out, 0, 1) = r0[5], MAT(out, 0, 2) = r0[6];
   MAT(out, 0, 3) = r0[7], MAT(out, 1, 0) = r1[4];
   MAT(out, 1, 1) = r1[5], MAT(out, 1, 2) = r1[6];
   MAT(out, 1, 3) = r1[7], MAT(out, 2, 0) = r2[4];
   MAT(out, 2, 1) = r2[5], MAT(out, 2, 2) = r2[6];
   MAT(out, 2, 3) = r2[7], MAT(out, 3, 0) = r3[4];
   MAT(out, 3, 1) = r3[5], MAT(out, 3, 2) = r3[6];
   MAT(out, 3, 3) = r3[7];

   return GL_TRUE;

#undef MAT
#undef SWAP_ROWS
}

//////////////////////////////////////////////////////////////////
//
//******************** NEW **********************
//
//               SKELETON ANIMATION
//
//******************** NEW **********************
//
//////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////
// SKELETON
//////////////////////////////////////////////////////////////////

// BONE

class Bone {
public:
  // bone ID
  char *id;
  // bone length (along y axis)
  float length;
  // parent bone if there's one
  Bone *daughterBone;
  Bone *parentBone;
  Bone *sisterBone;

  // translation
  Vector translation;
  // initial rotation
  Quaternion initialRotation;
  // animation rotation
  Quaternion animationRotation;

  // initial translation matrix computed from translation vector
  float      boneInitialTranslationMatrix[16];

  // initial and animation rotation computed from axis and angle
  float      boneAnimationRotationMatrix[16];
  float      boneInitialRotationMatrix[16];

  // joint Transformation Matrices (initial and current)
  float      initialJointTransformation[16];
  float      currentJointTransformation[16];

  // vertex animation Matrix expressed in the mesh local coordinates
  float      vertexAnimationMatrix[NBMAXMESH][16];

  Bone( void ) {
    id = new char[STRINGSIZE];
    length = 0;
    daughterBone = NULL;
    parentBone = NULL;
    sisterBone = NULL;

    // initial translation = Identity
    translation.init();
    memset( boneInitialTranslationMatrix , 0 , 16 * sizeof( float ) );
    boneInitialTranslationMatrix[0] = boneInitialTranslationMatrix[5] = boneInitialTranslationMatrix[10] = boneInitialTranslationMatrix[15] = 1.0;

    // initial rotation = Identity
    initialRotation.init();
    initialRotation.QuaternionToAngleAxis();
    memset( boneInitialRotationMatrix , 0 , 16 * sizeof( float ) );
    boneInitialRotationMatrix[0] = boneInitialRotationMatrix[5] = boneInitialRotationMatrix[10] = boneInitialRotationMatrix[15] = 1.0;

    // initial rotation animation = Identity
    animationRotation.init();
    animationRotation.QuaternionToAngleAxis();
    memset( boneAnimationRotationMatrix , 0 , 16 * sizeof( float ) );
    boneAnimationRotationMatrix[0] = boneAnimationRotationMatrix[5] = boneAnimationRotationMatrix[10] = boneAnimationRotationMatrix[15] = 1.0;

    for( int ind = 0 ; ind < NBMAXMESH ; ind++ ) {
      // initial vertex animation = Identity
      memset( vertexAnimationMatrix[ind] , 0 , 16 * sizeof( float ) );
      vertexAnimationMatrix[ind][0] = vertexAnimationMatrix[ind][5] 
	= vertexAnimationMatrix[ind][10] = vertexAnimationMatrix[ind][15] 
	= 1.0;
    }

    // initial joint transformation = Identity
    memset( initialJointTransformation , 0 , 16 * sizeof( float ) );
    initialJointTransformation[0] = initialJointTransformation[5] = initialJointTransformation[10] = initialJointTransformation[15] = 1.0;

    // current joint transformation = Identity
    memset( currentJointTransformation , 0 , 16 * sizeof( float ) );
    currentJointTransformation[0] = currentJointTransformation[5] = currentJointTransformation[10] = currentJointTransformation[15] = 1.0;
  }

  // graphical bone drawing: two orthogonal isocele triangles 
  void drawBone( void ) {
    glEnable( GL_COLOR_MATERIAL );
    glColor3f( .0 , .1  , 0.2 );
    glBegin( GL_QUADS );
    {
      glVertex3f( 0.02 , 0 , 0.02 );
      glVertex3f( -0.02 , 0 , -0.02 );
      glVertex3f( 0 , length , 0 );
      glVertex3f( 0 , length , 0 );
      glVertex3f( -0.02 , 0 , 0.02 );
      glVertex3f( 0.02 , 0 , -0.02 );
      glVertex3f( 0 , length , 0 );
      glVertex3f( 0 , length , 0 );
    }
    glEnd();
    glDisable( GL_COLOR_MATERIAL );
  }
  ~Bone( void ) {
    delete [] id;
  }
};

void render_one_bone( Bone *bone );
void compute_one_bone_transformation( Bone *bone );

//////////////////////////////////////////////////////////////////
// GOAL DATA
//////////////////////////////////////////////////////////////////

bool tirer = true;
map<string, int> map_bones;
float initRayon = 0.3;
float rayon = 0.3;
float delta_rayon = 0.0007;
bool resize_goal = false;
Vector Goal_Initial_Translation( -1.468, -37.702, -3 );
float interval_dir[2][2] = {
    {-6.94, 3.47}, // horizontal
    {1.17, 5.22} // vertical
};
float interval_but[2][2] = {
    {-5.96, 2.62}, // horizontal
    {1, 4.43} // vertical
};
float interval[2] = {15, 30};
float change_dir[2] = {0, 0};
float mat_Dir[2][3] = {
    { -1.468, -37.702, -3 },	// origine de la balle
    { -1.468, -.5, 2.05 }	// direction de tir
};
float goalTrans[3] = {0.00, 0.00, 0.00};
long change_speed = 1;
float speed_tir = 60;
float min_speed = -5.6;
float max_speed = 1.0;
float delta_speed = (max_speed - min_speed)/10;
float real_speed = min_speed;
int var_speed = 10;
float barre_dir[8][2] = {
    // vertical
    {4.5,-4.6},
    {4.1,-4.6},
    {4.1,-4.3},
    {4.5,-4.3},
    // horizontal
    {2.8,-5.3},
    {3.1,-5.3},
    {3.1,-5.7},
    {2.8,-5.7}
};
float coord_bar[4][2] = {
    // barre verticale
    {4.2,-3.5},
    {4.2,-5.3},
    // barre horizontale
    {4.0,-5.6},
    {1.8,-5.6}
};
long msg_index = -1;
bool draw_msg = false;
bool display_direction = true;
bool moving = false;

//////////////////////////////////////////////////////////////////
// SPHERE DISPLAY
//////////////////////////////////////////////////////////////////

class Sphere {
public:
  float R;
  int      SphereID;

  Sphere( void ) {
    R = rayon;
    SphereID =  glGenLists(1);
  }

  ~Sphere( void ) {
  }

  void sphere_position(float theta, float phi , float *p)
  {
    p[0] =   R * cos(phi) * cos(theta);
    p[1] =   R * sin(phi);
    p[2] =   R * cos(phi) * sin(theta);
  }

  void sphere_normal(float theta, float phi , float *n)
  {
    sphere_position(theta, phi , n);
  }

  void make_sphere( int stacks , int slices ) {
    int i, j;
    // sphere display list
    glNewList(SphereID, GL_COMPILE);
    {
      glEnable( GL_COLOR_MATERIAL );
      glColor3f( .2 , .05  , 0.03 );

      for(i = 0 ; i < stacks - 1 ; i++) {
	float t = i/(stacks-1.f);
	float t2;
	float phi = M_PI*t - M_PI/2;
	float phi2;

	if( i == stacks - 2 ) {
	  t2 = 1.0;
	  phi2 = M_PI/2;
	}
	else {
	  t2 = (i+1)/(stacks-1.f);
	  phi2 = M_PI*t2 - M_PI/2;
	}

	glBegin(GL_TRIANGLE_STRIP);
	{
	  for(j = 0 ; j <= slices-1 ; j++) {
	    float s;
	    float theta;

	    if( j == slices - 1 ) {
	      s = 1.0;
	      theta = 2 * M_PI;
	    }
	    else {
	      s = j/(slices-1.f);
	      theta = 2 * M_PI * s;
	    }
	  
	    float pos[3];
	    float norm[3];

	    sphere_normal (theta, phi, norm );
	    sphere_position (theta, phi, pos );
	    glNormal3fv( norm );
	    glVertex3fv( pos );

	    sphere_position (theta, phi2, pos );
	    sphere_normal (theta, phi2, norm );
	    glNormal3fv( norm );
	    glVertex3fv( pos );
	  }
	}
	glEnd();
      }
    }
    glEndList();
  }

  void draw_sphere( void ) {
    glCallList(SphereID);
  }
};

class Goal {
public:
  // goal ID
  char *id;
  // goal length (along y axis)
  float R;
  // translation
  Vector translation;
  // graphics
  Sphere *sphere;
  Goal( void ) {
    id = new char[STRINGSIZE];
    translation.init();
  }
  ~Goal( void ) {
    delete [] id;
    delete sphere;
  }
  void draw_goal( void ) {
    // Display List compiling
    sphere = new Sphere();
    sphere->make_sphere( 50 , 50 );
    glPushMatrix();
    {
      glTranslatef( translation.x , translation.y , translation.z );
      R = rayon;
      // calls rendering on the root bone (the first in the table of bones)
      sphere->draw_sphere();
    }
    glPopMatrix();
  }
};

Goal *goal;

void show_msg()
{
  glLoadIdentity();
  glEnable( GL_COLOR_MATERIAL );
  glPushMatrix();
  {
    if(msg_index == 0) {
      glColor3f(0,1,0);
      drawBitmapText("BUT !",-0.05,0.95,0);
    }
    if(msg_index == 1) {
      glColor3f(1,0,0);
      drawBitmapText("RATE",-0.05,0.95,0);
    }
    if(msg_index == 2) {
      glColor3f(1,0,0);
      drawBitmapText("BALLON TOUCHE",-0.32,0.88,0);
    }
    //glutSwapBuffers();
  }
  glPopMatrix();
  glDisable( GL_COLOR_MATERIAL );
  //glutPostRedisplay();
}

//////////////////////////////////////////////////////////////////
// SKELETON DATA
//////////////////////////////////////////////////////////////////
char     BoneFileName[STRINGSIZE];
Bone     *TabBones = NULL;
int      NbBones = 0;
int      CurrentActiveBone = 0;
//////////////////////////////////////////////////////////////////
// MESH TABLES
//////////////////////////////////////////////////////////////////

Vertex   *TabVertices = NULL;
Face     *TabFaces = NULL;
Vector   *TabNormals = NULL;
Mesh     *TabMeshes = NULL;
KP       *TabKPs = NULL;
int      NbMeshes = 0;
int      NbVertices = 0;
int      NbFaces = 0;
int      NbNormals = 0;
int      NbKPs = 0;

//////////////////// Déclaration de variables /////////////////////

bool afficher_bones = false;
int tailleTotale = sizeof(tab_Bones_KF_x);
int tailleLigne = sizeof(tab_Bones_KF_x[0]);
int nb_Bones = tailleTotale / tailleLigne;
int tailleElement = sizeof(tab_Bones_KF_x[0][0]);
int nbKF = tailleLigne / tailleElement;
float ** tab_trans;
int cte_fluid = 20;
int cmp_fin_bones = cte_fluid*nb_Bones;
int kf_index = 0;
int * tab_KF = NULL;
int nb_KF = 0;
clock_t startTime = clock();
clock_t endTime;

//////////////////////////////////////////////////////////////////
// Renvoyer des transitions élementaires de chaque bones
// pour un keyframe donné
//////////////////////////////////////////////////////////////////

float ** get_trans(int KF)
{
	float ** trans = new float * [nb_Bones];
	for(int i=0; i<nb_Bones; i++){
		trans[i] = new float[4];
	}
	for(int i=0; i<nb_Bones; i++){
		trans[i][0] = (tab_Bones_KF_a[i][KF] - TabBones[i].animationRotation.angle)/cte_fluid;
		trans[i][1] = (tab_Bones_KF_x[i][KF] - TabBones[i].animationRotation.axe.x)/cte_fluid;
		trans[i][2] = (tab_Bones_KF_y[i][KF] - TabBones[i].animationRotation.axe.y)/cte_fluid;
		trans[i][3] = (tab_Bones_KF_z[i][KF] - TabBones[i].animationRotation.axe.z)/cte_fluid;
	}
	return trans;
}

//////////////////////////////////////////////////////////////////
// Timer pour simuler la fluidité en passant
// d'une transition à l'autre
//////////////////////////////////////////////////////////////////

bool timer (float duree)
{
	endTime = clock();
	clock_t waitTime = endTime - startTime;
	waitTime = waitTime / 5000;
	if (waitTime >= duree)
	{
		startTime = endTime;
		return true;
	}
	return false;
}

//////////////////////////////////////////////////////////////////
// guardian moving data
//////////////////////////////////////////////////////////////////

float gard_move[1][4] =
	// intervals de mouvements verticals sans translation
	{2.5, 3.0, 3.94, 4.2};
float translations[1][2] = {0, 0};
int dir_trans = 0;
bool enable_tir = false;
bool ballon_arr = false;

//////////////////////////////////////////////////////////////////
// definir une direction de translation
//////////////////////////////////////////////////////////////////

void def_dir_trans() {
  switch(dir_trans) {
    // saut verticale
    case 1 : translations[0][0] = 0;
             translations[0][1] = .004;
             dir_trans = 2;
      break;
    // "tomber" après avoir sauter vers le haut
    case 2 : translations[0][0] = 0;
             translations[0][1] = -.003;
             dir_trans = 0;
      break;
    // sauter vers le coin haut gauche
    case 5 : translations[0][0] = .004;
             translations[0][1] = .003;
             dir_trans = 6;
      break;
    // "tomber" après avoir sauter vers le coin haut gauche
    case 6 : translations[0][0] = 0;
             translations[0][1] = -.006;
             dir_trans = 7;
      break;
    // se relever
    case 7 : translations[0][0] = -.002;
             translations[0][1] = .0005;
      break;
    // sauter vers l'état (7)
    case 8 : translations[0][0] = .000;
		     translations[0][1] = .004;
             dir_trans = 9;
      break;
    // "tomber" de l'état (7)
    case 9 : translations[0][0] = 0;
             translations[0][1] = -.006;
             dir_trans = 10;
      break;
    // se lever après avoir "tomber" de l'état (7)
    case 10 : translations[0][0] = -.0015;
              translations[0][1] = .0005;
              //dir_trans = 0;
      break;
    // sauter vers l'état (8)
    case 11 : translations[0][0] = .005;
		      translations[0][1] = -.002;
              dir_trans = 12;
      break;
    // se lever après avoir "tomber" de l'état (8)
    case 12 : translations[0][0] = -.0025;
		      translations[0][1] = .000;
      break;
    // sauter vers l'état (7) [droite]
    case 14 : translations[0][0] = .000;
		     translations[0][1] = .004;
             dir_trans = 15;
      break;
    // "tomber" de l'état (7) [droite]
    case 15 : translations[0][0] = 0;
             translations[0][1] = -.006;
             dir_trans = 16;
      break;
    // se lever après avoir "tomber" de l'état (7) [droite]
    case 16 : translations[0][0] = .0015;
              translations[0][1] = 0;
              dir_trans = 17;
      break;
    // ... [droite]
    case 17 : translations[0][0] = 0;
              translations[0][1] = .003;
              dir_trans = 0;
      break;
    // sauter vers le coin haut [droite]
    case 18 : translations[0][0] = -.004;
             translations[0][1] = .003;
             dir_trans = 19;
      break;
    // "tomber" après avoir sauter vers le coin haut [droite]
    case 19 : translations[0][0] = 0;
             translations[0][1] = -.006;
             dir_trans = 20;
      break;
    // se relever
    case 20 : translations[0][0] = .002;
             translations[0][1] = .002;
      break;
    // sauter vers l'état (8) [droite]
    case 21 : translations[0][0] = -.005;
		      translations[0][1] = -.002;
              dir_trans = 22;
      break;
    // se lever après avoir "tomber" de l'état (8) [droite]
    case 22 : translations[0][0] = .0015;
		      translations[0][1] = 0;
              dir_trans = 23;
      break;
    // se lever après avoir "tomber" de l'état (8) [droite]
    case 23 : translations[0][0] = .0050;
		      translations[0][1] = .003;
              dir_trans = 0;
      break;
    default : translations[0][0] = 0;
              translations[0][1] = 0;
  }
}

// initialiser le gardien

void init_gard() {
  cmp_fin_bones = cte_fluid*nb_Bones;
  tab_KF = new int[1];
  tab_KF[0] = 15;
  nb_KF = 1;
  enable_tir = true;
  tab_trans = get_trans(15);
  switch(dir_trans) {
    // se lever après (s'accroupir)
    case 3 : translations[0][0] = 0;
             translations[0][1] = .001;
             dir_trans = 0;
      break;
    // revenir après avoir s'incliner vers la gauche
    case 4 : translations[0][0] = -.0015;
             translations[0][1] = .001;
             dir_trans = 0;
      break;
    // revenir après avoir sauter vers le coin haut gauche
    case 5 : translations[0][0] = -.008;
             translations[0][1] = -.003;
             dir_trans = 0;
      break;
    // revenir après avoir sauter vers la gauche
    case 6 : translations[0][0] = -.0055;
             translations[0][1] = .004;
             dir_trans = 0;
      break;
    // revenir après avoir se relever
    case 7 : translations[0][0] = -.0035;
             translations[0][1] = .0035;
             dir_trans = 0;
      break;
    // se lever après avoir "tomber" de l'état (7)
    case 10 : translations[0][0] = 0;
              translations[0][1] = .0025;
              dir_trans = 0;
      break;
    // revenir de l'état 8 [droite]
    case 12 : translations[0][0] = -.004;
		      translations[0][1] = .003;
              dir_trans = 0;
      break;
    // revenir après avoir s'incliner vers la droite
    case 13 : translations[0][0] = .0015;
             translations[0][1] = .001;
             dir_trans = 0;
      break;
    // revenir de l'état 8 [gauche]
    case 20 : translations[0][0] = .0035;
             translations[0][1] = .002;
             dir_trans = 0;
      break;
  }
}

//////////////////////////////////////////////////////////////////
// Lecture des bones et application de la transition appropriée
//////////////////////////////////////////////////////////////////

void read_move()
{
	if(moving) {
		if(cmp_fin_bones>0){
			for (int i=0; i<nb_Bones; i++)
			{
				// chercher l'index du bone dans le TabBones
				int ind_Bones = map_bones.find(tab_Bones_ID[i])-> second;
				TabBones[ind_Bones].animationRotation.angle += tab_trans[i][0];
				TabBones[ind_Bones].animationRotation.axe.x += tab_trans[i][1];
				TabBones[ind_Bones].animationRotation.axe.y += tab_trans[i][2];
				TabBones[ind_Bones].animationRotation.axe.z += tab_trans[i][3];
				//// Translations ////
				TabBones[0].boneInitialTranslationMatrix[12] = TabBones[0].translation.x;
				TabBones[0].translation.x += translations[0][0];
				TabBones[10].boneInitialTranslationMatrix[12] = TabBones[10].translation.x;
				TabBones[10].translation.x += translations[0][0];
				TabBones[13].boneInitialTranslationMatrix[12] = TabBones[13].translation.x;
				TabBones[13].translation.x += translations[0][0];
				TabBones[0].boneInitialTranslationMatrix[14] = TabBones[0].translation.z;
				TabBones[0].translation.z += translations[0][1];
				TabBones[10].boneInitialTranslationMatrix[14] = TabBones[10].translation.z;
				TabBones[10].translation.z += translations[0][1];
				TabBones[13].boneInitialTranslationMatrix[14] = TabBones[13].translation.z;
				TabBones[13].translation.z += translations[0][1];
				for( int ind = 0 ; ind < NbVertices ; ind++ ) {
					TabVertices[ind].location.x += translations[0][0];
					TabVertices[ind].location.z += translations[0][1];
				}
				///////////////////////////////////////
				float vec[3];
				vec[0] = TabBones[i].animationRotation.axe.x;
				vec[1] = TabBones[i].animationRotation.axe.y;
				vec[2] = TabBones[i].animationRotation.axe.z;
				// transforms the rotation in a 4x4 matrix
				urot_about_axis_f( TabBones[i].boneAnimationRotationMatrix, TabBones[i].animationRotation.angle, vec);
				animate_vertices_in_mesh();
				//glDeleteLists(MeshID, 1);
				make_mesh();
				//glPopMatrix();
				//glutPostRedisplay();
				cmp_fin_bones--;
			}
		}
		// animation élementaire terminée
		else if (cmp_fin_bones==0)
		{
			if (enable_tir) {
				tirer = true;
				enable_tir = false;
				kf_index = 0;
			}
			def_dir_trans();
			if (kf_index < nb_KF-1) {
				kf_index++;
				cmp_fin_bones = cte_fluid*nb_Bones;
			}
			else {
				if (ballon_arr) {
					init_gard();
				}
				return;
			}
			tab_trans = get_trans(tab_KF[kf_index]);
		}
	}
}

void appel_continu()
{
	glutPostRedisplay();
	read_move();
	//window_idle();
}

//////////////////////////////////////////////////////////////////
// barre de vitesse de tir
//////////////////////////////////////////////////////////////////

void scroll_speed()
{
  if(change_speed == 0) return;
  if(speed_tir <= 30) return;
  if(change_speed == 1 && real_speed < max_speed)
  {
    real_speed += delta_speed;
    speed_tir -= var_speed;
  }
  if(change_speed == -1 && real_speed > min_speed)
  {
    real_speed -= delta_speed;
    speed_tir += var_speed;
  }
  change_speed = 0;
}

void show_scroll()
{
  scroll_speed();
  glBegin(GL_QUADS);
    glColor3d(1,0.1,0);
    glVertex3f(-7.1,-2.702,real_speed);
    glVertex3f(-6.8,-2.702,real_speed);
    glVertex3f(-6.8,-2.702,min_speed);
    glVertex3f(-7.1,-2.702,min_speed);
  glEnd();
}

//////////////////////////////////////////////////////////////////
// barres de directions de tir
//////////////////////////////////////////////////////////////////

void scroll_dir()
{
  // horizontal
  if(change_dir[0] != 0)
  {
    float delta = interval[0];
    float delta_vert = (coord_bar[2][0] - coord_bar[3][0])/delta;
    if(change_dir[0] == -1)
      delta_vert *= -1;
    if(barre_dir[4][0] <= (coord_bar[2][0]-2*delta_vert) && barre_dir[5][0] >= (coord_bar[3][0]-2*delta_vert))
    {
      barre_dir[4][0] += delta_vert;
      barre_dir[5][0] += delta_vert;
      barre_dir[6][0] += delta_vert;
      barre_dir[7][0] += delta_vert;
    }
    change_dir[0] = 0;
  }
  // vertical
  if(change_dir[1] != 0)
  {
    float delta = interval[1];
    float delta_hor = (coord_bar[0][1] - coord_bar[1][1])/delta;
    if(change_dir[1] == -1)
      delta_hor *= -1;
    if(barre_dir[2][1] <= (coord_bar[0][1]-2*delta_hor) && barre_dir[1][1] >= (coord_bar[1][1]-2*delta_hor))
    {
      barre_dir[0][1] += delta_hor;
      barre_dir[1][1] += delta_hor;
      barre_dir[2][1] += delta_hor;
      barre_dir[3][1] += delta_hor;
    }
    change_dir[1] = 0;
  }
}

void show_directions()
{
  scroll_dir();
  glBegin(GL_QUADS);
  {
    glColor3f(0,0.2,0.2);
    // barre de defilement verticale
    glVertex3f(coord_bar[0][0],-2.702,coord_bar[0][1]);
    glVertex3f(coord_bar[1][0],-2.702,coord_bar[1][1]);
    glVertex3f(coord_bar[1][0]+0.2,-2.702,coord_bar[1][1]);
    glVertex3f(coord_bar[0][0]+0.2,-2.702,coord_bar[0][1]);
  }
  glEnd();
  glBegin(GL_QUADS);
  {
    //glColor3f(0,0.2,0.2);
    // barre de defilement horizontale
    glVertex3f(coord_bar[2][0],-2.702,coord_bar[2][1]+0.2);
    glVertex3f(coord_bar[3][0],-2.702,coord_bar[3][1]+0.2);
    glVertex3f(coord_bar[3][0],-2.702,coord_bar[3][1]);
    glVertex3f(coord_bar[2][0],-2.702,coord_bar[2][1]);
  }
  glEnd();
  glBegin(GL_QUADS);
    glColor3d(0.3,1,1);
    // barre verticale
    glVertex3f(barre_dir[0][0],-2.702,barre_dir[0][1]);
    glVertex3f(barre_dir[1][0],-2.702,barre_dir[1][1]);
    glVertex3f(barre_dir[2][0],-2.702,barre_dir[2][1]);
    glVertex3f(barre_dir[3][0],-2.702,barre_dir[3][1]);
    // barre horizontale
    glVertex3f(barre_dir[4][0],-2.702,barre_dir[4][1]);
    glVertex3f(barre_dir[5][0],-2.702,barre_dir[5][1]);
    glVertex3f(barre_dir[6][0],-2.702,barre_dir[6][1]);
    glVertex3f(barre_dir[7][0],-2.702,barre_dir[7][1]);
  glEnd();
}

//////////////////////////////////////////////////////////////////
// MESH BONES AND KEYPOINT FILE PARSING
//////////////////////////////////////////////////////////////////

// OBJ file parsing (Alias Wavefront ASCII format)
void parse_mesh_obj( FILE *file )
{
  char    tag[512];
  char    line[512];
  char    ch;
  
  // Two comment lines
  // # Blender3D v244 OBJ File: Anime_Girl.blend
  // # www.blender3d.org
  fgets  ( line, 512, file );
  fgets  ( line, 512, file );
  
  // material file name
  fgets  ( line, 512, file );
  sscanf ( line, "%s %s", 
	   tag, 
	   MaterialFileName );
  // printf( "MaterialFileName %s %s\n" , tag , MaterialFileName );
  
  // bone file name
  fgets  ( line, 512, file );
  sscanf ( line, "%s %s", 
	   tag, 
	   BoneFileName );
  //printf( "\nParsing BoneFileName %s\n" , BoneFileName );
  
  // parses the bones
  FILE * fileBone = fopen( BoneFileName , "r" ); 
  if( !fileBone ) {
    //printf( "File %s not found, no bones defined for this mesh\n" , BoneFileName );
  }
  else {
    parse_Bone_obj( fileBone );
    fclose( fileBone );
  }
  //printf( "End of Parsing BoneFileName %s\n\n" , BoneFileName );

  //printf( "Parsing Mesh File %s\n" , MeshFileName );
  // mesh ID
  fgets  ( line, 512, file );
  sscanf ( line, "%s", tag );
  // printf( "Tag %s\n" , tag );
    
  while( strcmp( tag , "o" ) == 0 ) {
    if( NbMeshes >= NBMAXMESH ) {
      //printf( "Error: Excessive number of Meshes\n" );
      throw 0;
    }

    // mesh ID
    sscanf ( line, "%s %s", 
	     tag , TabMeshes[ NbMeshes ].id );
    //printf( "Mesh ID #%d %s\n" , NbMeshes , TabMeshes[ NbMeshes ].id );
    
    // next tag
    fgets  ( line, 512, file );
    sscanf ( line, "%s", tag );
  
    // Scan for Verts in this mesh
    int indVertexIni = NbVertices;
    while( strcmp( tag , "v" ) == 0 ) {
      if( NbVertices >= NBMAXVERTICES ) {
	//printf( "Error: Excessive number of vertices\n" );
	throw 0;
      }
      sscanf ( line, "%s %f %f %f", 
	       tag,
	       &(TabVertices[NbVertices].location.x),
	       &(TabVertices[NbVertices].location.y),
	       &(TabVertices[NbVertices].location.z) );
      TabVertices[NbVertices].curLocation = TabVertices[NbVertices].location;
      // printf( "vertex %f %f %f\n" , TabVertices[NbVertices].location.x,
      // 	      TabVertices[NbVertices].location.y,
      // 	      TabVertices[NbVertices].location.z );
      NbVertices++;
      fgets  ( line, 512, file );
      sscanf ( line, "%s", tag );
    }
    
    // Scan for Parent Bones in this mesh
    int nbBonesLoc = 0;
    int TabBonesLoc[NBMAXBONES];
    while( strcmp( tag , "bone" ) == 0 ) {
      char    tagAux[512];
      if( NbBones > 0 ) {
	char    tagAux[512];
	if( nbBonesLoc >= NBMAXBONES ) {
	  //printf( "Error: Excessive number of bones in object %s\n" , 
	//	  TabMeshes[ NbMeshes ].id  );
	  throw 0;
	}
	
	char    boneID[512];
	sscanf ( line, "%s %s", tag , boneID );
	bool bonefound = false;
	for( int indAux = 0 ; indAux < NbBones ; indAux++ ) {
	  if( strcmp( TabBones[indAux].id , boneID ) == 0 ) {
	    bonefound = true;
	    TabBonesLoc[nbBonesLoc] = indAux;
	    break;
	  }
	}
	if( !bonefound ) {
	  //printf( "Non bone parent group %s in object %s\n" , boneID , 
	//	  TabMeshes[ NbMeshes ].id );
	}
      }

      fgets  ( line, 512, file );
      sscanf ( line, "%s", tag );
      nbBonesLoc++;
    }

    // scans the weight vectors
    while( strcmp( tag , "vw" ) == 0 ) {
      int indVertex;
      int indBone;
      float w;

      sscanf ( line, "%s %d %d %f", tag , &indVertex , &indBone , &w );
      indVertex -= 1;
      indVertex += indVertexIni;
      indBone -= 1;

      // Scan for Bones in this mesh
      if( indBone >= NbBones ) {
	//printf( "Error: Incorrect bone index\nLine: %sBone index %d, Nb Bones %d, %dth weight for vertex %d\n" , 
	//	line, indBone , NbBones , nbBonesLoc + 2 , indVertex - indVertexIni + 1 );
	throw 0;
      }

      bool weightassigned = false;
      for( int indWeight = 0 ; indWeight < MAXBONEWEIGHTS ; indWeight++ ) {
	if( TabVertices[indVertex].wBones[indWeight] <= 0.0 ) {
	  TabVertices[indVertex].indBones[indWeight] 
	    = TabBonesLoc[indBone];
	  TabVertices[indVertex].wBones[indWeight] = w;
	  // printf( "vertex %d %s %f\n" , indVertex - indVertexIni + 1,
	  // 	  TabBones[TabVertices[indVertex].indBones[indWeight]].id,
	  //   TabVertices[indVertex].wBones[indWeight] );
	  weightassigned = true;
	  break;
	}
      }

      if( !weightassigned ) {
	//printf( "Error: Excessive number of bone weigths in object %s\n" , 
	//	TabMeshes[ NbMeshes ].id  );
	throw 0;
      }
      
      fgets  ( line, 512, file );
      sscanf ( line, "%s", tag );
    }
    
    // Scan for UV texture coordinates in this mesh (not used currently)
    while( strcmp( tag , "vt" ) == 0 ) {
      fgets  ( line, 512, file );
      sscanf ( line, "%s", tag );
    }

    // Scan for Norms in this mesh
    while( strcmp( tag , "vn" ) == 0 ) {
      if( NbNormals >= NBMAXVERTICES ) {
	//printf( "Error: Excessive number of vertices\n" );
	throw 0;
      }

      sscanf ( line, "%s %f %f %f", 
	       tag ,
	       &(TabNormals[NbNormals].x),
	       &(TabNormals[NbNormals].y),
	       &(TabNormals[NbNormals].z) );
      // printf( "normal %f %f %f\n" , TabNormals[NbNormals].x,
      // 	      TabNormals[NbNormals].y,
      // 	      TabNormals[NbNormals].z );
      NbNormals++;
      fgets  ( line, 512, file );
      sscanf ( line, "%s", tag );
    }
    
    // Scan for Mat in this mesh
    if( strcmp( tag , "usemtl" ) == 0 ) {
      sscanf ( line, "%s %s", 
	       tag , TabMeshes[ NbMeshes ].matId );
      // printf( "Mesh %d mat %s\n" , NbMeshes , TabMeshes[ NbMeshes ].matId );
      fgets  ( line, 512, file );
      sscanf ( line, "%s", tag );
    }
    
    TabMeshes[NbMeshes].indFaceIni = NbFaces;

    // Scan for Faces in this mesh
    while( strcmp( tag , "f" ) == 0 
	   || strcmp( tag , "usemtl" ) == 0
	   || strcmp( tag , "s" ) == 0 ) {
      if( NbFaces >= NBMAXFACES ) {
	//printf( "Error: Excessive number of faces\n" );
	throw 0;
      }

      // Scan for Mat in this mesh
      // currently only one mat per mesh
      if( strcmp( tag , "usemtl" ) == 0 ) {
	sscanf ( line, "%s", 
		 TabMeshes[ NbMeshes ].matId );
	// printf( "mat %s" , line );
      }
      // Scan for Smooth boolean in this mesh
      else if( strcmp( tag , "s" ) == 0 ) {
	sscanf ( line, "%s", tag );
      }
      // Scan for a Face in this mesh
      else {
	sscanf( line, "%s %d//%d %d//%d %d//%d", 
		tag,
		&(TabFaces[NbFaces].indVertex1),
		&(TabFaces[NbFaces].indNormal1),
		&(TabFaces[NbFaces].indVertex2),
		&(TabFaces[NbFaces].indNormal2),
		&(TabFaces[NbFaces].indVertex3),
		&(TabFaces[NbFaces].indNormal3) );
	// indices start from 1 in OBJ format
	// we make start from 0 for C++
	TabFaces[NbFaces].indVertex1--;
	TabFaces[NbFaces].indVertex2--;
	TabFaces[NbFaces].indVertex3--;
	TabFaces[NbFaces].indNormal1--;
	TabFaces[NbFaces].indNormal2--;
	TabFaces[NbFaces].indNormal3--;
	// printf( "face %d %d %d %d %d %d\n" , 
	//       TabFaces[NbFaces].indVertex1,
	//       TabFaces[NbFaces].indNormal1,
	//       TabFaces[NbFaces].indVertex2,
	//       TabFaces[NbFaces].indNormal2,
	//       TabFaces[NbFaces].indVertex3,
	//       TabFaces[NbFaces].indNormal3 );
	if( TabFaces[NbFaces].indVertex1 >= 0 
	    && TabFaces[NbFaces].indVertex2 >= 0 
	    && TabFaces[NbFaces].indVertex3 >= 0 
	    && TabFaces[NbFaces].indNormal1 >= 0 
	    && TabFaces[NbFaces].indNormal2 >= 0 
	    && TabFaces[NbFaces].indNormal3 >= 0 ) {
	  NbFaces++;
	}
      }

      if( !fgets( line, 512, file ) ) {
	TabMeshes[NbMeshes].indFaceEnd = NbFaces;
	//printf( "Mesh #%d %s Faces %d-%d\n" , NbMeshes , 
	//	TabMeshes[ NbMeshes ].id , 
	//	TabMeshes[ NbMeshes ].indFaceIni ,
	//	TabMeshes[ NbMeshes ].indFaceEnd );
	NbMeshes++;
	//printf( "End of Parsing Mesh File %s\n\n" , MeshFileName );
	return;
      }
      sscanf ( line, "%s", tag );
    }

    TabMeshes[NbMeshes].indFaceEnd = NbFaces;
    //printf( "Mesh #%d %s Faces %d-%d\n" , NbMeshes , 
	//    TabMeshes[ NbMeshes ].id , 
	 //   TabMeshes[ NbMeshes ].indFaceIni ,
	 //   TabMeshes[ NbMeshes ].indFaceEnd );
    NbMeshes++;
  }
  //printf( "End of Parsing Mesh File %s\n\n" , MeshFileName );
}
// Bone OBJ file parsing 
// (inhouse format inspired from the Alias Wavefront ASCII format)

void parse_one_Bone_obj( FILE *file , int level , char *nextTag ) {
  char    tag[512];
  char    id[512];
  char    line[512];

  strcpy( tag , nextTag );
  while( strcmp( tag , "transl" ) == 0 ) {
    // Scan for Bones in this mesh
    if( NbBones >= NBMAXBONES ) {
      //printf( "Error: Excessive number of Bones\n" );
      throw 0;
    }

    // has read the transl/ID line
    
    // stores the translation values
    fgets  ( line, 512, file );
    sscanf ( line, "%f %f %f", 
	     &(TabBones[NbBones].translation.x) , 
	     &(TabBones[NbBones].translation.y) , 
	     &(TabBones[NbBones].translation.z) );
    TabBones[NbBones].boneInitialTranslationMatrix[12] 
      = TabBones[NbBones].translation.x;
    TabBones[NbBones].boneInitialTranslationMatrix[13] 
      = TabBones[NbBones].translation.y;
    TabBones[NbBones].boneInitialTranslationMatrix[14] 
      = TabBones[NbBones].translation.z;

    // initialRotation tag
    fgets  ( line, 512, file );
    // stores the initialRotation values
    fgets  ( line, 512, file );
    sscanf ( line, "%f %f %f %f", 
	     &(TabBones[NbBones].initialRotation.w) , 
	     &(TabBones[NbBones].initialRotation.im.x) , 
	     &(TabBones[NbBones].initialRotation.im.y) , 
	     &(TabBones[NbBones].initialRotation.im.z) );
    TabBones[NbBones].initialRotation.QuaternionToAngleAxis();
    float vec[3];
    vec[0] = TabBones[NbBones].initialRotation.axe.x;
    vec[1] = TabBones[NbBones].initialRotation.axe.y;
    vec[2] = TabBones[NbBones].initialRotation.axe.z;
    urot_about_axis_f( TabBones[NbBones].boneInitialRotationMatrix , 
		       TabBones[NbBones].initialRotation.angle , vec );
  
    // bone
    fgets  ( line, 512, file );
    sscanf ( line, "%s %s", 
	     tag, TabBones[NbBones].id );
    // length
    fgets  ( line, 512, file );
    sscanf ( line, "%f", 
	     &(TabBones[NbBones].length) );
    // parent
    fgets  ( line, 512, file );
    sscanf ( line, "%s %s", 
	     tag, id );
    //printf( "Bone %s (parent: %s) (prof: %d)\n", 
	//    TabBones[NbBones].id , id , level );

    // empty line
    fgets  ( line, 512, file );

    // associates the parent bone with the current bone
    if( strcmp( id , "NULL" ) != 0 ) {
      bool parentfound = false;
      for( int ind = 0 ; ind < NbBones ; ind++ ) {
	if( strcmp( TabBones[ind].id , id ) == 0 ) {
	  TabBones[NbBones].parentBone = TabBones + ind;
	  if( !TabBones[ind].daughterBone ) {
	    TabBones[ind].daughterBone = TabBones + NbBones;
	  }
	  else {
	    Bone *currentBone = TabBones[ind].daughterBone;
	    while( currentBone->sisterBone ) {
	      currentBone = currentBone->sisterBone;
	    }
	    currentBone->sisterBone = TabBones + NbBones;
	  }
	  parentfound = true;
	  break;
	}
      }
      if( !parentfound ) {
	//printf( "Parent of bone %s (%s) not found!\n", 
	//	TabBones[NbBones].id , id  );
      }
    }
    // no parent chains with the root node
    else {
      // it is not the root node
      if( NbBones > 0 ) {
	Bone *currentBone = TabBones;
	while( currentBone->sisterBone ) {
	  currentBone = currentBone->sisterBone;
	}
	currentBone->sisterBone = TabBones + NbBones;
      }      
    }

    // next tag
    fgets  ( line, 512, file );
    sscanf ( line, "%s %s", tag, id );
    
    NbBones++;
    
    // daughter bone
    if( strcmp( tag , "transl" ) == 0 ) {
      parse_one_Bone_obj( file , level + 1 , tag );
      if( strcmp( tag , "end" ) == 0 ) { strcpy( nextTag , "end" ); return; }      

      // empty line: end of file
      fgets  ( line, 512, file );

      // if empty line: end of file
      if( !fgets  ( line, 512, file ) ) { strcpy( nextTag , "end" ); return; }      
      // non empty line: reads further (possible sister node)
      sscanf ( line, "%s %s", tag , id );
    }
    
    // end_bone tag
    else if( strcmp( tag , "bone_end" ) == 0 ) {
      // empty line
      fgets  ( line, 512, file );

      // if empty line: end of file
      if( !fgets  ( line, 512, file ) ) { strcpy( nextTag , "end" ); return; }      
      // non empty line: reads further (possible sister node)
      sscanf ( line, "%s %s", tag , id );
    }
  }
}

void parse_Bone_obj( FILE *file ) {
  char    tag[512];
  char    line[512];
  char    id[512];
  
  // One comment line
  // # Blender3D Bones File: Anime_Girl_Bones.blend
  fgets  ( line, 512, file );
  
  NbBones = 0;

  // next tag
  fgets  ( line, 512, file );
  sscanf ( line, "%s %s", tag , id );
  
  while( strcmp( tag , "transl" ) == 0 ) {
    parse_one_Bone_obj( file , 1 , tag );
  }
}
//////////////////////////////////////////////////////////////////
// INTIALIZATIONS
//////////////////////////////////////////////////////////////////

// scene initialization


void init_scene( void )
{
  TabVertices = new Vertex[ NBMAXVERTICES ];
  TabFaces = new Face[ NBMAXFACES ];
  TabNormals = new Vector[ NBMAXVERTICES ];
  TabMeshes = new Mesh[ NBMAXMESH ];
  TabKPs = new KP[ NBMAXKP ];

  //******************** NEW **********************
  TabBones = new Bone[ NBMAXBONES ];
  //******************** NEW **********************
  
  NbMeshes = 0;
  NbVertices = 0;
  NbFaces = 0;
  NbNormals = 0;
  //******************** NEW **********************
  NbBones = 0;
  //******************** NEW **********************

  // parses the mesh (obj format)
  FILE * fileMesh = fopen( MeshFileName , "r" ); 
  if( !fileMesh ) {
    //printf( "File %s not found\n" , MeshFileName );
    exit(0);
  }
  // calls bone file parsing
  parse_mesh_obj( fileMesh );
  fclose( fileMesh );

  // parses the keypoints
  FILE * fileKP = fopen( KPFileName , "r" ); 
  if( !fileKP ) {
    //printf( "File %s not found, no keypoint defined for this mesh\n" , KPFileName );
  }
  else {
    parse_KP_obj( fileKP );
    fclose( fileKP );
  }

  // Display List compiling
  MeshID =  glGenLists(1);
  TypeOfSurfaceDisplay = FULL;
  make_mesh();
  // locates the keypoints in the mesh
  locate_KP_in_mesh();
  // weights the vertices on these keypoints
  weight_vertices_on_KP_in_mesh( 0.03 , 0 , &linear_weight );
  // Goal initialization
  goal = new Goal;
  goal->translation = Goal_Initial_Translation;
  for(int i=0; i<NbBones; i++)
  {
    map_bones[TabBones[i].id] = i;
  }
}

int indexBone(char * id) {
	for (int i=0; i<NbBones; i++) {
		if (TabBones[i].id == id) {
			return i;
		}
	}
}

void productValMatrice(float v, float * matrix) {
	//for (int i=0; i<16; i++) {
		matrix[12] *= v;
		matrix[13] *= v;
		matrix[14] *= v;
	//}
}

//////////////////////////////////////////////////////////////////
// MESH ANIMATION
//////////////////////////////////////////////////////////////////

// moves each vertex according to the translation
// of the keypoints or the rotations of the bones
// to which this vertex is attached

void animate_one_vertex( int indMesh , int indVertex , Vertex *ptVertex ) {
  int        indKP;
  int        indBone;

  ptVertex->curLocation = ptVertex->location;
  if( ptVertex->indKP[0] >= 0 ) {
    for( int i = 0 ; i < MAXKPWEIGHTS ; i++ ) {
      // positive weight
      if( (indKP = ptVertex->indKP[i]) >= 0 ) {
	Vector t = TabKPs[indKP].translation;
	t *= ptVertex->wKP[i];
	ptVertex->curLocation += t;
	// 	printf( "KP %s Mesh %s Vertex %d Weight %.3f index(ini) %d\n" , 
	// 		TabKPs[indKP].id ,
	// 		TabMeshes[ TabKPs[indKP].indMesh ].id ,
	// 		indVertex , w , i );
      }
    }
  }

  // animation computation
  // each vertex is animated according to the current
  // transformations of the bones on which it is weighted
  // For each mesh, each bone stores in vertexAnimationMatrix
  // the transformation that must be applied to the vertices
  if( ptVertex->indBones[0] >= 0 ) {
    Point barycenter( 0 , 0 , 0 );
	float maxW = std::numeric_limits<float>::min();	// poid maximum du bone
	int indexB = 0;		// indice du bone ayant le poid maximum
	for( int i = 0 ; i < MAXBONEWEIGHTS ; i++ ) {
		// positive weight
		if( (indBone = ptVertex->indBones[i]) >= 0 ) {
		// chercher l'indice du bone "le plus lourd"
		//************************ TODO ***********************************
			if (maxW < ptVertex->wBones[i]) {
				maxW = ptVertex->wBones[i];
				indexB = indBone;
			}
		}
	}
	// computes the location of the transformed point
	barycenter.product(ptVertex->curLocation, TabBones[indexB].vertexAnimationMatrix[indMesh]);
	// stores the final transformation into ptVertex->curLocation 
	ptVertex->curLocation = barycenter;
  }
}

//////////////////////////////////////////////////////////////////
// SKELETON ANIMATION
//////////////////////////////////////////////////////////////////

void compute_one_bone_transformation( Bone *bone ) {
  if( !bone )
    return;
  // computes the initial and current joint transformation matrices 
  // the initial rotation and translation matrices
  float m1[ 16 ], m2[ 16 ], m3[ 16 ];
  // intermediary bone
  if( bone->parentBone ) {
    // initial joint transformation
    // could be computed only once in an independent function call
    //************************ TODO ***********************************
	product_matrix(m1,
		bone->parentBone->initialJointTransformation,
		bone->boneInitialTranslationMatrix);
	product_matrix(bone->initialJointTransformation,
		m1,
		bone->boneInitialRotationMatrix);
    // current joint transformation
    //************************ TODO ***********************************
	product_matrix(m2,
		bone->parentBone->currentJointTransformation,
		bone->boneInitialTranslationMatrix);
	product_matrix(m3,
		m2,
		bone->boneInitialRotationMatrix);
	product_matrix(bone->currentJointTransformation,
		m3,
		bone->boneAnimationRotationMatrix);
  }
  // root bone
  else {
    // initial joint transformation
    // could be computed only once in an independent function call
    //************************ TODO ***********************************
	product_matrix(bone->initialJointTransformation,
			bone->boneInitialTranslationMatrix,
			bone->boneInitialRotationMatrix);
    // current joint transformation
    //************************ TODO ***********************************
	product_matrix(m1,
		bone->boneInitialTranslationMatrix,
		bone->boneInitialRotationMatrix);
	product_matrix(bone->currentJointTransformation,
		m1,
		bone->boneAnimationRotationMatrix);
  }
  // The vertex coordinates are in the mesh local coordinate system
  // We compute the differential transformation from current joint
  // transformation and initial joint transformation for mesh skinning

  // invert the initial transformation of the joint 
  // initialJointTransformation
  //************************ TODO ***********************************
	float inverted_matrix [16];
	invert_matrix(bone->initialJointTransformation, inverted_matrix);
	for( int indMesh = 0 ; indMesh < NbMeshes ; indMesh++ ) {
		// currently the same computation is made for every mesh
		// because they are all given in the same coordinate system
		// combine with current joint transformation currentJointTransformation
		//************************ TODO ***********************************
		product_matrix(bone->vertexAnimationMatrix[indMesh], 
			bone->currentJointTransformation, 
			inverted_matrix);
	}
  // recursive call
  compute_one_bone_transformation( bone->daughterBone );

  // recursive call
  compute_one_bone_transformation( bone->sisterBone );
}

void compute_bone_transformations( void ) {
  // calls transformation computation on the root bone
  // (the first in the table of bones)
	compute_one_bone_transformation( TabBones );
}

//////////////////////////////////////////////////////////////////
// SKELETON DRAWING
//////////////////////////////////////////////////////////////////

void render_one_bone( Bone *bone ) {
  if( !bone ) 
    return;

  glPushMatrix();
  {
    //************************ TODO ***********************************
    // OpenGL transformations for drawing (initial position
    // + animation rotation)
    // utiliser bone->translation et bone->initialRotation pour
    // effectuer glTranslatef(x, y, z) et glRotatef (angle, x, y, z) 3 vecteurs de l'axe de rotation
    glTranslatef( bone->translation.x, bone->translation.y, bone->translation.z );
    glRotatef( bone->initialRotation.angle, bone->initialRotation.axe.x,
		bone->initialRotation.axe.y, bone->initialRotation.axe.z );
    glRotatef( bone->animationRotation.angle, bone->animationRotation.axe.x,
		bone->animationRotation.axe.y, bone->animationRotation.axe.z );
    // bone graphical rendering
    bone->drawBone();
    // recursive call
    render_one_bone( bone->daughterBone );
  }
  glPopMatrix();

  // recursive call
  render_one_bone( bone->sisterBone );
}

void render_bones( void ) {
  glPushMatrix();
  {
    // calls rendering on the root bone (the first in the table of bones)
    //glTranslatef( 0, 30, 0 );
    //glTranslatef(1.438, 2.702, -2.05);
    //glTranslatef(5, 0, 0);
    //glTranslatef(1.47, 0, -2.4);
    render_one_bone( TabBones );
  }
  glPopMatrix();
}

// voir si le ballon touche le gardien
bool gard_touche() {
  // les coordonnées du ballon
  float r = goal->sphere->R;
  float g_x[2] = {goal->translation.x-r, goal->translation.x+r};
  float g_y[2] = {goal->translation.y-r, goal->translation.y+r};
  float g_z[2] = {goal->translation.z-r, goal->translation.z+r};
  // parcourir tout les vertices du gardien
  bool flag = false;
  for(int ind = 0 ; ind < NbMeshes ; ind++) {
    for(int intFace = 0; intFace < NbFaces; intFace++) {
      Face item = TabFaces[intFace];
      if(((g_x[0] <= TabVertices[item.indVertex1].curLocation.x && g_x[1] >= TabVertices[item.indVertex1].curLocation.x) &&
         (g_y[0] <= TabVertices[item.indVertex1].curLocation.y && g_y[1] >= TabVertices[item.indVertex1].curLocation.y) &&
         (g_z[0] <= TabVertices[item.indVertex1].curLocation.z && g_z[1] >= TabVertices[item.indVertex1].curLocation.z)) ||
         ((g_x[0] <= TabVertices[item.indVertex2].curLocation.x && g_x[1] >= TabVertices[item.indVertex2].curLocation.x) &&
         (g_y[0] <= TabVertices[item.indVertex2].curLocation.y && g_y[1] >= TabVertices[item.indVertex2].curLocation.y) &&
         (g_z[0] <= TabVertices[item.indVertex2].curLocation.z && g_z[1] >= TabVertices[item.indVertex2].curLocation.z)) ||
         ((g_x[0] <= TabVertices[item.indVertex3].curLocation.x && g_x[1] >= TabVertices[item.indVertex3].curLocation.x) &&
         (g_y[0] <= TabVertices[item.indVertex3].curLocation.y && g_y[1] >= TabVertices[item.indVertex3].curLocation.y) &&
         (g_z[0] <= TabVertices[item.indVertex3].curLocation.z && g_z[1] >= TabVertices[item.indVertex3].curLocation.z)))
         flag = true;
    }
  }
  return flag;
}

// faire bouger le ballon
void move_goal() {
  // ballon prise par le gardien
  if(gard_touche()) {
    goalTrans[0] = 0.00;
    goalTrans[1] = 0.00;
    goalTrans[2] = 0.00;
    goal->translation = Goal_Initial_Translation;
    rayon = initRayon;
    resize_goal = false;
    real_speed = min_speed + delta_speed;
    speed_tir = 50;
    msg_index = 2;
    draw_msg = true;
    ballon_arr = true;
    return;
  }
  else if(goal->translation.y >= mat_Dir[1][1]) {
    // BUT
    if((goal->translation.x <= interval_but[0][1]) &&
      (goal->translation.x >= interval_but[0][0]) &&
      (goal->translation.z <= interval_but[1][1]) &&
      (goal->translation.z >= interval_but[1][0])) {
      msg_index = 0;
    }
    else {
      msg_index = 1;
    }
    real_speed = min_speed + delta_speed;
    speed_tir = 50;
    goalTrans[0] = 0.00;
    goalTrans[1] = 0.00;
    goalTrans[2] = 0.00;
    goal->translation = Goal_Initial_Translation;
    rayon = initRayon;
    resize_goal = false;
    draw_msg = true;
    ballon_arr = true;
    return;
  }
  // mettre à jour les coordonnées du ballon
  goal->translation.x += goalTrans[0];
  goal->translation.y += goalTrans[1];
  goal->translation.z += goalTrans[2];
  if(resize_goal) {
    float r = (real_speed - min_speed)/(delta_speed);
    rayon -= delta_rayon*(r/2+1);
  }
}

GLvoid release_key(int key, int x, int y)
{
  
}

// définir une direction de tir
void user_tir()
{
  float delta_x, delta_y, delta_z;
  float d_x, d_y, d_z;
  delta_x = mat_Dir[1][0] - mat_Dir[0][0];
  delta_y = mat_Dir[1][1] - mat_Dir[0][1];
  delta_z = mat_Dir[1][2] - mat_Dir[0][2];
  d_x = delta_x/speed_tir;
  d_y = delta_y/speed_tir;
  d_z = delta_z/speed_tir;
  goalTrans[0] = d_x;
  goalTrans[1] = d_y;
  goalTrans[2] = d_z;
  change_speed = 0;
}

GLvoid up_key(unsigned char key, int x, int y)
{
  // la touhe espace
  //if(key==(char)32)
  //{
    //user_tir();
  //}
}

//////////////////////////////////////////////////////////////////
// MAIN EVENT LOOP
//////////////////////////////////////////////////////////////////

int main(int argc, char **argv) 
{
#ifndef _WIN32
  int SocketToLocalServerFlags;
  /* socket creation */
  SocketToLocalServer = socket(AF_INET, SOCK_DGRAM, 0);

  if(SocketToLocalServer < 0) {
    //printf( "Error: udp server cannot open socket to local server!\n" ); 
    throw(0);
  }
  else {
    // Read the socket's flags
    SocketToLocalServerFlags = fcntl(SocketToLocalServer, F_GETFL, 0);
    // Sets the socket's flags to non-blocking
    SocketToLocalServerFlags |= O_NONBLOCK;
    int ret = fcntl(SocketToLocalServer, F_SETFL, SocketToLocalServerFlags );
    if(ret < 0) {
      ///printf(  "Error: udp server cannot set flag to non-blocking: %s!" , 
	//       strerror(errno) );
      throw(0);
    }
 }
#else
  // socket creation
  SocketToLocalServer = socket(AF_INET, SOCK_DGRAM, IPPROTO_UDP);

  if(SocketToLocalServer < 0) {
    //printf( "Error: udp server cannot open socket to local server!\n" ); 
    throw(0);
  }
  else {
    // Read the socket's flags
    unsigned long onoff=1;
    
    if (ioctlsocket(SocketToLocalServer, FIONBIO, &onoff) != 0){
      //printf(  "Error: udp server cannot set flag to non-blocking: %s!" , 
	//       strerror(errno) );
      throw(0);
    }
  }
#endif
    
  // bind local server port
  localServAddr.sin_family = AF_INET;
  localServAddr.sin_addr.s_addr = htonl(INADDR_ANY);
  localServAddr.sin_port = htons(Local_server_port);
  
  int rc = bind (SocketToLocalServer, 
		 (struct sockaddr *) &localServAddr,sizeof(localServAddr));
  if(rc < 0) {
    //printf( "Error: cannot bind locat port number %d!" , Local_server_port );
    throw(0);
  }
  //printf("udp server: waiting for data on port UDP %u\n", Local_server_port);

  if( argc >=2 ) {
    strcpy( MeshFileName , argv[1] );
  }
  else {    
    //printf( "Mesh file (football.obj):" );
    fflush( stdin );
    fgets( MeshFileName , STRINGSIZE , stdin );
    if( *MeshFileName == '\n' ) {
      strcpy( MeshFileName , "football.obj" ); ;
    }
    else {
      MeshFileName[ strlen( MeshFileName ) - 1 ] = 0; 
    }
  }
  strcpy( KPFileName , MeshFileName ); ;
  KPFileName[ strlen( KPFileName ) - 4 ] = 0; 
  strcat( KPFileName , "_KP.obj" );

  //printf( "Mesh file (%s)\n" , MeshFileName );
  //printf( "KP file (%s)\n" , KPFileName );

  // GLUT initialization
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);

  // window intialization
  glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);
  glutInitWindowPosition(0, 0);
  glutCreateWindow("...:::  Jeu videos de tir au but en 3D v1.0  :::...");

  CreateTexturePNG(&g_Texture[0], "images/grass1.png" );			
  // OpenGL and scene initialization
  initGL();
  init_scene();

  // GLUT callbacks
  // for drawing
  glutDisplayFunc(&window_display);
  // for window resize
  glutReshapeFunc(&window_reshape);
  // for keyboard events
  glutKeyboardFunc(&window_key);
  glutKeyboardUpFunc(&up_key);
  // for mouse clicks
  //glutMouseFunc(&window_mouseFunc);
  // for mouse drags
  //glutMotionFunc(&window_motionFunc);
  // for special keys
  glutSpecialFunc(&window_special_key);
  glutSpecialUpFunc(&release_key);

  //glutIdleFunc(&window_idle);

//// Itérer sur tout les bones et appliquer des transitions ////
  glutIdleFunc(&appel_continu);
  // main event loop
  glutMainLoop();

  return 1;
}

// key-point OBJ file parsing 
// (inhouse format inspired from the Alias Wavefront ASCII format)

void parse_KP_obj( FILE *file )
{
  char    tag[256];
  char    ignore[256];
  int     ignore1;
  int     ignore2;
  int     ignore3;
  char    line[256];
  char    meshID[256];
  int     indMesh;
  char    ch;
  
  // Two comment lines
  // # Anim_Girl Facial animation keypoints
  // # 
  fgets  ( line, 256, file );
  fgets  ( line, 256, file );
  
  // mesh ID
  fgets  ( line, 256, file );
  sscanf ( line, "%s", tag );
    
  NbKPs = 0;
  while( strcmp( tag , "o" ) == 0 ) {
    // mesh ID
    sscanf ( line, "%s %s", 
	     tag , meshID );

    // finds the index of the mesh associated with this KP
    indMesh = -1;
    for( int ind = 0 ; ind < NbMeshes ; ind++ ) {
      if( strcmp( meshID , TabMeshes[ ind ].id ) == 0 ) {
	indMesh = ind;
	// printf( "Mesh #%d ID %s\n" , ind , meshID );
      }
    }
    if( indMesh == -1 ) {
      //printf( "Error: KeyPoint Mesh ID [%s] not found\n" , meshID );
    }
    
    // next tag
    fgets  ( line, 256, file );
    sscanf ( line, "%s", tag );
  
    // Scan for KPs in this mesh
    int numberMeshKPs = 0;
    while( strcmp( tag , "kp" ) == 0 ) {
      if( NbKPs > NBMAXKP ) {
	//printf( "Error: Excessive number of KeyPoints\n" );
	throw 0;
      }
      sscanf ( line, "%s %s", 
	       tag, TabKPs[NbKPs].id );

      // stores the index of the mesh associated with this keypoint
      TabKPs[NbKPs].indMesh = indMesh;

      fgets  ( line, 256, file );
      sscanf ( line, "%s %f %f %f", 
	       tag,
	       &(TabKPs[NbKPs].location.x),
	       &(TabKPs[NbKPs].location.y),
	       &(TabKPs[NbKPs].location.z) );
      // printf( "vertex %f %f %f\n" , TabVertices[NbVertices].location.x,
      // 	      TabVertices[NbVertices].location.y,
      // 	      TabVertices[NbVertices].location.z );

      if( !fgets  ( line, 256, file ) ) {
	numberMeshKPs++;
	NbKPs++;
	if( indMesh >= 0 ) {
	  TabMeshes[indMesh].nbKPs = numberMeshKPs;
	}
	//printf( "Mesh #%d %s KPs %d\n" , indMesh , 
	//	TabMeshes[ indMesh ].id , 
	//	TabMeshes[ indMesh ].nbKPs );
	return;
      }

      sscanf ( line, "%s", tag );
      numberMeshKPs++;
      NbKPs++;
    }

    if( indMesh >= 0 ) {
      TabMeshes[indMesh].nbKPs = numberMeshKPs;
    }

    //printf( "Mesh #%d %s KPs %d\n" , indMesh , 
	//    TabMeshes[ indMesh ].id , 
	  //  TabMeshes[ indMesh ].nbKPs );
  }
}

//////////////////////////////////////////////////////////////////
// KEYPOINT BINDING AND VERTEX WEIGHTING
//////////////////////////////////////////////////////////////////

// for each keypoint, finds the nearest vertex in mesh
// not really used, just a check

void locate_KP_in_mesh( void ) {
  for( int indMesh = 0 ; indMesh < NbMeshes ; indMesh++ ) {
    for( int indKP = 0 ; indKP < NbKPs ; indKP++ ) {
      if( TabKPs[indKP].indMesh == indMesh ) {
	// accesses the vertices from a mesh and its faces
	float minDist = MAXFLOAT;
	int indVertexKP = -1;
	for (int indFace = TabMeshes[ indMesh ].indFaceIni ; 
	     indFace < TabMeshes[ indMesh ].indFaceEnd ; 
	     indFace++) {
	  float d;
	  if( (d = TabKPs[indKP].location.distance( 
                      TabVertices[ TabFaces[indFace].indVertex1 ].location))
	      < minDist ) {
	    indVertexKP = TabFaces[indFace].indVertex1;
	    minDist = d;
	  }
	  if( (d = TabKPs[indKP].location.distance( 
	              TabVertices[ TabFaces[indFace].indVertex2 ].location))
	      < minDist ) {
	    indVertexKP = TabFaces[indFace].indVertex2;
	    minDist = d;
	  }
	  if( (d = TabKPs[indKP].location.distance( 
	              TabVertices[ TabFaces[indFace].indVertex3 ].location))
	      < minDist ) {
	    indVertexKP = TabFaces[indFace].indVertex3;
	    minDist = d;
	  }
	}
	string kp_id = TabKPs[indKP].id;
	TabKPs[indKP].indVertex = indVertexKP;

	/*printf( "KP %s Mesh %s %f %f %f Vertex %d %f %f %f dist %f\n" , 
		TabKPs[indKP].id ,
		TabMeshes[ indMesh ].id ,
		TabKPs[indKP].location.x ,
		TabKPs[indKP].location.y ,
		TabKPs[indKP].location.z ,
		indVertexKP + 1 ,
		TabVertices[indVertexKP].location.x ,
		TabVertices[indVertexKP].location.y ,
		TabVertices[indVertexKP].location.z ,
		minDist );*/
      }
    }
  }
}

// a keypoint weighting scheme: linear weighting 

float linear_weight( float distance , float radius , int exponent ) {
  if( distance < radius ) {
    return 1 - distance / radius;
  }
  else {
    return 0.f;
  }
}

// inserer un point-clé et sont poids dans le tableau wKP[] d'un vertice

bool replace_or_insert_KP(int indVertex, int indKP, float weight) {
	bool replace = false;
	// est ce qu'il y a encore des cases vides
	bool vide = false;
	for (int i=0; i<4; i++) {
		if (TabVertices[indVertex].indKP[i] == -1) {
			vide = true;
			break;
		}
	}
	if (vide) { // inserer dans les cases vides
		for (int i=0; i<4; i++) {
			if (TabVertices[indVertex].indKP[i] == -1) {
				TabVertices[indVertex].indKP[i] = indKP;
				TabVertices[indVertex].wKP[i] = weight;
				break;
			}
		}
	}
	else { // toute les cases du tableau sont déjà renseignées
		// écraser les cases contenants les plus petites valeurs
		float minW = std::numeric_limits<float>::max();
		int index = 0;
		for (int i=0; i<4; i++) {
			if (TabVertices[indVertex].wKP[i] <= minW) {
				minW = TabVertices[indVertex].wKP[i];
				index = i;
			}
		}
		if (weight > TabVertices[indVertex].wKP[index]) {
			TabVertices[indVertex].indKP[index] = indKP;
			TabVertices[indVertex].wKP[index] = weight;
			replace = true;
		}
	}
	return replace;
}

// vertices weighting: each vertex is weighted on at most
// four keypoints
// if there are more than four keypoints with non null
// weight for this vertex, only the 4 KPs with the heighest
// weights are taken into account

WeightingType weight_one_vertex( int indVertex , int indKP ,
				 float radius , int exponent , 
				 float (*pt2Function)(float,float,int) ) {
  if( !TabVertices[ indVertex ].weighted ) {
	TabVertices[ indVertex ].weighted = true;
	// TODO
	// computes the distance and the weights
	float dx, dy, dz, distance, weight;
	dx = TabKPs[indKP].location.x - TabVertices[indVertex].location.x;
	dy = TabKPs[indKP].location.y - TabVertices[indVertex].location.y;
	dz = TabKPs[indKP].location.z - TabVertices[indVertex].location.z;
	distance = sqrt(dx*dx + dy*dy + dz*dz);
	weight = pt2Function(distance ,radius, 0);
	// inserer indKP et weight dans les tableaux indKP[] et wKP[]
	if (replace_or_insert_KP(indVertex, indKP, weight))
		return WeightSubstitution;
	else return Weighting;
    // return Weighting or WeightSubstitution
    // depending on whether the current weight has replaced
    // a previous one or not
  }
	else return NoWeighting;
}

// vertices weighting: weights all the vertices in a mesh

void weight_vertices_on_KP_in_mesh( float radius , int exponent , 
				    float (*pt2Function)(float,float,int) ) {
  for( int indMesh = 0 ; indMesh < NbMeshes ; indMesh++ ) {
    for( int indKP = 0 ; indKP < NbKPs ; indKP++ ) {
      if( TabKPs[indKP].indMesh == indMesh ) {
	int nbWeightedVertices = 0;

	// marks all the vertices as unprocessed for the current keypoint
	for( int indVertex = 0 ; indVertex < NbVertices ; indVertex++ ) {
	  TabVertices[ indVertex ].weighted = false;
	}	    

	// accesses the vertices from a mesh and its faces
	for (int indFace = TabMeshes[ indMesh ].indFaceIni ; 
	     indFace < TabMeshes[ indMesh ].indFaceEnd ; 
	     indFace++) {
	  if( weight_one_vertex( TabFaces[indFace].indVertex1 , indKP , 
				 radius , exponent , 
				 pt2Function ) == Weighting ) {
	    nbWeightedVertices++;
	  }
	  if( weight_one_vertex( TabFaces[indFace].indVertex2 , indKP , 
				 radius , exponent , 
				 pt2Function ) == Weighting ) {
	    nbWeightedVertices++;
	  }
	  if( weight_one_vertex( TabFaces[indFace].indVertex3 , indKP , 
				 radius , exponent , 
				 pt2Function ) == Weighting ) {
	    nbWeightedVertices++;
	  }
	}
	
	/*printf( "KP %s Mesh %s Nb weighted vertices %d\n" , 
		TabKPs[indKP].id ,
		TabMeshes[ TabKPs[indKP].indMesh ].id ,
		nbWeightedVertices );*/
      }
    }
  }
}

//////////////////////////////////////////////////////////////////
// MESH ANIMATION
//////////////////////////////////////////////////////////////////

// moves each vertex according to the translation
// of the keypoints to which this vertex is attached

// mesh animation mesh by mesh / face by face  / vertex by vertex

void animate_vertices_in_mesh( void ) {
  Vertex    *ptVertex = TabVertices;
  for( int indVertex = 0 ; indVertex < NbVertices ; 
       (ptVertex++ , indVertex++) ) {
    ptVertex->curLocation = ptVertex->location;
    ptVertex->updated = false;
  }
  // vertex update must be made mesh by mesh
  // because it depends on the mesh local coordinates
	compute_bone_transformations();
  for( int indMesh = 0 ; indMesh < NbMeshes ; indMesh++ ) {
    for (int i = TabMeshes[ indMesh ].indFaceIni ; 
	 i < TabMeshes[ indMesh ].indFaceEnd ; i++) {
      int indVertex = TabFaces[i].indVertex1;
      ptVertex = TabVertices + indVertex;
      if( ! ptVertex->updated ) {
	animate_one_vertex( indMesh , indVertex , ptVertex );
	ptVertex->updated = true;
      }
      indVertex = TabFaces[i].indVertex2;
      ptVertex = TabVertices + indVertex;
      if( ! ptVertex->updated ) {
	animate_one_vertex( indMesh , indVertex , ptVertex );
	ptVertex->updated = true;
      }
      indVertex = TabFaces[i].indVertex3;
      ptVertex = TabVertices + indVertex;
      if( ! ptVertex->updated ) {
	animate_one_vertex( indMesh , indVertex , ptVertex );
	ptVertex->updated = true;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////
// MESH DISPLAY
//////////////////////////////////////////////////////////////////

// should be taken from Mesh-display-modele.cpp

// face display

void make_mesh( void ) {
  glNewList(MeshID, GL_COMPILE);
  glPushMatrix();
  {
	//glTranslated(0,30,0);
    //glTranslatef(2.9, 5.402, -4.10);
    glScalef(1.6,1.6,1.6);
    //glTranslatef(1.47, 0, -2.4);
    for( int ind = 0 ; ind < NbMeshes ; ind++ ) {
      // displays wireframe + surface
	//255-239-213
	if( TypeOfSurfaceDisplay == FULL ) {
	// TODO
		glDisable(GL_COLOR_MATERIAL) ;
	glEnable(GL_LIGHTING) ;
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, Diffuse_silver) ;
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, Specular_silver) ;
	//glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, Ambient_silver) ;
	glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, &Shininess_silver) ;
	glBegin(GL_TRIANGLES) ;
	for (int intFace = 0; intFace < NbFaces; intFace++){
		Face item = TabFaces[intFace];
		if (intFace>4310 ||intFace <1530 ) 
		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, Ambient_silver) ;
		else {
			glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, Ambient_silver2) ;
		} 
		glNormal3f(TabNormals[item.indNormal1].x, TabNormals[item.indNormal1].y, TabNormals[item.indNormal1].z);
		glVertex3f(TabVertices[item.indVertex1].curLocation.x, TabVertices[item.indVertex1].curLocation.y, TabVertices[item.indVertex1].curLocation.z);
		glNormal3f(TabNormals[item.indNormal2].x, TabNormals[item.indNormal2].y, TabNormals[item.indNormal2].z);
		glVertex3f(TabVertices[item.indVertex2].curLocation.x, TabVertices[item.indVertex2].curLocation.y, TabVertices[item.indVertex2].curLocation.z);
		glNormal3f(TabNormals[item.indNormal3].x, TabNormals[item.indNormal3].y, TabNormals[item.indNormal3].z);
		glVertex3f(TabVertices[item.indVertex3].curLocation.x, TabVertices[item.indVertex3].curLocation.y, TabVertices[item.indVertex3].curLocation.z);
	}
	glEnd();
      }
      // displays wireframe
      else if( TypeOfSurfaceDisplay == MESH ) {
	// TODO
      }
      // displays surface
      else {
	// TODO
      }
    }
    }
  glPopMatrix();
  glEndList();
}

//////////////////////////////////////////////////////////////////
// INTIALIZATIONS
//////////////////////////////////////////////////////////////////

// OpenGL intialization

GLvoid initGL( void ) 
{
  // 3 light sources
  glLightfv(GL_LIGHT0, GL_AMBIENT, ambient_light0);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse_light0);
  glLightfv(GL_LIGHT0, GL_SPECULAR, specular_light0);
  glLightfv(GL_LIGHT0, GL_POSITION, light_position0);

  glLightfv(GL_LIGHT1, GL_AMBIENT, ambient_light1);
  glLightfv(GL_LIGHT1, GL_DIFFUSE, diffuse_light1);
  glLightfv(GL_LIGHT1, GL_SPECULAR, specular_light1);
  glLightfv(GL_LIGHT1, GL_POSITION, light_position1);

  glLightfv(GL_LIGHT2, GL_AMBIENT, ambient_light2);
  glLightfv(GL_LIGHT2, GL_DIFFUSE, diffuse_light2);
  glLightfv(GL_LIGHT2, GL_SPECULAR, specular_light2);
  glLightfv(GL_LIGHT2, GL_POSITION, light_position2);

  // Gouraud shading
  glShadeModel( GL_SMOOTH );

  // two side surface lighting
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE); 

  // BG color
  glClearColor( 1, 1, 1, 1 );  
      
  // z-buffer
  glEnable(GL_DEPTH_TEST);

  // lighting
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHT1);
  glEnable(GL_LIGHT2);

  // automatic normal normalization
  glEnable(GL_NORMALIZE);
}

//////////////////////////////////////////////////////////////////
// UDP CONNECTION
//////////////////////////////////////////////////////////////////

// MESSAGE PROCESSING

void processUDPMessages( void ) {

  // reads incoming UDP messages
  if( SocketToLocalServer >= 0 ) {
    char    message[1024];
    // init message buffer: Null values
    memset(message,0x0,1024);
    
    // receive message
    int n = recv(SocketToLocalServer, message, 1024, 0);
    if( n >= 0 ) {
      // scans the message and extract string & float values
      char MessID[256];
      float MessVector[4];
      int keyFrame;
      // printf( "Message size %d\n" , n );
      // printf( "Rec.: [%s]\n" , message );
      sscanf( message , "%s %d %f %f %f %f %f" , MessID , &keyFrame ,
	      MessVector , MessVector + 1 , MessVector + 2 , MessVector + 3 );
      //printf( "ID [%s] KF [%d] rot (%.3f,%.3f,%.3f,%.3f)\n" , MessID , keyFrame ,
	//      MessVector[0] , MessVector[1] , MessVector[2] , MessVector[3] );

      bool knownKP = false;
      /*for( int ind = 0 ; ind < NbKPs ; ind++ ) {
	if( strcmp( MessID , TabKPs[ind].id ) == 0 ) {
	  CurrentActiveKeyPoint = ind;
	  knownKP = true;
	  break;
	}
      }
      if( knownKP ) {
	TabKPs[CurrentActiveKeyPoint].translation.x  = MessVector[0];
	TabKPs[CurrentActiveKeyPoint].translation.y  = MessVector[1];
	TabKPs[CurrentActiveKeyPoint].translation.z  = MessVector[2]; //= boneInitialTranslationMatrix[14]
	animate_vertices_in_mesh();
	glDeleteLists(MeshID, 1);  
	make_mesh();
      }*/

      //******************** NEW **********************
      bool knownBone = false;
      for( int ind = 0 ; ind < NbBones ; ind++ ) {
	if( strcmp( MessID , TabBones[ind].id ) == 0 ) {
	  CurrentActiveBone = ind;
	  knownBone = true;
	  break;
	}
      }
      if( knownBone ) {
	TabBones[CurrentActiveBone].animationRotation.angle = MessVector[0];
	TabBones[CurrentActiveBone].animationRotation.axe.x = MessVector[1];
	TabBones[CurrentActiveBone].animationRotation.axe.y = MessVector[2];
	TabBones[CurrentActiveBone].animationRotation.axe.z = MessVector[3];
	// transforms the rotation in a 4x4 matrix
	urot_about_axis_f( TabBones[CurrentActiveBone].boneAnimationRotationMatrix ,
			   *MessVector , MessVector + 1 );

	// animates the vertices according to the bone displacements
	animate_vertices_in_mesh();
	glDeleteLists(MeshID, 1);  
	make_mesh();
      }
      //******************** NEW **********************

      if( !knownKP && !knownBone ) {
	//printf("UDP message with unknown Bone or KeyPoint (%s)\n" , MessID );
	return;
      }
	glutPostRedisplay();
    }
  }
}

//////////////////////////////////////////////////////////////////
// GLUT CALL-BACKS
//////////////////////////////////////////////////////////////////

// glut call-back: window display

GLvoid window_display( void )
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glLoadIdentity();
  glScalef( Zoom/5 , Zoom/5 , Zoom/5 );
  render_scene();

  glFlush();
}

// glut call-back: window resize

GLvoid window_reshape(GLsizei width, GLsizei height)
{  
  glViewport(0, 0, width, height);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(-1.2, 1.2, -1.2, 1.2, -10.2, 10.2);

  // toutes les transformations suivantes s\B4appliquent au mod\E8le de vue 
  glMatrixMode(GL_MODELVIEW);
}

// ballon tiré vers le milieu

void anim_gardien_0()
{
	// le gardien saute
	if(mat_Dir[1][2] > gard_move[0][3])
	{
		translations[0][0] = 0;
		translations[0][1] = -.001;
		tab_KF = new int[3];
		tab_KF[0] = 3;
		tab_KF[1] = 1;
		tab_KF[2] = 15;
		nb_KF = 3;
		// saut
		dir_trans = 1;
		tab_trans = get_trans(3);
		moving = true;
	}
	// lever les mains
	if(mat_Dir[1][2] <= gard_move[0][3] && mat_Dir[1][2] > gard_move[0][2])
	{
		tab_KF = new int[2];
		tab_KF[0] = 1;
		tab_KF[1] = 1;
		nb_KF = 2;
		tab_trans = get_trans(1);
		moving = true;
	}
	// les mains au niveau de la poitrine
	if(mat_Dir[1][2] <= gard_move[0][2] && mat_Dir[1][2] > gard_move[0][1])
	{
		tab_KF = new int[2];
		tab_KF[0] = 0;
		tab_KF[1] = 0;
		nb_KF = 2;
		tab_trans = get_trans(0);
		moving = true;
	}
	// baisser les mains
	if(mat_Dir[1][2] <= gard_move[0][1] && mat_Dir[1][2] > gard_move[0][0])
	{
		tab_KF = new int[2];
		tab_KF[0] = 2;
		tab_KF[1] = 2;
		nb_KF = 2;
		tab_trans = get_trans(2);
		moving = true;
	}
	// le gardien s'accroupie
	if(mat_Dir[1][2] <= gard_move[0][0])
	{
		translations[0][0] = 0;
		translations[0][1] = -.001;
		tab_KF = new int[2];
		tab_KF[0] = 4;
		tab_KF[1] = 4;
		nb_KF = 2;
		dir_trans = 3;
		tab_trans = get_trans(4);
		moving = true;
	}
}

// ballon tiré vers la gauche [niveau 1]

void anim_gardien_1()
{
	// le gardien saute
	if(mat_Dir[1][2] >= 2.7)
	{
		translations[0][0] = .0015;
		translations[0][1] = -.001;
		tab_KF = new int[4];
		tab_KF[0] = 5;
		tab_KF[1] = 7;
		tab_KF[2] = 8;
		tab_KF[3] = 9;
		nb_KF = 4;
		dir_trans = 8;
		tab_trans = get_trans(5);
		moving = true;
	}
	else
	{
		// s'incliner vers la gauche
		translations[0][0] = .0015;
		translations[0][1] = -.001;
		tab_KF = new int[1];
		tab_KF[0] = 5;
		nb_KF = 1;
		dir_trans = 4;
		tab_trans = get_trans(5);
		moving = true;
	}
}

// ballon tiré vers la gauche [niveau 2]

void anim_gardien_2()
{
	// le gardien saute vers le coin haut gauche
	if(mat_Dir[1][2] >= 2.7)
	{
		translations[0][0] = .0015;
		translations[0][1] = -.001;
		tab_KF = new int[4];
		tab_KF[0] = 5;
		tab_KF[1] = 6;
		tab_KF[2] = 8;
		tab_KF[3] = 9;
		nb_KF = 4;
		dir_trans = 5;
		tab_trans = get_trans(5);
		moving = true;
	}
	else
	{
		// sauter vers la gauche
		translations[0][0] = .0015;
		translations[0][1] = -.001;
		tab_KF = new int[3];
		tab_KF[0] = 5;
		tab_KF[1] = 8;
		tab_KF[2] = 9;
		nb_KF = 3;
		dir_trans = 11;
		tab_trans = get_trans(5);
		moving = true;
	}
}

// ballon tiré vers la droite [niveau 1]

void anim_gardien_left_1()
{
	// le gardien saute
	if(mat_Dir[1][2] >= 2.7)
	{
		translations[0][0] = -.0015;
		translations[0][1] = -.001;
		tab_KF = new int[4];
		tab_KF[0] = 10;
		tab_KF[1] = 11;
		tab_KF[2] = 12;
		tab_KF[3] = 13;
		nb_KF = 4;
		dir_trans = 14;
		tab_trans = get_trans(10);
		moving = true;
	}
	else
	{
		// s'incliner vers la droite
		translations[0][0] = -.0015;
		translations[0][1] = -.001;
		tab_KF = new int[1];
		tab_KF[0] = 10;
		nb_KF = 1;
		dir_trans = 13;
		tab_trans = get_trans(10);
		moving = true;
	}
}

// ballon tiré vers la droite [niveau 2]

void anim_gardien_left_2()
{
	// le gardien saute vers le coin haut droite
	if(mat_Dir[1][2] >= 2.7)
	{
		translations[0][0] = -.0015;
		translations[0][1] = -.001;
		tab_KF = new int[4];
		tab_KF[0] = 10;
		tab_KF[1] = 14;
		tab_KF[2] = 12;
		tab_KF[3] = 13;
		nb_KF = 4;
		dir_trans = 18;
		tab_trans = get_trans(10);
		moving = true;
	}
	else
	{
		// sauter vers la droite
		translations[0][0] = -.0015;
		translations[0][1] = -.001;
		tab_KF = new int[3];
		tab_KF[0] = 10;
		tab_KF[1] = 12;
		tab_KF[2] = 13;
		nb_KF = 3;
		dir_trans = 21;
		tab_trans = get_trans(10);
		moving = true;
	}
}

// FAIRE BOUGER LE GARDIEN VERS LE BALLON

void anim_gardien()
{
	cmp_fin_bones = cte_fluid*nb_Bones;
	// en fonction des intervales horizontaux
	//  -> ballon tiré au milieu
	if(mat_Dir[1][0] <= -1.3 && mat_Dir[1][0] >= -1.5) {
		anim_gardien_0();
	}
	//  -> ballon tiré vers la droite [niveau 1]
	if(mat_Dir[1][0] <= 0.1 && mat_Dir[1][0] >= -0.8) {
		anim_gardien_1();
	}
	//  -> ballon tiré vers la droite [niveau 2]
	if(mat_Dir[1][0] >= 0.5) {
		anim_gardien_2();
	}
	//  -> ballon tiré vers la gauche [niveau 1]
	if(mat_Dir[1][0] <= -2.0 && mat_Dir[1][0] >= -2.9) {
		anim_gardien_left_1();
	}
	//  -> ballon tiré vers la gauche [niveau 2]
	if(mat_Dir[1][0] <= -3.4) {
		anim_gardien_left_2();
	}
}

// glut call-back: keyboard events

GLvoid window_key(unsigned char key, int x, int y)
{
  // control of enter
  if(key==(char)13)
  {
    if(!tirer)
      return;
    tirer = false;
    ballon_arr = false;
    draw_msg = false;
    resize_goal = true;
    user_tir();
    // faire jouer le gardien
    anim_gardien();
  }
  // la touhe espace
  if(key==(char)32)
  {
    draw_msg = false;
    change_speed = 1;
  }
  switch (key) {
  case 'l': 
  case 'L':
    display_direction = !display_direction;
	break;
  case 'f': 
  case 'F':
    //if( TypeOfSurfaceDisplay != FULL ) {
      //TypeOfSurfaceDisplay = FULL; 
      //make_mesh();
      //glutPostRedisplay();
    //}
	break; 
  case '<':
    draw_msg = false;
    if(Zoom > 0.7) {
      Zoom /= 1.1;
      glutPostRedisplay();
    } else Zoom = 0.6;
    break; 
  case '>':
    draw_msg = false;
    Zoom *= 1.1;
    glutPostRedisplay();
    break;
  case '1':
    CurrentActiveKeyPoint = 0;
    break; 
  case '2':
    CurrentActiveKeyPoint = 1;
    break; 
  case '3':
    CurrentActiveKeyPoint = 2;
    break; 
  case '4':
    CurrentActiveKeyPoint = 3;
    break; 
  case '5':
    CurrentActiveKeyPoint = 4;
    break; 
  case '6':
    CurrentActiveKeyPoint = 5;
    break; 
  case '7':
    CurrentActiveKeyPoint = 6;
    break; 
  case '8':
    CurrentActiveKeyPoint = 7;
    break; 
  case '9':
    CurrentActiveKeyPoint = 8;
    break; 
  case KEY_ESC:  
    exit(1);                    
    break; 
  default:
    break;
  }     
}

// glut call-back: mouse click events

GLvoid window_mouseFunc(int button, int state, int x, int y)
{
  if (state == GLUT_DOWN && button == GLUT_LEFT_BUTTON) {
    mouse_pos_x = x;
    mouse_pos_y = y;
  }
}

// glut call-back: mouse drag events

GLvoid window_motionFunc(int x, int y)
{
  angle_x += y - mouse_pos_y;
  angle_y += x - mouse_pos_x;

  mouse_pos_x = x;
  mouse_pos_y = y;

  glutPostRedisplay();
}

// glut call-back: update function called at each frame

GLvoid window_idle( void ) 
{
   processUDPMessages();
}

// glut call-back: special key processing

#define ANIMATION_STEP    0.001
void window_special_key(int key, int x, int y)
{
  switch (key) {
  case GLUT_KEY_LEFT: // control of horizontal move
    draw_msg = false;
    if(mat_Dir[1][0] > interval_dir[0][0])
    {
      mat_Dir[1][0] -= (interval_dir[0][1]-interval_dir[0][0])/interval[0];
      change_dir[0] = -1;
    }
    break;
  case GLUT_KEY_RIGHT: // control of horizontal move
    draw_msg = false;
    if(mat_Dir[1][0] < interval_dir[0][1])
    {
      mat_Dir[1][0] += (interval_dir[0][1]-interval_dir[0][0])/interval[0];
      change_dir[0] = 1;
    }
    break;
  case GLUT_KEY_UP: // control of vertical move
    draw_msg = false;
    if(mat_Dir[1][2] < interval_dir[1][1])
    {
      mat_Dir[1][2] += (interval_dir[1][1]-interval_dir[1][0])/interval[1];
      change_dir[1] = 1;
    }
    break;
  case GLUT_KEY_DOWN: // control of vertical move
    draw_msg = false;
    if(mat_Dir[1][2] > interval_dir[1][0])
    {
      mat_Dir[1][2] -= (interval_dir[1][1]-interval_dir[1][0])/interval[1];
      change_dir[1] = -1;
    }
    break;
  case GLUT_KEY_HOME: // control of vertical move
    {
      TabKPs[CurrentActiveKeyPoint].translation.init();
      animate_vertices_in_mesh();
      glDeleteLists(MeshID, 1);  
      make_mesh();

      glutPostRedisplay();
    }
    break;
  case GLUT_KEY_END: // control of vertical move
    {
    }
    break;
  case GLUT_KEY_PAGE_UP: // control of view elevation
    {
      TabKPs[CurrentActiveKeyPoint].translation.z += ANIMATION_STEP;
      if (TabKPs[CurrentActiveKeyPoint].translation.z < 0)
	  TabKPs[CurrentActiveKeyPoint].translation.z = -TabKPs[CurrentActiveKeyPoint].translation.z;
      animate_vertices_in_mesh();
      glDeleteLists(MeshID, 1);  
      make_mesh();

      glutPostRedisplay();
    }
    break;
  case GLUT_KEY_PAGE_DOWN: // control of view elevation
    {
      TabKPs[CurrentActiveKeyPoint].translation.z -= ANIMATION_STEP;
      if (TabKPs[CurrentActiveKeyPoint].translation.z > 0)
	  TabKPs[CurrentActiveKeyPoint].translation.z = -TabKPs[CurrentActiveKeyPoint].translation.z;
      animate_vertices_in_mesh();
      glDeleteLists(MeshID, 1);  
      make_mesh();

      glutPostRedisplay();
    }
    break;
  default:
    //printf ("special key %d is not active.\n", key);
    break;
  }
}

//////////////////////////////////////////////////////////////////
// INTERACTIVE SCENE RENDERING
//////////////////////////////////////////////////////////////////

void render_scene( void )
{
  glPushMatrix();
  {
    glRotatef(angle_x, 1, 0, 0);
    glRotatef(angle_y, 0, 0, 1);
    //glTranslatef(0, 30, 0);
    dessinerCaisse();
    glEnable( GL_COLOR_MATERIAL );
    glEnable(GL_DEPTH_TEST);
    //dessinerGoal();
    // directions de tir
    glTranslatef(2.6, 2.702, -1.4);
    //glTranslatef(0, 0, 2.6);
    glCallList(MeshID);
    glScalef(1.6,1.6,1.6);
    glPushMatrix();
    {
      //glTranslatef(1.438, 2.702, -2.05);
      glBegin( GL_LINES );
      {
        if(display_direction)
        {
          // visualiser la direction choisie
          glColor3f(0,1,0);
          glVertex3f( mat_Dir[0][0], mat_Dir[0][1], mat_Dir[0][2] );
          glVertex3f( mat_Dir[1][0], mat_Dir[1][1], mat_Dir[1][2] );
        }
      }
      glEnd();
    }
    glPopMatrix();
    render_bones();
    glDisable(GL_DEPTH_TEST);
    //glCallList(BALLON);
    //glTranslatef(1.438, 2.702, -2.05);
    //******************** NEW **********************
    // no z-buffering for the bones so that they are visible whatever
    // the rendering mode
    make_mesh();
    //glutPostRedisplay();
    // dessiner une sphere
    goal->draw_goal();
    move_goal();
    glTranslatef(-.3, 0, 0.7);
    show_scroll();
    show_directions();
    if(draw_msg)
      show_msg();
    glEnable(GL_DEPTH_TEST);
  }
  glPopMatrix();

  glutSwapBuffers();
}
