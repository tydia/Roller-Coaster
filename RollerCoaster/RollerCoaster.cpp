/*
  CSCI 480
  Assignment 2
 */

#include <stdio.h>
#include <stdlib.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#include "pic.h"
#include <iostream>
#include <vector>
#include <math.h>
using namespace std;

/* for animation part. */
int animCounter = 0;
const int animDelay = 0;

/* used for using mouse click */
int g_iMenuId;
int g_iLeftMouseButton = 0;    /* 1 if pressed, 0 if not */
int g_iMiddleMouseButton = 0;
int g_iRightMouseButton = 0;

/* Define basis matrix for Catmull Rom spline*/
GLfloat s = 0.5;
GLfloat basis_mat[4][4] = {{ -s, 2-s,   s-2,  s}, 
                           {2*s, s-3, 3-2*s, -s}, 
                           { -s,   0,     s,  0},
                           {  0,   1,     0,  0}};

/* Define max line length for recuesive subdivision */
GLdouble maxlinelength = 0.001;

/* declare display lists */
GLuint splineDispList;
GLuint xSectionDispList;
GLuint groundDispList;
GLuint skyDispList;

/* for loading ground and sky textures*/
Pic * g_groundTex;
GLuint groundTexID;
Pic * g_skyTex;
GLuint skyTexID;

/* represents one control point along the spline */
struct point {
   double x;
   double y;
   double z;
};

/* declearations for subdivide */
vector <point> drawingPoints;
GLdouble umid;
vector<GLdouble> x0;
vector<GLdouble> x1;
point pt_1st;
point pt_2nd;
vector <point> tangentVecs;
vector<GLdouble> x0tangent;
vector<GLdouble> x1tangent;
point pt_1st_tangent;
point pt_2nd_tangent;

/* decalarations for 1) moving camera along spline, 2) create double rail*/
vector< vector<GLdouble> > Nvectors;
vector< vector<GLdouble> > Bvectors;
vector <GLdouble> T;
vector <GLdouble> N;
vector <GLdouble> B;
vector <GLdouble> V;

/* spline struct which contains how many control points, and an array of control points */
struct spline {
   int numControlPoints;
   struct point *points;
};

/* the spline array */
struct spline *g_Splines;

/* total number of splines */
int g_iNumOfSplines;

/* helper functions declartion */
vector<GLdouble> getCatmullRom (GLdouble u, GLfloat (*basis)[4][4], point *p0, point *p1, point *p2, point *p3);
vector<GLdouble> getUnitTangentVec (GLdouble u, GLfloat (*basis)[4][4], point *p0, point *p1, point *p2, point *p3);
void subdivide(GLdouble u0, GLdouble u1, GLdouble maxlinelength, point *p0, point *p1, point *p2, point *p3);
GLdouble dist(vector<GLdouble> * x0, vector<GLdouble> * x1);
pair< vector<GLdouble>, vector<GLdouble> > getNandB (vector<GLdouble>);
vector<GLdouble> unitCross(vector<GLdouble> a, vector<GLdouble> b);
vector<GLdouble> vecAddvec (vector<GLdouble> a, vector<GLdouble> b);
point ptAddvec (point a, vector<GLdouble> b);
vector<GLdouble> scalerMul (GLfloat s, vector<GLdouble> b);
vector<GLdouble> negVec (vector<GLdouble> a);


/* Write a screenshot to the specified filename */
void saveScreenshot (char *filename)
{
  int i, j;
  Pic *in = NULL;

  if (filename == NULL)
    return;

  /* Allocate a picture buffer */
  in = pic_alloc(640, 480, 3, NULL);

  printf("File to save to: %s\n", filename);

  for (i=479; i>=0; i--) {
    glReadPixels(0, 479-i, 640, 1, GL_RGB, GL_UNSIGNED_BYTE,
                 &in->pix[i*in->nx*in->bpp]);
  }

  if (jpeg_write(filename, in))
    printf("File saved Successfully\n");
  else
    printf("Error in Saving\n");

  pic_free(in);
}

/* 
  subdivide a spline segment recursively to
  1) Calculate location for each drawing point and store them.
  2) Calculate tangent vector of each drawing point and store them.
 */
void subdivide(GLdouble u0, GLdouble u1, GLdouble maxlinelength, point *p0, point *p1, point *p2, point *p3) {
  umid = (u0+u1)/2.0;
  x0 = getCatmullRom(u0, &basis_mat, p0, p1, p2, p3);
  x1 = getCatmullRom(u1, &basis_mat, p0, p1, p2, p3);

  if (dist(&x0, &x1) > maxlinelength) {
    subdivide(u0, umid, maxlinelength, p0, p1, p2, p3);
    subdivide(umid, u1, maxlinelength, p0, p1, p2, p3);
  } else {
    /* store x0, x1 as point structure to make semantic sense */
    pt_1st.x = x0.at(0);
    pt_1st.y = x0.at(1);
    pt_1st.z = x0.at(2);
    pt_2nd.x = x1.at(0);
    pt_2nd.y = x1.at(1);
    pt_2nd.z = x1.at(2);

    /* calculate x0 and x1's tangent vector */
    x0tangent = getUnitTangentVec(u0, &basis_mat, p0, p1, p2, p3);
    x1tangent = getUnitTangentVec(u1, &basis_mat, p0, p1, p2, p3);

    /* normalize */
    GLdouble x0TanVecLength = sqrt(x0tangent.at(0)*x0tangent.at(0) 
                                 + x0tangent.at(1)*x0tangent.at(1)
                                 + x0tangent.at(2)*x0tangent.at(2));
    GLdouble x1TanVecLength = sqrt(x1tangent.at(0)*x1tangent.at(0) 
                                 + x1tangent.at(1)*x1tangent.at(1)
                                 + x1tangent.at(2)*x1tangent.at(2));
    pt_1st_tangent.x = x0tangent.at(0)/x0TanVecLength;
    pt_1st_tangent.y = x0tangent.at(1)/x0TanVecLength;
    pt_1st_tangent.z = x0tangent.at(2)/x0TanVecLength;
    pt_2nd_tangent.x = x1tangent.at(0)/x1TanVecLength;
    pt_2nd_tangent.y = x1tangent.at(1)/x1TanVecLength;
    pt_2nd_tangent.z = x1tangent.at(2)/x1TanVecLength;

    /* store drawing points and their tangent vectors */
    drawingPoints.push_back(pt_1st);
    drawingPoints.push_back(pt_2nd);
    tangentVecs.push_back(pt_1st_tangent);
    tangentVecs.push_back(pt_2nd_tangent);
  }
}

/* calculate Euclidean distance between two points. Used in subdivide. */
GLdouble dist(vector<GLdouble> * x0, vector<GLdouble> * x1) {
  return sqrt (pow(x0->at(0) - x1->at(0), 2) 
             + pow(x0->at(1) - x1->at(1), 2) 
             + pow(x0->at(2) - x1->at(2), 2) );
}

/* get point locations using Catmull Rom spline */
vector<GLdouble> getCatmullRom (GLdouble u, GLfloat (*basis)[4][4], point *p0, point *p1, point *p2, point *p3) {
  vector<GLdouble> resultVec;

  /* Create [u^3 u^2 u 1] vector */
  GLdouble u2 = u*u;
  GLdouble u3 = u*u*u;
  GLdouble tmp[] = {u3, u2, u, 1};
  vector<GLdouble> u_vec(tmp, tmp+4);

  // [u^3 u^2 u 1] * basis_mat
  vector<GLdouble> leftMult;
  GLdouble sum;
  for (int m = 0; m < 4; m++) {
    sum = 0;
    for (int n = 0; n < 4; n++) {
      sum += u_vec.at(n) * (*basis)[n][m];
    }
    leftMult.push_back(sum);
  }

  // leftMult * control_points
  GLdouble tempx, tempy, tempz;
  tempx = leftMult.at(0) * p0->x 
        + leftMult.at(1) * p1->x
        + leftMult.at(2) * p2->x
        + leftMult.at(3) * p3->x;

  tempy = leftMult.at(0) * p0->y 
        + leftMult.at(1) * p1->y
        + leftMult.at(2) * p2->y
        + leftMult.at(3) * p3->y;

  tempz = leftMult.at(0) * p0->z 
        + leftMult.at(1) * p1->z
        + leftMult.at(2) * p2->z
        + leftMult.at(3) * p3->z;

  resultVec.push_back(tempx);
  resultVec.push_back(tempy);
  resultVec.push_back(tempz);

  return resultVec;
}

/* get normalized tangent vectors of points using Catmull Rom spline */
vector<GLdouble> getUnitTangentVec (GLdouble u, GLfloat (*basis)[4][4], point *p0, point *p1, point *p2, point *p3) {
  vector<GLdouble> resultVec;

  /* Create [3u^2 2u 1 0] vector */
  GLdouble tmp[] = {3*u*u, 2*u, 1, 0};
  vector<GLdouble> u_vec(tmp, tmp+4);

  // [3u^2 2u 1 0] * basis_mat
  vector<GLdouble> leftMult;
  GLdouble sum;
  for (int m = 0; m < 4; m++) {
    sum = 0;
    for (int n = 0; n < 4; n++) {
      sum += u_vec.at(n) * (*basis)[n][m];
    }
    leftMult.push_back(sum);
  }

  // leftMult * control_points
  GLdouble tempx, tempy, tempz;
  tempx = leftMult.at(0) * p0->x 
        + leftMult.at(1) * p1->x
        + leftMult.at(2) * p2->x
        + leftMult.at(3) * p3->x;

  tempy = leftMult.at(0) * p0->y 
        + leftMult.at(1) * p1->y
        + leftMult.at(2) * p2->y
        + leftMult.at(3) * p3->y;

  tempz = leftMult.at(0) * p0->z 
        + leftMult.at(1) * p1->z
        + leftMult.at(2) * p2->z
        + leftMult.at(3) * p3->z;

  /* normalize */
  GLdouble vecLength = sqrt(tempx*tempx + tempy*tempy + tempz*tempz);
  tempx = tempx / vecLength;
  tempy = tempy / vecLength;
  tempz = tempz / vecLength;

  resultVec.push_back(tempx);
  resultVec.push_back(tempy);
  resultVec.push_back(tempz);

  return resultVec;
}


void myinit() {
  glClearColor(0.0, 0.0, 0.0, 0.0);
  glEnable(GL_DEPTH_TEST);
  glShadeModel(GL_SMOOTH);

  // For all spline files
  for (int i = 0; i < g_iNumOfSplines; i++) {
    int currNumControlPts = g_Splines[i].numControlPoints;
    
    point Pstart = g_Splines[i].points[0];
    point Pend = g_Splines[i].points[currNumControlPts-1];

    /* Subdivide segments for current spline. */
    for (int j = 0; j < currNumControlPts - 3; j++) {
        subdivide(0, 1, maxlinelength, 
                  &g_Splines[i].points[j],
                  &g_Splines[i].points[j+1], 
                  &g_Splines[i].points[j+2],
                  &g_Splines[i].points[j+3]);
    }
  }

  groundDispList = glGenLists(1);
  glNewList(groundDispList, GL_COMPILE);
    glGenTextures(1, &groundTexID);
    g_groundTex = jpeg_read("ground.jpg", NULL);
    if (!g_groundTex) {
      printf ("error reading ground.jpg.\n");
      exit(1);
    }
    glBindTexture(GL_TEXTURE_2D, groundTexID);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, g_groundTex->nx, g_groundTex->ny, 0, GL_RGB, GL_UNSIGNED_BYTE, g_groundTex->pix);

    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
    glEnable(GL_TEXTURE_2D);
      glBegin(GL_QUADS);
        glTexCoord2f(0.0,0.0); glVertex3f(-20.0 , -10.0 , -5.0); 
        glTexCoord2f(0.0,1.0); glVertex3f(-20.0 , 50.0, -5.0); 
        glTexCoord2f(1.0,0.0); glVertex3f(40.0, 50.0, -5.0); 
        glTexCoord2f(1.0,1.0); glVertex3f(40.0,  -10.0, -5.0);
      glEnd();
    glDisable(GL_TEXTURE_2D);
  glEndList();


  skyDispList = glGenLists(1);
  glNewList(skyDispList, GL_COMPILE);
    glGenTextures(1, &skyTexID);
    g_skyTex = jpeg_read("sky.jpg", NULL);
    if (!g_skyTex) {
      printf ("error reading sky.jpg.\n");
      exit(1);
    }
    glBindTexture(GL_TEXTURE_2D, skyTexID);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, g_skyTex->nx, g_skyTex->ny, 0, GL_RGB, GL_UNSIGNED_BYTE, g_skyTex->pix);

    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
    glEnable(GL_TEXTURE_2D);
      glBegin(GL_QUADS);
        /* TOP */
        glTexCoord2f(0.0,0.0); glVertex3f(-20.0 , -10.0 , 35.0); 
        glTexCoord2f(0.0,1.0); glVertex3f(-20.0 , 50.0, 35.0); 
        glTexCoord2f(1.0,0.0); glVertex3f(40.0, 50.0, 35.0); 
        glTexCoord2f(1.0,1.0); glVertex3f(40.0,  -10.0, 35.0);

        /* LEFT */
        glTexCoord2f(0.0,0.0); glVertex3f(-20.0 , -10.0 , -5.0); 
        glTexCoord2f(0.0,1.0); glVertex3f(-20.0 , 50.0, -5.0); 
        glTexCoord2f(1.0,0.0); glVertex3f(-20.0, 50.0, 35.0); 
        glTexCoord2f(1.0,1.0); glVertex3f(-20.0,  -10.0, 35.0);

        /* RIGHT */
        glTexCoord2f(0.0,0.0); glVertex3f(40.0 , -10.0 , -5.0); 
        glTexCoord2f(0.0,1.0); glVertex3f(40.0 , 50.0, -5.0); 
        glTexCoord2f(1.0,0.0); glVertex3f(40.0, 50.0, 35.0); 
        glTexCoord2f(1.0,1.0); glVertex3f(40.0,  -10.0, 35.0);

        /* FRONT */
        glTexCoord2f(0.0,0.0); glVertex3f(-20.0 , 50.0 , -5.0); 
        glTexCoord2f(0.0,1.0); glVertex3f(40.0 , 50.0, -5.0); 
        glTexCoord2f(1.0,0.0); glVertex3f(40.0, 50.0, 35.0); 
        glTexCoord2f(1.0,1.0); glVertex3f(-20.0, 50.0, 35.0);


        /* BACK */
        glTexCoord2f(0.0,0.0); glVertex3f(-20.0 , -10.0 , -5.0); 
        glTexCoord2f(0.0,1.0); glVertex3f(40.0 , -10.0, -5.0); 
        glTexCoord2f(1.0,0.0); glVertex3f(40.0, -10.0, 35.0); 
        glTexCoord2f(1.0,1.0); glVertex3f(-20.0, -10.0, 35.0);
      glEnd();
    glDisable(GL_TEXTURE_2D);
  glEndList();

  /* calculate first local coordinate system */
  T.push_back(tangentVecs.at(0).x);
  T.push_back(tangentVecs.at(0).y);
  T.push_back(tangentVecs.at(0).z);
  V.push_back(1); V.push_back(1); V.push_back(1);

  N = unitCross(T,V);
  B = unitCross(T,N);
  Nvectors.push_back(N);
  Bvectors.push_back(B);
  /* calculate rest of local coordinate systems and store them */
  for (int i = 0; i < tangentVecs.size(); i++) {
    T.clear();  
    T.push_back(tangentVecs.at(i).x);
    T.push_back(tangentVecs.at(i).y);
    T.push_back(tangentVecs.at(i).z);

    N = unitCross(B, T);
    B = unitCross(T, N);

    Nvectors.push_back(N);
    Bvectors.push_back(B);
  }

  point v0u0, v1u0, v2u0, v3u0, 
        v0u1, v1u1, v2u1, v3u1;
  GLfloat alpha = 0.01;

  point v0,v1,v2,v3,                         // temp point for cross section positions
        leftV0, leftV1, leftV2, leftV3,      // temp point for left part of rail
        rightV0, rightV1, rightV2, rightV3;  // temp point for right part of rail
  vector<point> leftV0Data, leftV1Data, leftV2Data, leftV3Data,     // final data for left part of rail
                rightV0Data, rightV1Data, rightV2Data, rightV3Data; // final data for right part of rail
  vector<GLdouble> leftBotDirection,  // direction vector for left part of rail
                   rightBotDirection; // direction vector for right part of rail
  for (int i = 0; i <drawingPoints.size(); i++) {
    /* calculate four corners of a square shape cross section */
    v0 = ptAddvec(drawingPoints.at(i), scalerMul(alpha, vecAddvec(       Nvectors.at(i) ,        Bvectors.at(i))));
    v1 = ptAddvec(drawingPoints.at(i), scalerMul(alpha, vecAddvec(negVec(Nvectors.at(i)),        Bvectors.at(i))));
    v2 = ptAddvec(drawingPoints.at(i), scalerMul(alpha, vecAddvec(negVec(Nvectors.at(i)), negVec(Bvectors.at(i)))));
    v3 = ptAddvec(drawingPoints.at(i), scalerMul(alpha, vecAddvec(       Nvectors.at(i) , negVec(Bvectors.at(i)))));

    /* calculate the direction vector for moving v0,v1,v2,v3 to left bottom */
    leftBotDirection = scalerMul(0.1, vecAddvec(negVec(Nvectors.at(i)), negVec(Bvectors.at(i))));
    leftV0 = ptAddvec(v0, leftBotDirection);
    leftV1 = ptAddvec(v1, leftBotDirection);
    leftV2 = ptAddvec(v2, leftBotDirection);
    leftV3 = ptAddvec(v3, leftBotDirection);

    /* calculate the direction vector for moving v0,v1,v2,v3 to right bottom */
    rightBotDirection = scalerMul(0.1, vecAddvec(negVec(Nvectors.at(i)),        Bvectors.at(i)));
    rightV0 = ptAddvec(v0, rightBotDirection);
    rightV1 = ptAddvec(v1, rightBotDirection);
    rightV2 = ptAddvec(v2, rightBotDirection);
    rightV3 = ptAddvec(v3, rightBotDirection);

    /* store left cross section points and right cross section points */
    leftV0Data.push_back(leftV0);
    leftV1Data.push_back(leftV1);
    leftV2Data.push_back(leftV2);
    leftV3Data.push_back(leftV3);

    rightV0Data.push_back(rightV0);
    rightV1Data.push_back(rightV1);
    rightV2Data.push_back(rightV2);
    rightV3Data.push_back(rightV3);
  }

  xSectionDispList = glGenLists(1);
  glNewList(xSectionDispList, GL_COMPILE);
      for (int i=0; i<drawingPoints.size(); i+=2) {
        /* for left rail */
        v0u0 = leftV0Data.at(i);
        v1u0 = leftV1Data.at(i);
        v2u0 = leftV2Data.at(i);
        v3u0 = leftV3Data.at(i);
        v0u1 = leftV0Data.at(i+1);
        v1u1 = leftV1Data.at(i+1);
        v2u1 = leftV2Data.at(i+1);
        v3u1 = leftV3Data.at(i+1);

        /* connect 8 points with triangles (left rail) */
        glBegin(GL_TRIANGLE_STRIP);
          glVertex3d(v3u0.x, v3u0.y, v3u0.z); // Front-top-left
          glVertex3d(v0u0.x, v0u0.y, v0u0.z); // Front-top-right
          glVertex3d(v2u0.x, v2u0.y, v2u0.z); // Front-bottom-left
          glVertex3d(v1u0.x, v1u0.y, v1u0.z); // Front-bottom-right

          glVertex3d(v1u1.x, v1u1.y, v1u1.z); // Back-bottom-right
          glVertex3d(v0u0.x, v0u0.y, v0u0.z); // Front-top-right
          glVertex3d(v0u1.x, v0u1.y, v0u1.z); // Back-top-right
          glVertex3d(v3u0.x, v3u0.y, v3u0.z); // Front-top-left
          glVertex3d(v3u1.x, v3u1.y, v3u1.z); // Back-top-left
          glVertex3d(v2u0.x, v2u0.y, v2u0.z); // Front-bottom-left

          glVertex3d(v2u1.x, v2u1.y, v2u1.z); // Back-bottom-left
          glVertex3d(v1u1.x, v1u1.y, v1u1.z); // Back-bottom-right
          glVertex3d(v3u1.x, v3u1.y, v3u1.z); // Back-top-left
          glVertex3d(v0u1.x, v0u1.y, v0u1.z); // Back-top-right
        glEnd();

        /* for right rail */
        v0u0 = rightV0Data.at(i);
        v1u0 = rightV1Data.at(i);
        v2u0 = rightV2Data.at(i);
        v3u0 = rightV3Data.at(i);

        v0u1 = rightV0Data.at(i+1);
        v1u1 = rightV1Data.at(i+1);
        v2u1 = rightV2Data.at(i+1);
        v3u1 = rightV3Data.at(i+1);

        /* connect 8 points with triangles (right rail) */
        glBegin(GL_TRIANGLE_STRIP);
          glVertex3d(v3u0.x, v3u0.y, v3u0.z); // Front-top-left
          glVertex3d(v0u0.x, v0u0.y, v0u0.z); // Front-top-right
          glVertex3d(v2u0.x, v2u0.y, v2u0.z); // Front-bottom-left
          glVertex3d(v1u0.x, v1u0.y, v1u0.z); // Front-bottom-right

          glVertex3d(v1u1.x, v1u1.y, v1u1.z); // Back-bottom-right
          glVertex3d(v0u0.x, v0u0.y, v0u0.z); // Front-top-right
          glVertex3d(v0u1.x, v0u1.y, v0u1.z); // Back-top-right
          glVertex3d(v3u0.x, v3u0.y, v3u0.z); // Front-top-left
          glVertex3d(v3u1.x, v3u1.y, v3u1.z); // Back-top-left
          glVertex3d(v2u0.x, v2u0.y, v2u0.z); // Front-bottom-left

          glVertex3d(v2u1.x, v2u1.y, v2u1.z); // Back-bottom-left
          glVertex3d(v1u1.x, v1u1.y, v1u1.z); // Back-bottom-right
          glVertex3d(v3u1.x, v3u1.y, v3u1.z); // Back-top-left
          glVertex3d(v0u1.x, v0u1.y, v0u1.z); // Back-top-right
        glEnd();

        /* draw crossbars in every 2000 points inteval */
        if (i % 2000 == 0) {
          leftV1 = leftV1Data.at(i);
          rightV1 = rightV1Data.at(i);
          glLineWidth(80);
          glBegin(GL_LINES);
            glColor4f(0.894,0.741,0.043,1);
            glVertex3d(leftV1.x, leftV1.y, leftV1.z);
            glVertex3d(rightV1.x, rightV1.y, rightV1.z);
          glEnd();
        }
      }
  glEndList();
}

/* adds two vectors, returns vector */
vector<GLdouble> vecAddvec (vector<GLdouble> a, vector<GLdouble> b) {
  vector<GLdouble> c;
  c.push_back(a.at(0) + b.at(0));
  c.push_back(a.at(1) + b.at(1));
  c.push_back(a.at(2) + b.at(2));
  return c;
}

/* adds a point and a vector, returns point */
point ptAddvec (point a, vector<GLdouble> b) {
  point c;
  c.x = a.x + b.at(0);
  c.y = a.y + b.at(1);
  c.z = a.z + b.at(2);
  return c;
}

/* multiply input vector b by scaler s */
vector<GLdouble> scalerMul (GLfloat s, vector<GLdouble> b) {
  vector<GLdouble> c;
  c.push_back(s * b.at(0));
  c.push_back(s * b.at(1));
  c.push_back(s * b.at(2));
  return c;
}

/* neagate a vector */
vector<GLdouble> negVec (vector<GLdouble> a) {
  vector<GLdouble> b;
  b.push_back( - a.at(0));
  b.push_back( - a.at(1));
  b.push_back( - a.at(2));
  return b;
}

/* calculate cross product of two vectors, returned vector is normalized */
vector<GLdouble> unitCross(vector<GLdouble> a, vector<GLdouble> b) {
  vector<GLdouble> retVec;
  GLdouble x,y,z;
  x = a.at(1)*b.at(2) - a.at(2)*b.at(1);
  y = a.at(2)*b.at(0) - a.at(0)*b.at(2);
  z = a.at(0)*b.at(1) - a.at(1)*b.at(0);

  GLdouble vecLength = sqrt(x*x+y*y+z*z);

  retVec.push_back(x/vecLength);
  retVec.push_back(y/vecLength);
  retVec.push_back(z/vecLength);

  return retVec;
}


int fCnt = 0; /* counts frame for moving camera */
void display() {
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glLoadIdentity();

  /* setup modelview matrix by loading the current model-view matrix
     and do transformation on this current matrix */
  glMatrixMode(GL_MODELVIEW);

  /* move the camera along the spline */
  if (fCnt == 0) {
    gluLookAt(drawingPoints.at(0).x, 
              drawingPoints.at(0).y, 
              drawingPoints.at(0).z,
              drawingPoints.at(0).x + tangentVecs.at(0).x, 
              drawingPoints.at(0).y + tangentVecs.at(0).y, 
              drawingPoints.at(0).z + tangentVecs.at(0).z,
              Nvectors.at(0).at(0), 
              Nvectors.at(0).at(1), 
              Nvectors.at(0).at(2));
  } 
  else {
    /* move 800 points in each frame */
    for (int i = fCnt*4000; i < fCnt*4000+4000; i++) {
      /* if reached last point of spline, look at last pt's T,N */
      if (i >= tangentVecs.size()) {
        int lastPtInd = tangentVecs.size() - 1;
        gluLookAt(drawingPoints.at(lastPtInd).x, 
                  drawingPoints.at(lastPtInd).y, 
                  drawingPoints.at(lastPtInd).z,
                  drawingPoints.at(lastPtInd).x + tangentVecs.at(lastPtInd).x, 
                  drawingPoints.at(lastPtInd).y + tangentVecs.at(lastPtInd).y, 
                  drawingPoints.at(lastPtInd).z + tangentVecs.at(lastPtInd).z,
                  Nvectors.at(lastPtInd).at(0), 
                  Nvectors.at(lastPtInd).at(1), 
                  Nvectors.at(lastPtInd).at(2));
        break;
      }
      /* for rest of points, move the camera to the corresponding point */
      glLoadIdentity();
      gluLookAt(drawingPoints.at(i).x, 
                drawingPoints.at(i).y, 
                drawingPoints.at(i).z,
                drawingPoints.at(i).x + tangentVecs.at(i).x, 
                drawingPoints.at(i).y + tangentVecs.at(i).y, 
                drawingPoints.at(i).z + tangentVecs.at(i).z,
                Nvectors.at(i).at(0), 
                Nvectors.at(i).at(1), 
                Nvectors.at(i).at(2));
    }
  }

  fCnt++;

  glCallList(xSectionDispList);
  glCallList(groundDispList);
  glCallList(skyDispList);

  glutSwapBuffers();
}


void reshape(int w, int h) {
  /* maintain image aspect ratio */
  glViewport(0, 0, w, h);

  /* setup projection */
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  // glFrustum(-0.1, 0.1, 
  //     -float(h)/(10.0*float(w)), 
  //     float(h)/(10.0*float(w)), 0.03, 1000.0);

  gluPerspective(45.0, float(w)/float(h),0.1,1000.0);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}


void doIdle() { //TODO: last animation
  if (animCounter < animDelay) {
    animCounter++;
  }
  
  /* uncomment code inside if to take screenshots */
  char myFilenm[2048];
  if (animCounter >= animDelay && animCounter < animDelay+500)
  {
    // sprintf(myFilenm, "anim.%04d.jpg", animCounter-animDelay);
    // saveScreenshot(myFilenm);
    // animCounter++;
  }

  glutPostRedisplay();
}

void menufunc(int value)
{
  switch (value)
  {
    case 0:
      exit(0);
      break;
  }
}


void mousebutton(int button, int state, int x, int y)
{
  switch (button)
  {
    case GLUT_LEFT_BUTTON:
      g_iLeftMouseButton = (state==GLUT_DOWN);
      break;
    case GLUT_MIDDLE_BUTTON:
      g_iMiddleMouseButton = (state==GLUT_DOWN);
      break;
    case GLUT_RIGHT_BUTTON:
      g_iRightMouseButton = (state==GLUT_DOWN);
      break;
  }
}

int loadSplines(char *argv) {
  char *cName = (char *)malloc(128 * sizeof(char));
  FILE *fileList;
  FILE *fileSpline;
  int iType, i = 0, j, iLength;

  // printf("arg %s: ",argv);
  //RESULT IS track.txt

  /* load the track file */
  fileList = fopen(argv, "r");
  if (fileList == NULL) {
    printf ("can't open file\n");
    exit(1);
  }
  
  /* stores the number of splines in a global variable */
  fscanf(fileList, "%d", &g_iNumOfSplines);

  g_Splines = (struct spline *)malloc(g_iNumOfSplines * sizeof(struct spline));

  /* reads through the spline files */
  for (j = 0; j < g_iNumOfSplines; j++) {
    i = 0;
    fscanf(fileList, "%s", cName);
    fileSpline = fopen(cName, "r");

    if (fileSpline == NULL) {
      printf ("can't open file\n");
      exit(1);
    }

    /* gets length for spline file */
    fscanf(fileSpline, "%d %d", &iLength, &iType);

    /* allocate memory for all the points */
    g_Splines[j].points = (struct point *)malloc(iLength * sizeof(struct point));
    g_Splines[j].numControlPoints = iLength;

    /* saves the data to the struct */
    while (fscanf(fileSpline, "%lf %lf %lf", 
	   &g_Splines[j].points[i].x, 
	   &g_Splines[j].points[i].y, 
	   &g_Splines[j].points[i].z) != EOF) {
      i++;
    }
  }

  free(cName);

  return 0;
}




int main (int argc, char ** argv) {
  
  if (argc<2) {  
    printf ("usage: %s <trackfile>\n", argv[0]);
    exit(0);
  }

  loadSplines(argv[1]);

  /* OpenGL setup */
  glutInit(&argc,argv);

  /* set up double buffering and depth test*/
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_DEPTH);

  /* create window*/
  glutInitWindowSize(640, 480);
  glutInitWindowPosition(0, 0);
  glutCreateWindow("Roller Coaster");

  /* tells glut to use a particular display function to redraw */
  glutDisplayFunc(display);

  /* allow the user to quit using the right mouse button menu */
  g_iMenuId = glutCreateMenu(menufunc);
  glutSetMenu(g_iMenuId);
  glutAddMenuEntry("Quit",0);
  glutAttachMenu(GLUT_MIDDLE_BUTTON);

  glutReshapeFunc(reshape);

  glutIdleFunc(doIdle);
  glutMouseFunc(mousebutton);

  myinit();

  
  glutMainLoop();

  return(0);
}
