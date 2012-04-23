//============================================================================
// Name        : Media_y_PD.cpp
// Author      : Angelo M. Calvão
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <math.h>
#include <fstream>
#include <iostream>
#include <GL/glut.h>    // Header File For The GLUT Library
#include <GL/gl.h>	// Header File For The OpenGL32 Library
#include <GL/glu.h>	// Header File For The GLu32 Library
#include <GL/glx.h>     // Header file fot the glx libraries.
#include "Models.h"

int window, rrate, opt = 0, i;
GLuint base; // base display list for the font set.
int cont = 0;

using namespace std;

double_pendulum* Pendulo;


/*

GLvoid BuildFont(GLvoid)
{
   Display *dpy;
   XFontStruct *fontInfo;  // storage for our font.

   base = glGenLists(96);                      // storage for 96 characters.

   // load the font.  what fonts any of you have is going
   // to be system dependent, but on my system they are
   // in /usr/X11R6/lib/X11/fonts/*, with fonts.alias and
   // fonts.dir explaining what fonts the .pcf.gz files
   // are.  in any case, one of these 2 fonts should be
   // on your system...or you won't see any text.

   // get the current display.  This opens a second
   // connection to the display in the DISPLAY environment
   // value, and will be around only long enough to load
   // the font.
   dpy = XOpenDisplay(NULL); // default to DISPLAY env.

   fontInfo = XLoadQueryFont(dpy, "-adobe-helvetica-medium-r-normal--10-*-*-*-p-*-iso8859-1");
   if (fontInfo == NULL) {
       fontInfo = XLoadQueryFont(dpy, "fixed");
       if (fontInfo == NULL) {
           printf("no X font available?\n");
       }
   }

   // after loading this font info, this would probably be the time
   // to rotate, scale, or otherwise twink your fonts.

   // start at character 32 (space), get 96 characters (a few characters past z), and
   // store them starting at base.
   glXUseXFont(fontInfo->fid, 32, 96, base);

   // free that font's info now that we've got the
   // display lists.
   XFreeFont(dpy, fontInfo);

   // close down the 2nd display connection.
   XCloseDisplay(dpy);
}

GLvoid KillFont(GLvoid)                         // delete the font.
{
   glDeleteLists(base, 96);                    // delete all 96 characters.
}

GLvoid glPrint(char *text)                      // custom gl print routine.
{
   if (text == NULL) {                         // if there's no text, do nothing.
       return;
   }

   glPushAttrib(GL_LIST_BIT);                  // alert that we're about to offset the display lists with glListBase
   glListBase(base - 32);                      // sets the base character to 32.

   glCallLists(strlen(text), GL_UNSIGNED_BYTE, text); // draws the display list text.
   glPopAttrib();                              // undoes the glPushAttrib(GL_LIST_BIT);
}

 */

/* A general OpenGL initialization function.  Sets all of the initial parameters. */
void InitGL(int Width, int Height) // We call this right after our OpenGL window is created.
{
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f); // This Will Clear The Background Color To Black
    glClearDepth(1.0); // Enables Clearing Of The Depth Buffer
    glDepthFunc(GL_LESS); // The Type Of Depth Test To Do
    glEnable(GL_DEPTH_TEST); // Enables Depth Testing
    glShadeModel(GL_SMOOTH); // Enables Smooth Color Shading

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity(); // Reset The Projection Matrix

    gluPerspective(45.0f, (GLfloat) Width / (GLfloat) Height, 0.1f, 100.0f); // Calculate The Aspect Ratio Of The Window

    glMatrixMode(GL_MODELVIEW);

    //BuildFont();
}

/* The function called when our window is resized (which shouldn't happen, because we're fullscreen) */
void ReSizeGLScene(int Width, int Height) {
    if (Height == 0) // Prevent A Divide By Zero If The Window Is Too Small
        Height = 1;

    glViewport(0, 0, Width, Height); // Reset The Current Viewport And Perspective Transformation

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    gluPerspective(45.0f, (GLfloat) Width / (GLfloat) Height, 0.1f, 100.0f);
    glMatrixMode(GL_MODELVIEW);
}

/* The main drawing function. */



void DrawGLScene() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    glTranslatef(0.0f, 0.0f, -30.0f);
    //glRotatef(ax, 0, 1, 0);
    //glRotatef(ay, 1, 0, 0);
    glPointSize(1.5);
    float x1, y1, x, y, r = 0.2;

    GLUquadric* esfera;
    esfera = gluNewQuadric();



    /*
    double t=argument[0];
    double omega=argument[1];
    double theta=argument[2];
    double vel=argument[3];
    double pos=argument[4];
    double F=argument[5];
     */

    x = 2*sin((*Pendulo)[0]);
    y = -2*cos((*Pendulo)[0]);
    x1 = x+2*sin((*Pendulo)[1]);
    y1 = y-2*cos((*Pendulo)[1]);
    //Pendulo.set_F(+1*Pendulo.get_theta()-5*Pendulo.get_omega());

 
    (*Pendulo).next();
    
    //**************************************************************************************
    //Codigo que gera o pendulo
    glColor3f(0.3f, 0.8f, 0.3f);
    //WireSphere(r,25,25,x,y,0);


    glBegin(GL_LINES);
    glColor3f(1.0f, 0.0f, 0.0f);
    glVertex3f(0.0f, 0, 0.0f);
    glColor3f(0.0f, 1.0f, 0.0f);
    glVertex3f(x, y, 0.0f);
    glEnd();
    glColor3f(0.3f, 0.8f, 0.3f);
    //WireSphere(r,25,25,x1,y1,0);
    glBegin(GL_LINES);
    glColor3f(0.0f, 1.0f, 0.0f);
    glVertex3f(x, y, 0.0f);
    glColor3f(0.0f, 1.0f, 0.0f);
    glVertex3f(x1, y1, 0.0f);
    glEnd();
    //fim

    //Codigo que gera os eixos
    glBegin(GL_LINES);
    glColor3f(0.0f, 0.0f, 1.0f);
    glVertex3f(10.0f, 0.0, 0.0f);
    glVertex3f(-10.0f, 0.0, 0.0f);
    glEnd();

    glPushMatrix();
    //glTranslatef(2 * x, 2 * y, -30.0f);
    //gluSphere(esfera, 2 * r, 10, 10);
    //glTranslatef(-2*x,-2*y,-30);
    //glTranslated(0,0,-30);
    glPopMatrix();

    //glTranslatef(2 * x1, 2 * y1, -30.0f);
    //gluSphere(esfera, 2 * r, 10, 10);
    glPopMatrix();
    //glLoadIdentity();
    glutSwapBuffers();
}

int Pendulum_GL(int argc , char * argv[],double angle1,double angle2) {

	vector<double> variable(4),parameter(5);
        variable[V_THETA1] = angle1;
        variable[V_THETA2] = angle2;
        variable[V_OMEGA1] = 0.0;
        variable[V_OMEGA2] = 0.0;
        parameter[P_L1]= 0.30;
        parameter[P_L2]= 0.10;
        parameter[P_M1]= 0.25;
        parameter[P_M2]= 0.10;
        parameter[P_G]= 9.8;

        double_pendulum dp_attractor(variable,parameter,0.001,"rk");
 
	Pendulo = &dp_attractor;
    //srand(time(0));
    glutInit(&argc, argv);

    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
    glutInitWindowSize(800, 600);
    glutInitWindowPosition(60, 60);

    window = glutCreateWindow("Pendulo");

    glutDisplayFunc(&DrawGLScene);
    glutIdleFunc(&DrawGLScene);
    glutReshapeFunc(&ReSizeGLScene);


    InitGL(800, 600);
    glutMainLoop();
    delete[] Pendulo;
    return 1;
}

