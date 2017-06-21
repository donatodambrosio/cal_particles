#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include <model.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <utils_io.h>

struct ViewPort{
  int w; // width
  int h; // height
};
struct ViewPort wp;

struct ModelView{
  GLfloat x_rot;
  GLfloat y_rot;
  GLfloat z_trans;
};

struct ModelView model_view;
time_t start_time, end_time;

void drawAxes()
{
  glBegin(GL_LINES);

  glColor3f(1.0, 0.0, 0.0);
  glVertex3f(0.0, 0.0, 0.0);
  glVertex3f(1.0, 0.0, 0.0);

  glColor3f(0.0, 1.0, 0.0);
  glVertex3f(0.0, 0.0, 0.0);
  glVertex3f(0.0, 1.0, 0.0);

  glColor3f(0.0, 0.0, 1.0);
  glVertex3f(0.0, 0.0, 0.0);
  glVertex3f(0.0, 0.0, 1.0);

  glEnd();
}

void drawParticles()
{
  //float particle_size = 1.0/MAX_NUMBER_OF_PARTICLES_PER_CELL;
  float particle_radius = PARTICLE_RADIUS;

  CALint cell_x, cell_y, cell_z;
  CALreal px, py, pz;

  // Box
  glColor3f(1,1,1);
  glPushMatrix();
  float scale_x = X_CELLS-2;
  float scale_y = Y_CELLS-2;
  float scale_z = Z_CELLS-2;

  particle_radius *= scale_x*scale_y*scale_z;

  glScalef(scale_x,scale_y,scale_z);
  glPushAttrib(GL_LIGHTING_BIT);
  glDisable(GL_LIGHTING);
  glutWireCube(1.0);
  glPopAttrib();
  glPopMatrix();

  // Particles
  glColor3f(1,0,0);
  for (cell_x=0; cell_x<X_CELLS; cell_x++)
    for (cell_y=0; cell_y<Y_CELLS; cell_y++)
      for (cell_z=0; cell_z<Z_CELLS; cell_z++)
        if (calGet3Di(u_modellu,Q.imove[0],cell_x,cell_y,cell_z) != PARTICLE_BORDER)
          {
            for(int slot=0;slot<MAX_NUMBER_OF_PARTICLES_PER_CELL;slot++)
              if(calGet3Di(u_modellu,Q.imove[slot],cell_x,cell_y,cell_z) != PARTICLE_ABSENT)
                {
                  px = calGet3Dr(u_modellu,Q.rx[slot],cell_x,cell_y,cell_z) / CELL_SIDE;
                  py = calGet3Dr(u_modellu,Q.ry[slot],cell_x,cell_y,cell_z) / CELL_SIDE;
                  pz = calGet3Dr(u_modellu,Q.rz[slot],cell_x,cell_y,cell_z) / CELL_SIDE;

                  glPushMatrix();
                  glTranslatef(-X_CELLS/2, -Y_CELLS/2 , -Z_CELLS/2);
                  glTranslated(px,py,pz);
                  glutSolidSphere(particle_radius,10,10);
                  //glutWireSphere(particle_radius,20,20);
                  glPopMatrix();
                }
          }
}

void display(void)
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // particles
  GLfloat	 lightPos[]	= { 0.0f, 0.0f, 100.0f, 1.0f };
  int MAX = u_modellu->rows;

  if (MAX < u_modellu->columns)
    MAX = u_modellu->columns;
  if (MAX < u_modellu->slices)
    MAX = u_modellu->slices;

  glViewport (0, 0, (GLsizei)wp.w, (GLsizei)wp.h);
  glMatrixMode (GL_PROJECTION);
  glLoadIdentity ();
  gluPerspective(45.0, (GLfloat) wp.w/(GLfloat) wp.h, 1.0, 4*MAX);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt (0.0, 0.0, 2*MAX, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
  glPushMatrix();
    glTranslatef(0, 0, model_view.z_trans);
    glRotatef(model_view.x_rot, 1, 0, 0);
    glRotatef(model_view.y_rot, 0, 1, 0);
    glRotatef(-90,1,0,0);
    drawParticles();
  glPopMatrix();

  // axes
  glViewport (0, 0,200, 200);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(45.0, 1.0, 1.0, 10);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt (0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);

  lightPos[2] = 2*MAX;
  glLightfv(GL_LIGHT0, GL_POSITION, lightPos);

  glPushMatrix();
  glPushAttrib(GL_LIGHTING_BIT);
    glDisable(GL_LIGHTING);
    glRotatef(model_view.x_rot, 1, 0, 0);
    glRotatef(model_view.y_rot, 0, 1, 0);
    glRotatef(-90,1,0,0);
    drawAxes();
  glPopAttrib();
  glPopMatrix();

  glutSwapBuffers();
}

void simulationRun(void)
{
  CALbyte again;

  //exectutes the global transition function, the steering function and check for the stop condition.
  again = calRunCAStep3D(a_simulazioni);

  //simulation main loop
  a_simulazioni->step++;

#ifdef VERBOSE
  //graphic rendering
  printf("step: %d; \tactive cells: %d\r", a_simulazioni->step, a_simulazioni->ca3D->A.size_current);
  glutPostRedisplay();
#endif

  //check for the stop condition
  if (!again)
    {
      //breaking the simulation
      end_time = time(NULL);
      glutIdleFunc(NULL);
      printf("\n");
      printf("Simulation terminated\n");
      printf("Elapsed time: %lds\n", end_time - start_time);

      //graphic rendering
      printf("step: %d; \tactive cells: %d\r", a_simulazioni->step, a_simulazioni->ca3D->A.size_current);
      glutPostRedisplay();
      return;
    }

  char winwow_title[256];
  char steps_cstr[64];
  strcpy(winwow_title, "cal_DEM - step ");
  sprintf(steps_cstr, "%d, elapsed_time %.3f s", a_simulazioni->step, elapsed_time);
  strcat(winwow_title, steps_cstr);
  glutSetWindowTitle(winwow_title);
}

void init(void)
{
  GLfloat  ambientLight[] = { 0.2f, 0.2f, 0.2f, 1.0f };
  GLfloat  diffuseLight[] = { 0.75f, 0.75f, 0.75f, 1.0f };

  glEnable(GL_LIGHTING);
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambientLight);
  glLightfv(GL_LIGHT0, GL_DIFFUSE,  diffuseLight);
  glEnable(GL_LIGHT0);
  glEnable(GL_COLOR_MATERIAL);
  glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);

  glClearColor (0.0, 0.0, 0.0, 0.0);
  glShadeModel (GL_FLAT);

  glEnable (GL_DEPTH_TEST);

  model_view.x_rot = 0.0;
  model_view.y_rot = 0.0;
  model_view.z_trans = 0.0;

  printf("The 3D particles computational model\n");
  printf("Left click on the graphic window to start the simulation\n");
  printf("Right click on the graphic window to stop the simulation\n");
}

void reshape(int w, int h)
{
  wp.w = w;
  wp.h = h;
}

void mouse(int button, int state, int x, int y)
{
  switch (button) {
    case GLUT_LEFT_BUTTON:
      if (state == GLUT_DOWN)
        {
          start_time = time(NULL);
          glutIdleFunc(simulationRun);
        }
      break;
    case GLUT_MIDDLE_BUTTON:
    case GLUT_RIGHT_BUTTON:
      if (state == GLUT_DOWN)
        glutIdleFunc(NULL);
      break;
    default:
      break;
    }
}

void keyboard(unsigned char key, int x, int y)
{
  switch(key)
    {
    case ' ':
      mouse(GLUT_LEFT_BUTTON, GLUT_DOWN, x, y);
      break;
    case 'p':
      mouse(GLUT_RIGHT_BUTTON, GLUT_DOWN, x, y);
      break;
    case 27: //ESC
      exit(EXIT_SUCCESS);
    }

}

void specialKeys(int key, int x, int y){

  //GLubyte specialKey = glutGetModifiers();
  const GLfloat x_rot = 5.0, y_rot = 5.0, z_trans = 5.0;

  if (key==GLUT_KEY_F9)
    mouse(GLUT_LEFT_BUTTON, GLUT_DOWN, x, y);

  if(key==GLUT_KEY_DOWN){
      model_view.x_rot += x_rot;
    }
  if(key==GLUT_KEY_UP){
      model_view.x_rot -= x_rot;
    }
  if(key==GLUT_KEY_LEFT){
      model_view.y_rot -= y_rot;
    }
  if(key==GLUT_KEY_RIGHT){
      model_view.y_rot += y_rot;
    }
  if(key == GLUT_KEY_PAGE_UP){
      model_view.z_trans += z_trans;
    }
  if(key == GLUT_KEY_PAGE_DOWN){
      model_view.z_trans -= z_trans;
    }

  glutPostRedisplay();
}

int main(int argc, char** argv)
{
  partilu();

  if (argc > 1 && !strcmp(argv[1], "gl") )
    {
      glutInit(&argc, argv);
      glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
      glutInitWindowSize(1200, 800);
      glutInitWindowPosition(100, 100);
      glutCreateWindow("cal_DEM");
      init();
      glutDisplayFunc(display);
      glutReshapeFunc(reshape);
      glutSpecialFunc(specialKeys);
      glutKeyboardFunc(keyboard);
      glutMouseFunc(mouse);
      glutMainLoop();
    }
  else
    {

      clock_t begin = clock();
      clock_t end = begin;
      double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
      saveParticles(u_modellu, a_simulazioni->step, elapsed_time, time_spent, "./data/particles_t0.txt");

      calRun3D(a_simulazioni);

      end =  clock();
      time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
      saveParticles(u_modellu, a_simulazioni->step, elapsed_time, time_spent, "./data/particles_tf.txt");
    }

  return 0;
}
