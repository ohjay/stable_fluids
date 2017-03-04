#include "Fluid.h"

int window_width = 800;
int window_height = 800;

Fluid fluid;
float F[2];
float Ssource;

void init(void) {
    glClearColor(0.0, 0.0, 0.0, 0.0);
    F[0] = 5.0f; F[1] = 5.0f;
    Ssource = 150.0f;
}

// everything we need in order to redraw the scene
void display(void) {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glPushMatrix();
    
    // draw the density grid
    float alpha = 0.03;
    float x, y, c00, c01, c10, c11;
    
    int h = fluid.num_cells_y();
    int w = fluid.num_cells_x();
    float vsize = fluid.grid_spacing();
    
    glBegin(GL_QUADS);
    
    for (int r = 0; r < h; ++r) {
        y = (r - 0.5f) * vsize;
        for (int c = 0; c < w; ++c) {
            x = (c - 0.5f) * vsize;
            
            c00 = fluid.S_at(r, c);
            c01 = fluid.S_at(r, c + 1);
            c10 = fluid.S_at(r + 1, c);
            c11 = fluid.S_at(r + 1, c + 1);
            
            glColor4f(c11, c11, 0.1, alpha); glVertex2f(x + vsize, y + vsize);
            glColor4f(c01, c01, 0.2, alpha); glVertex2f(x, y + vsize);
            glColor4f(c00, c00, 0.3, alpha); glVertex2f(x, y);
            glColor4f(c10, c10, 0.4, alpha); glVertex2f(x + vsize, y);
        }
    }

	glEnd ();
	glFlush();
	glutSwapBuffers();
}

void reshape(int w, int h) {
    window_width = w;
    window_height = h;
    glViewport(0, 0, (GLsizei) w, (GLsizei) h);
}

void mouse(int button, int state, int x, int y) {
    switch (button) {
        case GLUT_LEFT_BUTTON:
            if (state == GLUT_DOWN) {
                // do something here, like create a new object
                F[0] += 5.0f;
                F[1] -= 1.0f;
            }
            break;
        default:
            break;
    }
}

void motion(int x, int y) {
    // do something here, like apply a force
}

void keyboard(unsigned char key, int x, int y) {
    switch (key) {
        case 'q':
            // do something here, like quit
            break;
        default:
            break;
    }
}

// render the next frame of our simulation
void idle(void) {
    fluid.step(F, Ssource);
    glutPostRedisplay();
}

int main(int argc, char* argv[]) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(window_width, window_height);
    glutInitWindowPosition(100, 100);
    glutCreateWindow("Stable Fluids");
    
    init();
    fluid.init(0.2f, 0.3f, 0.4f, 0.5f);
    
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutMouseFunc(mouse);
    glutMotionFunc(motion);
    glutKeyboardFunc(keyboard);
    glutIdleFunc(idle);
    glutMainLoop();
    
    return 0;
}
