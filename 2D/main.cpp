#include "Fluid.h"

int window_width = 800;
int window_height = 800;

Fluid fluid;
float F[2];
float Ssource;

void init(void) {
    glClearColor(0.0, 0.0, 0.0, 0.0);
    F[0] = -2.0f; F[1] = 3.0f;
    Ssource = 200.0f;
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
            
            // (TODO) these 1500s should not be here; they're just here for debugging help
            c00 = fluid.S_at(r, c)  / 1500.0f;
            c01 = fluid.S_at(r, c + 1) / 1500.0f;
            c10 = fluid.S_at(r + 1, c) / 1500.0f;
            c11 = fluid.S_at(r + 1, c + 1) / 1500.0f;
            
            glColor4f(c11, c11, 0.1, alpha); glVertex2f((x + vsize) / 25.0f - 0.5f, (y + vsize) / 25.0f - 0.5f);
            glColor4f(c01, c01, 0.2, alpha); glVertex2f(x / 25.0f - 0.5f, (y + vsize) / 25.0f - 0.5f);
            glColor4f(c00, c00, 0.3, alpha); glVertex2f(x / 25.0f - 0.5f, y / 25.0f - 0.5f);
            glColor4f(c10, c10, 0.4, alpha); glVertex2f((x + vsize) / 25.0f - 0.5f, y / 25.0f - 0.5f);
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
                F[0] += 12.0f;
                F[1] += 15.0f;
                Ssource += 120.0f;
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
    if (F[0] != 0) { F[0] = 0; }
    if (F[1] != 0) { F[1] = 0; }
    if (Ssource != 0) { Ssource = 0; }
    glutPostRedisplay();
}

int main(int argc, char* argv[]) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(window_width, window_height);
    glutInitWindowPosition(100, 100);
    glutCreateWindow("Stable Fluids");
    
    init();
    fluid.init(0.1f, 0.2f, 0.3f, 0.4f);
    
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutMouseFunc(mouse);
    glutMotionFunc(motion);
    glutKeyboardFunc(keyboard);
    glutIdleFunc(idle);
    glutMainLoop();
    
    return 0;
}
