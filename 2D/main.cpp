#include "Fluid.h"

int window_width = 800;
int window_height = 800;

Fluid fluid;

void display(void) {
    // everything we need in order to redraw the scene
}

void reshape(int w, int h) {
    window_width = w;
    window_height = h;
    
    // continue function here
}

void mouse(int button, int state, int x, int y) {
    switch (button) {
        case GLUT_LEFT_BUTTON:
            if (state == GLUT_DOWN) {
                // do something here, like create a new object
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

void idle(void) {
    // render the next frame of our simulation
    fluid.step();
    glutPostRedisplay();
}

int main(int argc, char* argv[]) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(window_width, window_height);
    glutInitWindowPosition(100, 100);
    glutCreateWindow("Stable Fluids");
    
    fluid.init();
    
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutMouseFunc(mouse);
    glutMotionFunc(motion);
    glutKeyboardFunc(keyboard);
    glutIdleFunc(idle);
    glutMainLoop();
    
    return 0;
}
