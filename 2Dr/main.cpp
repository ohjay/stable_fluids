#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "Fluid.h"

using namespace std;

int window_height = WINDOW_HEIGHT;
int window_width = WINDOW_WIDTH;
bool left_mouse_down, paused;
int mouse_y, mouse_x;

Fluid fluid;
float source, add_amount, force_y, force_x;
double cr, cg, cb, alpha;

void init(void) {
    glClearColor(0.0, 0.0, 0.0, 0.0);

    left_mouse_down = false;
    paused = false;

    cr = 0.7;
    cg = 0.9;
    cb = 0.3;
    alpha = 0.03;

    add_amount = 0.1f * max(CELLS_X, CELLS_Y);
    source = 0.0f;
    force_y = -5.0f;
    force_x = 5.0f;

    fluid.init();
}

void display(void) {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glPushMatrix();

    int cells_y = (DISPLAY_KEY == 1) ? CELLS_Y + 1 : CELLS_Y;
    int cells_x = (DISPLAY_KEY == 2) ? CELLS_X + 1 : CELLS_X;

    float color;
    for (int y = 0; y < cells_y; ++y) {
        for (int x = 0; x < cells_x; ++x) {
            if (DISPLAY_KEY == 0) {
                color = fluid.S_at(y, x);
            } else if (DISPLAY_KEY == 1) {
                color = fabs(fluid.Uy_at(y, x));
            } else if (DISPLAY_KEY == 2) {
                color = fabs(fluid.Ux_at(y, x));
            }

            glColor4f(cr * color, cg * color, cb * color, alpha);
            glRectf((x - 1.0f) * 2.0f / (cells_x - 2) - 1.0f, (y - 0.5f) * 2.0f / (cells_y - 2) - 1.0f,
                    (x + 1.0f) * 2.0f / (cells_x - 2) - 1.0f, (y + 0.5f) * 2.0f / (cells_y - 2) - 1.0f);
        }
    }

    glFlush();
    glutSwapBuffers();
}

void reshape(int w, int h) {
    window_height = h;
    window_width = w;

    glViewport(0, 0, (GLsizei) w, (GLsizei) h);
}

void motion(int x, int y) {
    int cells_y = (DISPLAY_KEY == 1) ? CELLS_Y + 1 : CELLS_Y;
    int cells_x = (DISPLAY_KEY == 2) ? CELLS_X + 1 : CELLS_X;

    mouse_y = (int) (((float) (window_height - y) / window_height) * (cells_y));
    mouse_x = (int) (((float) x / window_width) * (cells_x));
}

void mouse(int button, int state, int x, int y) {
    switch (button) {
        case GLUT_LEFT_BUTTON:
            left_mouse_down = state == GLUT_DOWN;
            motion(x, y);
            break;
        default:
            break;
    }
}

void keyboard(unsigned char key, int x, int y) {
    switch (key) {
        case 'q':
            throw "exit";
            break;
        case ' ':
            paused = !paused;
            break;
        case '=':
            add_amount += 0.005f * max(CELLS_X, CELLS_Y);
            break;
        case '-':
            add_amount = fmax(0.0f, add_amount - 0.005f * max(CELLS_X, CELLS_Y));
            break;
        default:
            break;
    }
}

void special_keyboard(int key, int x, int y) {
    switch (key) {
        case GLUT_KEY_UP:
            force_y += 1.0f;
            break;
        case GLUT_KEY_DOWN:
            force_y -= 1.0f;
            break;
        case GLUT_KEY_LEFT:
            force_x -= 1.0f;
            break;
        case GLUT_KEY_RIGHT:
            force_x += 1.0f;
            break;
        default:
            break;
    }
}

void idle(void) {
    if (left_mouse_down) {
        fluid.add_source_at(mouse_y, mouse_x, add_amount);
    }
    if (!paused) {
        fluid.step(force_y, force_x, source);
    }
    if (force_y != 0) { force_y = 0; }
    if (force_x != 0) { force_x = 0; }
    if (source != 0) { source = 0; }
    glutPostRedisplay();
}

int main(int argc, char* argv[]) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(window_width, window_height);
    glutInitWindowPosition(WINDOW_X, WINDOW_Y);
    glutCreateWindow("Unstable Fluids");

    init();

    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutMouseFunc(mouse);
    glutPassiveMotionFunc(motion);
    glutMotionFunc(motion);
    glutKeyboardFunc(keyboard);
    glutSpecialFunc(special_keyboard);
    glutIdleFunc(idle);

    try {
        glutMainLoop();
    } catch (const char* msg) {
        if (CLEANUP) { fluid.cleanup(); }
        cout << "[-] Program terminated." << endl;
    }

    return 0;
}
