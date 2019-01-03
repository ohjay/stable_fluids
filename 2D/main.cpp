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
int prev_mouse_y, prev_mouse_x;
int current_fluid;

Fluid fluid;
float add_amount;
double cr, cg, cb, alpha;
float fluid_colors[NUM_FLUIDS][3];

bool target_initialized;

void init(void) {
    glClearColor(0.0, 0.0, 0.0, 0.0);

    prev_mouse_y = 0;
    prev_mouse_x = 0;
    left_mouse_down = false;
    paused = false;

    cr = 0.7;
    cg = 0.9;
    cb = 0.3;
    alpha = 0.03;

    float colors[7][3] = ALL_COLORS; // {RED, GREEN, BLUE, YELLOW, CYAN, MAGENTA, WHITE};
    for (int i = 0; i < NUM_FLUIDS; ++i)
        memcpy(fluid_colors[i], colors[i], sizeof(colors[i]));

    add_amount = ADD_AMT_INIT * max(CELLS_X, CELLS_Y);
    current_fluid = 0;

    target_initialized = false;

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
                cr = 0.0f; cg = 0.0f; cb = 0.0f;
                for (int i = 0; i < NUM_FLUIDS; i++) {
                    cr += fluid_colors[i][0] * fluid.S_at(y, x, i);
                    cg += fluid_colors[i][1] * fluid.S_at(y, x, i);
                    cb += fluid_colors[i][2] * fluid.S_at(y, x, i);
                }
            } else {
                if (DISPLAY_KEY == 1) {
                    color = fabs(fluid.Uy_at(y, x));
                } else if (DISPLAY_KEY == 2) {
                    color = fabs(fluid.Ux_at(y, x));
                }
                cr = fluid_colors[0][0] * color;
                cg = fluid_colors[0][1] * color;
                cb = fluid_colors[0][2] * color;
            }

            glColor4f(cr, cg, cb, alpha);
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
    mouse_y = (int) (((float) (window_height - y) / window_height) * CELLS_Y);
    mouse_x = (int) (((float) x / window_width) * CELLS_X);
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
        case '[':
            current_fluid = max(0, current_fluid - 1);
            break;
        case ']':
            current_fluid = min(NUM_FLUIDS - 1, current_fluid + 1);
            break;
        case 't':
            if (!target_initialized) {
                fluid.save_density(current_fluid);
                target_initialized = true;
            }
            fluid.toggle_target_driven();
            break;
        default:
            break;
    }
}

void special_keyboard(int key, int x, int y) {
    switch (key) {
        case GLUT_KEY_UP:
            // force_y += 1.0f;
            break;
        case GLUT_KEY_DOWN:
            // force_y -= 1.0f;
            break;
        case GLUT_KEY_LEFT:
            // force_x -= 1.0f;
            break;
        case GLUT_KEY_RIGHT:
            // force_x += 1.0f;
            break;
        default:
            break;
    }
}

void idle(void) {
    if (left_mouse_down) {
        fluid.add_U_y_force_at(mouse_y, mouse_x, FORCE_SCALE * (mouse_y - prev_mouse_y));
        fluid.add_U_x_force_at(mouse_y, mouse_x, FORCE_SCALE * (mouse_x - prev_mouse_x));
        fluid.add_source_at(mouse_y, mouse_x, current_fluid, add_amount);
        prev_mouse_y = mouse_y;
        prev_mouse_x = mouse_x;
    }
    if (!paused) {
        fluid.step();
    }
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
