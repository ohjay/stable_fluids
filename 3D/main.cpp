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
int current_fluid;

Fluid fluid;
float source, add_amount, force_z, force_y, force_x;
int rx, ry, rz;
float zoom;
double cr, cg, cb, alpha;
float fluid_colors[NUM_FLUIDS][3];

void init(void) {
    glClearColor(0.0, 0.0, 0.0, 0.0);

    left_mouse_down = false;
    paused = false;

    cr = 0.7;
    cg = 0.9;
    cb = 0.3;
    alpha = 0.03;

    rx = 30; ry = 45; rz = 0;
    zoom = 1.0f;

    float colors[7][3] = ALL_COLORS; // {RED, GREEN, BLUE, YELLOW, CYAN, MAGENTA, WHITE};
    for (int i = 0; i < NUM_FLUIDS; ++i)
        memcpy(fluid_colors[i], colors[i], sizeof(colors[i]));

    add_amount = 0.1f * min(CELLS_X, min(CELLS_Y, CELLS_Z));
    current_fluid = 0;
    source = 10.0f;
    force_z = 5.0f;
    force_y = -5.0f;
    force_x = 5.0f;

    fluid.init();
}

void display(void) {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    int cells_z = (DISPLAY_KEY == 1) ? CELLS_Z + 1 : CELLS_Z;
    int cells_y = (DISPLAY_KEY == 2) ? CELLS_Y + 1 : CELLS_Y;
    int cells_x = (DISPLAY_KEY == 3) ? CELLS_X + 1 : CELLS_X;
    float color;

    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable( GL_BLEND );
    glPushMatrix();
    glScalef(zoom, zoom, zoom);
    glRotatef(rx, 1.0, 0.0, 0.0);
    glRotatef(ry, 0.0, 1.0, 0.0);
    glRotatef(rz, 0.0, 0.0, 1.0);
    glScalef(2.0f / cells_x, 2.0f / cells_y, 2.0f / cells_z);
    glTranslatef(-1.0f + 2.0f / cells_x / 2, -1.0f + 2.0f / cells_y / 2, -1.0f + 2.0f / cells_z / 2);
    // glTranslatef(-1, -1, -1);

    for (int z = 0; z < cells_z; ++z) {
        for (int y = 0; y < cells_y; ++y) {
            for (int x = 0; x < cells_x; ++x) {
                if (DISPLAY_KEY == 0) {
                    cr = 0.0f; cg = 0.0f; cb = 0.0f;
                    for (int i = 0; i < NUM_FLUIDS; i++) {
                        cr += fluid_colors[i][0] * fluid.S_at(z, y, x, i);
                        cg += fluid_colors[i][1] * fluid.S_at(z, y, x, i);
                        cb += fluid_colors[i][2] * fluid.S_at(z, y, x, i);
                    }
                } else {
                    if (DISPLAY_KEY == 0) {
                        color = fabs(fluid.Uz_at(z, y, x));
                    } else if (DISPLAY_KEY == 1) {
                        color = fabs(fluid.Uy_at(z, y, x));
                    } else if (DISPLAY_KEY == 2) {
                        color = fabs(fluid.Ux_at(z, y, x));
                    }
                    cr = fluid_colors[0][0] * color;
                    cg = fluid_colors[0][1] * color;
                    cb = fluid_colors[0][2] * color;
                }
                cout << fluid.S_at(z, y, x, 0) << endl;
                // glColor4f((x % 2) * 1.0f, 1.0f, 1.0f, 0.03f);
                // if (x == 0 || y == 0 || z == 0 || x == cells_x - 1 || y == cells_y - 1 || z == cells_z - 1)
                //     glColor4f(1.0f, 0.0f, 0.0f, 1.0f);
                // cout << cr << " " << cg << " " << cb << endl;
                // glColor4f(0.1 * x, 0.1 * x, 0.2 * x, 0.3);
                glColor4f(cr * 20, cg * 20, cb * 20, alpha);
                glutSolidCube(1.0f);  // scaled earlier in x, y, z

                glTranslatef(2.0f / cells_x, 0.0f, 0.0f);
            }
            glTranslatef(-2.0f, 0.0f, 0.0f);
            glTranslatef(0.0f, 2.0f / cells_y, 0.0f);
        }
        glTranslatef(0.0f, -2.0f, 0.0f);
        glTranslatef(0.0f, 0.0f, 2.0f / cells_z);
    }
    glPopMatrix();
    glPushMatrix();
    // glTranslatef(1.0f, 1.0f, -2.0f);
    glScalef(zoom, zoom, zoom);
    glRotatef(rx, 1.0, 0.0, 0.0);
    glRotatef(ry, 0.0, 1.0, 0.0);
    glRotatef(rz, 0.0, 0.0, 1.0);
    glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
    glutWireCube(0.7f);                 // where did this constant come from? D: dunno but it works
    glPopMatrix();

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
            // add_amount += 0.005f * max(CELLS_X, CELLS_Y);
            zoom += 0.1f;
            break;
        case '-':
            // add_amount = fmax(0.0f, add_amount - 0.005f * max(CELLS_X, CELLS_Y));
            zoom -= 0.1f;
            break;
        case '[':
            current_fluid = max(0, current_fluid - 1);
            break;
        case ']':
            current_fluid = min(NUM_FLUIDS - 1, current_fluid + 1);
            break;
        case ',':
            ry -= 5;
            break;
        case '.':
            ry += 5;
            break;
        default:
            break;
    }
}

void special_keyboard(int key, int x, int y) {
    switch (key) {
        case GLUT_KEY_UP:
            // force_y += 1.0f;
            rx += 5;
            break;
        case GLUT_KEY_DOWN:
            // force_y -= 1.0f;
            rx -= 5;
            break;
        case GLUT_KEY_LEFT:
            // force_x -= 1.0f;
            rz -= 5;
            break;
        case GLUT_KEY_RIGHT:
            // force_x += 1.0f;
            rz += 5;
            break;
        default:
            break;
    }
}

void idle(void) {
    if (source != 0) {
        for (int z = 0; z < CELLS_Z; ++z)
            for (int y = 0; y < CELLS_Y; ++y)
                for (int x = 0; x < CELLS_X; ++x)
                    fluid.add_source_at(z, y, x, current_fluid, add_amount);
        source = 0;
    }
    // glutPostRedisplay();
    // return;

    // if (left_mouse_down) {
    //     fluid.add_source_at(mouse_y, mouse_x, current_fluid, add_amount);
    // }
    if (!paused) {
        fluid.step(force_z, force_y, force_x, source);
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
