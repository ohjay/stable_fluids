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
float U_z_force;
int z_add_position;
int rx, ry, rz;
float zoom;
double cr, cg, cb, alpha;
float fluid_colors[NUM_FLUIDS][3];

void init(void) {
    glClearColor(0.0, 0.0, 0.0, 0.0);

    prev_mouse_y = 0;
    prev_mouse_x = 0;
    left_mouse_down = false;
    paused = false;

    cr = 0.7;
    cg = 0.9;
    cb = 0.3;
    alpha = 0.1;

    rx = 30; ry = 45; rz = 0;
    zoom = 1.0f;

    float colors[7][3] = ALL_COLORS;
    for (int i = 0; i < NUM_FLUIDS; ++i)
        memcpy(fluid_colors[i], colors[i], sizeof(colors[i]));

    add_amount = ADD_AMT_INIT * max(max(CELLS_X, CELLS_Y), CELLS_Z);
    current_fluid = 0;
    U_z_force = 0.0f;
    z_add_position = (int) (((float) CELLS_Z) / 2.0f);

    fluid.init();
}

void display(void) {
    float cellstep = 10.0f;
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_BLEND);
    glPushMatrix();
    glScalef(zoom, zoom, zoom);
    glRotatef(rx, 1.0, 0.0, 0.0);
    glRotatef(ry, 0.0, 1.0, 0.0);
    glRotatef(rz, 0.0, 0.0, 1.0);
    glScalef(2.0f / CELLS_X, 2.0f / CELLS_Y, 2.0f / CELLS_Z);
    glTranslatef(-5.0f, -5.0f, -5.0f);

    // draw colored voxels (3D grid)
    float color, _alpha, total_S;
    for (int z = 0; z < CELLS_Z; ++z) {
        for (int y = 0; y < CELLS_Y; ++y) {
            for (int x = 0; x < CELLS_X; ++x) {
                if (DISPLAY_KEY == 0) {
                    total_S = 0.0f;
                    cr = 0.0f; cg = 0.0f; cb = 0.0f;
                    for (int i = 0; i < NUM_FLUIDS; i++) {
                        cr += fluid_colors[i][0] * fluid.S_at(z, y, x, i);
                        cg += fluid_colors[i][1] * fluid.S_at(z, y, x, i);
                        cb += fluid_colors[i][2] * fluid.S_at(z, y, x, i);
                        total_S += fluid.S_at(z, y, x, i);
                    }
                } else {
                    total_S = 1.0f;
                    if (DISPLAY_KEY == 1) {
                        color = fabs(fluid.Uz_at(z, y, x));
                    } else if (DISPLAY_KEY == 2) {
                        color = fabs(fluid.Uy_at(z, y, x));
                    } else if (DISPLAY_KEY == 3) {
                        color = fabs(fluid.Ux_at(z, y, x));
                    }
                    cr = fluid_colors[current_fluid][0] * color;
                    cg = fluid_colors[current_fluid][1] * color;
                    cb = fluid_colors[current_fluid][2] * color;
                }
                // without smart alpha blending,
                // black voxels in the front cover colored voxels in the back
                if (ALPHA_OPTION == 2) {
                    _alpha = fmin(alpha, alpha * pow(total_S, 2) * 100);
                } else {
                    _alpha = fmin(alpha, alpha * pow(total_S, 3) * 1e4);
                }
                glColor4f(cr * COLOR_SCALE, cg * COLOR_SCALE, cb * COLOR_SCALE, _alpha);
                glutSolidCube(1.0f);  // scaled earlier in x, y, z

                glTranslatef(cellstep / CELLS_X, 0.0f, 0.0f);
            }
            glTranslatef(-cellstep, 0.0f, 0.0f);
            glTranslatef(0.0f, cellstep / CELLS_Y, 0.0f);
        }
        glTranslatef(0.0f, -cellstep, 0.0f);
        glTranslatef(0.0f, 0.0f, cellstep / CELLS_Z);
    }
    glPopMatrix();

    // draw wire cube
    glPushMatrix();
    glScalef(zoom, zoom, zoom);
    glRotatef(rx, 1.0, 0.0, 0.0);
    glRotatef(ry, 0.0, 1.0, 0.0);
    glRotatef(rz, 0.0, 0.0, 1.0);
    glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
    glutWireCube(0.7f);
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
    mouse_y = (int) (((float) (window_height - y) / window_height) * CELLS_Y);
    mouse_x = (int) (((float) x / window_width) * CELLS_X);
}

void mouse(int button, int state, int x, int y) {
    switch (button) {
        case GLUT_LEFT_BUTTON:
            left_mouse_down = state == GLUT_DOWN;
            motion(x, y);
            if (left_mouse_down) {
                prev_mouse_y = mouse_y;
                prev_mouse_x = mouse_x;
            }
            break;
        default:
            break;
    }
}

void keyboard(unsigned char key, int x, int y) {
    switch (key) {
        case 'a':
            z_add_position = max(1, z_add_position - 1);
            break;
        case 'd':
            z_add_position = min(CELLS_Z - 2, z_add_position + 1);
            break;
        case 'q':
            throw "exit";
            break;
        case ' ':
            paused = !paused;
            break;
        case '=':
            zoom += 0.1f;
            break;
        case '-':
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
            rx += 5;
            break;
        case GLUT_KEY_DOWN:
            rx -= 5;
            break;
        case GLUT_KEY_LEFT:
            rz -= 5;
            break;
        case GLUT_KEY_RIGHT:
            rz += 5;
            break;
        default:
            break;
    }
}

void idle(void) {
    if (left_mouse_down) {
        U_z_force = (mouse_y - prev_mouse_y) * (mouse_x - prev_mouse_x);
        fluid.add_U_z_force_at(z_add_position, mouse_y, mouse_x, FORCE_SCALE * U_z_force);
        fluid.add_U_y_force_at(z_add_position, mouse_y, mouse_x, FORCE_SCALE * (mouse_y - prev_mouse_y));
        fluid.add_U_x_force_at(z_add_position, mouse_y, mouse_x, FORCE_SCALE * (mouse_x - prev_mouse_x));
        fluid.add_source_at(z_add_position, mouse_y, mouse_x, current_fluid, add_amount);
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
