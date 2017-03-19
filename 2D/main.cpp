#include "Fluid.h"

int window_width = 600;
int window_height = 600;

Fluid fluid;
float F[2];
float Ssource, add_amount;
int Fy, Fx;
bool left_mouse_down, paused;
float cr, cg, cb, alpha;

void init(void) {
    glClearColor(0.0, 0.0, 0.0, 0.0);
    // F[0] = -6.0f; F[1] = 6.0f;
    F[0] = 2.0f; F[1] = 0.0f;
    Fy = CELLS_PER_SIDE / 2;
    Fx = CELLS_PER_SIDE / 2;
    Ssource = 0.0f;
    cr = 0.7, cg = 0.2, cb = 0.9;
    alpha = 0.03;
    paused = false;
    add_amount = 200.0f;
}

// everything we need in order to redraw the scene
void display(void) {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glPushMatrix();

    // draw the density grid
    float x, y, c00, c01, c10, c11;

    int h = CELLS_PER_SIDE;
    int w = CELLS_PER_SIDE;
    float vsize = fluid.grid_spacing();
    float half_side = CELLS_PER_SIDE / 2.0f;

    glBegin(GL_QUADS);

    for (int r = 0; r < h; ++r) {
        y = (r + 0.5f) * vsize;
        for (int c = 0; c < w; ++c) {
            x = (c + 0.5f) * vsize;

            // (TODO) these 1500s should not be here; they're just here for debugging help
            c00 = fluid.S_at(r, c) / 500.0f;
            c01 = c10 = c11 = c00;
            c00 = c01 = fluid.U10_at(r, c);
            c10 = c11 = fluid.U11_at(r, c);
            // c01 = fluid.S_at(r, c + 1) / 500.0f;
            // c10 = fluid.S_at(r + 1, c) / 500.0f;
            // c11 = fluid.S_at(r + 1, c + 1) / 500.0f;

            glColor4f(cr*c11, cg*c11, cb*c11, alpha); glVertex2f((x + vsize) / half_side - 1.0f, (y + vsize) / half_side - 1.0f);
            glColor4f(cr*c01, cg*c01, cb*c01, alpha); glVertex2f(x / half_side - 1.0f, (y + vsize) / half_side - 1.0f);
            glColor4f(cr*c00, cg*c00, cb*c00, alpha); glVertex2f(x / half_side - 1.0f, y / half_side - 1.0f);
            glColor4f(cr*c10, cg*c10, cb*c10, alpha); glVertex2f((x + vsize) / half_side - 1.0f, y / half_side - 1.0f);
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
            left_mouse_down = state == GLUT_DOWN;
            // if (state == GLUT_DOWN) {
            //     // do something here, like create a new object
            //     // F[0] += 12.0f;
            //     // F[1] += 15.0f;
            //     // Ssource += 120.0f;
            //
            //     // (TODO) should really use SIDE_LEN
            //     Fy = (int) (((float) (window_width - y) / window_width) * CELLS_PER_SIDE);
            //     Fx = (int) (((float) x / window_width) * CELLS_PER_SIDE);
            //     fluid.add_S_at(Fy, Fx, 500.0f);
            // }
            break;
        default:
            break;
    }
}

void motion(int x, int y) {
    // do something here, like apply a force
    Fy = (int) (((float) (window_width - y) / window_width) * CELLS_PER_SIDE);
    Fx = (int) (((float) x / window_width) * CELLS_PER_SIDE);
}

void keyboard(unsigned char key, int x, int y) {
    switch (key) {
        case 'q':
            // do something here, like quit
            throw "exit";
            break;
        case ' ':
            paused = !paused;
            break;
        case '=':
            add_amount += 100.0f;
            break;
        case '-':
            add_amount = max(0.0f, add_amount - 100.0f);
            break;
        default:
            break;
    }
}

// render the next frame of our simulation
void idle(void) {
    if (left_mouse_down)
        fluid.add_S_at(Fy, Fx, add_amount);
    if (!paused) fluid.step(F, Ssource, Fy, Fx);
    sleep(1);
    if (F[0] != 0) { F[0] = 0; }
    if (F[1] != 0) { F[1] = 0; }
    if (Ssource != 0) { Ssource = 0; }
    glutPostRedisplay();
}

int main(int argc, char* argv[]) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(window_width, window_height);
    glutInitWindowPosition(400, 100);
    glutCreateWindow("Stable Fluids");

    init();
    fluid.init(VISC, KS, AS, DT);

    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutMouseFunc(mouse);
    glutPassiveMotionFunc(motion);
    glutMotionFunc(motion);
    glutKeyboardFunc(keyboard);
    glutIdleFunc(idle);

    try {
        glutMainLoop();
    } catch (const char* msg) {
        fluid.cleanup();
        std::cout << "[-] Program terminated." << std::endl;
    }

    return 0;
}
