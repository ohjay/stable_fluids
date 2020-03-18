# Stable Fluids (CS 284B Project #1)
Reimplementation of Jos Stam's [Stable Fluids](https://d2f99xq7vri1nk.cloudfront.net/legacy_app_files/pdf/ns.pdf)
and Fattal et al.'s [Target-Driven Smoke Animation](http://www.cs.huji.ac.il/labs/cglab/projects/tdsmoke/tdsmoke.pdf).

![stable_2d](https://user-images.githubusercontent.com/8358648/50723766-d10deb80-1096-11e9-9ea2-509681ff63a9.gif)
![stable_3d](https://user-images.githubusercontent.com/8358648/50723767-d2d7af00-1096-11e9-935e-06322229eee8.gif)
![stable_face1](https://user-images.githubusercontent.com/8358648/50723768-d4a17280-1096-11e9-8ecb-68ae34f848b9.gif)
![stable_face2](https://user-images.githubusercontent.com/8358648/50723769-d66b3600-1096-11e9-931a-a6a85a6425ac.gif)
![stable_s1](https://user-images.githubusercontent.com/8358648/50723771-d9febd00-1096-11e9-956a-ff5c9198bf03.gif)
![stable_s2](https://user-images.githubusercontent.com/8358648/50723772-dbc88080-1096-11e9-88be-88fa2ccaf1b2.gif)

## Implementation and Usage Details (2D)
### Usage
```
make fresh; ./sim
```

### Controls
Click and drag to add mass-density and manipulate the velocity field.

| Input     | Effect                                                        |
| :-------: | :------------------------------------------------------------ |
| [         | Switch to previous color.                                     |
| ]         | Switch to next color.                                         |
| r         | Toggle rainbow mode (automatic color switching).              |
| t         | Switch to target-driven mode, with current density as target. |
| -         | Decrease source add amount.                                   |
| +         | Increase source add amount.                                   |
| Space     | Pause simulation.                                             |
| q, Ctrl-Q | Quit simulation.                                              |

### Parameters
All of our tunable parameters are stored within `params.h`. In particular, one
can select whether to visualize the density grid or a velocity grid
by modifying `DISPLAY_KEY` as follows:

```cpp
#define DISPLAY_KEY 0 // visualize density grid
#define DISPLAY_KEY 1 // visualize vertical velocity grid
#define DISPLAY_KEY 2 // visualize horizontal velocity grid
```

### Axis Directions
In our fluid representation, the y-axis points downward and the x-axis points
to the right. Unfortunately, the OpenGL / GLUT functions involve an
upward-pointing y-axis... which makes for some fun reversals in `main.cpp`.

### Grid
We represent our velocity fields using colocated grids. We also tried using
staggered MAC grids, but it made indexing annoying and didn't seem to improve
the visual quality of results.

**Previous staggered grid implementation:**
_take `(y, x)` to be the center of the cell.
Interpret horizontal velocity stored at `(y, x)` as velocity at `(y, x - 0.5)`.
Interpret vertical velocity stored at `(y, x)` as velocity at `(y - 0.5, x)`._

### Indexing
All of our grids are stored as linear arrays in row-major order. Therefore, to
index into one of our grids at position `(y, x)` we would actually access
position `y * CELLS_X + x`. When indexing, we will always involve the vertical
dimension before the horizontal dimension (à la NumPy).

We include the extra parameter `key` with many of our functions. This can take
on any of three values; the choice of which depends on the grid currently in use.

```cpp
0 // a scalar field
1 // a velocity grid containing vertical components
2 // a velocity grid containing horizontal components
```

## References
1. [Stable Fluids](https://d2f99xq7vri1nk.cloudfront.net/legacy_app_files/pdf/ns.pdf) (Stam)
2. [Target-Driven Smoke Animation](http://www.cs.huji.ac.il/labs/cglab/projects/tdsmoke/tdsmoke.pdf) (Fattal, Lischinski)
3. [Real-Time Fluid Dynamics for Games](https://pdfs.semanticscholar.org/847f/819a4ea14bd789aca8bc88e85e906cfc657c.pdf) (Stam)
