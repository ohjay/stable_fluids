# Stable Fluids (CS 284B Project #1)
Reimplementation of Jos Stam's [Stable Fluids](http://www.dgp.toronto.edu/people/stam/reality/Research/pdf/ns.pdf)
and Fattal et al.'s [Target-Driven Smoke Animation](http://www.cs.huji.ac.il/labs/cglab/projects/tdsmoke/tdsmoke.pdf).

## Implementation and Usage Details (2D)
### Usage
```
make fresh; ./sim
```

### Controls
| Input     | Effect                                                        |
| :-------: | :------------------------------------------------------------ |
| [         | Switch to previous color.                                     |
| ]         | Switch to next color.                                         |
| t         | Switch to target-driven mode, with current density as target. |
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
