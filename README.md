# Stable Fluids (CS 284B Project #1)

## Implementation and Usage Details (2D)
### Parameters
All of our tunable parameters are stored within `params.h`. In particular, one
can select whether to visualize the density grid or the velocity grid
by modifying `DISPLAY_KEY` as follows:

```cpp
#define DISPLAY_KEY 0 // visualize density grid
#define DISPLAY_KEY 1 // visualize vertical velocity grid
#define DISPLAY_KEY 2 // visualize horizontal velocity grid
```

### Axis Directions
The y-axis points downward and the x-axis points to the right.

### Velocity Fields
We represent our velocity fields using staggered MAC grids. **If we take `(y, x)`
to mean the center of the cell**, then the horizontal velocity stored at `(y, x)`
is really the velocity at `(y, x - 0.5)`. Likewise, the vertical velocity stored
at `(y, x)` is really the velocity at `(y - 0.5, x)`.

### Indexing
All of our grids are stored as linear arrays in row-major order. Therefore, to
index into one of our grids at position `(y, x)` we would actually access
position `y * CELLS_X + x`. When indexing, we will always involve the vertical
dimension before the horizontal dimension (à la NumPy).

Also, since indexing is different for each of our grids (e.g. the staggered
horizontal grid has one extra value per row), we include the extra parameter
`key` with many of our functions. `key` can take on any of three values; the
choice of which depends on the grid currently in use.

```cpp
0 // a scalar field
1 // a velocity grid containing vertical components
2 // a velocity grid containing horizontal components
```
