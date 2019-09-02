#ifndef OUTPUT_PNG_HPP
#define OUTPUT_PNG_HPP

// Code for outputting the grid to PNG format.

struct xy_model_grid;

void xy_grid_to_png(const xy_model_grid &g, const char *fname);


#endif // OUTPUT_PNG_HPP
