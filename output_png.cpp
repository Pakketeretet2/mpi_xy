#include "output_png.hpp"

#include <cstdlib>
#include <iostream>
#include <png.h>


#include "xy_grid.hpp"


void xy_grid_to_png(const xy_model_grid &g, const char *fname)
{
	// TODO: Check that these fit in an int?
	int img_width = g.Nx;
	int img_height = g.Ny;
	png_byte color_type, bit_depth;
	png_structp png_ptr;
	png_infop info_ptr;
	png_bytep row_ptr;

	color_type = PNG_COLOR_TYPE_RGB;
	bit_depth = 8;

	// Use C-style file manipulation:
	FILE *fp = fopen(fname, "wb");
	if (!fp) {
		std::cerr << "Failed to open file \"" << fname
		          << "\"for writing png!\n";
		return;
	}
	png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,
	                                  NULL, NULL, NULL);
	if (!png_ptr) {
		std::cerr << "Failed to create PNG r/w struct!\n";
		goto CLOSE_FILE;
		return;
	}
	info_ptr = png_create_info_struct(png_ptr);
	if (!info_ptr) {
		std::cerr << "Failed to create PNG info struct!\n";
		goto CLOSE_FILE;
		return;
	}
	if (setjmp(png_jmpbuf(png_ptr))) {
		std::cerr << "Error during init IO!\n";
		goto DESTROY_PNG_STRUCTS;
		return;
	}
	png_init_io(png_ptr, fp);

	
	if (setjmp(png_jmpbuf(png_ptr))) {
		std::cerr << "Error writing header!\n";
		goto DESTROY_PNG_STRUCTS;
		return;
	}
	png_set_IHDR(png_ptr, info_ptr, img_width, img_height,
	             bit_depth, color_type, PNG_INTERLACE_NONE,
	             PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
	png_write_info(png_ptr, info_ptr);

	// Right now the color scheme is simple.

	row_ptr = new png_byte[img_width*3];
	if (!row_ptr) {
		std::cerr << "Allocation of row buffer failed!\n";
		goto DESTROY_PNG_STRUCTS;
	}

	struct byte3 {
		double r,g,b;
	};
	
	for (long iy = 0; iy < g.Ny; ++iy) {
		for (long ix = 0; ix < g.Nx; ++ix) {
			long i = ix + g.Nx*iy;
			double val = g.s[i];
			if (val > 0) {
				*(row_ptr + 3*ix  ) = 255;
				*(row_ptr + 3*ix+1) = 0;
				*(row_ptr + 3*ix+2) = 0;
			} else {
				*(row_ptr + 3*ix  ) = 0;
				*(row_ptr + 3*ix+1) = 0;
				*(row_ptr + 3*ix+2) = 255;
			}
		}
		png_write_row(png_ptr, row_ptr);
	}
	png_write_end(png_ptr, NULL);

	if (row_ptr) delete [] row_ptr;
	
 DESTROY_PNG_STRUCTS:
	if (png_ptr && info_ptr) {
		png_destroy_write_struct(&png_ptr, &info_ptr);
	}

 CLOSE_FILE:
	fclose(fp);
}
