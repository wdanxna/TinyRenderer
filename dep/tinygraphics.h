#ifndef _my_graphics_
#define _my_graphics_

#include <dep/tgaimage.h>
namespace tinygraphics {
    
    const TGAColor white(255,255,255,255);
    const TGAColor red(255, 0, 0, 255);

    void line(int x0, int y0, int x1, int y1, TGAImage &image, const TGAColor& color);
}

#endif