#ifndef _my_graphics_
#define _my_graphics_

#include <dep/tgaimage.h>
#include <dep/geometry2.h>
namespace tinygraphics {

    const TGAColor white(255,255,255,255);
    const TGAColor red(255, 0, 0, 255);
    const TGAColor green(0, 255, 0, 255);
    const TGAColor blue(0, 0, 255, 255);
    const TGAColor yellow(255, 255, 153, 255);

    void line(int x0, int y0, int x1, int y1, TGAImage &image, const TGAColor& color);
    void line(Vec2i t0, Vec2i t1, TGAImage &image, const TGAColor& color);
    Vec3f barycentric(Vec2i *pts, Vec2i P);
    void triangle(Vec2i* t, TGAImage& image, const TGAColor& color);
    void triangle(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage& image, const TGAColor& color);
}

#endif