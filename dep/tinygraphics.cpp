#include <dep/tinygraphics.h>
namespace tinygraphics {
void line(int x0, int y0, int x1, int y1, TGAImage &image, const TGAColor& color) {
    bool steep = false;
    if (abs(x1 - x0) < abs(y1 - y0)) {
        //steep line, tanspose
        steep = true;
        std::swap(x0, y0);
        std::swap(x1, y1);
    }

    if (x0 > x1) {
        std::swap(x0, x1);
        std::swap(y0, y1);
    }

    int dx = x1 - x0;
    int dy = y1 - y0;
    int derr = std::abs(dy << 1);
    int error = 0;
    int y = y0;
    for (int x = x0; x <= x1; x++) {
        if (steep) {
            image.set(y, x, color);
        } else {
            image.set(x, y, color);
        }

        error += derr;
        while (error > dx) {
            y += (y0 < y1) ? 1 : -1;
            error -= (dx<<1);
        }
    }
}
}