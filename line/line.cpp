#include <dep/tgaimage.h>

const TGAColor white = TGAColor(255,255,255,255);
const TGAColor red = TGAColor(255,0,0,255);

void line(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor color) {
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

    for (int x = x0; x <= x1; x++) {
        float f = (x-x0) / (float)(x1-x0);
        int y = y0*(1.0-f) + y1*f;
        if (steep) {
            image.set(y, x, color);
        } else {
            image.set(x, y, color);
        }
    }
}

//optimized version 1:
//take the division out of the loop. by defineing a slope and an error variable
//and increment the error by slope in every iteration, the error variable keeps
//track of the difference between the current y value and what's the value it should be.
//since the y is an integer, we must decide when that increment happen. It is reasonable
//to choose error > 0.5 as the threashold, its just like round up operation.
//keep mind of the slope calculation which is wrapped by std::abs since y1 < y0 could happen.
void line2(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor color) {
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
    float slope = std::abs(dy / (float)dx);
    float error = .0f;
    int y = y0;
    for (int x = x0; x <= x1; x++) {
        if (steep) {
            image.set(y, x, color);
        } else {
            image.set(x, y, color);
        }

        error += slope;
        while (error > 0.5) {
            y += (y0 < y1) ? 1 : -1;
            error -= 1.0;
        }
    }
}

//optimized version 2:
//remove floating point and division.
//in the previous version, the error = n * dy/dx (n is the iteration count)
//the threashold check: error > 0.5 
//can be rewrite to: n*dy/dx > 1/2
//transform a little bit: n * 2 * dy/dx > 1
//then: n * 2*dy > dx
//if we defined slop' = 2*dy, and theashold check at dx, we can eliminate division and floating point use.
//error' = n*2*dy = 2*dx*error. when error = 1, error' should be 2*dx*1 = 2*dx. this 2*dx should be substract
//from error' whenever threashold is reached.
void line3(int x0, int y0, int x1, int y1, TGAImage &image, const TGAColor& color) {
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

#define LINE line3

int main () {
    TGAImage image(100, 100, TGAImage::RGB);
    LINE(13, 20, 80, 40, image, white);
    LINE(20, 13, 40, 80, image, red);
    LINE(80, 40, 13, 20, image, red);//<-- test for symmetry

    image.flip_vertically();
    image.write_tga_file("line3.tga");
    return 0;
}