#include <iostream>
#include <dep/tgaimage.h>
#include <dep/model.h>
#include <dep/tinygraphics.h>

using namespace tinygraphics;
int main() {
    int width = 500;
    int height = 500;
    TGAImage img(width, height, TGAImage::RGB);

    Model model("../res/african_head.obj");

    for (int i = 0; i < model.nfaces(); ++i) {
        auto face = model.face(i);
        for (int j = 0; j < 3; j++) {
            Vec3f v0 = model.vert(face[j]);
            Vec3f v1 = model.vert(face[(j+1)%3]);
            int x0 = (v0.x+1.0f) * width/2.0f;
            int y0 = (v0.y+1.0f) * height/2.0f;
            int x1 = (v1.x+1.0f) * width/2.0f;
            int y1 = (v1.y+1.0f) * height/2.0f;
            line(x0, y0, x1, y1, img, white);
        }
    }

    img.flip_vertically();
    img.write_tga_file("face.tga");

    return 0;
}