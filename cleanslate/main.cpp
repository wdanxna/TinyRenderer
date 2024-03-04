#include <dep/tgaimage.h>
#include <dep/model.h>
const int width = 500;
const int height = 500;

Vec3f barycentric(Vec3f *v, Vec3f p) {
    Vec3f p0p1 = v[1] - v[0];
    Vec3f p0p2 = v[2] - v[0];
    Vec3f pp0 = v[0] - p;

    auto i = cross(Vec3f(p0p1.x, p0p2.x, pp0.x), 
                   Vec3f(p0p1.y, p0p2.y, pp0.y));
    i = i / i.z;
    return Vec3f(1.0-i.x-i.y, i.x, i.y);
}

void triangle(Vec3f* v, const TGAColor& color, TGAImage& out) {
    int left = INT_MAX;
    int right = 0;
    int bottom = INT_MAX;
    int up = 0;
    for (int i = 0; i < 3; i++) {
        left = std::min(width-1, std::max(0, std::min(left, int(v[i].x + 0.5))));
        right = std::min(width-1, std::max(0, std::max(right, int(v[i].x + 0.5))));
        up = std::min(height-1, std::max(0, std::max(up, int(v[i].y + 0.5))));
        bottom = std::min(height-1, std::max(0, std::min(bottom, int(v[i].y + 0.5))));
    }
    for (int x = left; x <= right; x++) {
        for (int y = bottom; y <= up; y++) {
            auto u = barycentric(v, {(float)x, (float)y, (float)0});
            if (u.x < 0 || u.y < 0 || u.z < 0) continue;

            mat<3,3,float> verts;
            verts.set_col(0, v[0]);
            verts.set_col(1, v[1]);
            verts.set_col(2, v[2]);
            auto frag_vert = verts * u;

            out.set(x, y, color);
        }
    }
}

int main() {

    TGAImage framebuffer(width, height, TGAImage::RGBA);

    Vec3f verts[3] = {
        {0, 0, 0},
        {100, 100, 0},
        {50, 200, 0}
    };

    triangle(verts, TGAColor(255, 0, 0, 255), framebuffer);

    framebuffer.flip_vertically();//make (0,0) at bottom left, x going right, y going up
    framebuffer.write_tga_file("frame.tga");
    return 0;
}