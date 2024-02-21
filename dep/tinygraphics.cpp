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

    void line(Vec2i t0, Vec2i t1, TGAImage &image, const TGAColor& color) {
        line(t0.x, t0.y, t1.x, t1.y, image, color);
    }

   
Vec3f barycentric(Vec2i *pts, Vec2i P) {
    // A = pts[0], B = pts[1], C = pts[2]
    // we can treat AB, AC as 2 base vectors
    // then we can express AP by the libear combination of these 2 basis:
    // u*AB + v*AC = AP
    // transform a little bit to yield:
    // u*AB + v*AC + PA = 0
    //
    // that is:
    // (u, v, 1) dot (ABx, ACx, PAx) = 0
    // (u, v, 1) dot (ABy, ACy, PAy) = 0
    // that is equivalent to say that the (u, v, 1) is perpendicular to both vectors
    // so we can do a cross product of these 2 vectors
    // then normalize the last component to yield (u, v, 1)
    auto& A = pts[0]; 
    auto& B = pts[1]; 
    auto& C = pts[2];
    Vec3f u = Vec3f(B.x-A.x, C.x-A.x, A.x-P.x)^Vec3f(B.y-A.y, C.y-A.y, A.y-P.y);
    // since u as coordinate should be integer, therefore u.z < 1.0 means that AP does not stick
    // out (enough) to make an eligible triangle (degenerate triangle), we ignore this case
    if (std::abs(u.z) < 1.0f) return Vec3f{-1.0f, -1.0f, -1.0f};
    auto ret = Vec3f{
        1.0f - (u.x + u.y)/u.z,
        u.x/u.z,
        u.y/u.z
    };
    return ret;
}

const int max_width = 1000;
const int max_height = 1000;
void triangle(Vec2i* t, TGAImage& image, const TGAColor& color) {
    // find bounding box first
    auto bottomLeft = Vec2i(INT_MAX, INT_MAX);
    auto topRight = Vec2i(0, 0);
    for (int i = 0; i < 3; i++) {
        auto& v = t[i];
        bottomLeft.x = std::max(std::min(bottomLeft.x, v.x), 0);
        bottomLeft.y = std::max(std::min(bottomLeft.y, v.y), 0);
        topRight.x = std::min(std::max(topRight.x, v.x), max_width);
        topRight.y = std::min(std::max(topRight.y, v.y), max_height);
    }

    for (int x = bottomLeft.x; x <= topRight.x; x++) {
        for (int y = bottomLeft.y; y <= topRight.y; y++) {
            auto u = barycentric(t, Vec2i{x, y});
            // since u shuold be treated as integers, then u.z < 1.0 means that AP does not stick
            // out (enough) to make an eligible triangle (degenerate triangle), we ignore this case
            //if (u.z < 1.0f) continue;
            //check if this point is outside the triangle
            if (u.x < 0.0f || u.y < 0.0f || u.z < 0.0f) continue;
            image.set(x, y, color);
        }
    }
}

void triangle(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage& image, const TGAColor& color) {
    Vec2i ts[3] = {t0, t1, t2};
    triangle(ts, image, color);
}


}