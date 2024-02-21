#include <dep/tgaimage.h>
#include <dep/model.h>
#include <dep/tinygraphics.h>
using namespace tinygraphics;

void outline(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage& image, const TGAColor& color) {
    line(t0, t1, image, color);
    line(t1, t2, image, color);
    line(t2, t0, image, color);
}


void triangle_refImp(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage &image, TGAColor color) { 
    if (t0.y==t1.y && t0.y==t2.y) return; // I dont care about degenerate triangles 
    // sort the vertices, t0, t1, t2 lower−to−upper (bubblesort yay!) 
    if (t0.y>t1.y) std::swap(t0, t1); 
    if (t0.y>t2.y) std::swap(t0, t2); 
    if (t1.y>t2.y) std::swap(t1, t2); 
    int total_height = t2.y-t0.y; 
    for (int i=0; i<total_height; i++) { 
        bool second_half = i>t1.y-t0.y || t1.y==t0.y; 
        int segment_height = second_half ? t2.y-t1.y : t1.y-t0.y; 
        float alpha = (float)i/total_height; 
        float beta  = (float)(i-(second_half ? t1.y-t0.y : 0))/segment_height; // be careful: with above conditions no division by zero here 
        Vec2i A =               t0 + (t2-t0)*alpha; 
        Vec2i B = second_half ? t1 + (t2-t1)*beta : t0 + (t1-t0)*beta; 
        if (A.x>B.x) std::swap(A, B); 
        for (int j=A.x; j<=B.x; j++) { 
            image.set(j, t0.y+i, color); // attention, due to int casts t0.y+i != A.y 
        } 
    } 
}

void triangle(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage& image, const TGAColor& color) {
    //bubble sort, make sure all vertices are sorted in ascending order by y value.
    if (t0.y > t1.y) std::swap(t0, t1);
    if (t1.y > t2.y) std::swap(t1, t2);
    if (t0.y > t1.y) std::swap(t0, t1);
    //draw bottom half 
    for (int y = t0.y; y < t1.y; y++) {
        float lr = (y - t0.y) / (float)(t2.y - t0.y);
        int lx = (1.0-lr)*t0.x + lr*t2.x;
        Vec2i left{lx, y};

        float rr = (y - t0.y) / (float)(t1.y - t0.y);
        int rx = (1.0-rr)*t0.x + rr*t1.x;
        Vec2i right(rx, y);

        line(left, right, image, color);
    }
    //draw upper half
    for (int y = t1.y; y <= t2.y; y++) {
        float lr = (y - t0.y) / (float)(t2.y - t0.y);
        int lx = (1.0-lr)*t0.x + lr*t2.x;
        Vec2i left{lx, y};

        float rr = (y - t1.y) / (float)(t2.y - t1.y + 1);//to prevent dividing by zero, maybe 1.0-alhpa is good enough?
        int rx = (1.0-rr)*t1.x + rr*t2.x;
        Vec2i right(rx, y);

        line(left, right, image, color);
    }
}


//refactor
void triangle2(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage& image, const TGAColor& color) {
    if (t0.y == t1.y && t1.y == t2.y) return;//skip degenerate triangle
    if (t0.y > t1.y) std::swap(t0, t1);
    if (t1.y > t2.y) std::swap(t1, t2);
    if (t0.y > t1.y) std::swap(t0, t1);

    int total_height = t2.y - t0.y;
    for (int i = 0; i <= total_height; ++i) {
        bool bottom_half = i < (t1.y - t0.y) || (t1.y == t2.y);
        float alpha = (float)(i - (bottom_half?0:(t1.y-t0.y))) / (bottom_half ? (t1.y - t0.y) : (t2.y - t1.y));
        float beta = (float)i / total_height;
        /*
            It would leave holes if we interpolate points like this, probabaly because floating
            point rounding issue? Im not sure.

            Vec2i pa = (bottom_half ? t0 : t1)*(1.0f - alpha) + (bottom_half ? t1 : t2)*alpha;
            Vec2i pb = t0*(1.0f - beta) + t2*beta;
        */
        Vec2i pa = {
            (int)((1.0f-alpha)*(bottom_half ? t0.x : t1.x) + alpha*(bottom_half?t1.x : t2.x)),
            t0.y + i
        };

        Vec2i pb = {
            (int)((1.0f-beta)*t0.x + beta*t2.x),
            t0.y + i
        };
        // image.set(pa.x, pa.y, yellow);
        // image.set(pb.x, pb.y, blue);
        line(pa, pb, image, color);
    }
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
    // since u shuold be treated as integers, therefore u.z < 1.0 means that AP does not stick
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
void triangle3(Vec2i* t, TGAImage& image, const TGAColor& color) {
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
           
            //if (u.z < 1.0f) continue;
            //check if this point is outside the triangle
            if (u.x < 0.0f || u.y < 0.0f || u.z < 0.0f) continue;
            image.set(x, y, color);
        }
    }
}

void triangle3(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage& image, const TGAColor& color) {
    Vec2i ts[3] = {t0, t1, t2};
    triangle3(ts, image, color);
}

#define TRIANGLE triangle3

int main() {
    TGAImage img(200, 200, TGAImage::RGB);

    Vec2i t0[3] = {{10, 70},{50, 160},{70, 80}};
    Vec2i t1[3] = {{180, 50},{150, 1},{70, 180}};
    Vec2i t2[3] = {{180, 150}, {120, 160}, {130, 180}};
    Vec2i t3[3] = {{50, 50}, {30, 90}, {70, 90}};
    
    TRIANGLE(t0[0], t0[1], t0[2], img, red);
    outline(t0[0], t0[1], t0[2], img, green);
   
    TRIANGLE(t1[0], t1[1], t1[2], img, white);
    outline(t1[0], t1[1], t1[2], img, white);
    
    TRIANGLE(t2[0], t2[1], t2[2], img, green);
    outline(t2[0], t2[1], t2[2], img, green);

    // outline(t3[0], t3[1], t3[2], img, red);
    // TRIANGLE(t3[0], t3[1], t3[2], img, yellow);
    

    img.flip_vertically();
    img.write_tga_file("tri.tga");
    return 0;
}