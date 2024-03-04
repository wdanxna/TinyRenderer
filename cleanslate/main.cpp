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

void triangle(Vec3f* v, const TGAColor& color, float* zbuffer, TGAImage& out) {
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
            
            if (zbuffer[x + y * width] - frag_vert.z < 1e-5) {
                zbuffer[x + y * width] = frag_vert.z;
                out.set(x, y, color);
            }
        }
    }
}

mat<4,4,float> viewport(int width, int height) {
    //viewport = scale * T
    mat<4,4,float> m;
    m[0] = {width * 0.5f, 0.f, 0.f,  width * 0.5f};
    m[1] = {0.f, height * 0.5f, 0.f,  height * 0.5f};
    m[2] = {0.f, 0.f,  0.5f, 0.5f};
    m[3] = {0.f, 0.f,  0.f,  1.f};
    return m;
}

int main() {

    TGAImage framebuffer(width, height, TGAImage::RGBA);
    TGAImage zimg(width, height, TGAImage::RGBA);

    mat<4,4,float> vp = viewport(width, height);

    float zbuffer[width*height];
    for (int i = 0; i < width*height;i++) {
        zbuffer[i] = std::numeric_limits<float>::min();
    }

    Model model("../res/african_head.obj");
    for (int i = 0; i < model.nfaces(); i++) {
        auto vert_ids = model.face(i);
        Vec3f object_verts[3];
        Vec3f screen_verts[3];
        for (int j = 0; j < 3; j++) {
            object_verts[j] = model.vert(vert_ids[j]);
            screen_verts[j] = proj<3>(vp * embed<4>(object_verts[j]));
        }

        //calculate normal
        Vec3f n = cross((object_verts[1] - object_verts[0]), (object_verts[2] - object_verts[0])).normalize();
        float shade = std::max(0.0f, n.z);
        triangle(screen_verts, TGAColor(255 * shade, 255* shade, 255* shade, 255), zbuffer, framebuffer);
    }

    //zbuffer visualization
    for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            zimg.set(x, y, TGAColor(255,255,255,255) * zbuffer[x + y*width]);
        }
    }
    zimg.flip_vertically();
    zimg.write_tga_file("z.tga");


    framebuffer.flip_vertically();//make (0,0) at bottom left, x going right, y going up
    framebuffer.write_tga_file("frame.tga");
    return 0;
}