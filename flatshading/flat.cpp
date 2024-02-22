#include <iostream>
#include <dep/tgaimage.h>
#include <dep/model.h>
#include <dep/tinygraphics.h>

using namespace tinygraphics;


const int max_width = 1000;
const int max_height = 1000;
void triangle(Vec3f* vert, float* zbuffer, Vec3f* texcoords, TGAImage& tex, float intensity, TGAImage& image) {
    // find bounding box first
    auto bottomLeft = Vec2i(INT_MAX, INT_MAX);
    auto topRight = Vec2i(0, 0);
    for (int i = 0; i < 3; i++) {
        auto& v = vert[i];
        bottomLeft.x = std::max(std::min(bottomLeft.x, (int)v.x), 0);
        bottomLeft.y = std::max(std::min(bottomLeft.y, (int)v.y), 0);
        topRight.x = std::min(std::max(topRight.x, (int)v.x), max_width);
        topRight.y = std::min(std::max(topRight.y, (int)v.y), max_height);
    }

    for (int x = bottomLeft.x; x <= topRight.x; x++) {
        for (int y = bottomLeft.y; y <= topRight.y; y++) {
            Vec2i ti[3] = {
                {(int)vert[0].x, (int)vert[0].y},
                {(int)vert[1].x, (int)vert[1].y},
                {(int)vert[2].x, (int)vert[2].y}
            };

            Vec3f u = barycentric(ti, Vec2i{x, y});
            //check if this point is outside the triangle
            if (u.x < 0.0f || u.y < 0.0f || u.z < 0.0f) continue;

            float z = 0.0;
            Vec3f texcoord;
            for (int i = 0; i < 3; i++) {
                //interpolate z
                z +=  u.raw[i] * vert[i].z;
                //interpolate tex coords
                for (int j = 0; j < 3; j++) {
                    texcoord.raw[i] += u.raw[j] * texcoords[j].raw[i];
                }
            }

            TGAColor color = tex.get(
                tex.get_width() * texcoord.x,
                tex.get_height() * texcoord.y
            );

            if (zbuffer[y*image.get_width() + x] < z) {
                zbuffer[y*image.get_width() + x] = z;
                image.set(x, y, TGAColor(
                    color.r * intensity,
                    color.g * intensity,
                    color.b * intensity,
                    color.a
                ));
            }
        }
    }
}

int main() {
    int width = 500;
    int height = 500;
    TGAImage img(width, height, TGAImage::RGB);
    TGAImage zimg(width, height, TGAImage::RGB);
    TGAImage texture;
    texture.read_tga_file("../res/african_head_diffuse.tga");
    texture.flip_vertically();

    float zbuffer[width * height];
    for (int j = 0; j < width*height; j++) {
        zbuffer[j] = std::numeric_limits<float>::min();
    }

    Model model("../res/african_head.obj");
    Vec3f light{0, 0, -1};

    for (int i = 0; i < model.nfaces(); ++i) {
        auto face = model.face(i);
        Vec3f screen_coords[3];
        Vec3f world_coords[3];
        Vec3f tex_coords[3];
       
        // foreach vertex in a triangle
        for (int j = 0; j < 3; j++) {
            const auto [vert_idx, tex_idx] = face[j];
            world_coords[j] = model.vert(vert_idx);
            tex_coords[j] = model.tex(tex_idx);
            //do model space to screen space transform
            screen_coords[j] = Vec3f(
                (world_coords[j].x + 1.0) * width / 2.0f,
                (world_coords[j].y + 1.0) * height / 2.0f,
                world_coords[j].z
            );
        }
        //Note: The direction of the normal vector matters!
        Vec3f n = (world_coords[2] - world_coords[0]) ^ (world_coords[1] - world_coords[0]);
        n.normalize();
        float shade = n*light;
        
        ///----- Not so good back surface removeal
        /// if (shade > 0)
        ///    triangle(screen_coords, img, color);
        ///------

        triangle(screen_coords, zbuffer, tex_coords, texture, shade, img);
    }

    for (int j = 0; j < width*height; j++) {
       TGAColor d((int)255*zbuffer[j], (int)255*zbuffer[j], (int)255*zbuffer[j], 255);
       zimg.set(j % width, j / width, d);
    }
    zimg.flip_vertically();
    zimg.write_tga_file("z.tga");

    img.flip_vertically();
    img.write_tga_file("face.tga");

    return 0;
}