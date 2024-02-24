#include <iostream>
#include <dep/tgaimage.h>
#include <dep/model.h>
#include <dep/tinygraphics.h>

using namespace tinygraphics;


const int max_width = 500;
const int max_height = 500;
void triangle(Vec3f* vert, float* zbuffer, Vec3f* texcoords, TGAImage& tex, float intensity, TGAImage& image) {
    // find bounding box first
    auto bottomLeft = Vec2i(INT_MAX, INT_MAX);
    auto topRight = Vec2i(0, 0);
    for (int i = 0; i < 3; i++) {
        auto& v = vert[i];
        bottomLeft.x = std::max(std::min(bottomLeft.x, (int)v.x), 0);
        bottomLeft.x = std::min(bottomLeft.x, max_width);
        bottomLeft.y = std::max(std::min(bottomLeft.y, (int)v.y), 0);
        bottomLeft.y = std::min(bottomLeft.y, max_height);

        topRight.x = std::min(std::max(topRight.x, (int)v.x), max_width);
        topRight.x = std::max(topRight.x, 0);
        topRight.y = std::min(std::max(topRight.y, (int)v.y), max_height);
        topRight.y = std::max(topRight.y, 0);
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

Matrix lookAt(const Vec3f& center, const Vec3f& eye, const Vec3f& up) {
    //construct an orthognal coordinate system which origin is at 'eye'
    // and z is eye-center
    Vec3f z = (eye - center).normalize();
    Vec3f x = (up ^ z).normalize();
    Vec3f y = (z ^ x).normalize();

    //note that the matrix representation in this project is row major, that 
    //is, m[i] is a row, m[0..4][0] is the first column.
    //to construct a camera transformation matrix which transform any point p
    //in original system to the camera system. we can consider the fact that the
    //new coordinate of p' is equivalent to project the p onto each basis of the new coord system.
    //to project any vector onto another vector, we can just do a simple dot product, that is
    //equivalent to make the x^,y^,z^ of camera coordinate system as the first 3 row of the matrix. 
    Matrix view = Matrix::identity(4);
    view[0][0] = x.x;
    view[0][1] = x.y;
    view[0][2] = x.z;

    view[1][0] = y.x;
    view[1][1] = y.y;
    view[1][2] = y.z;

    view[2][0] = z.x;
    view[2][1] = z.y;
    view[2][2] = z.z;

    Matrix translate = Matrix::identity(4);
    translate[0][3] = -eye.x;
    translate[1][3] = -eye.y;
    translate[2][3] = -eye.z;

    return view * translate;
}

int main() {
    int width = 501;
    int height = 501;
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

    //c is equivalent to the focal length (the distance between the pinhole and the projection plane)
    //that is, the smaller the c is the larger the FOV will be and strongger the perspective effect.
    float c = 0.5f;

    Matrix projection;
    projection[0][0] = 1.0f;
    projection[1][1] = 1.0f;
    projection[2][2] = 1.0f;
    projection[3][3] = 1.0f;
    projection[3][2] = -1.0f/c;

    Vec3f camera{0.15f, 0.05f, 0.6f};
    Matrix view = lookAt({0.0f, 0.0f, -0.01f}, camera, {0.0f, 1.0f, 0.0f});

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
            
            //project the world coordinates onto the z=0 plane
            Matrix c;
            c[0][0] = world_coords[j].raw[0];
            c[1][0] = world_coords[j].raw[1];
            c[2][0] = world_coords[j].raw[2];
            c[3][0] = 1.0f;

            Matrix p = projection * view * c;
            world_coords[j].x = p[0][0] / p[3][0];
            world_coords[j].y = p[1][0] / p[3][0];

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