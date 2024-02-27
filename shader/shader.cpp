#include <iostream>
#include <dep/tgaimage.h>
#include <dep/model.h>
#include <dep/tinygraphics.h>

using namespace tinygraphics;

Matrix embed(const Vec3f& v, float fill = 1.0f) {
    Matrix ret;
    ret[0][0] = v.x;
    ret[1][0] = v.y;
    ret[2][0] = v.z;
    ret[3][0] = fill;
    return ret;
}

Vec3f project(Matrix m) {
    Vec3f ret;
    ret.x = m[0][0];
    ret.y = m[1][0];
    ret.z = m[2][0];
    return ret;
}

Vec3f project2(Matrix m) {
    Vec3f ret;
    ret.x = m[0][0] / m[3][0];
    ret.y = m[1][0]/ m[3][0];
    ret.z = m[2][0]/ m[3][0];
    return ret;
}


struct Shader {
    Model* model;
    virtual Matrix vertex(int iface, int nthvert) = 0;
    virtual bool fragment(const Vec3f& barycentric, TGAColor& color) = 0;
};

struct GouraudShader : public Shader {
    Vec3f lightDir;
    Matrix modelView;
    Matrix modelView_invtrans;
    Matrix proj;

    Vec3f varying_intensity;//write by vertex shader, read by fragment shader

    Matrix vertex(int iface, int nthvert) override {
        const auto verts = model->face(iface);
        const auto [v, t, n] = verts[nthvert];

        Vec3f vertex = model->vert(v);
        Vec3f tex = model->tex(t);
        Vec3f norm = model->norm(n);
        //calculate lighting per vertex
        Vec3f normal_view = project(modelView_invtrans * embed(norm));
        Vec3f l = project(modelView * embed(lightDir*-1, 0.0f));
        normal_view.normalize();
        l.normalize();//turns out normalization is very important
        varying_intensity.raw[nthvert] = std::max(0.1f, normal_view * l);

        //do vertex transformation
        return proj * modelView * embed(vertex);
    }

    bool fragment(const Vec3f& barycentric, /*out*/TGAColor& color) override {
        //interpolate intensity passed by vertex shader
        float intensity = barycentric * varying_intensity;
        color = TGAColor(255, 255, 255, 255) * intensity;
        return true;
    }
};

const int max_width = 500;
const int max_height = 500;
void triangle(
    const Model& model, 
    int iface, 
    Shader& shader, 
    float* zbuffer,
    TGAImage& image) {

    //do vertex shader
    auto face = model.face(iface);
    Vec3f screen_coords[3];
    Vec3f tex_coords[3];
    Vec3f norm[3];
    for (int i = 0; i < 3; i++) {
        auto v = shader.vertex(iface, i);
        screen_coords[i] = Vec3f(
            (v.x() / v.w() + 1.0) * 501 / 2.0f,
            (v.y() / v.w() + 1.0) * 501 / 2.0f,
            (v.z() / v.w() + 1.0) * 255 / 2.0f
        );
    }


    // find bounding box first
    auto bottomLeft = Vec2i(INT_MAX, INT_MAX);
    auto topRight = Vec2i(0, 0);
    for (int i = 0; i < 3; i++) {
        auto& v = screen_coords[i];
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
                {(int)screen_coords[0].x, (int)screen_coords[0].y},
                {(int)screen_coords[1].x, (int)screen_coords[1].y},
                {(int)screen_coords[2].x, (int)screen_coords[2].y}
            };

            Vec3f u = barycentric(ti, Vec2i{x, y});
            //check if this point is outside the triangle
            if (u.x < 0.0f || u.y < 0.0f || u.z < 0.0f) continue;

            float z = 0.0;
            Vec3f texcoord, norm;
            for (int i = 0; i < 3; i++) {
                //interpolate z
                z +=  u.raw[i] * screen_coords[i].z;
            }

            if (zbuffer[y*image.get_width() + x] < z) {
                TGAColor color;
                if (shader.fragment(u, color)) {
                    zbuffer[y*image.get_width() + x] = z;
                    image.set(x, y, color);
                }
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

    //c is equivalent to the focal length (the distance between the pinhole and the projection plane)
    //that is, the smaller the c is the larger the FOV will be and strongger the perspective effect.
    float c = 0.5f;
    Matrix projection = Matrix::identity(4);
    projection[3][2] = -1.0f/c;

    // Vec3f camera{0.15f / 2, 0.05f / 2, 0.6f / 2};
    Vec3f camera{.3f, .4f, .5f};
    Matrix view = lookAt({0.0f, 0.0f, 0.001f}, camera, {0.0f, 1.0f, 0.0f});

    Matrix model2world = Matrix::identity(4); // transform from model space to world space

    Vec3f light{0,1, -1};

    GouraudShader shader;
    shader.model = &model;
    shader.lightDir = light;
    shader.modelView = view * model2world;
    shader.modelView_invtrans = (view * model2world).inverse().transpose();
    shader.proj = projection;

    for (int i = 0; i < model.nfaces(); ++i) {
        triangle(model, i, shader, zbuffer, img);
    }

    for (int j = 0; j < width*height; j++) {
       TGAColor d((int)zbuffer[j], (int)zbuffer[j], (int)zbuffer[j], 255);
       zimg.set(j % width, j / width, d);
    }
    zimg.flip_vertically();
    zimg.write_tga_file("z.tga");

    img.flip_vertically();
    img.write_tga_file("face.tga");

    return 0;
}