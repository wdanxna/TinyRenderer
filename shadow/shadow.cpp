#include <iostream>
#include <dep/tgaimage.h>
#include <dep/model.h>
#include <dep/tinygraphics.h>

using namespace tinygraphics;

Matrix lookAt(const Vec3f& center, const Vec3f& eye, const Vec3f& up) {
    //construct an orthognal coordinate system which origin is at 'eye'
    // and z is eye-center
    Vec3f z = (eye - center).normalize();
    Vec3f x = cross(up, z).normalize();
    Vec3f y = cross(z, x).normalize();

    //note that the matrix representation in this project is row major, that 
    //is, m[i] is a row, m[0..4][0] is the first column.
    //to construct a camera transformation matrix which transform any point p
    //in original system to the camera system. we can consider the fact that the
    //new coordinate of p' is equivalent to project the p onto each basis of the new coord system.
    //to project any vector onto another vector, we can just do a simple dot product, that is
    //equivalent to make the x^,y^,z^ of camera coordinate system as the first 3 row of the matrix. 
    Matrix view = Matrix::identity();
    view[0][0] = x.x;
    view[0][1] = x.y;
    view[0][2] = x.z;

    view[1][0] = y.x;
    view[1][1] = y.y;
    view[1][2] = y.z;

    view[2][0] = z.x;
    view[2][1] = z.y;
    view[2][2] = z.z;

    Matrix translate = Matrix::identity();
    translate[0][3] = -eye.x;
    translate[1][3] = -eye.y;
    translate[2][3] = -eye.z;

    return view * translate;
}

struct Shader {
    Model* model;
    virtual Vec4f vertex(int iface, int nthvert) = 0;
    virtual bool fragment(const Vec3f& barycentric, TGAColor& color) = 0;
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
            (v[0] / v[3] + 1.0) * 501 / 2.0f,
            (v[1] / v[3] + 1.0) * 501 / 2.0f,
            (v[2] / v[3] + 1.0) * 255 / 2.0f
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
                z +=  u[i] * screen_coords[i].z;
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



struct PhongShader : public Shader {
    Vec3f lightDir;
    Matrix modelView;
    Matrix modelView_invtrans;
    Matrix projection;
    Matrix shadowM;

    TGAImage* texture;
    TGAImage* normalmap;
    TGAImage* specularmap;
    TGAImage* shadowmap;

    Vec3f varying_intensity;//write by vertex shader, read by fragment shader
    mat<2, 3, float> varying_uv;
    mat<3, 3, float> varying_norm, varying_pos;

    Vec4f vertex(int iface, int nthvert) override {
        const auto verts = model->face(iface);
        const auto [v, t, n] = verts[nthvert];

        Vec3f vertex = model->vert(v);
        Vec3f tex = model->tex(t);
        Vec3f norm = proj<3>(modelView_invtrans * embed<4>(model->norm(n)));

        varying_uv.set_col(nthvert, proj<2>(tex));
        varying_norm.set_col(nthvert, norm);

        Vec4f mv = modelView * embed<4>(vertex);
        varying_pos.set_col(nthvert, proj<3>(mv));
        //do vertex transformation
        Vec4f clip = projection * mv;
        return clip;
    }

    bool fragment(const Vec3f& barycentric, /*out*/TGAColor& color) override {
        //interpolate intensity passed by vertex shader
        Vec3f l = proj<3>(modelView * embed<4>(lightDir*-1, 0.0f)).normalize();
        //interpolate uv
        Vec2f uv = varying_uv * barycentric;

        //find tangent space
        Vec3f bn = (varying_norm * barycentric).normalize();
        mat<3,3,float> A;
        A[0] = varying_pos.col(1) - varying_pos.col(0);
        A[1] = varying_pos.col(2) - varying_pos.col(0);
        A[2] = bn;
        mat<3,3,float> AI = A.invert();

        Vec3f i = AI * Vec3f(varying_uv[0][1] - varying_uv[0][0], varying_uv[0][2] - varying_uv[0][0], 0.0f);
        Vec3f j = AI * Vec3f(varying_uv[1][1] - varying_uv[1][0], varying_uv[1][2] - varying_uv[1][0], 0.0f);
        mat<3,3,float> tangentM;
        tangentM.set_col(0, i.normalize());
        tangentM.set_col(1, j.normalize());
        tangentM.set_col(2, bn);

        TGAColor n_sample = normalmap->get(
            int(uv.x * normalmap->get_width()), 
            int(uv.y * normalmap->get_height()));
        Vec3f n;
        for (int i=0; i<3; i++)
            n[2-i] = (float)n_sample.raw[i]/255.f*2.f - 1.f;
        n = (tangentM * n).normalize();

        //specular map
        TGAColor spec_sample = specularmap->get(
            int(uv.x * specularmap->get_width()), 
            int(uv.y * specularmap->get_height()));
        float pow = spec_sample.raw[0] / 1.0f;
        //calculate reflect vector
        Vec3f r = (n*(n*l*2.0f)-l).normalize();
        //dot product between r and view vector (0, 0, 1), why (0,0,1)?
        //since after camera transform, we always located at (0,0,c) and looking
        //at (0, 0, -1), the view vector is just -(0, 0, -1) -> (0, 0, 1)
        float spec = powf(std::max(r.z, 0.0f), std::max(1.0f, pow));

        //calc shadow
        Vec4f pos = embed<4>(varying_pos * barycentric);//<-- view space vertex pos
        Vec4f p = shadowM * pos;
        int shadow_val = shadowmap->get((p[0]/p[3]+1.0)/2.0f * shadowmap->get_width(), 
                           (p[1]/p[3]+1.0)/2.0f * shadowmap->get_height()).r;
        int z_val = std::max(0, std::min(255, (int)(255*(p[2]+1.0f)/2.0f)));
        float shadow = 0.3f + 0.7 * (shadow_val < (2+z_val));

        TGAColor diffuse = texture->get(
            int(uv.x * texture->get_width()), 
            int(uv.y * texture->get_height())) * std::max(0.0f, n*l) * shadow * 1.2;

        int ambient = 20;

        color = TGAColor(
            (int)std::min<float>(ambient + diffuse.r + diffuse.r * spec, 255),
            (int)std::min<float>(ambient + diffuse.g + diffuse.g * spec, 255),
            (int)std::min<float>(ambient + diffuse.b + diffuse.b * spec, 255),
            255
        );
        return true;
    }
};

class ShadowShader : public Shader {
public:
    //uniforms
    Matrix mvpM;
    //varyings
    mat<3,3,float> varying_pos;

    ShadowShader(Matrix &mvp) : mvpM{mvp} {}
    virtual Vec4f vertex(int iface, int nthvert) override {
        const auto& face = model->face(iface);
        auto [vi, ti, ni] = face[nthvert];
        auto vertex = model->vert(vi);
        auto texcoord = model->tex(ti);
        auto norm = model->norm(ni);
        Vec4f clip_pos = mvpM * embed<4>(vertex);

        varying_pos.set_col(nthvert, proj<3>(clip_pos));
        return clip_pos;
    }

    virtual bool fragment(const Vec3f& barycentric, TGAColor& color) override {
        //interpolate the ndc vertex coordinates
        Vec3f v = varying_pos * barycentric;
        color = TGAColor(255, 255, 255, 255) * std::max(0.0f, ((v.z + 1.f)/2.0f));
        return true;
    }
};


int main() {
    int width = 501;
    int height = 501;
    TGAImage img(width, height, TGAImage::RGB);
    TGAImage zimg(width, height, TGAImage::RGB);

    TGAImage texture;
    texture.read_tga_file("../res/diablo3_pose_diffuse.tga");
    texture.flip_vertically();

    TGAImage normalmap;
    normalmap.read_tga_file("../res/diablo3_pose_nm_tangent.tga");
    normalmap.flip_vertically();

    TGAImage specularmap;
    specularmap.read_tga_file("../res/diablo3_pose_spec.tga");
    specularmap.flip_vertically();

    float zbuffer[width * height];
    for (int j = 0; j < width*height; j++) {
        zbuffer[j] = std::numeric_limits<float>::min();
    }

    Model model("../res/diablo3_pose.obj");

    //c is equivalent to the focal length (the distance between the pinhole and the projection plane)
    //that is, the smaller the c is the larger the FOV will be and strongger the perspective effect.
    float c = 1.0f;
    Matrix projection = Matrix::identity();
    projection[3][2] = -1.0f/c;

    // Vec3f camera{0.25f, 0.05f, 0.6f};
    // Vec3f camera{.3f, .2f, 0.2f};
    Vec3f camera{0.0, 0.0f, 0.3f};
    // Vec3f camera{-.3f, .2f, 0.2f};
    Matrix view = lookAt({0.0f, 0.0f, 0.01f}, camera, {0.0f, 1.0f, 0.0f});

    Matrix model2world = Matrix::identity(); // transform from model space to world space


    Vec3f light_pos{.3f, .2f, 0.2f};
    Vec3f light_dir{0, 0, -0.1};
    //shadow pass
    Matrix shadow_mv = lookAt(light_dir, light_pos, {0, 1.0, 0});
    Matrix shadow_mvp = projection * shadow_mv;
    ShadowShader shadowS(shadow_mvp);
    shadowS.model = &model;
    TGAImage shadowmap(width, height, TGAImage::RGB);
    for (int i = 0; i < model.nfaces(); ++i) {
        triangle(model, i, shadowS, zbuffer, shadowmap);
    }


    for (int j = 0; j < width*height; j++) {
        zbuffer[j] = std::numeric_limits<float>::min();
    }
    //main pass
    PhongShader shader;
    shader.model = &model;
    shader.lightDir = light_dir - light_pos;
    shader.modelView = view * model2world;
    shader.modelView_invtrans = (view * model2world).invert_transpose();
    shader.projection = projection;
    shader.shadowM = shadow_mvp * (view * model2world).invert();//transform from view space to shadow clip space
    shader.texture = &texture;
    shader.normalmap = &normalmap;
    shader.specularmap = &specularmap;
    shader.shadowmap = &shadowmap;
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
    img.write_tga_file("shadow.tga");

    shadowmap.flip_vertically();
    shadowmap.write_tga_file("shadowmap.tga");

    return 0;
}