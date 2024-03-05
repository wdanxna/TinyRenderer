#include <dep/tgaimage.h>
#include <dep/model.h>
#include "main.h"

const int width = 500;
const int height = 500;


mat<4,4,float> viewport(int width, int height) {
    //viewport = scale * T
    mat<4,4,float> m;
    m[0] = {width * 0.5f, 0.f, 0.f,  width * 0.5f};
    m[1] = {0.f, height * 0.5f, 0.f,  height * 0.5f};
    m[2] = {0.f, 0.f,  0.5f, 0.5f};
    m[3] = {0.f, 0.f,  0.f,  1.f};
    return m;
}

mat<4,4,float> perspective(float c) {
    mat<4, 4, float> m;
    m[0] = {1.f,  0.0f, 0.0f, 0.0f};
    m[1] = {0.0f, 1.f,  0.0f, 0.0f};
    m[2] = {0.0f, 0.0f, 1.0f, 0.0f};
    m[3] = {0.0f, 0.0f, -1.0f/c, 1.f};
    return m;
}

mat<4,4,float> lookAt(const Vec3f& eye, const Vec3f& target, const Vec3f& up) {
    mat<4,4,float> t = mat<4,4,float>::identity();
    t[0][3] = -eye.x;
    t[1][3] = -eye.y;
    t[2][3] = -eye.z;

    Vec3f z = (eye - target).normalize();
    Vec3f x = cross(up, z).normalize();
    Vec3f y = cross(z, x).normalize();

    mat<4,4,float> rot;
    rot.set_col(0, embed<4>(x, 0.0f));
    rot.set_col(1, embed<4>(y, 0.0f));
    rot.set_col(2, embed<4>(z, 0.0f));
    rot.set_col(3, {0.0f, 0.0f, 0.0f, 1.0f});

    return rot.invert() * t;
}


Vec3f barycentric(mat<3,2,float> v, Vec2f p) {
    Vec2f p0p1 = v[1] - v[0];
    Vec2f p0p2 = v[2] - v[0];
    Vec2f pp0 = v[0] - p;

    auto i = cross(Vec3f(p0p1.x, p0p2.x, pp0.x), 
                   Vec3f(p0p1.y, p0p2.y, pp0.y));
    i = i / i.z;
    return Vec3f(1.0-i.x-i.y, i.x, i.y);
}

class Shader {
public:
    Model* _model;
    virtual Vec4f vertex(int iface, int nthvert) = 0;
    virtual void fragment(Vec3f clip_bary, TGAColor& color) = 0;
};

void triangle(
    Shader& shader,
    mat<4,3,float> clip_verts,
    mat<4,4,float> viewport_matrix,
    float* zbuffer, TGAImage& out) {
    mat<3,4,float> cv = clip_verts.transpose();
    mat<3,2,float> sv;
    for (int i = 0; i < 3; i++) {
        sv[i] = proj<2>(viewport_matrix * (cv[i]/cv[i][3]));
    }

    int left = INT_MAX;
    int right = 0;
    int bottom = INT_MAX;
    int up = 0;
    for (int i = 0; i < 3; i++) {
        left = std::min(width-1, std::max(0, std::min(left, int(sv[i].x + 0.5))));
        right = std::min(width-1, std::max(0, std::max(right, int(sv[i].x + 0.5))));
        up = std::min(height-1, std::max(0, std::max(up, int(sv[i].y + 0.5))));
        bottom = std::min(height-1, std::max(0, std::min(bottom, int(sv[i].y + 0.5))));
    }
    for (int x = left; x <= right; x++) {
        for (int y = bottom; y <= up; y++) {
            auto screen_bary = barycentric(sv, {(float)x, (float)y});
            if (screen_bary.x < 0 || screen_bary.y < 0 || screen_bary.z < 0) continue;
            //derive clip space barycentric from screen space barycentric
            auto clip_bary = Vec3f{ screen_bary.x/cv[0][3], screen_bary.y/cv[1][3], screen_bary.z/cv[2][3]};
            clip_bary = clip_bary / (clip_bary[0] + clip_bary[1] + clip_bary[2]);

            TGAColor color;
            shader.fragment(clip_bary, color);

            //shift z to prevent cutoff
            float frag_depth = (clip_verts * clip_bary)[2] / 4.0f + 1.0f;

            if (zbuffer[x + y * width] < frag_depth) {
                zbuffer[x + y * width] = frag_depth;
                out.set(x, y, color);
            }
        }
    }
}

void rasterize(
    Shader& shader,
    Matrix &projection, 
    Matrix &view, 
    Matrix &modelworld, 
    Matrix &vp, 
    float* zbuffer, 
    TGAImage &framebuffer)
{
    for (int i = 0; i < shader._model->nfaces(); i++)
    {
        // auto vert_ids = shader._model->face(i);
        // Vec3f object_verts[3];
        mat<3,4,float> clip_verts;
        Vec3f ndc_verts[3];
        Vec3f screen_verts[3];
        
        for (int j = 0; j < 3; j++)
        {
            clip_verts[j] = shader.vertex(i, j);
            ndc_verts[j] = proj<3>(clip_verts[j] / clip_verts[j][3]);
            screen_verts[j] = proj<3>(vp * embed<4>(ndc_verts[j]));
        }

        triangle(shader, clip_verts.transpose(), vp, zbuffer, framebuffer);
    }
}


class Gouraud : public Shader {
public:
    //uniforms
    mat<4,4,float> u_mvp;
    Vec3f u_lightDir;

    //varyings
    mat<3,3,float> varying_normals;

    Gouraud(mat<4,4,float> mvp, Vec3f lightDir) 
    : u_mvp{mvp}, u_lightDir{lightDir} {};

    virtual Vec4f vertex(int iface, int nthvert) override {
        const auto& face = _model->face(iface);
        Vec3f obj_v = _model->vert(face[nthvert]);
        Vec3f obj_n = _model->normal(iface, nthvert);
        varying_normals.set_col(nthvert, obj_n);

        return u_mvp * embed<4>(obj_v);
    }

    virtual void fragment(Vec3f clip_bary, TGAColor& color) override {
        Vec3f n = (varying_normals * clip_bary).normalize();
        float intensity = std::max(0.0f, n*(u_lightDir*-1.0f).normalize());
        color = TGAColor(255 * intensity, 255 * intensity, 255 * intensity, 255);
    }
};

int main() {

    TGAImage framebuffer(width, height, TGAImage::RGBA);
    TGAImage zimg(width, height, TGAImage::RGBA);

    mat<4,4,float> modelworld = mat<4,4,float>::identity();
    modelworld[0][3] = 0;//tx
    modelworld[1][3] = 0.1;//ty
    modelworld[2][3] = -0.3f;//tz

    Vec3f eye {0.2f/1.0f, 0.1f/1.0f, 0.5/1.0f};
    mat<4,4,float> view = lookAt(eye, {0.0f, 0.0f, 0.001f}, {0.0f, 1.0f, 0.0f});

    mat<4,4,float> vp = viewport(width, height);
    mat<4,4,float> projection = perspective(1.0f);

    float zbuffer[width*height];
    for (int i = 0; i < width*height;i++) {
        zbuffer[i] = std::numeric_limits<float>::min();
    }

    Model head("../res/african_head.obj");
    Gouraud gouraud(projection*view*modelworld, {0.2f, -0.3f, -1.0f});
    gouraud._model = &head;
    rasterize(gouraud, projection, view, modelworld, vp, zbuffer, framebuffer);

    Model floor("../res/floor.obj");
    gouraud._model = &floor;
    rasterize(gouraud, projection, view, modelworld, vp, zbuffer, framebuffer);

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