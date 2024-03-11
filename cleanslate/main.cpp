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
    if (i.z < 1e-2) return {-1,1,1};
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
    const int width = out.get_width();
    const int height = out.get_height();
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
            float frag_depth = (clip_verts * clip_bary)[2];

            if (zbuffer[x + y * width] < frag_depth) {
                zbuffer[x + y * width] = frag_depth;
                TGAColor color;
                shader.fragment(clip_bary, color);
                out.set(x, y, color);
            }
        }
    }
}

void rasterize(
    Shader& shader,
    const Matrix &vp, 
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
    mat<4,4,float> u_norm;
    Vec3f u_lightDir;

    //varyings
    mat<3,3,float> varying_normals;

    Gouraud(mat<4,4,float> mvp, mat<4,4,float> norm_m, Vec3f lightDir) 
    : u_mvp{mvp}, u_norm{norm_m}, u_lightDir{lightDir} {};

    virtual Vec4f vertex(int iface, int nthvert) override {
        const auto& face = _model->face(iface);
        Vec3f obj_v = _model->vert(face[nthvert]);
        Vec3f obj_n = _model->normal(iface, nthvert);
        varying_normals.set_col(nthvert, proj<3>(u_norm * embed<4>(obj_n)));

        return u_mvp * embed<4>(obj_v);
    }

    virtual void fragment(Vec3f clip_bary, TGAColor& color) override {
        Vec3f n = (varying_normals * clip_bary).normalize();
        float intensity = std::max(0.0f, n*(u_lightDir*-1.0f).normalize());
        color = TGAColor(255 * intensity, 255 * intensity, 255 * intensity, 255);
    }
};

class Phong : public Shader {
public:
    //uniforms
    mat<4,4,float> u_mvp;
    mat<4,4,float> u_nm;//normal matrix
    mat<4,4,float> u_sm;//shadow matrix, transform from clip space to shadow space
    Vec3f u_lightPos, u_lightLookat;

    float* u_shadowmap;
    TGAImage* u_aomap;

    //varyings
    mat<3,3, float> varying_normals;
    mat<2,3, float> varying_uvs;
    mat<3,3,float> varying_pos;

    Phong(mat<4,4,float> mvp, mat<4,4,float> nm, mat<4,4,float> sm, Vec3f lightPos, Vec3f lightLookat)
    : u_mvp{mvp}, u_nm{nm}, u_sm{sm}, u_lightPos{lightPos}, u_lightLookat{lightLookat} {}

    virtual Vec4f vertex(int iface, int nthvert) override {
        auto face = _model->face(iface);
        auto obj_v = _model->vert(face[nthvert]);
        auto obj_n = _model->normal(iface, nthvert);
        auto uv = _model->uv(iface, nthvert);
        auto clip_v = u_mvp * embed<4>(obj_v);
        varying_normals.set_col(nthvert, proj<3>(u_nm * embed<4>(obj_n)));
        varying_uvs.set_col(nthvert, uv);
        varying_pos[nthvert] = proj<3>(clip_v);
        return clip_v;
    }

    virtual void fragment(Vec3f clip_bary, TGAColor& color) override {
        Vec2f uv = varying_uvs * clip_bary;
        Vec3f n = (varying_normals * clip_bary).normalize();
        Vec3f l = proj<3>(u_mvp * embed<4>(u_lightPos) - u_mvp * embed<4>(u_lightLookat)).normalize();
        Vec3f v = varying_pos.transpose() * clip_bary;
        
        //calculate tangent space normal
        mat<3,3,float> A;
        A[0] = varying_pos[1] - varying_pos[0];
        A[1] = varying_pos[2] - varying_pos[0];
        A[2] = n;
        
        auto AI = A.invert();
        Vec3f i = AI * Vec3f{varying_uvs[0][1] - varying_uvs[0][0],
                              varying_uvs[0][2] - varying_uvs[0][0],
                              0.0f};
        Vec3f j = AI * Vec3f{varying_uvs[1][1] - varying_uvs[1][0],
                              varying_uvs[1][2] - varying_uvs[1][0],
                              0.0f};

        mat<3,3,float> tangent;
        tangent.set_col(0, i.normalize());
        tangent.set_col(1, j.normalize());
        tangent.set_col(2, n);

        n = (tangent * _model->normal(uv));

        auto sv = u_sm * embed<4>(v);
        sv = sv / sv[3];
        auto sx = int((sv[0]+1.0f)/2.0f * width + 0.5f);
        auto sy = int((sv[1]+1.0f)/2.0f * height + 0.5f);
        float shadow = 1.0f;
        if (sx >= 0 && sx < width && sy >= 0 && sy < height) {
            float shadow_val = u_shadowmap[sx + sy * width];
            float shadow_z = sv[2];
            shadow = 0.3f + 0.7f * (shadow_val < shadow_z+0.02);
        }
        // shadow = 1.0f;
        
        //ambient term
        float ambient = 0.5f;
        float occlu = 1.0f;
        if (u_aomap) {
            occlu = u_aomap->get(uv.x*u_aomap->get_width(), uv.y*u_aomap->get_height())[0] / 255.0f;
            ambient *= occlu;
        }
        //diffuse term
        TGAColor c = _model->diffuse(uv);
        float diffuse = std::max(0.0f, n*l) * 1.8f;
        //specular term
        Vec3f r = (n*(2.0f*(l*n))-l).normalize();
        float spec = std::powf(std::max(0.0f, r.z), std::max(1.0f, _model->specular(uv))) * 1.0f;

        color = TGAColor(
            std::min(255, int(c[2]*ambient) + int(c[2]*diffuse*shadow) + int(c[2]*spec*shadow)),
            std::min(255, int(c[1]*ambient) + int(c[1]*diffuse*shadow) + int(c[1]*spec*shadow)),
            std::min(255, int(c[0]*ambient) + int(c[0]*diffuse*shadow) + int(c[0]*spec*shadow)),
            255
        );
    }

};

class Shadow : public Shader {
public:
    mat<4,4,float> u_mvp;

    mat<4,3,float> varying_pos;
    Shadow(mat<4,4,float> mvp) : u_mvp{mvp} {}

    virtual Vec4f vertex(int iface, int nthvert) override {
        auto clip_v = u_mvp * embed<4>(_model->vert(iface, nthvert));
        varying_pos.set_col(nthvert, clip_v);
        return clip_v;
    }

    virtual void fragment(Vec3f clip_bary, TGAColor& color) override {
        auto v = varying_pos * clip_bary;

        color = TGAColor(255,255,255,255) * ((v[2]/v[3] + 3.0f)/6.0f);
        color[3] = 255;
    }
};

class AoImage : public Shader {
public:
    mat<4,4,float> u_mvp;
    float* u_zbuffer;
    TGAImage* u_occlu;

    mat<4,3,float> varying_pos;
    mat<2,3, float> varying_uvs;

    AoImage(mat<4,4,float> mvp, float* zbuffer, TGAImage* occlu) 
        : u_mvp{mvp}, u_zbuffer{zbuffer}, u_occlu{occlu} {}

    virtual Vec4f vertex(int iface, int nthvert) override {
        auto clip_v = u_mvp * embed<4>(_model->vert(iface, nthvert));
        varying_pos.set_col(nthvert, clip_v);
        varying_uvs.set_col(nthvert, _model->uv(iface, nthvert));
        return clip_v;
    }

    virtual void fragment(Vec3f clip_bary, TGAColor& color) override {
        auto v = varying_pos * clip_bary;
        auto uv = varying_uvs * clip_bary;
        
        auto sx = int((v[0]/v[3]+1.0f)/2.0f * width + 0.5f);
        auto sy = int((v[1]/v[3]+1.0f)/2.0f * height + 0.5f);
        if (sx >= 0 && sx < width && sy >= 0 && sy < height) {
            float z = u_zbuffer[sx + sy * width];
            if ((v[2]/v[3]) - z >= 0.0f) {
                u_occlu->set(
                    uv.x*u_occlu->get_width(), 
                    uv.y*u_occlu->get_height(), TGAColor(255));
            }
        }
        // color = TGAColor(255,255,255,255) * ((v[2]/v[3] + 3.0f)/6.0f);
        // color[3] = 255;
    }
};

Vec3f rand_point_on_unit_sphere() {
    float u = (float)rand()/(float)RAND_MAX;
    float v = (float)rand()/(float)RAND_MAX;
    float theta = 2.f*M_PI*u;
    float phi   = acos(2.f*v - 1.f);
    return Vec3f(sin(phi)*cos(theta), sin(phi)*sin(theta), cos(phi));
}

void make_aoimage(Model& model, TGAImage& aoimage) {
    const int ao_width = aoimage.get_width();
    const int ao_height = aoimage.get_height();
    const int samples = 1000;
    float* zbuff = new float[ao_width*ao_height];
    float* shadow = new float[ao_width*ao_height];
    TGAImage tmp(width, height, TGAImage::RGB);
    TGAImage ao_tmp(ao_width, ao_height, TGAImage::RGB);
    TGAImage ao_fb(ao_width, ao_height, TGAImage::RGB);
    srand(time(0));
    for (int iter = 1; iter <= samples; iter++) {
        tmp.clear();
        auto p = rand_point_on_unit_sphere();
        p.y = std::abs(p.y);
        auto view = lookAt(p, {0,0,0.01}, {0.0, 1.0, 0.0});
        auto vp = viewport(width, height);
        std::fill(zbuff, zbuff + (width * height), std::numeric_limits<float>::lowest());
        std::fill(zbuff, shadow + (width * height), std::numeric_limits<float>::lowest());
        //first pass, get z buffer
        Shadow zshader(view);
        zshader._model = &model;
        rasterize(zshader, vp, shadow, tmp);
        //visual debug
        // tmp.flip_vertically();
        // tmp.write_tga_file("tmp.tga");

        //second pass, generate aoimage
        vp = viewport(ao_width, ao_height);
        AoImage aoshader(view, shadow, &ao_tmp);
        aoshader._model = &model;
        std::fill(zbuff, zbuff + (ao_width * ao_height), std::numeric_limits<float>::lowest());
        rasterize(aoshader, vp, zbuff, ao_fb);

        for (int i = 0; i < ao_width; i++) {
            for (int j = 0; j < ao_height; j++) {
                float prev = aoimage.get(i, j)[0];
                float cur = ao_tmp.get(i, j)[0];
                float t = (prev*(iter-1)+cur*0.5)/(float)iter;
                TGAColor c(t,t,t);
                aoimage.set(i,j, c);
            }
        }
    }

    delete[] zbuff;
    delete[] shadow;
}

float max_elevation_angle(float *zbuffer, Vec2f p, Vec2f dir) {
    float maxangle = 0;
    for (float step = 0.0f; step < 100; step+=1.0f) {
        Vec2f sample = p + dir*step;
        if (sample.x >= width || sample.y >= height || sample.x < 0 || sample.y < 0) break;
        float dist = (sample-p).norm();
        if (dist < 1.0) continue;
        float elevation = zbuffer[int(sample.x) + int(sample.y)*width] - zbuffer[int(p.x) + int(p.y)*width];
        if (elevation > 0.22) elevation = 0;
        maxangle = std::max(maxangle, atanf(elevation / dist));
    }
    return maxangle;
}

int main() {

    TGAImage framebuffer(width, height, TGAImage::RGBA);
    TGAImage zimg(width, height, TGAImage::RGBA);

    mat<4,4,float> modelworld = mat<4,4,float>::identity();
    modelworld[0][0] = 0.9f;
    modelworld[1][1] = 0.9f;
    modelworld[2][2] = 0.9f;
    modelworld[0][3] = 0.2;//tx
    modelworld[1][3] = 0.0;//ty
    modelworld[2][3] = 0.2f;//tz

    mat<4,4,float> floor_modelworld = mat<4,4,float>::identity();
    floor_modelworld[0][3] = 0;//tx
    floor_modelworld[1][3] = 0.15;//ty
    floor_modelworld[2][3] = -0.5f;//tz

    Vec3f eye {0.07, 0.1f, 0.25f};
    // Vec3f eye {0.5, 0.0f, 0.0f};
    mat<4,4,float> view = lookAt(eye, {0.0f, 0.0f, 0.001f}, {0.0f, 1.0f, 0.0f});
    mat<4,4,float> vp = viewport(width, height);
    mat<4,4,float> projection = perspective(1.5f);

    float zbuffer[width*height];
    std::fill(zbuffer, zbuffer + (width * height), std::numeric_limits<float>::lowest());

    Model head("../res/diablo3_pose.obj");
    Model floor("../res/floor.obj");

    // TGAImage *aoimage = new TGAImage(1024, 1024, TGAImage::RGB);
    // make_aoimage(head, *aoimage);
    // aoimage->flip_vertically();
    // aoimage->write_tga_file("ao2.tga");

    //shadow pass
    TGAImage shadowimg(width, height, TGAImage::RGBA);
    float shadowmap[width*height];
    std::fill(shadowmap, shadowmap + (width * height), std::numeric_limits<float>::lowest());
    Vec3f light_pos {0.5f, 0.5f, 0.5f};
    Vec3f light_lookat {0.0f, 0.0f, 0.0f};
    mat<4,4,float> light_view = lookAt(light_pos, light_lookat, {0.0f, 1.0f, 0.0f});
    Shadow shadow(light_view * modelworld);
    shadow._model = &head;
    rasterize(shadow, vp, shadowmap, shadowimg);

    shadow.u_mvp = light_view * floor_modelworld;
    shadow._model = &floor;
    rasterize(shadow, vp, shadowmap, shadowimg);
    shadowimg.flip_vertically();
    shadowimg.write_tga_file("shadowmap.tga");


    //SSAO
    TGAImage ssao_frame(width, height, TGAImage::RGB);
    {
        //depth pass
        Shadow zshader(projection*view*modelworld);
        zshader._model = &head;
        rasterize(zshader, vp, zbuffer, shadowimg);
        zshader.u_mvp = projection*view*floor_modelworld;
        zshader._model = &floor;
        rasterize(zshader, vp, zbuffer, shadowimg);

        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {
                if (zbuffer[x + y * width] < -1e5) continue;
                float total = 0;
                for (float a = 0; a < M_PI*2 - 1e-4; a += M_PI_4) {
                    //8 direction
                    float u = M_PI*2 + (float)rand()/(float)RAND_MAX;
                    total += M_PI_2 - max_elevation_angle(zbuffer, {(float)x, (float)y}, {cos(a+u), sin(a+u)});
                }
                total /= M_PI_2 * 8;
                total = powf(total, 100.f);
                ssao_frame.set(x, y, TGAColor(total * 255, total * 255, total * 255));
            }
        }
        ssao_frame.flip_vertically();
        ssao_frame.write_tga_file("ssao.tga");
    }

    std::fill(zbuffer, zbuffer + (width * height), std::numeric_limits<float>::lowest());
    //light pass
    Phong shader(
        projection*view*modelworld, 
        (view*modelworld).invert_transpose(), 
        light_view * modelworld * (view*modelworld).invert(), 
        light_pos,
        light_lookat);

    shader.u_shadowmap = shadowmap;
    shader._model = &head;
    TGAImage *ao = new TGAImage();
    ao->read_tga_file("./ssao.tga");
    ao->flip_vertically();
    shader.u_aomap = ao;
    rasterize(shader, vp, zbuffer, framebuffer);
    delete ao;

    shader.u_mvp = projection*view*floor_modelworld;
    shader.u_nm = (view*floor_modelworld).invert_transpose();
    shader.u_sm = light_view * floor_modelworld * (view * floor_modelworld).invert();
    shader._model = &floor;
    shader.u_aomap = nullptr;
    rasterize(shader, vp, zbuffer, framebuffer);

    //zbuffer visualization
    for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            zimg.set(x, y, TGAColor(255,255,255,255) * ((zbuffer[x + y*width] + 3.0f)/6.0f));
        }
    }
    zimg.flip_vertically();
    zimg.write_tga_file("z.tga");


    framebuffer.flip_vertically();//make (0,0) at bottom left, x going right, y going up
    framebuffer.write_tga_file("frame.tga");
    return 0;
}