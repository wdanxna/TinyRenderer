#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include "model.h"

Model::Model(const char *filename) : verts_(), faces_() {
    std::ifstream in;
    in.open (filename, std::ifstream::in);
    if (in.fail()) return;
    std::string line;
    while (!in.eof()) {
        std::getline(in, line);
        std::istringstream iss(line.c_str());
        char trash;
        if (!line.compare(0, 2, "v ")) {
            iss >> trash;
            Vec3f v;
            for (int i=0;i<3;i++) iss >> v.raw[i];
            verts_.push_back(v);
        } else if (!line.compare(0, 2, "vt")) {
            //parse texture coordinates
            iss >> trash >> trash;
            Vec3f uv;
            for (int i = 0; i < 3; i++) iss >> uv.raw[i];
            texcoords_.push_back(uv);
        } else if (!line.compare(0, 2, "f ")) {
            std::vector<std::tuple<int/*vert idx*/, int /*tex idx*/>> f;
            int vert_idx, tex_idx, itrash;
            iss >> trash;
            while (iss >> vert_idx >> trash >> tex_idx >> trash >> itrash) {
                vert_idx--; // in wavefront obj all indices start at 1, not zero
                tex_idx--;
                f.push_back({vert_idx, tex_idx});
            }
            faces_.push_back(f);
        }
    }
    std::cerr << "# v# " << verts_.size() << " f# "  << faces_.size() << std::endl;
}

Model::~Model() {
}

int Model::nverts() {
    return (int)verts_.size();
}

int Model::nfaces() {
    return (int)faces_.size();
}

std::vector<std::tuple<int/*vert idx*/, int /*tex idx*/>> Model::face(int idx) {
    return faces_[idx];
}

Vec3f Model::vert(int i) {
    return verts_[i];
}

Vec3f Model::tex(int i) {
    return texcoords_[i];
}

