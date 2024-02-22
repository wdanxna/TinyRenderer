#ifndef __MODEL_H__
#define __MODEL_H__

#include <vector>
#include <tuple>
#include "geometry.h"

class Model {
private:
	std::vector<Vec3f> verts_;
	std::vector<Vec3f> texcoords_;
	std::vector<
		//[(vert, tex), (vert, tex), (vert, tex)]
		std::vector<std::tuple<int/*vert idx*/, int /*tex idx*/>>
	> faces_;
public:
	Model(const char *filename);
	~Model();
	int nverts();
	int nfaces();
	Vec3f vert(int i);
	Vec3f tex(int i);
	std::vector<std::tuple<int/*vert idx*/, int /*tex idx*/>> face(int idx);
};

#endif //__MODEL_H__
