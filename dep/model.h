#ifndef __MODEL_H__
#define __MODEL_H__

#include <vector>
#include <tuple>
#include "geometry.h"

class Model {
private:
	std::vector<Vec3f> verts_;
	std::vector<Vec3f> texcoords_;
	std::vector<Vec3f> norm_;
	std::vector<
		//[(vert, tex, norm), (vert, tex, norm), (vert, tex, norm)]
		std::vector<std::tuple<
			int/*vert idx*/, int /*tex idx*/, int /*normal idx*/>>
	> faces_;
public:
	Model(const char *filename);
	~Model();
	int nverts() const;
	int nfaces() const;
	Vec3f vert(int i) const;
	Vec3f tex(int i) const;
	Vec3f norm(int i) const;
	std::vector<std::tuple<int/*vert idx*/, int /*tex idx*/, int /*normal idx*/>> face(int idx) const;
};

#endif //__MODEL_H__
