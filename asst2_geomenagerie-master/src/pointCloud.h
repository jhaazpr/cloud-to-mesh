#ifndef CGL_POINT_CLOUD_H
#define CGL_POINT_CLOUD_H

#include "scene.h"
#include "material.h"

namespace CGL {

    /**
     * A PointCloud is just a storage element for vertices and potentially normals,
     * it doesn't have any polygon rendering abilites. These should be transformed
     * and refactored as meshes if you want to do anything besides rendering points.
     */
    struct PointCloud : Instance {
        std::string id;
        std::string name;

        std::vector<Vector3D> vertices; ///< polygon vertex array
        std::vector<Vector3D> normals; ///< polygon normal array

    }; // struct Polymesh

    // std::ostream& operator<<(std::ostream& os, const Polymesh& polymesh);

} // namespace CGL

#endif // CGL_POINT_CLOUD_H
