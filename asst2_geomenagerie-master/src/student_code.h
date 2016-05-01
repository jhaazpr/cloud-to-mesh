#ifndef STUDENT_CODE_H
#define STUDENT_CODE_H

#include "halfEdgeMesh.h"
#include "bezierPatch.h"

using namespace std;

namespace CGL {

    class MeshResampler {
    public:

        MeshResampler() {
        };

        ~MeshResampler() {
        }

        void upsample(HalfedgeMesh& mesh);
        void simplify(HalfedgeMesh& mesh);
    };

    // Mesh Construction
    void BPA(std::vector<Vector3D>& vertices);


}

#endif // STUDENT_CODE_H
