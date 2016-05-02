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
    Polymesh BPA(std::vector<Vector3D>& vertices);
    class BPAEdge;
    class BPALoop;
    class BPAFront;

    class BPAEdge {
    public:
      BPAEdge();
      BPAEdge( Index i, Index j, BPAEdge& prev_edge, BPAEdge& next_edge,
               BPALoop& my_loop );

      // Either ACTIVE or BOUNDARY (not active)
      bool is_active;
      Index i;
      Index j;
      BPAEdge& prev_edge;
      BPAEdge& next_edge;
      BPALoop& my_loop;

      /**
       * To do
       */
      void join(Index index);

      /**
       * To do
       */
      void glue(BPAEdge other_edge, BPAFront front);

      /**
       * To do
       */
      Vector3D ball_pivot(double rho);

      /**
       * To do
       */
      void mark_boundary(void);

    };

    class BPALoop {
    public:
      /**
       * A BPALoop contains a pointer to the first of the linked list of edges,
       * and to the BPAFront that owns it.
       */
      BPALoop( BPAEdge& start_edge, BPAFront& my_front );

      BPAEdge& start_edge;
      BPAFront& my_front;
    };

    class BPAFront {
    public:
      BPAFront( std::vector<Vector3D> vertices , double rho );

      /**
       * Core fields
       */
      std::vector<BPALoop> loops;
      std::vector<Vector3D> vertices;
      double rho;

      /**
       * Keep track of our own lil polymesh that we will of course
       * keep updated the entire time.
       */
       Polymesh& pm;

      /**
       * For each vertex index, and array of booleans to convey the state
       * of the vertex. These need to be updated upon:
       * --- insterting triangles into the mesh
       * --- join and glue operations
       */
      std::vector<bool> vertices_on_front;
      std::vector<bool> vertices_used;

      /**
       * Pulls any active edge from the front.
       */
      BPAEdge& get_active_edge(void);

      /**
       * Add an edge as a new loop in the front. Don't forget to update stuff.
       */
      void insert_edge(BPAEdge& edge);

      /**
       * Grab three vertices from a seed triangle that the ball rolls onto.
       */
      std::vector<Vector3D> find_seed_triangle(void);
    };


}

#endif // STUDENT_CODE_H
