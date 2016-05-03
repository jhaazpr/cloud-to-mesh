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
    class BPAEdge;
    class BPALoop;
    class BPAFront;

    class BPAEdge {
    public:
      BPAEdge();
      BPAEdge( Index i, Index j, Index o, BPAEdge *prev_edge, BPAEdge *next_edge,
               BPALoop *my_loop );

      // Either ACTIVE or BOUNDARY (not active)
      bool is_active;
      Index i;
      Index j;
      Index o;
      BPAEdge* prev_edge;
      BPAEdge* next_edge;
      BPALoop* my_loop;


      /**
       * If other_edge utilizes the same points, as this edge does,
       * resolve the conflict by removing one of the edges.
       */
      void glue(BPAEdge *other_edge, BPAFront *front);

      /**
       *  Rotates the ball by considering all points within a 2*rho radius
       *  and checking whether the center lies on the circle gamma. Returns
       *  True if a vertex k was found.
       */
      bool ball_pivot(double rho, Index *k);

      /**
       * Marks an edge as boundary if already fully explored.
       */
      void mark_not_active(void);

    };

    class BPALoop {
    public:
      /**
       * A BPALoop contains a pointer to the first of the linked list of edges,
       * and to the BPAFront that owns it.
       */
      BPALoop( BPAEdge *start_edge, BPAFront *my_front );

      BPAEdge *start_edge;
      BPAFront *my_front;
    };

    class BPAFront {
    public:
      BPAFront( std::vector<Vector3D> vertices, Polymesh* pm);
      void BPA(double rho);
      /**
       * Core fields
       */
      std::vector<BPALoop> loops;
      std::vector<Vector3D> vertices;
      /**
       * Takes in the index of a vertex and inserts edges ik and kj
       * by removing the edge ij from the loop the edge belongs too.
       */
      void join(BPAEdge* e_ij, Index k);

     /**
      * Adds the triangle i,j,k to the polymesh
      */
      void output_triangle(Index i, Index j, Index k);
      /**
       * Keep track of our own lil polymesh that we will of course
       * keep updated the entire time.
       */
       Polymesh *pm;

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
      bool get_active_edge(BPAEdge * e);

      /**
       * Add an edge as a new loop in the front. Don't forget to update stuff.
       * Returns the loop in which the method inserted the edge.
       */
      BPALoop *insert_edge(BPAEdge *edge);

      /**
       * Grab three vertices from a seed triangle that the ball rolls onto.
       */
      bool find_seed_triangle(std::vector<Index> *indices, double rho);
    private:
      bool find_seed_trangle_indices(std::vector<Index> * indices, double rho);
    };


}

#endif // STUDENT_CODE_H
