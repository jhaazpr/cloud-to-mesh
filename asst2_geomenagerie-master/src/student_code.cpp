/*
 * Student solution for UC Berkeley Project 2
 *
 * Implemented by ____ on ____.
 *
 */

#include "student_code.h"
#include "mutablePriorityQueue.h"

namespace CGL {
    // In case you want a pyramid NOW.
    Polymesh BPA_hardcode(std::vector<Vector3D>& vertices);

    /**
     * Because of the interface in meshEdit that already exists, we need to return
     * a Polymesh which then gets converted into a HalfedgeMesh. But, it's difficult
     * to construct a Polymesh if we're doing a lot of edge operations. So, we become
     * kludgy and create a HalfedgeMesh, using all the methods and trappings of a
     * HalfedgeMesh, then rip the HalfedgeMesh's vertices and Polygons, put them into
     * a Polymesh, which then gets refactored into a HalfedgeMesh anyway in the
     * init_polymesh() function in meshEdit.cpp
     */
    Polymesh BPA(std::vector<Vector3D>& vertices) {
      cout << "BPA yo, voyteces: " << vertices.size() << endl;

      Polymesh pm = BPA_hardcode(vertices);
      return pm;
    }

    //FIXME: issues initializing edges
    // BPAEdge::BPAEdge( void )
    //   : i(0), j(0), prev_edge(nullptr), next_edge(nullptr), my_loop(nullptr)
    // {
    //   i = 0;
    //   j = 0;
    //   prev_edge = this;
    //   next_edge = this;
    //   my_loop = this;
    // }

    BPAEdge::BPAEdge( Index i, Index j, BPAEdge& prev_edge, BPAEdge& next_edge,
             BPALoop& my_loop )
             : i(i), j(j), prev_edge(prev_edge), next_edge(next_edge), my_loop(my_loop)
    {
    }

    /**
     * Documentation goes here
     */
    void BPAEdge::join(Index index) {
      //TODO
    }

    /**
     * Documentation goes here
     */
    void BPAEdge::glue(BPAEdge other_edge, BPAFront front) {
      //TODO
    }

    /**
     * Documentation goes here
     */
    Vector3D BPAEdge::ball_pivot(double rho) {
      //TODO
      Vector3D new_vertex;
      return new_vertex;
    }

    /**
     * Documentation goes here
     */
    void BPAEdge::mark_boundary(void) {
      //TODO
    }

    BPALoop::BPALoop( BPAEdge& start_edge, BPAFront& my_front )
      : start_edge(start_edge), my_front(my_front)
    {
    }

    //FIXME: needs to initalize polymesh
    // BPAFront::BPAFront( std::vector<Vector3D> vertices , double rho )
    //   : vertices(vertices), rho(rho)
    // {
    // }

    /**
     * Pulls any active edge from the front.
     */
     //FIXME: need to actually initialize edge
    // BPAEdge& BPAFront::get_active_edge(void) {
    //   BPAEdge& edge();
    //   return edge;
    // }

    /**
     * Add an edge as a new loop in the front. Don't forget to update stuff.
     */
    void BPAFront::insert_edge(BPAEdge& edge) {
      //TODO
    }

    /**
     * Grab three vertices from a seed triangle that the ball rolls onto.
     */
    std::vector<Vector3D> BPAFront::find_seed_triangle(void) {
      //TODO
      std::vector<Vector3D> triangle;
      return triangle;
    }

    ////////////////////////////////////////////////////////////////////

    void BezierPatch::preprocess() {
        // TODO Part 1.
        // TODO If you use the matrix form for Bezier patch evaluation, you will need to
        // TODO compute your matrices based on the 16 control points here.
        // TODO You will also need to define your matrices
        // TODO as member variables in the "BezierPatch" class.
        // TODO If you use De Casteljau's recursive algorithm, you will not need to do anything here.
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                Cx[i][j] = controlPoints[i][j].x;
                Cy[i][j] = controlPoints[i][j].y;
                Cz[i][j] = controlPoints[i][j].z;
            }
        }

    }

    Vector3D BezierPatch::evaluate(double u, double v) const {
        // TODO Part 1.
        // TODO Returns the 3D point whose parametric coordinates are (u, v) on the Bezier patch.
        // TODO Note that both u and v are within [0, 1].
        Vector3D point;
        double datau[] = {pow(u,3), pow(u,2), u, 1,0,0,0,0,0,0,0,0,0,0,0,0};
        Matrix4x4 U = Matrix4x4(datau);
        double datav[] = {pow(v,3), pow(v,2), v, 1,0,0,0,0,0,0,0,0,0,0,0,0};
        Matrix4x4 V = Matrix4x4(datav).T();
        double datam[] = {-1, 3, -3, 1,3,-6,3,0,-3,3,0,0,1,0,0,0};
        Matrix4x4 M = Matrix4x4(datam);
        point.x = (U*M*Cx*M.T()*V)[0][0];
        point.y = (U*M*Cy*M.T()*V)[0][0];
        point.z = (U*M*Cz*M.T()*V)[0][0];

        return point;
    }

    void BezierPatch::add2mesh(Polymesh* mesh) const {
        // TODO Part 1.
        // TODO Tessellate the given Bezier patch into triangles uniformly on a 8x8 grid(8x8x2=128 triangles) in parameter space.
        // TODO You will call your own evaluate function here to compute vertex positions of the tessellated triangles.
        // TODO The "addTriangle" function inherited from the "BezierPatchLoader" class may help you add triangles to the output mesh.
        double t = 1.0/8;
        for (double i=0; i<1; i+=t) {
            for(double j=0; j<1; j+=t) {
                Vector3D p1 = evaluate(i,j);
                Vector3D p2 = evaluate(i+t,j);
                Vector3D p3 = evaluate(i,j+t);
                Vector3D p4 = evaluate(i+t,j+t);
                addTriangle(mesh,p1,p2,p3);
                addTriangle(mesh,p3,p2,p4);
            }
        }


    }

    Vector3D Vertex::normal(void) const
    // TODO Part 2.
    // TODO Returns an approximate unit normal at this vertex, computed by
    // TODO taking the area-weighted average of the normals of neighboring
    // TODO triangles, then normalizing.
    {
        // TODO Compute and return the area-weighted unit normal.
        Vector3D n(0,0,0);
        HalfedgeCIter h = halfedge();
        Vector3D p0 = halfedge()->vertex()->position;
          do {
            HalfedgeCIter h_twin = h->twin();
            Vector3D p1 = h_twin->vertex()->position;
            Vector3D edge1 = p0-p1;

            h = h_twin->next();
            Vector3D p2 = h->twin()->vertex()->position;
            Vector3D edge2 = p2-p0;

            n = n+cross(edge1,edge2);
          } while(h != halfedge());
        return n.unit();
    }

    EdgeIter HalfedgeMesh::flipEdge(EdgeIter e0) {
        // TODO Part 3.
        // TODO This method should flip the given edge and return an iterator to the flipped edge.
        if (e0->isBoundary()) {
            return EdgeIter();
        }
        // Get elements
        HalfedgeIter he0 = e0->halfedge();
        EdgeIter e1 = e0;
        FaceIter f0 = he0->face();
        VertexIter c = he0->vertex();
        HalfedgeIter he1 = he0->next();
        EdgeIter e2 = he1->edge();
        HalfedgeIter he2 = he1->next();
        EdgeIter e3 = he2->edge();
        VertexIter d = he2->vertex();
        HalfedgeIter he3 = he2->next()->twin();
        FaceIter f1 = he3->face();
        VertexIter b = he3->vertex();
        HalfedgeIter he4 = he3->next();
        EdgeIter e4 = he4->edge();
        HalfedgeIter he5 = he4->next();
        EdgeIter e5 = he5->edge();
        VertexIter a = he5->vertex();

        // Set elements
        a->halfedge() = he5;
        b->halfedge() = he1;
        c->halfedge() = he4;
        d->halfedge() = he2;
        f0->halfedge() = he0;
        f1->halfedge() = he3;
        e1->halfedge() = he0;
        e2->halfedge() = he1;
        e3->halfedge() = he2;
        e4->halfedge() = he4;
        e5->halfedge() = he5;

        he0->next() = he2;
        he0->vertex() = a;
        he0->edge() = e1;
        he0->face() = f0;

        he1->next() = he3;
        he1->vertex() = b;
        he1->edge() = e2;
        he1->face() = f1;

        he2->next() = he4;
        he2->vertex() = d;
        he2->edge() = e3;
        he2->face() = f0;

        he3->next() = he5;
        he3->vertex() = d;
        he3->edge() = e1;
        he3->face() = f1;

        he4->next() = he0;
        he4->vertex() = c;
        he4->edge() = e4;
        he4->face() = f0;

        he5->next() = he1;
        he5->vertex() = a;
        he5->edge() = e5;
        he5->face() = f1;


        return e0;
    }

    VertexIter HalfedgeMesh::splitEdge(EdgeIter e0) {
        // TODO Part 4.
        // TODO This method should split the given edge and return an iterator to the newly inserted vertex.
        // TODO The halfedge of this vertex should point along the edge that was split, rather than the new edges.
        if (e0->isBoundary()) {
            return VertexIter();
        }
        HalfedgeIter he0 = e0->halfedge();
        EdgeIter e1 = e0;
        FaceIter f0 = he0->face();
        VertexIter c = he0->vertex();
        HalfedgeIter he1 = he0->next();
        EdgeIter e2 = he1->edge();
        HalfedgeIter he2 = he1->next();
        EdgeIter e3 = he2->edge();
        VertexIter d = he2->vertex();
        HalfedgeIter he3 = he2->next()->twin();
        FaceIter f1 = he3->face();
        VertexIter b = he3->vertex();
        HalfedgeIter he4 = he3->next();
        EdgeIter e4 = he4->edge();
        HalfedgeIter he5 = he4->next();
        EdgeIter e5 = he5->edge();
        VertexIter a = he5->vertex();

        Vector3D epos = (c->position+b->position)/2;
        VertexIter e = newVertex();
        e->position = epos;

        // new elements
        FaceIter f2 = newFace();
        FaceIter f3 = newFace();
        EdgeIter e6 = newEdge();
        EdgeIter e7 = newEdge();
        EdgeIter e8 = newEdge();
        HalfedgeIter he6 = newHalfedge();
        HalfedgeIter he7 = newHalfedge();
        HalfedgeIter he8 = newHalfedge();
        HalfedgeIter he9 = newHalfedge();
        HalfedgeIter he10 = newHalfedge();
        HalfedgeIter he11 = newHalfedge();

        // Set elements
        f0->halfedge() = he0;
        f1->halfedge() = he3;
        f2->halfedge() = he7;
        f3->halfedge() = he6;

        e1->halfedge() = he0;
        e2->halfedge() = he1;
        e3->halfedge() = he2;
        e4->halfedge() = he4;
        e5->halfedge() = he5;
        e6->halfedge() = he6;
        e7->halfedge() = he8;
        e8->halfedge() = he10;

        a->halfedge() = he5;
        b->halfedge() = he1;
        c->halfedge() = he4;
        d->halfedge() = he2;
        e->halfedge() = he6;

        he0->next() = he8;
        he0->vertex() = c;
        he0->edge() = e1;
        he0->face() = f0;

        he1->next() = he9;
        he1->vertex() = b;
        he1->edge() = e2;
        he1->face() = f3;

        he2->next() = he0;
        he2->vertex() = d;
        he2->edge() = e3;
        he2->face() = f0;

        he3->next() = he4;
        he3->vertex() = e;
        he3->edge() = e1;
        he3->face() = f1;

        he4->next() = he10;
        he4->vertex() = c;
        he4->edge() = e4;
        he4->face() = f1;

        he5->next() = he7;
        he5->vertex() = a;
        he5->edge() = e5;
        he5->face() = f2;

        he6->setNeighbors(he1,he7,e,e6,f3);
        he7->setNeighbors(he11,he6,b,e6,f2);
        he8->setNeighbors(he2,he9,e,e7,f0);
        he9->setNeighbors(he6,he8,d,e7,f3);
        he10->setNeighbors(he3,he11,a,e8,f1);
        he11->setNeighbors(he5,he10,e,e8,f2);

        return e;
    }

    void MeshResampler::simplify(HalfedgeMesh& mesh){

    }

    void HalfedgeMesh::face_quadric_error() {
        // Itertate through all of the faces in the mesh and set their quadric_error values.
        for (FaceIter f = this->facesBegin(); f != this->facesEnd(); f++) {
            double a = f->normal().x;
            double b = f->normal().y;
            double c = f->normal().z;
            double xo = f->halfedge()->vertex()->position.x;
            double yo = f->halfedge()->vertex()->position.y;
            double zo = f->halfedge()->vertex()->position.z;
            double d = -a*xo - b*yo - c*zo;
            double data[] = {pow(a,2.0),a*b,a*c,a*d,a*b,pow(b,2.0),b*c,b*d,a*c,b*c,pow(c,2.0),c*d,a*d,b*d,c*d,pow(d,2.0)};
          f->quadric = Matrix4x4(data);
        }
    }

    void HalfedgeMesh::vertex_quadric_error() {
        for (VertexIter v = this->verticesBegin(); v != this->verticesEnd(); v++) {
            Matrix4x4 sum;
            HalfedgeCIter h = v->halfedge();    // get one of the outgoing halfedges of the vertex
              do {
                HalfedgeCIter h_twin = h->twin(); // get the vertex of the current halfedge.
                FaceCIter f = h_twin->face();       // get face of twin halfedge.
                sum += f->quadric;               //  Add to sum.
                h = h_twin->next();               // move to the next outgoing halfedge of the vertex.
              } while(h != v->halfedge());        // keep going until we're back at the beginning
            v->quadric = sum;
        }
    }

    double HalfedgeMesh::cost(VertexIter v0, VertexIter v1) {
        // Vector4D v = Vector4D(v0.x,v0.y,v0.z,1)-Vector4D(v1.x,v1.y,v1.z,1);
        Matrix4x4 Q = (v0->quadric + v1->quadric);
        double data[] = {Q(0,0),Q(0,1),Q(0,2),Q(0,3),Q(1,0),Q(1,1),Q(1,2),Q(1,3),Q(2,0),Q(2,1),Q(2,2),Q(2,3),0,0,0,1};
        Matrix4x4 Qu = Matrix4x4(data);
        Vector4D v = Q.inv()*Vector4D(0,0,0,1);
        return dot(v,(Q*v));
    }

    void MeshResampler::upsample(HalfedgeMesh& mesh)
    // TODO Part 5.
    // This routine should increase the number of triangles in the mesh using Loop subdivision.
    {
        // Each vertex and edge of the original surface can be associated with a vertex in the new (subdivided) surface.
        // Therefore, our strategy for computing the subdivided vertex locations is to *first* compute the new positions
        // using the connectity of the original (coarse) mesh; navigating this mesh will be much easier than navigating
        // the new subdivided (fine) mesh, which has more elements to traverse.  We will then assign vertex positions in
        // the new mesh based on the values we computed for the original mesh.


        // TODO Compute new positions for all the vertices in the input mesh, using the Loop subdivision rule,
        // TODO and store them in Vertex::newPosition. At this point, we also want to mark each vertex as being
        // TODO a vertex of the original mesh.
        for(VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++) {
            v->isNew = false;
            Vector3D neighbor_position_sum;
            Vector3D original_position = v->position;
            double n = 0.0;
            double u;
            VertexIter move = v;
            HalfedgeIter h = move->halfedge();    // get one of the outgoing halfedges of the vertex
                do {
                    HalfedgeIter h_twin = h->twin(); // get the vertex of the current halfedge
                    VertexIter move = h_twin->vertex(); // vertex is 'source' of the half edge.
                                                  // so h->vertex() is v,
                                                  // whereas h_twin->vertex() is the neighbor vertex.
                    neighbor_position_sum += move->position;
                    n+= 1;
                    h = h_twin->next();               // move to the next outgoing halfedge of the vertex.
                } while(h != move->halfedge());        // keep going until we're back at the beginning
            if (n == 3) {
                u = 3.0/16.0;
            } else {
                u = 3.0/(8.0*n);
            }
            v->newPosition = (1.0 - n*u) * original_position + u * neighbor_position_sum;
        }

        // TODO Next, compute the updated vertex positions associated with edges, and store it in Edge::newPosition.
        for(EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
            e->isNew = false;
            if (!e->isBoundary()) {
                HalfedgeIter h = e->halfedge();
                Vector3D A = h->vertex()->position;
                h = h->next();
                Vector3D B = h->vertex()->position;
                h = h->next();
                Vector3D D = h->vertex()->position;
                h = h->next()->twin()->next()->next();
                Vector3D C = h->vertex()->position;
                e->newPosition = ((3.0/8.0) * (A + B)) + ((1.0/8.0) * (C + D));
            }
        }


        // TODO Next, we're going to split every edge in the mesh, in any order.  For future
        // TODO reference, we're also going to store some information about which subdivided
        // TODO edges come from splitting an edge in the original mesh, and which edges are new,
        // TODO by setting the flat Edge::isNew.  Note that in this loop, we only want to iterate
        // TODO over edges of the original mesh---otherwise, we'll end up splitting edges that we
        // TODO just split (and the loop will never end!)
        EdgeIter e = mesh.edgesBegin();
        while (e != mesh.edgesEnd()) {

          // get the next edge NOW!
          EdgeIter nextEdge = e;
          nextEdge++;
          VertexIter v1 = e->halfedge()->vertex();
          VertexIter v2 = e->halfedge()->next()->vertex();

          // now, even if splitting the edge deletes it...
          if ((!e->isNew) && (!v1->isNew && !v2->isNew)) {
            VertexIter v = mesh.splitEdge(e);
            v->isNew = true;
            v->position = e->newPosition;
            HalfedgeIter s = v->halfedge();
            s = s->next()->next();
            s->edge()->isNew = true;
            s = s->next()->twin()->next();
            s->edge()->isNew = true;
          }

          // ...we still have a valid reference to the next edge.
          e = nextEdge;
        }


        // TODO Now flip any new edge that connects an old and new vertex.
        EdgeIter ef = mesh.edgesBegin();
        while (ef != mesh.edgesEnd()) {

          // get the next edge NOW!
          EdgeIter nextEdgef = ef;
          nextEdgef++;
          VertexIter v1 = ef->halfedge()->vertex();
          VertexIter v2 = ef->halfedge()->next()->vertex();

          // now, even if splitting the edge deletes it...
          if ((ef->isNew) &&  ((v1->isNew && !v2->isNew) || (v2->isNew && !v1->isNew))) {
            mesh.flipEdge(ef);
          }

          // ...we still have a valid reference to the next edge.
          ef = nextEdgef;
        }


        // TODO Finally, copy the new vertex positions into final Vertex::position.
        for(VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++) {
            if (!v->isNew) {
                v->position = v->newPosition;
            }
        }

    }

    // TODO Part 6.
    // TODO There's also some code you'll need to complete in "Shader/frag" file.

    Polymesh BPA_hardcode(std::vector<Vector3D>& vertices) {
      cout << "BPA yo, voyteces: " << vertices.size() << endl;


      // Build the HalfedgeMesh based on what we got
      Polymesh pm;
      PolyList polygons; // a std::vector of Polygon structs
      std::vector<Vector3D> hc_vertices;


      // Let's hard code to test!!
      hc_vertices.push_back(Vector3D(1.0, 0.0, 0.0));
      hc_vertices.push_back(Vector3D(0.0, 1.0, 0.0));
      hc_vertices.push_back(Vector3D(0.0, 0.0, 1.0));
      hc_vertices.push_back(Vector3D(-1.0, 0.0, 0.0));
      hc_vertices.push_back(Vector3D(0.0, -1.0, 0.0));

      Polygon polygon1;
      polygon1.vertex_indices.push_back(0);
      polygon1.vertex_indices.push_back(1);
      polygon1.vertex_indices.push_back(2);

      Polygon polygon2;
      polygon2.vertex_indices.push_back(1);
      polygon2.vertex_indices.push_back(3);
      polygon2.vertex_indices.push_back(2);

      Polygon polygon3;
      polygon3.vertex_indices.push_back(3);
      polygon3.vertex_indices.push_back(4);
      polygon3.vertex_indices.push_back(2);

      Polygon polygon4;
      polygon4.vertex_indices.push_back(4);
      polygon4.vertex_indices.push_back(0);
      polygon4.vertex_indices.push_back(2);

      pm.vertices = hc_vertices;
      pm.polygons.push_back(polygon1);
      pm.polygons.push_back(polygon2);
      pm.polygons.push_back(polygon3);
      pm.polygons.push_back(polygon4);
      return pm;
    }

}
