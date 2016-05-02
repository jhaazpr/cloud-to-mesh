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

    VertexIter HalfedgeMesh::collapseEdge(EdgeIter e0){
        if (e0->isBoundary()) {
            return VertexIter();
        }
        // HalfedgeIter h1 = e0->halfedge();
        // HalfedgeIter h2 = h1->next();
        // HalfedgeIter h3 = h2->next();
        // HalfedgeIter h4 = h3->twin();
        // HalfedgeIter h5 = h4->next();
        // HalfedgeIter h6 = h5->next();
        // HalfedgeIter h7 = h6->twin();
        // HalfedgeIter h8 = h7->next();
        // HalfedgeIter h9 = h8->next();
        // HalfedgeIter h10 = h9->twin();
        // HalfedgeIter h11 = h10->next();
        // HalfedgeIter h12 = h11->next();
        // HalfedgeIter h13 = h12->twin();
        // HalfedgeIter h14 = h13->next();
        // HalfedgeIter h15 = h14->next();
        // HalfedgeIter h16 = h15->twin();
        // HalfedgeIter h17 = h16->next();
        // HalfedgeIter h18 = h17->next();
        // HalfedgeIter h19 = h17->twin();
        // HalfedgeIter h20 = h19->next();
        // HalfedgeIter h21 = h20->next();
        // HalfedgeIter h22 = h21->twin();
        // HalfedgeIter h23 = h22->next();
        // HalfedgeIter h24 = h23->next();
        // HalfedgeIter h25 = h24->twin();
        // HalfedgeIter h26 = h25->next();
        // HalfedgeIter h27 = h26->next();
        // HalfedgeIter h28 = h27->twin();
        // HalfedgeIter h29 = h28->next();
        // HalfedgeIter h30 = h29->next();


        // FaceIter f1 = h1->face();
        // FaceIter f2 = h18->face();
        // FaceIter f3 = h30->face();
        // FaceIter f4 = h4->face();
        // FaceIter f5 = h14->face();
        // FaceIter f6 = h19->face();
        // FaceIter f7 = h25->face();
        // FaceIter f8 = h24->face();
        // FaceIter f9 = h7->face();
        // FaceIter f10 = h10->face();

        // EdgeIter e1 = h2->edge();
        // EdgeIter e2 = h3->edge();
        // EdgeIter e3 = h17->edge();
        // EdgeIter e4 = h16->edge();
        // EdgeIter e5 = h14->edge();
        // EdgeIter e6 = h12->edge();
        // EdgeIter e7 = h5->edge();
        // EdgeIter e8 = h6->edge();
        // EdgeIter e9 = h21->edge();
        // EdgeIter e10 = h20->edge();
        // EdgeIter e11 = h28->edge();
        // EdgeIter e12 = h29->edge();
        // EdgeIter e13 = h26->edge();
        // EdgeIter e14 = h23->edge();
        // EdgeIter e15 = h11->edge();
        // EdgeIter e16 = h8->edge();
        // EdgeIter e17 = h9->edge();
        // EdgeIter e18 = h25->edge();

        // VertexIter v1 = h4->vertex();
        // VertexIter v2 = h18->vertex();
        // VertexIter v3 = h3->vertex();
        // VertexIter v4 = h17->vertex();
        // VertexIter v5 = h12->vertex();
        // VertexIter v6 = h27->vertex();
        // VertexIter v7 = h23->vertex();
        // VertexIter v8 = h6->vertex();
        // VertexIter v9 = h24->vertex();
        // VertexIter v10 = h9->vertex();

        HalfedgeIter h1 = e0->halfedge();
        HalfedgeIter h2 = h1->next();
        HalfedgeIter h3 = h2->next();
        HalfedgeIter h4 = h3->twin();
        HalfedgeIter h7 = h4->next()->next()->twin();
        HalfedgeIter h10 = h7->next()->next()->twin();
        HalfedgeIter h13 = h10->next()->next()->twin();
        HalfedgeIter h14 = h13->next();
        HalfedgeIter h15 = h14->next();
        HalfedgeIter h18 = h1->twin();
        HalfedgeIter h16 = h18->next();
        HalfedgeIter h17 = h16->next();
        HalfedgeIter h19 = h17->twin();
        HalfedgeIter h22 = h19->next()->next()->twin();
        HalfedgeIter h25 = h22->next()->next()->twin();
        HalfedgeIter h28 = h25->next()->next()->twin();
        HalfedgeIter h29 = h28->next();
        HalfedgeIter h30 = h29->next();

        FaceIter f1 = h1->face();
        FaceIter f2 = h18->face();
        FaceIter f3 = h30->face();
        FaceIter f5 = h15->face();

        EdgeIter e1 = h2->edge();
        EdgeIter e4 = h16->edge();

        VertexIter v1 = h1->vertex();
        VertexIter v2 = h18->vertex();
        VertexIter v3 = h3->vertex();
        VertexIter v4 = h17->vertex();


        h4->vertex() = v1;
        h7->vertex() = v1;
        h10->vertex() = v1;
        h13->vertex() = v1;
        h19->vertex() = v1;
        h22->vertex() = v1;
        h25->vertex() = v1;
        h28->vertex() = v1;

        h3->next() = h28;
        h29->next() = h3;

        h17->next() = h13;
        h14->next() = h17;

        h3->face() = f3;
        h17->face() = f5;

        f3->halfedge() = h3;
        f5->halfedge() = h17;

        v1->halfedge() = h4;
        v3->halfedge() = h3;
        v4->halfedge() = h17;

        Vector3D p = (v1->position+v2->position)/2.0;
        v1->position = p;
        v2->position = p;

        deleteVertex(v2);

        deleteFace(f1);
        deleteFace(f2);

        deleteEdge(e0);
        deleteEdge(e1);
        deleteEdge(e4);

        deleteHalfedge(h2);
        deleteHalfedge(h30);
        deleteHalfedge(h16);
        deleteHalfedge(h15);
        // deleteHalfedge(h18);
        // deleteHalfedge(h1);
        
        
        return v1;
        // return VertexIter();
    }

    EdgeRecord::EdgeRecord( EdgeIter& _edge )
   : edge( _edge )
   {
        VertexIter v1 = _edge->halfedge()->vertex();
        VertexIter v2 = _edge->halfedge()->twin()->vertex();
        Matrix4x4 K = v1->quadric+v2->quadric;
        double data[] = {K(0,0),K(0,1),K(0,2),K(1,0),K(1,1),K(1,2),K(2,0),K(2,1),K(2,2)};
        Matrix3x3 A = Matrix3x3(data);
        Vector3D b = Vector3D(-K(0,3),-K(1,3),-K(2,3));
        Vector3D x = A.inv() * b;
        optimalPoint = x;
        score = dot(x,(K*x));
        edge = _edge;
    }

    void MeshResampler::simplify( HalfedgeMesh& mesh )
   {
      // TODO Compute initial quadrics for each face by simply writing the plane
      // equation for the face in homogeneous coordinates.  These quadrics should
      // be stored in Face::quadric
        for (FaceIter f = mesh.facesBegin(); f != mesh.facesEnd(); f++) {
            double a = f->normal().x;
            double b = f->normal().y;
            double c = f->normal().z;
            double xo = f->halfedge()->vertex()->position.x;
            double yo = f->halfedge()->vertex()->position.y;
            double zo = f->halfedge()->vertex()->position.z;
            double d = -a*xo - b*yo - c*zo;
            Vector4D v = Vector4D(a,b,c,d);
            f->quadric = outer(v,v);
        }


      // TODO Compute an initial quadric for each vertex as the sum of the quadrics
      // associated with the incident faces, storing it in Vertex::quadric
        for (VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++) {
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


      // TODO Build a priority queue of edges according to their quadric error cost,
      // TODO i.e., by building an EdgeRecord for each edge and sticking it in the queue.

        MutablePriorityQueue<EdgeRecord> queue;
        for(EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
            EdgeRecord myRecord(e);
            queue.insert(myRecord);
        }
        
        Size Start = mesh.nFaces();
        // while(mesh.nFaces() >= 10) {
            cout << Start << endl;
            // Get the cheapest edge from the queue.
            EdgeRecord bestEdge = queue.top();
            // Remove the cheapest edge from the queue by calling pop().
            queue.pop();
            // Compute the new quadric by summing the quadrics at its two endpoints.
            VertexIter v1 = bestEdge.edge->halfedge()->vertex();
            VertexIter v2 = bestEdge.edge->halfedge()->twin()->vertex();
            Matrix4x4 K = v1->quadric+v2->quadric; 
            VertexIter verts[] = {v1,v2};
            // Remove any edge touching either of its endpoints from the queue.
            for (VertexIter vr: verts) {
                HalfedgeIter h1 = vr->halfedge();    // get one of the outgoing halfedges of the vertex
                do {
                    HalfedgeIter h1_twin = h1->twin(); // get the vertex of the current halfedge
                    VertexIter vr = h1_twin->vertex(); // vertex is 'source' of the half edge.
                    queue.remove(h1_twin->edge());
                    h1 = h1_twin->next();               // move to the next outgoing halfedge of the vertex.
                  } while(h1 != vr->halfedge());        // keep going until we're back at the beginning
            }
            // Collapse the edge.
            VertexIter v = mesh.collapseEdge(bestEdge.edge);
            mesh.edgesBegin();
            cout << Start << endl;
            v->position = bestEdge.optimalPoint;
            // Set the quadric of the new vertex to the quadric computed in Step 2.
            v->quadric = K;
            // Insert any edge touching the new vertex into the queue, creating new edge records for each of them.
            HalfedgeIter h = v->halfedge();    // get one of the outgoing halfedges of the vertex
            do {
                HalfedgeIter h_twin = h->twin(); // get the vertex of the current halfedge
                VertexIter v = h_twin->vertex(); // vertex is 'source' of the half edge.
                EdgeRecord myRecord(h_twin->edge());
                queue.insert(myRecord);
                h = h_twin->next();               // move to the next outgoing halfedge of the vertex.
              } while(h != v->halfedge());        // keep going until we're back at the beginning
        // }

      // TODO Until we reach the target edge budget, collapse the best edge.  Remember
      // TODO to remove from the queue any edge that touches the collapsing edge BEFORE
      // TODO it gets collapsed, and add back into the queue any edge touching the collapsed
      // TODO vertex AFTER it's been collapsed.  Also remember to assign a quadric to the
      // TODO collapsed vertex, and to pop the collapsed edge off the top of the queue.
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
