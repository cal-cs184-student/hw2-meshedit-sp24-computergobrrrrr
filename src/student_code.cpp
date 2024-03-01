#include "student_code.h"
#include "CGL/vector3D.h"
#include "halfEdgeMesh.h"
#include "mutablePriorityQueue.h"

using namespace std;

namespace CGL
{

  /**
   * Evaluates one step of the de Casteljau's algorithm using the given points and
   * the scalar parameter t (class member).
   *
   * @param points A vector of points in 2D
   * @return A vector containing intermediate points or the final interpolated vector
   */
  std::vector<Vector2D> BezierCurve::evaluateStep(std::vector<Vector2D> const &points)
  { 
    // DONE Part 1.

    int numPoints = points.size();
    if (numPoints == 1) return points;

    std::vector<Vector2D> resPoints;
    for (int i = 0; i < numPoints - 1; i++) {
      resPoints.push_back(points[i] * (1 - t) + points[i + 1] * t);
    }

    return resPoints;
  }

  /**
   * Evaluates one step of the de Casteljau's algorithm using the given points and
   * the scalar parameter t (function parameter).
   *
   * @param points    A vector of points in 3D
   * @param t         Scalar interpolation parameter
   * @return A vector containing intermediate points or the final interpolated vector
   */
  std::vector<Vector3D> BezierPatch::evaluateStep(std::vector<Vector3D> const &points, double t) const
  {
    // DONE Part 2.

    int numPoints = points.size();
    if (numPoints == 1) return points;

    std::vector<Vector3D> resPoints;
    for (int i = 0; i < numPoints - 1; i++) {
      resPoints.push_back(points[i] * (1 - t) + points[i + 1] * t);
    }

    return resPoints;
  }

  /**
   * Fully evaluates de Casteljau's algorithm for a vector of points at scalar parameter t
   *
   * @param points    A vector of points in 3D
   * @param t         Scalar interpolation parameter
   * @return Final interpolated vector
   */
  Vector3D BezierPatch::evaluate1D(std::vector<Vector3D> const &points, double t) const
  {
    // DONE Part 2.

    int numPoints = points.size();
    if (numPoints == 1) return points[0];

    std::vector<Vector3D> resPoints = evaluateStep(points, t);
    return evaluate1D(resPoints, t);
  }

  /**
   * Evaluates the Bezier patch at parameter (u, v)
   *
   * @param u         Scalar interpolation parameter
   * @param v         Scalar interpolation parameter (along the other axis)
   * @return Final interpolated vector
   */
  Vector3D BezierPatch::evaluate(double u, double v) const 
  {  
    // DONE Part 2.

    std::vector<Vector3D> resPoints;
    for (int i = 0; i < controlPoints.size() ; i++) {
      CGL::Vector3D tempPoints = evaluate1D(controlPoints[i], u);
      resPoints.push_back(tempPoints);      
    }

    return evaluate1D(resPoints, v);
  }

  double calculateArea(Vector3D a, Vector3D b, Vector3D c) {
    Vector3D ba = b - a;
    Vector3D ca = c - a;
    return cross(ba, ca).norm() / 2;
  }

  Vector3D Vertex::normal( void ) const
  {
    // DONE Part 3.
    // Returns an approximate unit normal at this vertex, computed by
    // taking the area-weighted average of the normals of neighboring
    // triangles, then normalizing.
    HalfedgeCIter endHalfEdge = halfedge();
    HalfedgeCIter h = halfedge();
    Vector3D resVec(0, 0, 0);

    double area = calculateArea(
      h->vertex()->position,
      h->next()->next()->vertex()->position, // Reverse order of normals (from spec?)
      h->next()->vertex()->position
    );
    Vector3D normal = h->face()->normal();
    resVec += normal * area;

    while (h->twin()->next() != endHalfEdge) {
      h = h->twin()->next();
      if (h->face()->isBoundary()) continue;
      double area = calculateArea(
        h->vertex()->position,
        h->next()->next()->vertex()->position, // Reverse order of normals (from spec?)
        h->next()->vertex()->position
      );
      Vector3D normal = h->face()->normal();
      resVec += normal * area; // This is the weighted area calc
    }

    return resVec.unit();
  }

  EdgeIter HalfedgeMesh::flipEdge( EdgeIter e0 )
  {
    // TODO Part 4.
    // This method should flip the given edge and return an iterator to the flipped edge.

    if (e0->isBoundary()) return e0;

    // First group inner halfedges
    HalfedgeIter inner_h0 = e0->halfedge();
    HalfedgeIter inner_h1 = inner_h0->next();
    HalfedgeIter inner_h2 = inner_h1->next();

    // Second group inner halfedges
    HalfedgeIter inner_h3 = inner_h0->twin();
    HalfedgeIter inner_h4 = inner_h3->next();
    HalfedgeIter inner_h5 = inner_h4->next();

    // Outer halfedges
    HalfedgeIter outer_h0 = inner_h1->twin();
    HalfedgeIter outer_h1 = inner_h2->twin();
    HalfedgeIter outer_h2 = inner_h4->twin();
    HalfedgeIter outer_h3 = inner_h5->twin();

    // Full edges
    EdgeIter e1 = inner_h1->edge();
    EdgeIter e2 = inner_h2->edge();
    EdgeIter e3 = inner_h4->edge();
    EdgeIter e4 = inner_h5->edge();

    // Vertices
    VertexIter v0 = inner_h0->vertex();
    VertexIter v1 = inner_h2->vertex();
    VertexIter v2 = inner_h3->vertex();
    VertexIter v3 = inner_h5->vertex();

    FaceIter f0 = inner_h0->face();
    FaceIter f1 = inner_h3->face();

    inner_h0->setNeighbors(inner_h1,inner_h3,v1,e0,f0);
    inner_h1->setNeighbors(inner_h2,outer_h3,v3,e4,f0);
    inner_h2->setNeighbors(inner_h0,outer_h0,v2,e1,f0);
    inner_h3->setNeighbors(inner_h4,inner_h0,v3,e0,f1);
    inner_h4->setNeighbors(inner_h5,outer_h1,v1,e2,f1);
    inner_h5->setNeighbors(inner_h3,outer_h2,v0,e3,f1);

    FaceIter outer_h0_face = outer_h0->face();
    FaceIter outer_h1_face = outer_h1->face();
    FaceIter outer_h2_face = outer_h2->face();
    FaceIter outer_h3_face = outer_h3->face();
    outer_h0->setNeighbors(outer_h0->next(),inner_h2,v1,e1,outer_h0_face);
    outer_h1->setNeighbors(outer_h1->next(),inner_h4,v0,e2,outer_h1_face);
    outer_h2->setNeighbors(outer_h2->next(),inner_h5,v3,e3,outer_h2_face);
    outer_h3->setNeighbors(outer_h3->next(),inner_h1,v2,e4,outer_h3_face);

    e0->halfedge() = inner_h0;
    e1->halfedge() = inner_h2;
    e2->halfedge() = inner_h4;
    e3->halfedge() = inner_h5;
    e4->halfedge() = inner_h1;

    v0->halfedge() = inner_h5;
    v1->halfedge() = inner_h4;
    v2->halfedge() = inner_h2;
    v3->halfedge() = inner_h1;

    f0->halfedge() = inner_h0;
    f1->halfedge() = inner_h3;

    return e0;
  }

  VertexIter HalfedgeMesh::splitEdge( EdgeIter e0 )
  {
    // DONE Part 5. LETS GOOOOOOOO
    // This method should split the given edge and return an iterator to the newly inserted vertex.
    // The halfedge of this vertex should point along the edge that was split, rather than the new edges.

    if (e0->isBoundary()) return e0->halfedge()->vertex();

    // New elements created
    HalfedgeIter new_h0 = newHalfedge();
    HalfedgeIter new_h1 = newHalfedge();
    HalfedgeIter new_h2 = newHalfedge();
    HalfedgeIter new_h3 = newHalfedge();
    HalfedgeIter new_h4 = newHalfedge();
    HalfedgeIter new_h5 = newHalfedge();

    EdgeIter new_e0 = newEdge();
    EdgeIter new_e1 = newEdge();
    EdgeIter new_e2 = newEdge();
    new_e0->isNew = true;
    new_e1->isNew = true;
    new_e2->isNew = true;

    VertexIter new_v0 = newVertex();
    new_v0->isNew = true;

    FaceIter new_f0 = newFace();
    FaceIter new_f1 = newFace();

    // First group inner halfedges
    HalfedgeIter inner_h0 = e0->halfedge();
    HalfedgeIter inner_h1 = inner_h0->next();
    HalfedgeIter inner_h2 = inner_h1->next();

    // Second group inner halfedges
    HalfedgeIter inner_h3 = inner_h0->twin();
    HalfedgeIter inner_h4 = inner_h3->next();
    HalfedgeIter inner_h5 = inner_h4->next();

    // Outer halfedges
    HalfedgeIter outer_h0 = inner_h1->twin();
    HalfedgeIter outer_h1 = inner_h2->twin();
    HalfedgeIter outer_h2 = inner_h4->twin();
    HalfedgeIter outer_h3 = inner_h5->twin();

    // Full edges
    EdgeIter e1 = inner_h1->edge();
    EdgeIter e2 = inner_h2->edge();
    EdgeIter e3 = inner_h4->edge();
    EdgeIter e4 = inner_h5->edge();

    // Vertices
    VertexIter v0 = inner_h0->vertex();
    VertexIter v1 = inner_h2->vertex();
    VertexIter v2 = inner_h3->vertex();
    VertexIter v3 = inner_h5->vertex();

    FaceIter f0 = inner_h0->face();
    FaceIter f1 = inner_h3->face();

    inner_h0->setNeighbors(inner_h1, inner_h3, new_v0, e0, f0);
    inner_h1->setNeighbors(inner_h2, outer_h0, v2, e1, f0);
    inner_h2->setNeighbors(inner_h0, new_h1, v1, new_e0, f0);

    inner_h3->setNeighbors(inner_h4,inner_h0,v2,e0,f1);
    inner_h4->setNeighbors(inner_h5,new_h5,new_v0,new_e2,f1);
    inner_h5->setNeighbors(inner_h3,outer_h3,v3,e4,f1);

    FaceIter outer_h0_face = outer_h0->face();
    FaceIter outer_h1_face = outer_h1->face();
    FaceIter outer_h2_face = outer_h2->face();
    FaceIter outer_h3_face = outer_h3->face();
    outer_h0->setNeighbors(outer_h0->next(),inner_h1,v1,e1,outer_h0_face);
    outer_h1->setNeighbors(outer_h1->next(),new_h2,v0,e2,outer_h1_face);
    outer_h2->setNeighbors(outer_h2->next(),new_h4,v3,e3,outer_h2_face);
    outer_h3->setNeighbors(outer_h3->next(),inner_h5,v2,e4,outer_h3_face);

    e0->halfedge() = inner_h0;
    e1->halfedge() = inner_h1;
    e2->halfedge() = new_h2;
    e3->halfedge() = new_h4;
    e4->halfedge() = inner_h5;

    v0->halfedge() = new_h0;
    v1->halfedge() = new_h2;
    v2->halfedge() = inner_h1;
    v3->halfedge() = inner_h5;

    f0->halfedge() = inner_h0;
    f1->halfedge() = inner_h3;

    // Assigning new elements
    new_h0->setNeighbors(new_h1,new_h3,v0,new_e1,new_f0);
    new_h1->setNeighbors(new_h2,inner_h2,new_v0,new_e0,new_f0);
    new_h2->setNeighbors(new_h0,outer_h1,v1,e2,new_f0);

    new_h3->setNeighbors(new_h4,new_h0,new_v0,new_e1,new_f1);
    new_h4->setNeighbors(new_h5,outer_h2,v0,e3,new_f1);
    new_h5->setNeighbors(new_h3,inner_h4,v3,new_e2,new_f1);

    new_e0->halfedge() = inner_h2;
    new_e1->halfedge() = new_h0;
    new_e2->halfedge() = inner_h4;

    new_v0->halfedge() = inner_h0;
    new_v0->position = (v0->position + v2->position) / 2;

    new_f0->halfedge() = new_h0;
    new_f1->halfedge() = new_h3;

    // Handling isNew
    e0->isNew = false;
    new_e0->isNew = true;
    new_e1->isNew = false;
    new_e2->isNew = true;

    return new_v0;
  }



  void MeshResampler::upsample( HalfedgeMesh& mesh )
  {
    // TODO Part 6.
    // This routine should increase the number of triangles in the mesh using Loop subdivision.
    // One possible solution is to break up the method as listed below.

    VertexIter v;
    EdgeIter e;

    // 1. Compute new positions for all the vertices in the input mesh, using the Loop subdivision rule,
    // and store them in Vertex::newPosition. At this point, we also want to mark each vertex as being
    // a vertex of the original mesh.
    for (e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
      HalfedgeIter h0 = e->halfedge();
      HalfedgeIter h1 = h0->next();
      HalfedgeIter h2 = h1->next();
      HalfedgeIter h3 = h0->twin();
      HalfedgeIter h4 = h3->next();
      HalfedgeIter h5 = h4->next();

      VertexIter v0 = h0->vertex();
      VertexIter v1 = h2->vertex();
      VertexIter v2 = h3->vertex();
      VertexIter v3 = h5->vertex();

      e->newPosition = 0.375 * (v0->position + v2->position) + 0.125 * (v1->position + v3->position);
      e->isNew = false;
    }
    
    // 2. Compute the updated vertex positions associated with edges, and store it in Edge::newPosition.
    for (v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++) { 
      Vector3D original_neighbor_position_sum(0,0,0);
      HalfedgeIter h = v->halfedge();

      original_neighbor_position_sum += h->twin()->vertex()->position;
      while (h->twin()->next() != v->halfedge()) {
        h = h->twin()->next();
        original_neighbor_position_sum += h->twin()->vertex()->position;
      }
      
      float n = (float) v->degree();
      float u;
      if (n == 3.0) u = 0.1875;
      else u = 0.375 / n;
      
      v->newPosition = (1 - n * u) * v->position + u * original_neighbor_position_sum;  
      v->isNew = false;
    }

    // 3. Split every edge in the mesh, in any order. For future reference, we're also going to store some
    // information about which subdivide edges come from splitting an edge in the original mesh, and which edges
    // are new, by setting the flat Edge::isNew. Note that in this loop, we only want to iterate over edges of
    // the original mesh---otherwise, we'll end up splitting edges that we just split (and the loop will never end!)
    for (e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
      VertexIter v0 = e->halfedge()->vertex();
      VertexIter v1 = e->halfedge()->twin()->vertex();

      if (v0->isNew || v1->isNew) continue;

      VertexIter v = mesh.splitEdge(e);
      v->newPosition = e->newPosition;
    }

    // 4. Flip any new edge that connects an old and new vertex.
    for (e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
      VertexIter v0 = e->halfedge()->vertex();
      VertexIter v1 = e->halfedge()->twin()->vertex();
      if (e->isNew) {
        if (v0->isNew != v1->isNew) mesh.flipEdge(e);
      }
    }

    // 5. Copy the new vertex positions into final Vertex::position.
    for (v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++) {
      v->position = v->newPosition;
      v->isNew = false;
    }
  }
}
