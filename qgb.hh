#pragma once

#include <geometry.hh>

// 1. Set up interpolants
// 2. Set up midpoint
// 3. Evaluate
// The patch can be directly re-evaluated after modification.

// Boundary curve CPs are assumed to be in the correct order and orientation.

class QGB {
public:
  using Boundary = std::array<Geometry::Point3D, 3>;
  // Constructor
  QGB(size_t n);

  // Getters & setters
  size_t size() const;
  Boundary boundary(size_t i) const;
  void setBoundary(size_t i, const Boundary &b);
  Geometry::Point3D midpoint() const;
  void setMidpoint(const Geometry::Point3D &p);
  void resetMidpoint();

  // Evaluation
  Geometry::Point3D eval(const Geometry::Point2D &uv, double *deficiency = nullptr) const;
  Geometry::TriMesh eval(size_t resolution) const;

private:
  void updateCentralControlPoint();

  size_t n_;
  Geometry::Point3D central_cp_, midpoint_;
  std::vector<Boundary> boundaries_;

  class Domain {
  public:
    Domain(size_t n);
    Geometry::Point2D center() const;
    Geometry::DoubleVector barycentric(const Geometry::Point2D &uv) const;
    Geometry::Point2DVector parameters(size_t resolution) const;
    Geometry::TriMesh meshTopology(size_t resolution) const;
    bool onEdge(size_t resolution, size_t index) const;

  private:
    size_t n_;
    Geometry::Point2DVector points_;
  } domain_;
};
