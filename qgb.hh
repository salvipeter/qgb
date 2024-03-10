#pragma once

#include "domain.hh"

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
  Boundary boundary(size_t i) const;
  void setBoundary(size_t i, const Boundary &b);
  Geometry::Point3D midpoint() const;
  void setMidpoint(const Geometry::Point3D &p);
  void resetMidpoint();
  const Domain *domain() const;

  // Evaluation
  void updateDomain();
  Geometry::Point3D eval(const Geometry::Point2D &uv, double *deficiency = nullptr) const;
  Geometry::TriMesh eval(size_t resolution) const;

private:
  void updateCentralControlPoint();

  size_t n_;
  Geometry::Point3D central_cp_, midpoint_;
  std::vector<Boundary> boundaries_;
  std::unique_ptr<Domain> domain_;
};
