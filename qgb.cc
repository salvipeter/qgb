#include "qgb.hh"
#include "regular-domain.hh"

#include <algorithm>
#include <cmath>
#include <numeric>

using namespace Geometry;


// Constructor

QGB::QGB(size_t n) : n_(n)
{
  boundaries_.resize(n_);
  domain_ = std::make_unique<RegularDomain>(n_);
}


// Constraint modifications

QGB::Boundary
QGB::boundary(size_t i) const {
  return boundaries_[i];
}

void
QGB::setBoundary(size_t i, const Boundary &b) {
  boundaries_[i] = b;
}

void
QGB::updateCentralControlPoint() {
  Point2D center = domain_->center();
  central_cp_ = { 0, 0, 0 };
  double def;
  auto s = eval(center, &def);
  if (std::abs(def) < epsilon)
    central_cp_ = midpoint_;  // as good as anything else
  else
    central_cp_ = (midpoint_ - s) / def;
}

Geometry::Point3D
QGB::midpoint() const {
  return midpoint_;
}

void
QGB::setMidpoint(const Point3D &p) {
  midpoint_ = p;
  updateCentralControlPoint();
}

void
QGB::resetMidpoint() {
  midpoint_ = Point3D(0,0,0);
  for (size_t i = 0; i < n_; ++i)
    midpoint_ += boundaries_[i][1];
  midpoint_ /= n_;
  updateCentralControlPoint();
}

const Domain *
QGB::domain() const {
  return domain_.get();
}


// Evaluation

static Point2DVector localParameters(const DoubleVector &bc) {
  size_t n = bc.size();
  Point2DVector sds;
  for (size_t i = 0; i < n; ++i) {
    size_t im = (i + n - 1) % n;
    double s = 0.5;
    double denom = bc[im] + bc[i];
    if (denom > epsilon)
      s = bc[i] / denom;
    sds.emplace_back(s, 1 - denom);
  }
  return sds;
}

Point3D
QGB::eval(const Point2D &uv, double *deficiency) const {
  auto bc = domain_->barycentric(uv);
  auto sds = localParameters(bc);
  double bsum = 0.0;
  Point3D p(0, 0, 0);
  for (size_t i = 0; i < n_; ++i) {
    auto d2 = std::pow(1.0 - sds[i][1], 2);
    auto s = sds[i][0];
    auto blend = 0.5 * d2 * (1.0 - s) * (1.0 - s);
    bsum += blend;
    p += boundaries_[i][0] * blend;
    blend = 2.0 * d2 * s * (1.0 - s);
    bsum += blend;
    p += boundaries_[i][1] * blend;
    blend = 0.5 * d2 * s * s;
    bsum += blend;
    p += boundaries_[i][2] * blend;
  }
  double def = 1.0 - bsum;
  p += central_cp_ * def;
  if (deficiency)
    *deficiency = def;
  return p;
}

TriMesh
QGB::eval(size_t resolution) const {
  TriMesh mesh = domain_->meshTopology(resolution);
  Point2DVector uvs = domain_->parameters(resolution);
  PointVector points; points.reserve(uvs.size());
  std::transform(uvs.begin(), uvs.end(), std::back_inserter(points),
                 [&](const Point2D &uv) { return eval(uv); });
  mesh.setPoints(points);
  return mesh;
}
