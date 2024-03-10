#include "qgb.hh"

#include <algorithm>
#include <cmath>
#include <numeric>

using namespace Geometry;


// Constructor

QGB::QGB(size_t n) : n_(n), domain_(n)
{
  boundaries_.resize(n_);
}


// Constraint modifications

size_t
QGB::size() const {
  return n_;
}

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
  Point2D center = domain_.center();
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
  auto bc = domain_.barycentric(uv);
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
  TriMesh mesh = domain_.meshTopology(resolution);
  Point2DVector uvs = domain_.parameters(resolution);
  PointVector points; points.reserve(uvs.size());
  std::transform(uvs.begin(), uvs.end(), std::back_inserter(points),
                 [&](const Point2D &uv) { return eval(uv); });
  mesh.setPoints(points);
  return mesh;
}

QGB::Domain::Domain(size_t n) : n_(n) {
  if (n_ == 4)
    points_ = { {1,1}, {-1,1}, {-1,-1}, {1,-1} };
  else {
    double alpha = 2.0 * M_PI / n_;
    for (size_t i = 0; i < n_; ++i)
      points_.emplace_back(std::cos(alpha * i), std::sin(alpha * i));
  }
}

Point2D
QGB::Domain::center() const {
  return { 0, 0 };
}

DoubleVector
QGB::Domain::barycentric(const Point2D &uv) const {
  Vector2DVector vectors; vectors.reserve(n_);
  std::transform(points_.begin(), points_.end(), std::back_inserter(vectors),
                 [uv](const Point2D &p) { return uv - p; });

  DoubleVector areas; areas.reserve(n_);
  for (size_t i = 0; i < n_; ++i) {
    const Vector2D &si = vectors[i];
    const Vector2D &si1 = vectors[(i+1)%n_];
    areas.push_back((si[0] * si1[1] - si[1] * si1[0]) / 2.0);
  }

  DoubleVector l; l.reserve(n_);

  for (size_t i = 0; i < n_; ++i) {
    size_t i_1 = (i + n_ - 1) % n_, i1 = (i + 1) % n_;
    double Ai = 1.0, Ai_1 = 1.0, Ai_1i = 1.0;
    for (size_t j = 0; j < n_; ++j) {
      if (j == i)
        Ai_1 *= areas[j];
      else if (j == i_1)
        Ai *= areas[j];
      else {
        Ai_1 *= areas[j];
        Ai *= areas[j];
        Ai_1i *= areas[j];
      }
    }
    const Vector2D &si_1 = vectors[i_1];
    const Vector2D &si1 = vectors[i1];
    double Bi = (si_1[0] * si1[1] - si_1[1] * si1[0]) / 2.0;
    l.push_back(Ai_1 + Ai - Bi * Ai_1i);
  }

  double sum = std::accumulate(l.begin(), l.end(), 0.0);
  std::transform(l.begin(), l.end(), l.begin(), [sum](double x) { return x / sum; });
  return l;
}

static size_t meshSize(size_t n, size_t resolution) {
  if (n == 3)
    return (resolution + 1) * (resolution + 2) / 2;
  if (n == 4)
    return (resolution + 1) * (resolution + 1);
  return 1 + n * resolution * (resolution + 1) / 2;
}

Point2DVector
QGB::Domain::parameters(size_t resolution) const {
  size_t size = meshSize(n_, resolution);
  Point2DVector result;
  result.reserve(size);

  if (n_ == 3) {
    for (size_t j = 0; j <= resolution; ++j) {
      double u = (double)j / resolution;
      auto p = points_[0] * u + points_[2] * (1 - u);
      auto q = points_[1] * u + points_[2] * (1 - u);
      for (size_t k = 0; k <= j; ++k) {
        double v = j == 0 ? 1.0 : (double)k / j;
        result.push_back(p * (1 - v) + q * v);
      }
    }
  } else if (n_ == 4) {
    for (size_t j = 0; j <= resolution; ++j) {
      double u = (double)j / resolution;
      auto p = points_[0] * (1 - u) + points_[1] * u;
      auto q = points_[3] * (1 - u) + points_[2] * u;
      for (size_t k = 0; k <= resolution; ++k) {
        double v = (double)k / resolution;
        result.push_back(p * (1 - v) + q * v);
      }
    }
  } else { // n_ > 4
    Point2D center(0.0, 0.0);
    result.push_back(center);
    for (size_t j = 1; j <= resolution; ++j) {
      double u = (double)j / (double)resolution;
      for (size_t k = 0; k < n_; ++k)
        for (size_t i = 0; i < j; ++i) {
          double v = (double)i / (double)j;
          Point2D ep = points_[(k+n_-1)%n_] * (1.0 - v) + points_[k] * v;
          Point2D p = center * (1.0 - u) + ep * u;
          result.push_back(p);
        }
    }
  }
  return result;
}

TriMesh
QGB::Domain::meshTopology(size_t resolution) const {
  TriMesh mesh;
  mesh.resizePoints(meshSize(n_, resolution));

  if (n_ == 3) {
    size_t prev = 0, current = 1;
    for (size_t i = 0; i < resolution; ++i) {
      for (size_t j = 0; j < i; ++j) {
        mesh.addTriangle(current + j, current + j + 1, prev + j);
        mesh.addTriangle(current + j + 1, prev + j + 1, prev + j);
      }
      mesh.addTriangle(current + i, current + i + 1, prev + i);
      prev = current;
      current += i + 2;
    }
  } else if (n_ == 4) {
    for (size_t i = 0; i < resolution; ++i)
      for (size_t j = 0; j < resolution; ++j) {
        size_t index = i * (resolution + 1) + j;
        mesh.addTriangle(index, index + resolution + 1, index + 1);
        mesh.addTriangle(index + 1, index + resolution + 1, index + resolution + 2);
      }
  } else { // n_ > 4
    size_t inner_start = 0, outer_vert = 1;
    for (size_t layer = 1; layer <= resolution; ++layer) {
      size_t inner_vert = inner_start, outer_start = outer_vert;
      for (size_t side = 0; side < n_; ++side) {
        size_t vert = 0;
        while(true) {
          size_t next_vert = (side == n_ - 1 && vert == layer - 1) ? outer_start : (outer_vert + 1);
          mesh.addTriangle(inner_vert, outer_vert, next_vert);
          ++outer_vert;
          if (++vert == layer)
            break;
          size_t inner_next = (side == n_ - 1 && vert == layer - 1) ? inner_start : (inner_vert + 1);
          mesh.addTriangle(inner_vert, next_vert, inner_next);
          inner_vert = inner_next;
        }
      }
      inner_start = outer_start;
    }
  }
  return mesh;
}

bool
QGB::Domain::onEdge(size_t resolution, size_t index) const {
  if (n_ == 3) {
    if (index >= meshSize(3, resolution) - resolution - 1)
      return true;
    auto issquare = [](size_t n) {
                      size_t root = std::round(std::sqrt(n));
                      return root * root == n;
                    };
    size_t n = index * 8 + 1;
    return issquare(n) || issquare(n + 8);
  }
  if (n_ == 4) {
    return index <= resolution || index >= (resolution + 1) * resolution ||
      index % (resolution + 1) == 0 || index % (resolution + 1) == resolution;
  }
  return index >= meshSize(n_, resolution) - n_ * resolution;
}
