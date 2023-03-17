/// @ref geom_polygon
/// @file geom/polygon.hpp
///
/// @defgroup geom_polygon Polygon
/// @ingroup geom
///

#pragma once

#include <vector>

#include "setup.hpp"
#include "line.hpp"
#include "segment.hpp"
#include "aabb.hpp"
#include "plane.hpp"
#include "triangle.hpp"

#include "../vector_extensions.hpp"
#include "../matrix_extensions.hpp"
#if GLM_MESSAGES == GLM_ENABLE && !defined(GLM_EXT_INCLUDED)
  #pragma message("GLM: GLM_EXT_GEOM_polygon extension included")
#endif

GLM_GEOM_BEGIN_NAMESPACE
/// @addtogroup geom_polygon
/// @{

/// <summary>
/// Describes the thickness of the polygon, i.e., how the third dimension
/// relates to the plane for 'contains' operations.
/// </summary>
enum class ePolyContains {
  Positive,  // Boundary extends in the positive direction: [0, +dist]
  Negative,  // Boundary extends in the negative direction: [-dist, 0]
  Unidirectional,  // Boundary extends in both directions: [-0.5*dist, 0.5*dist]
};

/// <summary>
/// A (N - 1) dimension closed surface in N dimension space.
/// </summary>
template<length_t L, typename T, qualifier Q>
struct Polygon {

  // -- Implementation detail --

  typedef T value_type;
  typedef Polygon<L, T, Q> type;
  typedef vec<L, T, Q> point_type;
  using List = std::vector<point_type>;  // @TODO: parameterize

  // -- Data --

  List p;  // Vertices of this polygon

  Polygon()
    : p() {
  }

  Polygon(Polygon<L, T, Q> const &poly)
    : p(poly.p) {
  }

  Polygon<L, T, Q> &operator=(Polygon<L, T, Q> const &poly) {
    p = poly.p;
    return *this;
  }

  GLM_FUNC_QUALIFIER size_t size() const {
    return p.size();
  }

  GLM_FUNC_QUALIFIER point_type const &back() const {
    return p.back();
  }

  GLM_FUNC_QUALIFIER point_type &operator[](size_t i) {
    return p.operator[](i);
  }

  GLM_FUNC_QUALIFIER point_type const &operator[](size_t i) const {
    return p.operator[](i);
  }

  GLM_FUNC_QUALIFIER typename List::const_iterator begin() const {
    return p.begin();
  }

  GLM_FUNC_QUALIFIER typename List::const_iterator cbegin() const {
    return p.cbegin();
  }

  GLM_FUNC_QUALIFIER typename List::const_iterator end() const {
    return p.end();
  }

  GLM_FUNC_QUALIFIER typename List::const_iterator cend() const {
    return p.cend();
  }
};

template<length_t L, typename T, qualifier Q>
static Polygon<L, T, Q> operator-(Polygon<L, T, Q> const &polygon) {
  Polygon<L, T, Q> p(polygon);
  for (size_t i = 0; i < p.size(); ++i)
    p[i] = operator-(p[i]);
  return p;
}

template<length_t L, typename T, qualifier Q>
static bool operator==(Polygon<L, T, Q> const &p1, Polygon<L, T, Q> const &p2) {
  if (p1.size() != p2.size())
    return false;
  for (size_t i = 0; i < p1.size(); ++i) {
    if (p1[i] != p2[i])
      return false;
  }
  return true;
}

template<length_t L, typename T, qualifier Q>
static Polygon<L, T, Q> operator+(Polygon<L, T, Q> const &polygon, vec<L, T, Q> const &offset) {
  Polygon<L, T, Q> p(polygon);
  for (size_t i = 0; i < p.size(); ++i)
    p[i] = p[i] + offset;
  return p;
}

template<length_t L, typename T, qualifier Q>
static Polygon<L, T, Q> operator-(Polygon<L, T, Q> const &polygon, vec<L, T, Q> const &offset) {
  Polygon<L, T, Q> p(polygon);
  for (size_t i = 0; i < p.size(); ++i)
    p[i] = p[i] - offset;
  return p;
}

template<length_t L, typename T, qualifier Q>
static Polygon<L, T, Q> operator*(mat<3, 3, T, Q> const &transform, Polygon<L, T, Q> const &polygon) {
  Polygon<L, T, Q> p(polygon);
  for (size_t i = 0; i < p.size(); ++i)
    p[i] = transform * p[i];
  return p;
}

template<length_t L, typename T, qualifier Q>
static Polygon<L, T, Q> operator*(mat<3, 4, T, Q> const &transform, Polygon<L, T, Q> const &polygon) {
  Polygon<L, T, Q> p(polygon);
  for (size_t i = 0; i < p.size(); ++i)
    p[i] = transform * p[i];
  return p;
}

template<length_t L, typename T, qualifier Q>
static Polygon<L, T, Q> operator*(mat<4, 3, T, Q> const &transform, Polygon<L, T, Q> const &polygon) {
  Polygon<L, T, Q> p(polygon);
  for (size_t i = 0; i < p.size(); ++i)
    p[i] = transformPos(transform, p[i]);
  return p;
}

template<length_t L, typename T, qualifier Q>
static Polygon<L, T, Q> operator*(mat<4, 4, T, Q> const &transform, Polygon<L, T, Q> const &polygon) {
  Polygon<L, T, Q> p(polygon);
  for (size_t i = 0; i < p.size(); ++i)
    p[i] = transformPos(transform, p[i]);
  return p;
}

template<length_t L, typename T, qualifier Q>
static Polygon<L, T, Q> operator*(qua<T, Q> const &transform, Polygon<L, T, Q> const &polygon) {
  Polygon<L, T, Q> p(polygon);
  for (size_t i = 0; i < p.size(); ++i)
    p[i] = transform * p[i];
  return p;
}

/// Return the number of points in the polygon
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER size_t length(Polygon<L, T, Q> const &polygon) {
  return polygon.size();
}

/// Return a vertex of this polygon: [0, length(polygon) - 1]
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> vertex(Polygon<L, T, Q> const &polygon, size_t i) {
  const size_t size = polygon.size();
  return (size == 0 || i >= size) ? vec<L, T, Q>(T(0)) : polygon[i];
}

/// Return a Segment between two adjacent vertices of the polygon.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER Segment<L, T, Q> edge(Polygon<L, T, Q> const &polygon, size_t i) {
  const size_t size = polygon.size();
  if (size == 0 || i >= size)
    return Segment<L, T, Q>(T(0));
  else if (size == 1)
    return Segment<L, T, Q>(polygon[0], polygon[0]);
  else {
    return Segment<L, T, Q>(polygon[i], polygon[(i + 1) % size]);
  }
}

/// Return a Segment between two adjacent vertices of the polygon, in the local space of the polygon.
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER Segment<2, T, Q> edge2d(Polygon<3, T, Q> const &polygon, size_t i) {
  const size_t size = polygon.size();
  if (size == 0 || i >= size)
    return Segment<2, T, Q>(T(0));
  else if (size == 1)
    return Segment<2, T, Q>(vec<2, T, Q>(T(0)), vec<2, T, Q>(T(0)));
  else {
    return Segment<2, T, Q>(mapTo2D(polygon, i), mapTo2D(polygon, (i + 1) % size));
  }
}

/// <summary>
/// Return the normal of the given edge, i.e., the vector perpendicular to the
/// plane the polygon exists in.
/// </summary>
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> edgeNormal(Polygon<L, T, Q> const &polygon, size_t idx) {
  return normalize(cross(edge(polygon, idx).dir(), normalCCW(polygon)));
}

/// Return the normal plane of the given edge.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER Plane<L, T, Q> edgePlane(Polygon<L, T, Q> const &polygon, size_t idx) {
  return Plane<L, T, Q>(edge(polygon, idx).a, edgeNormal(polygon, idx));
}

/// Return an extreme point along the polygon, i.e., the furthest point in a given direction.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> extremePoint(Polygon<L, T, Q> const &polygon, vec<L, T, Q> const &direction, T &projectionDistance) {
  projectionDistance = -std::numeric_limits<T>::infinity();

  vec<L, T, Q> mostExtreme(T(0));
  for (size_t i = 0; i < polygon.size(); ++i) {
    const T d = dot(direction, polygon[i]);
    if (d > projectionDistance) {
      projectionDistance = d;
      mostExtreme = polygon[i];
    }
  }
  return mostExtreme;
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> extremePoint(Polygon<L, T, Q> const &polygon, vec<L, T, Q> const &direction) {
  T projectionDistance(0);
  return extremePoint(polygon, direction, projectionDistance);
}

/// Project the Polygon onto the provided axis.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER void projectToAxis(Polygon<L, T, Q> const &polygon, vec<L, T, Q> const &direction, T &outMin, T &outMax) {
  outMin = dot(extremePoint(polygon, -direction), direction);
  outMax = dot(extremePoint(polygon, direction), direction);
}

/// <summary>
/// Test whether the diagonal that joins the two vertices lies inside the
/// polygon and is not intersected by edges of the polygon.
/// </summary>
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER_OUTLINE bool diagonalExists(Polygon<L, T, Q> const &polygon, size_t i, size_t j) {
  const size_t size = polygon.size();
  if (i > j) {  // Ensure "i" is the minimal index.
    const size_t tmp = i;
    i = j;
    j = tmp;
  }

  if (size < 3 || i == j)  // Degenerate if i == j.
    return false;
  else if (i >= size || j >= size)
    return false;
  else if (i + 1 == j)  // Is this Segment an edge of this polygon?
    return false;

  GLM_GEOM_ASSUME(isPlanar(polygon), false);
  const Plane<L, T, Q> polygonPlane = planeCCW(polygon);
  const Segment<L, T, Q> diag = project(polygonPlane, Segment<L, T, Q>(polygon[i], polygon[j]));

  // First check that this diagonal line is not intersected by an edge of this polygon.
  for (size_t k = 0; k < size; ++k) {
    if (!(k == i || k + 1 == i || k == j)) {
      const Segment<L, T, Q> d = project(polygonPlane, Segment<L, T, Q>(polygon[k], polygon[k + 1]));
      if (intersects(d, diag))
        return false;
    }
  }

  return isConvex(polygon);
}

/// <summary>
/// Return the diagonal that joins the two given vertices of the polygon.
/// If |i - j| == 1, then an edge of the polygon is returned.
/// </summary>
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER Segment<L, T, Q> diagonal(Polygon<L, T, Q> const &polygon, size_t i, size_t j) {
  const size_t size = polygon.size();
  const vec<L, T, Q> a = i < size ? polygon[i] : vec<L, T, Q>(T(0));
  const vec<L, T, Q> b = j < size ? polygon[j] : vec<L, T, Q>(T(0));
  return Segment<L, T, Q>(a, b);
}

/// Generates the U-vector (i.e., local space x axis) of the polygon.
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<3, T, Q> basisU(Polygon<3, T, Q> const &polygon) {
  if (polygon.size() < 2)
    return right<T, Q>();
  return normalize(polygon[1] - polygon[0]);
}

/// Generates the V-vector (i.e., local-space y axis) of the polygon.
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<3, T, Q> basisV(Polygon<3, T, Q> const &polygon) {
  if (polygon.size() < 2)
    return up<T, Q>();
  return normalize(cross(normalCCW(polygon), basisU(polygon)));
}

/// Maps the given (world) space point to the local 2D space of the polygon.
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<2, T, Q> mapTo2D(Polygon<3, T, Q> const &polygon, vec<3, T, Q> const &point) {
  const vec<3, T, Q> bu = basisU(polygon);
  const vec<3, T, Q> bv = basisV(polygon);
  const vec<3, T, Q> pt = point - ((polygon.size() == 0) ? vec<3, T, Q>(T(0)) : polygon[0]);
  return vec<2, T, Q>(dot(pt, bu), dot(pt, bv));
}

/// Map the given vertex to the local 2D space of the polygon.
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<2, T, Q> mapTo2D(Polygon<3, T, Q> const &polygon, size_t i) {
  return (i < polygon.size()) ? mapTo2D(polygon, polygon[i]) : vec<2, T, Q>(T(0));
}

/// Map the given local 2D space coordinate to a 3D point world space coordinate.
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<3, T, Q> mapFrom2D(Polygon<3, T, Q> const &polygon, vec<2, T, Q> const &point) {
  if (polygon.size() == 0)
    return vec<3, T, Q>(T(0));
  return polygon[0] + point.x * basisU(polygon) + point.y * basisV(polygon);
}

/// Return the surface area of the polygon.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T area(Polygon<L, T, Q> const &polygon) {
  const size_t size = polygon.size();
  if (size <= 2)
    return T(0);

  GLM_GEOM_ASSUME(isPlanar(polygon), T(0));
  vec<L, T, Q> area(0);
  size_t i = size - 1;
  for (size_t j = 0; j < size; ++j) {
    area += cross(polygon[i], polygon[j]);
    i = j;
  }
  return abs(dot(normalCCW(polygon), area)) * T(0.5);
}

/// Return the total edge length of the polygon.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T perimeter(Polygon<L, T, Q> const &polygon) {
  T perimeter = T(0);
  for (size_t i = 0; i < polygon.size(); ++i)
    perimeter += length(edge(polygon, i));
  return perimeter;
}

/// Return the center of mass of the polygon.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> centroid(Polygon<L, T, Q> const &polygon) {
  const size_t size = polygon.size();
  if (size == 0)
    return vec<L, T, Q>(T(0));

  vec<L, T, Q> centroid(T(0));
  for (const typename Polygon<L, T, Q>::point_type &p : polygon)
    centroid += p;
  return centroid / static_cast<T>(size);
}

/// Test if the polygon is planar, i.e., all of its vertices lie on the same plane.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER_OUTLINE bool isPlanar(Polygon<L, T, Q> const &polygon, T epsSq = epsilon2<T>()) {
  const size_t size = polygon.size();
  if (size == 0)
    return false;
  else if (size <= 3)
    return true;

  const vec<L, T, Q> normal = cross(polygon[1] - polygon[0], polygon[2] - polygon[0]);
  const T lenSq = length2(normal);
  for (size_t i = 3; i < size; ++i) {
    const T d = dot(normal, polygon[i] - polygon[0]);
    if (d * d > epsSq * lenSq) {
      return false;
    }
  }
  return true;
}

/// Test if the polygon is simple, i.e., no two non-consecutive edges have a point in common.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER_OUTLINE bool isSimple(Polygon<L, T, Q> const &polygon) {
  GLM_GEOM_ASSUME(isPlanar(polygon), false);

  const size_t size = polygon.size();
  const Plane<L, T, Q> plane = planeCCW(polygon);
  for (size_t i = 0; i < size; ++i) {
    const Segment<L, T, Q> si = project(plane, edge(polygon, i));
    for (size_t j = i + 2; j < size; ++j) {
      if (i != 0 || j != (size - 1)) {
        const Segment<L, T, Q> sj = project(plane, edge(polygon, j));
        if (intersects(si, sj)) {
          return false;
        }
      }
    }
  }
  return true;
}

/// Test if the polygon is null, i.e., has no vertices.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool isNull(Polygon<L, T, Q> const &polygon) {
  return polygon.size() == 0;
}

/// Test if all vertices of the polygon are finite
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool isfinite(Polygon<L, T, Q> const &polygon) {
  if (polygon.size() == 0)
    return true;

  for (const typename Polygon<L, T, Q>::point_type &p : polygon) {
    if (!all(isfinite(p)))
      return false;
  }
  return true;
}

/// <summary>
/// Return true if the polygon is degenerate:
///   1. It has two-or-less vertices;
///   2. its surface area is less or equal than a given epsilon
/// </summary>
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool isDegenerate(Polygon<L, T, Q> const &polygon, T eps = epsilon<T>()) {
  return polygon.size() < 3 || area(polygon) <= eps;
}

/// @private
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool orientedCCW(vec<2, T, Q> const &a, vec<2, T, Q> const &b, vec<2, T, Q> const &c) {
  return (a.x - c.x) * (b.y - c.y) - (a.y - c.y) * (b.x - c.x) >= T(0);
}

/// <summary>
/// Test whether the polygon is convex: For each pair of points inside the
/// polygon, the segment joining those points is also completely inside the
/// polygon.
/// </summary>
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER_OUTLINE bool isConvex(Polygon<3, T, Q> const &polygon) {
  const size_t size = polygon.size();
  if (size == 0)
    return false;
  else if (size <= 3)
    return true;

  GLM_GEOM_ASSUME(isPlanar(polygon), false);
  size_t i = size - 2;
  size_t j = size - 1;
  size_t k = 0;
  while (k < size) {
    const vec<2, T, Q> a = mapTo2D(polygon, i);
    const vec<2, T, Q> b = mapTo2D(polygon, j);
    const vec<2, T, Q> c = mapTo2D(polygon, k);
    if (!orientedCCW(a, b, c))
      return false;

    i = j;
    j = k;
    ++k;
  }
  return true;
}

/// Return a point on the perimeter of this polygon
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER_OUTLINE vec<L, T, Q> pointOnEdge(Polygon<L, T, Q> const &polygon, T dist) {
  assert(dist >= zero<T>() && dist <= one<T>());  // Only defined in [0, 1]
  const size_t size = polygon.size();
  if (size == 0)
    return vec<L, T, Q>(T(0));
  else if (size < 2)
    return polygon[0];

  T d = perimeter(polygon) * (dist - floor(dist));
  for (size_t i = 0; i < size; ++i) {
    const Segment<L, T, Q> e = edge(polygon, i);
    const T len = length(e);
    if (detail::approx_zero(len))
      return vec<L, T, Q>(T(0));  // degenerate polygon

    if (d <= len)
      return getPoint(e, d / len);
    d -= len;
  }
  return polygon[0];  // Throw an assertion failure, should not reach this.
}

/// <summary>
/// Return the Plane the Polygon is contained in.
///
/// The normal of the plane points to the direction from which the vertices wind
/// in counter-clockwise order.
/// </summary>
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER_OUTLINE Plane<3, T, Q> planeCCW(Polygon<3, T, Q> const &polygon) {
  const size_t size = polygon.size();
  const vec<3, T, Q> hint = forward<T, Q>();
  const vec<3, T, Q> hint2 = up<T, Q>();
  switch (size) {
    case 0: return Plane<3, T, Q>(T(0));
    case 1: return planeFrom(polygon[0], up<T, Q>());
    case 2: {
      const vec<3, T, Q> d = normalize(polygon[1] - polygon[0]);
      return planeFrom(Line<3, T, Q>(polygon[0], d), perpendicular(d, hint, hint2));
    }
    case 3: return planeFrom(polygon[0], polygon[1], polygon[2]);
    default: {
      break;
    }
  }

  for (size_t i = 0; i < size - 2; ++i) {
    for (size_t j = i + 1; j < size - 1; ++j) {
      const vec<3, T, Q> pij = polygon[j] - polygon[i];
      for (size_t k = j + 1; k < size; ++k) {
        vec<3, T, Q> normal = cross(pij, polygon[k] - polygon[i]);
        const T lenSq = length2(normal);
        if (lenSq > epsilon<T>()) {
          normal /= sqrt(lenSq);
          return Plane<3, T, Q>(normal, dot(normal, polygon[i]));
        }
      }
    }
  }

  // Collinear points cannot form a plane.
  const vec<3, T, Q> d = normalize(polygon[1] - polygon[0]);
  return planeFrom(Line<3, T, Q>(polygon[0], d), perpendicular(d, hint, hint2));
}

/// Return the normal of the polygon in the counter-clockwise direction.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> normalCCW(Polygon<L, T, Q> const &polygon) {
  return planeCCW(polygon).normal;
}

/// Return the (clockwise, i.e., the normal points in the clockwise direction) plane this polygon is contained in.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER Plane<3, T, Q> planeCW(Polygon<L, T, Q> const &polygon) {
  return reverseNormal(planeCCW(polygon));
}

/// Return the normal of the polygon in the clockwise direction.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> normalCW(Polygon<L, T, Q> const &polygon) {
  return planeCW(polygon).normal;
}

/// Return the smallest AABB that encloses the polygon.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER AABB<L, T, Q> minimalEnclosingAABB(Polygon<L, T, Q> const &polygon) {
  AABB<L, T, Q> aabb(T(0));
  if (polygon.size() > 0) {
    aabb.setNegativeInfinity();
    for (const typename Polygon<L, T, Q>::point_type &p : polygon) {
      aabb.enclose(p);
    }
  }
  return aabb;
}

#if 0
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER OBB<T, Q> minimalEnclosingOBB(Polygon<3, T, Q> const &polygon) { }
#endif

/* Test if the given object (worldspace) is fully contained inside the polygon */

/// @private epsilon value for containment tests.
template<typename genType>
GLM_FUNC_QUALIFIER GLM_CONSTEXPR genType thickness_eps() {
  return intersect_eps<genType>();
}

template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER_OUTLINE bool contains(Polygon<3, T, Q> const &polygon, vec<3, T, Q> const &point, ePolyContains type, T thickness = thickness_eps<T>()) {
  const size_t size = polygon.size();
  if (size < 3)
    return false;

  const vec<3, T, Q> bu = basisU(polygon);
  const vec<3, T, Q> bv = basisV(polygon);
  if (!areOrthonormal(bu, bv, epsilon<T>()))
    return false;

  // Check if the point is within the plane of the polygon
  const T thick_eps = thickness + thickness_eps<T>();
  const T p_diff = dot(cross(bu, bv), polygon[0] - point);
  bool contains = true;
  switch (type) {
    case ePolyContains::Positive: contains = p_diff >= T(0) && p_diff <= thick_eps; break;
    case ePolyContains::Negative: contains = p_diff <= T(0) && -p_diff <= thick_eps; break;
    case ePolyContains::Unidirectional: {
      contains = abs(p_diff) <= T(0.5) * thick_eps;
      break;
    }
  }

  if (!contains)
    return false;

  // Crossings Test
  const T eps = epsilon<T>();
  vec<3, T, Q> vt = polygon.back() - point;
  vec<2, T, Q> p0 = vec<2, T, Q>(dot(vt, bu), dot(vt, bv));
  if (abs(p0.y) < eps)
    p0.y = -eps;

  size_t numIntersections = 0;
  for (size_t i = 0; i < size; ++i) {
    vt = polygon[i] - point;

    vec<2, T, Q> p1 = vec<2, T, Q>(dot(vt, bu), dot(vt, bv));
    if (abs(p1.y) < eps)
      p1.y = -eps;

    if (p0.y * p1.y < T(0)) {
      if (min(p0.x, p1.x) > T(0))
        ++numIntersections;
      else if (max(p0.x, p1.x) > T(0)) {
        const vec<2, T, Q> delta = p1 - p0;
        if (!detail::exactly_zero(delta.y)) {
          const T t = -p0.y / delta.y;
          const T x = p0.x + t * delta.x;
          if (t >= T(0) && t <= T(1) && x > T(0)) {
            ++numIntersections;
          }
        }
      }
    }
    p0 = p1;
  }

  return (numIntersections % 2) == 1;
}

/// Polygon/Point containment test, expanding 0.5 * thickness in each halfspace
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool contains(Polygon<3, T, Q> const &polygon, vec<3, T, Q> const &point, T thickness = thickness_eps<T>()) {
  return contains(polygon, point, ePolyContains::Unidirectional, thickness);
}

/// Polygon/Point containment test, expanding thickness in the positive halfspace
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool containsPositive(Polygon<3, T, Q> const &polygon, vec<3, T, Q> const &point, T thickness = thickness_eps<T>()) {
  return contains(polygon, point, ePolyContains::Positive, thickness);
}

/// Polygon/Point containment test, expanding thickness in the negative halfspace
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool containsNegative(Polygon<3, T, Q> const &polygon, vec<3, T, Q> const &point, T thickness = thickness_eps<T>()) {
  return contains(polygon, point, ePolyContains::Negative, thickness);
}

template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER_OUTLINE bool contains2D(Polygon<3, T, Q> const &polygon, Segment<3, T, Q> const &localSegment) {
  if (polygon.size() < 3)
    return false;

  const vec<3, T, Q> bu = basisU(polygon);
  const vec<3, T, Q> bv = basisV(polygon);

  Segment<3, T, Q> edge(T(0));
  edge.a = vec<3, T, Q>(dot(polygon.back(), bu), dot(polygon.back(), bv), T(0));
  for (const typename Polygon<3, T, Q>::point_type &p : polygon) {
    edge.b = vec<3, T, Q>(dot(p, bu), dot(p, bv), T(0));
    if (intersects(edge, localSegment))
      return false;

    edge.a = edge.b;
  }

  // segment is fully inside or outside: determine which
  return contains(polygon, mapFrom2D(polygon, vec<2, T, Q>(localSegment.a.x, localSegment.a.y)));
}

// Polygon/Polygon containment test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool contains(Polygon<L, T, Q> const &polygon, Polygon<L, T, Q> const &worldSpacePolygon, T thickness = thickness_eps<T>()) {
  if (polygon.size() == 0)
    return false;

  for (const typename Polygon<L, T, Q>::point_type &p : worldSpacePolygon) {
    if (!contains(polygon, p, thickness)) {
      return false;
    }
  }
  return true;
}

// Polygon/Segment containment test
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool contains(Polygon<3, T, Q> const &polygon, Segment<3, T, Q> const &worldSpaceSegment, T thickness = thickness_eps<T>()) {
  const size_t size = polygon.size();
  if (size < 3)
    return false;

  const Plane<3, T, Q> plane = planeCCW(polygon);
  if (distance(plane, worldSpaceSegment.a) > thickness || distance(plane, worldSpaceSegment.b) > thickness)
    return false;

  const Segment<3, T, Q> l = project(plane, worldSpaceSegment);
  if (!contains(polygon, l.a) || !contains(polygon, l.b))
    return false;

  for (size_t i = 0; i < size; ++i) {
    if (intersects(project(plane, edge(polygon, i)), l)) {
      return false;
    }
  }
  return true;
}

// Polygon/Triangle containment test
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool contains(Polygon<3, T, Q> const &polygon, Triangle<3, T, Q> const &worldSpaceTriangle, T thickness = thickness_eps<T>()) {
  return contains(polygon, edge(worldSpaceTriangle, 0), thickness)
         && contains(polygon, edge(worldSpaceTriangle, 1), thickness)
         && contains(polygon, edge(worldSpaceTriangle, 2), thickness);
}

// Polygon/Capsule containment test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool contains(Polygon<L, T, Q> const &polygon, Capsule<L, T, Q> const &capsule) {
  return contains(capsule, polygon);
}

#if 0
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool contains(Polygon<3, T, Q> const &polygon, Disc<3, T, Q> const &disc) { }
#endif

/* Return the closest point between an AABB and the given object */

#if 0
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Polygon<L, T, Q> const &polygon, vec<L, T, Q> const &point) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Polygon<L, T, Q> const &polygon, Line<L, T, Q> const &line, T &d) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Polygon<L, T, Q> const &polygon, Segment<L, T, Q> const &segment, T &d) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Polygon<L, T, Q> const &polygon, Ray<L, T, Q> const &ray, T &d) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Polygon<L, T, Q> const &polygon, Sphere<L, T, Q> const &sphere) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Polygon<L, T, Q> const &polygon, Capsule<L, T, Q> const &capsule) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Polygon<L, T, Q> const &polygon, Disc<L, T, Q> const &disc) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Polygon<L, T, Q> const &polygon, Plane<L, T, Q> const &plane) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Polygon<L, T, Q> const &polygon, Triangle<L, T, Q> const &triangle) { }
#endif

/* Test whether the polygon and the given object intersect */

/// Polygon/Line intersection test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Polygon<L, T, Q> const &polygon, Line<L, T, Q> const &line) {
  T d(0);
  const Plane<L, T, Q> plane = planeCCW(polygon);
  return intersects(plane, line, d) ? contains(polygon, getPoint(line, d)) : false;
}

/// Polygon/Ray intersection test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Polygon<L, T, Q> const &polygon, Ray<L, T, Q> const &ray) {
  T d(0);
  const Plane<L, T, Q> plane = planeCCW(polygon);
  return intersects(plane, ray, d) ? contains(polygon, getPoint(ray, d)) : false;
}

template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER_OUTLINE bool intersects2D(Polygon<3, T, Q> const &polygon, Segment<3, T, Q> const &localSpaceSegment) {
  const size_t size = polygon.size();
  if (size < 3)
    return false;

  const vec<3, T, Q> bu = basisU(polygon);
  const vec<3, T, Q> bv = basisV(polygon);

  Segment<3, T, Q> edge(T(0));
  edge.a = vec<3, T, Q>(dot(polygon.back(), bu), dot(polygon.back(), bv), T(0));
  for (size_t i = 0; i < size; ++i) {
    edge.b = vec<3, T, Q>(dot(polygon[i], bu), dot(polygon[i], bv), 0);
    if (intersects(edge, localSpaceSegment))
      return true;
    edge.a = edge.b;
  }

  // segment is fully inside or outside: determine which
  return contains(polygon, mapFrom2D(polygon, vec<2, T, Q>(localSpaceSegment.a.x, localSpaceSegment.a.y)));
}

/// Polygon/Segment intersection test
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER_OUTLINE bool intersects(Polygon<3, T, Q> const &polygon, Segment<3, T, Q> const &line) {
  const Plane<3, T, Q> plane = planeCCW(polygon);
  const T denom = dot(plane.normal, line.b - line.a);  // Compute line-plane intersection
  if (abs(denom) < intersect_eps<T>()) {  // They are planar
    const vec<2, T, Q> a = mapTo2D(polygon, line.a);
    const vec<2, T, Q> b = mapTo2D(polygon, line.b);
    const Segment<3, T, Q> segment(vec<3, T, Q>(a.x, a.y, T(0)), vec<3, T, Q>(b.x, b.y, T(0)));
    return intersects2D(polygon, segment);
  }

  const T t = (plane.d - dot(plane.normal, line.a)) / denom;
  return (t >= T(0) && t <= T(1)) ? contains(polygon, getPoint(line, t)) : false;
}

/// Polygon/Plane intersection test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Polygon<L, T, Q> const &polygon, Plane<L, T, Q> const &plane) {
  if (polygon.size() == 0)
    return false;

  // Project this polygon onto the plane. If there are points on both sides of
  // the plane, then the polygon intersects it.
  T minD = std::numeric_limits<T>::infinity();
  T maxD = -std::numeric_limits<T>::infinity();
  for (const typename Polygon<L, T, Q>::point_type &p : polygon) {
    const T d = signedDistance(plane, p);
    minD = min(minD, d);
    maxD = max(maxD, d);
  }

  return minD <= intersect_eps<T>() && maxD >= -intersect_eps<T>();
}

#if 0
template<length_t L, typename T, qualifier Q, typename Object>
static bool intersectsObject(Polygon<L, T, Q> const &polygon, Object const &obj, T thickness = thickness_eps<T>()) {
  return false;
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Polygon<L, T, Q> const &polygon, Polygon<L, T, Q> const &other) {
  return intersectsObject(polygon, other);
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Polygon<L, T, Q> const &polygon, Triangle<L, T, Q> const &triangle) {
  return intersectsObject(polygon, triangle);
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Polygon<L, T, Q> const &polygon, Sphere<L, T, Q> const &sphere) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Polygon<L, T, Q> const &polygon, Capsule<L, T, Q> const &capsule) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Polygon<L, T, Q> const &polygon, Disc<L, T, Q> const &disc) { }
#endif

/// @}
GLM_GEOM_END_NAMESPACE
#if GLM_GEOM_TOSTRING
namespace glm {
  namespace detail {
    template<length_t L, typename T, qualifier Q>
    struct compute_to_string<Polygon<L, T, Q>> {
      GLM_GEOM_QUALIFIER std::string call(Polygon<L, T, Q> const &polygon) {
        ((void)polygon);
        return std::string("Polygon");
      }
    };
  }
}
#endif
