/// @ref geom_segment
/// @file geom/segment.hpp
///
/// @defgroup geom_segment Segment
/// @ingroup geom
///

#pragma once

#include "setup.hpp"
#include "line.hpp"
#include "triangle.hpp"
#if GLM_MESSAGES == GLM_ENABLE && !defined(GLM_EXT_INCLUDED)
  #pragma message("GLM: GLM_EXT_GEOM_segment extension included")
#endif

GLM_GEOM_BEGIN_NAMESPACE
/// @addtogroup geom_segment
/// @{

/// <summary>
/// A line represented by (1 - T) * A + T * B, where A and B are the finite
/// endpoints of the segment and 0 <= T <= 1.
/// </summary>
template<length_t L, typename T, qualifier Q>
struct Segment {

  // -- Implementation detail --

  typedef T value_type;
  typedef Segment<L, T, Q> type;
  typedef vec<L, T, Q> point_type;

  // -- Data --

  point_type a;  // Starting point of this segment
  point_type b;  // End point of this segment

#if GLM_CONFIG_DEFAULTED_DEFAULT_CTOR == GLM_ENABLE
  Segment() GLM_DEFAULT_CTOR;
#else
  Segment()
  #if GLM_CONFIG_CTOR_INIT != GLM_CTOR_INIT_DISABLE
    : a(T(0)), b(T(0))
  #endif
  {
  }
#endif

  Segment(value_type scalar)
    : a(scalar), b(scalar) {
  }

  Segment(point_type const &begin, point_type const &end)
    : a(begin), b(end) {
  }

  Segment(Segment<L, T, Q> const &segment)
    : a(segment.a), b(segment.b) {
  }

  Segment<L, T, Q> &operator=(Segment<L, T, Q> const &segment) {
    a = segment.a;
    b = segment.b;
    return *this;
  }

  GLM_FUNC_QUALIFIER point_type dir() const {
    return normalize(b - a);
  }

  GLM_FUNC_QUALIFIER point_type dir2() const {
    return b - a;
  }
};

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER Line<L, T, Q> toLine(Segment<L, T, Q> const &segment) {
  return Line<L, T, Q>(segment.a, segment.dir());
}

template<length_t L, typename T, qualifier Q>
static Segment<L, T, Q> operator-(Segment<L, T, Q> const &segment) {
  return Segment<L, T, Q>(-segment.b, -segment.a);
}

template<length_t L, typename T, qualifier Q>
static bool operator==(Segment<L, T, Q> const &l1, Segment<L, T, Q> const &l2) {
  return l1.a == l2.a && l1.b == l2.b;
}

template<length_t L, typename T, qualifier Q>
static Segment<L, T, Q> operator+(Segment<L, T, Q> const &segment, vec<L, T, Q> const &offset) {
  return Segment<L, T, Q>(segment.a + offset, segment.b + offset);
}

template<length_t L, typename T, qualifier Q>
static Segment<L, T, Q> operator-(Segment<L, T, Q> const &segment, vec<L, T, Q> const &offset) {
  return Segment<L, T, Q>(segment.a - offset, segment.b - offset);
}

template<typename T, qualifier Q>
static Segment<3, T, Q> operator*(mat<3, 3, T, Q> const &m, Segment<3, T, Q> const &line) {
  return Segment<3, T, Q>(m * line.a, m * line.b);
}

template<typename T, qualifier Q>
static Segment<3, T, Q> operator*(mat<3, 4, T, Q> const &m, Segment<3, T, Q> const &line) {
  return Segment<3, T, Q>(m * line.a, m * line.b);
}

template<typename T, qualifier Q>
static Segment<3, T, Q> operator*(mat<4, 3, T, Q> const &m, Segment<3, T, Q> const &line) {
  return Segment<3, T, Q>(transformPos(m, line.a), transformPos(m, line.b));
}

template<typename T, qualifier Q>
static Segment<3, T, Q> operator*(mat<4, 4, T, Q> const &m, Segment<3, T, Q> const &line) {
  return Segment<3, T, Q>(transformPos(m, line.a), transformPos(m, line.b));
}

template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER Segment<3, T, Q> operator*(qua<T, Q> const &q, Segment<3, T, Q> const &line) {
  return Segment<3, T, Q>(q * line.a, q * line.b);
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool equal(Segment<L, T, Q> const &x, Segment<L, T, Q> const &y, T eps = epsilon<T>()) {
  return all(equal(x.a, y.a, eps)) && all(equal(x.b, y.b, eps));
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool equal(Segment<L, T, Q> const &x, Segment<L, T, Q> const &y, vec<L, T, Q> const &eps) {
  return all(equal(x.a, y.a, eps)) && all(equal(x.b, y.b, eps));
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool equal(Segment<L, T, Q> const &x, Segment<L, T, Q> const &y, int MaxULPs) {
  return all(equal(x.a, y.a, MaxULPs)) && all(equal(x.b, y.b, MaxULPs));
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool equal(Segment<L, T, Q> const &x, Segment<L, T, Q> const &y, vec<L, int, Q> const &MaxULPs) {
  return all(equal(x.a, y.a, MaxULPs)) && all(equal(x.b, y.b, MaxULPs));
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool notEqual(Segment<L, T, Q> const &x, Segment<L, T, Q> const &y, T eps = epsilon<T>()) {
  return any(notEqual(x.a, y.a, eps)) || any(notEqual(x.b, y.b, eps));
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool notEqual(Segment<L, T, Q> const &x, Segment<L, T, Q> const &y, vec<L, T, Q> const &eps) {
  return any(notEqual(x.a, y.a, eps)) || any(notEqual(x.b, y.b, eps));
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool notEqual(Segment<L, T, Q> const &x, Segment<L, T, Q> const &y, int MaxULPs) {
  return any(notEqual(x.a, y.a, MaxULPs)) || any(notEqual(x.b, y.b, MaxULPs));
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool notEqual(Segment<L, T, Q> const &x, Segment<L, T, Q> const &y, vec<L, int, Q> const &MaxULPs) {
  return any(notEqual(x.a, y.a, MaxULPs)) || any(notEqual(x.b, y.b, MaxULPs));
}

/// Return the length of the segment
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T length(Segment<L, T, Q> const &segment) {
  return distance(segment.a, segment.b);
}

/// Return the squared length of the segment
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T length2(Segment<L, T, Q> const &segment) {
  return distance2(segment.a, segment.b);
}

/// Test if any component of the segment is infinite
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool isinf(Segment<L, T, Q> const &segment) {
  return any(isinf(segment.a)) || any(isinf(segment.b));
}

/// Test if any component of the segment is NaN
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool isnan(Segment<L, T, Q> const &segment) {
  return any(isnan(segment.a)) || any(isnan(segment.b));
}

/// Test if all components of the segment are finite
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool isfinite(Segment<L, T, Q> const &segment) {
  return all(isfinite(segment.a)) && all(isfinite(segment.b));
}

/// Get a point along the line at a given distance (parametric point)
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> getPoint(Segment<L, T, Q> const &segment, T d) {
  return (T(1) - d) * segment.a + d * segment.b;
}

/// Return the center point of this segment; getPoint(line, T(0.5))
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> centerPoint(Segment<L, T, Q> const &segment) {
  return (segment.a + segment.b) * T(0.5);
}

/// Reverse the direction of the segment.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER Segment<L, T, Q> reverse(Segment<L, T, Q> const &segment) {
  return Segment<L, T, Q>(segment.b, segment.a);
}

/// Return the (normalized) direction vector that points from line.a to line.b
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> dir(Segment<L, T, Q> const &segment) {
  return normalize(segment.dir2());
}

/// Return an extreme point along the segment, i.e., the furthest point in a given direction.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> extremePoint(Segment<L, T, Q> const &segment, vec<L, T, Q> const &direction) {
  return dot(direction, segment.dir2()) >= T(0) ? segment.b : segment.a;
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> extremePoint(Segment<L, T, Q> const &segment, vec<L, T, Q> const &direction, T &projectionDistance) {
  vec<L, T, Q> point = extremePoint(segment, direction);
  projectionDistance = dot(point, direction);
  return point;
}

/// Project the Segment onto the given axis (direction), i.e., collapse the segment onto an axis.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER void projectToAxis(Segment<L, T, Q> const &segment, vec<L, T, Q> const &direction, T &outMin, T &outMax) {
  outMin = dot(direction, segment.a);
  outMax = dot(direction, segment.b);
  if (outMax < outMin) {
    const T temp = outMin;
    outMin = outMax;
    outMax = temp;
  }
}

/* Return the closest point on this segment to the given object */

/// Return the closest point between the Segment and given point; including parametric segment coordinate.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Segment<L, T, Q> const &segment, vec<L, T, Q> const &point, T &d) {
  const vec<L, T, Q> dir = segment.dir2();
  d = glm::clamp(dot(point - segment.a, dir) / length2(dir), T(0), T(1));
  return segment.a + d * dir;
}

/// <summary>
/// Return the closest point between the Segment and given Ray; including the
/// parametric coordinates for segment (d1) and ray (d2)
/// </summary>
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Segment<L, T, Q> const &segment, Ray<L, T, Q> const &ray, T &d1, T &d2) {
  closestPoint(ray, segment, d2, d1);
  return getPoint(segment, d1);
}

/// <summary>
/// Return the closest point between the Segment and given line; including the
/// parametric coordinates for segment (d1) and line (d2)
/// </summary>
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Segment<L, T, Q> const &segment, Line<L, T, Q> const &line, T &d1, T &d2) {
  closestPointLineLine(line.pos, line.dir, segment.a, segment.dir2(), d2, d1);
  if (d1 < T(0)) {
    d1 = T(0);
    closestPoint(line, segment.a, d2);
    return segment.a;
  }
  else if (d1 > T(1)) {
    d1 = T(1);
    closestPoint(line, segment.b, d2);
    return segment.b;
  }
  return getPoint(segment, d1);
}

/// Return the closest point between two Segments; including parametric coordinates both segments
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER_OUTLINE vec<L, T, Q> closestPoint(Segment<L, T, Q> const &segment, Segment<L, T, Q> const &other, T &d1, T &d2) {
  closestPointLineLine(segment.a, segment.dir2(), other.a, other.dir2(), d1, d2);

  if (d1 >= T(0) && d1 <= T(1) && d2 >= T(0) && d2 <= T(1))
    return segment.a + d1 * segment.dir2();
  else if (d1 >= T(0) && d1 <= T(1)) {
    const vec<L, T, Q> p = (d2 < T(0)) ? other.a : other.b;
    d2 = (d2 < T(0)) ? T(0) : T(1);
    return closestPoint(segment, p, d1);
  }
  else {
    const vec<L, T, Q> p = (d1 < T(0)) ? segment.a : segment.b;
    const vec<L, T, Q> p2 = (d2 < T(0)) ? other.a : other.b;
    d1 = (d1 < T(0)) ? T(0) : T(1);
    d2 = (d2 < T(0)) ? T(0) : T(1);

    T dt(0), dt2(0);
    const vec<L, T, Q> pt = closestPoint(segment, p2, dt);
    const vec<L, T, Q> pt2 = closestPoint(other, p, dt2);
    if (distance2(pt, p2) <= distance2(pt2, p)) {
      d1 = dt;
      return pt;
    }
    else {
      d2 = dt2;
      return p;
    }
  }
}

/// @private see triangle.hpp
template<length_t L, typename T, qualifier Q>
GLM_GEOM_DECL vec<L, T, Q> closestPointTriangleSegment(Triangle<L, T, Q> const &t, Segment<L, T, Q> const &segment, T &outU, T &outV, T &outD);

/// Segment/Triangle intersection test, includes the parametric and barycentric intersection coordinates.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Segment<L, T, Q> const &segment, Triangle<L, T, Q> const &triangle, T &d, T &u, T &v) {
  closestPointTriangleSegment<L, T, Q, Segment<L, T, Q>>(triangle, segment, u, v, d);
  return getPoint(segment, d);
}

/// Return the closest point between the Segment and given point
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Segment<L, T, Q> const &segment, vec<L, T, Q> const &point) {
  T d(0);
  return closestPoint(segment, point, d);
}

/// Return the closest point between the Segment and given Ray
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Segment<L, T, Q> const &segment, Ray<L, T, Q> const &ray) {
  T d1(0), d2(0);
  return closestPoint(segment, ray, d1, d2);
}

/// Return the closest point between the Segment and given Line
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Segment<L, T, Q> const &segment, Line<L, T, Q> const &line) {
  T d1(0), d2(0);
  return closestPoint(segment, line, d1, d2);
}

/// Return the closest point between two Segments
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Segment<L, T, Q> const &segment, Segment<L, T, Q> const &other) {
  T d1(0), d2(0);
  return closestPoint(segment, other, d1, d2);
}

/// Return the closest point between the Segment and given Triangle
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<3, T, Q> closestPoint(Segment<L, T, Q> const &segment, Triangle<L, T, Q> const &triangle) {
  T u, v, d(0);
  return closestPoint(segment, triangle, d, u, v);
}

/* Test if the given object is fully contained on the segment */

/// Return true if the point is on the segment.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool contains(Segment<L, T, Q> const &segment, vec<L, T, Q> const &point, T distanceThreshold = contains_eps<T>()) {
  T d(0);
  return distance2(closestPoint(segment, point, d), point) <= distanceThreshold;
}

/// Return true if 'rhs' lays on the segment.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool contains(Segment<L, T, Q> const &segment, Segment<L, T, Q> const &rhs, T distanceThreshold = contains_eps<T>()) {
  return contains(segment, rhs.a, distanceThreshold) && contains(segment, rhs.b, distanceThreshold);
}

/* Return the distance between the segment and the given object. */

/// Return the distance between the Segment and given point; including parametric segment coordinate.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Segment<L, T, Q> const &segment, vec<L, T, Q> const &point, T &d) {
  return distance(closestPoint(segment, point, d), point);
}

/// <summary>
/// Return the closest point between the Segment and given Ray; including the
/// parametric coordinates for segment (d1) and ray (d2)
/// </summary>
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Segment<L, T, Q> const &segment, Ray<L, T, Q> const &ray, T &d1, T &d2) {
  const vec<L, T, Q> point = closestPoint(segment, ray, d1, d2);
  return distance(point, getPoint(ray, d2));
}

/// <summary>
/// Return the closest point between the Segment and given Line; including the
/// parametric coordinates for segment (d1) and line (d2)
/// </summary>
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Segment<L, T, Q> const &segment, Line<L, T, Q> const &line, T &d1, T &d2) {
  const vec<L, T, Q> point = closestPoint(segment, line, d1, d2);
  return distance(point, getPoint(line, d2));
}

/// Return the closest point between two Segments; including the parametric coordinates.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Segment<L, T, Q> const &segment, Segment<L, T, Q> const &other, T &d1, T &d2) {
  const vec<L, T, Q> point = closestPoint(segment, other, d1, d2);
  return distance(point, getPoint(other, d2));
}

/// Return the distance between the Segment and given point
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Segment<L, T, Q> const &segment, vec<L, T, Q> const &point) {
  T d(0);
  return distance(segment, point, d);
}

/// Return the distance between the Segment and given Ray
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Segment<L, T, Q> const &segment, Ray<L, T, Q> const &ray) {
  T d1(0), d2(0);
  return distance(segment, ray, d1, d2);
}

/// Return the distance between the Segment and given Line
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Segment<L, T, Q> const &segment, Line<L, T, Q> const &line) {
  T d1(0), d2(0);
  return distance(segment, line, d1, d2);
}

/// Return the distance between two Segments
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Segment<L, T, Q> const &segment, Segment<L, T, Q> const &other) {
  T d1(0), d2(0);
  return distance(segment, other, d1, d2);
}

/// Return the squared distance between the Segment and given point; including parametric segment coordinate.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance2(Segment<L, T, Q> const &segment, vec<L, T, Q> const &point, T &d) {
  return distance2(closestPoint(segment, point, d), point);
}

/// Return the squared distance between two Segments; including the parametric coordinates.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance2(Segment<L, T, Q> const &segment, Segment<L, T, Q> const &other, T &d1, T &d2) {
  const vec<L, T, Q> point = closestPoint(segment, other, d1, d2);
  return distance2(point, getPoint(other, d2));
}

/// Return the squared distance between the Segment and given point
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance2(Segment<L, T, Q> const &segment, vec<L, T, Q> const &point) {
  T d(0);
  return distance2(segment, point, d);
}

/// Return the distance between two Segments
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance2(Segment<L, T, Q> const &segment, Segment<L, T, Q> const &other) {
  T d1(0), d2(0);
  return distance2(segment, other, d1, d2);
}

/// Return the distance between the Segment and given Sphere
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Segment<L, T, Q> const &segment, Sphere<L, T, Q> const &sphere) {
  T d(0);
  return max(T(0), distance(segment, sphere.pos, d) - sphere.r);
}

/// Return the distance between the Segment and given Capsule
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Segment<L, T, Q> const &segment, Capsule<L, T, Q> const &capsule) {
  return max(T(0), distance(segment, capsule.l) - capsule.r);
}

/// Return the distance between the Segment and given Plane
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Segment<L, T, Q> const &segment, Plane<L, T, Q> const &plane) {
  const T aDist = signedDistance(plane, segment.a);
  const T bDist = signedDistance(plane, segment.b);
  if (aDist * bDist <= T(0))
    return T(0);  // Was an intersection
  return min(abs(aDist), abs(bDist));
}

/* Test whether the segment and the given object intersect */

// Segment/AABB intersection test, includes parametric intersection coordinates
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Segment<L, T, Q> const &segment, AABB<L, T, Q> const &aabb, T &d1, T &d2) {
  return intersects(aabb, segment, d1, d2);
}

// Segment/OBB intersection test, includes parametric intersection coordinates
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Segment<3, T, Q> const &segment, OBB<T, Q> const &obb, T &dNear, T &dFar) {
  return intersects(obb, segment, dNear, dFar);
}

// Segment/AABB intersection test, includes parametric intersection coordinates
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER int intersects(Segment<L, T, Q> const &segment, Sphere<L, T, Q> const &sphere, T &d1, T &d2) {
  return intersects(sphere, segment, d1, d2);
}

// Segment/Capsule intersection test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Segment<L, T, Q> const &segment, Capsule<L, T, Q> const &capsule) {
  return intersects(capsule, segment);
}

/// Segment/Plane intersection test, includes the parametric intersection coordinate.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Segment<L, T, Q> const &segment, Plane<L, T, Q> const &plane, T &d) {
  return intersects(plane, segment, d);
}

/// Segment/Disc intersection test, includes the parametric intersection coordinate.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Segment<L, T, Q> const &segment, Disc<L, T, Q> const &disc, T &d) {
  return intersects(disc, segment, d);
}

/// Segment/Triangle intersection test, includes the parametric and barycentric intersection coordinates.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Segment<L, T, Q> const &segment, Triangle<L, T, Q> const &triangle, T &d, T &u, T &v) {
  return intersects(triangle, segment, u, v, d);
}

// Segment/Ray intersection test, includes parametric intersection coordinates for segment (d1) and ray (d2)
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Segment<L, T, Q> const &segment, Ray<L, T, Q> const &ray, T &d1, T &d2, T eps = intersect_eps<T>()) {
  return distance(segment, ray, d1, d2) <= eps;
}

// Segment/Line intersection test, includes parametric intersection coordinates for segment (d1) and line (d2)
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Segment<L, T, Q> const &segment, Line<L, T, Q> const &line, T &d1, T &d2, T eps = intersect_eps<T>()) {
  return distance(segment, line, d1, d2) <= eps;
}

// Segment/Segment intersection test, includes parametric intersection coordinates for segment (d1) and other (d2)
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Segment<L, T, Q> const &segment, Segment<L, T, Q> const &other, T &d1, T &d2, T eps = intersect_eps<T>()) {
  return distance2(segment, other, d1, d2) <= eps;
}

/// Segment/AABB intersection test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Segment<L, T, Q> const &segment, AABB<L, T, Q> const &aabb) {
  return intersects(segment, aabb);
}

/// Segment/OBB intersection test
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Segment<3, T, Q> const &segment, OBB<T, Q> const &obb) {
  return interscts(obb, segment);
}

/// Segment/Sphere intersection test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Segment<L, T, Q> const &segment, Sphere<L, T, Q> const &sphere) {
  return intersects(segment, sphere);
}

/// Segment/Plane intersection test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Segment<L, T, Q> const &segment, Plane<L, T, Q> const &plane) {
  const T aDist = signedDistance(plane, segment.a);
  const T bDist = signedDistance(plane, segment.b);
  return aDist * bDist <= T(0);
}

/// Segment/Disc intersection test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Segment<L, T, Q> const &segment, Disc<L, T, Q> const &disc) {
  T d(0);
  return intersects(disc, segment, d);
}

/// Segment/Segment intersection test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Segment<L, T, Q> const &segment, Segment<L, T, Q> const &other, T eps = intersect_eps<T>()) {
  T d1(0), d2(0);
  return intersects(segment, other, d1, d2, eps);
}

/// Segment/Triangle intersection test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Segment<L, T, Q> const &segment, Triangle<L, T, Q> const &triangle) {
  T d(0), u(0), v(0);
  return intersects(triangle, segment, u, v, d);
}

/// @}
GLM_GEOM_END_NAMESPACE
#if GLM_GEOM_TOSTRING
namespace glm {
  namespace detail {
    template<glm::length_t L, typename T, qualifier Q>
    struct compute_to_string<Segment<L, T, Q>> {
      GLM_GEOM_QUALIFIER std::string call(Segment<L, T, Q> const &segment) {
        return detail::format("segment(%s, %s)",
          glm::to_string(segment.a).c_str(),
          glm::to_string(segment.b).c_str()
        );
      }
    };
  }
}
#endif
