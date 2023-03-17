/// @ref geom_line
/// @file geom/line.hpp
///
/// @defgroup geom_line Line
/// @ingroup geom
///

#pragma once

#include "setup.hpp"

#if GLM_MESSAGES == GLM_ENABLE && !defined(GLM_EXT_INCLUDED)
  #pragma message("GLM: GLM_EXT_GEOM_line extension included")
#endif

GLM_GEOM_BEGIN_NAMESPACE
/// @addtogroup geom_line
/// @{

/// <summary>
/// A line represented by P + T * N, where P is the line origin, N is its
/// direction vector, and T is any real number.
/// </summary>
template<length_t L, typename T, qualifier Q>
struct Line {

  // -- Implementation detail --

  typedef T value_type;
  typedef Line<L, T, Q> type;
  typedef vec<L, T, Q> point_type;

  // -- Data --

  point_type pos;  // Origin of this line.
  point_type dir;  // Normalized direction of this line.

#if GLM_CONFIG_DEFAULTED_DEFAULT_CTOR == GLM_ENABLE
  Line() GLM_DEFAULT_CTOR;
#else
  Line()
  #if GLM_CONFIG_CTOR_INIT != GLM_CTOR_INIT_DISABLE
    : pos(T(0)), dir(T(0))
  #endif
  {
  }
#endif

  Line(value_type scalar)
    : pos(scalar), dir(scalar) {
  }

  Line(point_type const &position, point_type const &direction)
    : pos(position), dir(normalize(direction)) {
  }

  Line(Line<L, T, Q> const &line)
    : pos(line.pos), dir(line.dir) {
  }

  Line<L, T, Q> &operator=(Line<L, T, Q> const &line) {
    pos = line.pos;
    dir = line.dir;
    return *this;
  }
};

template<length_t L, typename T, qualifier Q>
static Line<L, T, Q> operator-(Line<L, T, Q> const &line) {
  return Line<L, T, Q>(line.pos, -line.dir);
}

template<length_t L, typename T, qualifier Q>
static bool operator==(Line<L, T, Q> const &l1, Line<L, T, Q> const &l2) {
  return l1.pos == l2.pos && l1.dir == l2.dir;
}

template<length_t L, typename T, qualifier Q>
static Line<L, T, Q> operator+(Line<L, T, Q> const &ray, vec<L, T, Q> const &offset) {
  return Line<L, T, Q>(ray.pos + offset, ray.dir);
}

template<length_t L, typename T, qualifier Q>
static Line<L, T, Q> operator-(Line<L, T, Q> const &ray, vec<L, T, Q> const &offset) {
  return Line<L, T, Q>(ray.pos - offset, ray.dir);
}

template<typename T, qualifier Q>
static Line<3, T, Q> operator*(mat<3, 3, T, Q> const &m, Line<3, T, Q> const &line) {
  return Line<3, T, Q>(m * line.pos, m * line.dir);
}

template<typename T, qualifier Q>
static Line<3, T, Q> operator*(mat<3, 4, T, Q> const &m, Line<3, T, Q> const &line) {
  return Line<3, T, Q>(m * line.pos, m * line.dir);
}

template<typename T, qualifier Q>
static Line<3, T, Q> operator*(mat<4, 3, T, Q> const &m, Line<3, T, Q> const &line) {
  return Line<3, T, Q>(transformPos(m, line.pos), transformDir(m, line.dir));
}

template<typename T, qualifier Q>
static Line<3, T, Q> operator*(mat<4, 4, T, Q> const &m, Line<3, T, Q> const &line) {
  return Line<3, T, Q>(transformPos(m, line.pos), transformDir(m, line.dir));
}

template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER Line<3, T, Q> operator*(qua<T, Q> const &q, Line<3, T, Q> const &line) {
  return Line<3, T, Q>(q * line.pos, q * line.dir);
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool equal(Line<L, T, Q> const &x, Line<L, T, Q> const &y, T eps = epsilon<T>()) {
  return all(equal(x.pos, y.pos, eps)) && all(equal(x.dir, y.dir, eps));
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool equal(Line<L, T, Q> const &x, Line<L, T, Q> const &y, vec<L, T, Q> const &eps) {
  return all(equal(x.pos, y.pos, eps)) && all(equal(x.dir, y.dir, eps));
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool equal(Line<L, T, Q> const &x, Line<L, T, Q> const &y, int MaxULPs) {
  return all(equal(x.pos, y.pos, MaxULPs)) && all(equal(x.dir, y.dir, MaxULPs));
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool equal(Line<L, T, Q> const &x, Line<L, T, Q> const &y, vec<L, int, Q> const &MaxULPs) {
  return all(equal(x.pos, y.pos, MaxULPs)) && all(equal(x.dir, y.dir, MaxULPs));
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool notEqual(Line<L, T, Q> const &x, Line<L, T, Q> const &y, T eps = epsilon<T>()) {
  return any(notEqual(x.pos, y.pos, eps)) || any(notEqual(x.dir, y.dir, eps));
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool notEqual(Line<L, T, Q> const &x, Line<L, T, Q> const &y, vec<L, T, Q> const &eps) {
  return any(notEqual(x.pos, y.pos, eps)) || any(notEqual(x.dir, y.dir, eps));
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool notEqual(Line<L, T, Q> const &x, Line<L, T, Q> const &y, int MaxULPs) {
  return any(notEqual(x.pos, y.pos, MaxULPs)) || any(notEqual(x.dir, y.dir, MaxULPs));
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool notEqual(Line<L, T, Q> const &x, Line<L, T, Q> const &y, vec<L, int, Q> const &MaxULPs) {
  return any(notEqual(x.pos, y.pos, MaxULPs)) || any(notEqual(x.dir, y.dir, MaxULPs));
}

/* Forward declarations */

/// Return the closest point pair on two Lines.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool closestPointLineLine(vec<L, T, Q> const &v0, vec<L, T, Q> const &v10, vec<L, T, Q> const &v2, vec<L, T, Q> const &v32, T &d1, T &d2);

/// Return the closest point (along an edge) between the Triangle and provided Line.
template<length_t L, typename T, qualifier Q, typename Line>
GLM_GEOM_QUALIFIER_OUTLINE vec<L, T, Q> closestPointTriangleLine(Triangle<L, T, Q> const &t, Line const &line, T &outU, T &outV, T &outD);

/// Test if any component of the line is infinite.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool isinf(Line<L, T, Q> const &line) {
  return any(isinf(line.pos)) || any(isinf(line.dir));
}

/// Test if any component of the line is NaN
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool isnan(Line<L, T, Q> const &line) {
  return any(isnan(line.pos)) || any(isnan(line.dir));
}

/// Test if all components of the line are finite.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool isfinite(Line<L, T, Q> const &line) {
  return all(isfinite(line.pos)) && all(isfinite(line.dir));
}

/// Get a point along the Line at a given distance (parametric point).
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> getPoint(Line<L, T, Q> const &line, T d) {
  return line.pos + d * line.dir;
}

/* Return the closest point on this line to the given object */

/// Return the closest point between the Line and given point; including parametric line coordinate.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Line<L, T, Q> const &line, vec<L, T, Q> const &point, T &d) {
  d = dot(point - line.pos, line.dir);
  return getPoint(line, d);
}

/// Return the closest point between two Line's; including parametric line coordinate.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Line<L, T, Q> const &line, Line<L, T, Q> const &other, T &d1, T &d2) {
  closestPointLineLine(line.pos, line.dir, other.pos, other.dir, d1, d2);
  return getPoint(line, d1);
}

/// Return the closest point between a Line and Segment; including parametric coordinates.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Line<L, T, Q> const &line, Segment<L, T, Q> const &segment, T &d1, T &d2) {
  closestPointLineLine(line.pos, line.dir, segment.a, segment.dir2(), d1, d2);
  if (d2 < T(0)) {
    d2 = T(0);
    return closestPoint(line, segment.a, d1);
  }
  else if (d2 > T(1)) {
    d2 = T(1);
    return closestPoint(line, segment.b, d1);
  }
  return getPoint(line, d1);
}

/// Return the closest point between a Line and Ray; including parametric coordinates.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Line<L, T, Q> const &line, Ray<L, T, Q> const &ray, T &d1, T &d2) {
  closestPointLineLine(line.pos, line.dir, ray.pos, ray.dir, d1, d2);
  if (d2 >= T(0))
    return getPoint(line, d1);
  else {
    d2 = T(0);
    return closestPoint(line, ray.pos, d1);
  }
}

/// Return the closest point between a Line and Triangle; including parametric Line and Triangle coordinates.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Line<L, T, Q> const &line, Triangle<L, T, Q> const &triangle, T &d, T &u, T &v) {
  d = intersectTriangleLine(triangle, line.pos, line.dir, u, v);  // Compute distance along line
  if (glm::isinf(d)) {
    closestPointTriangleLine<L, T, Q, Line<L, T, Q>>(triangle, line, u, v, d);
  }
  return getPoint(line, d);
}

/// Return the closest point between a Line and point
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Line<L, T, Q> const &line, vec<L, T, Q> const &point) {
  T d(0);
  return closestPoint(line, point, d);
}

/// Return the closest point between two Lines
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Line<L, T, Q> const &line, Line<L, T, Q> const &other) {
  T d1(0), d2(0);
  return closestPoint(line, other, d1, d2);
}

/// Return the closest point between a Line and Segment
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Line<L, T, Q> const &line, Segment<L, T, Q> const &segment) {
  T d1(0), d2(0);
  return closestPoint(line, segment, d1, d2);
}

/// Return the closest point between a Line and Ray
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Line<L, T, Q> const &line, Ray<L, T, Q> const &ray) {
  T d1(0), d2(0);
  return closestPoint(line, ray, d1, d2);
}

/// Return the closest point between a Line and Triangle
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Line<L, T, Q> const &line, Triangle<L, T, Q> const &triangle) {
  T u, v, d(0);
  return closestPoint(line, triangle, d, u, v);
}

/* Test if the given object is fully contained on the line */

/// Return true if the point is on the Line.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool contains(Line<L, T, Q> const &line, vec<L, T, Q> const &point, T distanceThreshold = contains_eps<T>()) {
  T d(0);
  return distance2(closestPoint(line, point, d), point) <= distanceThreshold;
}

/// Return true if the Line and Ray are equal
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool contains(Line<L, T, Q> const &line, Ray<L, T, Q> const &ray, T distanceThreshold = contains_eps<T>()) {
  return contains(line, ray.pos, distanceThreshold) && all(epsilonEqual(line.dir, ray.dir, distanceThreshold));
}

/// Return true if the Segment lays on the Line
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool contains(Line<L, T, Q> const &line, Segment<L, T, Q> const &segment, T distanceThreshold = contains_eps<T>()) {
  return contains(line, segment.a, distanceThreshold) && contains(line, segment.b, distanceThreshold);
}

/* Return the distance between the line and the given object. */

/// Return the distance between the Line and given point; including parametric line coordinate.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Line<L, T, Q> const &line, vec<L, T, Q> const &point, T &d) {
  return distance(closestPoint(line, point, d), point);
}

/// Return the closest point between the Line and given Ray; including parametric coordinates for line (d1) and ray (d2)
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Line<L, T, Q> const &line, Ray<L, T, Q> const &ray, T &d1, T &d2) {
  const vec<L, T, Q> point = closestPoint(line, ray, d1, d2);
  return distance(point, getPoint(ray, d2));
}

/// Return the distance between two Lines; including parametric coordinates both lines
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Line<L, T, Q> const &line, Line<L, T, Q> const &other, T &d1, T &d2) {
  const vec<L, T, Q> point = closestPoint(line, other, d1, d2);
  return distance(point, getPoint(other, d2));
}

/// Return the closest point between the Line and given Segment; including parametric coordinates for line (d1) and line (d2)
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Line<L, T, Q> const &line, Segment<L, T, Q> const &segment, T &d1, T &d2) {
  const vec<L, T, Q> point = closestPoint(line, segment, d1, d2);
  return (d2 >= T(0) && d2 <= T(1)) ? distance(point, getPoint(segment, d2)) : T(-1);  // Invalid
}

/// Return the distance between the Line and given Sphere
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Line<L, T, Q> const &line, Sphere<L, T, Q> const &sphere) {
  T ignore = T(0);
  return max(T(0), distance(line, sphere.pos, ignore) - sphere.r);
}

/// Return the distance between the Line and given Capsule
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Line<L, T, Q> const &line, Capsule<L, T, Q> const &capsule) {
  return max(T(0), distance(line, capsule.l) - capsule.r);
}

/// Return the distance between the Line and given point
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Line<L, T, Q> const &line, vec<L, T, Q> const &point) {
  T d(0);
  return distance(line, point, d);
}

/// Return the distance between the Line and given Ray
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Line<L, T, Q> const &line, Ray<L, T, Q> const &ray) {
  T d1(0), d2(0);
  return distance(line, ray, d1, d2);
}

/// Return the distance between two Lines
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Line<L, T, Q> const &line, Line<L, T, Q> const &other) {
  T d1(0), d2(0);
  return distance(line, other, d1, d2);
}

/// Return the distance between the Line and given Segment
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Line<L, T, Q> const &line, Segment<L, T, Q> const &segment) {
  T d1(0), d2(0);
  return distance(line, segment, d1, d2);
}

/* Test whether the line and the given object intersect */

// Line/AABB intersection test, includes parametric intersection coordinates
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Line<L, T, Q> const &line, AABB<L, T, Q> const &aabb, T &dNear, T &dFar) {
  return intersects(aabb, line, dNear, dFar);
}

// Line/Ray intersection test, includes parametric intersection coordinates for line (d1) and ray (d2)
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Line<L, T, Q> const &line, Ray<L, T, Q> const &ray, T &d1, T &d2, T eps = intersect_eps<T>()) {
  return distance(line, ray, d1, d2) <= eps;
}

// Line/Line intersection test, includes parametric intersection coordinates for line (d1) and other (d2)
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Line<L, T, Q> const &line, Line<L, T, Q> const &other, T &d1, T &d2, T eps = intersect_eps<T>()) {
  return distance(line, other, d1, d2) <= eps;
}

// Line/Segment intersection test, includes parametric intersection coordinates for line (d1) and segment (d2)
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Line<L, T, Q> const &line, Segment<L, T, Q> const &segment, T &d1, T &d2, T eps = intersect_eps<T>()) {
  return distance(line, segment, d1, d2) <= eps;
}

// Line/Sphere intersection test, includes parametric intersection coordinates
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER int intersects(Line<L, T, Q> const &line, Sphere<L, T, Q> const &s, T &d1, T &d2) {
  return intersects(s, line, d1, d2);
}

// Line/Capsule intersection test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Line<L, T, Q> const &line, Capsule<L, T, Q> const &capsule) {
  return intersects(capsule, line);
}

// Line/Plane intersection test, includes parametric intersection coordinate
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Line<L, T, Q> const &line, Plane<L, T, Q> const &plane, T &d) {
  return intersects(plane, line, d);
}

/// Line/Triangle intersection test, includes the parametric and barycentric intersection coordinates.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Line<L, T, Q> const &line, Triangle<L, T, Q> const &triangle, T &d, T &u, T &v) {
  return intersects(triangle, line, u, v, d);
}

/// Line/Disc intersection test, includes the parametric intersection coordinate.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Line<L, T, Q> const &line, Disc<L, T, Q> const &disc, T &d) {
  return intersects(disc, line, d);
}

// Ray/OBB intersection test, includes parametric intersection coordinates
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Line<3, T, Q> const &line, OBB<T, Q> const &obb, T &dNear, T &dFar) {
  return intersects(obb, line, dNear, dFar);
}

/// Line/AABB intersection test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Line<L, T, Q> const &line, AABB<L, T, Q> const &aabb) {
  return intersects(aabb, line);
}

/// Line/Sphere intersection test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Line<L, T, Q> const &line, Sphere<L, T, Q> const &s) {
  T d1(0), d2(0);
  return intersects(line, s, d1, d2);
}

/// Line/Plane intersection test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Line<L, T, Q> const &line, Plane<L, T, Q> const &plane) {
  T d(0);
  return intersects(line, plane, d);
}

/// Line/Triangle intersection test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Line<L, T, Q> const &line, Triangle<L, T, Q> const &triangle) {
  T u, v, d(0);
  return intersects(line, triangle, d, u, v);
}

/// Line/Disc intersection test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Line<L, T, Q> const &line, Disc<L, T, Q> const &disc) {
  T d(0);
  return intersects(disc, line, d);
}

/// Line/OBB intersection test
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Line<3, T, Q> const &line, OBB<T, Q> const &obb) {
  return intersects(obb, line);
}

/// Convert the Line to a Segment.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER Segment<L, T, Q> toSegment(Line<L, T, Q> const &line, T end) {
  return Segment<L, T, Q>(line.pos, getPoint(line, end));
}

/// Convert the Line to a Segment.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER Segment<L, T, Q> toSegment(Line<L, T, Q> const &line, T start, T end) {
  return Segment<L, T, Q>(getPoint(line, start), getPoint(line, end));
}

/// Project the Line onto the given axis (direction), i.e., collapse the line onto an axis.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER void projectToAxis(Line<L, T, Q> const &line, vec<L, T, Q> const &direction, T &outMin, T &outMax) {
  if (areOrthonormal2(line.dir, direction, epsilon<T>()))
    outMin = outMax = dot(direction, line.pos);
  else {
    outMin = -std::numeric_limits<T>::infinity();
    outMax = std::numeric_limits<T>::infinity();
  }
}

/// @private
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool closestPointLineLine(vec<L, T, Q> const &v0, vec<L, T, Q> const &v1, vec<L, T, Q> const &v2, vec<L, T, Q> const &v3, T &d1, T &d2) {
  const bool isNullV1 = isNull(v1, epsilon<T>());
  const bool isNullV3 = isNull(v3, epsilon<T>());
  if (isNullV1 || isNullV3) {  // at least one line is degenerate
    d1 = isNullV1 ? T(0) : T(0.5);
    d2 = isNullV3 ? T(0) : T(0.5);
    return false;
  }

  const T d33 = dot(v3, v3);
  if (!detail::exactly_zero(d33)) {
    const vec<L, T, Q> v4 = v0 - v2;
    const T d43 = dot(v4, v3);
    const T d31 = dot(v3, v1);
    const T denom = dot(v1, v1) * d33 - d31 * d31;
    d1 = detail::exactly_zero(denom) ? T(0) : (d43 * d31 - dot(v4, v1) * d33) / denom;
    d2 = (d43 + d1 * d31) / d33;
    return true;
  }
  else {
    d1 = d2 = T(0);
    return false;
  }
}

/// @}
GLM_GEOM_END_NAMESPACE
#if GLM_GEOM_TOSTRING
namespace glm {
  namespace detail {
    template<glm::length_t L, typename T, qualifier Q>
    struct compute_to_string<Line<L, T, Q>> {
      GLM_GEOM_QUALIFIER std::string call(Line<L, T, Q> const &line) {
        return detail::format("line(%s, %s)",
          glm::to_string(line.pos).c_str(),
          glm::to_string(line.dir).c_str()
        );
      }
    };
  }
}
#endif
