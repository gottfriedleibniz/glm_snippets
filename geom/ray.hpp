/// @ref geom_ray
/// @file geom/ray.hpp
///
/// @defgroup geom_ray Ray
/// @ingroup geom
///

#pragma once

#include "setup.hpp"
#include "line.hpp"
#include "segment.hpp"
#include "triangle.hpp"

#include <glm/gtx/compatibility.hpp>
#if GLM_MESSAGES == GLM_ENABLE && !defined(GLM_EXT_INCLUDED)
  #pragma message("GLM: GLM_EXT_GEOM_ray extension included")
#endif

GLM_GEOM_BEGIN_NAMESPACE
/// @addtogroup geom_ray
/// @{

/// <summary>
/// A line represented by P + T * N, where P is the ray origin, N is its
/// direction vector, and T >= 0 (extending to infinity).
/// </summary>
template<length_t L, typename T, qualifier Q>
struct Ray : Line<L, T, Q> {

  // -- Implementation detail --

  typedef T value_type;
  typedef Ray<L, T, Q> type;
  typedef vec<L, T, Q> point_type;

  // -- Data --

  using Line<L, T, Q>::pos;  // Origin of this ray
  using Line<L, T, Q>::dir;  // Normalized direction vector of this ray

#if GLM_CONFIG_DEFAULTED_DEFAULT_CTOR == GLM_ENABLE
  Ray() GLM_DEFAULT_CTOR;
#else
  Ray()
  #if GLM_CONFIG_CTOR_INIT != GLM_CTOR_INIT_DISABLE
    : Line<L, T, Q>()
  #endif
  {
  }
#endif

  Ray(T scalar)
    : Line<L, T, Q>(scalar) {
  }

  Ray(point_type const &position, point_type const &direction)
    : Line<L, T, Q>(position, direction) {
  }

  Ray(Line<L, T, Q> const &line)
    : Line<L, T, Q>(line) {
  }

  Ray(Ray<L, T, Q> const &ray)
    : Line<L, T, Q>(ray) {
  }

  Ray<L, T, Q> &operator=(Ray<L, T, Q> const &ray) {
    pos = ray.pos;
    dir = ray.dir;
    return *this;
  }
};

template<length_t L, typename T, qualifier Q>
static Ray<L, T, Q> operator-(Ray<L, T, Q> const &ray) {
  return Ray<L, T, Q>(ray.pos, -ray.dir);
}

template<length_t L, typename T, qualifier Q>
static bool operator==(Ray<L, T, Q> const &r1, Ray<L, T, Q> const &r2) {
  return r1.pos == r2.pos && r1.dir == r2.dir;
}

template<length_t L, typename T, qualifier Q>
static Ray<L, T, Q> operator+(Ray<L, T, Q> const &ray, vec<L, T, Q> const &offset) {
  return Ray<L, T, Q>(ray.pos + offset, ray.dir);
}

template<length_t L, typename T, qualifier Q>
static Ray<L, T, Q> operator-(Ray<L, T, Q> const &ray, vec<L, T, Q> const &offset) {
  return Ray<L, T, Q>(ray.pos - offset, ray.dir);
}

template<typename T, qualifier Q>
static Ray<3, T, Q> operator*(mat<3, 3, T, Q> const &m, Ray<3, T, Q> const &ray) {
  return Ray<3, T, Q>(m * ray.pos, m * ray.dir);
}

template<typename T, qualifier Q>
static Ray<3, T, Q> operator*(mat<3, 4, T, Q> const &m, Ray<3, T, Q> const &ray) {
  return Ray<3, T, Q>(m * ray.pos, m * ray.dir);
}

template<typename T, qualifier Q>
static Ray<3, T, Q> operator*(mat<4, 3, T, Q> const &m, Ray<3, T, Q> const &ray) {
  return Ray<3, T, Q>(transformPos(m, ray.pos), transformDir(m, ray.dir));
}

template<typename T, qualifier Q>
static Ray<3, T, Q> operator*(mat<4, 4, T, Q> const &m, Ray<3, T, Q> const &ray) {
  return Ray<3, T, Q>(transformPos(m, ray.pos), transformDir(m, ray.dir));
}

template<typename T, qualifier Q>
static Ray<3, T, Q> operator*(qua<T, Q> const &q, Ray<3, T, Q> const &ray) {
  return Ray<3, T, Q>(q * ray.pos, q * ray.dir);
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool equal(Ray<L, T, Q> const &x, Ray<L, T, Q> const &y, T eps = epsilon<T>()) {
  return all(equal(x.pos, y.pos, eps)) && all(equal(x.dir, y.dir, eps));
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool equal(Ray<L, T, Q> const &x, Ray<L, T, Q> const &y, vec<L, T, Q> const &eps) {
  return all(equal(x.pos, y.pos, eps)) && all(equal(x.dir, y.dir, eps));
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool equal(Ray<L, T, Q> const &x, Ray<L, T, Q> const &y, int MaxULPs) {
  return all(equal(x.pos, y.pos, MaxULPs)) && all(equal(x.dir, y.dir, MaxULPs));
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool equal(Ray<L, T, Q> const &x, Ray<L, T, Q> const &y, vec<L, int, Q> const &MaxULPs) {
  return all(equal(x.pos, y.pos, MaxULPs)) && all(equal(x.dir, y.dir, MaxULPs));
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool notEqual(Ray<L, T, Q> const &x, Ray<L, T, Q> const &y, T eps = epsilon<T>()) {
  return any(notEqual(x.pos, y.pos, eps)) || any(notEqual(x.dir, y.dir, eps));
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool notEqual(Ray<L, T, Q> const &x, Ray<L, T, Q> const &y, vec<L, T, Q> const &eps) {
  return any(notEqual(x.pos, y.pos, eps)) || any(notEqual(x.dir, y.dir, eps));
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool notEqual(Ray<L, T, Q> const &x, Ray<L, T, Q> const &y, int MaxULPs) {
  return any(notEqual(x.pos, y.pos, MaxULPs)) || any(notEqual(x.dir, y.dir, MaxULPs));
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool notEqual(Ray<L, T, Q> const &x, Ray<L, T, Q> const &y, vec<L, int, Q> const &MaxULPs) {
  return any(notEqual(x.pos, y.pos, MaxULPs)) || any(notEqual(x.dir, y.dir, MaxULPs));
}

/// Test if any component of the ray is infinite.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool isinf(Ray<L, T, Q> const &line) {
  return any(isinf(line.pos)) || any(isinf(line.dir));
}

/// Test if any component of the ray is NaN
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool isnan(Ray<L, T, Q> const &line) {
  return any(isnan(line.pos)) || any(isnan(line.dir));
}

/// Test if all components of the ray are finite
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool isfinite(Ray<L, T, Q> const &line) {
  return all(isfinite(line.pos)) && all(isfinite(line.dir));
}

/// <summary>
/// Get a point along the Ray at a given distance (parametric value). Passing
/// negative values to this function treats the ray as if it were a line.
/// </summary>
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> getPoint(Ray<L, T, Q> const &ray, T d) {
  return ray.pos + d * ray.dir;
}

/* Return the closest point on this ray to the given object */

/// Return the closest point between the Ray and given point; including parametric ray coordinate.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Ray<L, T, Q> const &ray, vec<L, T, Q> const &point, T &d) {
  d = max(T(0), dot(point - ray.pos, ray.dir));
  return getPoint(ray, d);
}

/// Return the closest point between the Ray and given line; including parametric coordinates for ray (d1) and line (d2)
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER_OUTLINE vec<L, T, Q> closestPoint(Ray<L, T, Q> const &ray, Line<L, T, Q> const &line, T &d1, T &d2) {
  closestPointLineLine(ray.pos, ray.dir, line.pos, line.dir, d1, d2);
  if (d1 < T(0)) {
    d1 = T(0);
    closestPoint(line, ray.pos, d2);
    return ray.pos;
  }
  return getPoint(ray, d1);
}

/// Return the closest point between two Rays; including parametric coordinates both rays
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER_OUTLINE vec<L, T, Q> closestPoint(Ray<L, T, Q> const &ray, Ray<L, T, Q> const &other, T &d1, T &d2) {
  closestPointLineLine(ray.pos, ray.dir, other.pos, other.dir, d1, d2);
  if (d1 < T(0) && d2 < T(0)) {
    const vec<L, T, Q> pt = closestPoint(ray, other.pos, d1);
    const vec<L, T, Q> pt2 = closestPoint(other, ray.pos, d2);
    if (distance2(pt, other.pos) <= distance2(pt2, ray.pos)) {
      d2 = T(0);
      return pt;
    }
    else {
      d1 = T(0);
      return ray.pos;
    }
  }
  else if (d1 < T(0)) {
    closestPoint(other, ray.pos, d2);
    d1 = T(0);
    d2 = max(T(0), d2);
    return ray.pos;
  }
  else if (d2 < T(0)) {
    const vec<L, T, Q> pt = closestPoint(ray, other.pos, d1);
    d1 = max(T(0), d1);
    d2 = T(0);
    return pt;
  }
  else {
    return getPoint(ray, d1);
  }
}

/// <summary>
/// Return the closest point between the Ray and given Segment including the
/// parametric coordinates for ray (d1) and segment (d2)
/// </summary>
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER_OUTLINE vec<L, T, Q> closestPoint(Ray<L, T, Q> const &ray, Segment<L, T, Q> const &segment, T &d1, T &d2) {
  closestPointLineLine(ray.pos, ray.dir, segment.a, segment.dir2(), d1, d2);
  if (d1 < T(0)) {
    d1 = T(0);
    if (d2 >= T(0) && d2 <= T(1)) {
      closestPoint(segment, ray.pos, d2);
      return ray.pos;
    }
    else {
      const T t2 = (d2 < T(0)) ? T(0) : T(1);
      const vec<L, T, Q> p = (d2 < T(0)) ? segment.a : segment.b;
      const vec<L, T, Q> pt = closestPoint(ray, p, d1);
      const vec<L, T, Q> pt2 = closestPoint(segment, ray.pos, d2);
      if (distance2(pt, p) <= distance2(pt2, ray.pos)) {
        d2 = t2;
        return pt;
      }
      else {
        d1 = T(0);
        return ray.pos;
      }
    }
  }
  else if (d2 < T(0)) {
    d2 = T(0);
    return closestPoint(ray, segment.a, d1);
  }
  else if (d2 > T(1)) {
    d2 = T(1);
    return closestPoint(ray, segment.b, d1);
  }
  return getPoint(ray, d1);
}

/// Return the closest point between a Ray and point
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Ray<L, T, Q> const &ray, vec<L, T, Q> const &point) {
  T d(0);
  return closestPoint(ray, point, d);
}

/// Return the closest point between a Ray and Line
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Ray<L, T, Q> const &ray, Line<L, T, Q> const &line) {
  T d1(0), d2(0);
  return closestPoint(ray, line, d1, d2);
}

/// Return the closest point between a Ray and Segment
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Ray<L, T, Q> const &ray, Segment<L, T, Q> const &segment) {
  T d1(0), d2(0);
  return closestPoint(ray, segment, d1, d2);
}

/// Return the closest point between two Rays
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Ray<L, T, Q> const &ray, Ray<L, T, Q> const &other) {
  T d1(0), d2(0);
  return closestPoint(ray, other, d1, d2);
}

/* Test if the given object is fully contained on the ray */

/// Return true if the point is on the Ray.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool contains(Ray<L, T, Q> const &ray, vec<L, T, Q> const &point, T distanceThreshold = contains_eps<T>()) {
  T ignore(0);
  return distance2(closestPoint(ray, point, ignore), point) <= distanceThreshold;
}

/// Return true if the Segment lays on the Ray
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool contains(Ray<L, T, Q> const &ray, Segment<L, T, Q> const &segment, T distanceThreshold = contains_eps<T>()) {
  return contains(ray, segment.a, distanceThreshold) && contains(ray, segment.b, distanceThreshold);
}

/* Return the distance between the ray and the given object */

/// Return the distance between the Ray and given point; including parametric ray coordinate.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Ray<L, T, Q> const &ray, vec<L, T, Q> const &point, T &d) {
  return distance(closestPoint(ray, point, d), point);
}

/// Return the distance between two Rays; including parametric coordinates both rays
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Ray<L, T, Q> const &ray, Ray<L, T, Q> const &other, T &d1, T &d2) {
  const vec<L, T, Q> point = closestPoint(ray, other, d1, d2);
  return distance(point, getPoint(other, d2));
}

/// Return the closest point between the Ray and given Line; including parametric coordinates for ray (d1) and line (d2)
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Ray<L, T, Q> const &ray, Line<L, T, Q> const &line, T &d1, T &d2) {
  const vec<L, T, Q> point = closestPoint(ray, line, d1, d2);
  return distance(point, getPoint(line, d2));
}

/// Return the closest point between the Ray and given Segment; including parametric coordinates for ray (d1) and line (d2)
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Ray<L, T, Q> const &ray, Segment<L, T, Q> const &segment, T &d1, T &d2) {
  const vec<L, T, Q> point = closestPoint(ray, segment, d1, d2);
  return distance(point, getPoint(segment, d2));
}

/// Return the distance between the Ray and given Sphere
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Ray<L, T, Q> const &ray, Sphere<L, T, Q> const &sphere) {
  return max(T(0), distance(ray, sphere.pos) - sphere.r);
}

/// Return the distance between the Ray and given Capsule
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Ray<L, T, Q> const &ray, Capsule<L, T, Q> const &capsule) {
  return max(T(0), distance(ray, capsule.l) - capsule.r);
}

/// Return the distance between the Ray and given point
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Ray<L, T, Q> const &ray, vec<L, T, Q> const &point) {
  T d(0);
  return distance(ray, point, d);
}

/// Return the distance between two Rays
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Ray<L, T, Q> const &ray, Ray<L, T, Q> const &other) {
  T d1(0), d2(0);
  return distance(ray, other, d1, d2);
}

/// Return the distance between the AABB and given Line
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Ray<L, T, Q> const &ray, Line<L, T, Q> const &line) {
  T d1(0), d2(0);
  return distance(ray, line, d1, d2);
}

/// Return the distance between the AABB and given Segment
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Ray<L, T, Q> const &ray, Segment<L, T, Q> const &segment) {
  T d1(0), d2(0);
  return distance(ray, segment, d1, d2);
}

/* Test whether the ray and the given object intersect */

// Ray/Sphere intersection test, includes parametric intersection coordinates
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER int intersects(Ray<L, T, Q> const &ray, Sphere<L, T, Q> const &sphere, T &d1, T &d2) {
  return intersects(sphere, ray, d1, d2);
}

// Ray/AABB intersection test, includes parametric intersection coordinates
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Ray<L, T, Q> const &ray, AABB<L, T, Q> const &aabb, T &d1, T &d2) {
  return intersects(aabb, ray, d1, d2);
}

// Ray/OBB intersection test, includes parametric intersection coordinates
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Ray<3, T, Q> const &ray, OBB<T, Q> const &obb, T &dNear, T &dFar) {
  return intersects(obb, ray, dNear, dFar);
}

// Ray/Ray intersection test, includes parametric intersection coordinates for ray (d1) and other (d2)
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Ray<L, T, Q> const &ray, Ray<L, T, Q> const &other, T &d1, T &d2, T eps = intersect_eps<T>()) {
  return distance(ray, other, d1, d2) <= eps;
}

// Ray/Line intersection test, includes parametric intersection coordinates for ray (d1) and line (d2)
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Ray<L, T, Q> const &ray, Line<L, T, Q> const &line, T &d1, T &d2, T eps = intersect_eps<T>()) {
  return distance(ray, line, d1, d2) <= eps;
}

// Ray/Segment intersection test, includes parametric intersection coordinates for ray (d1) and segment (d2)
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Ray<L, T, Q> const &ray, Segment<L, T, Q> const &segment, T &d1, T &d2, T eps = intersect_eps<T>()) {
  return distance(ray, segment, d1, d2) <= eps;
}

/// Ray/Plane intersection test, includes the parametric intersection coordinate.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Ray<L, T, Q> const &ray, Plane<L, T, Q> const &plane, T &d) {
  return intersects(plane, ray, d);
}

/// Ray/Disc intersection test, includes the parametric intersection coordinate.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Ray<L, T, Q> const &ray, Disc<L, T, Q> const &disc, T &d) {
  return intersects(disc, ray, d);
}

/// Ray/Triangle intersection test, includes the parametric and barycentric intersection coordinates.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Ray<L, T, Q> const &ray, Triangle<L, T, Q> const &triangle, T &d, T &u, T &v) {
  d = intersectTriangleLine(triangle, ray.pos, ray.dir, u, v);
  return glm::isfinite(d) && d >= T(0);
}

/// Ray/Capsule intersection test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Ray<L, T, Q> const &ray, Capsule<L, T, Q> const &capsule) {
  return intersects(capsule, ray);
}

/// Ray/Sphere intersection test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Ray<L, T, Q> const &ray, Sphere<L, T, Q> const &sphere) {
  T d1(0), d2(0);
  return intersects(ray, sphere, d1, d2) > 0;
}

/// Ray/AABB intersection test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Ray<L, T, Q> const &ray, AABB<L, T, Q> const &aabb) {
  T d1(0), d2(0);
  return intersects(ray, aabb, d1, d2) > 0;
}

/// Ray/OBB intersection test
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Ray<3, T, Q> const &ray, OBB<T, Q> const &obb) {
  return intersects(obb, ray);
}

/// Ray/Plane intersection test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Ray<L, T, Q> const &ray, Plane<L, T, Q> const &plane) {
  T d(0);
  return intersects(ray, plane, d);
}

/// Ray/Disc intersection test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Ray<L, T, Q> const &ray, Disc<L, T, Q> const &disc) {
  T d(0);
  return intersects(disc, ray, d);
}

/// Ray/Triangle intersection test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Ray<L, T, Q> const &ray, Triangle<L, T, Q> const &triangle) {
  T u, v, d(0);
  return intersects(ray, triangle, d, u, v);
}

/// Convert a Ray to a Segment.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER Segment<L, T, Q> toSegment(Ray<L, T, Q> const &ray, T dEnd) {
  return Segment<L, T, Q>(ray.pos, getPoint(ray, dEnd));
}

/// Convert a Ray to a Segment.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER Segment<L, T, Q> toSegment(Ray<L, T, Q> const &ray, T dStart, T dEnd) {
  return Segment<L, T, Q>(getPoint(ray, dStart), getPoint(ray, dEnd));
}

/// Project the Ray onto the given axis (direction), i.e., collapse the ray onto an axis.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER void projectToAxis(Ray<L, T, Q> const &ray, vec<L, T, Q> const &direction, T &outMin, T &outMax) {
  const T d = dot(direction, ray.dir);
  outMin = outMax = dot(direction, ray.pos);
  if (d > epsilon<T>())
    outMax = std::numeric_limits<T>::infinity();
  else if (d < -epsilon<T>())
    outMin = -std::numeric_limits<T>::infinity();
}

/// @}
GLM_GEOM_END_NAMESPACE
#if GLM_GEOM_TOSTRING
namespace glm {
  namespace detail {
    template<glm::length_t L, typename T, qualifier Q>
    struct compute_to_string<Ray<L, T, Q>> {
      GLM_GEOM_QUALIFIER std::string call(Ray<L, T, Q> const &ray) {
        return detail::format("ray(%s, %s)",
          glm::to_string(ray.pos).c_str(),
          glm::to_string(ray.dir).c_str()
        );
      }
    };
  }
}
#endif
