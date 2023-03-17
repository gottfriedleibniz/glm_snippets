/// @ref geom_sphere
/// @file geom/sphere.hpp
///
/// @defgroup geom_sphere Sphere
/// @ingroup geom
///

#pragma once

#include <algorithm>

#include "setup.hpp"
#include "line.hpp"
#include "segment.hpp"
#include "ray.hpp"
#include "triangle.hpp"

#if GLM_MESSAGES == GLM_ENABLE && !defined(GLM_EXT_INCLUDED)
  #pragma message("GLM: GLM_EXT_GEOM_sphere extension included")
#endif

GLM_GEOM_BEGIN_NAMESPACE
/// @addtogroup geom_sphere
/// @{

/// <summary>
/// A hypersphere represented by |X - P| = R where P is the sphere centroid and
/// R is its radius.
/// </summary>
template<length_t L, typename T, qualifier Q>
struct Sphere {

  // -- Implementation detail --

  typedef T value_type;
  typedef Sphere<L, T, Q> type;
  typedef vec<L, T, Q> point_type;

  // -- Data --

  point_type pos;  // Center point of this sphere
  value_type r;  // Radius of this sphere

#if GLM_CONFIG_DEFAULTED_DEFAULT_CTOR == GLM_ENABLE
  Sphere() GLM_DEFAULT_CTOR;
#else
  Sphere()
  #if GLM_CONFIG_CTOR_INIT != GLM_CTOR_INIT_DISABLE
    : pos(T(0)), r(T(0))
  #endif
  {
  }
#endif

  Sphere(T scalar)
    : pos(scalar), r(scalar) {
  }

  Sphere(point_type const &position, value_type radius)
    : pos(position), r(radius) {
  }

  Sphere(Sphere<L, T, Q> const &sphere)
    : pos(sphere.pos), r(sphere.r) {
  }

  Sphere<L, T, Q> &operator=(Sphere<L, T, Q> const &sphere) {
    pos = sphere.pos;
    r = sphere.r;
    return *this;
  }

  void setDegenerate() {
    pos = point_type(std::numeric_limits<T>::quiet_NaN());
    r = std::numeric_limits<T>::quiet_NaN();
  }

  void enclose(point_type const &point, T eps = T(0)) {
    const point_type d = point - pos;
    const T dist2 = length2(d);
    if (dist2 + eps > r * r) {
      const T dist = sqrt(dist2);
      const T halfDist = (dist - r) * T(0.5);
      pos += d * halfDist / dist;
      r += halfDist + epsilon<T>();  // see rationale in optimalEnclosingSphere
    }
  }
};

template<length_t L, typename T, qualifier Q>
static Sphere<L, T, Q> operator-(Sphere<L, T, Q> const &sphere) {
  return Sphere<L, T, Q>(-sphere.pos, sphere.r);
}

template<length_t L, typename T, qualifier Q>
static bool operator==(Sphere<L, T, Q> const &s1, Sphere<L, T, Q> const &s2) {
  return s1.pos == s2.pos && detail::equal_to(s1.r, s2.r);
}

template<length_t L, typename T, qualifier Q>
static Sphere<L, T, Q> operator+(Sphere<L, T, Q> const &sphere, vec<L, T, Q> const &offset) {
  return Sphere<L, T, Q>(sphere.pos + offset, sphere.r);
}

template<length_t L, typename T, qualifier Q>
static Sphere<L, T, Q> operator-(Sphere<L, T, Q> const &sphere, vec<L, T, Q> const &offset) {
  return Sphere<L, T, Q>(sphere.pos - offset, sphere.r);
}

template<typename T, qualifier Q>
static Sphere<3, T, Q> operator*(mat<3, 3, T, Q> const &m, Sphere<3, T, Q> const &sphere) {
  return Sphere<3, T, Q>(m * sphere.pos, length(m[0]) * sphere.r);
}

template<typename T, qualifier Q>
static Sphere<3, T, Q> operator*(mat<3, 4, T, Q> const &m, Sphere<3, T, Q> const &sphere) {
  return Sphere<3, T, Q>(m * sphere.pos, length(m[0]) * sphere.r);
}

template<typename T, qualifier Q>
static Sphere<3, T, Q> operator*(mat<4, 3, T, Q> const &m, Sphere<3, T, Q> sphere) {
  const T scale = length(vec<3, T, Q>(m[0].x, m[0].y, m[0].z));
  return Sphere<3, T, Q>(transformPos(m, sphere.pos), scale * sphere.r);
}

template<typename T, qualifier Q>
static Sphere<3, T, Q> operator*(mat<4, 4, T, Q> const &m, Sphere<3, T, Q> sphere) {
  const T scale = length(vec<3, T, Q>(m[0].x, m[0].y, m[0].z));
  return Sphere<3, T, Q>(transformPos(m, sphere.pos), scale * sphere.r);
}

template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER Sphere<3, T, Q> operator*(qua<T, Q> const &q, Sphere<3, T, Q> sphere) {
  return Sphere<3, T, Q>(q * sphere.pos, sphere.r);
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool equal(Sphere<L, T, Q> const &x, Sphere<L, T, Q> const &y, T eps = epsilon<T>()) {
  return all(equal(x.pos, y.pos, eps)) && glm::equal(x.r, y.r, eps);
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool equal(Sphere<L, T, Q> const &x, Sphere<L, T, Q> const &y, vec<L, T, Q> const &eps) {
  return all(equal(x.pos, y.pos, eps)) && glm::equal(x.r, y.r, eps[0]);
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool equal(Sphere<L, T, Q> const &x, Sphere<L, T, Q> const &y, int MaxULPs) {
  return all(equal(x.pos, y.pos, MaxULPs)) && glm::equal(x.r, y.r, MaxULPs);
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool equal(Sphere<L, T, Q> const &x, Sphere<L, T, Q> const &y, vec<L, int, Q> const &MaxULPs) {
  return all(equal(x.pos, y.pos, MaxULPs)) && glm::equal(x.r, y.r, MaxULPs[0]);
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool notEqual(Sphere<L, T, Q> const &x, Sphere<L, T, Q> const &y, T eps = epsilon<T>()) {
  return any(notEqual(x.pos, y.pos, eps)) || glm::notEqual(x.r, y.r, eps);
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool notEqual(Sphere<L, T, Q> const &x, Sphere<L, T, Q> const &y, vec<L, T, Q> const &eps) {
  return any(notEqual(x.pos, y.pos, eps)) || glm::notEqual(x.r, y.r, eps[0]);
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool notEqual(Sphere<L, T, Q> const &x, Sphere<L, T, Q> const &y, int MaxULPs) {
  return any(notEqual(x.pos, y.pos, MaxULPs)) || glm::notEqual(x.r, y.r, MaxULPs);
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool notEqual(Sphere<L, T, Q> const &x, Sphere<L, T, Q> const &y, vec<L, int, Q> const &MaxULPs) {
  return any(notEqual(x.pos, y.pos, MaxULPs)) || glm::notEqual(x.r, y.r, MaxULPs[0]);
}

/// Return the center of mass of the Sphere.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> centroid(Sphere<L, T, Q> const &sphere) {
  return sphere.pos;
}

/// Return the smallest AABB that encloses the Sphere.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER AABB<L, T, Q> maximalContainedAABB(Sphere<L, T, Q> const &sphere) {
  AABB<L, T, Q> aabb;
  aabb.setFromCenterAndSize(sphere.pos, vec<L, T, Q>(sphere.r * sqrt(T(3))));
  return aabb;
}

/// Test if any component of the Sphere is infinite
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool isinf(Sphere<L, T, Q> const &sphere) {
  return any(isinf(sphere.pos)) || glm::isinf(sphere.r);
}

/// Test if any component of the sphere is NaN
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool isnan(Sphere<L, T, Q> const &sphere) {
  return any(isnan(sphere.pos)) || glm::isnan(sphere.r);
}

/// Test if all components of the sphere are finite
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool isfinite(Sphere<L, T, Q> const &sphere) {
  return all(isfinite(sphere.pos)) && glm::isfinite(sphere.r);
}

/// Test whether the sphere is degenerate, i.e., not finite or the radius is less-or-equal to zero.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool isDegenerate(Sphere<L, T, Q> const &sphere) {
  return !(sphere.r > T(0) || !all(isfinite(sphere.pos)));
}

/// Return the volume of the Sphere.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T volume(Sphere<L, T, Q> const &sphere) {
  return T(4) * pi<T>() * sphere.r * sphere.r * sphere.r / T(3);
}

/// Return the surface area of the Sphere.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T surfaceArea(Sphere<L, T, Q> const &sphere) {
  return T(4) * pi<T>() * sphere.r * sphere.r;
}

/// Return the area of the Circle.
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER T area(Sphere<2, T, Q> const &sphere) {
  return pi<T>() * sphere.r * sphere.r;
}

/// Return an extreme point along the Sphere, i.e., the furthest point in a given direction.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> extremePoint(Sphere<L, T, Q> const &sphere, vec<L, T, Q> const &direction) {
  const T len = length(direction);
  if (detail::approx_zero(len))
    return sphere.pos;
  return sphere.pos + direction * (sphere.r / len);
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> extremePoint(Sphere<L, T, Q> const &sphere, vec<L, T, Q> const &direction, T &projectionDistance) {
  const vec<L, T, Q> point = extremePoint(sphere, direction);
  projectionDistance = dot(point, direction);
  return point;
}

/// Project the Sphere onto the given axis (direction), i.e., collapse the Sphere onto an axis.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER void projectToAxis(Sphere<L, T, Q> const &sphere, vec<L, T, Q> const &direction, T &outMin, T &outMax) {
  const T d = dot(direction, sphere.pos);
  outMin = d - sphere.r;
  outMax = d + sphere.r;
}

/* Return the point on the sphere closest to the given object */

/// Return point on the Sphere closest to the given point.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Sphere<L, T, Q> const &sphere, vec<L, T, Q> const &point) {
  const T d = distance(sphere.pos, point);
  const T t = (d >= sphere.r ? sphere.r : d);
  return sphere.pos + (point - sphere.pos) * (t / d);
}

#if 0
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Sphere<L, T, Q> const &sphere, AABB<L, T, Q> const &aabb) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Sphere<L, T, Q> const &sphere, Line<L, T, Q> const &line) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Sphere<L, T, Q> const &sphere, Segment<L, T, Q> const &segment) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Sphere<L, T, Q> const &sphere, Ray<L, T, Q> const &ray) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Sphere<L, T, Q> const &sphere, Plane<L, T, Q> const &plane) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Sphere<L, T, Q> const &sphere, Sphere<L, T, Q> const &other) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Sphere<L, T, Q> const &sphere, Capsule<L, T, Q> const &capsule) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Sphere<L, T, Q> const &sphere, Disc<L, T, Q> const &disc) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Sphere<L, T, Q> const &sphere, Triangle<L, T, Q> const &triangle) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Sphere<L, T, Q> const &sphere, Polygon<L, T, Q> const &polygon) { }

template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<3, T, Q> closestPoint(Sphere<3, T, Q> const &sphere, OBB<T, Q> const &obb) { }
#endif

/* Test if the given object is fully contained within the sphere */

/// Sphere/Point containment test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool contains(Sphere<L, T, Q> const &sphere, vec<L, T, Q> const &point, T eps = epsilon<T>()) {
  const T r = sphere.r + eps;
  return distance2(sphere.pos, point) <= r * r;
}

/// Sphere/Segment containment test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool contains(Sphere<L, T, Q> const &sphere, Segment<L, T, Q> const &segment, T eps = epsilon<T>()) {
  return contains(sphere, segment.a, eps) && contains(sphere, segment.b, eps);
}

/// Sphere/Sphere containment test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool contains(Sphere<L, T, Q> const &a, Sphere<L, T, Q> const &b, T eps = epsilon<T>()) {
  return distance2(a.pos, b.pos) + b.r - a.r <= eps;
}

/// Sphere/Capsule containment test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool contains(Sphere<L, T, Q> const &sphere, Capsule<L, T, Q> const &capsule, T eps = epsilon<T>()) {
  return distance(sphere.pos, capsule.l.a) + capsule.r <= (sphere.r + eps)
         && distance(sphere.pos, capsule.l.b) + capsule.r <= (sphere.r + eps);
}

/// Sphere/AABB containment test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool contains(Sphere<L, T, Q> const &sphere, AABB<L, T, Q> const &aabb) {
  bool result = true;
  for (length_t i = 0; i < 8 && result; ++i)
    result &= contains(sphere, cornerPoint(aabb, i));
  return true;
}

/// Sphere/OBB containment test
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool contains(Sphere<3, T, Q> const &sphere, OBB<T, Q> const &obb) {
  bool result = true;
  for (length_t i = 0; i < 8 && result; ++i)
    result &= contains(sphere, cornerPoint(obb, i));
  return true;
}

/// Sphere/Triangle containment test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool contains(Sphere<L, T, Q> const &sphere, Triangle<L, T, Q> const &triangle, T eps = epsilon<T>()) {
  return contains(sphere, triangle.a, eps)
         && contains(sphere, triangle.b, eps)
         && contains(sphere, triangle.c, eps);
}

#if 0
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool contains(Sphere<L, T, Q> const &sphere, Disc<L, T, Q> const &disc) { }
#endif

/* Return the distance between the sphere and the given object */

/// Return the distance between the Sphere and given point
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Sphere<L, T, Q> const &sphere, vec<L, T, Q> const &point) {
  return max(T(0), distance(sphere.pos, point) - sphere.r);
}

/// Return the distance between two Spheres
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Sphere<L, T, Q> const &sphere, Sphere<L, T, Q> const &other) {
  return max(T(0), distance(sphere.pos, other.pos) - sphere.r - other.r);
}

/// Return the distance between the Sphere and given Capsule
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Sphere<L, T, Q> const &sphere, Capsule<L, T, Q> const &capsule) {
  return distance(capsule, sphere);
}

/// Return the distance between the Sphere and given AABB
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Sphere<L, T, Q> const &sphere, AABB<L, T, Q> const &aabb) {
  return distance(aabb, sphere);
}

/// Return the distance between the Sphere and given OBB
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Sphere<3, T, Q> const &sphere, OBB<T, Q> const &obb) {
  return distance(obb, sphere);
}

/// Return the distance between the Sphere and given Ray
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Sphere<L, T, Q> const &sphere, Ray<L, T, Q> const &ray) {
  return distance(ray, sphere);
}

/// Return the distance between the Sphere and given Segment
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Sphere<L, T, Q> const &sphere, Segment<L, T, Q> const &segment) {
  return distance(segment, sphere);
}

/// Return the distance between the Sphere and given Line
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Sphere<L, T, Q> const &sphere, Line<L, T, Q> const &line) {
  return distance(line, sphere);
}

/// Return the distance between the Sphere and given Plane
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Sphere<L, T, Q> const &sphere, Plane<L, T, Q> const &plane) {
  return distance(plane, sphere);
}

/// Return the distance between the Sphere and given Triangle
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Sphere<L, T, Q> const &sphere, Triangle<L, T, Q> const &triangle) {
  return distance(triangle, sphere);
}

#if 0
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Sphere<L, T, Q> const &sphere, Disc<L, T, Q> const &disc) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Sphere<L, T, Q> const &sphere, Polygon<L, T, Q> const &polygon) { }
#endif

/* Test whether the sphere and the given object intersect */

/// @private: Generic Line/Sphere Intersection
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER int intersectLine(Line<L, T, Q> const &line, Sphere<L, T, Q> const &sphere, T &d1, T &d2) {
  const vec<L, T, Q> a = line.pos - sphere.pos;
  const T C = dot(a, a) - (sphere.r * sphere.r);
  const T B = T(2) * dot(a, line.dir);

  T D = B * B - T(4) * C;
  if (D < T(0)) {  // No intersections.
    d1 = std::numeric_limits<T>::infinity();
    d2 = -std::numeric_limits<T>::infinity();
    return 0;
  }
  else if (D < epsilon<T>()) {  // tangent to sphere
    d1 = d2 = -B * T(0.5);
    return 1;
  }
  else {
    D = sqrt(D);
    d1 = (-B - D) * T(0.5);
    d2 = (-B + D) * T(0.5);
    return 2;
  }
}

/// Sphere/Sphere intersection test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Sphere<L, T, Q> const &sphere, Sphere<L, T, Q> const &other) {
  return distance2(sphere.pos, other.pos) <= ((sphere.r + other.r) * (sphere.r + other.r));
}

/// Sphere/Capsule intersection test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Sphere<L, T, Q> const &sphere, Capsule<L, T, Q> const &capsule) {
  return intersects(capsule, sphere);
}

/// Sphere/Line intersection test, includes parametric segment coordinates.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER int intersects(Sphere<L, T, Q> const &sphere, Line<L, T, Q> const &line, T &d1, T &d2) {
  return intersectLine(line, sphere, d1, d2);
}

/// Sphere/Segment intersection test, includes parametric liune coordinates.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER int intersects(Sphere<L, T, Q> const &sphere, Segment<L, T, Q> const &segment, T &d1, T &d2) {
  const int numIntersections = intersectLine(toLine(segment), sphere, d1, d2);
  if (numIntersections == 0)
    return 0;

  const T lineLength = length(segment);
  if (d2 < T(0) || d1 > lineLength)
    return 0;

  d1 /= lineLength;
  d2 /= lineLength;
  return numIntersections;
}

/// Sphere/Ray intersection test, includes parametric ray coordinates.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER int intersects(Sphere<L, T, Q> const &sphere, Ray<L, T, Q> const &ray, T &d1, T &d2) {
  int numIntersections = intersectLine(ray, sphere, d1, d2);
  if (d1 < T(0) && numIntersections == 2)  // behind the ray.
    d1 = d2;
  return (d1 >= T(0)) ? numIntersections : 0;  // otherwise, negative direction of the ray.
}

/// Sphere/AABB intersection test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Sphere<L, T, Q> const &sphere, AABB<L, T, Q> const &aabb) {
  return intersects(aabb, sphere);
}

/// Sphere/OBB intersection test
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Sphere<3, T, Q> const &sphere, OBB<T, Q> const &obb) {
  return intersects(obb, sphere);
}

/// Sphere/Plane intersection test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Sphere<L, T, Q> const &sphere, Plane<L, T, Q> const &plane) {
  return intersects(plane, sphere);
}

/// Sphere/Triangle intersection test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Sphere<L, T, Q> const &sphere, Triangle<L, T, Q> const &triangle) {
  vec<L, T, Q> intersectionPt;
  return intersects(triangle, sphere, intersectionPt);
}

/// Sphere/Line intersection test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER int intersects(Sphere<L, T, Q> const &sphere, Line<L, T, Q> const &line) {
  T d1(0), d2(0);
  return intersectLine(line, sphere, d1, d2);
}

/// Sphere/Ray intersection test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER int intersects(Sphere<L, T, Q> const &sphere, Ray<L, T, Q> const &ray) {
  T d1(0), d2(0);
  return intersects(sphere, ray, d1, d2);
}

/// Sphere/Segment intersection test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER int intersects(Sphere<L, T, Q> const &sphere, Segment<L, T, Q> const &segment) {
  T d1(0), d2(0);
  return intersects(sphere, segment, d1, d2);
}

/// Return the intersection between the Sphere and Plane.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER Disc<L, T, Q> clipPlane(Sphere<L, T, Q> const &sphere, Plane<L, T, Q> const &plane) {
  const vec<L, T, Q> pos = closestPoint(plane, sphere.pos);
  const vec<L, T, Q> normal = plane.normal;
  const T r = sqrt(sphere.r * sphere.r - distance2(pos, sphere.pos));
  return Disc<L, T, Q>(pos, normal, r);
}

/* Expand the sphere to enclose the given object */

/// <summary>
/// Generic enclosing algorithm. Fetch each corner point of the given object and
/// enclose points by the furthest distance from the sphere to shortest.
/// @private
/// </summary>
template<length_t L, typename T, qualifier Q, typename Object, length_t NumCorners>
static Sphere<L, T, Q> encloseObject(Sphere<L, T, Q> const &sphere, Object const &obj) {
  struct Tuple {
    length_t idx;
    vec<L, T, Q> point;
    T distance;
    bool operator<(Tuple const &rhs) const {
      return detail::equal_to(distance, rhs.distance) ? idx < rhs.idx : distance < rhs.distance;
    }
  };

  GLM_STATIC_ASSERT(NumCorners > 0, "NumCorners");
  Tuple corners[static_cast<size_t>(NumCorners)];
  for (length_t i = 0; i < NumCorners; ++i) {
    const vec<L, T, Q> point = cornerPoint(obj, i);
    corners[i] = { i, point, distance2(sphere.pos, point) };
  }

#if 0 && GLM_HAS_CXX11_STL
  std::sort(corners, corners + NumCorners);
#else
  qsort(corners, NumCorners, sizeof(Tuple), [](const void *a, const void *b) -> int {
    const Tuple &t_a = *static_cast<const Tuple *>(a);
    const Tuple &t_b = *static_cast<const Tuple *>(b);
    return (t_a < t_b) ? 1 : ((t_b < t_a) ? -1 : (t_a.idx - t_b.idx));
  });
#endif

  Sphere<L, T, Q> result(sphere);
  for (length_t i = NumCorners; i > 0; --i)
    result.enclose(corners[i - 1].point);
  return result;
}

/// Expand the Sphere to enclose the given point
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER Sphere<L, T, Q> enclose(Sphere<L, T, Q> const &sphere, vec<L, T, Q> const &point, T eps = epsilon<T>()) {
  Sphere<L, T, Q> result(sphere);
  result.enclose(point, eps);
  return result;
}

/// Expand the Sphere to enclose the given Segment
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER Sphere<L, T, Q> enclose(Sphere<L, T, Q> const &sphere, Segment<L, T, Q> const &segment) {
  Sphere<L, T, Q> result(sphere);
  if (distance2(sphere.pos, segment.a) > distance2(sphere.pos, segment.b)) {
    result.enclose(segment.a);
    result.enclose(segment.b);
  }
  else {
    result.enclose(segment.b);
    result.enclose(segment.a);
  }
  return result;
}

/// Expand the Sphere to enclose the given AABB
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER Sphere<3, T, Q> enclose(Sphere<3, T, Q> const &sphere, AABB<3, T, Q> const &aabb) {
  return encloseObject<3, T, Q, AABB<3, T, Q>, 8>(sphere, aabb);  // @TODO: Replace specialization w/ 2^L
}

/// Expand the Sphere to enclose the given AABB
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER Sphere<2, T, Q> enclose(Sphere<2, T, Q> const &sphere, AABB<2, T, Q> const &aabb) {
  return encloseObject<2, T, Q, AABB<2, T, Q>, 4>(sphere, aabb);
}

/// Expand the Sphere to enclose the given OBB
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER Sphere<3, T, Q> enclose(Sphere<3, T, Q> const &sphere, OBB<T, Q> const &obb) {
  return encloseObject<3, T, Q, OBB<T, Q>, 8>(sphere, obb);
}

/// Expand the Sphere to enclose the given Triangle
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER Sphere<L, T, Q> enclose(Sphere<L, T, Q> const &sphere, Triangle<L, T, Q> const &triangle) {
  return encloseObject<L, T, Q, Triangle<L, T, Q>, 3>(sphere, triangle);
}

/// Expand the Sphere to enclose the other Sphere
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER Sphere<L, T, Q> enclose(Sphere<L, T, Q> const &sphere, Sphere<L, T, Q> const &other) {
  const vec<L, T, Q> furthestPoint = scaleLength((other.pos - sphere.pos), other.r);
  Sphere<L, T, Q> result(sphere);
  result.enclose(other.pos + furthestPoint);
  result.enclose(other.pos - furthestPoint);
  return result;
}

/// Expand the Sphere to enclose the given Capsule
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER Sphere<L, T, Q> enclose(Sphere<L, T, Q> const &sphere, Capsule<L, T, Q> const &capsule) {
  const T dista = distance2(sphere.pos, capsule.l.a);
  const T distb = distance2(sphere.pos, capsule.l.b);
  Sphere<L, T, Q> result(sphere);
  if (dista > distb) {
    result = enclose(result, Sphere<L, T, Q>(capsule.l.a, capsule.r));
    result = enclose(result, Sphere<L, T, Q>(capsule.l.b, capsule.r));
  }
  else {
    result = enclose(result, Sphere<L, T, Q>(capsule.l.b, capsule.r));
    result = enclose(result, Sphere<L, T, Q>(capsule.l.a, capsule.r));
  }
  return result;
}

/// Expand the radius of the Sphere until it encloses the given point.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER Sphere<L, T, Q> extendRadiusToContain(Sphere<L, T, Q> const &sphere, vec<L, T, Q> const &point, T eps = epsilon<T>()) {
  const T radius = distance(sphere.pos, point) + eps;
  return Sphere<L, T, Q>(sphere.pos, max(sphere.r, radius));
}

/// Expand the radius of the Sphere until it encloses the given Sphere.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER Sphere<L, T, Q> extendRadiusToContain(Sphere<L, T, Q> const &sphere, Sphere<L, T, Q> const &other, T eps = epsilon<T>()) {
  const T radius = distance(sphere.pos, other.pos) + other.r + eps;
  return Sphere<L, T, Q>(sphere.pos, max(sphere.r, radius));
}

/// <summary>
/// Return the minimal bounding Sphere for two points, i.e., the smallest volume
/// Sphere that contains the provided points.
/// </summary>
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER Sphere<3, T, Q> optimalEnclosingSphere(vec<3, T, Q> const &a, vec<3, T, Q> const &b, T eps = epsilon<T>()) {
  const vec<3, T, Q> pos = (a + b) * T(0.5);
  if (all(isfinite(pos)))
    return Sphere<3, T, Q>(pos, length(b - pos) + eps);
  return Sphere<3, T, Q>(vec<3, T, Q>(T(0)), T(0));
}

/// <summary>
/// Return the minimal bounding sphere for three points, i.e., the smallest
/// volume Sphere that contains the provided points.
/// </summary>
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER_OUTLINE Sphere<3, T, Q> optimalEnclosingSphere(vec<3, T, Q> const &a, vec<3, T, Q> const &b, vec<3, T, Q> const &c, T eps = epsilon<T>()) {
  Sphere<3, T, Q> sphere;
  const vec<3, T, Q> ab = b - a;
  const vec<3, T, Q> ac = c - a;

  T s(0), t(0);
  const bool success = !areCollinear(ab, ac, epsilon<T>()) && fitSphereThroughPoints(ab, ac, s, t);

  // Points are either collinear or sufficiently far enough away from <a, b, c>
  if (!success || abs(s) > T(10000) || abs(t) > T(10000)) {
    const vec<3, T, Q> minPt = min(a, min(b, c));
    const vec<3, T, Q> maxPt = max(a, max(b, c));
    sphere.pos = (minPt + maxPt) * T(0.5);
    sphere.r = distance(sphere.pos, minPt);
  }
  else if (s < T(0)) {
    sphere.pos = (a + c) * T(0.5);
    sphere.r = max(distance(a, c) * T(0.5), distance(b, sphere.pos));
  }
  else if (t < T(0)) {
    sphere.pos = (a + b) * T(0.5);
    sphere.r = max(distance(a, b) * T(0.5), distance(c, sphere.pos));
  }
  else if (s + t > T(1)) {
    sphere.pos = (b + c) * T(0.5);
    sphere.r = max(distance(b, c) * T(0.5), distance(a, sphere.pos));
  }
  else {
    sphere.pos = a + s * ab + t * ac;
    sphere.r = sqrt(max(distance2(sphere.pos, a), max(distance2(sphere.pos, b), distance2(sphere.pos, c))));
  }

  sphere.r += T(2) * eps;
  return sphere;
}

/// <summary>
/// Return the minimal bounding Sphere for four points, i.e., the smallest
/// volume Sphere that contains the provided points.
/// </summary>
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER_OUTLINE Sphere<3, T, Q> optimalEnclosingSphere(vec<3, T, Q> const &a, vec<3, T, Q> const &b, vec<3, T, Q> const &c, vec<3, T, Q> const &d, T eps = epsilon<T>()) {
  T s(0), t(0), u(0);
  Sphere<3, T, Q> sphere;

  const vec<3, T, Q> ab = b - a;
  const vec<3, T, Q> ac = c - a;
  const vec<3, T, Q> ad = d - a;

  bool success = fitSphereThroughPoints(ab, ac, ad, s, t, u);
  if (!success || s < T(0) || t < T(0) || u < T(0) || s + t + u > T(1)) {
    sphere = optimalEnclosingSphere(a, b, c, eps);
    if (!contains(sphere, d)) {
      sphere = optimalEnclosingSphere(a, b, d, eps);
      if (!contains(sphere, c)) {
        sphere = optimalEnclosingSphere(a, c, d, eps);
        if (!contains(sphere, b)) {
          sphere = optimalEnclosingSphere(b, c, d, eps);
          sphere.r = max(sphere.r, distance(a, sphere.pos) + eps);
        }
      }
    }
  }
  else {
    sphere.pos = a + s * ab + t * ac + u * ad;
    sphere.r = sqrt(max(distance2(sphere.pos, a), max(
      distance2(sphere.pos, b),
      max(distance2(sphere.pos, c), distance2(sphere.pos, d))
    )));
  }

  sphere.r += T(2) * eps;
  return sphere;
}

/// <summary>
/// Return the minimal bounding Sphere for five points, i.e., the smallest
/// volume sphere that contains the points.
///
/// A minimal enclosing sphere can be defined by four points (or fewer). At
/// least one of the provided points is considered redundant.
/// </summary>
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER_OUTLINE Sphere<3, T, Q> optimalEnclosingSphere(vec<3, T, Q> const &a, vec<3, T, Q> const &b, vec<3, T, Q> const &c, vec<3, T, Q> const &d, vec<3, T, Q> const &e, size_t &redundantPoint, T eps = epsilon<T>()) {
  Sphere<3, T, Q> s = optimalEnclosingSphere(b, c, d, e, eps);
  if (contains(s, a, eps))
    return (redundantPoint = 0, s);
  else if (contains((s = optimalEnclosingSphere(a, c, d, e, eps)), b, eps))
    return (redundantPoint = 1, s);
  else if (contains((s = optimalEnclosingSphere(a, b, d, e, eps)), c, eps))
    return (redundantPoint = 2, s);
  else if (contains((s = optimalEnclosingSphere(a, b, c, e, eps)), d, eps))
    return (redundantPoint = 3, s);
  else {
    s = optimalEnclosingSphere(a, b, c, d, eps);
    return (redundantPoint = 4, s);
  }
}

template<typename T, qualifier Q, class Vector>
GLM_GEOM_QUALIFIER_OUTLINE Sphere<3, T, Q> optimalEnclosingSphere(Vector const &pts, T eps = epsilon<T>()) {
  switch (pts.size()) {
    case 0: return Sphere<3, T, Q>(T(0));
    case 1: return Sphere<3, T, Q>(pts[0], T(0));
    case 2: return optimalEnclosingSphere<T, Q>(pts[0], pts[1], eps);
    case 3: return optimalEnclosingSphere<T, Q>(pts[0], pts[1], pts[2], eps);
    case 4: return optimalEnclosingSphere<T, Q>(pts[0], pts[1], pts[2], pts[3], eps);
    default: {
      break;
    }
  }

  // The set of supporting points for the minimal sphere.
  typename Vector::size_type sp[4] = { 0, 1, 2, 3 };
  bool expendable[4] = { true, true, true, true };

  Sphere<3, T, Q> s = optimalEnclosingSphere(pts[sp[0]], pts[sp[1]], pts[sp[2]], pts[sp[3]], eps);
  T r2 = s.r * s.r + eps;
  for (typename Vector::size_type i = 4; i < pts.size(); ++i) {
    if (i == sp[0] || i == sp[1] || i == sp[2] || i == sp[3])
      continue;

    // If the next point does not fit inside the current sphere, compute the
    // minimal sphere that contains it.
    if (distance2(pts[i], s.pos) > r2) {
      size_t redundant = 0;
      s = optimalEnclosingSphere(pts[sp[0]], pts[sp[1]], pts[sp[2]], pts[sp[3]], pts[i], redundant, eps);
      r2 = s.r * s.r + eps;

      // One of the five points is now redundant and can be removed.
      if (redundant != 4 && (sp[redundant] < i || expendable[redundant])) {
        sp[redundant] = i;
        expendable[redundant] = false;
        if (sp[0] < i) expendable[0] = true;
        if (sp[1] < i) expendable[1] = true;
        if (sp[2] < i) expendable[2] = true;
        if (sp[3] < i) expendable[3] = true;
        i = 0;
      }
    }
  }

  return s;
}

/// <summary>
/// Return true if a sphere can be fit through zero, 'ab', and 'ac'. Storing its
/// barycentric coordinates in 's', and 't'.
/// </summary>
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool fitSphereThroughPoints(vec<3, T, Q> const &ab, vec<3, T, Q> const &ac, T &s, T &t) {
  const T bb = dot(ab, ab);
  const T cc = dot(ac, ac);
  const T bc = dot(ab, ac);

  T denom = bb * cc - bc * bc;
  if (!detail::approx_zero(denom)) {
    denom = T(0.5) / denom;
    s = (cc * bb - bc * cc) * denom;
    t = (cc * bb - bc * bb) * denom;
    return true;
  }
  return false;
}

/// <summary>
/// Return true if a sphere can be fit through zero, 'ab', 'ac', and 'ad'.
/// Storing its barycentric coordinates in 's', 't', and 'u'.
/// </summary>
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER_OUTLINE bool fitSphereThroughPoints(vec<3, T, Q> const &ab, vec<3, T, Q> const &ac, vec<3, T, Q> const &ad, T &s, T &t, T &u) {
  const T bb = dot(ab, ab);
  const T bc = dot(ab, ac);
  const T bd = dot(ab, ad);
  const T cc = dot(ac, ac);
  const T cd = dot(ac, ad);
  const T dd = dot(ad, ad);
  const mat<3, 3, T, Q> ms(bb, bc, bd, bc, cc, cd, bd, cd, dd);
  if (invertible(ms)) {
    const mat<3, 3, T, Q> m = inverse(ms);
    const vec<3, T, Q> v = m * vec<3, T, Q>(bb * T(0.5), cc * T(0.5), dd * T(0.5));
    s = v.x;
    t = v.y;
    u = v.z;
    return true;
  }
  return false;
}

/// Fit a Sphere through the two given points.
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER Sphere<3, T, Q> fitThroughPoints(vec<3, T, Q> const &a, vec<3, T, Q> const &b) {
  return optimalEnclosingSphere(a, b);
}

/// <summary>
/// Fit a Sphere through three points. Returning the Sphere that contains a, b,
/// and c while minimizing its volume.
/// </summary>
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER Sphere<3, T, Q> fitThroughPoints(vec<3, T, Q> const &a, vec<3, T, Q> const &b, vec<3, T, Q> const &c) {
  Sphere<3, T, Q> sphere;
  const vec<3, T, Q> ab = b - a;
  const vec<3, T, Q> ac = c - a;

  T s(0), t(0);
  if (fitSphereThroughPoints(ab, ac, s, t)) {
    const vec<3, T, Q> center = s * ab + t * ac;
    sphere.r = length(center);
    sphere.pos = a + center;
  }
  else {
    sphere.setDegenerate();
  }
  return sphere;
}

/// <summary>
/// Fit a Sphere through four points. Four non-coplanar points can uniquely
/// define a Sphere in three dimensions.
/// </summary>
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER Sphere<3, T, Q> fitThroughPoints(vec<3, T, Q> const &a, vec<3, T, Q> const &b, vec<3, T, Q> const &c, vec<3, T, Q> const &d) {
  Sphere<3, T, Q> sphere;
  T s(0), t(0), u(0);
  const vec<3, T, Q> ab = b - a;
  const vec<3, T, Q> ac = c - a;
  const vec<3, T, Q> ad = d - a;
  if (fitSphereThroughPoints(ab, ac, ad, s, t, u)) {
    const vec<3, T, Q> center = s * ab + t * ac + u * ad;
    sphere.r = length(center);
    sphere.pos = a + center;
  }
  else {
    sphere.setDegenerate();
  }
  return sphere;
}

/// @}
GLM_GEOM_END_NAMESPACE
#if GLM_GEOM_TOSTRING
namespace glm {
  namespace detail {
    template<glm::length_t L, typename T, qualifier Q>
    struct compute_to_string<Sphere<L, T, Q>> {
      GLM_GEOM_QUALIFIER std::string call(Sphere<L, T, Q> const &sphere) {
        char const *LiteralStr = literal<T, std::numeric_limits<T>::is_iec559>::value();
        std::string FormatStr(detail::format("sphere(%%s, %s)", LiteralStr));
        return detail::format(FormatStr.c_str(),
          glm::to_string(sphere.pos).c_str(),
          static_cast<typename cast<T>::value_type>(sphere.r)
        );
      }
    };
  }
}
#endif
