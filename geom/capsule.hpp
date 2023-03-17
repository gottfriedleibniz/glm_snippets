/// @ref geom_capsule
/// @file geom/capsule.hpp
///
/// @defgroup geom_capsule Capsule
/// @ingroup geom
///

#pragma once

#include "setup.hpp"
#include "aabb.hpp"
#include "disc.hpp"
#include "segment.hpp"

#include "../vector_extensions.hpp"
#include "../quat_extensions.hpp"
#if GLM_MESSAGES == GLM_ENABLE && !defined(GLM_EXT_INCLUDED)
  #pragma message("GLM: GLM_EXT_GEOM_aabb extension included")
#endif

GLM_GEOM_BEGIN_NAMESPACE
/// @addtogroup geom_capsule
/// @{

/// <summary>
/// A set of points equidistant (by some radius) from a segment.
/// </summary>
template<length_t L, typename T, qualifier Q>
struct Capsule {

  // -- Implementation detail --

  typedef T value_type;
  typedef Capsule<L, T, Q> type;
  typedef vec<L, T, Q> point_type;
  typedef Segment<L, T, Q> line_type;

  // -- Data --

  line_type l;  // Coordinates of this cylinder: { bottom, top }
  value_type r;  // Radius of this capsule

#if GLM_CONFIG_DEFAULTED_DEFAULT_CTOR == GLM_ENABLE
  Capsule() GLM_DEFAULT_CTOR;
#else
  Capsule()
  #if GLM_CONFIG_CTOR_INIT != GLM_CTOR_INIT_DISABLE
    : l(T(0)), r(T(0))
  #endif
  {
  }
#endif

  Capsule(T scalar)
    : l(scalar), r(scalar) {
  }

  Capsule(line_type const &line, value_type radius)
    : l(line), r(radius) {
  }

  Capsule(point_type const &bottom, point_type const &top, value_type radius)
    : l(bottom, top), r(radius) {
  }

  Capsule(Capsule<L, T, Q> const &capsule)
    : l(capsule.l), r(capsule.r) {
  }

  Capsule<L, T, Q> &operator=(Capsule<L, T, Q> const &capsule) {
    l = capsule.l;
    r = capsule.r;
    return *this;
  }
};

template<length_t L, typename T, qualifier Q>
static Capsule<L, T, Q> operator-(Capsule<L, T, Q> const &capsule) {
  return Capsule<L, T, Q>(capsule.l, capsule.r);
}

template<length_t L, typename T, qualifier Q>
static bool operator==(Capsule<L, T, Q> const &s1, Capsule<L, T, Q> const &s2) {
  return operator==(s1.l, s2.l) && detail::equal_to(s1.r, s2.r);
}

template<length_t L, typename T, qualifier Q>
static Capsule<L, T, Q> operator+(Capsule<L, T, Q> const &capsule, vec<L, T, Q> const &offset) {
  return Capsule<L, T, Q>(capsule.l + offset, capsule.r);
}

template<length_t L, typename T, qualifier Q>
static Capsule<L, T, Q> operator-(Capsule<L, T, Q> const &capsule, vec<L, T, Q> const &offset) {
  return Capsule<L, T, Q>(capsule.l - offset, capsule.r);
}

template<typename T, qualifier Q>
static Capsule<3, T, Q> operator*(mat<3, 3, T, Q> const &m, Capsule<3, T, Q> const &capsule) {
  return Capsule<3, T, Q>(m * capsule.l, length(m[0]) * capsule.r);
}

template<typename T, qualifier Q>
static Capsule<3, T, Q> operator*(mat<3, 4, T, Q> const &m, Capsule<3, T, Q> const &capsule) {
  return Capsule<3, T, Q>(m * capsule.l, length(m[0]) * capsule.r);
}

template<typename T, qualifier Q>
static Capsule<3, T, Q> operator*(mat<4, 3, T, Q> const &m, Capsule<3, T, Q> capsule) {
  const T scale = length(vec<3, T, Q>(m[0].x, m[0].y, m[0].z));
  return Capsule<3, T, Q>(m * capsule.l, scale * capsule.r);
}

template<typename T, qualifier Q>
static Capsule<3, T, Q> operator*(mat<4, 4, T, Q> const &m, Capsule<3, T, Q> capsule) {
  const T scale = length(vec<3, T, Q>(m[0].x, m[0].y, m[0].z));
  return Capsule<3, T, Q>(m * capsule.l, scale * capsule.r);
}

template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER Capsule<3, T, Q> operator*(qua<T, Q> const &q, Capsule<3, T, Q> capsule) {
  return Capsule<3, T, Q>(q * capsule.l, capsule.r);
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool equal(Capsule<L, T, Q> const &x, Capsule<L, T, Q> const &y, T eps = epsilon<T>()) {
  return equal(x.l, y.l, eps) && glm::equal(x.r, y.r, eps);
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool equal(Capsule<L, T, Q> const &x, Capsule<L, T, Q> const &y, vec<L, T, Q> const &eps) {
  return equal(x.l, y.l, eps) && glm::equal(x.r, y.r, eps[0]);
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool equal(Capsule<L, T, Q> const &x, Capsule<L, T, Q> const &y, int MaxULPs) {
  return equal(x.l, y.l, MaxULPs) && glm::equal(x.r, y.r, MaxULPs);
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool equal(Capsule<L, T, Q> const &x, Capsule<L, T, Q> const &y, vec<L, int, Q> const &MaxULPs) {
  return equal(x.l, y.l, MaxULPs) && glm::equal(x.r, y.r, MaxULPs[0]);
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool notEqual(Capsule<L, T, Q> const &x, Capsule<L, T, Q> const &y, T eps = epsilon<T>()) {
  return notEqual(x.l, y.l, eps) || glm::notEqual(x.r, y.r, eps);
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool notEqual(Capsule<L, T, Q> const &x, Capsule<L, T, Q> const &y, vec<L, T, Q> const &eps) {
  return notEqual(x.l, y.l, eps) || glm::notEqual(x.r, y.r, eps[0]);
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool notEqual(Capsule<L, T, Q> const &x, Capsule<L, T, Q> const &y, int MaxULPs) {
  return notEqual(x.l, y.l, MaxULPs) || glm::notEqual(x.r, y.r, MaxULPs);
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool notEqual(Capsule<L, T, Q> const &x, Capsule<L, T, Q> const &y, vec<L, int, Q> const &MaxULPs) {
  return notEqual(x.l, y.l, MaxULPs) || glm::notEqual(x.r, y.r, MaxULPs[0]);
}

/// Test if any component of the capsule is infinite
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool isinf(Capsule<L, T, Q> const &capsule) {
  return isinf(capsule.l) || glm::isinf(capsule.r);
}

/// Test if any component of the capsule is NaN
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool isnan(Capsule<L, T, Q> const &capsule) {
  return isnan(capsule.l) || glm::isnan(capsule.r);
}

/// Test if all components of the capsule are finite
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool isfinite(Capsule<L, T, Q> const &capsule) {
  return isfinite(capsule.l) && glm::isfinite(capsule.r);
}

/// Test whether the capsule is degenerate, i.e., not finite or if the normal is null
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool isDegenerate(Capsule<L, T, Q> const &capsule) {
  return !isfinite(capsule) || capsule.r <= T(0);
}

/// Return the center of mass of the Capsule.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> centroid(Capsule<L, T, Q> const &capsule) {
  return centerPoint(capsule.l);
}

/// Return the diameter of the Capsule
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T diameter(Capsule<L, T, Q> const &capsule) {
  return T(2) * capsule.r;
}

/// Return the total height of the Capsule.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T height(Capsule<L, T, Q> const &capsule) {
  return length(capsule.l) + diameter(capsule);
}

/// Return the height of only the cylindrical section.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T heightCylinder(Capsule<L, T, Q> const &capsule) {
  return length(capsule.l);
}

/// Return the volume of the Capsule
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T volume(Capsule<L, T, Q> const &capsule) {
  const T r = capsule.r;
  return (pi<T>() * r * r * length(capsule.l)) + (T(4) * pi<T>() * r * r * r / T(3));
}

/// Return the surface area of the Capsule
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T surfaceArea(Capsule<L, T, Q> const &capsule) {
  const T r = capsule.r;
  return (T(2) * pi<T>() * r * length(capsule.l)) + (T(4) * pi<T>() * r * r);
}

/// Return the top-most point of this Capsule, i.e., the extreme point along its normal.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> top(Capsule<L, T, Q> const &capsule) {
  return capsule.l.b + capsule.l.dir() * capsule.r;
}

/// Return the bottom-most point of this Capsule, i.e., the extreme point along its reverse normal.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> bottom(Capsule<L, T, Q> const &capsule) {
  return capsule.l.a - capsule.l.dir() * capsule.r;
}

/// Return the Segment corresponding to the top/bottom coordinates of the Capsule
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER Segment<L, T, Q> segment(Capsule<L, T, Q> const &capsule) {
  return Segment<L, T, Q>(bottom(capsule), top(capsule));
}

/// Return an extreme point along the Capsule, i.e., the furthest point in a given direction
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> extremePoint(Capsule<L, T, Q> const &capsule, vec<L, T, Q> const &direction) {
  const T len = length(direction);
  if (!detail::exactly_zero(len)) {
    const vec<L, T, Q> pos = dot(direction, capsule.l.b - capsule.l.a) >= T(0) ? capsule.l.b : capsule.l.a;
    return pos + direction * (capsule.r / len);
  }
  return capsule.l.b;
}

/// Project the Capsule onto the provided axis
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER void projectToAxis(Capsule<L, T, Q> const &capsule, vec<L, T, Q> const &direction, T &outMin, T &outMax) {
  const T la = dot(direction, capsule.l.a);
  const T lb = dot(direction, capsule.l.b);
  outMin = min(la, lb) - capsule.r;
  outMax = max(la, lb) + capsule.r;
#if 0
  outMin = min(la, lb) - capsule.r * length(direction);
  outMax = max(la, lb) + capsule.r * length(direction);
#endif
}

/// Return the cross-section of this Capsule at a relative, i.e., between 0 and 1 inclusive, height/offset.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER Disc<L, T, Q> crossSection(Capsule<L, T, Q> const &capsule, T offset) {
  const T r = capsule.r;
  const T pos = offset * height(capsule);
  const vec<L, T, Q> up = capsule.l.dir();
  const vec<L, T, Q> centerPos = bottom(capsule) + up * pos;
  if (offset < r)
    return Disc<L, T, Q>(centerPos, up, sqrt(r * r - (r - offset) * (r - offset)));
  else if (offset < length(capsule.l) + r)
    return Disc<L, T, Q>(centerPos, up, r);
  else {
    const T d = offset - r - length(capsule.l);
    return Disc<L, T, Q>(centerPos, up, sqrt(r * r - d * d));
  }
}

/// Generate a point inside this Capsule at a given height, angle, and distance.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> pointInside(Capsule<L, T, Q> const &capsule, T height, T angle, T dist) {
  return getPoint(crossSection(capsule, height), T(2) * pi<T>() * angle, dist);
}

/// Return the smallest AABB that encloses the Capsule
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER AABB<L, T, Q> minimalEnclosingAABB(Capsule<L, T, Q> const &capsule) {
  const vec<L, T, Q> vec(capsule.r);
  return AABB<L, T, Q>(min(capsule.l.a, capsule.l.b) - vec, max(capsule.l.a, capsule.l.b) + vec);
}

/// Return the smallest OBB that encloses the Capsule
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER OBB<T, Q> minimalEnclosingOBB(Capsule<3, T, Q> const &capsule) {
  const vec<3, T, Q> basisX = capsule.l.dir();
  vec<3, T, Q> basisY, basisZ;
  perpendicularBasis(basisX, basisY, basisZ);
  return OBB<T, Q>(
    centroid(capsule),
    fromBasis(basisX, basisY, basisZ),
    vec<3, T, Q>(T(0.5) * height(capsule))
  );
}

/* Return the closest point between an capsule and the given object */

/// Return the closest point between an Capsule and the given object
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Capsule<L, T, Q> const &capsule, vec<L, T, Q> const &point) {
  const vec<L, T, Q> closest = closestPoint(capsule.l, point);
  if (distance2(closest, point) <= capsule.r * capsule.r)
    return point;
  return closest + scaleLength(point - closest, capsule.r);
}

#if 0
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Capsule<L, T, Q> const &capsule, AABB<L, T, Q> const &aabb) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Capsule<L, T, Q> const &capsule, Line<L, T, Q> const &line) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Capsule<L, T, Q> const &capsule, Segment<L, T, Q> const &segment) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Capsule<L, T, Q> const &capsule, Ray<L, T, Q> const &ray) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Capsule<L, T, Q> const &capsule, Plane<L, T, Q> const &plane) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Capsule<L, T, Q> const &capsule, Triangle<L, T, Q> const &triangle) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Capsule<L, T, Q> const &capsule, Sphere<L, T, Q> const &sphere) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Capsule<L, T, Q> const &capsule, Capsule<L, T, Q> const &other) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Capsule<L, T, Q> const &capsule, Disc<L, T, Q> const &disc) { }

template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<3, T, Q> closestPoint(Capsule<3, T, Q> const &capsule, OBB<T, Q> const &obb) { }
#endif

/* Return the distance between the capsule and the given object. */

/// Return the distance between the Capsule and given point
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Capsule<L, T, Q> const &capsule, vec<L, T, Q> const &point) {
  return max(T(0), distance(capsule.l, point) - capsule.r);
}

/// Return the distance between the Capsule and given P lane
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Capsule<L, T, Q> const &capsule, Plane<L, T, Q> const &plane) {
  return distance(plane, capsule);
}

/// Return the distance between the Capsule and given Sphere
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Capsule<L, T, Q> const &capsule, Sphere<L, T, Q> const &sphere) {
  return max(T(0), distance(capsule, sphere.pos) - sphere.r);
}

/// Return the distance between two Capsules
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Capsule<L, T, Q> const &capsule, Capsule<L, T, Q> const &other) {
  return max(T(0), distance(capsule.l, other.l) - capsule.r - other.r);
}

/// Return the distance between the Capsule and given Line
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Capsule<L, T, Q> const &capsule, Line<L, T, Q> const &line) {
  return distance(line, capsule);
}

/// Return the distance between the Capsule and given Segment
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Capsule<L, T, Q> const &capsule, Segment<L, T, Q> const &segment) {
  return distance(segment, capsule);
}

/// Return the distance between the Capsule and given Ray
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Capsule<L, T, Q> const &capsule, Ray<L, T, Q> const &ray) {
  return distance(ray, capsule);
}

#if 0
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Capsule<L, T, Q> const &capsule, Disc<L, T, Q> const &disc) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Capsule<L, T, Q> const &capsule, AABB<L, T, Q> const &aabb) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Capsule<L, T, Q> const &capsule, Triangle<L, T, Q> const &triangle) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Capsule<L, T, Q> const &capsule, Polygon<L, T, Q> const &polygon) { }
#endif

/* Test if the given objects are fully contained inside the capsule */

/// Capsule/Point containment test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool contains(Capsule<L, T, Q> const &capsule, vec<L, T, Q> const &point, T eps = epsilon<T>()) {
  return distance(capsule.l, point) <= capsule.r + eps;
}

/// Capsule/AABB containment test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool contains(Capsule<L, T, Q> const &capsule, AABB<L, T, Q> const &aabb) {
  bool result = true;
  for (length_t i = 0; i < 8 && result; ++i)
    result &= contains(capsule, cornerPoint(aabb, i));
  return result;
}

/// Capsule/OBB containment test
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool contains(Capsule<3, T, Q> const &capsule, OBB<T, Q> const &obb) {
  bool result = true;
  for (length_t i = 0; i < 8 && result; ++i)
    result &= contains(capsule, cornerPoint(obb, i));
  return result;
}

/// Capsule/Segment containment test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool contains(Capsule<L, T, Q> const &capsule, Segment<L, T, Q> const &segment) {
  return contains(capsule, segment.a) && contains(capsule, segment.b);
}

/// Capsule/Triangle containment test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool contains(Capsule<L, T, Q> const &capsule, Triangle<L, T, Q> const &triangle) {
  return contains(capsule, triangle.a) && contains(capsule, triangle.b) && contains(capsule, triangle.c);
}

/// Capsule/Polygon containment test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool contains(Capsule<L, T, Q> const &capsule, Polygon<L, T, Q> const &polygon) {
  bool result = true;
  const size_t size = polygon.size();
  for (size_t i = 0; i < size && result; ++i)
    result &= contains(capsule, polygon[i]);
  return result;
}

/* Test whether the capsule and the given object intersect */

/// Capsule/Line intersection test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Capsule<L, T, Q> const &capsule, Line<L, T, Q> const &line, T eps = epsilon<T>()) {
  return distance(capsule.l, line) <= capsule.r + eps;
}

/// Capsule/Segment intersection test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Capsule<L, T, Q> const &capsule, Segment<L, T, Q> const &segment, T eps = epsilon<T>()) {
  return distance(capsule.l, segment) <= capsule.r + eps;
}

/// Capsule/Ray intersection test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Capsule<L, T, Q> const &capsule, Ray<L, T, Q> const &ray, T eps = epsilon<T>()) {
  return distance(capsule.l, ray) <= capsule.r + eps;
}

/// Capsule/Sphere intersection test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Capsule<L, T, Q> const &capsule, Sphere<L, T, Q> const &sphere, T eps = epsilon<T>()) {
  const T r = capsule.r + sphere.r + eps;
  return distance2(capsule.l, sphere.pos) <= r * r;
}

/// Capsule/Plane intersection test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Capsule<L, T, Q> const &capsule, Plane<L, T, Q> const &plane, T eps = epsilon<T>()) {
  return distance(capsule.l, plane) <= capsule.r + eps;
}

/// Capsule/Capsule intersection test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Capsule<L, T, Q> const &capsule, Capsule<L, T, Q> const &other, T eps = epsilon<T>()) {
  const T r = capsule.r + other.r + eps;
  return distance2(capsule.l, other.l) <= r * r;
}

/// Capsule/Triangle intersection test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Capsule<L, T, Q> const &capsule, Triangle<L, T, Q> const &triangle, T eps = epsilon<T>()) {
  T d(0);
  const T r = capsule.r + eps;
  const vec<L, T, Q> point = closestPoint(triangle, capsule.l, d);
  return distance2(point, getPoint(capsule.l, d)) <= r * r;
}

#if 0
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Capsule<L, T, Q> const &capsule, AABB<L, T, Q> const &aabb) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Capsule<L, T, Q> const &capsule, Polygon<L, T, Q> const &polygon) { }

template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Capsule<3, T, Q> const &capsule, OBB<T, Q> const &obb) { }
#endif

/// @}
GLM_GEOM_END_NAMESPACE
#if GLM_GEOM_TOSTRING
namespace glm {
  namespace detail {
    template<glm::length_t L, typename T, qualifier Q>
    struct compute_to_string<Capsule<L, T, Q>> {
      GLM_GEOM_QUALIFIER std::string call(Capsule<L, T, Q> const &capsule) {
        char const *LiteralStr = literal<T, std::numeric_limits<T>::is_iec559>::value();
        std::string FormatStr(detail::format("capsule(%%s, %s)", LiteralStr));
        return detail::format(FormatStr.c_str(),
          glm::to_string(capsule.l).c_str(),
          static_cast<typename cast<T>::value_type>(capsule.r)
        );
      }
    };
  }
}
#endif
