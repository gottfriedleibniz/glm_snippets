/// @ref geom_plane
/// @file geom/plane.hpp
///
/// @defgroup geom_plane Plane
/// @ingroup geom
///

#pragma once

#include "setup.hpp"
#include "line.hpp"
#include "segment.hpp"
#include "ray.hpp"
#include "capsule.hpp"
#include "disc.hpp"

#include <glm/gtc/matrix_inverse.hpp>
#include <glm/gtx/projection.hpp>
#include <glm/gtx/quaternion.hpp>

#include "../matrix_extensions.hpp"
#if GLM_MESSAGES == GLM_ENABLE && !defined(GLM_EXT_INCLUDED)
  #pragma message("GLM: GLM_EXT_GEOM_plane extension included")
#endif

GLM_GEOM_BEGIN_NAMESPACE
/// @addtogroup geom_plane
/// @{

/// <summary>
/// A hyperplane represented by dot(N, X) == C, where N is a normal vector, C is
/// the hyperplane constant, and X is any point on the plane. A halfspace is
/// represented by dot(N, X) >= C or <= C.
/// </summary>
template<length_t L, typename T, qualifier Q>
struct Plane {

  // -- Implementation detail --

  typedef T value_type;
  typedef Plane<L, T, Q> type;
  typedef vec<L, T, Q> point_type;

  // -- Data --

  point_type normal;  // Direction of this plane
  T d;  // Offset of this plane

#if GLM_CONFIG_DEFAULTED_DEFAULT_CTOR == GLM_ENABLE
  Plane() GLM_DEFAULT_CTOR;
#else
  Plane()
  #if GLM_CONFIG_CTOR_INIT != GLM_CTOR_INIT_DISABLE
    : normal(T(0)), d(T(0))
  #endif
  {
  }
#endif

  Plane(T scalar)
    : normal(scalar), d(scalar) {
  }

  Plane(point_type const &direction, T offset)
    : normal(direction), d(offset) {
  }

  Plane(point_type const &point, point_type const &normal_)
    : normal(normal_), d(dot(point, normal)) {
  }

  Plane(Plane<L, T, Q> const &plane)
    : normal(plane.normal), d(plane.d) {
  }

  Plane<L, T, Q> &operator=(Plane<L, T, Q> const &plane) {
    normal = plane.normal;
    d = plane.d;
    return *this;
  }
};

template<length_t L, typename T, qualifier Q>
static Plane<L, T, Q> operator-(Plane<L, T, Q> const &plane) {
  return Plane<L, T, Q>(-plane.normal, plane.d);
}

template<length_t L, typename T, qualifier Q>
static bool operator==(Plane<L, T, Q> const &p1, Plane<L, T, Q> const &p2) {
  return p1.normal == p2.normal && detail::equal_to(p1.d, p2.d);
}

template<length_t L, typename T, qualifier Q>
static Plane<L, T, Q> operator+(Plane<L, T, Q> const &plane, vec<L, T, Q> const &offset) {
  return Plane<L, T, Q>(plane.normal, plane.d - dot(plane.normal, offset));
}

template<length_t L, typename T, qualifier Q>
static Plane<L, T, Q> operator-(Plane<L, T, Q> const &plane, vec<L, T, Q> const &offset) {
  return Plane<L, T, Q>(plane.normal, plane.d + dot(plane.normal, offset));
}

template<typename T, qualifier Q>
static Plane<3, T, Q> operator*(mat<3, 3, T, Q> const &m, Plane<3, T, Q> const &plane) {
  const mat<3, 3, T, Q> r = inverseTranspose(m);
  return Plane<3, T, Q>(plane.normal * r, plane.d);
}

template<typename T, qualifier Q>
static Plane<3, T, Q> operator*(mat<3, 4, T, Q> const &m, Plane<3, T, Q> const &plane) {
  const mat<3, 3, T, Q> r(inverse(mat<3, 3, T, Q>(m)));
  return Plane<3, T, Q>(plane.normal * r, plane.d);
}

template<typename T, qualifier Q>
static Plane<3, T, Q> operator*(mat<4, 3, T, Q> const &m, Plane<3, T, Q> const &plane) {
  const mat<3, 3, T, Q> r(inverse(mat<3, 3, T, Q>(m)));
  return Plane<3, T, Q>(plane.normal * r, plane.d + dot(plane.normal, r * m[3]));
}

template<typename T, qualifier Q>
static Plane<3, T, Q> operator*(mat<4, 4, T, Q> const &m, Plane<3, T, Q> const &plane) {
  const mat<3, 3, T, Q> r(inverse(mat<3, 3, T, Q>(m)));
  return Plane<3, T, Q>(plane.normal * r, plane.d + dot(plane.normal, r * m[3]));
}

template<typename T, qualifier Q>
static Plane<3, T, Q> operator*(qua<T, Q> const &q, Plane<3, T, Q> const &plane) {
  return operator*(mat3_cast(q), plane);
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool equal(Plane<L, T, Q> const &x, Plane<L, T, Q> const &y, T eps = epsilon<T>()) {
  return all(equal(x.normal, y.normal, eps)) && glm::equal(x.d, y.d, eps);
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool equal(Plane<L, T, Q> const &x, Plane<L, T, Q> const &y, vec<L, T, Q> const &eps) {
  return all(equal(x.normal, y.normal, eps)) && glm::equal(x.d, y.d, eps[0]);
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool equal(Plane<L, T, Q> const &x, Plane<L, T, Q> const &y, int MaxULPs) {
  return all(equal(x.normal, y.normal, MaxULPs)) && glm::equal(x.d, y.d, MaxULPs);
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool equal(Plane<L, T, Q> const &x, Plane<L, T, Q> const &y, vec<L, int, Q> const &MaxULPs) {
  return all(equal(x.normal, y.normal, MaxULPs)) && glm::equal(x.d, y.d, MaxULPs[0]);
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool notEqual(Plane<L, T, Q> const &x, Plane<L, T, Q> const &y, T eps = epsilon<T>()) {
  return any(notEqual(x.normal, y.normal, eps)) || glm::notEqual(x.d, y.d, eps);
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool notEqual(Plane<L, T, Q> const &x, Plane<L, T, Q> const &y, vec<L, T, Q> const &eps) {
  return any(notEqual(x.normal, y.normal, eps)) || glm::notEqual(x.d, y.d, eps[0]);
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool notEqual(Plane<L, T, Q> const &x, Plane<L, T, Q> const &y, int MaxULPs) {
  return any(notEqual(x.normal, y.normal, MaxULPs)) || glm::notEqual(x.d, y.d, MaxULPs);
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool notEqual(Plane<L, T, Q> const &x, Plane<L, T, Q> const &y, vec<L, int, Q> const &MaxULPs) {
  return any(notEqual(x.normal, y.normal, MaxULPs)) || glm::notEqual(x.d, y.d, MaxULPs[0]);
}

/// Construct a Plane by specifying a Ray that exists along the Plane and its normal
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER Plane<L, T, Q> planeFrom(Ray<L, T, Q> const &ray, vec<L, T, Q> const &normal) {
  const vec<L, T, Q> perpNormal = normal - proj(normal, ray.dir);
  return Plane<L, T, Q>(ray.pos, normalize(perpNormal));
}

/// Construct a Plane by specifying a line that exists along the Plane and its normal
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER Plane<L, T, Q> planeFrom(Line<L, T, Q> const &line, vec<L, T, Q> const &normal) {
  const vec<L, T, Q> perpNormal = normal - proj(normal, line.dir);
  return Plane<L, T, Q>(line.pos, normalize(perpNormal));
}

/// Construct a Plane by specifying a segment that exists along the Plane and its normal
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER Plane<L, T, Q> planeFrom(Segment<L, T, Q> const &segment, vec<L, T, Q> const &normal) {
  const vec<L, T, Q> perpNormal = normal - proj(normal, segment.b - segment.a);
  return Plane<L, T, Q>(segment.a, normalize(perpNormal));
}

/// Construct a Plane by specifying a point on the Plane and its normal.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER Plane<3, T, Q> planeFrom(vec<L, T, Q> const &point, vec<L, T, Q> const &normal) {
  return Plane<L, T, Q>(point, normal);
}

/// Construct a Plane by specifying three points on the Plane.
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER Plane<3, T, Q> planeFrom(vec<3, T, Q> const &v1, vec<3, T, Q> const &v2, vec<3, T, Q> const &v3) {
  vec<3, T, Q> normal = cross(v2 - v1, v3 - v1);
  const T len = length(normal);
  if (len > epsilon<T>()) {
    normal /= len;
    return Plane<3, T, Q>(normal, dot(normal, v1));
  }
  return Plane<3, T, Q>(vec<3, T, Q>{ T(0), T(0), T(1) }, T(0));
}

/// Test if any component of the plane is infinite.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool isinf(Plane<L, T, Q> const &plane) {
  return any(isinf(plane.normal)) || glm::isinf(plane.d);
}

/// Test if any component of the plane is isnan.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool isnan(Plane<L, T, Q> const &plane) {
  return any(isnan(plane.normal)) || glm::isnan(plane.d);
}

/// Test if all components of the plane are finite.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool isfinite(Plane<L, T, Q> const &plane) {
  return all(isfinite(plane.normal)) && glm::isfinite(plane.d);
}

// Test whether the plane is degenerate, i.e., not finite or if the normal is null
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool isDegenerate(Plane<L, T, Q> const &plane) {
  return !all(isfinite(plane.normal)) || isNull(plane.normal, epsilon<T>()) || !glm::isfinite(plane.d);
}

/// Return true if two Planes are parallel.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool isParallel(Plane<L, T, Q> const &plane, Plane<L, T, Q> const &other, T eps = epsilon<T>()) {
  return all(epsilonEqual(plane.normal, other.normal, eps));
}

/// Return true if the Plane contains/passes-through the origin (i.e., T(0)).
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool passesThroughOrigin(Plane<L, T, Q> const &plane, T eps = epsilon<T>()) {
  return abs(plane.d) <= eps;
}

/// Return the angle (radians) of intersection between two Planes.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T angle(Plane<L, T, Q> const &plane, Plane<L, T, Q> const &other) {
  return dot(plane.normal, other.normal);
}

/// Reverse the direction of the Plane normal, while still representing the same set of points.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER Plane<L, T, Q> reverseNormal(Plane<L, T, Q> const &plane) {
  return Plane<L, T, Q>(-plane.normal, -plane.d);
}

/// <summary>
/// Return a point on this Plane. The returned point has the property that the
/// line passing through "it" and (0,0,0) is perpendicular to this Plane.
/// </summary>
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> pointOnPlane(Plane<L, T, Q> const &plane) {
  return plane.normal * plane.d;
}

/// Return a point on the Plane at the given parameterized (u, v) coordinates.
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<3, T, Q> point(Plane<3, T, Q> const &plane, T u, T v) {
  vec<3, T, Q> b1, b2;
  perpendicularBasis(plane.normal, b1, b2);
  return pointOnPlane(plane) + u * b1 + v * b2;
}

template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<3, T, Q> point(Plane<3, T, Q> const &plane, T u, T v, vec<3, T, Q> const &referenceOrigin) {
  vec<3, T, Q> b1, b2;
  perpendicularBasis(plane.normal, b1, b2);
  return project(plane, referenceOrigin) + u * b1 + v * b2;
}

/// Refract the given incident vector along the Plane.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> refract(Plane<L, T, Q> const &plane, vec<L, T, Q> const &vec, T eta) {
  return refract(vec, plane.normal, eta);
}

/// Refract the given incident vector along the Plane.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> refract(Plane<L, T, Q> const &plane, vec<L, T, Q> const &vec, T negativeSideRefractionIndex, T positiveSideRefractionIndex) {
  return refract(vec, plane.normal, negativeSideRefractionIndex, positiveSideRefractionIndex);
}

/* Orthographically project the given object onto the plane. */

/// Orthographically project the point onto the Plane
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> project(Plane<L, T, Q> const &plane, vec<L, T, Q> const &point) {
  return (point - (dot(plane.normal, point) - plane.d) * plane.normal);
}

/// Orthographically project the Segment onto the Plane
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER Segment<L, T, Q> project(Plane<L, T, Q> const &plane, Segment<L, T, Q> const &segment) {
  return Segment<L, T, Q>(project(plane, segment.a), project(plane, segment.b));
}

/// Orthographically project the Line onto the Plane
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER Line<L, T, Q> project(Plane<L, T, Q> const &plane, Line<L, T, Q> const &line, bool *nonDegenerate = GLM_NULLPTR) {
  Line<L, T, Q> l;
  l.pos = project(plane, line.pos);
  l.dir = normalize(line.dir - proj(line.dir, plane.normal));
  if (nonDegenerate)
    *nonDegenerate = (length(l.dir) > T(0));
  return l;
}

/// Orthographically project the Ray onto the Plane
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER Line<L, T, Q> project(Plane<L, T, Q> const &plane, Ray<L, T, Q> const &ray, bool *nonDegenerate = GLM_NULLPTR) {
  Line<L, T, Q> l;
  l.pos = project(plane, ray.pos);
  l.dir = normalize(ray.dir - proj(ray.dir, plane.normal));
  if (nonDegenerate)
    *nonDegenerate = (length(l.dir) > T(0));
  return l;
}

/// Orthographically project the Triangle onto the Plane
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER Triangle<L, T, Q> project(Plane<L, T, Q> const &plane, Triangle<L, T, Q> const &triangle) {
  return Triangle<L, T, Q>(project(plane, triangle.a), project(plane, triangle.b), project(plane, triangle.c));
}

/// Projects the given point to the negative halfspace of the Plane
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> projectToNegativeHalf(Plane<L, T, Q> const &plane, vec<L, T, Q> const &point) {
  return point - max(T(0), (dot(plane.normal, point) - plane.d)) * plane.normal;
}

/// Projects the given point to the positive halfspace of the Plane
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> projectToPositiveHalf(Plane<L, T, Q> const &plane, vec<L, T, Q> const &point) {
  return point - min(T(0), (dot(plane.normal, point) - plane.d)) * plane.normal;
}

/* Return the signed distance between the plane and the given object. */

/// Return the signed distance of this Plane to the given point
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T signedDistance(Plane<L, T, Q> const &plane, vec<L, T, Q> const &point) {
  return dot(plane.normal, point) - plane.d;
}

/// Return the signed distance of this Plane to the given Object
template<length_t L, typename T, qualifier Q, typename Object>
GLM_GEOM_QUALIFIER T signedDistance(Plane<L, T, Q> const &plane, Object const &object) {
  T pMin(0), pMax(0);
  projectToAxis(object, plane.normal, pMin, pMax);
  pMin -= plane.d;
  pMax -= plane.d;
  if (pMin * pMax <= T(0))
    return T(0);
  return abs(pMin) < abs(pMax) ? pMin : pMax;
}

/// Return true if two points are on the same side of this Plane.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool areOnSameSide(Plane<L, T, Q> const &plane, vec<L, T, Q> const &p1, vec<L, T, Q> const &p2) {
  return signedDistance(plane, p1) * signedDistance(plane, p2) >= T(0);
}

/// Test if the given direction vector points towards the positive side of this Plane.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool isInPositiveDirection(Plane<L, T, Q> const &plane, vec<L, T, Q> const &direction) {
  return dot(plane.normal, direction) >= T(0);
}

/// Test if the given point exists on the positive side of this Plane.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool isOnPositiveSide(Plane<L, T, Q> const &plane, vec<L, T, Q> const &point) {
  return signedDistance(plane, point) >= T(0);
}

/// <summary>
/// Triangle/Plane intersection test. Returning:
///   1 - If the Triangle is completely in the positive halfspace of the Plane;
///  -1 - If the Triangle is completely in the negative halfspace of the Plane;
///   0 - If the Triangle intersects the Plane.
/// </summary>
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER int examineSide(Plane<L, T, Q> const &plane, Triangle<L, T, Q> const &triangle, T eps = epsilon<T>()) {
  const T a = signedDistance(plane, triangle.a);
  const T b = signedDistance(plane, triangle.b);
  const T c = signedDistance(plane, triangle.c);
  if (a >= -eps && b >= -eps && c >= -eps)
    return 1;
  if (a <= eps && b <= eps && c <= eps)
    return -1;
  return 0;
}

/* Return the distance between the plane and the given object(s). */

/// Return the distance between the Plane and given point
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Plane<L, T, Q> const &plane, vec<L, T, Q> const &point) {
  return abs(signedDistance(plane, point));
}

/// Return the distance between the Plane and given Segment
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Plane<L, T, Q> const &plane, Segment<L, T, Q> const &segment) {
  return distance(segment, plane);
}

/// Return the distance between the Plane and given Sphere
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Plane<L, T, Q> const &plane, Sphere<L, T, Q> const &sphere) {
  return max(T(0), distance(plane, sphere.pos) - sphere.r);
}

/// Return the distance between the Plane and given Capsule
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Plane<L, T, Q> const &plane, Capsule<L, T, Q> const &capsule) {
  return max(T(0), distance(plane, capsule.l) - capsule.r);
}

/// Return the signed distance between the Plane and given AABB
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T signedDistance(Plane<L, T, Q> const &plane, AABB<L, T, Q> const &aabb) {
  return signedDistance<L, T, Q, AABB<L, T, Q>>(plane, aabb);
}

/// Return the signed distance between the Plane and given OBB
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER T signedDistance(Plane<3, T, Q> const &plane, OBB<T, Q> const &obb) {
  return signedDistance<3, T, Q, OBB<T, Q>>(plane, obb);
}

/// Return the signed distance between the Plane and given Line
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T signedDistance(Plane<L, T, Q> const &plane, Line<L, T, Q> const &line) {
  return signedDistance<L, T, Q, Line<L, T, Q>>(plane, line);
}

/// Return the signed distance between the Plane and given Segment
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T signedDistance(Plane<L, T, Q> const &plane, Segment<L, T, Q> const &segment) {
  return signedDistance<L, T, Q, Segment<L, T, Q>>(plane, segment);
}

/// Return the signed distance between the Plane and given Ray
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T signedDistance(Plane<L, T, Q> const &plane, Ray<L, T, Q> const &ray) {
  return signedDistance<L, T, Q, Ray<L, T, Q>>(plane, ray);
}

/// Return the signed distance between the Plane and given Sphere
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T signedDistance(Plane<L, T, Q> const &plane, Sphere<L, T, Q> const &sphere) {
  return signedDistance<L, T, Q, Sphere<L, T, Q>>(plane, sphere);
}

/// Return the signed distance between the Plane and given Capsule
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T signedDistance(Plane<L, T, Q> const &plane, Capsule<L, T, Q> const &capsule) {
  return signedDistance<L, T, Q, Capsule<L, T, Q>>(plane, capsule);
}

#if 0
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T signedDistance(Plane<L, T, Q> const &plane, Disc<L, T, Q> const &disc) { }
#endif

/// Return the signed distance between the Plane and given Triangle
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T signedDistance(Plane<L, T, Q> const &plane, Triangle<L, T, Q> const &triangle) {
  return signedDistance<L, T, Q, Triangle<L, T, Q>>(plane, triangle);
}

/// Return an affine transformation matrix that projects orthographically onto the Plane
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER mat<4, 3, T, Q> orthoProjection(Plane<3, T, Q> const &plane) {
  return orthoProjection<4, 3, T, Q>(plane.normal, plane.d);
}

/// Mirrors the given point with respect to the Plane.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> mirror(Plane<L, T, Q> const &plane, vec<L, T, Q> const &point) {
  return point - T(2) * (dot(point, plane.normal) - plane.d) * plane.normal;
}

/// Return a transformation matrix that mirrors objects along the Plane.
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER mat<4, 3, T, Q> mirrorMatrix(Plane<3, T, Q> const &plane) {
  return planeMirror<4, 3, T, Q>(plane.normal, plane.d);
}

/* Return the closest point on this plane to the given object */

/// Return the closest point between the Plane and given point.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Plane<L, T, Q> const &plane, vec<L, T, Q> const &point) {
  return project(plane, point);
}

/// Return the closest point between the Plane and given Ray
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Plane<L, T, Q> const &plane, Ray<L, T, Q> const &ray) {
  const T denom = dot(plane.normal, ray.dir);
  if (detail::approx_zero(denom))
    return project(plane, ray.pos);

  const T t = (plane.d - dot(plane.normal, ray.pos)) / denom;
  if (t >= T(0))
    return getPoint(ray, t);
  else
    return project(plane, ray.pos);
}

/// Return the closest point between the Plane and given Segment
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Plane<L, T, Q> const &plane, Segment<L, T, Q> const &segment) {
  const T aDist = dot(plane.normal, segment.a);
  const T bDist = dot(plane.normal, segment.b);
  const T denom = bDist - aDist;
  if (detail::approx_zero(denom))
    return project(plane, abs(aDist) < abs(bDist) ? segment.a : segment.b);

  const T t = glm::clamp((plane.d - dot(plane.normal, segment.a)) / denom, T(0), T(1));
  return project(plane, getPoint(segment, t));
}

/* Test if this plane contains the given object(s) */

/// Plane/Point containment test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool contains(Plane<L, T, Q> const &plane, vec<L, T, Q> const &point, T distanceThreshold = contains_eps<T>()) {
  return distance(plane, point) <= distanceThreshold;
}

/// Plane/Line containment test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool contains(Plane<L, T, Q> const &plane, Line<L, T, Q> const &line, T distanceThreshold = contains_eps<T>()) {
  return contains(plane, line.pos, distanceThreshold) && areOrthonormal2(line.dir, plane.normal, distanceThreshold);
}

/// Plane/Ray containment test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool contains(Plane<L, T, Q> const &plane, Ray<L, T, Q> const &ray, T distanceThreshold = contains_eps<T>()) {
  return contains(plane, ray.pos, distanceThreshold) && areOrthonormal2(ray.dir, plane.normal, distanceThreshold);
}

/// Plane/Segment containment test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool contains(Plane<L, T, Q> const &plane, Segment<L, T, Q> const &segment, T distanceThreshold = contains_eps<T>()) {
  return contains(plane, segment.a, distanceThreshold) && contains(plane, segment.b, distanceThreshold);
}

/// Plane/Triangle containment test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool contains(Plane<L, T, Q> const &plane, Triangle<L, T, Q> const &triangle, T distanceThreshold = contains_eps<T>()) {
  return contains(plane, triangle.a, distanceThreshold)
         && contains(plane, triangle.b, distanceThreshold)
         && contains(plane, triangle.c, distanceThreshold);
}

/// Plane/Disc containment test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool contains(Plane<L, T, Q> const &plane, Disc<L, T, Q> const &disc, T eps = epsilon<T>()) {
  return contains(disc, plane, eps);
}

/* Test whether the plane and the given object intersect */

/// This function attempts to improve stability with lines that are almost parallel with the Plane.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersectLinePlane(vec<L, T, Q> const &planeNormal, T planeD, vec<L, T, Q> const &linePos, vec<L, T, Q> const &lineDir, T &d, T eps = intersect_eps<T>()) {
  const T denom = dot(planeNormal, lineDir);
  if (abs(denom) >= eps) {  // line starting point to point of intersection
    d = (planeD - dot(planeNormal, linePos)) / denom;
    return true;
  }

  d = T(0);
  return epsilonEqual(dot(planeNormal, linePos), planeD, eps);
}

/// Plane/Ray intersection test, includes parametric ray coordinate.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Plane<L, T, Q> const &plane, Ray<L, T, Q> const &ray, T &d) {
  if (intersectLinePlane(plane.normal, plane.d, ray.pos, ray.dir, d))
    return d >= T(0);
  return false;
}

/// Plane/Line intersection test, includes parametric line coordinate.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Plane<L, T, Q> const &plane, Line<L, T, Q> const &line, T &d) {
  return intersectLinePlane(plane.normal, plane.d, line.pos, line.dir, d);
}

/// Plane/Segment intersection test, includes parametric segment coordinate.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Plane<L, T, Q> const &plane, Segment<L, T, Q> const &segment, T &d) {
  if (intersectLinePlane(plane.normal, plane.d, segment.a, segment.dir(), d)) {
    d /= length(segment);
    return d >= T(0) && d <= T(1);
  }
  return false;
}

/// Plane/Sphere intersection test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Plane<L, T, Q> const &plane, Sphere<L, T, Q> const &sphere, T eps = epsilon<T>()) {
  return distance(plane, sphere.pos) <= sphere.r + eps;
}

/// Plane/AABB intersection test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Plane<L, T, Q> const &plane, AABB<L, T, Q> const &aabb) {
  const vec<L, T, Q> c = centroid(aabb);
  const vec<L, T, Q> e = halfSize(aabb);

  T r(0);  // Compute projection interval radius; aabb.center + t * plane.normal
  for (length_t i = 0; i < L; ++i)
    r += e[i] * abs(plane.normal[i]);
  return abs(dot(plane.normal, c) - plane.d) <= r;
}

/// Plane/OBB intersection test
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Plane<3, T, Q> const &plane, OBB<T, Q> const &obb) {
  return intersects(obb, plane);
}

/// Plane/Ray intersection test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Plane<L, T, Q> const &plane, Ray<L, T, Q> const &ray) {
  T d(0);
  return intersects(plane, ray, d);
}

/// Plane/Line intersection test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Plane<L, T, Q> const &plane, Line<L, T, Q> const &line) {
  T d(0);
  return intersects(plane, line, d);
}

/// Plane/Segment intersection test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Plane<L, T, Q> const &plane, Segment<L, T, Q> const &segment) {
  T d(0);
  return intersects(plane, segment, d);
}

/// Plane/Plane intersection test; includes line of intersection
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Plane<L, T, Q> const &plane, Plane<L, T, Q> const &other, Line<L, T, Q> &line) {
  const vec<L, T, Q> direction = cross(plane.normal, other.normal);
  const T denominator = length2(direction);
  if (detail::exactly_zero(denominator)) {
    line = Line<L, T, Q>(T(0));
    return false;  // If direction is zero, the planes are parallel/coincident.
  }
  else {
    const vec<L, T, Q> temp = plane.d * other.normal - other.d * plane.normal;
    line = Line<L, T, Q>(cross(temp, direction), normalize(direction));
    return true;
  }
}

/// Plane/Plane/Plane intersection tests; includes point of intersection
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Plane<L, T, Q> const &a, Plane<L, T, Q> const &b, Plane<L, T, Q> const &c, vec<L, T, Q> &result) {
  const T denom = dot(cross(a.normal, b.normal), c.normal);
  if (denom >= intersect_eps<T>()) {
    result = ((cross(b.normal, c.normal) * a.d) + (cross(c.normal, a.normal) * b.d) + (cross(a.normal, b.normal) * c.d)) / denom;
    return true;
  }
  return false;
}

/// Plane/Capsule intersection test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Plane<L, T, Q> const &plane, Capsule<L, T, Q> const &capsule) {
  return intersects(capsule, plane);
}

/// Plane/Disc intersection test, includes points of intersection
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER int intersects(Plane<L, T, Q> const &plane, Disc<L, T, Q> const &disc, vec<L, T, Q> &pt1, vec<L, T, Q> &pt2, T eps = intersect_eps<T>()) {
  return intersects(disc, plane, pt1, pt2, eps);
}

/// Plane/Triangle intersection test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Plane<L, T, Q> const &plane, Triangle<L, T, Q> const &triangle) {
  const T a = signedDistance(plane, triangle.a);
  const T b = signedDistance(plane, triangle.b);
  const T c = signedDistance(plane, triangle.c);
  return (a * b <= T(0) || a * c <= T(0));
}

/* Plane/Object clipping */

/// @private In-place clipping
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool clip(Plane<L, T, Q> const &plane, vec<L, T, Q> &a, vec<L, T, Q> &b) {
  T t(0);
  bool intersects = intersectLinePlane(plane.normal, plane.d, a, b - a, t);
  if (!intersects || t <= T(0) || t >= T(1))
    return signedDistance(plane, a) > T(0);  // Within the positive/negative halfspace

  const vec<L, T, Q> pt = a + (b - a) * t;  // Point of intersection
  if (isOnPositiveSide(plane, a))
    b = pt;
  else
    a = pt;
  return true;
}

/// <summary>
/// Clips a Segment against the Plane, i.e., remove part of the line that exists
/// in the negative halfspace of the Plane.
/// </summary>
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER Segment<L, T, Q> clip(Plane<L, T, Q> const &plane, Segment<L, T, Q> const &segment) {
  Segment<L, T, Q> result(segment);
  return clip(plane, result.a, result.b) ? result : segment;
}

/// <summary>
/// Clips a line against the Plane, i.e., remove part of the line that exists in
/// the negative halfspace of the Plane.
/// </summary>
/// <returns>
///   0 - If clipping removes the entire Line (the Line exists entirely in the negative halfspace).
///   1 - If clipping results in a Ray (clipped at the point of intersection).
///   2 - If clipping keeps the entire Line (the Line exists entirely in the positive halfspace).
/// </returns>
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER int clip(Plane<L, T, Q> const &plane, Line<L, T, Q> const &line, Ray<L, T, Q> &outRay) {
  T t(0);
  if (!intersectLinePlane(plane.normal, plane.d, line.pos, line.dir, t)) {
    outRay.pos = line.pos;
    outRay.dir = line.dir;
    return signedDistance(plane, line.pos) <= T(0) ? 0 : 2;  // Completely within the positive/negative halfspace
  }

  outRay.pos = line.pos + line.dir * t;
  if (dot(line.dir, plane.normal) >= T(0))
    outRay.dir = line.dir;
  else
    outRay.dir = -line.dir;
  return 1;
}

#if 0
/// <summary>
/// Clip a polygon against a Plane, i.e., remove part(s) of the polygon that lie
/// in the negative halfspace of the Plane and return a new polygon.
/// </summary>
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER Polygon<L, T, Q> clip(Plane<L, T, Q> const &plane, Polygon<L, T, Q> const &polygon) {
}
#endif

/// <summary>
/// Clip a Triangle against the Plane, i.e., create one or more Triangles that
/// reside strictly in the positive/negative halfspaces of the Plane.
/// </summary>
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER int clip(Plane<L, T, Q> const &plane, Triangle<L, T, Q> const &triangle, Triangle<L, T, Q> &t1, Triangle<L, T, Q> &t2) {
  const bool aSide = isOnPositiveSide(plane, triangle.a);
  const bool bSide = isOnPositiveSide(plane, triangle.b);
  const bool cSide = isOnPositiveSide(plane, triangle.c);
  switch ((aSide ? 1 : 0) + (bSide ? 1 : 0) + (cSide ? 1 : 0)) {
    case 1: {
      if (bSide)
        t1 = Triangle<L, T, Q>(triangle.b, triangle.c, triangle.a);
      else if (cSide)
        t1 = Triangle<L, T, Q>(triangle.c, triangle.b, triangle.a);
      else
        t1 = triangle;

      T t, r;
      intersects(plane, Segment<L, T, Q>(t1.a, t1.b), t);
      intersects(plane, Segment<L, T, Q>(t1.a, t1.c), r);
      t1.b = t1.a + (t1.b - t1.a) * t;
      t1.c = t1.a + (t1.c - t1.a) * r;
      return 1;
    }
    case 2: {
      if (!bSide)
        t1 = Triangle<L, T, Q>(triangle.b, triangle.c, triangle.a);
      else if (!cSide)
        t1 = Triangle<L, T, Q>(triangle.c, triangle.b, triangle.a);
      else
        t1 = triangle;

      T t, r;
      intersects(plane, Segment<L, T, Q>(t1.a, t1.b), t);
      intersects(plane, Segment<L, T, Q>(t1.a, t1.c), r);
      t2 = Triangle<L, T, Q>(t1.c, t1.a + (t1.c - t1.a) * r, t1.a + (t1.b - t1.a) * t);
      t1.a = t2.c;
      return 2;
    }
    case 3: {  // All vertices are on the positive side.
      t1 = triangle;
      return 1;
    }
    default: {
      break;
    }
  }
  return 0;
}

/// @}
GLM_GEOM_END_NAMESPACE
#if GLM_GEOM_TOSTRING
namespace glm {
  namespace detail {
    template<glm::length_t L, typename T, qualifier Q>
    struct compute_to_string<Plane<L, T, Q>> {
      GLM_GEOM_QUALIFIER std::string call(Plane<L, T, Q> const &plane) {
        char const *LiteralStr = literal<T, std::numeric_limits<T>::is_iec559>::value();
        std::string FormatStr(detail::format("plane(%%s, %s)", LiteralStr));
        return detail::format(FormatStr.c_str(),
          glm::to_string(plane.normal).c_str(),
          static_cast<typename cast<T>::value_type>(plane.d)
        );
      }
    };
  }
}
#endif
