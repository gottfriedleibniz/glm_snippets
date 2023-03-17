/// @ref geom_obb
/// @file geom/obb.hpp
///
/// @defgroup geom_obb OBB
/// @ingroup geom
///

#pragma once

#include "setup.hpp"
#include "aabb.hpp"
#include "segment.hpp"
#include "sphere.hpp"
#include "polygon.hpp"
#include "vector_extensions.hpp"

#include <glm/gtc/matrix_inverse.hpp>
#include <glm/gtx/component_wise.hpp>
#include <glm/gtx/quaternion.hpp>
#include <glm/gtx/scalar_relational.hpp>
#include <glm/ext/quaternion_common.hpp>

#include "../quat_extensions.hpp"
#if GLM_MESSAGES == GLM_ENABLE && !defined(GLM_EXT_INCLUDED)
  #pragma message("GLM: GLM_EXT_GEOM_obb extension included")
#endif

GLM_GEOM_BEGIN_NAMESPACE
/// @addtogroup geom_obb
/// @{

/// <summary>
/// An oriented bounding box.
/// </summary>
template<typename T, qualifier Q>
struct OBB {

  // -- Implementation detail --

  typedef T value_type;
  typedef OBB<T, Q> type;
  typedef vec<3, T, Q> point_type;
  typedef qua<T, Q> rotation_type;

  // -- Data --

  point_type position;
  rotation_type rotation;
  point_type half;  // half lengths

#if GLM_CONFIG_DEFAULTED_DEFAULT_CTOR == GLM_ENABLE
  OBB() GLM_DEFAULT_CTOR;
#else
  OBB()
  #if GLM_CONFIG_CTOR_INIT != GLM_CTOR_INIT_DISABLE
    : position(T(0)), rotation(), half(T(0))
  #endif
  {
  }
#endif

  OBB(point_type const &p, rotation_type const &r, point_type const &h)
    : position(p), rotation(r), half(h) {
  }

  OBB(AABB<3, T, Q> const &aabb)
    : position(centroid(aabb)),
      rotation(quat_identity<T, Q>()),
      half(halfSize(aabb)) {
  }

  OBB(Sphere<3, T, Q> const &sphere)
    : position(sphere.pos),
      rotation(quat_identity<T, Q>()),
      half(sphere.r) {
  }

  vec<3, T, Q> Localize(vec<3, T, Q> const &point) const {
    return conjugate(rotation) * (point - position);
  }

  mat<3, 3, T, Q> Axes() const {
    return mat3_cast(rotation);
  }

  vec<3, T, Q> XAxis() const {
    return rotation * vec<3, T, Q>(T(1), T(0), T(0));
  }
  vec<3, T, Q> YAxis() const {
    return rotation * vec<3, T, Q>(T(0), T(1), T(0));
  }
  vec<3, T, Q> ZAxis() const {
    return rotation * vec<3, T, Q>(T(0), T(0), T(1));
  }

  vec<3, T, Q> Up() const {
    return rotation * glm::up<T, Q>();
  }
  vec<3, T, Q> Right() const {
    return rotation * glm::right<T, Q>();
  }
  vec<3, T, Q> Forward() const {
    return rotation * glm::forward<T, Q>();
  }
};

template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER OBB<T, Q> obbFromAABB(AABB<3, T, Q> const &aabb) {
  return OBB<T, Q>(aabb);
}

template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER OBB<T, Q> obbFromSphere(Sphere<3, T, Q> const &sphere) {
  return OBB<T, Q>(sphere);
}

/// Apply a generic matrix transformation to the OBB.
template<typename T, qualifier Q, class Matrix>
GLM_GEOM_QUALIFIER_OUTLINE OBB<T, Q> transformAsOBB(OBB<T, Q> const &obb, Matrix const &m) {
  OBB<T, Q> Result;
  Result.position = transformPos(m, obb.position);
  Result.rotation = quat_cast(m) * obb.rotation;
  Result.half = extractScale(m) * obb.half;
  return Result;
}

/// Apply a generic matrix transformation to the AABB, returning an OBB.
template<typename T, qualifier Q, class Matrix>
GLM_GEOM_QUALIFIER_OUTLINE OBB<T, Q> transformAsOBB(AABB<3, T, Q> const &aabb, Matrix const &m) {
  OBB<T, Q> Result;
  Result.position = centroid(aabb);
  Result.rotation = quat_cast(m);
  Result.half = extractScale(m) * halfSize(aabb);
  return Result;
}

template<typename T, qualifier Q>
static OBB<T, Q> operator-(OBB<T, Q> const &obb) {
  return OBB<T, Q>(obb.position, glm::inverse(obb.rotation), obb.half);
}

template<typename T, qualifier Q>
static bool operator==(OBB<T, Q> const &o1, OBB<T, Q> const &o2) {
  return o1.position == o2.position && o1.rotation == o2.rotation && o1.half == o2.half;
}

template<typename T, qualifier Q>
static OBB<T, Q> operator+(OBB<T, Q> const &obb, vec<3, T, Q> const &point) {
  return OBB<T, Q>(obb.position + point, obb.rotation, obb.half);
}

template<typename T, qualifier Q>
static OBB<T, Q> operator-(OBB<T, Q> const &obb, vec<3, T, Q> const &point) {
  return OBB<T, Q>(obb.position - point, obb.rotation, obb.half);
}

template<typename T, qualifier Q>
static OBB<T, Q> operator*(mat<3, 3, T, Q> const &m, OBB<T, Q> const &obb) {
  return transformAsOBB<T, Q, mat<3, 3, T, Q>>(obb, m);
}

template<typename T, qualifier Q>
static OBB<T, Q> operator*(mat<3, 4, T, Q> const &m, OBB<T, Q> const &obb) {
  return transformAsOBB<T, Q, mat<3, 4, T, Q>>(obb, m);
}

template<typename T, qualifier Q>
static OBB<T, Q> operator*(mat<4, 3, T, Q> const &m, OBB<T, Q> const &obb) {
  return transformAsOBB<T, Q, mat<4, 3, T, Q>>(obb, m);
}

template<typename T, qualifier Q>
static OBB<T, Q> operator*(mat<4, 4, T, Q> const &m, OBB<T, Q> const &obb) {
  return transformAsOBB<T, Q, mat<4, 4, T, Q>>(obb, m);
}

template<typename T, qualifier Q>
static OBB<T, Q> operator*(qua<T, Q> const &q, OBB<T, Q> const &obb) {
  return transformAsOBB<T, Q, mat<3, 3, T, Q>>(obb, mat3_cast(q));
}

template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool equal(OBB<T, Q> const &x, OBB<T, Q> const &y, T eps = epsilon<T>()) {
  return all(equal(x.position, y.position, eps))
         && all(equal(x.rotation, y.rotation, eps))
         && all(equal(x.half, y.half, eps));
}

template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool equal(OBB<T, Q> const &x, OBB<T, Q> const &y, vec<3, T, Q> const &eps) {
  return all(equal(x.position, y.position, eps))
         && all(equal(x.rotation, y.rotation, vec<4, T, Q>(eps.x)))
         && all(equal(x.half, y.half, eps));
}

template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool equal(OBB<T, Q> const &x, OBB<T, Q> const &y, int MaxULPs) {
  return all(equal(x.position, y.position, MaxULPs))
         && all(equal(x.rotation, y.rotation, MaxULPs))
         && all(equal(x.half, y.half, MaxULPs));
}

template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool equal(OBB<T, Q> const &x, OBB<T, Q> const &y, vec<3, int, Q> const &MaxULPs) {
  return all(equal(x.position, y.position, MaxULPs))
         && all(equal(x.rotation, y.rotation, vec<4, int, Q>(MaxULPs.x)))
         && all(equal(x.half, y.half, MaxULPs));
}

template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool notEqual(OBB<T, Q> const &x, OBB<T, Q> const &y, T eps = epsilon<T>()) {
  return any(notEqual(x.position, y.position, eps))
         || any(notEqual(x.rotation, y.rotation, eps))
         || any(notEqual(x.half, y.half, eps));
}

template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool notEqual(OBB<T, Q> const &x, OBB<T, Q> const &y, vec<3, T, Q> const &eps) {
  return any(notEqual(x.position, y.position, eps))
         || any(notEqual(x.rotation, y.rotation, vec<4, T, Q>(eps, eps.x)))
         || any(notEqual(x.half, y.half, eps));
}

template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool notEqual(OBB<T, Q> const &x, OBB<T, Q> const &y, int MaxULPs) {
  return any(notEqual(x.position, y.position, MaxULPs))
         || any(notEqual(x.rotation, y.rotation, MaxULPs))
         || any(notEqual(x.half, y.half, MaxULPs));
}

template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool notEqual(OBB<T, Q> const &x, OBB<T, Q> const &y, vec<3, int, Q> const &MaxULPs) {
  return any(notEqual(x.position, y.position, MaxULPs))
         || any(notEqual(x.rotation, y.rotation, vec<4, T, Q>(MaxULPs, MaxULPs.x)))
         || any(notEqual(x.half, y.half, MaxULPs));
}

/// Test if any component of the OBB is infinite.
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool isinf(OBB<T, Q> const &obb) {
  return any(isinf(obb.position))
         || any(isinf(obb.rotation))
         || any(isinf(obb.half));
}

/// Test if any component of the OBB is NaN
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool isnan(OBB<T, Q> const &obb) {
  return any(isnan(obb.position))
         || any(isnan(obb.rotation))
         || any(isnan(obb.half));
}

/// Test if all components of the OBB are finite.
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool isfinite(OBB<T, Q> const &obb) {
  return all(isfinite(obb.position))
         && all(isfinite(obb.rotation))
         && all(isfinite(obb.half));
}

/// Return true if the OBB is degenerate, i.e., does not span in a strictly positive volume.
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool isDegenerate(OBB<T, Q> const &obb) {
  return obb.half.x <= T(0) || obb.half.y <= T(0) || obb.half.z <= T(0);
}

/// Return obb.position
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool centroid(OBB<T, Q> const &obb) {
  return obb.position;
}

/// Return the length of the OBB along each dimension: half * 2
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<3, T, Q> size(OBB<T, Q> const &obb) {
  return T(2) * obb.half;
}

/// <summary>
/// Return a diagonal vector of the OBB corresponding to <0,0,0> to <x,y,z> in
/// the world space of the OBB
/// </summary>
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<3, T, Q> halfDiagonal(OBB<T, Q> const &obb) {
  const mat<3, 3, T, Q> axis = obb.Axes();
  return axis[0] * obb.half.x + axis[1] * obb.half.y + axis[2] * obb.half.z;
}

/// <summary>
/// Return a diagonal vector of the OBB corresponding to <-x,-y,-z> to <x,y,z>
/// in the world space of the OBB
/// </summary>
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<3, T, Q> diagonal(OBB<T, Q> const &obb) {
  return T(2) * halfDiagonal(obb);
}

/// Return the volume of the OBB
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER T volume(OBB<T, Q> const &obb) {
  const vec<3, T, Q> s = size(obb);
  return s.x * s.y * s.z;
}

/// Return the surface area of the faces of the OBB
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER T surfaceArea(OBB<T, Q> const &obb) {
  const vec<3, T, Q> s = size(obb);
  return T(2) * (s.x * s.y + s.x * s.z + s.y * s.z);
}

/// Return the smallest Sphere that contains the OBB.
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER Sphere<3, T, Q> minimalEnclosingSphere(OBB<T, Q> const &obb) {
  return Sphere<3, T, Q>(obb.position, length(halfDiagonal(obb)));
}

/// Return the largest Sphere that can fit inside the OBB.
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER Sphere<3, T, Q> maximalContainedSphere(OBB<T, Q> const &obb) {
  return Sphere<3, T, Q>(obb.position, compMin(obb.half));
}

/// Return the transformation matrix that maps from OBB local space to world space.
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER mat<4, 4, T, Q> localToWorld(OBB<T, Q> const &obb) {
  const mat<3, 3, T, Q> axis = obb.Axes();
  const vec<3, T, Q> halfdiag = axis[0] * obb.half.x - axis[1] * obb.half.y - axis[2] * obb.half.z;
  return mat<4, 4, T, Q>(
    vec<4, T, Q>(axis[0], T(0)),
    vec<4, T, Q>(axis[1], T(0)),
    vec<4, T, Q>(axis[2], T(0)),
    vec<4, T, Q>(obb.position - halfdiag, T(1))
  );
}

/// Return the transformation matrix that maps from world space to OBB local space.
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER mat<4, 4, T, Q> worldToLocal(OBB<T, Q> const &obb) {
  return affineInverse(localToWorld(obb));
}

/// Return a corner point of the OBB: [0, 7].
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER_OUTLINE vec<3, T, Q> cornerPoint(OBB<T, Q> const &obb, length_t index) {
  const mat<3, 3, T, Q> axis = obb.Axes();
  const vec<3, T, Q> offset(
    (index & 4) ? obb.half.x : -obb.half.x,
    (index & 2) ? obb.half.y : -obb.half.y,
    (index & 1) ? obb.half.z : -obb.half.z
  );
  return obb.position + offset.x * axis[0] + offset.y * axis[1] + offset.z * axis[2];
}

/// Return an edge of the OBB: [0, 11].
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER_OUTLINE Segment<3, T, Q> edge(OBB<T, Q> const &obb, length_t edgeIndex) {
  switch (edgeIndex) {
    case 1: return Segment<3, T, Q>(cornerPoint(obb, 0), cornerPoint(obb, 2));
    case 2: return Segment<3, T, Q>(cornerPoint(obb, 0), cornerPoint(obb, 4));
    case 3: return Segment<3, T, Q>(cornerPoint(obb, 1), cornerPoint(obb, 3));
    case 4: return Segment<3, T, Q>(cornerPoint(obb, 1), cornerPoint(obb, 5));
    case 5: return Segment<3, T, Q>(cornerPoint(obb, 2), cornerPoint(obb, 3));
    case 6: return Segment<3, T, Q>(cornerPoint(obb, 2), cornerPoint(obb, 6));
    case 7: return Segment<3, T, Q>(cornerPoint(obb, 3), cornerPoint(obb, 7));
    case 8: return Segment<3, T, Q>(cornerPoint(obb, 4), cornerPoint(obb, 5));
    case 9: return Segment<3, T, Q>(cornerPoint(obb, 4), cornerPoint(obb, 6));
    case 10: return Segment<3, T, Q>(cornerPoint(obb, 5), cornerPoint(obb, 7));
    case 11: return Segment<3, T, Q>(cornerPoint(obb, 6), cornerPoint(obb, 7));
    case 0:
    default: {
      return Segment<3, T, Q>(cornerPoint(obb, 0), cornerPoint(obb, 1));
    }
  }
}

/// Return an extreme point along the OBB, i.e., the farthest point in a given direction.
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<3, T, Q> extremePoint(OBB<T, Q> const &obb, vec<3, T, Q> const &direction) {
  const mat<3, 3, T, Q> axis = obb.Axes();
  vec<3, T, Q> pt = obb.position;
  pt += axis[0] * (dot(direction, axis[0]) >= T(0) ? obb.half.x : -obb.half.x);
  pt += axis[1] * (dot(direction, axis[1]) >= T(0) ? obb.half.y : -obb.half.y);
  pt += axis[2] * (dot(direction, axis[2]) >= T(0) ? obb.half.z : -obb.half.z);
  return pt;
}

template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<3, T, Q> extremePoint(OBB<T, Q> const &obb, vec<3, T, Q> const &direction, T &distance) {
  const vec<3, T, Q> ep = extremePoint(obb, direction);
  distance = dot(ep, direction);
  return ep;
}

/// Project the OBB onto the provided axis.
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER void projectToAxis(OBB<T, Q> const &obb, vec<3, T, Q> const &axis, T &dMin, T &dMax) {
  const mat<3, 3, T, Q> obbAxis = obb.Axes();
  const T x = abs(dot(axis, obbAxis[0]) * obb.half.x);
  const T y = abs(dot(axis, obbAxis[1]) * obb.half.y);
  const T z = abs(dot(axis, obbAxis[2]) * obb.half.z);
  const T pos = dot(axis, obb.position);
  dMin = pos - x - y - z;
  dMax = pos + x + y + z;
}

/// Return a point along an edge of the OBB.
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER_OUTLINE vec<3, T, Q> pointOnEdge(OBB<T, Q> const &obb, length_t edgeIndex, T u) {
  const length_t eid = glm::clamp(edgeIndex, length_t(0), length_t(11));
  const mat<3, 3, T, Q> axis = obb.Axes();
  const vec<3, T, Q> d = axis[eid / 4] * (T(2) * u - T(1)) * obb.half[eid / 4];
  switch (eid) {
    case 1: return obb.position - obb.half.y * axis[1] + obb.half.z * axis[2] + d;
    case 2: return obb.position + obb.half.y * axis[1] - obb.half.z * axis[2] + d;
    case 3: return obb.position + obb.half.y * axis[1] + obb.half.z * axis[2] + d;
    case 4: return obb.position - obb.half.x * axis[0] - obb.half.z * axis[2] + d;
    case 5: return obb.position - obb.half.x * axis[0] + obb.half.z * axis[2] + d;
    case 6: return obb.position + obb.half.x * axis[0] - obb.half.z * axis[2] + d;
    case 7: return obb.position + obb.half.x * axis[0] + obb.half.z * axis[2] + d;
    case 8: return obb.position - obb.half.x * axis[0] - obb.half.y * axis[1] + d;
    case 9: return obb.position - obb.half.x * axis[0] + obb.half.y * axis[1] + d;
    case 10: return obb.position + obb.half.x * axis[0] - obb.half.y * axis[1] + d;
    case 11: return obb.position + obb.half.x * axis[0] + obb.half.y * axis[1] + d;
    case 0:
    default: {
      return obb.position - obb.half.y * axis[1] - obb.half.z * axis[2] + d;
    }
  }
}

/// Return the point at the center of the given face, [0, 5], of the OBB.
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER_OUTLINE vec<3, T, Q> faceCentroid(OBB<T, Q> const &obb, length_t faceIndex) {
  switch (faceIndex) {
    case 1: return obb.position + obb.half.x * obb.XAxis();
    case 2: return obb.position - obb.half.y * obb.YAxis();
    case 3: return obb.position + obb.half.y * obb.YAxis();
    case 4: return obb.position - obb.half.z * obb.ZAxis();
    case 5: return obb.position + obb.half.z * obb.ZAxis();
    case 0:
    default: {
      return obb.position - obb.half.x * obb.XAxis();
    }
  }
}

/// Return the surface normal of the given face of the OBB.
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER_OUTLINE vec<3, T, Q> faceNormalOBB(OBB<T, Q> const &obb, length_t faceIndex) {
  switch (faceIndex) {
    case 1: return obb.XAxis();
    case 2: return -obb.YAxis();
    case 3: return obb.YAxis();
    case 4: return -obb.ZAxis();
    case 5: return obb.ZAxis();
    case 0:
    default: {
      return -obb.XAxis();
    }
  }
}

/// Generate a point on the surface of the given face of the OBB.
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER_OUTLINE vec<3, T, Q> facePoint(OBB<T, Q> const &obb, length_t faceIndex, T u, T v) {
  const length_t fid = glm::clamp(faceIndex, length_t(0), length_t(5));
  const length_t uid = fid / 2;
  const length_t vid = (fid / 2 + 1) % 3;
  const mat<3, 3, T, Q> axis = obb.Axes();
  const vec<3, T, Q> u3 = axis[uid] * (T(2) * u - T(1)) * obb.half[uid];
  const vec<3, T, Q> v3 = axis[vid] * (T(2) * v - T(1)) * obb.half[vid];
  switch (fid) {
    case 1: return obb.position + obb.half.z * axis[2] + u3 + v3;
    case 2: return obb.position - obb.half.x * axis[0] + u3 + v3;
    case 3: return obb.position + obb.half.x * axis[0] + u3 + v3;
    case 4: return obb.position - obb.half.y * axis[1] + u3 + v3;
    case 5: return obb.position + obb.half.y * axis[1] + u3 + v3;
    case 0:
    default: {
      return obb.position - obb.half.z * axis[2] + u3 + v3;
    }
  }
}

/// Generate a Plane (point and normal) for the given face of the OBB.
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER Plane<3, T, Q> facePlane(OBB<T, Q> const &obb, length_t faceIndex) {
  return Plane<3, T, Q>(faceCentroid(obb, faceIndex), faceNormalOBB(obb, faceIndex));
}

/* Return the distance between the OBB and the given object. */

/// Return the distance between the OBB and given point
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(OBB<T, Q> const &obb, vec<3, T, Q> const &point) {
  const vec<3, T, Q> local = obb.Localize(point);
  return distance(local, clamp(local, -obb.half, obb.half));
}

/// Return the distance between the OBB and given Sphere
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(OBB<T, Q> const &obb, Sphere<3, T, Q> const &sphere) {
  return max(T(0), distance(obb, sphere.pos) - sphere.r);
}

/// Return the squared distance between the OBB and given point
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance2(OBB<T, Q> const &obb, vec<3, T, Q> const &point) {
  const vec<3, T, Q> local = obb.Localize(point);
  return distance2(local, clamp(local, -obb.half, obb.half));
}

/* Test if the given objects are fully contained inside the OBB */

/// OBB/Point containment test
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool contains(OBB<T, Q> const &obb, vec<3, T, Q> const &point, T distanceThreshold = contains_eps<T>()) {
  const vec<3, T, Q> local = obb.Localize(point);
  const vec<3, T, Q> half = obb.half + distanceThreshold;
  return all(greaterThanEqual(local, -half))
         && all(lessThanEqual(local, half));
}

/// OBB/AABB containment test
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool contains(OBB<T, Q> const &obb, AABB<3, T, Q> const &aabb, T distanceThreshold = contains_eps<T>()) {
  bool result = true;
  for (length_t i = 0; i < 8 && result; ++i)
    result &= contains(obb, cornerPoint(aabb, i), distanceThreshold);
  return result;
}

/// OBB/Segment containment test
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool contains(OBB<T, Q> const &obb, Segment<3, T, Q> const &segment, T distanceThreshold = contains_eps<T>()) {
  return contains(obb, segment.a, distanceThreshold)
         && contains(obb, segment.b, distanceThreshold);
}

/// OBB/Triangle containment test
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool contains(OBB<T, Q> const &obb, Triangle<3, T, Q> const &triangle, T distanceThreshold = contains_eps<T>()) {
  return contains(obb, triangle.a, distanceThreshold)
         && contains(obb, triangle.b, distanceThreshold)
         && contains(obb, triangle.c, distanceThreshold);
}

/// OBB/OBB containment test
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool contains(OBB<T, Q> const &obb, OBB<T, Q> const &other, T distanceThreshold = contains_eps<T>()) {
  bool result = true;
  for (length_t i = 0; i < 8 && result; ++i)
    result &= contains(obb, cornerPoint(other, i), distanceThreshold);
  return result;
}

#if 0
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool contains(OBB<T, Q> const &obb, Sphere<3, T, Q> const &sphere) { }
#endif

/// OBB/Polygon containment test
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool contains(OBB<T, Q> const &obb, Polygon<3, T, Q> const &polygon, T distanceThreshold = contains_eps<T>()) {
  bool result = true;
  const size_t size = polygon.size();
  for (size_t i = 0; i < size && result; ++i)
    result &= contains(obb, polygon[i], distanceThreshold);
  return true;
}

/// Return the closest point between the OBB and given point.
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<3, T, Q> closestPoint(OBB<T, Q> const &obb, vec<3, T, Q> const &point) {
  const vec<3, T, Q> dir = point - obb.position;
  const mat<3, 3, T, Q> axis = obb.Axes();
  vec<3, T, Q> result = obb.position;
  for (length_t i = 0; i < 3; ++i) {
    const vec<3, T, Q> a = axis[i];
    result += (glm::clamp(dot(dir, a), -obb.half[i], obb.half[i]) * a);
  }
  return result;
}

/* Functions to expand the OBB to enclose the given objects */

template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER OBB<T, Q> enclose(OBB<T, Q> const &obb, vec<3, T, Q> const &point) {
  const mat<3, 3, T, Q> axis = obb.Axes();
  OBB<T, Q> result = obb;
  vec<3, T, Q> p = point - result.position;
  for (length_t i = 0; i < 3; ++i) {
    const T daxis = dot(p, axis[i]);
    const T dobb = abs(daxis) - obb.half[i];
    if (dobb > T(0)) {
      const T expand = dobb * T(0.5);
      result.half[i] += expand;
      result.position += (daxis > T(0) ? T(1) : T(-1)) * (axis[i] * expand);
    }
  }
  return result;
}

#if 0
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER OBB<T, Q> enclose(OBB<T, Q> const &obb, Segment<3, T, Q> const &segment) { }

template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER OBB<T, Q> enclose(OBB<T, Q> const &obb, Triangle<3, T, Q> const &triangle) { }

template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER OBB<T, Q> enclose(OBB<T, Q> const &obb, Sphere<3, T, Q> const &sphere) { }

template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER OBB<T, Q> enclose(OBB<T, Q> const &obb, Capsule<3, T, Q> const &capsule) { }

template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER OBB<T, Q> enclose(OBB<T, Q> const &obb, AABB<3, T, Q> const &aabb) { }

template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER OBB<T, Q> enclose(OBB<T, Q> const &obb, Polygon<3, T, Q> const &polygon) { }

template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER OBB<T, Q> enclose(OBB<T, Q> const &obb, OBB<T, Q> const &other) { }
#endif

/* SAT Test helper functions */

template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER void intervalSAT(OBB<T, Q> const &obb, vec<3, T, Q> const &axis, T &min, T &max, T eps = T(0)) {
  const vec<3, T, Q> obbAxis = conjugate(obb.rotation) * axis;
  const T obbProjPos = dot(obb.position, axis);
  const T obbProjHalf = dot(obb.half, abs(obbAxis));
  min = obbProjPos - obbProjHalf - eps;
  max = obbProjPos + obbProjHalf + eps;
}

template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER void intervalSAT(Triangle<3, T, Q> const &triangle, vec<3, T, Q> const &axis, T &tmin, T &tmax, T eps = T(0)) {
  const T da = dot(axis, triangle.a);
  const T db = dot(axis, triangle.b);
  const T dc = dot(axis, triangle.c);
  tmin = min(da, min(db, dc)) - eps;
  tmax = max(da, max(db, dc)) + eps;
}

/* Test whether the OBB and the given object intersect */

/// OBB/Line intersection test, includes parametric line coordinates.
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(OBB<T, Q> const &obb, Line<3, T, Q> const &line, T &dNear, T &dFar) {
  GLM_GEOM_ASSUME(lessThanEqual(dNear, dFar), false);
  const Line<3, T, Q> localLine = worldToLocal(obb) * line;
  const AABB<3, T, Q> localAABB(vec<3, T, Q>(T(0)), size(obb));
  return intersects(localAABB, localLine, dNear, dFar);
}

/// OBB/Segment intersection test, includes parametric segment coordinates.
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(OBB<T, Q> const &obb, Segment<3, T, Q> const &line, T &dNear, T &dFar) {
  GLM_GEOM_ASSUME(lessThanEqual(dNear, dFar), false);
  GLM_GEOM_ASSUME(greaterThanEqual(dNear, -epsilon<T>()) && lessThanEqual(dFar, T(1) + epsilon<T>()), false);
  const Segment<3, T, Q> localLine = worldToLocal(obb) * line;
  const AABB<3, T, Q> localAABB(vec<3, T, Q>(T(0)), size(obb));
  return intersects(localAABB, localLine, dNear, dFar);
}

/// OBB/Ray intersection test, includes parametric ray coordinates.
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(OBB<T, Q> const &obb, Ray<3, T, Q> const &ray, T &dNear, T &dFar) {
  GLM_GEOM_ASSUME(lessThanEqual(dNear, dFar), false);
  GLM_GEOM_ASSUME(greaterThanEqual(dNear, -epsilon<T>()), false);
  const Ray<3, T, Q> localRay = worldToLocal(obb) * ray;
  const AABB<3, T, Q> localAABB(vec<3, T, Q>(T(0)), size(obb));
  return intersects(localAABB, localRay, dNear, dFar);
}

/// OBB/Sphere intersection test
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(OBB<T, Q> const &obb, Sphere<3, T, Q> const &sphere, T eps = epsilon<T>()) {
  const T r = sphere.r + eps;
  const vec<3, T, Q> pt = closestPoint(obb, sphere.pos);
  return distance2(sphere.pos, pt) <= r * r;
}

/// OBB/Plane intersection test
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(OBB<T, Q> const &obb, Plane<3, T, Q> const &plane) {
  const mat<3, 3, T, Q> axis = obb.Axes();
  const T len = obb.half.x * abs(dot(plane.normal, axis[0]))
                + obb.half.y * abs(dot(plane.normal, axis[1]))
                + obb.half.z * abs(dot(plane.normal, axis[2]));
  return distance(plane, obb.position) <= len;
}

/// OBB/OBB intersection test
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER_OUTLINE bool intersects(OBB<T, Q> const &obb, OBB<T, Q> const &other, T eps = T(0)) {
  auto OverlapOnAxis = [&obb, &other, &eps](vec<3, T, Q> const &axis) -> bool {
    T obbProjMin, obbProjMax, otherProjMin, otherProjMax;
    intervalSAT(obb, axis, obbProjMin, obbProjMax);
    intervalSAT(other, axis, otherProjMin, otherProjMax, eps);
    return obbProjMax >= otherProjMin && obbProjMin <= otherProjMax;
  };

  const mat<3, 3, T, Q> obbAxis = obb.Axes();
  const mat<3, 3, T, Q> otherAxis = other.Axes();
  return OverlapOnAxis(obbAxis[0])
         && OverlapOnAxis(obbAxis[1])
         && OverlapOnAxis(obbAxis[2])
         && OverlapOnAxis(otherAxis[0])
         && OverlapOnAxis(otherAxis[1])
         && OverlapOnAxis(otherAxis[2])
         && OverlapOnAxis(cross(obbAxis[0], otherAxis[0]))
         && OverlapOnAxis(cross(obbAxis[0], otherAxis[1]))
         && OverlapOnAxis(cross(obbAxis[0], otherAxis[2]))
         && OverlapOnAxis(cross(obbAxis[1], otherAxis[0]))
         && OverlapOnAxis(cross(obbAxis[1], otherAxis[1]))
         && OverlapOnAxis(cross(obbAxis[1], otherAxis[2]))
         && OverlapOnAxis(cross(obbAxis[2], otherAxis[0]))
         && OverlapOnAxis(cross(obbAxis[2], otherAxis[1]))
         && OverlapOnAxis(cross(obbAxis[2], otherAxis[2]));
}

/// OBB/Triangle intersection test
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(OBB<T, Q> const &obb, Triangle<3, T, Q> const &triangle, T eps = T(0)) {
  auto OverlapOnAxis = [&obb, &triangle, &eps](vec<3, T, Q> const &axis) -> bool {
    T obbProjMin, obbProjMax, triProjMin, triProjMax;
    intervalSAT(obb, axis, obbProjMin, obbProjMax);
    intervalSAT(triangle, axis, triProjMin, triProjMax, eps);
    return obbProjMax >= triProjMin && obbProjMin <= triProjMax;
  };

  const mat<3, 3, T, Q> obbAxis = obb.Axes();
  const vec<3, T, Q> t0 = triangle.b - triangle.a;
  const vec<3, T, Q> t1 = triangle.c - triangle.b;
  const vec<3, T, Q> t2 = triangle.a - triangle.c;
  return OverlapOnAxis(obbAxis[0])
         && OverlapOnAxis(obbAxis[1])
         && OverlapOnAxis(obbAxis[2])
         && OverlapOnAxis(cross(t0, t1))
         && OverlapOnAxis(cross(obbAxis[0], t0))
         && OverlapOnAxis(cross(obbAxis[0], t1))
         && OverlapOnAxis(cross(obbAxis[0], t2))
         && OverlapOnAxis(cross(obbAxis[1], t0))
         && OverlapOnAxis(cross(obbAxis[1], t1))
         && OverlapOnAxis(cross(obbAxis[1], t2))
         && OverlapOnAxis(cross(obbAxis[2], t0))
         && OverlapOnAxis(cross(obbAxis[2], t1))
         && OverlapOnAxis(cross(obbAxis[2], t2));
}

/// OBB/AABB intersection test
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER_OUTLINE bool intersects(OBB<T, Q> const &obb, AABB<3, T, Q> const &aabb) {
  return intersects(obb, OBB<T, Q>(aabb));
}

/// OBB/Line intersection test
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(OBB<T, Q> const &obb, Line<3, T, Q> const &line) {
  T dNear = -std::numeric_limits<T>::infinity();
  T dFar = std::numeric_limits<T>::infinity();
  return intersects(obb, line, dNear, dFar);
}

/// OBB/Segment intersection test
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(OBB<T, Q> const &obb, Segment<3, T, Q> const &line) {
  T dNear(0), dFar(1);
  return intersects(obb, line, dNear, dFar);
}

/// OBB/Ray intersection test
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(OBB<T, Q> const &obb, Ray<3, T, Q> const &ray) {
  T dNear = T(0);
  T dFar = std::numeric_limits<T>::infinity();
  return intersects(obb, ray, dNear, dFar);
}

#if 0
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(OBB<T, Q> const &obb, Capsule<3, T, Q> const &capsule) { }

template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(OBB<T, Q> const &obb, Polygon<3, T, Q> const &polygon) { }

template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(OBB<T, Q> const &obb, Triangle<3, T, Q> const &triangle) {
  const Triangle<3, T, Q> localTriangle = worldToLocal(obb) * triangle;
  const AABB<3, T, Q> localAABB(vec<3, T, Q>(T(0)), size(obb));
  return intersects(localAABB, localTriangle);
}
#endif

/// @}
GLM_GEOM_END_NAMESPACE
#if GLM_GEOM_TOSTRING
namespace glm {
  namespace detail {
    template<typename T, qualifier Q>
    struct compute_to_string<OBB<T, Q>> {
      GLM_GEOM_QUALIFIER std::string call(OBB<T, Q> const &obb) {
        return detail::format("OBB(%s, %s, %s)",
          glm::to_string(obb.position).c_str(),
          glm::to_string(obb.rotation).c_str(),
          glm::to_string(obb.half).c_str()
        );
      }
    };
  }
}
#endif
