/// @ref geom_aabb
/// @file geom/aabb.hpp
///
/// @defgroup geom_aabb AABB
/// @ingroup geom
///

#pragma once

#include "setup.hpp"
#include "line.hpp"
#include "segment.hpp"
#include "ray.hpp"
#include "triangle.hpp"
#include "capsule.hpp"

#include <glm/gtx/compatibility.hpp>
#include <glm/gtx/norm.hpp>
#include <glm/gtx/quaternion.hpp>

#include "../matrix_extensions.hpp"
#if GLM_MESSAGES == GLM_ENABLE && !defined(GLM_EXT_INCLUDED)
  #pragma message("GLM: GLM_EXT_GEOM_aabb extension included")
#endif

GLM_GEOM_BEGIN_NAMESPACE
/// @addtogroup geom_aabb
/// @{

/// <summary>
/// An axis-aligned bounding box.
/// </summary>
template<length_t L, typename T, qualifier Q>
struct AABB {

  // -- Implementation detail --

  typedef T value_type;
  typedef AABB<L, T, Q> type;
  typedef vec<L, T, Q> point_type;

  // -- Data --

  point_type minPoint;  // Minimum extent of this AABB
  point_type maxPoint;  // Maximum extent of this AABB

#if GLM_CONFIG_DEFAULTED_DEFAULT_CTOR == GLM_ENABLE
  AABB() GLM_DEFAULT_CTOR;
#else
  AABB()
  #if GLM_CONFIG_CTOR_INIT != GLM_CTOR_INIT_DISABLE
    : minPoint(T(-1)), maxPoint(T(1))
  #endif
  {
  }
#endif

  AABB(T scalar)
    : minPoint(scalar), maxPoint(scalar) {
  }

  AABB(point_type const &min, point_type const &max)
    : minPoint(min), maxPoint(max) {
  }

  AABB(AABB<L, T, Q> const &aabb)
    : minPoint(aabb.minPoint), maxPoint(aabb.maxPoint) {
  }

  AABB<L, T, Q> &operator=(AABB<L, T, Q> const &aabb) {
    minPoint = aabb.minPoint;
    maxPoint = aabb.maxPoint;
    return *this;
  }

  GLM_FUNC_QUALIFIER void setNegativeInfinity() {
    minPoint = vec<L, T, Q>(std::numeric_limits<T>::infinity());
    maxPoint = vec<L, T, Q>(-std::numeric_limits<T>::infinity());
  }

  GLM_FUNC_QUALIFIER void setFromCenterAndSize(point_type const &center, point_type const &size) {
    vec<L, T, Q> halfSize = T(0.5) * size;
    minPoint = center - halfSize;
    maxPoint = center + halfSize;
  }

  GLM_FUNC_QUALIFIER void enclose(point_type const &point) {
    minPoint = min(minPoint, point);
    maxPoint = max(maxPoint, point);
  }
};

template<typename T, qualifier Q, class Matrix>
GLM_GEOM_DECL AABB<3, T, Q> transformAsAABB(AABB<3, T, Q> const &aabb, Matrix const &m);

template<typename T, qualifier Q, class Matrix>
GLM_GEOM_DECL AABB<2, T, Q> transformAsAABB(AABB<2, T, Q> const &aabb, Matrix const &m);

template<length_t L, typename T, qualifier Q>
GLM_GEOM_DECL vec<L, T, Q> centroid(AABB<L, T, Q> const &aabb);

template<length_t L, typename T, qualifier Q>
static AABB<L, T, Q> operator-(AABB<L, T, Q> const &aabb) {
  return AABB<L, T, Q>(-aabb.maxPoint, -aabb.minPoint);
}

template<length_t L, typename T, qualifier Q>
static bool operator==(AABB<L, T, Q> const &a1, AABB<L, T, Q> const &a2) {
  return a1.minPoint == a2.minPoint && a1.maxPoint == a2.maxPoint;
}

template<length_t L, typename T, qualifier Q>
static AABB<L, T, Q> operator+(AABB<L, T, Q> const &aabb, vec<L, T, Q> const &point) {
  return AABB<L, T, Q>(aabb.minPoint + point, aabb.maxPoint + point);
}

template<length_t L, typename T, qualifier Q>
static AABB<L, T, Q> operator-(AABB<L, T, Q> const &aabb, vec<L, T, Q> const &point) {
  return AABB<L, T, Q>(aabb.minPoint - point, aabb.maxPoint - point);
}

template<length_t L, typename T, qualifier Q>
static AABB<L, T, Q> operator*(mat<3, 3, T, Q> const &m, AABB<L, T, Q> const &aabb) {
  return transformAsAABB<T, Q, mat<3, 3, T, Q>>(aabb, m);
}

template<length_t L, typename T, qualifier Q>
static AABB<L, T, Q> operator*(mat<3, 4, T, Q> const &m, AABB<L, T, Q> const &aabb) {
  return transformAsAABB<T, Q, mat<3, 4, T, Q>>(aabb, m);
}

template<length_t L, typename T, qualifier Q>
static AABB<L, T, Q> operator*(mat<4, 3, T, Q> const &m, AABB<L, T, Q> const &aabb) {
  return transformAsAABB<T, Q, mat<4, 3, T, Q>>(aabb, m);
}

template<length_t L, typename T, qualifier Q>
static AABB<L, T, Q> operator*(mat<4, 4, T, Q> const &m, AABB<L, T, Q> const &aabb) {
  return transformAsAABB<T, Q, mat<4, 4, T, Q>>(aabb, m);
}

template<length_t L, typename T, qualifier Q>
static AABB<L, T, Q> operator*(qua<T, Q> const &q, AABB<L, T, Q> const &aabb) {
  const vec<L, T, Q> center = q * centroid(aabb);
  const vec<L, T, Q> newDir = abs((q * size(aabb)) * T(0.5));
  return AABB<L, T, Q>(center - newDir, center + newDir);
}

template<typename T, qualifier Q>
static AABB<2, T, Q> operator*(qua<T, Q> const &q, AABB<2, T, Q> const &aabb) {
  return operator*(mat3_cast(q), aabb);
}

/// Returns the component-wise comparison of: |x.min - y.min| < eps && |x.max - y.max| < eps
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool equal(AABB<L, T, Q> const &x, AABB<L, T, Q> const &y, T eps = epsilon<T>()) {
  return all(equal(x.minPoint, y.minPoint, eps)) && all(equal(x.maxPoint, y.maxPoint, eps));
}

/// Returns the component-wise comparison of: |x.min - y.min| < eps && |x.max - y.max| < eps
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool equal(AABB<L, T, Q> const &x, AABB<L, T, Q> const &y, vec<L, T, Q> const &eps) {
  return all(equal(x.minPoint, y.minPoint, eps)) && all(equal(x.maxPoint, y.maxPoint, eps));
}

/// Returns the component-wise comparison between two AABBs in term of ULPs.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool equal(AABB<L, T, Q> const &x, AABB<L, T, Q> const &y, int MaxULPs) {
  return all(equal(x.minPoint, y.minPoint, MaxULPs)) && all(equal(x.maxPoint, y.maxPoint, MaxULPs));
}

/// Returns the component-wise comparison between two AABBs in term of ULPs.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool equal(AABB<L, T, Q> const &x, AABB<L, T, Q> const &y, vec<L, int, Q> const &MaxULPs) {
  return all(equal(x.minPoint, y.minPoint, MaxULPs)) && all(equal(x.maxPoint, y.maxPoint, MaxULPs));
}

/// Returns the component-wise comparison of: |x.min - y.min| >= eps || |x.max - y.max| >= eps
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool notEqual(AABB<L, T, Q> const &x, AABB<L, T, Q> const &y, T eps = epsilon<T>()) {
  return any(notEqual(x.minPoint, y.minPoint, eps)) || any(notEqual(x.maxPoint, y.maxPoint, eps));
}

/// Returns the component-wise comparison of: |x.min - y.min| >= eps || |x.max - y.max| >= eps
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool notEqual(AABB<L, T, Q> const &x, AABB<L, T, Q> const &y, vec<L, T, Q> const &eps) {
  return any(notEqual(x.minPoint, y.minPoint, eps)) || any(notEqual(x.maxPoint, y.maxPoint, eps));
}

/// Returns the component-wise comparison between two AABBs in term of ULPs.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool notEqual(AABB<L, T, Q> const &x, AABB<L, T, Q> const &y, int MaxULPs) {
  return any(notEqual(x.minPoint, y.minPoint, MaxULPs)) || any(notEqual(x.maxPoint, y.maxPoint, MaxULPs));
}

/// Returns the component-wise comparison between two AABBs in term of ULPs.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool notEqual(AABB<L, T, Q> const &x, AABB<L, T, Q> const &y, vec<L, int, Q> const &MaxULPs) {
  return any(notEqual(x.minPoint, y.minPoint, MaxULPs)) || any(notEqual(x.maxPoint, y.maxPoint, MaxULPs));
}

/// Create an AABB by specifying its center and size (along each dimension).
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER AABB<L, T, Q> aabbFromCenterAndSize(vec<L, T, Q> const &center, vec<L, T, Q> const &size) {
  const vec<L, T, Q> halfSize = T(0.5) * size;
  return AABB<L, T, Q>(center - halfSize, center + halfSize);
}

/// Create an AABB by specifying its center and size (uniform on each dimension).
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER AABB<L, T, Q> aabbFromCenterAndSize(vec<L, T, Q> const &center, T size) {
  const vec<L, T, Q> halfSize = vec<L, T, Q>{ T(0.5) * size };
  return AABB<L, T, Q>(center - halfSize, center + halfSize);
}

/// Create the smallest possible AABB, in terms of volume, that contains the provided sphere.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER AABB<L, T, Q> aabbFromSphere(Sphere<L, T, Q> const &sphere) {
  const vec<L, T, Q> d(sphere.r);
  return AABB<L, T, Q>(sphere.pos - d, sphere.pos + d);
}

/// Create the smallest possible AABB, in terms of volume, that contains the provided OBB.
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER AABB<3, T, Q> aabbFromOBB(OBB<T, Q> const &obb) {
  return aabbFromCenterAndSize<3, T, Q>(obb.position, diagonal(obb));
}

/// Test if any component of the AABB is infinite.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool isinf(AABB<L, T, Q> const &aabb) {
  return any(isinf(aabb.minPoint)) || any(isinf(aabb.maxPoint));
}

/// Test if any component of the AABB is NaN.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool isnan(AABB<L, T, Q> const &aabb) {
  return any(isnan(aabb.minPoint)) || any(isnan(aabb.maxPoint));
}

/// Test if all components of the AABB are finite.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool isfinite(AABB<L, T, Q> const &aabb) {
  return all(isfinite(aabb.minPoint)) && all(isfinite(aabb.maxPoint));
}

template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER T width(AABB<2, T, Q> const &aabb) {
  return aabb.maxPoint.x - aabb.minPoint.x;
}

template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER T height(AABB<2, T, Q> const &aabb) {
  return aabb.maxPoint.y - aabb.minPoint.y;
}

/// Return true if the AABB is degenerate, i.e., does not span in a strictly positive volume
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool isDegenerate(AABB<L, T, Q> const &aabb) {
  bool result = false;
  for (length_t i = 0; i < L; ++i)
    result |= !aabb.minPoint[i] >= aabb.maxPoint[i];
  return result;
}

template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool isDegenerate(AABB<3, T, Q> const &aabb) {
  return aabb.minPoint.x >= aabb.maxPoint.x || aabb.minPoint.y >= aabb.maxPoint.y || aabb.minPoint.z >= aabb.maxPoint.z;
}

template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool isDegenerate(AABB<2, T, Q> const &aabb) {
  return aabb.minPoint.x >= aabb.maxPoint.x || aabb.minPoint.y >= aabb.maxPoint.y;
}

template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool hasNegativeVolume(AABB<2, T, Q> const &aabb) {
  return aabb.maxPoint.x < aabb.minPoint.x || aabb.maxPoint.y < aabb.minPoint.y;
}

/// Return the center point of the AABB
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> centroid(AABB<L, T, Q> const &aabb) {
  return (aabb.minPoint + aabb.maxPoint) * T(0.5);
}

/// <summary>
/// Generates a point inside the AABB. "p" is a vector of normalized values
/// (i.e., between [0, 1]) along each axis, relative to the minpoint.
/// </summary>
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> pointInside(AABB<L, T, Q> const &aabb, vec<L, T, Q> const &p) {
  const vec<L, T, Q> d = aabb.maxPoint - aabb.minPoint;
  return aabb.minPoint + d * p;
}

/// Return the smallest Sphere that contains the AABB.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER Sphere<L, T, Q> minimalEnclosingSphere(AABB<L, T, Q> const &aabb) {
  return Sphere<L, T, Q>(centroid(aabb), length(aabb.maxPoint - aabb.minPoint) * T(0.5));
}

/// Return the largest Sphere that can fit inside the AABB.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER Sphere<L, T, Q> maximalContainedSphere(AABB<L, T, Q> const &aabb) {
  const vec<L, T, Q> hsize = halfSize(aabb);
  return Sphere<L, T, Q>(centroid(aabb), min(hsize.x, min(hsize.y, hsize.z)));
}

/// Return an edge of the AABB: [0, 11].
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER_OUTLINE Segment<3, T, Q> edge(AABB<3, T, Q> const &aabb, length_t edgeIndex) {
  switch (edgeIndex) {
    case 1: return Segment<3, T, Q>(aabb.minPoint, vec<3, T, Q>(aabb.minPoint.x, aabb.maxPoint.y, aabb.minPoint.z));
    case 2: return Segment<3, T, Q>(aabb.minPoint, vec<3, T, Q>(aabb.maxPoint.x, aabb.minPoint.y, aabb.minPoint.z));
    case 3: return Segment<3, T, Q>(vec<3, T, Q>(aabb.minPoint.x, aabb.minPoint.y, aabb.maxPoint.z), vec<3, T, Q>(aabb.minPoint.x, aabb.maxPoint.y, aabb.maxPoint.z));
    case 4: return Segment<3, T, Q>(vec<3, T, Q>(aabb.minPoint.x, aabb.minPoint.y, aabb.maxPoint.z), vec<3, T, Q>(aabb.maxPoint.x, aabb.minPoint.y, aabb.maxPoint.z));
    case 5: return Segment<3, T, Q>(vec<3, T, Q>(aabb.minPoint.x, aabb.maxPoint.y, aabb.minPoint.z), vec<3, T, Q>(aabb.minPoint.x, aabb.maxPoint.y, aabb.maxPoint.z));
    case 6: return Segment<3, T, Q>(vec<3, T, Q>(aabb.minPoint.x, aabb.maxPoint.y, aabb.minPoint.z), vec<3, T, Q>(aabb.maxPoint.x, aabb.maxPoint.y, aabb.minPoint.z));
    case 7: return Segment<3, T, Q>(vec<3, T, Q>(aabb.minPoint.x, aabb.maxPoint.y, aabb.maxPoint.z), aabb.maxPoint);
    case 8: return Segment<3, T, Q>(vec<3, T, Q>(aabb.maxPoint.x, aabb.minPoint.y, aabb.minPoint.z), vec<3, T, Q>(aabb.maxPoint.x, aabb.minPoint.y, aabb.maxPoint.z));
    case 9: return Segment<3, T, Q>(vec<3, T, Q>(aabb.maxPoint.x, aabb.minPoint.y, aabb.minPoint.z), vec<3, T, Q>(aabb.maxPoint.x, aabb.maxPoint.y, aabb.minPoint.z));
    case 10: return Segment<3, T, Q>(vec<3, T, Q>(aabb.maxPoint.x, aabb.minPoint.y, aabb.maxPoint.z), aabb.maxPoint);
    case 11: return Segment<3, T, Q>(vec<3, T, Q>(aabb.maxPoint.x, aabb.maxPoint.y, aabb.minPoint.z), aabb.maxPoint);
    case 0:
    default: {  // First Point
      return Segment<3, T, Q>(aabb.minPoint, vec<3, T, Q>(aabb.minPoint.x, aabb.minPoint.y, aabb.maxPoint.z));
    }
  }
}

template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER_OUTLINE Segment<2, T, Q> edge(AABB<2, T, Q> const &aabb, length_t edgeIndex) {
  switch (edgeIndex) {
    case 1: return Segment<2, T, Q>(vec<2, T, Q>(aabb.maxPoint.x, aabb.minPoint.y), aabb.maxPoint);
    case 2: return Segment<2, T, Q>(aabb.maxPoint, vec<2, T, Q>(aabb.minPoint.x, aabb.maxPoint.y));
    case 3: return Segment<2, T, Q>(vec<2, T, Q>(aabb.minPoint.x, aabb.maxPoint.y), aabb.minPoint);
    case 0:
    default: {  // First Point
      return Segment<2, T, Q>(aabb.minPoint, vec<2, T, Q>(aabb.maxPoint.x, aabb.minPoint.y));
    }
  }
}

/// Return a corner point of the AABB: [0, 7].
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER_OUTLINE vec<3, T, Q> cornerPoint(AABB<3, T, Q> const &aabb, length_t index) {
  return vec<3, T, Q>(
    (index & 4) ? aabb.maxPoint.x : aabb.minPoint.x,
    (index & 2) ? aabb.maxPoint.y : aabb.minPoint.y,
    (index & 1) ? aabb.maxPoint.z : aabb.minPoint.z
  );
}

/// Return a corner point of the AABB: [0, 3].
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER_OUTLINE vec<2, T, Q> cornerPoint(AABB<2, T, Q> const &aabb, length_t index) {
  return vec<2, T, Q>(
    (index & 2) ? aabb.maxPoint.x : aabb.minPoint.x,
    (index & 1) ? aabb.maxPoint.y : aabb.minPoint.y
  );
}

/// Return an extreme point along the AABB, i.e., the farthest point in a given direction.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> extremePoint(AABB<L, T, Q> const &aabb, vec<L, T, Q> const &direction) {
  vec<L, T, Q> result(T(0));
  for (length_t i = 0; i < L; ++i)
    result[i] = direction[i] >= T(0) ? aabb.maxPoint[i] : aabb.minPoint[i];
  return result;
}

template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<3, T, Q> extremePoint(AABB<3, T, Q> const &aabb, vec<3, T, Q> const &direction) {
  return vec<3, T, Q>(
    (direction.x >= T(0) ? aabb.maxPoint.x : aabb.minPoint.x),
    (direction.y >= T(0) ? aabb.maxPoint.y : aabb.minPoint.y),
    (direction.z >= T(0) ? aabb.maxPoint.z : aabb.minPoint.z)
  );
}

/// Return a point along an edge of the AABB.
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER_OUTLINE vec<3, T, Q> pointOnEdge(AABB<3, T, Q> const &aabb, length_t edgeIndex, T u) {
  const vec<3, T, Q> d = aabb.maxPoint - aabb.minPoint;
  switch (edgeIndex) {
    case 1: return vec<3, T, Q>(aabb.minPoint.x, aabb.maxPoint.y, aabb.minPoint.z + u * d.z);
    case 2: return vec<3, T, Q>(aabb.maxPoint.x, aabb.minPoint.y, aabb.minPoint.z + u * d.z);
    case 3: return vec<3, T, Q>(aabb.maxPoint.x, aabb.maxPoint.y, aabb.minPoint.z + u * d.z);
    case 4: return vec<3, T, Q>(aabb.minPoint.x, aabb.minPoint.y + u * d.y, aabb.minPoint.z);
    case 5: return vec<3, T, Q>(aabb.maxPoint.x, aabb.minPoint.y + u * d.y, aabb.minPoint.z);
    case 6: return vec<3, T, Q>(aabb.minPoint.x, aabb.minPoint.y + u * d.y, aabb.maxPoint.z);
    case 7: return vec<3, T, Q>(aabb.maxPoint.x, aabb.minPoint.y + u * d.y, aabb.maxPoint.z);
    case 8: return vec<3, T, Q>(aabb.minPoint.x + u * d.x, aabb.minPoint.y, aabb.minPoint.z);
    case 9: return vec<3, T, Q>(aabb.minPoint.x + u * d.x, aabb.minPoint.y, aabb.maxPoint.z);
    case 10: return vec<3, T, Q>(aabb.minPoint.x + u * d.x, aabb.maxPoint.y, aabb.minPoint.z);
    case 11: return vec<3, T, Q>(aabb.minPoint.x + u * d.x, aabb.maxPoint.y, aabb.maxPoint.z);
    case 0:  // First point
    default: {
      return vec<3, T, Q>(aabb.minPoint.x, aabb.minPoint.y, aabb.minPoint.z + u * d.z);
    }
  }
}

/// Return the point at the center of the given face, [0, 5], of the AABB.
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER_OUTLINE vec<3, T, Q> faceCentroid(AABB<3, T, Q> const &aabb, length_t faceIndex) {
  const vec<3, T, Q> center = (aabb.minPoint + aabb.maxPoint) * T(0.5);
  switch (faceIndex) {
    case 1: return vec<3, T, Q>(aabb.maxPoint.x, center.y, center.z);
    case 2: return vec<3, T, Q>(center.x, aabb.minPoint.y, center.z);
    case 3: return vec<3, T, Q>(center.x, aabb.maxPoint.y, center.z);
    case 4: return vec<3, T, Q>(center.x, center.y, aabb.minPoint.z);
    case 5: return vec<3, T, Q>(center.x, center.y, aabb.maxPoint.z);
    case 0:
    default: {
      return vec<3, T, Q>(aabb.minPoint.x, center.y, center.z);
    }
  }
}

/// Generate a point on the surface of the given face of the AABB.
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER_OUTLINE vec<3, T, Q> facePoint(AABB<3, T, Q> const &aabb, length_t faceIndex, T u, T v) {
  const vec<3, T, Q> d = aabb.maxPoint - aabb.minPoint;
  switch (faceIndex) {
    case 1: return vec<3, T, Q>(aabb.maxPoint.x, aabb.minPoint.y + u * d.y, aabb.minPoint.z + v * d.z);
    case 2: return vec<3, T, Q>(aabb.minPoint.x + u * d.x, aabb.minPoint.y, aabb.minPoint.z + v * d.z);
    case 3: return vec<3, T, Q>(aabb.minPoint.x + u * d.x, aabb.maxPoint.y, aabb.minPoint.z + v * d.z);
    case 4: return vec<3, T, Q>(aabb.minPoint.x + u * d.x, aabb.minPoint.y + v * d.y, aabb.minPoint.z);
    case 5: return vec<3, T, Q>(aabb.minPoint.x + u * d.x, aabb.minPoint.y + v * d.y, aabb.maxPoint.z);
    case 0:
    default: {
      return vec<3, T, Q>(aabb.minPoint.x, aabb.minPoint.y + u * d.y, aabb.minPoint.z + v * d.z);
    }
  }
}

/// Return the surface normal of the given face of the AABB.
template<typename T, qualifier Q = glm::defaultp>
GLM_GEOM_QUALIFIER_OUTLINE vec<3, T, Q> faceNormalAABB(length_t faceIndex) {
  switch (faceIndex) {
    case 1: return vec<3, T, Q>(T(1), T(0), T(0));
    case 2: return vec<3, T, Q>(T(0), T(-1), T(0));
    case 3: return vec<3, T, Q>(T(0), T(1), T(0));
    case 4: return vec<3, T, Q>(T(0), T(0), T(-1));
    case 5: return vec<3, T, Q>(T(0), T(0), T(1));
    case 0:
    default: {
      return vec<3, T, Q>(T(-1), T(0), T(0));
    }
  }
}

/// Generate a Plane (point and normal) for the given face of the AABB.
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER Plane<3, T, Q> facePlane(AABB<3, T, Q> const &aabb, length_t faceIndex) {
  return Plane<3, T, Q>(faceCentroid(aabb, faceIndex), faceNormalAABB<T, Q>(faceIndex));
}

/// Generate an AABB that encloses the given set of points.
template<length_t L, typename T, qualifier Q, class Vector>
GLM_GEOM_QUALIFIER AABB<L, T, Q> minimalEnclosingAABB(Vector const &points) {
  AABB<L, T, Q> aabb(T(0));
  aabb.setNegativeInfinity();
  for (auto it = points.begin(); it != points.end(); ++it)
    aabb.enclose(*it);
  return aabb;
}

/// Generate an AABB that encloses the given point iterators.
template<class Iterator, length_t L, typename T, qualifier Q = glm::defaultp>
GLM_GEOM_QUALIFIER AABB<L, T, Q> minimalEnclosingAABB(Iterator begin, Iterator const &end) {
  AABB<L, T, Q> aabb(T(0));
  aabb.setNegativeInfinity();
  for (; begin != end; ++begin)
    aabb.enclose(*begin);
  return aabb;
}

/// Return the length of the AABB along each dimension.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> size(AABB<L, T, Q> const &aabb) {
  return aabb.maxPoint - aabb.minPoint;
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> halfSize(AABB<L, T, Q> const &aabb) {
  return size(aabb) * T(0.5);
}

/// Return the volume of the AABB.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T volume(AABB<L, T, Q> const &aabb) {
  const vec<L, T, Q> s = size(aabb);

  T result(1);
  for (length_t i = 0; i < L; ++i)
    result *= s[i];
  return result;
}

template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER T volume(AABB<3, T, Q> const &aabb) {
  const vec<3, T, Q> s = size(aabb);
  return s.x * s.y * s.z;
}

/// Return the surface area of the faces of the AABB.
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER T surfaceArea(AABB<3, T, Q> const &aabb) {
  const vec<3, T, Q> s = size(aabb);
  return T(2) * (s.x * s.y + s.x * s.z + s.y * s.z);
}

/// Apply a uniform scale to the AABB
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER AABB<L, T, Q> scale(AABB<L, T, Q> const &aabb, vec<L, T, Q> const &centerPoint, T scaleFactor) {
  return AABB<L, T, Q>(
    (aabb.minPoint - centerPoint) * scaleFactor + centerPoint,
    (aabb.maxPoint - centerPoint) * scaleFactor + centerPoint
  );
}

/// Grow an AABB by the given size.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER AABB<L, T, Q> grow(AABB<L, T, Q> const &aabb, T amount) {
  return AABB<L, T, Q>(aabb.minPoint - (T(0.5) * amount), aabb.maxPoint + (T(0.5) * amount));
}

/// Project the AABB onto the provided axis.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER void projectToAxis(AABB<L, T, Q> const &aabb, vec<L, T, Q> const &axis, T &dMin, T &dMax) {
  const vec<L, T, Q> c = (aabb.minPoint + aabb.maxPoint) * T(0.5);
  const vec<L, T, Q> e = aabb.maxPoint - c;
  const T r = abs(dot(e, abs(axis)));
  const T s = dot(axis, c);  // distance center/plane.
  dMin = s - r;
  dMax = s + r;
}

/// Compute a 30-bit morton code for a point within an AABB
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER unsigned int morton3D(AABB<L, T, Q> const &aabb, vec<L, T, Q> const &point) {
  return morton3D((point - aabb.minPoint) / size(aabb));
}

/// Apply a generic matrix transformation to the AABB.
template<typename T, qualifier Q, class Matrix>
GLM_GEOM_QUALIFIER_OUTLINE AABB<3, T, Q> transformAsAABB(AABB<3, T, Q> const &aabb, Matrix const &m) {
  const vec<3, T, Q> centroid = (aabb.minPoint + aabb.maxPoint) * T(0.5);
  const vec<3, T, Q> halfSize = centroid - aabb.minPoint;
  const vec<3, T, Q> newCenter = transformPos(m, centroid);
  const vec<3, T, Q> newDir = vec<3, T, Q>(
    abs(m[0][0] * halfSize.x) + abs(m[1][0] * halfSize.y) + abs(m[2][0] * halfSize.z),
    abs(m[0][1] * halfSize.x) + abs(m[1][1] * halfSize.y) + abs(m[2][1] * halfSize.z),
    abs(m[0][2] * halfSize.x) + abs(m[1][2] * halfSize.y) + abs(m[2][2] * halfSize.z)
  );
  return AABB<3, T, Q>(newCenter - newDir, newCenter + newDir);
}

/// Apply a generic matrix transformation to the AABB.
template<typename T, qualifier Q, class Matrix>
GLM_GEOM_QUALIFIER_OUTLINE AABB<2, T, Q> transformAsAABB(AABB<2, T, Q> const &aabb, Matrix const &m) {
  const T ax = m[0][0] * aabb.minPoint.x;
  const T bx = m[0][0] * aabb.maxPoint.x;
  const T ay = m[1][0] * aabb.minPoint.y;
  const T by = m[1][0] * aabb.maxPoint.y;
  const T ax2 = m[0][1] * aabb.minPoint.x;
  const T bx2 = m[0][1] * aabb.maxPoint.x;
  const T ay2 = m[1][1] * aabb.minPoint.y;
  const T by2 = m[1][1] * aabb.maxPoint.y;

  T ox = T(0), oy = T(0);
  GLM_IF_CONSTEXPR(m.length() > 3) {
    ox = m[3][0];
    oy = m[3][1];
  }
  return AABB<2, T, Q>(
    vec<2, T, Q>(min(ax, bx) + min(ay, by) + ox, min(ax2, bx2) + min(ay2, by2) + oy),
    vec<2, T, Q>(max(ax, bx) + max(ay, by) + ox, max(ax2, bx2) + max(ay2, by2) + oy)
  );
}

/* Return the closest point between an AABB and the given object */

/// Return the closest point between the AABB and given point.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(AABB<L, T, Q> const &aabb, vec<L, T, Q> const &target) {
  return clamp(target, aabb.minPoint, aabb.maxPoint);
}

#if 0
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(AABB<L, T, Q> const &aabb, AABB<L, T, Q> const &other) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(AABB<L, T, Q> const &aabb, Line<L, T, Q> const &line) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(AABB<L, T, Q> const &aabb, Segment<L, T, Q> const &segment) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(AABB<L, T, Q> const &aabb, Ray<L, T, Q> const &ray) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(AABB<L, T, Q> const &aabb, Plane<L, T, Q> const &plane) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(AABB<L, T, Q> const &aabb, Triangle<L, T, Q> const &triangle) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(AABB<L, T, Q> const &aabb, Sphere<L, T, Q> const &sphere) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(AABB<L, T, Q> const &aabb, Capsule<L, T, Q> const &capsule) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(AABB<L, T, Q> const &aabb, Disc<L, T, Q> const &disc) { }

template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<3, T, Q> closestPoint(AABB<3, T, Q> const &aabb, OBB<T, Q> const &obb) { }
#endif

/* Return the distance between the AABB and the given object. */

/// Return the distance between the AABB and given point
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(AABB<L, T, Q> const &aabb, vec<L, T, Q> const &point) {
  return distance(closestPoint(aabb, point), point);
}

/// Return the distance between the AABB and given Sphere
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(AABB<L, T, Q> const &aabb, Sphere<L, T, Q> const &sphere) {
  return max(T(0), distance(aabb, sphere.pos) - sphere.r);
}

#if 0
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(AABB<L, T, Q> const &aabb, AABB<L, T, Q> const &other) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(AABB<L, T, Q> const &aabb, Line<L, T, Q> const &line) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(AABB<L, T, Q> const &aabb, Segment<L, T, Q> const &segment) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(AABB<L, T, Q> const &aabb, Ray<L, T, Q> const &ray) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(AABB<L, T, Q> const &aabb, Plane<L, T, Q> const &plane) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(AABB<L, T, Q> const &aabb, Triangle<L, T, Q> const &triangle) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(AABB<L, T, Q> const &aabb, Capsule<L, T, Q> const &capsule) { }

template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(AABB<3, T, Q> const &aabb, OBB<T, Q> const &obb) { }
#endif

/* Test if the given objects are fully contained inside the AABB */

/// AABB/Point containment test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool contains(AABB<L, T, Q> const &aabb, vec<L, T, Q> const &target) {
  bool result = true;
  for (length_t i = 0; i < L; ++i)
    result &= (aabb.minPoint[i] <= target[i] && target[i] <= aabb.maxPoint[i]);
  return result;
}

/// AABB/Span containment test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool contains(AABB<L, T, Q> const &aabb, vec<L, T, Q> const &minPoint, vec<L, T, Q> const &maxPoint) {
  bool result = true;
  for (length_t i = 0; i < L; ++i)
    result &= (aabb.minPoint[i] <= minPoint[i] && maxPoint[i] <= aabb.maxPoint[i]);
  return result;
}

template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool contains(AABB<3, T, Q> const &aabb, vec<3, T, Q> const &point) {
  return aabb.minPoint.x <= point.x && point.x <= aabb.maxPoint.x
         && aabb.minPoint.y <= point.y && point.y <= aabb.maxPoint.y
         && aabb.minPoint.z <= point.z && point.z <= aabb.maxPoint.z;
}

template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool contains(AABB<3, T, Q> const &aabb, vec<3, T, Q> const &minPoint, vec<3, T, Q> const &maxPoint) {
  return aabb.minPoint.x <= minPoint.x && maxPoint.x <= aabb.maxPoint.x
         && aabb.minPoint.y <= minPoint.y && maxPoint.y <= aabb.maxPoint.y
         && aabb.minPoint.z <= minPoint.z && maxPoint.z <= aabb.maxPoint.z;
}

template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool contains(AABB<2, T, Q> const &aabb, vec<2, T, Q> const &point) {
  return aabb.minPoint.x <= point.x && point.x <= aabb.maxPoint.x
         && aabb.minPoint.y <= point.y && point.y <= aabb.maxPoint.y;
}

template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool contains(AABB<2, T, Q> const &aabb, vec<2, T, Q> const &minPoint, vec<2, T, Q> const &maxPoint) {
  return aabb.minPoint.x <= minPoint.x && maxPoint.x <= aabb.maxPoint.x
         && aabb.minPoint.y <= minPoint.y && maxPoint.y <= aabb.maxPoint.y;
}

/// AABB/AABB containment test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool contains(AABB<L, T, Q> const &aabb, AABB<L, T, Q> const &otheraabb) {
  return contains(aabb, otheraabb.minPoint, otheraabb.maxPoint);
}

/// AABB/Segment containment test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool contains(AABB<L, T, Q> const &aabb, Segment<L, T, Q> const &segment) {
  return contains(aabb, min(segment.a, segment.b), max(segment.a, segment.b));
}

/// AABB/Triangle containment test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool contains(AABB<L, T, Q> const &aabb, Triangle<L, T, Q> const &triangle) {
  return contains(aabb, boundingAABB(triangle));
}

/// AABB/Sphere containment test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool contains(AABB<L, T, Q> const &aabb, Sphere<L, T, Q> const &sphere) {
  const vec<L, T, Q> dir(sphere.r);
  return contains(aabb, sphere.pos - dir, sphere.pos + dir);
}

#if 0
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool contains(AABB<L, T, Q> const &aabb, Disc<L, T, Q> const &disc) { }
#endif

/// AABB/Capsule containment test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool contains(AABB<L, T, Q> const &aabb, Capsule<L, T, Q> const &capsule) {
  return contains(aabb, minimalEnclosingAABB(capsule));
}

/// AABB/Polygon containment test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool contains(AABB<L, T, Q> const &aabb, Polygon<L, T, Q> const &polygon) {
  return contains(aabb, minimalEnclosingAABB(polygon));
}

#if 0
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool contains(AABB<3, T, Q> const &aabb, OBB<T, Q> const &obb) { }
#endif

/* Functions to expand the AABB to enclose the given objects. */

/// Expand the AABB to enclose the given point
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER AABB<L, T, Q> enclose(AABB<L, T, Q> const &aabb, vec<L, T, Q> const &point) {
  AABB<L, T, Q> result(aabb);
  result.enclose(point);
  return result;
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER AABB<L, T, Q> enclose(AABB<L, T, Q> const &aabb, vec<L, T, Q> const &aabbMinPoint, vec<L, T, Q> const &aabbMaxPoint) {
  AABB<L, T, Q> result(aabb);
  result.enclose(aabbMinPoint);
  result.enclose(aabbMaxPoint);
  return result;
}

/// Expand the AABB to enclose the given Segment
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER AABB<L, T, Q> enclose(AABB<L, T, Q> const &aabb, Segment<L, T, Q> const &segment) {
  return enclose(aabb, min(segment.a, segment.b), max(segment.a, segment.b));
}

/// Expand the AABB to enclose the given Triangle
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER AABB<L, T, Q> enclose(AABB<L, T, Q> const &aabb, Triangle<L, T, Q> const &triangle) {
  return enclose(aabb, min(triangle.a, triangle.b, triangle.c), max(triangle.a, triangle.b, triangle.c));
}

/// Expand the AABB to enclose the given Sphere
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER AABB<L, T, Q> enclose(AABB<L, T, Q> const &aabb, Sphere<L, T, Q> const &sphere) {
  const vec<L, T, Q> d(sphere.r);
  return enclose(aabb, sphere.pos - d, sphere.pos + d);
}

/// Expand the AABB to enclose the given Capsule
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER AABB<L, T, Q> enclose(AABB<L, T, Q> const &aabb, Capsule<L, T, Q> const &capsule) {
  const vec<L, T, Q> d(capsule.r);
  return enclose(aabb, min(capsule.l.a, capsule.l.b) - d, max(capsule.l.a, capsule.l.b) + d);
}

/// Expand the AABB to enclose the given AABB
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER AABB<L, T, Q> enclose(AABB<L, T, Q> const &aabb, AABB<L, T, Q> const &other) {
  return enclose(aabb, other.minPoint, other.maxPoint);
}

/// Expand the AABB to enclose the given Polygon
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER AABB<L, T, Q> enclose(AABB<L, T, Q> const &aabb, Polygon<L, T, Q> const &polygon) {
  return enclose(aabb, minimalEnclosingAABB(polygon));
}

#if 0
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER AABB<3, T, Q> enclose(AABB<3, T, Q> const &aabb, OBB<T, Q> const &obb) { }
#endif

/// Expand the AABB in both directions by the given delta
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER AABB<L, T, Q> expand(AABB<L, T, Q> const &aabb, vec<L, T, Q> const &delta) {
  return AABB<L, T, Q>(aabb.minPoint - delta, aabb.maxPoint + delta);
}

/// Clamp the AABB to be contained within the specified (other) AABB.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER AABB<L, T, Q> clamp(AABB<L, T, Q> const &aabb, AABB<L, T, Q> const &other) {
  return AABB<L, T, Q>(
    clamp(aabb.minPoint, other.minPoint, other.maxPoint),
    clamp(aabb.maxPoint, other.minPoint, other.maxPoint)
  );
}

/// Generalized compute the intersection of a Line (or Ray) and the AABB.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersectLineAABB(AABB<L, T, Q> const &aabb, Line<L, T, Q> const &line, T &dNear, T &tFar) {
  for (length_t i = 0; i < L; ++i) {  // test each cardinal plane
    if (!detail::approx_zero(line.dir[i])) {
      const T recipDir = T(1) / line.dir[i];
      const T t1 = (aabb.minPoint[i] - line.pos[i]) * recipDir;
      const T t2 = (aabb.maxPoint[i] - line.pos[i]) * recipDir;
      if (t1 < t2)
        dNear = max(t1, dNear), tFar = min(t2, tFar);
      else  // Swap t1 and t2.
        dNear = max(t2, dNear), tFar = min(t1, tFar);

      if (dNear > tFar) {
        return false;  // exit before entering; AABB missed.
      }
    }
    else if (line.pos[i] < aabb.minPoint[i] || line.pos[i] > aabb.maxPoint[i]) {
      return false;  // AABB completely missed.
    }
  }
  return dNear <= tFar;
}

/// Return the intersection of a Line (or Ray) and the AABB.
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersectLineAABB(AABB<3, T, Q> const &aabb, Line<3, T, Q> const &line, T &dNear, T &tFar) {
  if (!detail::approx_zero(line.dir.x)) {  // test each cardinal plane: X, Y and Z.
    const T recipDir = T(1) / line.dir.x;
    const T t1 = (aabb.minPoint.x - line.pos.x) * recipDir;
    const T t2 = (aabb.maxPoint.x - line.pos.x) * recipDir;
    if (t1 < t2)
      dNear = max(t1, dNear), tFar = min(t2, tFar);
    else  // Swap t1 and t2.
      dNear = max(t2, dNear), tFar = min(t1, tFar);

    if (dNear > tFar)
      return false;  // exit before entering; AABB missed.
  }
  else if (line.pos.x < aabb.minPoint.x || line.pos.x > aabb.maxPoint.x)
    return false;  // AABB completely missed.

  if (!detail::approx_zero(line.dir.y)) {
    const T recipDir = T(1) / line.dir.y;
    const T t1 = (aabb.minPoint.y - line.pos.y) * recipDir;
    const T t2 = (aabb.maxPoint.y - line.pos.y) * recipDir;
    if (t1 < t2)
      dNear = max(t1, dNear), tFar = min(t2, tFar);
    else
      dNear = max(t2, dNear), tFar = min(t1, tFar);

    if (dNear > tFar)
      return false;
  }
  else if (line.pos.y < aabb.minPoint.y || line.pos.y > aabb.maxPoint.y)
    return false;

  if (!detail::approx_zero(line.dir.z)) {
    const T recipDir = T(1) / line.dir.z;
    const T t1 = (aabb.minPoint.z - line.pos.z) * recipDir;
    const T t2 = (aabb.maxPoint.z - line.pos.z) * recipDir;
    if (t1 < t2)
      dNear = max(t1, dNear), tFar = min(t2, tFar);
    else
      dNear = max(t2, dNear), tFar = min(t1, tFar);
  }
  else if (line.pos.z < aabb.minPoint.z || line.pos.z > aabb.maxPoint.z)
    return false;

  return dNear <= tFar;
}

/* Test whether the AABB and the given object intersect */

/// AABB/Line intersection test, includes parametric line coordinates.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(AABB<L, T, Q> const &aabb, Line<L, T, Q> const &line, T &dNear, T &dFar) {
  GLM_GEOM_ASSUME(lessThanEqual(dNear, dFar), false);
  return intersectLineAABB(aabb, line, dNear, dFar);
}

/// AABB/Line intersection test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(AABB<L, T, Q> const &aabb, Line<L, T, Q> const &line) {
  T dNear = -std::numeric_limits<T>::infinity();
  T dFar = std::numeric_limits<T>::infinity();
  return intersects(aabb, line, dNear, dFar);
}

/// AABB/Ray intersection test, includes parametric ray coordinate.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(AABB<L, T, Q> const &aabb, Ray<L, T, Q> const &ray, T &dNear, T &dFar) {
  GLM_GEOM_ASSUME(lessThanEqual(dNear, dFar), false);
  GLM_GEOM_ASSUME(greaterThanEqual(dNear, -epsilon<T>()), false);
  return intersectLineAABB(aabb, ray, dNear, dFar);
}

/// AABB/Ray intersection test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(AABB<L, T, Q> const &aabb, Ray<L, T, Q> const &ray) {
  T dNear = T(0);
  T dFar = std::numeric_limits<T>::infinity();
  return intersects(aabb, ray, dNear, dFar);
}

/// AABB/Segment intersection test, includes parametric segment coordinates.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(AABB<L, T, Q> const &aabb, Segment<L, T, Q> const &segment, T &dNear, T &dFar) {
  GLM_GEOM_ASSUME(lessThanEqual(dNear, dFar), false);
  GLM_GEOM_ASSUME(greaterThanEqual(dNear, -epsilon<T>()) && lessThanEqual(dFar, T(1) + epsilon<T>()), false);
  const vec<L, T, Q> dir = segment.dir2();
  const T len = length(dir);
  if (len <= epsilon<T>()) {  // Degenerate segment
    dNear = T(0), dFar = T(1);
    return contains(aabb, segment.a);
  }

  T sNear = dNear * len, sFar = dFar * len;
  if (intersectLineAABB(aabb, Line<L, T, Q>(segment.a, dir / len), sNear, sFar)) {
    dNear = sNear / len;
    dFar = sFar / len;
    return true;
  }
  return false;
}

/// AABB/Segment intersection test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(AABB<L, T, Q> const &aabb, Segment<L, T, Q> const &segment) {
  T dNear(0), dFar(1);
  return intersects(aabb, segment, dNear, dFar);
}

/// AABB/AABB intersection test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(AABB<L, T, Q> const &aabb, AABB<L, T, Q> const &other) {
  bool result = true;
  for (length_t i = 0; i < L; ++i)
    result &= aabb.minPoint[i] < other.maxPoint[i] && other.minPoint[i] < aabb.maxPoint[i];
  return result;
}

/// AABB/Sphere intersection test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(AABB<L, T, Q> const &aabb, Sphere<L, T, Q> const &sphere, T eps = epsilon<T>()) {
  const T r = sphere.r + eps;
  const vec<L, T, Q> pt = closestPoint(aabb, sphere.pos);
  return distance2(sphere.pos, pt) <= r * r;
}

/// AABB/OBB intersection test
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(AABB<3, T, Q> const &aabb, OBB<T, Q> const &obb) {
  return intersects(obb, aabb);
}

#if 0
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(AABB<L, T, Q> const &aabb, Disc<L, T, Q> const &disc) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects((AABB<L, T, Q> const &aabb, Capsule<L, T, Q> const &capsule) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(AABB<L, T, Q> const &aabb, Triangle<L, T, Q> const &triangle) { }
#endif

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(AABB<L, T, Q> const &aabb, Plane<L, T, Q> const &plane) {
  return intersects(plane, aabb);
}

/// AABB/Ray intersection test using the slab method. For performance, the ray
/// direction must be its element-wise inverse to its actual direction.
/// @see <a href="https://tavianator.com/2011/ray_box.html">Fast, Branchless Ray/Bounding Box Intersections</a>
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool slabs(AABB<3, T, Q> const &aabb, Ray<3, T, Q> const &ray, T &dNear, T &dFar) {
  const vec<3, T, Q> tMin = (aabb.minPoint - ray.pos) * /* INV */ ray.dir;
  const vec<3, T, Q> tMax = (aabb.maxPoint - ray.pos) * /* INV */ ray.dir;
  const vec<3, T, Q> t1 = min(tMin, tMax);
  const vec<3, T, Q> t2 = max(tMin, tMax);
  dNear = max(max(t1.x, t1.y), t1.z);
  dFar = min(min(t2.x, t2.y), t2.z);
  return dNear <= dFar;
}

/// AABB/Ray intersection test using the slab method
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool slabs(AABB<3, T, Q> const &aabb, Ray<3, T, Q> const &ray) {
  T dNear = T(0), dFar = std::numeric_limits<T>::infinity();
  return slabs(aabb, ray, dNear, dFar);
}

/* Specializations */

template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(AABB<3, T, Q> const &aabb, AABB<3, T, Q> const &other) {
  return aabb.minPoint.x < other.maxPoint.x
         && aabb.minPoint.y < other.maxPoint.y
         && aabb.minPoint.z < other.maxPoint.z
         && other.minPoint.x < aabb.maxPoint.x
         && other.minPoint.y < aabb.maxPoint.y
         && other.minPoint.z < aabb.maxPoint.z;
}

template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(AABB<2, T, Q> const &aabb, AABB<2, T, Q> const &other) {
  return aabb.minPoint.x < other.maxPoint.x
         && aabb.minPoint.y < other.maxPoint.y
         && other.minPoint.x < aabb.maxPoint.x
         && other.minPoint.y < aabb.maxPoint.y;
}

/// @}
GLM_GEOM_END_NAMESPACE
#if GLM_GEOM_TOSTRING
namespace glm {
  namespace detail {
    template<glm::length_t L, typename T, qualifier Q>
    struct compute_to_string<AABB<L, T, Q>> {
      GLM_GEOM_QUALIFIER std::string call(AABB<L, T, Q> const &aabb) {
        return detail::format("aabb(%s, %s)",
          glm::to_string(aabb.minPoint).c_str(),
          glm::to_string(aabb.maxPoint).c_str()
        );
      }
    };
  }
}
#endif
