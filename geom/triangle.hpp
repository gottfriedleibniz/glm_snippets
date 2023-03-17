/// @ref geom_triangle
/// @file geom/triangle.hpp
///
/// @defgroup geom_triangle Triangle
/// @ingroup geom
///

#pragma once

#include "setup.hpp"
#include "line.hpp"
#include "segment.hpp"
#include "plane.hpp"
#include "ray.hpp"

#include <glm/gtx/compatibility.hpp>
#include <glm/ext/scalar_relational.hpp>
#if GLM_MESSAGES == GLM_ENABLE && !defined(GLM_EXT_INCLUDED)
  #pragma message("GLM: GLM_EXT_GEOM_triangle extension included")
#endif

GLM_GEOM_BEGIN_NAMESPACE
/// @addtogroup geom_triangle
/// @{

/// <summary>
/// Three vertices in N-dimension space (N >= 2).
/// </summary>
template<length_t L, typename T, qualifier Q>
struct Triangle {

  // -- Implementation detail --

  typedef T value_type;
  typedef Triangle<L, T, Q> type;
  typedef vec<L, T, Q> point_type;

  // -- Data --

  point_type a;
  point_type b;
  point_type c;

#if GLM_CONFIG_DEFAULTED_DEFAULT_CTOR == GLM_ENABLE
  Triangle() GLM_DEFAULT_CTOR;
#else
  Triangle()
  #if GLM_CONFIG_CTOR_INIT != GLM_CTOR_INIT_DISABLE
    : a(T(0)), b(T(0)), c(T(0))
  #endif
  {
  }
#endif

  Triangle(T scalar)
    : a(scalar), b(scalar), c(scalar) {
  }

  Triangle(Triangle<L, T, Q> const &tri)
    : a(tri.a), b(tri.b), c(tri.c) {
  }

  Triangle(point_type const &a_, point_type const &b_, point_type const &c_)
    : a(a_), b(b_), c(c_) {
  }

  Triangle<L, T, Q> &operator=(Triangle<L, T, Q> const &t) {
    a = t.a;
    b = t.b;
    c = t.c;
    return *this;
  }
};

// operators

template<length_t L, typename T, qualifier Q>
static Triangle<L, T, Q> operator-(Triangle<L, T, Q> const &t) {
  return Triangle<L, T, Q>(t.a, t.c, t.b);
}

template<length_t L, typename T, qualifier Q>
static bool operator==(Triangle<L, T, Q> const &t1, Triangle<L, T, Q> const &t2) {
  return t1.a == t1.a && t1.b == t2.b && t1.c == t2.c;
}

template<length_t L, typename T, qualifier Q>
static Triangle<L, T, Q> operator+(Triangle<L, T, Q> const &t, vec<L, T, Q> const &offset) {
  return Triangle<L, T, Q>(t.a + offset, t.b + offset, t.c + offset);
}

template<length_t L, typename T, qualifier Q>
static Triangle<L, T, Q> operator-(Triangle<L, T, Q> const &t, vec<L, T, Q> const &offset) {
  return Triangle<L, T, Q>(t.a - offset, t.b - offset, t.c - offset);
}

template<typename T, qualifier Q>
static Triangle<3, T, Q> operator*(mat<3, 3, T, Q> const &m, Triangle<3, T, Q> const &t) {
  return Triangle<3, T, Q>(m * t.a, m * t.b, m * t.c);
}

template<typename T, qualifier Q>
static Triangle<3, T, Q> operator*(mat<3, 4, T, Q> const &m, Triangle<3, T, Q> const &t) {
  return Triangle<3, T, Q>(m * t.a, m * t.b, m * t.c);
}

template<typename T, qualifier Q>
static Triangle<3, T, Q> operator*(mat<4, 3, T, Q> const &m, Triangle<3, T, Q> const &t) {
  return Triangle<3, T, Q>(transformPos(m, t.a), transformPos(m, t.b), transformPos(m, t.c));
}

template<typename T, qualifier Q>
static Triangle<3, T, Q> operator*(mat<4, 4, T, Q> const &m, Triangle<3, T, Q> const &t) {
  return Triangle<3, T, Q>(transformPos(m, t.a), transformPos(m, t.b), transformPos(m, t.c));
}

template<typename T, qualifier Q>
static Triangle<3, T, Q> operator*(qua<T, Q> const &q, Triangle<3, T, Q> const &t) {
  return Triangle<3, T, Q>(q * t.a, q * t.b, q * t.c);
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool equal(Triangle<L, T, Q> const &x, Triangle<L, T, Q> const &y, T eps = epsilon<T>()) {
  return all(equal(x.a, y.a, eps)) && all(equal(x.b, y.b, eps)) && all(equal(x.c, y.c, eps));
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool equal(Triangle<L, T, Q> const &x, Triangle<L, T, Q> const &y, vec<L, T, Q> const &eps) {
  return all(equal(x.a, y.a, eps)) && all(equal(x.b, y.b, eps)) && all(equal(x.c, y.c, eps));
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool equal(Triangle<L, T, Q> const &x, Triangle<L, T, Q> const &y, int MaxULPs) {
  return all(equal(x.a, y.a, MaxULPs)) && all(equal(x.b, y.b, MaxULPs)) && all(equal(x.c, y.c, MaxULPs));
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool equal(Triangle<L, T, Q> const &x, Triangle<L, T, Q> const &y, vec<L, int, Q> const &MaxULPs) {
  return all(equal(x.a, y.a, MaxULPs)) && all(equal(x.b, y.b, MaxULPs)) && all(equal(x.c, y.c, MaxULPs));
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool notEqual(Triangle<L, T, Q> const &x, Triangle<L, T, Q> const &y, T eps = epsilon<T>()) {
  return any(notEqual(x.a, y.a, eps)) || any(notEqual(x.b, y.b, eps)) || any(notEqual(x.c, y.c, eps));
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool notEqual(Triangle<L, T, Q> const &x, Triangle<L, T, Q> const &y, vec<L, T, Q> const &eps) {
  return any(notEqual(x.a, y.a, eps)) || any(notEqual(x.b, y.b, eps)) || any(notEqual(x.c, y.c, eps));
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool notEqual(Triangle<L, T, Q> const &x, Triangle<L, T, Q> const &y, int MaxULPs) {
  return any(notEqual(x.a, y.a, MaxULPs)) || any(notEqual(x.b, y.b, MaxULPs)) || any(notEqual(x.c, y.c, MaxULPs));
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool notEqual(Triangle<L, T, Q> const &x, Triangle<L, T, Q> const &y, vec<L, int, Q> const &MaxULPs) {
  return any(notEqual(x.a, y.a, MaxULPs)) || any(notEqual(x.b, y.b, MaxULPs)) || any(notEqual(x.c, y.c, MaxULPs));
}

/// Test if any component of the triangle is infinite
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool isinf(Triangle<L, T, Q> const &t) {
  return any(isinf(t.a)) || any(isinf(t.b)) || any(isinf(t.c));
}

/// Test if any component of the triangle is NaN
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool isnan(Triangle<L, T, Q> const &t) {
  return any(isnan(t.a)) || any(isnan(t.b)) || any(isnan(t.c));
}

/// Test if all components of the triangle are finite
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool isfinite(Triangle<L, T, Q> const &t) {
  return all(isfinite(t.a)) && all(isfinite(t.b)) && all(isfinite(t.c));
}

/// Test whether the triangle is degenerate: infinite or with zero surface area.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool isDegenerate(Triangle<L, T, Q> const &t, T eps = epsilon<T>()) {
  return all(equal(t.a, t.b, eps)) || all(equal(t.a, t.c, eps)) || all(equal(t.b, t.c, eps));
}

/// Return the centroid of the Triangle.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> centroid(Triangle<L, T, Q> const &t) {
  return (t.a + t.b + t.c) * (T(1) / T(3));
}

/// Return the surface area of the Triangle.
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER T area(Triangle<3, T, Q> const &t) {
  return T(0.5) * length(cross(t.b - t.a, t.c - t.a));
}

/// Return the surface area of the Triangle.
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER T area(Triangle<2, T, Q> const &t) {
  return (t.a.x - t.b.x) * (t.b.y - t.c.y) - (t.b.x - t.c.x) * (t.a.y - t.b.y);
}

/// Return the barycentric U-coordinate at a given point on the triangle.
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER T signedArea(Triangle<3, T, Q> const &t, vec<3, T, Q> const &pt) {
  return dot(cross(t.b - pt, t.c - pt), normalize(cross(t.b - t.a, t.c - t.a)));
}

/// Return the total edge length of the Triangle.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T perimeter(Triangle<L, T, Q> const &t) {
  return distance(t.a, t.b) + distance(t.b, t.c) + distance(t.c, t.a);
}

/// Return an edge of the Triangle.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER Segment<L, T, Q> edge(Triangle<L, T, Q> const &t, length_t index) {
  switch (index) {
    case 1: return Segment<L, T, Q>(t.b, t.c);
    case 2: return Segment<L, T, Q>(t.c, t.a);
    default: {
      return Segment<L, T, Q>(t.a, t.b);
    }
  }
}

/// Return a vertex of the triangle.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> vertex(Triangle<L, T, Q> const &t, length_t index) {
  switch (index) {
    case 1: return t.b;
    case 2: return t.c;
    default: {
      return t.a;
    }
  }
}

/// Return a vertex of the triangle.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> cornerPoint(Triangle<L, T, Q> const &t, length_t index) {
  return vertex(t, index);
}

#if 0
/// Return all face normals of the Triangle.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER void faceNormals(Triangle<L, T, Q> const &t, vec<L, T, Q> (&normals)[4]) {
  normals[0] = cross(t.b - t.a, t.c - t.a);
  normals[1] = cross(normals[0], t.b - t.a);
  normals[2] = cross(normals[0], t.c - t.a);
  normals[3] = cross(normals[0], t.c - t.b);
}

/// Return the directions of each edge.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER void faceNormals(Triangle<L, T, Q> const &t, vec<L, T, Q> (&directions)[3]) {
  directions[0] = normalize(t.b - t.a);
  directions[1] = normalize(t.c - t.a);
  directions[2] = normalize(t.c - t.b);
}
#endif

/// Project the Triangle onto the provided axis.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER void projectToAxis(Triangle<L, T, Q> const &t, vec<L, T, Q> const &axis, T &dMin, T &dMax) {
  const T db = dot(axis, t.b);
  const T dc = dot(axis, t.c);

  dMin = dMax = dot(axis, t.a);
  dMin = min(dMin, db);
  dMax = max(dMax, db);
  dMin = min(dMin, dc);
  dMax = max(dMax, dc);
}

/// Return an extreme point (furthest point) in the given direction.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> extremePoint(Triangle<L, T, Q> const &t, vec<L, T, Q> const &direction) {
  vec<L, T, Q> extremePt(T(0));
  T extremeDist = std::numeric_limits<T>::min();
  for (length_t i = 0; i < 3; ++i) {
    const vec<L, T, Q> &pt = vertex(t, i);
    const T d = dot(direction, pt);
    if (d > extremeDist) {
      extremeDist = d;
      extremePt = pt;
    }
  }
  return extremePt;
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> extremePoint(Triangle<3, T, Q> const &t, vec<3, T, Q> const &direction, T &distance) {
  const vec<L, T, Q> extremePt = extremePoint(t, direction);
  distance = dot(extremePt, direction);
  return extremePt;
}

/// Return the minimal AABB that encloses the Triangle
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER AABB<L, T, Q> boundingAABB(Triangle<L, T, Q> const &t) {
  return AABB<L, T, Q>(min(t.a, t.b, t.c), max(t.a, t.b, t.c));
}

/// @private
template<typename T>
GLM_GEOM_QUALIFIER T triangleArea2D(T x1, T y1, T x2, T y2, T x3, T y3) {
  return (x1 - x2) * (y2 - y3) - (x2 - x3) * (y1 - y2);
}

/// <summary>
/// Express the given point (a position on the plane formed by the triangle) in
/// terms of barycentric u, v, w coordinates.
///
/// To map to (u, v) coordinates use: (v, w).
/// </summary>
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<3, T, Q> barycentricUVW(Triangle<3, T, Q> const &t, vec<3, T, Q> const &point) {
  const vec<3, T, Q> m = cross(t.b - t.a, t.c - t.a);  // Unnormalized triangle normal.
  const vec<3, T, Q> m_abs = abs(m);

  T nu, nv, d;
  if (m_abs.x >= m_abs.y && m_abs.x >= m_abs.z) {  // YZ plane projection
    nu = triangleArea2D(point.y, point.z, t.b.y, t.b.z, t.c.y, t.c.z);  // P-B-C Area
    nv = triangleArea2D(point.y, point.z, t.c.y, t.c.z, t.a.y, t.a.z);  // P-C-A Area
    d = T(1) / m.x;
  }
  else if (m_abs.y >= m_abs.z) {  // XZ plane projection
    nu = triangleArea2D(point.x, point.z, t.b.x, t.b.z, t.c.x, t.c.z);
    nv = triangleArea2D(point.x, point.z, t.c.x, t.c.z, t.a.x, t.a.z);
    d = T(1) / -m.y;
  }
  else {  // XY plane projection
    nu = triangleArea2D(point.x, point.y, t.b.x, t.b.y, t.c.x, t.c.y);
    nv = triangleArea2D(point.x, point.y, t.c.x, t.c.y, t.a.x, t.a.y);
    d = T(1) / m.z;
  }

  const T u = nu * d;
  const T v = nv * d;
  return vec<3, T, Q>(u, v, T(1) - u - v);
}

/// <summary>
/// Express the given point in terms of barycentric u, v coordinates. To map to
/// (u, v, w) coordinates use: (1 - u - v, u, v).
/// </summary>
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<2, T, Q> barycentricUV(Triangle<3, T, Q> const &t, vec<3, T, Q> const &point) {
  const vec<3, T, Q> uvw = barycentricUVW(t, point);
  return vec<2, T, Q>(uvw.y, uvw.z);
}

/// Return the point at the given barycentric (UV) coordinates.
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<3, T, Q> barycentricPoint(Triangle<3, T, Q> const &t, T u, T v) {
  return t.a + ((t.b - t.a) * u + (t.c - t.a) * v);
}

/// barycentricPoint with vec2 uv coordinates
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<3, T, Q> barycentricPoint(Triangle<3, T, Q> const &t, vec<2, T, Q> const &uv) {
  return barycentricPoint(t, uv.x, uv.y);
}

/// Return the point at the given barycentric (UVW) coordinates.
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<3, T, Q> barycentricPoint(Triangle<3, T, Q> const &t, T u, T v, T w) {
  return u * t.a + v * t.b + w * t.c;
}

/// barycentricPoint with vec2 uvw coordinates
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<3, T, Q> barycentricPoint(Triangle<3, T, Q> const &t, vec<3, T, Q> const &uvw) {
  return barycentricPoint(t, uvw.x, uvw.y, uvw.z);
}

/// Return the Plane of the Triangle with counter-clockwise orientation.
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER Plane<3, T, Q> planeCCW(Triangle<3, T, Q> const &t) {
  return planeFrom(t.a, t.b, t.c);
}

/// Return an unnormalized counter-clockwise oriented normal for the given Triangle.
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<3, T, Q> unnormalizedNormalCCW(Triangle<3, T, Q> const &t) {
  return cross(t.b - t.a, t.c - t.a);
}

/// Return a counter-clockwise oriented normal for the given Triangle.
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<3, T, Q> normalCCW(Triangle<3, T, Q> const &t) {
  return normalize(unnormalizedNormalCCW(t));
}

/// Return the Plane of the Triangle with clockwise orientation.
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER Plane<3, T, Q> planeCW(Triangle<3, T, Q> const &t) {
  return planeFrom(t.a, t.c, t.b);
}

/// Return an unnormalized clockwise oriented normal for the given triangle.
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<3, T, Q> unnormalizedNormalCW(Triangle<3, T, Q> const &t) {
  return cross(t.c - t.a, t.b - t.a);
}

/// Return a clockwise oriented normal for the given triangle.
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<3, T, Q> normalCW(Triangle<3, T, Q> const &t) {
  return normalize(unnormalizedNormalCW(t));
}

/// Return true if the given point is contained within the Triangle.
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool contains(Triangle<3, T, Q> const &t, vec<3, T, Q> const &point, T thicknessSq = epsilon2<T>()) {
  const vec<3, T, Q> normal = cross(t.b - t.a, t.c - t.a);
  const T d = dot(normal, t.b - point);
  if (d * d <= thicknessSq * length2(normal)) {
    const vec<3, T, Q> br = barycentricUVW(t, point);
    return br.x >= -epsilon<T>() && br.y >= -epsilon<T>() && br.z >= -epsilon<T>();
  }
  return false;
}

/// Return true if the given Segment is fully contained within the Triangle.
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool contains(Triangle<3, T, Q> const &t, Segment<3, T, Q> const &segment, T thicknessSq = epsilon2<T>()) {
  return contains(t, segment.a, thicknessSq) && contains(t, segment.b, thicknessSq);
}

/// Return true if the given triangle is fully contained within the Triangle.
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool contains(Triangle<3, T, Q> const &t, Triangle<3, T, Q> const &other, T thicknessSq = epsilon2<T>()) {
  return contains(t, other.a, thicknessSq)
         && contains(t, other.b, thicknessSq)
         && contains(t, other.c, thicknessSq);
}

#if 0
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool contains(Triangle<3, T, Q> const &t, Disc<3, T, Q> const &disc, T thicknessSq = epsilon2<T>()) { }
#endif

/// <summary>
/// Line/Triangle intersection, includes the barycentric and parametric
/// intersection coordinates (infinity otherwise).
/// @private
/// </summary>
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER_OUTLINE T intersectTriangleLine(Triangle<L, T, Q> const &t, vec<L, T, Q> const &linePos, vec<L, T, Q> const &lineDir, T &u, T &v) {
  const T eps = intersect_eps<T>();
  const vec<L, T, Q> e1 = t.b - t.a;
  const vec<L, T, Q> e2 = t.c - t.a;
  const vec<L, T, Q> vt = linePos - t.a;
  const vec<L, T, Q> vp = cross(lineDir, e2);
  const vec<L, T, Q> vq = cross(vt, e1);

  const T det = dot(e1, vp);
  if (abs(det) <= eps) {  // determinant zero: line on plane of triangle
    u = v = std::numeric_limits<T>::infinity();
    return std::numeric_limits<T>::infinity();
  }

  const T invDet = T(1) / det;  // compute barycentric coordinates.
  u = dot(vt, vp) * invDet;
  v = dot(lineDir, vq) * invDet;
  if (u < -eps || u > (T(1) + eps))
    return std::numeric_limits<T>::infinity();
  if (v < -eps || (u + v) > (T(1) + eps))
    return std::numeric_limits<T>::infinity();

  return dot(e2, vq) * invDet;
}

/// <summary>
/// Return the closest point between the Triangle (along an edge) and provided
/// point. Also storing the closest Triangle point in the U,V,W barycentric.
/// </summary>
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER_OUTLINE vec<L, T, Q> closestPointTriangle(Triangle<L, T, Q> const &t, vec<L, T, Q> const &p, T &u, T &v, T &w) {
  const vec<L, T, Q> ba = t.b - t.a;
  const vec<L, T, Q> ca = t.c - t.a;
  const vec<L, T, Q> pa = p - t.a;
  const vec<L, T, Q> bp = p - t.b;
  const vec<L, T, Q> cp = p - t.c;

  const T d1 = dot(ba, pa);
  const T d2 = dot(ca, pa);
  const T d3 = dot(ba, bp);
  const T d4 = dot(ca, bp);
  const T d5 = dot(ba, cp);
  const T d6 = dot(ca, cp);

  const T vc = d1 * d4 - d3 * d2;
  const T vb = d5 * d2 - d1 * d6;
  const T va = d3 * d6 - d5 * d4;

  const T zero = T(0);
  if (d1 <= zero && d2 <= zero) {  // P is in the vertex region outside A.
    v = w = zero;
    u = T(1);
    return t.a;  // (1, 0, 0)
  }
  else if (d3 >= zero && d4 <= d3) {  // P is in the vertex region outside B.
    u = w = zero;
    v = T(1);
    return t.b;  // (0, 1, 0)
  }
  else if (vc <= zero && d1 >= zero && d3 <= zero) {  // P is in edge region of AB
    w = zero;
    v = d1 / (d1 - d3);
    u = T(1) - v;
    return t.a + v * ba;  // (1 - v, v, 0)
  }
  else if (d6 >= zero && d5 <= d6) {  // P is in the vertex region outside C.
    u = v = zero;
    w = T(1);
    return t.c;  // (0, 0, 1)
  }
  else if (vb <= zero && d2 >= zero && d6 <= zero) {  // P is in edge region of AC
    v = zero;
    w = d2 / (d2 - d6);
    u = T(1) - w;
    return t.a + v * ca;  // (1 - w, 0, w)
  }
  else if (va <= zero && d4 - d3 >= zero && d5 - d6 >= zero) {  // P is in edge region of BC
    u = zero;
    w = (d4 - d3) / (d4 - d3 + d5 - d6);
    v = T(1) - w;
    return t.b + v * (t.c - t.b);  // (0, 1 - w, w)
  }

  const T denom = T(1) / (va + vb + vc);  // P must be inside the face
  v = vb * denom;
  w = vc * denom;
  u = T(1) - v - w;  // va * denom
  return t.a + ba * v + ca * w;
}

/// <summary>
/// Compute the closest distance between a Triangle (along an edge) and a Line.
///
/// Returning a point on the Triangles edge 'closest' to the line, its
/// barycentric coordinates (outU, outV), and distance along the line (outD).
/// </summary>
template<length_t L, typename T, qualifier Q, typename Line>
GLM_GEOM_QUALIFIER_OUTLINE vec<L, T, Q> closestPointTriangleLine(Triangle<L, T, Q> const &t, Line const &line, T &outU, T &outV, T &outD) {
  T d1, d2, d3, d_unused;

  const vec<L, T, Q> pt1 = closestPoint(edge(t, 0), line, d_unused, d1);
  const vec<L, T, Q> pt2 = closestPoint(edge(t, 1), line, d_unused, d2);
  const vec<L, T, Q> pt3 = closestPoint(edge(t, 2), line, d_unused, d3);
  const T dist1 = distance2(pt1, getPoint(line, d1));
  const T dist2 = distance2(pt2, getPoint(line, d2));
  const T dist3 = distance2(pt3, getPoint(line, d3));

  T resultDist = d3;
  vec<L, T, Q> result = pt3;
  if (dist1 <= dist2 && dist1 <= dist3) {
    result = pt1;
    resultDist = d1;
  }
  else if (dist2 <= dist3) {
    result = pt2;
    resultDist = d2;
  }

  outD = resultDist;
  outU = barycentricUV(t, result).x;
  outV = barycentricUV(t, result).y;
  return result;
}

/// <summary>
/// Compute the closest distance between a Triangle (along an edge) and a Segment.
///
/// Returning a point on the Triangles edge 'closest' to the Line, its
/// barycentric coordinates (outU, outV), and distance along the line (outD).
/// </summary>
template<length_t L, typename T, qualifier Q, typename Line>
GLM_GEOM_QUALIFIER_OUTLINE vec<L, T, Q> closestPointTriangleSegment(Triangle<L, T, Q> const &t, Segment<L, T, Q> const &segment, T &outU, T &outV, T &outD) {
  outD = intersectTriangleLine(t, segment.a, segment.b - segment.a, outU, outV);
  if (outD >= T(0.0) && outD <= T(1.0)) {
    return barycentricPoint(t, outU, outV);
  }

  T d1(0);
  const vec<L, T, Q> pt1 = closestPointTriangleLine<L, T, Q, Segment<L, T, Q>>(t, segment, outU, outV, d1);
  const vec<L, T, Q> pt2 = closestPoint(t, segment.a);
  const vec<L, T, Q> pt3 = closestPoint(t, segment.b);
  const T l1 = distance2(pt1, getPoint(segment, d1));
  const T l2 = distance2(pt2, segment.a);
  const T l3 = distance2(pt3, segment.b);
  if (l1 <= l2 && l1 <= l3) {
    outD = d1;
    return pt1;
  }
  else if (l2 <= l3) {
    outD = T(0);
    return pt2;
  }
  else {
    outD = T(1);
    return pt3;
  }
}

/* Return the closest point on the triangle to the given object */

/// Return the closest point between the Triangle and given point, includes the
/// barycentric coordinate of the closest point on the triangle
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Triangle<L, T, Q> const &t, vec<L, T, Q> const &p, T &u, T &v, T &w) {
  return closestPointTriangle(t, p, u, v, w);
}

/// Return the closest point between the Triangle and given point
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Triangle<L, T, Q> const &t, vec<L, T, Q> const &p) {
  T u, v, w;
  return closestPointTriangle(t, p, u, v, w);
}

/// Return the closest point between the Triangle and given Segment, includes the
/// parametric intersection coordinate.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Triangle<L, T, Q> const &t, Segment<L, T, Q> const &segment, T &d) {
  T u, v;
  return closestPointTriangleSegment<L, T, Q, Segment<L, T, Q>>(t, segment, u, v, d);
}

/// Return the closest point between the Triangle and given Line, includes the
/// parametric intersection coordinate.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Triangle<L, T, Q> const &t, Line<L, T, Q> const &line, T &d) {
  T u, v;
  d = intersectTriangleLine(t, line.pos, line.dir, u, v);
  if (glm::isfinite(d))
    return barycentricPoint(t, u, v);
  return closestPointTriangleLine<L, T, Q, Line<L, T, Q>>(t, line, u, v, d);
}

#if 0
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Triangle<L, T, Q> const &t, AABB<L, T, Q> const &aabb) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Triangle<L, T, Q> const &t, Ray<L, T, Q> const &ray, T &d) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Triangle<L, T, Q> const &t, Plane<L, T, Q> const &plane) { }

template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<3, T, Q> closestPoint(Triangle<3, T, Q> const &t, Triangle<3, T, Q> const &other, vec<3, T, Q> &otherPt) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Triangle<L, T, Q> const &t, Sphere<L, T, Q> const &sphere) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Triangle<L, T, Q> const &t, Capsule<L, T, Q> const &capsule) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Triangle<L, T, Q> const &t, Disc<L, T, Q> const &disc) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPoint(Triangle<L, T, Q> const &t, Polygon<L, T, Q> const &polygon) { }

template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<3, T, Q> closestPoint(Triangle<3, T, Q> const &t, OBB<T, Q> const &obb) { }
#endif

/* Test whether the triangle and the given object intersect */

/// Triangle/Segment intersection test, includes the barycentric and parametric intersection coordinates.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Triangle<L, T, Q> const &t, Segment<L, T, Q> const &segment, T &u, T &v, T &d) {
  d = intersectTriangleLine(t, segment.a, segment.b - segment.a, u, v);
  return d >= T(0.0) && d <= T(1.0);
}

/// Triangle/Line intersection test, includes the barycentric and parametric intersection coordinates.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Triangle<L, T, Q> const &t, Line<L, T, Q> const &line, T &u, T &v, T &d) {
  d = intersectTriangleLine(t, line.pos, line.dir, u, v);
  return glm::isfinite(d);
}

/// Triangle/Ray intersection test, includes the barycentric and parametric intersection coordinates.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Triangle<L, T, Q> const &t, Ray<L, T, Q> const &ray, T &u, T &v, T &d) {
  d = intersectTriangleLine(t, ray.pos, ray.dir, u, v);
  return glm::isfinite(d) && d >= T(0);
}

/// Triangle/Plane intersection test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Triangle<L, T, Q> const &t, Plane<L, T, Q> const &plane) {
  return intersects(plane, t);
}

/// Triangle/Sphere intersection test, includes intersection coordinate.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Triangle<L, T, Q> const &t, Sphere<L, T, Q> const &sphere, vec<L, T, Q> &intersection, T eps = epsilon<T>()) {
  const T r = sphere.r + eps;
  intersection = closestPoint(t, sphere.pos);
  return distance2(intersection, sphere.pos) <= r * r;
}

/// Triangle/Segment intersection test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Triangle<L, T, Q> const &t, Segment<L, T, Q> const &segment) {
  T u, v;
  const T d = intersectTriangleLine(t, segment.a, segment.b - segment.a, u, v);
  return d >= T(0.0) && d <= T(1.0);
}

/// Triangle/Line intersection test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Triangle<L, T, Q> const &t, Line<L, T, Q> const &line) {
  T u, v;
  const T d = intersectTriangleLine(t, line.pos, line.dir, u, v);
  return glm::isfinite(d);
}

/// Triangle/Ray intersection test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Triangle<L, T, Q> const &t, Ray<L, T, Q> const &ray) {
  T u, v;
  const T d = intersectTriangleLine(t, ray.pos, ray.dir, u, v);
  return glm::isfinite(d) && d >= T(0);
}

/// Triangle/Sphere intersection test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Triangle<L, T, Q> const &t, Sphere<L, T, Q> const &sphere) {
  vec<L, T, Q> unused;
  return intersects(t, sphere, unused);
}

/// Triangle/Capsule intersection test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Triangle<L, T, Q> const &triangle, Capsule<L, T, Q> const &capsule) {
  return intersects(capsule, triangle);
}

/// Triangle/OBB intersection test
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Triangle<3, T, Q> const &triangle, OBB<T, Q> const &obb) {
  return intersects(obb, triangle);
}

#if 0
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Triangle<L, T, Q> const &t, Triangle<L, T, Q> const &other, Segment<L, T, Q> &intersection) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Triangle<L, T, Q> const &t, AABB<L, T, Q> const &aabb) { }
#endif

/* Return the distance between the triangle and the given object */

/// Return the distance between the Triangle and given point.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Triangle<L, T, Q> const &t, vec<L, T, Q> const &p) {
  return distance(closestPoint(t, p), p);
}

/// Return the distance between a Triangle and Sphere.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Triangle<L, T, Q> const &t, Sphere<L, T, Q> const &s) {
  return max(T(0), distance(t, s.pos) - s.r);
}

/// Return the distance between a Triangle and Capsule.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Triangle<L, T, Q> const &t, Capsule<L, T, Q> const &capsule) {
  T d(0);
  const vec<L, T, Q> point = closestPoint(t, capsule.l, d);
  return max(T(0), distance(point, getPoint(capsule.l, d)) - capsule.r);
}

#if 0
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Triangle<L, T, Q> const &t, AABB<L, T, Q> const &aabb) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Triangle<L, T, Q> const &t, Line<L, T, Q> const &line) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Triangle<L, T, Q> const &t, Segment<L, T, Q> const &segment) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Triangle<L, T, Q> const &t, Ray<L, T, Q> const &ray) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Triangle<L, T, Q> const &t, Disc<L, T, Q> const &disc) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Triangle<L, T, Q> const &t, Plane<L, T, Q> const &plane) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Triangle<L, T, Q> const &t, Triangle<L, T, Q> const &triangle) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Triangle<L, T, Q> const &t, Polygon<L, T, Q> const &polygon) { }

template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distance(Triangle<3, T, Q> const &t, OBB<T, Q> const &obb) { }
#endif

/// @}
GLM_GEOM_END_NAMESPACE
#if GLM_GEOM_TOSTRING
namespace glm {
  namespace detail {
    template<glm::length_t L, typename T, qualifier Q>
    struct compute_to_string<Triangle<L, T, Q>> {
      GLM_GEOM_QUALIFIER std::string call(Triangle<L, T, Q> const &t) {
        return detail::format("triangle(%s, %s, %s)",
          glm::to_string(t.a).c_str(),
          glm::to_string(t.b).c_str(),
          glm::to_string(t.c).c_str()
        );
      }
    };
  }
}
#endif
