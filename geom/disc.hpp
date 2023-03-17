/// @ref geom_disc
/// @file geom/disc.hpp
///
/// @defgroup geom_disc Disc
/// @ingroup geom
///

#pragma once

#include "setup.hpp"
#include "plane.hpp"

#include "../matrix_extensions.hpp"
#if GLM_MESSAGES == GLM_ENABLE && !defined(GLM_EXT_INCLUDED)
  #pragma message("GLM: GLM_EXT_GEOM_disc extension included")
#endif

GLM_GEOM_BEGIN_NAMESPACE
/// @addtogroup geom_disc
/// @{

/// <summary>
/// The intersection of a sphere |X - P|^2 == R^2 and plane dot(N, X - P) == 0,
/// where P is the discs center, R its radius, and N is its normal vector.
/// </summary>
template<length_t L, typename T, qualifier Q>
struct Disc {

  // -- Implementation detail --

  typedef T value_type;
  typedef Disc<L, T, Q> type;
  typedef vec<L, T, Q> point_type;

  // -- Data --

  point_type pos;  // Center point of this disc
  point_type normal;  // Direction of this disc
  value_type r;  // Radius of this disc

#if GLM_CONFIG_DEFAULTED_DEFAULT_CTOR == GLM_ENABLE
  Disc() GLM_DEFAULT_CTOR;
#else
  Disc()
  #if GLM_CONFIG_CTOR_INIT != GLM_CTOR_INIT_DISABLE
    : pos(T(0)), normal(T(0)), r(T(0))
  #endif
  {
  }
#endif

  Disc(T scalar)
    : pos(scalar), normal(scalar), r(scalar) {
  }

  Disc(point_type const &position, point_type const &normal_, value_type radius)
    : pos(position), normal(normal_), r(radius) {
  }

  Disc(Disc<L, T, Q> const &disc)
    : pos(disc.pos), normal(disc.normal), r(disc.r) {
  }

  Disc<L, T, Q> &operator=(Disc<L, T, Q> const &disc) {
    pos = disc.pos;
    normal = disc.normal;
    r = disc.r;
    return *this;
  }
};

template<length_t L, typename T, qualifier Q>
static Disc<L, T, Q> operator-(Disc<L, T, Q> const &disc) {
  return Disc<L, T, Q>(disc.pos, -disc.normal, disc.r);
}

template<length_t L, typename T, qualifier Q>
static bool operator==(Disc<L, T, Q> const &s1, Disc<L, T, Q> const &s2) {
  return s1.pos == s2.pos && s1.normal == s2.normal && detail::equal_to(s1.r, s2.r);
}

template<length_t L, typename T, qualifier Q>
static Disc<L, T, Q> operator+(Disc<L, T, Q> const &disc, vec<L, T, Q> const &offset) {
  return Disc<L, T, Q>(disc.pos + offset, disc.normal, disc.r);
}

template<length_t L, typename T, qualifier Q>
static Disc<L, T, Q> operator-(Disc<L, T, Q> const &disc, vec<L, T, Q> const &offset) {
  return Disc<L, T, Q>(disc.pos - offset, disc.normal, disc.r);
}

template<typename T, qualifier Q>
static Disc<3, T, Q> operator*(mat<3, 3, T, Q> const &m, Disc<3, T, Q> const &disc) {
  return Disc<3, T, Q>(m * disc.pos, normalize(m * disc.normal), length(m[0]) * disc.r);
}

template<typename T, qualifier Q>
static Disc<3, T, Q> operator*(mat<3, 4, T, Q> const &m, Disc<3, T, Q> const &disc) {
  return Disc<3, T, Q>(transformPos(m, disc.pos), normalize(transformDir(m, disc.normal)), length(m[0]) * disc.r);
}

template<typename T, qualifier Q>
static Disc<3, T, Q> operator*(mat<4, 3, T, Q> const &m, Disc<3, T, Q> disc) {
  const T scale = length(vec<3, T, Q>(m[0].x, m[0].y, m[0].z));
  return Disc<3, T, Q>(transformPos(m, disc.pos), transformDir(m, disc.normal), scale * disc.r);
}

template<typename T, qualifier Q>
static Disc<3, T, Q> operator*(mat<4, 4, T, Q> const &m, Disc<3, T, Q> disc) {
  const T scale = length(vec<3, T, Q>(m[0].x, m[0].y, m[0].z));
  return Disc<3, T, Q>(transformPos(m, disc.pos), transformDir(m, disc.normal), scale * disc.r);
}

template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER Disc<3, T, Q> operator*(qua<T, Q> const &q, Disc<3, T, Q> disc) {
  return Disc<3, T, Q>(q * disc.pos, q * disc.normal, disc.r);
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool equal(Disc<L, T, Q> const &x, Disc<L, T, Q> const &y, T eps = epsilon<T>()) {
  return all(equal(x.pos, y.pos, eps)) && all(equal(x.normal, y.normal, eps)) && glm::equal(x.r, y.r, eps);
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool equal(Disc<L, T, Q> const &x, Disc<L, T, Q> const &y, vec<L, T, Q> const &eps) {
  return all(equal(x.pos, y.pos, eps)) && all(equal(x.normal, y.normal, eps)) && glm::equal(x.r, y.r, eps[0]);
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool equal(Disc<L, T, Q> const &x, Disc<L, T, Q> const &y, int MaxULPs) {
  return all(equal(x.pos, y.pos, MaxULPs)) && all(equal(x.normal, y.normal, MaxULPs)) && glm::equal(x.r, y.r, MaxULPs);
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool equal(Disc<L, T, Q> const &x, Disc<L, T, Q> const &y, vec<L, int, Q> const &MaxULPs) {
  return all(equal(x.pos, y.pos, MaxULPs)) && all(equal(x.normal, y.normal, MaxULPs)) && glm::equal(x.r, y.r, MaxULPs[0]);
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool notEqual(Disc<L, T, Q> const &x, Disc<L, T, Q> const &y, T eps = epsilon<T>()) {
  return any(notEqual(x.pos, y.pos, eps)) || any(notEqual(x.normal, y.normal, eps)) || glm::notEqual(x.r, y.r, eps);
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool notEqual(Disc<L, T, Q> const &x, Disc<L, T, Q> const &y, vec<L, T, Q> const &eps) {
  return any(notEqual(x.pos, y.pos, eps)) || any(notEqual(x.normal, y.normal, eps)) || glm::notEqual(x.r, y.r, eps[0]);
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool notEqual(Disc<L, T, Q> const &x, Disc<L, T, Q> const &y, int MaxULPs) {
  return any(notEqual(x.pos, y.pos, MaxULPs)) || any(notEqual(x.normal, y.normal, MaxULPs)) || glm::notEqual(x.r, y.r, MaxULPs);
}

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER GLM_CONSTEXPR bool notEqual(Disc<L, T, Q> const &x, Disc<L, T, Q> const &y, vec<L, int, Q> const &MaxULPs) {
  return any(notEqual(x.pos, y.pos, MaxULPs)) || any(notEqual(x.normal, y.normal, MaxULPs)) || glm::notEqual(x.r, y.r, MaxULPs[0]);
}

/// Test if any component of the disc is infinite.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool isinf(Disc<L, T, Q> const &disc) {
  return any(isinf(disc.pos)) || any(isinf(disc.normal)) || glm::isinf(disc.r);
}

/// Test if any component of the disc is NaN
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool isnan(Disc<L, T, Q> const &disc) {
  return any(isnan(disc.pos)) || any(isnan(disc.normal)) || glm::isnan(disc.r);
}

/// Test if all components of the disc are finite.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool isfinite(Disc<L, T, Q> const &disc) {
  return all(isfinite(disc.pos)) && all(isfinite(disc.normal)) && glm::isfinite(disc.r);
}

/// Test whether the disc is degenerate, i.e., not finite or if the normal is null
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool isDegenerate(Disc<L, T, Q> const &disc) {
  return !all(isfinite(disc.normal)) || isNull(disc.normal, epsilon<T>()) || !glm::isfinite(disc.r);
}

/// Return the U-vector (i.e., local space x axis) of the Disc
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<3, T, Q> basisU(Disc<3, T, Q> const &disc) {
  return perpendicular(disc.normal);
}

/// Return the V-vector (i.e., local-space y axis) of the Disc
template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<3, T, Q> basisV(Disc<3, T, Q> const &disc) {
  return perpendicular2(disc.normal);
}

/// <summary>
/// Return a point along this Disc, relative to a normalized offset, where 'd'
/// is a value between zero (center of Disc) and one (edge of Disc).
/// </summary>
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<3, T, Q> getPoint(Disc<L, T, Q> const &disc, T radians, T d = T(1)) {
  T sin, cos;
  sincos(radians, sin, cos);
  return disc.pos + disc.r * d * (sin * basisU(disc) + cos * basisV(disc));
}

/// Return an extreme point along the Disc edge (Circle), i.e., the furthest point in a given direction
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> extremePoint(Disc<L, T, Q> const &disc, vec<L, T, Q> const &dir) {
  const vec<L, T, Q> d = projNorm(disc.normal, dir);
  return isNull(d, T(0)) ? disc.pos : disc.pos + scaleLength(d, disc.r);
}

/// Return the Plane that defines the Disc, i.e., all points of the disc lie on the returned Plane.
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER Plane<L, T, Q> containingPlane(Disc<L, T, Q> const &disc) {
  return Plane<L, T, Q>(disc.pos, disc.normal);
}

/// Return the point on the Circle that is closest the given position
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPointToEdge(Disc<L, T, Q> const &disc, vec<L, T, Q> const &point) {
  const vec<L, T, Q> pointOnPlane = project(containingPlane(disc), point);
  const vec<L, T, Q> diff = pointOnPlane - disc.pos;
  return isNull(diff, T(0)) ? getPoint(disc, T(0)) : disc.pos + scaleLength(diff, disc.r);
}

/// Return the point on the Disc that is closest the given position
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPointToDisc(Disc<L, T, Q> const &disc, vec<L, T, Q> const &point) {
  const vec<L, T, Q> pointOnPlane = project(containingPlane(disc), point);
  vec<L, T, Q> diff = pointOnPlane - disc.pos;
  const T dist = length2(diff);
  if (dist > disc.r * disc.r)
    diff *= (disc.r / sqrt(dist));
  return disc.pos + diff;
}

/// Return the distance between the Circle and given position
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distanceToEdge(Disc<L, T, Q> const &disc, vec<L, T, Q> const &point) {
  return distance(closestPointToEdge(disc, point), point);
}

/// Return the distance between the Disc and given position
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distanceToDisc(Disc<L, T, Q> const &disc, vec<L, T, Q> const &point) {
  return distance(closestPointToDisc(disc, point), point);
}

/// Return true if the given point lies on the circle
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool edgeContains(Disc<L, T, Q> const &disc, vec<L, T, Q> const &point, T maxDistance) {
  return distanceToEdge(disc, point) <= maxDistance;
}

#if 0
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPointToEdge(Disc<L, T, Q> const &disc, Line<L, T, Q> const &line, T &d) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPointToEdge(Disc<L, T, Q> const &disc, Segment<L, T, Q> const &segment, T &d) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPointToEdge(Disc<L, T, Q> const &disc, Ray<L, T, Q> const &ray, T &d) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPointToDisc(Disc<L, T, Q> const &disc, Line<L, T, Q> const &line, T &d) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPointToDisc(Disc<L, T, Q> const &disc, Segment<L, T, Q> const &segment, T &d) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER vec<L, T, Q> closestPointToDisc(Disc<L, T, Q> const &disc, Ray<L, T, Q> const &ray, T &d) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distanceToEdge(Disc<L, T, Q> const &disc, Line<L, T, Q> const &line, T &d) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distanceToEdge(Disc<L, T, Q> const &disc, Segment<L, T, Q> const &segment, T &d) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distanceToEdge(Disc<L, T, Q> const &disc, Ray<L, T, Q> const &ray, T &d) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distanceToDisc(Disc<L, T, Q> const &disc, Line<L, T, Q> const &line, T &d) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distanceToDisc(Disc<L, T, Q> const &disc, Segment<L, T, Q> const &segment, T &d) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER T distanceToDisc(Disc<L, T, Q> const &disc, Ray<L, T, Q> const &ray, T &d) { }

template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool discContains(Disc<L, T, Q> const &disc, vec<L, T, Q> const &point) { }
#endif

/// Disc/Plane containment test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool contains(Disc<L, T, Q> const &disc, Plane<L, T, Q> const &plane, T eps = epsilon<T>()) {
  const bool same_dir = disc.r <= eps || glm::equal(abs(dot(plane.normal, disc.normal)), T(1), eps);
  return contains(plane, disc.pos, eps) && same_dir;
}

/// Disc/Line intersection test, includes parametric line coordinate
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Disc<L, T, Q> const &disc, Line<L, T, Q> const &line, T &d, T eps = epsilon<T>()) {
  const T r = disc.r + eps;
  return !intersects(line, containingPlane(disc), d)
         && distance2(disc.pos, getPoint(line, d)) <= (r * r);
}

/// Disc/Ray intersection test, includes parametric line coordinate
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Disc<L, T, Q> const &disc, Ray<L, T, Q> const &ray, T &d, T eps = epsilon<T>()) {
  const T r = disc.r + eps;
  return !intersects(ray, containingPlane(disc), d)
         && distance2(disc.pos, getPoint(ray, d)) <= (r * r);
}

/// Disc/Segment intersection test, includes parametric line coordinate
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Disc<L, T, Q> const &disc, Segment<L, T, Q> const &segment, T &d, T eps = epsilon<T>()) {
  const T r = disc.r + eps;
  return !intersects(segment, containingPlane(disc), d)
         && distance2(disc.pos, getPoint(segment, d)) <= (r * r);
}

/// Disc/Plane intersection test, includes points of intersection
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER int intersects(Disc<L, T, Q> const &disc, Plane<L, T, Q> const &plane, vec<L, T, Q> &pt1, vec<L, T, Q> &pt2, T eps = epsilon<T>()) {
  Line<L, T, Q> line(T(0));
  pt1 = pt2 = vec<L, T, Q>(T(0));
  if (intersects(plane, containingPlane(disc), line)) {
    line.pos -= disc.pos;  // Shift to origin
    const T a = T(1);
    const T b = T(2) * dot(line.pos, line.dir);
    const T c = length2(line.pos) - disc.r * disc.r;
    const T discrim = b * b - T(4) * a * c;
    if (discrim < -eps)
      return 0;
    else if (discrim < eps) {
      const T denom = T(1) / (T(2) * a);
      pt1 = pt2 = disc.pos + getPoint(line, -b * denom);
      return 1;
    }
    else {
      const T denom = T(1) / (T(2) * a);
      const T sqrtdisc = sqrt(discrim);
      pt1 = disc.pos + getPoint(line, (-b + sqrtdisc) * denom);
      pt2 = disc.pos + getPoint(line, (-b - sqrtdisc) * denom);
      return 2;
    }
  }
  return 0;
}

/// Disc/Line intersection test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Disc<L, T, Q> const &disc, Line<L, T, Q> const &line) {
  T d;
  return intersects(disc, line, d);
}

/// Disc/Ray intersection test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Disc<L, T, Q> const &disc, Ray<L, T, Q> const &ray) {
  T d;
  return intersects(disc, ray, d);
}

/// Disc/Segment intersection test
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER bool intersects(Disc<L, T, Q> const &disc, Segment<L, T, Q> const &segment) {
  T d;
  return intersects(disc, segment, d);
}

#if 0
template<length_t L, typename T, qualifier Q>
GLM_GEOM_QUALIFIER int intersects(Disc<L, T, Q> const &disc, AABB<L, T, Q> const &aabb) { }

template<typename T, qualifier Q>
GLM_GEOM_QUALIFIER int intersects(Disc<3, T, Q> const &disc, OBB<T, Q> const &obb) { }
#endif

/// @}
GLM_GEOM_END_NAMESPACE
#if GLM_GEOM_TOSTRING
namespace glm {
  namespace detail {
    template<glm::length_t L, typename T, qualifier Q>
    struct compute_to_string<Disc<L, T, Q>> {
      GLM_GEOM_QUALIFIER std::string call(Disc<L, T, Q> const &disc) {
        char const *LiteralStr = literal<T, std::numeric_limits<T>::is_iec559>::value();
        std::string FormatStr(detail::format("disc(%%s, %%s, %s)", LiteralStr));
        return detail::format(FormatStr.c_str(),
          glm::to_string(disc.pos).c_str(),
          glm::to_string(disc.normal).c_str(),
          static_cast<typename cast<T>::value_type>(disc.r)
        );
      }
    };
  }
}
#endif
