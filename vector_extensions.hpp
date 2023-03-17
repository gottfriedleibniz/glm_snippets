/// @ref gtx_vector_extensions
/// @file vector_extensions.hpp
///
/// @defgroup gtx_vector_extensions GLM_GTX_vector_extensions
/// @ingroup gtx
///
/// GLM vector extensions
///  1. API unifying functions, usually handling functions without genType or vec<1, genType> declarations
///  2. Support for C99/C++11 math functions
///  3. Functions emulated/ported from other popular vector-math libraries
///
/// @see core_func_common
/// @see ext_vector_common

#pragma once
#if !defined(GLM_ENABLE_EXPERIMENTAL)
  #define GLM_ENABLE_EXPERIMENTAL
#endif

#include <limits>
#include <type_traits>

#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>
#include <glm/gtc/random.hpp>
#include <glm/gtx/compatibility.hpp>
#include <glm/gtx/extended_min_max.hpp>
#include <glm/gtx/orthonormalize.hpp>
#include <glm/gtx/projection.hpp>
#include <glm/gtx/vector_angle.hpp>
#include <glm/gtx/vector_query.hpp>
#include <glm/ext/scalar_common.hpp>
#include <glm/ext/scalar_constants.hpp>
#include <glm/ext/quaternion_trigonometric.hpp>

#include "scalar_extensions.hpp"
#if GLM_MESSAGES == GLM_ENABLE && !defined(GLM_EXT_INCLUDED)
  #pragma message("GLM: GLM_GTX_vector_ext extension included")
#endif

namespace glm {

  /// Controlled by compilation flags:
  ///   - GLM_FORCE_LEFT_HANDED
  ///   - GLM_FORCE_Z_UP
  ///
  /// @addtogroup gtx_handed_coordinate_space
  /// @{

  /// Return the direction vector representing the "right" to an actor in world space.
  template<typename T, qualifier Q = defaultp>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR vec<3, T, Q> right() {
    return vec<3, T, Q>(T(1), T(0), T(0));
  }

  /// Return the direction vector representing "up" to an actor in world space.
  template<typename T, qualifier Q = defaultp>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR vec<3, T, Q> up() {
    return vec<3, T, Q>(
#if defined(GLM_FORCE_Z_UP)
      T(0), T(0), T(1)
#else
      T(0), T(1), T(0)
#endif
    );
  }

  /// Return the direction vector representing "forward" to an actor in left-handed world space.
  template<typename T, qualifier Q = defaultp>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR vec<3, T, Q> forwardLH() {
    return vec<3, T, Q>(
#if defined(GLM_FORCE_Z_UP)
      T(0), T(-1), T(0)
#else
      T(0), T(0), T(1)
#endif
    );
  }

  /// Return the direction vector representing "forward" to an actor in right-handed world space.
  template<typename T, qualifier Q = defaultp>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR vec<3, T, Q> forwardRH() {
    return vec<3, T, Q>(
#if defined(GLM_FORCE_Z_UP)
      T(0), T(1), T(0)
#else
      T(0), T(0), T(-1)
#endif
    );
  }

  /// Return the direction vector representing "forward" to an actor in world space.
  template<typename T, qualifier Q = defaultp>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR vec<3, T, Q> forward() {
#if defined(GLM_FORCE_LEFT_HANDED)
    return forwardLH<T, Q>();
#else
    return forwardRH<T, Q>();
#endif
  }

  /// @}

  /// @addtogroup gtx_vector_extensions
  /// @{

  /// <summary>
  /// all(equal(x, y)) shorthand
  /// @private @see core_func_vector_relational
  /// </summary>
  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR bool all_equal(vec<L, T, Q> const &x, vec<L, T, Q> const &y) {
    return all(equal(x, y));
  }

  /// <summary>
  /// all(equal(x, y, eps)) shorthand
  /// @private @see core_func_vector_relational
  /// </summary>
  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR bool all_equal(vec<L, T, Q> const &x, vec<L, T, Q> const &y, T eps) {
    return all(equal(x, y, eps));
  }

  /// <summary>
  /// all(equal(x, y, MaxULPs)) shorthand
  /// @private @see core_func_vector_relational
  /// </summary>
  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR bool all_equal(vec<L, T, Q> const &x, vec<L, T, Q> const &y, int MaxULPs) {
    return all(equal(x, y, MaxULPs));
  }

  /// <summary>
  /// all(equal(x, y, eps)) shorthand
  /// @private @see core_func_vector_relational
  /// </summary>
  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR bool all_equal(vec<L, T, Q> const &x, vec<L, T, Q> const &y, vec<L, T, Q> const &eps) {
    return all(equal(x, y, eps));
  }

  /// <summary>
  /// all(equal(x, y, MaxULPs)) shorthand
  /// @private @see core_func_vector_relational
  /// </summary>
  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR bool all_equal(vec<L, T, Q> const &x, vec<L, T, Q> const &y, vec<L, int, Q> const &MaxULPs) {
    return all(equal(x, y, MaxULPs));
  }

  /// <summary>
  /// any(notEqual(x, y)) shorthand
  /// @private @see core_func_vector_relational
  /// </summary>
  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR bool any_notEqual(vec<L, T, Q> const &x, vec<L, T, Q> const &y) {
    return any(notEqual(x, y));
  }

  /// <summary>
  /// any(notEqual(x, y, eps)) shorthand
  /// @private @see core_func_vector_relational
  /// </summary>
  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR bool any_notEqual(vec<L, T, Q> const &x, vec<L, T, Q> const &y, T eps) {
    return any(notEqual(x, y, eps));
  }

  /// <summary>
  /// any(notEqual(x, y, MaxULPs)) shorthand
  /// @private @see core_func_vector_relational
  /// </summary>
  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR bool any_notEqual(vec<L, T, Q> const &x, vec<L, T, Q> const &y, int MaxULPs) {
    return any(notEqual(x, y, MaxULPs));
  }

  /// <summary>
  /// any(notEqual(x, y, eps)) shorthand
  /// @private @see core_func_vector_relational
  /// </summary>
  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR bool any_notEqual(vec<L, T, Q> const &x, vec<L, T, Q> const &y, vec<L, T, Q> const &eps) {
    return any(notEqual(x, y, eps));
  }

  /// <summary>
  /// any(notEqual(x, y, MaxULPs)) shorthand
  /// @private @see core_func_vector_relational
  /// </summary>
  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR bool any_notEqual(vec<L, T, Q> const &x, vec<L, T, Q> const &y, vec<L, int, Q> const &MaxULPs) {
    return any(notEqual(x, y, MaxULPs));
  }

  /// <summary>
  /// all(lessThan(x, y)) shorthand
  /// @private @see core_func_vector_relational
  /// </summary>
  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR bool all_lessThan(vec<L, T, Q> const &x, vec<L, T, Q> const &y) {
    return all(lessThan(x, y));
  }

  /// <summary>
  /// all(lessThanEqual(x, y)) shorthand
  /// @private @see core_func_vector_relational
  /// </summary>
  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR bool all_lessThanEqual(vec<L, T, Q> const &x, vec<L, T, Q> const &y) {
    return all(lessThanEqual(x, y));
  }

  /// <summary>
  /// all(greaterThan(x, y)) shorthand
  /// @private @see core_func_vector_relational
  /// </summary>
  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR bool all_greaterThan(vec<L, T, Q> const &x, vec<L, T, Q> const &y) {
    return all(greaterThan(x, y));
  }

  /// <summary>
  /// all(greaterThanEqual(x, y)) shorthand
  /// @private @see core_func_vector_relational
  /// </summary>
  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR bool all_greaterThanEqual(vec<L, T, Q> const &x, vec<L, T, Q> const &y) {
    return all(greaterThanEqual(x, y));
  }

  /// <summary>
  /// any(isinf(x)) shorthand
  /// @private @see core_func_common
  /// </summary>
  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR bool any_isinf(vec<L, T, Q> const &x) {
    return any(isinf(x));
  }

  /// <summary>
  /// all(isfinite(x)) shorthand
  /// @private @see core_func_common
  /// </summary>
  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR bool all_isfinite(vec<L, T, Q> const &x) {
    return all(isfinite(x));
  }

  /// <summary>
  /// glm::any(glm::isnan(...)) shorthand
  /// @private @see core_func_common
  /// </summary>
  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR bool any_isnan(vec<L, T, Q> const &x) {
    return any(isnan(x));
  }

  /// <summary>
  /// Returns 1.0 if >= 0, or –1.0 if x < 0
  /// @see core_func_common
  /// </summary>
  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<L, T, Q> signP(vec<L, T, Q> const &x) {
    return vec<L, T, Q>(lessThanEqual(vec<L, T, Q>(0), x)) - vec<L, T, Q>(lessThan(x, vec<L, T, Q>(0)));
  }

  /// <summary>
  /// Returns 1.0 if x > 0, or –1.0 if x <= 0
  /// @see core_func_common
  /// </summary>
  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<L, T, Q> signN(vec<L, T, Q> const &x) {
    return vec<L, T, Q>(lessThan(vec<L, T, Q>(0), x)) - vec<L, T, Q>(lessThanEqual(x, vec<L, T, Q>(0)));
  }

  /// <summary>
  /// ceilMultiple that accepts a scalar multiple parameter.
  ///
  /// @tparam T Floating-point or integer scalar types
  /// @param Source Source values to which is applied the function
  /// @param Multiple Must be a null or positive value
  /// @see gtc_round
  /// </summary>
  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<L, T, Q> ceilMultiple(vec<L, T, Q> const &Source, T Multiple) {
    return ceilMultiple(Source, vec<L, T, Q>(Multiple));
  }

  /// <summary>
  /// floorMultiple that accepts a scalar multiple parameter.
  ///
  /// @tparam T Floating-point or integer scalar types
  /// @param Source Source values to which is applied the function
  /// @param Multiple Must be a null or positive value
  /// @see gtc_round
  /// </summary>
  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<L, T, Q> floorMultiple(vec<L, T, Q> const &Source, T Multiple) {
    return floorMultiple(Source, vec<L, T, Q>(Multiple));
  }

  /// <summary>
  /// roundMultiple that accepts a scalar multiple parameter.
  ///
  /// @tparam T Floating-point or integer scalar types
  /// @param Source Source values to which is applied the function
  /// @param Multiple Must be a null or positive value
  /// @see gtc_round
  /// </summary>
  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<L, T, Q> roundMultiple(vec<L, T, Q> const &Source, T Multiple) {
    return roundMultiple(Source, vec<L, T, Q>(Multiple));
  }

  /* Numeric extensions */

  /// <summary>
  /// Return true if all vector elements are equal (or within eps of each other).
  /// @see gtx_vector_query
  /// </summary>
  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER bool isUniform(vec<L, T, Q> const &v, T eps = epsilon<T>()) {
    bool result = true;
    for (length_t i = 1; i < L; ++i)
      result &= equal(v[i], v[0], eps);
    return result;
  }

  /// <summary>
  /// Reverse the elements of a vector
  /// @see gtx_vector_query
  /// </summary>
  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<L, T, Q> reverse(vec<L, T, Q> const &v) {
    vec<L, T, Q> result;
    for (length_t i = 0; i < L; ++i)
      result[i] = v[L - i - 1];
    return result;
  }

  /// <summary>
  /// Calculate sin and cos simultaneously.
  /// @param[in] v vector representing angles in radians
  /// @param[out] s the sine of value in the range [-1 ; +1]
  /// @param[out] c the cosine of value in the range [-1 ; +1]
  /// @see core_func_trigonometric
  /// </summary>
  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER void sincos(vec<L, T, Q> const &v, vec<L, T, Q> &s, vec<L, T, Q> &c) {
    s = sin(v);
    c = cos(v);
  }

  /// <summary>
  /// Cardinal sine
  /// @see core_func_trigonometric
  /// </summary>
  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<L, T, Q> sinc(vec<L, T, Q> const &v) {
    return detail::functor1<vec, L, T, T, Q>::call(sinc, v);
  }

  /// <summary>
  /// Normalized cardinal sine
  /// @see core_func_trigonometric
  /// </summary>
  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<L, T, Q> sincn(vec<L, T, Q> const &v) {
    return detail::functor1<vec, L, T, T, Q>::call(sincn, v);
  }

  /// <summary>
  /// Logistic function with basic overflow handling; underflow to-be-determined.
  /// @see core_func_exponential
  /// </summary>
  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<L, T, Q> logistic(vec<L, T, Q> const &v) {
    GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'logistic' only accept floating-point inputs.");
    return detail::functor1<vec, L, T, T, Q>::call(logistic, v);
  }

  /// <summary>
  /// Create a normalized vector2 from an angle (in radians).
  /// @see core_func_trigonometric
  /// </summary>
  template<typename T, qualifier Q = glm::defaultp>
  GLM_FUNC_QUALIFIER vec<2, T, Q> fromAngle(T angle) {
    T sin, cos;
    sincos(angle, sin, cos);
    return vec<2, T, Q>(sin, cos);
  }

  /// <summary>
  /// Return the direction vector given spherical coordinates.
  /// @see core_func_trigonometric
  /// </summary>
  template<typename T, qualifier Q = defaultp>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR vec<3, T, Q> spherical(T phi, T theta) {
    T sinphi, cosphi, sintheta, costheta;  // see sphericalRand
    sincos(phi, sinphi, cosphi);
    sincos(theta, sintheta, costheta);
    return vec<3, T, Q>(sinphi * costheta, sinphi * sintheta, cosphi);
  }

  /// <summary>
  /// Return a copy of the vector 'v' with its length clamped to 'maxLength'
  /// @see core_func_geometric
  /// </summary>
  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<L, T, Q> clampLength(vec<L, T, Q> const &v, T maxLength) {
    return (length2(v) > (maxLength * maxLength)) ? (normalize(v) * maxLength) : v;
  }

  template<typename genType>
  GLM_FUNC_QUALIFIER genType clampLength(genType x, genType maxLength) {
    return clampLength(vec<1, genType>(x), maxLength).x;
  }

  /// <summary>
  /// Scales the length of vector "v" to "newLength".
  /// @see core_func_geometric
  /// </summary>
  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<L, T, Q> scaleLength(vec<L, T, Q> const &v, T newLength) {
    GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'scaleLength' only accept floating-point inputs");
    const T sqlen = length2(v);
    if (sqlen < epsilon<T>()) {
      vec<L, T, Q> result(T(0));
      result[0] = newLength;
      return result;
    }
    return v * (newLength / sqrt(sqlen));
  }

  template<typename genType>
  GLM_FUNC_QUALIFIER genType scaleLength(genType x, genType newLength) {
    return scaleLength(vec<1, genType>(x), newLength).x;
  }

  /// <summary>
  /// Returns the homogenized vector: divides all components by w
  /// @see core_func_geometric
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<3, T, Q> homogenize(vec<4, T, Q> const &v) {
    return vec<3, T, Q>(v.x / v.w, v.y / v.w, v.z / v.w);
  }

  /// <summary>
  /// Returns the cross product of v and {1,0,0}.
  /// @tparam T Floating-point scalar types.
  /// @see core_func_geometric
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<3, T, Q> crossXAxis(vec<3, T, Q> const &v) {
    return vec<3, T, Q>(T(0), v.z, -v.y);
  }

  /// <summary>
  /// Returns the cross product of v and {0,1,0}.
  /// @tparam T Floating-point scalar types.
  /// @see core_func_geometric
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<3, T, Q> crossYAxis(vec<3, T, Q> const &v) {
    return vec<3, T, Q>(-v.z, T(0), v.x);
  }

  /// <summary>
  /// Returns the cross product of v and {0,0,1}.
  /// @tparam T Floating-point scalar types.
  /// @see core_func_geometric
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<3, T, Q> crossZAxis(vec<3, T, Q> const &v) {
    return vec<3, T, Q>(v.y, -v.x, T(0));
  }

  /// <summary>
  /// Returns the cross product of {1,0,0} and v.
  /// @tparam T Floating-point scalar types.
  /// @see core_func_geometric
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<3, T, Q> xAxisCross(vec<3, T, Q> const &v) {
    return vec<3, T, Q>(T(0), -v.z, v.y);
  }

  /// <summary>
  /// Returns the cross product of {0,1,0} and v.
  /// @tparam T Floating-point scalar types.
  /// @see core_func_geometric
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<3, T, Q> yAxisCross(vec<3, T, Q> const &v) {
    return vec<3, T, Q>(v.z, T(0), -v.x);
  }

  /// <summary>
  /// Returns the cross product of {0,0,1} and v.
  /// @tparam T Floating-point scalar types.
  /// @see core_func_geometric
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<3, T, Q> zAxisCross(vec<3, T, Q> const &v) {
    return vec<3, T, Q>(-v.y, v.x, T(0));
  }

  /// <summary>
  /// areOrthonormal that assumes the vectors are normalized
  /// @see gtx_vector_query
  /// @private
  /// </summary>
  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER bool areOrthonormal2(vec<L, T, Q> const &v0, vec<L, T, Q> const &v1, T const &epsilon) {
    // assert(isNormalized(v0, epsilon));
    // assert(isNormalized(v1, epsilon));
    return abs(dot(v0, v1)) <= epsilon;
  }

  template<typename genType>
  GLM_FUNC_QUALIFIER bool areOrthonormal2(genType v0, genType v1, genType eps = epsilon<genType>()) {
    return abs(dot(v0, v1)) <= eps;
  }

  /// <summary>
  /// Create a 'hint' axis for perpendicular/basis calculations.
  /// @see gtx_perpendicular
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<3, T, Q> hint(vec<3, T, Q> const &v) {
    return ((v.x * v.x) < T(0.5) * length2(v)) ? right<T, Q>() : forward<T, Q>();
  }

  /// <summary>
  /// Return a direction vector perpendicular to 'v'
  /// @param[in] v direction to compute the perpendicular of.
  /// @param[in] hint reference axis/vector to computer perpendicular.
  /// @param[in] hint2 Alternative reference axis if 'v' points towards 'hint'
  /// @see gtx_perpendicular
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<3, T, Q> perpendicular(vec<3, T, Q> const &v, vec<3, T, Q> const &hint = forward<T, Q>(), vec<3, T, Q> const &hint2 = up<T, Q>()) {
    const vec<3, T, Q> v2 = cross(v, hint);
    return detail::approx_zero(dot(v2, v2)) ? hint2 : normalize(v2);
  }

  /// <summary>
  /// Return a vector that is perpendicular to 'v' and the vector returned by perpendicular.
  /// @see gtx_perpendicular
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<3, T, Q> perpendicular2(vec<3, T, Q> const &v, vec<3, T, Q> const &hint = forward<T, Q>(), vec<3, T, Q> const &hint2 = up<T, Q>()) {
    return normalize(cross(v, perpendicular(v, hint, hint2)));
  }

  /// <summary>
  /// Update vectors 'out' and 'out2' to be orthogonal to 'v' and each other.
  /// @see gtx_perpendicular
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER void perpendicularBasis(vec<3, T, Q> const &v, vec<3, T, Q> &out, vec<3, T, Q> &out2) {
    GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'perpendicularBasis' only accept floating-point inputs");
    const T s = v.z >= T(0) ? T(1) : T(-1);
    const T a = T(-1) / (s + v.z);
    const T b = v.x * v.y * a;
    out = vec<3, T, Q>(T(1) + s * v.x * v.x * a, s * b, -s * v.x);
    out2 = vec<3, T, Q>(b, s + v.y * v.y * a, -v.y);
  }

  /// <summary>
  /// A mutable glm::orthonormalize implementation.
  /// @see gtx_orthonormalize
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER void orthonormalize2(vec<3, T, Q> &x, vec<3, T, Q> &y) {
    x = normalize(x);
    y = orthonormalize(y, x);
  }

  /// <summary>
  /// Normalize the provided vectors and make them orthogonal to each other.
  /// @see gtx_orthonormalize
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER void orthonormalize3(vec<3, T, Q> &x, vec<3, T, Q> &y, vec<3, T, Q> &z) {
    x = normalize(x);
    y = orthonormalize(y, x);
    const T dot0 = dot(x, z);
    const T dot1 = dot(y, z);
    z = normalize(z - (y * dot1 + x * dot0));
  }

  /// <summary>
  /// glm::proj with the assumption 'Normal' is already normalized.
  /// @see gtx_projection
  /// </summary>
  template<typename genType>
  GLM_FUNC_QUALIFIER genType projNorm(genType const &x, genType const &Normal) {
    return dot(x, Normal) * Normal;
  }

  /// <summary>
  /// Project a vector onto this plane defined by its normal orthogonal
  /// @see gtx_projection
  /// </summary>
  template<typename genType>
  GLM_FUNC_QUALIFIER genType projPlane(genType const &x, genType const &Normal) {
    return x - proj(x, Normal);
  }

  /// <summary>
  /// Breaks this vector down into parallel and perpendicular components with respect to the given direction
  /// @see gtx_projection
  /// </summary>
  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER void projDecompose(vec<L, T, Q> const &v, vec<L, T, Q> const &direction, vec<L, T, Q> &outPara, vec<L, T, Q> &outPerp) {
    outPara = proj(v, direction);
    outPerp = v - outPara;
  }

  template<typename genType>
  GLM_FUNC_QUALIFIER void projDecompose(genType v, genType direction, genType &outPara, genType &outPerp) {
    vec<1, genType> vParallel, vPerpendicular;
    projDecompose(vec<1, genType>(v), vec<1, genType>(direction), vParallel, vPerpendicular);
    outPara = vParallel.x;
    outPerp = vPerpendicular.x;
  }

  /// <summary>
  /// Return true if the three given points are collinear, i.e., lie on the same line.
  /// @see gtx_vector_query
  /// </summary>
  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER bool areCollinear(vec<L, T, Q> const &p1, vec<L, T, Q> const &p2, vec<L, T, Q> const &p3, T eps = epsilon<T>()) {
    return length2(cross(p2 - p1, p3 - p1)) <= eps;
  }

  /// <summary>
  /// Encode a normal using a spherical coordinate system
  /// @see gtc_packing
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<2, T, Q> sphericalEncode(vec<3, T, Q> const &v) {
    const vec<2, T, Q> Result(atan2(v.y, v.x) * one_over_pi<T>(), v.z);
    return Result * T(0.5) + T(0.5);
  }

  /// <summary>
  /// Decode a vector from a spherical coordinate system
  /// @see gtc_packing
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<3, T, Q> sphericalDecode(vec<2, T, Q> const &v) {
    const vec<2, T, Q> ang = v * T(2) - T(1);
    const vec<2, T, Q> sc(sin(ang.x * pi<T>()), cos(ang.x * pi<T>()));
    const vec<2, T, Q> phi(sqrt(T(1) - ang.y * ang.y), ang.y);
    return vec<3, T, Q>(sc.y * phi.x, sc.x * phi.x, phi.y);
  }

  /// <summary>
  /// Encode a normal using a octahedron coordinate system
  /// @see gtc_packing
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<2, T, Q> octahedronEncode(vec<3, T, Q> const &v) {
    const vec<3, T, Q> n = v / (abs(v.x) + abs(v.y) + abs(v.z));
    vec<2, T, Q> Result(T(0));
    if (n.z >= T(0)) {
      Result.x = n.x;
      Result.y = n.y;
    }
    else {
      Result.x = (T(1) - abs(n.y)) * (n.x >= T(0) ? T(1) : -T(1));
      Result.y = (T(1) - abs(n.x)) * (n.y >= T(0) ? T(1) : -T(1));
    }
    Result.x = Result.x * T(0.5) + T(0.5);
    Result.y = Result.y * T(0.5) + T(0.5);
    return Result;
  }

  /// <summary>
  /// Decode a vector from a octahedron coordinate system
  /// @see gtc_packing
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<3, T, Q> octahedronDecode(vec<2, T, Q> const &v) {
    const vec<2, T, Q> f(v.x * T(2) - T(1), v.y * T(2) - T(1));
    const vec<3, T, Q> n(f.x, f.y, T(1) - abs(f.x) - abs(f.y));
    const T t = saturate(-n.z);
    return normalize(n + vec<3, T, Q>(n.x >= T(0) ? -t : t, n.y >= T(0) ? -t : t, T(0)));
  }

  /// <summary>
  /// Return a reflection vector according to a incident and surface normal.
  ///
  /// @param[in] I Incident vector
  /// @param[in] N Surface normal
  /// @param[in] negativeSideRefractionIndex Refraction index of material exiting
  /// @param[in] positiveSideRefractionIndex Refraction index of material entering
  /// @see core_func_geometric
  /// </summary>
  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<L, T, Q> refract(vec<L, T, Q> const &I, vec<L, T, Q> const &N, T negativeSideRefractionIndex, T positiveSideRefractionIndex) {
    return refract(I, N, negativeSideRefractionIndex / positiveSideRefractionIndex);
  }

  /// <summary>
  /// Return a vector containing the Cartesian coordinates of a point specified
  /// in barycentric (relative to a N-Dimension triangle) coordinates.
  /// </summary>
  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<L, T, Q> barycentric(vec<L, T, Q> const &value1, vec<L, T, Q> const &value2, vec<L, T, Q> const &value3, T amount1, T amount2) {
    return (value1 + (amount1 * (value2 - value1))) + (amount2 * (value3 - value1));
  }

  template<typename genType>
  GLM_FUNC_QUALIFIER genType barycentric(genType value1, genType value2, genType value3, genType amount1, genType amount2) {
    return barycentric(vec<1, genType>(value1), vec<1, genType>(value2), vec<1, genType>(value3), amount1, amount2).x;
  }

  /// <summary>
  /// Wraps x between [0, maxValue]
  /// @see gtx_wrap
  /// </summary>
  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<L, T, Q> wrap(vec<L, T, Q> const &x, vec<L, T, Q> const &maxValue) {
    vec<L, T, Q> Result(T(0));
    for (length_t i = 0; i < L; ++i)
      Result[i] = wrap<T>(x[i], maxValue[i]);
    return Result;
  }

  /// <summary>
  /// Wraps x between [0, maxValue]
  /// @see gtx_wrap
  /// </summary>
  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<L, T, Q> wrap(vec<L, T, Q> const &x, T maxValue) {
    return wrap(x, vec<L, T, Q>(maxValue));
  }

  /// <summary>
  /// Loops 't' so that it is never greater than 'length' and less than zero.
  /// This function is an emulation of: Unity.Mathf.Repeat
  /// @see gtx_wrap
  /// </summary>
  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<L, T, Q> loopRepeat(vec<L, T, Q> const &t, vec<L, T, Q> const &length) {
    return detail::functor2<vec, L, T, Q>::call(loopRepeat, t, length);
  }

  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<L, T, Q> loopRepeat(vec<L, T, Q> const &t, T length) {
    if (detail::exactly_zero(length))
      return vec<L, T, Q>(T(0));
    return loopRepeat(t, vec<L, T, Q>(length));
  }

  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<L, T, Q> pingPong(vec<L, T, Q> const &v, vec<L, T, Q> const &length) {
    const vec<L, T, Q> t = loopRepeat(v, length * vec<L, T, Q>(2));
    return length - abs(t - length);
  }

  /// <summary>
  /// @see gtx_fast_trigonometry
  /// </summary>
  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<L, T, Q> wrapAngleSigned(vec<L, T, Q> const &x) {
    return detail::functor1<vec, L, T, T, Q>::call(wrapAngleSigned, x);
  }

  /// <summary>
  /// A lerp implementation that ensures values interpolate correctly when wrapped around two-pi.
  /// This function emulates Unity.Mathf.LerpAngle
  /// @see gtx_fast_trigonometry
  /// </summary>
  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<L, T, Q> lerpAngle(vec<L, T, Q> const &x, vec<L, T, Q> const &y, T t) {
    GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'lerpAngle' only accept floating-point inputs");
    vec<L, T, Q> Result(T(0));
    for (length_t i = 0; i < L; ++i)
      Result[i] = lerpAngle<T>(x[i], y[i], t);
    return Result;
  }

  /// <summary>
  /// A lerp implementation that ensures values interpolate correctly when wrapped around two-pi.
  /// @see gtx_fast_trigonometry
  /// </summary>
  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<L, T, Q> lerpAngle(vec<L, T, Q> const &x, vec<L, T, Q> const &y, vec<L, T, Q> const &t) {
    GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'lerpAngle' only accept floating-point inputs");
    vec<L, T, Q> Result(T(0));
    for (length_t i = 0; i < L; ++i)
      Result[i] = lerpAngle<T>(x[i], y[i], t[i]);
    return Result;
  }

  /// <summary>
  /// Return a position between two points, moving no further than maxDist.
  /// This function emulates Unity.Vector3.MoveTowards
  /// @see gtx_functions
  /// </summary>
  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<L, T, Q> moveTowards(vec<L, T, Q> const &current, vec<L, T, Q> const &target, T maxDist) {
    GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'moveTowards' only accept floating-point inputs");
    const vec<L, T, Q> delta = target - current;
    const T sqdist = dot(delta, delta);
    if (detail::approx_zero(sqdist) || (maxDist >= T(0) && sqdist <= maxDist * maxDist))
      return target;
    return current + (delta / (sqrt(sqdist) * maxDist));
  }

  /// <summary>
  /// Round a 'value' to a specified grid size.
  /// @see gtc_round
  /// @see gtx_scalar_extensions
  /// </summary>
  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<L, T, Q> snap(vec<L, T, Q> const &x, vec<L, T, Q> const &y) {
    return detail::functor2<vec, L, T, Q>::call(snap, x, y);
  }

  /// <summary>
  /// Round a 'value' to a specified grid size.
  /// @see gtc_round
  /// @see gtx_scalar_extensions
  /// </summary>
  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<L, T, Q> snap(vec<L, T, Q> const &x, T y) {
    return snap(x, vec<L, T, Q>(y));
  }

  /// <summary>
  /// Returns the normalized vector pointing from 'x' to 'y'.
  /// @see gtx_vector_common
  /// </summary>
  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR vec<L, T, Q> direction(vec<L, T, Q> const &x, vec<L, T, Q> const &y) {
    return normalize(y - x);
  }

  /// <summary>
  /// Returns a vector 't' such that lerp(x, y, t) == value (or 0 if x == y)
  /// @see gtx_compatibility
  /// </summary>
  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<L, T, Q> lerpInverse(vec<L, T, Q> const &x, vec<L, T, Q> const &y, vec<L, T, Q> const &value) {
    vec<L, T, Q> result(T(0));
    for (length_t i = 0; i < L; ++i)
      result[i] = lerpInverse(x[i], y[i], value[i]);
    return result;
  }

  /// <summary>
  /// Returns a vector 't' such that lerp(x, y, t) == value (or 0 if x == y)
  /// @see gtx_compatibility
  /// </summary>
  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<L, T, Q> lerpInverse(vec<L, T, Q> const &x, vec<L, T, Q> const &y, T value) {
    vec<L, T, Q> result(T(0));
    for (length_t i = 0; i < L; ++i)
      result[i] = lerpInverse(x[i], y[i], value);
    return result;
  }

  /// <summary>
  /// Normalize Lerp
  /// @see gtx_compatibility
  /// </summary>
  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<L, T, Q> nlerp(vec<L, T, Q> const &x, vec<L, T, Q> const &y, vec<L, T, Q> const &t) {
    return normalize(lerp(x, y, t));
  }

  /// <summary>
  /// Normalize Lerp
  /// @see gtx_compatibility
  /// </summary>
  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<L, T, Q> nlerp(vec<L, T, Q> const &x, vec<L, T, Q> const &y, T t) {
    return normalize(lerp(x, y, t));
  }

  /* Missing implicit genType support & API unification */

  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<L, T, Q> fclamp(vec<L, T, Q> const &x) {
    return fclamp(x, T(0), T(1));
  }

  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<1, T, Q> lerp(vec<1, T, Q> const &x, vec<1, T, Q> const &y, T a) {
    return mix(x, y, a);
  }

  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<1, T, Q> lerp(vec<1, T, Q> const &x, vec<1, T, Q> const &y, vec<1, T, Q> const &a) {
    return mix(x, y, a);
  }

  /* Functions with additional integral type support. */

  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER typename std::enable_if<std::is_integral<T>::value, vec<L, T, Q>>::type iceil(vec<L, T, Q> const &x) {
    return x;
  }

  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER typename std::enable_if<std::is_integral<T>::value, vec<L, T, Q>>::type ifloor(vec<L, T, Q> const &x) {
    return x;
  }

  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER typename std::enable_if<!std::is_integral<T>::value, vec<L, T, Q>>::type iceil(vec<L, T, Q> const &x) {
    return ceil(x);
  }

  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER typename std::enable_if<!std::is_integral<T>::value, vec<L, T, Q>>::type ifloor(vec<L, T, Q> const &x) {
    return floor(x);
  }

  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER typename std::enable_if<!std::is_integral<T>::value, vec<L, T, Q>>::type imod(vec<L, T, Q> const &x, T y) {
    return mod(x, y);
  }

  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER typename std::enable_if<!std::is_integral<T>::value, vec<L, T, Q>>::type imod(vec<L, T, Q> const &x, vec<L, T, Q> const &y) {
    return mod(x, y);
  }

  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<L, T, Q> pow(vec<L, T, Q> const &base, T exponent) {
    return pow(base, vec<L, T, Q>(exponent));
  }

  /* glm/gtx/associated.hpp extensions */

  /// @private
  namespace detail {
    template<typename T, bool Aligned>
    struct compute_associated {};

    template<length_t L, typename T, qualifier Q, bool Aligned>
    struct compute_associated<vec<L, T, Q>, Aligned> {
      static GLM_FUNC_QUALIFIER GLM_CONSTEXPR vec<L, T, Q> eq(vec<L, T, Q> const &cA, vec<L, T, Q> const &cB, vec<L, T, Q> const &vA, vec<L, T, Q> const &vB) {
        vec<L, T, Q> Result;
        for (length_t i = 0; i < L; ++i)
          Result[i] = equal_to(cA[i], cB[i]) ? vA[i] : vB[i];
        return Result;
      }

      static GLM_FUNC_QUALIFIER GLM_CONSTEXPR vec<L, T, Q> gt(vec<L, T, Q> const &cA, vec<L, T, Q> const &cB, vec<L, T, Q> const &vA, vec<L, T, Q> const &vB) {
        vec<L, T, Q> Result;
        for (length_t i = 0; i < L; ++i)
          Result[i] = (cA[i] > cB[i]) ? vA[i] : vB[i];
        return Result;
      }

      static GLM_FUNC_QUALIFIER GLM_CONSTEXPR vec<L, T, Q> gte(vec<L, T, Q> const &cA, vec<L, T, Q> const &cB, vec<L, T, Q> const &vA, vec<L, T, Q> const &vB) {
        vec<L, T, Q> Result;
        for (length_t i = 0; i < L; ++i)
          Result[i] = (cA[i] >= cB[i]) ? vA[i] : vB[i];
        return Result;
      }
    };
  }

  /// <summary>
  /// Equal comparison between 2 variables and returns 2 associated variable values
  /// @see gtx_associated_min_max
  /// </summary>
  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<L, T, Q> associatedEqual(vec<L, T, Q> const &cA, vec<L, T, Q> const &cB, vec<L, T, Q> const &vA, vec<L, T, Q> const &vB) {
    return detail::compute_associated<vec<L, T, Q>, detail::is_aligned<Q>::value>::eq(cA, cB, vA, vB);
  }

  /// <summary>
  /// Greater-than comparison between 2 variables and returns 2 associated variable values
  /// @see gtx_associated_min_max
  /// </summary>
  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<L, T, Q> associatedGreater(vec<L, T, Q> const &cA, vec<L, T, Q> const &cB, vec<L, T, Q> const &vA, vec<L, T, Q> const &vB) {
    return detail::compute_associated<vec<L, T, Q>, detail::is_aligned<Q>::value>::gt(cA, cB, vA, vB);
  }

  /// <summary>
  /// Greater-than-or-equal-to comparison between 2 variables and returns 2 associated variable values
  /// @see gtx_associated_min_max
  /// </summary>
  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<L, T, Q> associatedGreaterEqual(vec<L, T, Q> const &cA, vec<L, T, Q> const &cB, vec<L, T, Q> const &vA, vec<L, T, Q> const &vB) {
    return detail::compute_associated<vec<L, T, Q>, detail::is_aligned<Q>::value>::gte(cA, cB, vA, vB);
  }

  /*
  ** {======================================================
  ** Patches
  ** =======================================================
  */

#if GLM_CONFIG_SIMD == GLM_ENABLE
  #if 0 && !(GLM_ARCH & GLM_ARCH_SSE41_BIT)
  /// <summary>
  /// @private
  /// @GLMFix: Avoid imprecise _mm_floor_ps emulation when GLM_ARCH_SSE41_BIT is
  ///          not enabled: +/-2^23 (8388608.0).
  /// </summary>
  namespace detail {
  }
  #endif

  #if (GLM_ARCH & GLM_ARCH_SSE2_BIT)
  /// <summary>
  /// @private
  /// @GLMFix: glm/glm/detail/func_integer_simd.inl:42:32: error: could not convert ‘add0’ from ‘const __m128i’ to ‘glm::vec<4, unsigned int, glm::aligned_highp>’
  /// </summary>
  namespace detail {
    template<>
    struct compute_bitfieldBitCountStep<4, uint, glm::aligned_highp, true, true> {
      GLM_FUNC_QUALIFIER static vec<4, uint, glm::aligned_highp> call(vec<4, uint, glm::aligned_highp> const &v, uint Mask, uint Shift) {
        __m128i const set0 = v.data;

        __m128i const set1 = _mm_set1_epi32(static_cast<int>(Mask));
        __m128i const and0 = _mm_and_si128(set0, set1);
        __m128i const sft0 = _mm_slli_epi32(set0, static_cast<int>(Shift));
        __m128i const and1 = _mm_and_si128(sft0, set1);
        __m128i const add0 = _mm_add_epi32(and0, and1);

        vec<4, uint, glm::aligned_highp> Result;
        Result.data = add0;
        return Result;
      }
    };

    template<>
    struct compute_bitfieldReverseStep<4, uint, glm::aligned_highp, true, true> {
      GLM_FUNC_QUALIFIER static vec<4, uint, glm::aligned_highp> call(vec<4, uint, glm::aligned_highp> const &v, uint Mask, uint Shift) {
        __m128i const set0 = v.data;

        __m128i const set1 = _mm_set1_epi32(static_cast<int>(Mask));
        __m128i const and1 = _mm_and_si128(set0, set1);
        __m128i const sft1 = _mm_slli_epi32(and1, static_cast<int>(Shift));
        __m128i const set2 = _mm_andnot_si128(set0, _mm_set1_epi32(-1));
        __m128i const and2 = _mm_and_si128(set0, set2);
        __m128i const sft2 = _mm_srai_epi32(and2, static_cast<int>(Shift));
        __m128i const or0 = _mm_or_si128(sft1, sft2);

        vec<4, uint, glm::aligned_highp> Result;
        Result.data = or0;
        return Result;
      }
    };
  }
  #endif
#endif

  /// <summary>
  /// @private
  /// @GLMFix: A implementation of glm::angle that is numerically stable at all angles.
  /// </summary>
  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER T _angle(vec<L, T, Q> const &x, vec<L, T, Q> const &y) {
    GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'angle' only accept floating-point inputs");
    const vec<L, T, Q> xyl = x * length(y);
    const vec<L, T, Q> yxl = y * length(x);
    const T n = length(xyl - yxl);
    return detail::approx_zero(n) ? T(0) : T(2) * atan2(n, length(xyl + yxl));
  }

  template<typename genType>
  GLM_FUNC_QUALIFIER genType _angle(genType const &x, genType const &y) {
    return angle<genType>(x, y);
  }

  /// <summary>
  /// orientedAngle implementation that uses _angle
  /// @see gtx_vector_angle
  /// @private
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER T _orientedAngle(vec<2, T, Q> const &x, vec<2, T, Q> const &y) {
    GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'orientedAngle' only accept floating-point inputs");
    const T Angle = _angle(x, y);
    const T partialCross = x.x * y.y - y.x * x.y;
    return (partialCross > T(0)) ? Angle : -Angle;
  }

  /// <summary>
  /// orientedAngle implementation that uses _angle
  /// @see gtx_vector_angle
  /// @private
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER T _orientedAngle(vec<3, T, Q> const &x, vec<3, T, Q> const &y, vec<3, T, Q> const &ref) {
    GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'orientedAngle' only accept floating-point inputs");
    const T Angle = _angle(x, y);
    return mix(Angle, -Angle, dot(ref, cross(x, y)) < T(0));
  }

  /// <summary>
  /// @private
  /// @GLMFix: Generalized slerp implementation; glm/ext/quaternion_common.inl
  /// </summary>
  template<length_t L, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<L, T, Q> _slerp(vec<L, T, Q> const &x, vec<L, T, Q> const &y, T const &a) {
    GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'slerp' only accept floating-point inputs");
    const T CosAlpha = dot(x, y);
    if (CosAlpha > static_cast<T>(1) - epsilon<T>())
      return mix(x, y, a);

    const T Alpha = acos(CosAlpha);  // get angle (0 -> pi)
    const T SinAlpha = sin(Alpha);  // get sine of angle between vectors (0 -> 1)
    const T t1 = sin((static_cast<T>(1) - a) * Alpha) / SinAlpha;
    const T t2 = sin(a * Alpha) / SinAlpha;
    return x * t1 + y * t2;
  }

  template<typename genType>
  GLM_FUNC_QUALIFIER genType _slerp(genType x, genType y, genType a) {
    return _slerp(vec<1, genType>(x), vec<1, genType>(y), a).x;
  }

  /* }====================================================== */

  /// @}
}
