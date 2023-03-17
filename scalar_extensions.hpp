/// @ref gtx_scalar_extensions
/// @file scalar_extensions.hpp
///
/// @defgroup gtx_scalar_extensions GLM_GTX_scalar_extensions
/// @ingroup gtx
///
/// GLM scalar extensions
///
/// @see core_func_common
/// @see ext_scalar_common

#pragma once
#if !defined(GLM_ENABLE_EXPERIMENTAL)
  #define GLM_ENABLE_EXPERIMENTAL
#endif

#include <cmath>
#include <limits>
#include <type_traits>

#include <glm/glm.hpp>
#include <glm/fwd.hpp>

#include <glm/detail/compute_vector_relational.hpp>
#include <glm/gtc/bitfield.hpp>
#include <glm/gtc/color_space.hpp>
#include <glm/gtc/epsilon.hpp>
#include <glm/gtc/packing.hpp>
#include <glm/gtx/common.hpp>
#include <glm/gtx/compatibility.hpp>
#include <glm/gtx/spline.hpp>
#include <glm/ext/scalar_common.hpp>

#include <glm/ext/scalar_constants.hpp>
#if GLM_HAS_CXX11_STL && (GLM_LANG & GLM_LANG_CXX14_FLAG)
  #include <functional>  // std::equal_to/std::not_equal_to
#endif

#if GLM_MESSAGES == GLM_ENABLE && !defined(GLM_EXT_INCLUDED)
  #pragma message("GLM: GLM_GTX_scalar_ext extension included")
#endif

// GLM_HAS_BUILTIN __has_builtin wrapper
#if !defined(GLM_HAS_BUILTIN)
#if defined(__has_builtin)
  #define GLM_HAS_BUILTIN(x) __has_builtin(x)
#else
  #define GLM_HAS_BUILTIN(x) 0
#endif
#endif

namespace glm {

  /// <summary>
  /// @private
  /// Functions to suppress float-equal warnings and extension to compute_vector_relational.hpp
  /// </summary>
  namespace detail {
    template<typename genType>
    GLM_FUNC_QUALIFIER GLM_CONSTEXPR bool equal_to(genType x, genType y) {
#if GLM_HAS_CXX11_STL && (GLM_LANG & GLM_LANG_CXX14_FLAG)
      return std::equal_to<genType>()(x, y);
#else
      return x == y;
#endif
    }

    template<typename genType>
    GLM_FUNC_QUALIFIER GLM_CONSTEXPR bool not_equal_to(genType x, genType y) {
#if GLM_HAS_CXX11_STL && (GLM_LANG & GLM_LANG_CXX14_FLAG)
      return std::not_equal_to<genType>()(x, y);
#else
      return x != y;
#endif
    }

    template<typename genType>
    GLM_FUNC_QUALIFIER GLM_CONSTEXPR bool exactly_zero(genType x) {
      return equal_to(x, genType(0));
    }

    template<typename genType>
    GLM_FUNC_QUALIFIER GLM_CONSTEXPR bool exactly_one(genType x) {
      return equal_to(x, genType(1));
    }

    template<typename T>
    struct compute_equal<T, true> {
      GLM_FUNC_QUALIFIER GLM_CONSTEXPR static bool call(T a, T b) {
        return equal_to(a, b);
      }
    };

    /* Approximately ... */

    template<typename genType>
    GLM_FUNC_QUALIFIER GLM_CONSTEXPR bool approx_zero(genType x, genType eps = epsilon<genType>()) {
      return epsilonEqual(x, genType(0), eps);
    }
  }

  /// @addtogroup gtx_scalar_extensions
  /// @{

  /* Constants */
#if 0
  /// @GLMEpsilon: Reduce glm::epsilon<double>()
  ///
  /// Some glm functions (extensions) are numerically unstable at the extremes
  /// and a small (and non-configurable) tolerance is unforgiving to any sort of
  /// precision-loss/drift (e.g., ffast-math, polynomial approximations, etc).
  /// Reference:
  ///   glm/ext/quaternion_common.inl
  ///   glm/ext/quaternion_exponential.inl
  ///   glm/gtc/quaternion.inl
  ///   glm/gtx/quaternion.inl
  ///   glm/gtx/rotate_vector.inl
  /// @private
  template<>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR double epsilon() {
    return static_cast<double>(glm::epsilon<float>());
  }
#endif

  /* Increase epsilon tolerances for functions that tend to drift */

  /// <summary>
  /// Precomputed sqrt(epsilon<genType>()) value.
  /// @see ext_scalar_constants
  /// </summary>
  template<typename genType>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR genType epsilon2() {
    GLM_STATIC_ASSERT(std::numeric_limits<genType>::is_iec559, "'epsilon2' only accepts floating-point inputs");
    return genType(1E-4);  // sqrt(epsilon<genType>());
  }

  /// <summary>
  /// Epsilon value for intersection tests.
  /// @see ext_scalar_constants
  /// </summary>
  template<typename genType>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR genType intersect_eps() {
    GLM_STATIC_ASSERT(std::numeric_limits<genType>::is_iec559, "'intersection_epsilon' only accepts floating-point inputs");
    return genType(1E-4);
  }

  /// <summary>
  /// Epsilon value for containment tests.
  /// @see ext_scalar_constants
  /// </summary>
  template<typename genType>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR genType contains_eps() {
    GLM_STATIC_ASSERT(std::numeric_limits<genType>::is_iec559, "'contains_eps' only accepts floating-point inputs");
    return genType(1E-3);
  }

  /// <summary>
  /// all(equal(x, y)) shorthand
  /// @private @see gtx_scalar_relational
  /// </summary>
  template<typename genType>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR bool all_equal(genType x, genType y) {
    return all(equal(x, y));
  }

  /// <summary>
  /// all(equal(x, y, eps)) shorthand
  /// @private @see gtx_scalar_relational
  /// </summary>
  template<typename genType>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR bool all_equal(genType x, genType y, genType eps) {
    return all(equal(x, y, eps));
  }

  /// <summary>
  /// all(equal(x, y, MaxULPs)) shorthand
  /// @private @see gtx_scalar_relational
  /// </summary>
  template<typename genType>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR bool all_equal(genType x, genType y, int MaxULPs) {
    return all(equal(x, y, MaxULPs));
  }

  /// <summary>
  /// any(notEqual(x, y)) shorthand
  /// @private @see gtx_scalar_relational
  /// </summary>
  template<typename genType>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR bool any_notEqual(genType x, genType y) {
    return any(notEqual(x, y));
  }

  /// <summary>
  /// any(notEqual(x, y, eps)) shorthand
  /// @private @see gtx_scalar_relational
  /// </summary>
  template<typename genType>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR bool any_notEqual(genType x, genType y, genType eps) {
    return any(notEqual(x, y, eps));
  }

  /// <summary>
  /// any(notEqual(x, y, eps, MaxULPs)) shorthand
  /// @private @see gtx_scalar_relational
  /// </summary>
  template<typename genType>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR bool any_notEqual(genType x, genType y, int MaxULPs) {
    return any(notEqual(x, y, MaxULPs));
  }

  /// <summary>
  /// all(lessThan(x, y)) shorthand
  /// @private @see gtx_scalar_relational
  /// </summary>
  template<typename genType>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR bool all_lessThan(genType x, genType y) {
    return all(lessThan(x, y));
  }

  /// <summary>
  /// all(lessThanEqual(x, y)) shorthand
  /// @private @see gtx_scalar_relational
  /// </summary>
  template<typename genType>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR bool all_lessThanEqual(genType x, genType y) {
    return all(lessThanEqual(x, y));
  }

  /// <summary>
  /// all(greaterThan(x, y)) shorthand
  /// @private @see gtx_scalar_relational
  /// </summary>
  template<typename genType>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR bool all_greaterThan(genType x, genType y) {
    return all(greaterThan(x, y));
  }

  /// <summary>
  /// all(greaterThanEqual(x, y)) shorthand
  /// @private @see gtx_scalar_relational
  /// </summary>
  template<typename genType>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR bool all_greaterThanEqual(genType x, genType y) {
    return all(greaterThanEqual(x, y));
  }

  /// <summary>
  /// any(isinf(x)) shorthand
  /// @private @see core_func_common
  /// </summary>
  template<typename genType>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR bool any_isinf(genType x) {
    return any(isinf(x));
  }

  /// <summary>
  /// all(isfinite(x)) shorthand
  /// @private @private @see core_func_common
  /// </summary>
  template<typename genType>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR bool all_isfinite(genType x) {
    return all(isfinite(x));
  }

  /// <summary>
  /// glm::any(glm::isnan(...)) shorthand
  /// @see core_func_common
  /// </summary>
  template<typename genType>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR bool any_isnan(genType x) {
    return any(isnan(x));
  }

  /// <summary>
  /// Returns 1.0 if >= 0, or –1.0 if x < 0.
  /// @see core_func_common
  /// </summary>
  template<typename genType>
  GLM_FUNC_QUALIFIER genType signP(genType x) {
    GLM_STATIC_ASSERT(std::numeric_limits<genType>::is_iec559 || (std::numeric_limits<genType>::is_signed && std::numeric_limits<genType>::is_integer), "'sign' only accept signed inputs");
    return (x >= genType(0)) ? genType(1) : genType(-1);
  }

  /// <summary>
  /// Returns 1.0 if x > 0, or –1.0 if x <= 0.
  /// @see core_func_common
  /// </summary>
  template<typename genType>
  GLM_FUNC_QUALIFIER genType signN(genType x) {
    GLM_STATIC_ASSERT(std::numeric_limits<genType>::is_iec559 || (std::numeric_limits<genType>::is_signed && std::numeric_limits<genType>::is_integer), "'sign' only accept signed inputs");
    return (x > genType(0)) ? genType(1) : genType(-1);
  }

  /// <summary>
  /// Calculate sin and cos simultaneously.
  ///
  /// @param[in] v value representing angle in radians
  /// @param[out] s the sine of value in the range [-1 ; +1]
  /// @param[out] c the cosine of value in the range [-1 ; +1]
  /// @see core_func_trigonometric
  /// </summary>
  template<typename genType>
  GLM_FUNC_QUALIFIER void sincos(genType v, genType &s, genType &c) {
    s = sin(v);
    c = cos(v);
  }

  /// <summary>
  /// Cardinal sine
  /// @see core_func_trigonometric
  /// </summary>
  template<typename genType>
  GLM_FUNC_QUALIFIER genType sinc(genType x) {
    return detail::exactly_zero(x) ? genType(1) : glm::sin(x) / x;
  }

  /// <summary>
  /// Normalized cardinal sine
  /// @see core_func_trigonometric
  /// </summary>
  template<typename genType>
  GLM_FUNC_QUALIFIER genType sincn(genType x) {
    return sinc(glm::pi<genType>() * x);
  }

  /// <summary>
  /// Logistic function with basic overflow handling; underflow to-be-determined.
  /// @see core_func_exponential
  /// </summary>
  template<typename genType>
  GLM_FUNC_QUALIFIER genType logistic(genType x) {
    GLM_STATIC_ASSERT(std::numeric_limits<genType>::is_iec559, "'logistic' only accept floating-point inputs.");
    const genType e = exp(min(x, genType(44.3614196)));  // exp(-44.3614196) ~ 2^{−64}
    return e / (genType(1.0) + e);
  }

  template<>
  GLM_FUNC_QUALIFIER float logistic(float x) {
    const float e = exp(min(x, 16.6355324f));  // exp(−16.6355324) ~ 2^{−24}
    return e / (1.0f + e);
  }

  /// <summary>
  /// Wraps 'value' between [0, maxValue]
  /// @see gtx_wrap
  /// </summary>
  template<typename genType>
  GLM_FUNC_QUALIFIER genType wrap(genType value, genType maxValue) {
    return fmod(value, maxValue) + ((value < genType(0)) ? maxValue : genType(0));
  }

  /// <summary>
  /// Wraps [minValue, maxValue].
  /// @see gtx_wrap
  /// </summary>
  // template<typename genType>
  // GLM_FUNC_QUALIFIER genType wrap(genType value, genType minValue, genType maxValue) {
  //   return wrap(value - minValue, maxValue - minValue) + minValue;
  // }

  /// <summary>
  /// Loops 't' so that it is never greater than 'length' and less than zero.
  /// @see gtx_wrap
  /// </summary>
  template<typename genType>
  GLM_FUNC_QUALIFIER genType loopRepeat(genType t, genType length) {
    if (detail::exactly_zero(length))
      return 0;
    return clamp(t - floor(t / length) * length, genType(0), length);
  }

  /// <summary>
  /// Returns a value that will increment and decrement 't' between 0 and
  /// 'length'. This function emulates Unity.Mathf.PingPong.
  /// @see gtx_wrap
  /// </summary>
  template<typename genType>
  GLM_FUNC_QUALIFIER genType pingPong(genType t, genType length) {
    t = loopRepeat(t, length * genType(2));
    return length - abs(t - length);
  }

  /// <summary>
  /// glm::wrapAngle defined over [-pi, pi]
  /// @see gtx_fast_trigonometry
  /// </summary>
  template<typename genType>
  GLM_FUNC_QUALIFIER genType wrapAngleSigned(genType value) {
    if (value >= genType(0))
      return fmod(value + pi<genType>(), two_pi<genType>()) - pi<genType>();
    return fmod(value - pi<genType>(), two_pi<genType>()) + pi<genType>();
  }

  /// <summary>
  /// Return the shortest difference between two angles (in radians).
  /// @see ext_scalar_common
  /// </summary>
  template<typename genType>
  GLM_FUNC_QUALIFIER genType deltaAngle(genType a, genType b) {
    const genType dt = loopRepeat((b - a), two_pi<genType>());
    return min(two_pi<genType>() - dt, dt);
  }

  /// <summary>
  /// A lerp implementation that ensures values interpolate correctly when
  /// wrapped around two-pi
  /// @see gtx_fast_trigonometry
  /// </summary>
  template<typename genType>
  GLM_FUNC_QUALIFIER genType lerpAngle(genType a, genType b, genType t) {
    GLM_STATIC_ASSERT(std::numeric_limits<genType>::is_iec559, "'lerpAngle' only accept floating-point inputs");
    const genType dt = loopRepeat((b - a), two_pi<genType>());
    return lerp(a, a + (dt > pi<genType>() ? dt - two_pi<genType>() : dt), t);
  }

  /// <summary>
  /// Return a position between two points, moving no further than maxDist.
  /// @see gtx_functions
  /// </summary>
  template<typename genType>
  GLM_FUNC_QUALIFIER genType moveTowards(genType current, genType target, genType maxDist) {
    GLM_STATIC_ASSERT(std::numeric_limits<genType>::is_iec559, "'moveTowards' only accept floating-point inputs");
    if (abs(target - current) <= maxDist)
      return target;
    return current + sign(target - current) * maxDist;
  }

  /// <summary>
  /// Round a 'value' to a specified grid size.
  /// @see gtc_round
  /// </summary>
  template<typename genType>
  GLM_FUNC_QUALIFIER genType snap(genType value, genType step) {
    if (!detail::exactly_zero(step))
      return floor((value / step) + genType(0.5)) * step;
    return value;
  }

  /// <summary>
  /// Returns the normalized vector pointing from "x" to "y"
  /// @see ext_vector_common
  /// </summary>
  template<typename genType>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR genType direction(genType x, genType y) {
    return normalize(vec<1, genType>(y - x)).x;
  }

  /// <summary>
  /// Returns a value 't' such that lerp(x, y, t) == value (or 0 if x == y)
  /// @see gtx_compatibility
  /// </summary>
  template<typename genType>
  GLM_FUNC_QUALIFIER genType lerpInverse(genType x, genType y, genType value) {
    return equal(x, y, epsilon<genType>()) ? genType(0) : (value - x) / (y - x);
  }

  /// <summary>
  /// Expands a 10-bit integer into 30 bits by inserting 2 zeros after each bit.
  ///
  /// @see <a href="https://developer.nvidia.com/blog/thinking-parallel-part-iii-tree-construction-gpu/">Thinking Parallel, Part III: Tree Construction on the GPU</a>
  /// </summary>
  GLM_FUNC_QUALIFIER unsigned int expandBits(unsigned int v) {
    v = (v * 0x00010001u) & 0xFF0000FFu;
    v = (v * 0x00000101u) & 0x0F00F00Fu;
    v = (v * 0x00000011u) & 0xC30C30C3u;
    v = (v * 0x00000005u) & 0x49249249u;
    return v;
  }

  /// <summary>
  /// Calculates a 30-bit Morton code for the given 3D point located within the unit cube [0,1].
  /// @see <a href="https://developer.nvidia.com/blog/thinking-parallel-part-iii-tree-construction-gpu/">Thinking Parallel, Part III: Tree Construction on the GPU</a>
  /// </summary>
  GLM_FUNC_QUALIFIER unsigned int morton3D(float x, float y, float z) {
    x = min(max(x * 1024.0f, 0.0f), 1023.0f);
    y = min(max(y * 1024.0f, 0.0f), 1023.0f);
    z = min(max(z * 1024.0f, 0.0f), 1023.0f);
    unsigned int xx = expandBits(static_cast<unsigned int>(x));
    unsigned int yy = expandBits(static_cast<unsigned int>(y));
    unsigned int zz = expandBits(static_cast<unsigned int>(z));
    return xx * 4 + yy * 2 + zz;
  }

  /// <summary>
  /// Calculates a 30-bit Morton code for the given 3D vector located within the unit cube.
  /// </summary>
  template<qualifier Q>
  GLM_FUNC_QUALIFIER unsigned int morton3D(vec<3, float, Q> const &x) {
    return morton3D(x.x, x.y, x.z);
  }

  /* glm/gtx/associated.hpp extensions */

  /// <summary>
  /// Equal comparison between 2 variables and returns 2 associated variable values
  /// @see gtx_associated_min_max
  /// </summary>
  template<typename genType>
  GLM_FUNC_QUALIFIER genType associatedEqual(genType cA, genType cB, genType vA, genType vB) {
    return detail::equal_to(cA, cB) ? vA : vB;
  }

  /// <summary>
  /// Greater-than comparison between 2 variables and returns 2 associated variable values
  /// @see gtx_associated_min_max
  /// </summary>
  template<typename genType>
  GLM_FUNC_QUALIFIER genType associatedGreater(genType cA, genType cB, genType vA, genType vB) {
    return (cA > cB) ? vA : vB;
  }

  /// <summary>
  /// Greater-than-or-equal-to comparison between 2 variables and returns 2 associated variable values
  /// @see gtx_associated_min_max
  /// </summary>
  template<typename genType>
  GLM_FUNC_QUALIFIER genType associatedGreaterEqual(genType cA, genType cB, genType vA, genType vB) {
    return (cA >= cB) ? vA : vB;
  }

  /* Functions with additional integral type support */

  template<typename T>
  GLM_FUNC_QUALIFIER typename std::enable_if<std::is_integral<T>::value, T>::type iceil(T x) {
    return x;
  }

  template<typename T>
  GLM_FUNC_QUALIFIER typename std::enable_if<std::is_floating_point<T>::value, T>::type iceil(T x) {
    return ceil(x);
  }

  template<typename T>
  GLM_FUNC_QUALIFIER typename std::enable_if<std::is_integral<T>::value, T>::type ifloor(T x) {
    return x;
  }

  template<typename T>
  GLM_FUNC_QUALIFIER typename std::enable_if<std::is_floating_point<T>::value, T>::type ifloor(T x) {
    return floor(x);
  }

  template<typename T>
  GLM_FUNC_QUALIFIER typename std::enable_if<std::is_integral<T>::value && std::is_signed<T>::value, T>::type imod(T x, T y) {
    return detail::exactly_zero(y) ? T(0) : ((x % y) + y) % y;
  }

  template<typename T>
  GLM_FUNC_QUALIFIER typename std::enable_if<std::is_integral<T>::value && !std::is_signed<T>::value, T>::type imod(T x, T y) {
    return x - y * (x / y);
  }

  template<typename T>
  GLM_FUNC_QUALIFIER typename std::enable_if<std::is_floating_point<T>::value, T>::type imod(T x, T y) {
    return mod(x, y);
  }

  template<typename T>
  GLM_FUNC_QUALIFIER typename std::enable_if<std::is_integral<T>::value, T>::type pow(T x, unsigned int y) {
    if (detail::exactly_zero(y))
      return x >= T(0) ? T(1) : T(-1);
    else {
      T result = x;
      for (unsigned int i = 1; i < y; ++i)
        result *= x;
      return result;
    }
  }

  /* Missing implicit genType support & API unification */

  template<typename genType>
  GLM_FUNC_QUALIFIER genType fclamp(genType x) {
    return fclamp(x, genType(0), genType(1));
  }

  template<typename genType>
  GLM_FUNC_QUALIFIER genType nlerp(genType x, genType y, genType t) {
    return lerp(x, y, t);
  }

  template<typename genType>
  GLM_FUNC_QUALIFIER bool isUniform(genType v, genType eps = epsilon<genType>()) {
    ((void)v);
    ((void)eps);
    return true;
  }

  template<typename genType>
  GLM_FUNC_QUALIFIER genType reverse(genType v) {
    return v;
  }

  template<typename genType>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR genType compAdd(genType v) {
    return v;
  }

  template<typename genType>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR genType compMul(genType v) {
    return v;
  }

  template<typename genType>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR genType compMin(genType v) {
    return v;
  }

  template<typename genType>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR genType compMax(genType v) {
    return v;
  }

  template<typename genType>
  GLM_FUNC_QUALIFIER genType normalize(genType x) {
    GLM_STATIC_ASSERT(std::numeric_limits<genType>::is_iec559, "'normalize' accepts only floating-point inputs");
    return x < genType(0) ? genType(-1) : genType(1);
  }

  template<typename genType>
  GLM_FUNC_QUALIFIER bool isNormalized(genType x, genType eps = epsilon<genType>()) {
    GLM_STATIC_ASSERT(std::numeric_limits<genType>::is_iec559, "'isNormalized' only accept floating-point inputs");
    return abs(x - static_cast<genType>(1)) <= static_cast<genType>(2) * eps;
  }

  template<typename genType>
  GLM_FUNC_QUALIFIER bool isNull(genType x, genType eps = epsilon<genType>()) {
    GLM_STATIC_ASSERT(std::numeric_limits<genType>::is_iec559, "'isNull' only accept floating-point inputs");
    return x <= eps;
  }

  template<typename genType>
  GLM_FUNC_QUALIFIER bool isCompNull(genType v, genType eps = epsilon<genType>()) {
    return abs(v) < eps;
  }

  template<typename genType>
  GLM_FUNC_QUALIFIER bool areOrthonormal(genType v0, genType v1, genType eps = epsilon<genType>()) {
    return isNormalized(v0, eps) && isNormalized(v1, eps) && (abs(dot(v0, v1)) <= eps);
  }

  template<typename genType>
  GLM_FUNC_QUALIFIER bool areOrthogonal(genType v0, genType v1, genType eps = epsilon<genType>()) {
    return areOrthogonal(vec<1, genType>(v0), vec<1, genType>(v1), eps);
  }

  template<typename genType>
  GLM_FUNC_QUALIFIER genType normalizeDot(genType x, genType y) {
    return normalizeDot(vec<1, genType>(x), vec<1, genType>(y));
  }

  template<typename genType>
  GLM_FUNC_QUALIFIER genType fastNormalizeDot(genType x, genType y) {
    return fastNormalizeDot(vec<1, genType>(x), vec<1, genType>(y));
  }

  template<typename genType>
  GLM_FUNC_QUALIFIER bool openBounded(genType Value, genType Min, genType Max) {
    return openBounded(vec<1, genType>(Value), vec<1, genType>(Min), vec<1, genType>(Max)).x;
  }

  template<typename genType>
  GLM_FUNC_QUALIFIER bool closeBounded(genType Value, genType Min, genType Max) {
    return closeBounded(vec<1, genType>(Value), vec<1, genType>(Min), vec<1, genType>(Max)).x;
  }

  template<typename floatType, typename T>
  GLM_FUNC_QUALIFIER floatType compNormalize(T x) {
    return compNormalize<floatType>(vec<1, T>(x)).x;
  }

  template<typename T, typename floatType>
  GLM_FUNC_QUALIFIER T compScale(floatType x) {
    return compScale<T>(vec<1, floatType>(x)).x;
  }

  GLM_FUNC_QUALIFIER uint16 packHalf(float v) {
    return packHalf<1, defaultp>(vec<1, float>(v)).x;
  }

  GLM_FUNC_QUALIFIER float unpackHalf(uint16 v) {
    return unpackHalf<1, defaultp>(vec<1, uint16>(v)).x;
  }

  template<typename uintType, typename floatType>
  GLM_FUNC_QUALIFIER uintType packUnorm(floatType v) {
    return packUnorm<uintType, 1, floatType, defaultp>(vec<1, floatType>(v)).x;
  }

  template<typename floatType, typename uintType>
  GLM_FUNC_QUALIFIER floatType unpackUnorm(uintType v) {
    return unpackUnorm<floatType, 1, uintType, defaultp>(vec<1, uintType>(v)).x;
  }

  template<typename intType, typename floatType>
  GLM_FUNC_QUALIFIER intType packSnorm(floatType v) {
    return packSnorm<intType, 1, floatType, defaultp>(vec<1, floatType>(v)).x;
  }

  template<typename floatType, typename intType>
  GLM_FUNC_QUALIFIER floatType unpackSnorm(intType v) {
    return unpackSnorm<floatType, 1, intType, defaultp>(vec<1, intType>(v)).x;
  }

  template<typename genType>
  GLM_FUNC_QUALIFIER genType catmullRom(genType v1, genType v2, genType v3, genType v4, genType s) {
    return catmullRom<vec<1, genType>>(vec<1, genType>(v1), vec<1, genType>(v2), vec<1, genType>(v3), vec<1, genType>(v4), s).x;
  }

  template<typename genType>
  GLM_FUNC_QUALIFIER genType hermite(genType v1, genType t1, genType v2, genType t2, genType s) {
    return hermite<vec<1, genType>>(vec<1, genType>(v1), vec<1, genType>(t1), vec<1, genType>(v2), vec<1, genType>(t2), s).x;
  }

  template<typename genType>
  GLM_FUNC_QUALIFIER genType cubic(genType v1, genType v2, genType v3, genType v4, genType s) {
    return cubic<vec<1, genType>>(vec<1, genType>(v1), vec<1, genType>(v2), vec<1, genType>(v3), vec<1, genType>(v4), s).x;
  }

  template<typename genType>
  GLM_FUNC_QUALIFIER genType convertLinearToSRGB(genType ColorLinear) {
    return convertLinearToSRGB(vec<1, genType>(ColorLinear)).x;
  }

  template<typename genType>
  GLM_FUNC_QUALIFIER genType convertLinearToSRGB(genType ColorLinear, genType Gamma) {
    return convertLinearToSRGB(vec<1, genType>(ColorLinear), Gamma).x;
  }

  template<typename genType>
  GLM_FUNC_QUALIFIER genType convertSRGBToLinear(genType ColorSRGB) {
    return convertSRGBToLinear(vec<1, genType>(ColorSRGB)).x;
  }

  template<typename genType>
  GLM_FUNC_QUALIFIER genType convertSRGBToLinear(genType ColorSRGB, genType Gamma) {
    return convertSRGBToLinear(vec<1, genType>(ColorSRGB), Gamma).x;
  }

  /*
  ** {======================================================
  ** Patches
  ** =======================================================
  */

  /// <summary>
  /// @GLMFix: -Werror when using bitfieldFillOne and bitfieldFillZero:
  /// libs/glm/glm/./gtc/bitfield.inl:229:29: warning: comparison of integer expressions of different signedness: ‘int’ and ‘long unsigned int’ [-Wsign-compare]
  ///   229 | return Bits >= sizeof(genIUType) * 8 ? ~static_cast<genIUType>(0) : (static_cast<genIUType>(1) << Bits) - static_cast<genIUType>(1);
  ///
  /// @private copied from glm/detail/func_integer.inl
  /// </summary>
  template<>
  GLM_FUNC_QUALIFIER int mask(int Bits) {
    return Bits >= static_cast<int>(sizeof(int) * 8) ? ~static_cast<int>(0) : (static_cast<int>(1) << Bits) - static_cast<int>(1);
  }

  /// <summary>
  /// @private
  /// @GLMFix: missing long double support
  /// </summary>
  template<typename genType>
  GLM_FUNC_QUALIFIER bool epsilonEqual(genType const &x, genType const &y, genType const &epsilon) {
    return abs(x - y) < epsilon;
  }

  /// <summary>
  /// @private
  /// @GLMFix: missing long double support
  /// </summary>
  template<typename genType>
  GLM_FUNC_QUALIFIER bool epsilonNotEqual(genType const &x, genType const &y, genType const &epsilon) {
    return abs(x - y) >= epsilon;
  }

  /* }====================================================== */

  /// @}
}
