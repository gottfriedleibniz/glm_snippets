/// @ref gtx_quat_extensions
/// @file quat_extensions.hpp
///
/// @defgroup gtx_quat_extensions GLM_GTX_quat_extensions
/// @ingroup gtx
///
/// GLM quaternion extensions
///  1. API unifying functions
///  2. Functions that exist for rotation matrices but not for quaternions
///  3. Functions emulated/ported from other popular vector-math libraries
///
/// @see core_func_common
/// @see ext_quaternion_common

#pragma once
#if !defined(GLM_ENABLE_EXPERIMENTAL)
  #define GLM_ENABLE_EXPERIMENTAL
#endif

#include <limits>

#include <glm/glm.hpp>
#include <glm/ext/quaternion_common.hpp>
#include <glm/gtx/quaternion.hpp>
#include <glm/gtx/euler_angles.hpp>
#include <glm/gtx/fast_square_root.hpp>
#include <glm/ext/matrix_transform.hpp>
#include <glm/ext/quaternion_common.hpp>
#include <glm/ext/quaternion_trigonometric.hpp>
#include <glm/ext/quaternion_geometric.hpp>
#include <glm/ext/scalar_constants.hpp>

#include "scalar_extensions.hpp"
#include "vector_extensions.hpp"
#include "matrix_extensions.hpp"

#if GLM_MESSAGES == GLM_ENABLE && !defined(GLM_EXT_INCLUDED)
  #pragma message("GLM: GLM_GTX_quat_ext extension included")
#endif

namespace glm {

  /// @addtogroup gtx_quat_extensions
  /// @{

  /* EulerAngles -> Quaternion; @TODO: Optimize */

  template<typename T, qualifier Q = defaultp> GLM_FUNC_QUALIFIER qua<T, Q> quatEulerAngleX(T angleX) { return toQuat(eulerAngleX(angleX)); }
  template<typename T, qualifier Q = defaultp> GLM_FUNC_QUALIFIER qua<T, Q> quatEulerAngleY(T angleY) { return toQuat(eulerAngleY(angleY)); }
  template<typename T, qualifier Q = defaultp> GLM_FUNC_QUALIFIER qua<T, Q> quatEulerAngleZ(T angleZ) { return toQuat(eulerAngleZ(angleZ)); }
  template<typename T, qualifier Q = defaultp> GLM_FUNC_QUALIFIER qua<T, Q> quatEulerAngleXY(T angleX, T angleY) { return toQuat(eulerAngleXY(angleX, angleY)); }
  template<typename T, qualifier Q = defaultp> GLM_FUNC_QUALIFIER qua<T, Q> quatEulerAngleXZ(T angleX, T angleZ) { return toQuat(eulerAngleXZ(angleX, angleZ)); }
  template<typename T, qualifier Q = defaultp> GLM_FUNC_QUALIFIER qua<T, Q> quatEulerAngleYX(T angleY, T angleX) { return toQuat(eulerAngleYX(angleY, angleX)); }
  template<typename T, qualifier Q = defaultp> GLM_FUNC_QUALIFIER qua<T, Q> quatEulerAngleYZ(T angleY, T angleZ) { return toQuat(eulerAngleYZ(angleY, angleZ)); }
  template<typename T, qualifier Q = defaultp> GLM_FUNC_QUALIFIER qua<T, Q> quatEulerAngleZX(T angleZ, T angleX) { return toQuat(eulerAngleZX(angleZ, angleX)); }
  template<typename T, qualifier Q = defaultp> GLM_FUNC_QUALIFIER qua<T, Q> quatEulerAngleZY(T angleZ, T angleY) { return toQuat(eulerAngleZY(angleZ, angleY)); }
  template<typename T, qualifier Q = defaultp> GLM_FUNC_QUALIFIER qua<T, Q> quatEulerAngleXYX(T t1, T t2, T t3) { return toQuat(eulerAngleXYX(t1, t2, t3)); }
  template<typename T, qualifier Q = defaultp> GLM_FUNC_QUALIFIER qua<T, Q> quatEulerAngleXZX(T t1, T t2, T t3) { return toQuat(eulerAngleXZX(t1, t2, t3)); }
  template<typename T, qualifier Q = defaultp> GLM_FUNC_QUALIFIER qua<T, Q> quatEulerAngleYXY(T t1, T t2, T t3) { return toQuat(eulerAngleYXY(t1, t2, t3)); }
  template<typename T, qualifier Q = defaultp> GLM_FUNC_QUALIFIER qua<T, Q> quatEulerAngleYZY(T t1, T t2, T t3) { return toQuat(eulerAngleYZY(t1, t2, t3)); }
  template<typename T, qualifier Q = defaultp> GLM_FUNC_QUALIFIER qua<T, Q> quatEulerAngleZXZ(T t1, T t2, T t3) { return toQuat(eulerAngleZXZ(t1, t2, t3)); }
  template<typename T, qualifier Q = defaultp> GLM_FUNC_QUALIFIER qua<T, Q> quatEulerAngleZYZ(T t1, T t2, T t3) { return toQuat(eulerAngleZYZ(t1, t2, t3)); }

  /* Quaternion -> EulerAngles; @TODO: Optimize */

  template<typename T, qualifier Q> GLM_FUNC_QUALIFIER void extractEulerAngleXYX(qua<T, Q> const &q, T &t1, T &t2, T &t3) { extractEulerAngleXYX(toMat3(q), t1, t2, t3); }
  template<typename T, qualifier Q> GLM_FUNC_QUALIFIER void extractEulerAngleXYZ(qua<T, Q> const &q, T &t1, T &t2, T &t3) { extractEulerAngleXYZ(toMat3(q), t1, t2, t3); }
  template<typename T, qualifier Q> GLM_FUNC_QUALIFIER void extractEulerAngleXZX(qua<T, Q> const &q, T &t1, T &t2, T &t3) { extractEulerAngleXZX(toMat3(q), t1, t2, t3); }
  template<typename T, qualifier Q> GLM_FUNC_QUALIFIER void extractEulerAngleXZY(qua<T, Q> const &q, T &t1, T &t2, T &t3) { extractEulerAngleXZY(toMat3(q), t1, t2, t3); }
  template<typename T, qualifier Q> GLM_FUNC_QUALIFIER void extractEulerAngleYXY(qua<T, Q> const &q, T &t1, T &t2, T &t3) { extractEulerAngleYXY(toMat3(q), t1, t2, t3); }
  template<typename T, qualifier Q> GLM_FUNC_QUALIFIER void extractEulerAngleYXZ(qua<T, Q> const &q, T &t1, T &t2, T &t3) { extractEulerAngleYXZ(toMat3(q), t1, t2, t3); }
  template<typename T, qualifier Q> GLM_FUNC_QUALIFIER void extractEulerAngleYZX(qua<T, Q> const &q, T &t1, T &t2, T &t3) { extractEulerAngleYZX(toMat3(q), t1, t2, t3); }
  template<typename T, qualifier Q> GLM_FUNC_QUALIFIER void extractEulerAngleYZY(qua<T, Q> const &q, T &t1, T &t2, T &t3) { extractEulerAngleYZY(toMat3(q), t1, t2, t3); }
  template<typename T, qualifier Q> GLM_FUNC_QUALIFIER void extractEulerAngleZXY(qua<T, Q> const &q, T &t1, T &t2, T &t3) { extractEulerAngleZXY(toMat3(q), t1, t2, t3); }
  template<typename T, qualifier Q> GLM_FUNC_QUALIFIER void extractEulerAngleZXZ(qua<T, Q> const &q, T &t1, T &t2, T &t3) { extractEulerAngleZXZ(toMat3(q), t1, t2, t3); }
  template<typename T, qualifier Q> GLM_FUNC_QUALIFIER void extractEulerAngleZYX(qua<T, Q> const &q, T &t1, T &t2, T &t3) { extractEulerAngleZYX(toMat3(q), t1, t2, t3); }
  template<typename T, qualifier Q> GLM_FUNC_QUALIFIER void extractEulerAngleZYZ(qua<T, Q> const &q, T &t1, T &t2, T &t3) { extractEulerAngleZYZ(toMat3(q), t1, t2, t3); }

  /* EulerAngles -> Quaternion (Optimized) */

  /// <summary>
  /// Creates a quaternion from euler angles (X * Y * Z).
  /// @see gtx_euler_angles
  /// </summary>
  template<typename T, qualifier Q = defaultp>
  GLM_FUNC_QUALIFIER qua<T, Q> quatEulerAngleXYZ(T t1, T t2, T t3) {
    GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'quatEulerAngle' only accept floating-point inputs");

    vec<3, T, Q> s, c;
    sincos(T(0.5) * vec<3, T, Q>(t1, t2, t3), s, c);
    return qua<T, Q>(
      c.x * c.y * c.z - s.x * s.y * s.z,
      s.x * c.y * c.z + c.x * s.y * s.z,
      c.x * s.y * c.z - s.x * c.y * s.z,
      c.x * c.y * s.z + s.x * s.y * c.z
    );
  }

  /// <summary>
  /// Creates a quaternion from euler angles (X * Z * Y).
  /// @see gtx_euler_angles
  /// </summary>
  template<typename T, qualifier Q = defaultp>
  GLM_FUNC_QUALIFIER qua<T, Q> quatEulerAngleXZY(T t1, T t2, T t3) {
    GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'quatEulerAngle' only accept floating-point inputs");

    vec<3, T, Q> s, c;
    sincos(T(0.5) * vec<3, T, Q>(t1, t3, t2), s, c);
    return qua<T, Q>(
      c.x * c.y * c.z + s.x * s.y * s.z,
      s.x * c.y * c.z - c.x * s.y * s.z,
      c.x * s.y * c.z - s.x * c.y * s.z,
      c.x * c.y * s.z + s.x * s.y * c.z
    );
  }

  /// <summary>
  /// Creates a quaternion from euler angles (Y * X * Z).
  /// @see gtx_euler_angles
  /// </summary>
  template<typename T, qualifier Q = defaultp>
  GLM_FUNC_QUALIFIER qua<T, Q> quatEulerAngleYXZ(T t1, T t2, T t3) {
    GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'quatEulerAngle' only accept floating-point inputs");

    vec<3, T, Q> s, c;
    sincos(T(0.5) * vec<3, T, Q>(t2, t1, t3), s, c);
    return qua<T, Q>(
      c.x * c.y * c.z + s.x * s.y * s.z,
      s.x * c.y * c.z + c.x * s.y * s.z,
      c.x * s.y * c.z - s.x * c.y * s.z,
      c.x * c.y * s.z - s.x * s.y * c.z
    );
  }

  /// <summary>
  /// Creates a quaternion from euler angles (Y * Z * X).
  /// @see gtx_euler_angles
  /// </summary>
  template<typename T, qualifier Q = defaultp>
  GLM_FUNC_QUALIFIER qua<T, Q> quatEulerAngleYZX(T t1, T t2, T t3) {
    GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'quatEulerAngle' only accept floating-point inputs");

    vec<3, T, Q> s, c;
    sincos(T(0.5) * vec<3, T, Q>(t3, t1, t2), s, c);
    return qua<T, Q>(
      c.x * c.y * c.z - s.x * s.y * s.z,
      s.x * c.y * c.z + c.x * s.y * s.z,
      c.x * s.y * c.z + s.x * c.y * s.z,
      c.x * c.y * s.z - s.x * s.y * c.z
    );
  }

  /// <summary>
  /// Creates a quaternion from euler angles (Z * X * Y).
  /// @see gtx_euler_angles
  /// </summary>
  template<typename T, qualifier Q = defaultp>
  GLM_FUNC_QUALIFIER qua<T, Q> quatEulerAngleZXY(T t1, T t2, T t3) {
    GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'quatEulerAngle' only accept floating-point inputs");

    vec<3, T, Q> s, c;
    sincos(T(0.5) * vec<3, T, Q>(t2, t3, t1), s, c);
    return qua<T, Q>(
      c.x * c.y * c.z - s.x * s.y * s.z,
      s.x * c.y * c.z - c.x * s.y * s.z,
      c.x * s.y * c.z + s.x * c.y * s.z,
      c.x * c.y * s.z + s.x * s.y * c.z
    );
  }

  /// <summary>
  /// Creates a quaternion from euler angles (Z * Y * X).
  /// @see gtx_euler_angles
  /// </summary>
  template<typename T, qualifier Q = defaultp>
  GLM_FUNC_QUALIFIER qua<T, Q> quatEulerAngleZYX(T t1, T t2, T t3) {
    GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'quatEulerAngle' only accept floating-point inputs");

    vec<3, T, Q> s, c;
    sincos(T(0.5) * vec<3, T, Q>(t3, t2, t1), s, c);
    return qua<T, Q>(
      c.x * c.y * c.z + s.x * s.y * s.z,
      s.x * c.y * c.z - c.x * s.y * s.z,
      c.x * s.y * c.z + s.x * c.y * s.z,
      c.x * c.y * s.z - s.x * s.y * c.z
    );
  }

  /* quaternion-as-vector4 operations */

  /// <summary>
  /// Returns the component-wise comparison of result x == y.
  /// @see ext_quaternion_relational
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR vec<4, bool, Q> equal(qua<T, Q> const &x, qua<T, Q> const &y, int MaxULPs) {
    return equal(x, y, vec<4, int, Q>(MaxULPs));
  }

  /// <summary>
  /// Returns the component-wise comparison of result x == y.
  /// @see ext_quaternion_relational
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<4, bool, Q> equal(qua<T, Q> const &x, qua<T, Q> const &y, vec<4, T, Q> const &eps) {
    const vec<4, T, Q> v(x.x - y.x, x.y - y.y, x.z - y.z, x.w - y.w);
    return lessThan(abs(v), eps);
  }

  /// <summary>
  /// Returns the component-wise comparison of result x == y.
  /// @see ext_quaternion_relational
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR vec<4, bool, Q> equal(qua<T, Q> const &x, qua<T, Q> const &y, vec<4, int, Q> const &MaxULPs) {
    return equal(vec<4, T, Q>(x.x, x.y, x.z, x.w), vec<4, T, Q>(y.x, y.y, y.z, y.w), MaxULPs);
  }

  /// <summary>
  /// Returns the component-wise comparison of result x != y.
  /// @see ext_quaternion_relational
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR vec<4, bool, Q> notEqual(qua<T, Q> const &x, qua<T, Q> const &y, int MaxULPs) {
    return notEqual(x, y, vec<4, int, Q>(MaxULPs));
  }

  /// <summary>
  /// Returns the component-wise comparison of result x != y.
  /// @see ext_quaternion_relational
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR vec<4, bool, Q> notEqual(qua<T, Q> const &x, qua<T, Q> const &y, vec<4, int, Q> const &MaxULPs) {
    return not_(equal(x, y, MaxULPs));
  }

  /// <summary>
  /// Returns the component-wise comparison of result x != y.
  /// @see ext_quaternion_relational
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<4, bool, Q> notEqual(qua<T, Q> const &x, qua<T, Q> const &y, vec<4, T, Q> const &eps) {
    const vec<4, T, Q> v(x.x - y.x, x.y - y.y, x.z - y.z, x.w - y.w);
    return greaterThanEqual(abs(v), eps);
  }

  /// <summary>
  /// all(equal(x, y)) shorthand
  /// @private @see ext_quaternion_relational
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR bool all_equal(qua<T, Q> const &x, qua<T, Q> const &y) {
    return all(equal(x, y));
  }

  /// <summary>
  /// all(equal(x, y, eps)) shorthand
  /// @private @see ext_quaternion_relational
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR bool all_equal(qua<T, Q> const &x, qua<T, Q> const &y, T eps) {
    return all(equal(x, y, eps));
  }

  /// <summary>
  /// all(equal(x, y, MaxULPs)) shorthand
  /// @private @see ext_quaternion_relational
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR bool all_equal(qua<T, Q> const &x, qua<T, Q> const &y, int MaxULPs) {
    return all(equal(x, y, MaxULPs));
  }

  /// <summary>
  /// all(equal(x, y, eps)) shorthand
  /// @private @see ext_quaternion_relational
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR bool all_equal(qua<T, Q> const &x, qua<T, Q> const &y, vec<4, T, Q> const &eps) {
    return all(equal(x, y, eps));
  }

  /// <summary>
  /// all(equal(x, y, MaxULPs)) shorthand
  /// @private @see ext_quaternion_relational
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR bool all_equal(qua<T, Q> const &x, qua<T, Q> const &y, vec<4, int, Q> const &MaxULPs) {
    return all(equal(x, y, MaxULPs));
  }

  /// <summary>
  /// any(notEqual(x, y)) shorthand
  /// @private @see ext_quaternion_relational
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR bool any_notEqual(qua<T, Q> const &x, qua<T, Q> const &y) {
    return any(notEqual(x, y));
  }

  /// <summary>
  /// any(notEqual(x, y, eps)) shorthand
  /// @private @see ext_quaternion_relational
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR bool any_notEqual(qua<T, Q> const &x, qua<T, Q> const &y, T eps) {
    return any(notEqual(x, y, eps));
  }

  /// <summary>
  /// any(notEqual(x, y, MaxULPs)) shorthand
  /// @private @see ext_quaternion_relational
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR bool any_notEqual(qua<T, Q> const &x, qua<T, Q> const &y, int MaxULPs) {
    return any(notEqual(x, y, MaxULPs));
  }

  /// <summary>
  /// any(notEqual(x, y, eps)) shorthand
  /// @private @see ext_quaternion_relational
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR bool any_notEqual(qua<T, Q> const &x, qua<T, Q> const &y, vec<4, T, Q> const &eps) {
    return any(notEqual(x, y, eps));
  }

  /// <summary>
  /// any(notEqual(x, y, MaxULPs)) shorthand
  /// @private @see ext_quaternion_relational
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR bool any_notEqual(qua<T, Q> const &x, qua<T, Q> const &y, vec<4, int, Q> const &MaxULPs) {
    return any(notEqual(x, y, MaxULPs));
  }

  /// <summary>
  /// all(lessThan(x, y)) shorthand
  /// @private @see gtc_quaternion
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR bool all_lessThan(qua<T, Q> const &x, qua<T, Q> const &y) {
    return all(lessThan(x, y));
  }

  /// <summary>
  /// all(lessThanEqual(x, y)) shorthand
  /// @private @see gtc_quaternion
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR bool all_lessThanEqual(qua<T, Q> const &x, qua<T, Q> const &y) {
    return all(lessThanEqual(x, y));
  }

  /// <summary>
  /// all(greaterThan(x, y)) shorthand
  /// @private @see gtc_quaternion
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR bool all_greaterThan(qua<T, Q> const &x, qua<T, Q> const &y) {
    return all(greaterThan(x, y));
  }

  /// <summary>
  /// all(greaterThanEqual(x, y)) shorthand
  /// @private @see gtc_quaternion
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR bool all_greaterThanEqual(qua<T, Q> const &x, qua<T, Q> const &y) {
    return all(greaterThanEqual(x, y));
  }

  /// <summary>
  /// Test whether each quaternion component is a finite value
  /// @see gtx_compatibility
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<4, bool, Q> isfinite(qua<T, Q> const &x) {
    return vec<4, bool, Q>(isfinite(x.x), isfinite(x.y), isfinite(x.z), isfinite(x.w));
  }

  /// <summary>
  /// Get the shortest equivalent of the rotation.
  /// @see ext_quaternion_common
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER qua<T, Q> shortestEquivalent(qua<T, Q> const &q) {
    return (q.w < T(0.0)) ? -q : q;
  }

  /// <summary>
  /// Return true if the quaternion is invertible, i.e., is non-zero and finite.
  /// @see ext_quaternion_common
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER bool invertible(qua<T, Q> const &q, T eps = epsilon<T>()) {
    return all(isfinite(q)) && length2(q) > eps;
  }

  /// <summary>
  /// normalized lerp that does not sanitize 't'.
  /// @see ext_quaternion_common
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER qua<T, Q> nlerp(qua<T, Q> const &x, qua<T, Q> const &y, T a) {
    GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'nlerp' only accept floating-point inputs");
    return normalize(x * (static_cast<T>(1) - a) + (y * a));
  }

  /// <summary>
  /// Create a quaternion from barycentric coordinates.
  /// @see ext_quaternion_common
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER qua<T, Q> barycentric(qua<T, Q> const &q1, qua<T, Q> const &q2, qua<T, Q> const &q3, T u, T v) {
    GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'barycentric' only accept floating-point inputs");
    const qua<T, Q> start = slerp(q1, q2, u + v);
    const qua<T, Q> end = slerp(q1, q3, u + v);
    return slerp(start, end, v / (u + v));
  }

  /// <summary>
  /// Return the absolute angle between two quaternions.
  /// @see ext_quaternion_trigonometric
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER T angle(qua<T, Q> const &x, qua<T, Q> const &y) {
    return deltaAngle(T(0), angle(y * conjugate(x)));
  }

  /// <summary>
  /// Return the oriented angle between two quaternions based on a reference axis.
  /// @see ext_quaternion_trigonometric
  /// @see gtx_vector_angle
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER T orientedAngle(qua<T, Q> const &x, qua<T, Q> const &y, vec<3, T, Q> const &ref) {
    const qua<T, Q> rot = y * conjugate(x);
    return deltaAngle(T(0), angle(rot)) * sign(dot(ref, axis(rot)));
  }

  /// <summary>
  /// Rotate a quaternion 'x' towards a given target 'y', rotating no more than 'maxRadians'.
  /// @see gtx_functions
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER qua<T, Q> rotateTowards(qua<T, Q> const &x, qua<T, Q> const &y, T maxRadians) {
    const T q_angle = angle(x, y);
    return detail::approx_zero(q_angle) ? y : slerp(x, y, min(T(1), maxRadians / q_angle));
  }

  /// <summary>
  /// quatLookAt alternative (from O3DE).
  ///
  /// @see gtc_quaternion
  /// @see gtx_handed_coordinate_space
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER qua<T, Q> fromBasis(vec<3, T, Q> const &basisX, vec<3, T, Q> const &basisY, vec<3, T, Q> const &basisZ) {
    T trace(0);
    qua<T, Q> result = identity<qua<T, Q>>();
    if (basisZ.z < T(0)) {
      if (basisX.x > basisY.y) {
        trace = T(1.0) + basisX.x - basisY.y - basisZ.z;
        result = qua<T, Q>(basisY.z - basisZ.y, trace, basisX.y + basisY.x, basisZ.x + basisX.z);
      }
      else {
        trace = T(1.0) - basisX.x + basisY.y - basisZ.z;
        result = qua<T, Q>(basisZ.x - basisX.z, basisX.y + basisY.x, trace, basisY.z + basisZ.y);
      }
    }
    else {
      if (basisX.x < -basisY.y) {
        trace = T(1.0) - basisX.x - basisY.y + basisZ.z;
        result = qua<T, Q>(basisX.y - basisY.x, basisZ.x + basisX.z, basisY.z + basisZ.y, trace);
      }
      else {
        trace = T(1.0) + basisX.x + basisY.y + basisZ.z;
        result = qua<T, Q>(trace, basisY.z - basisZ.y, basisZ.x - basisX.z, basisX.y - basisY.x);
      }
    }
    return result * (T(0.5) * inversesqrt(trace));
  }

  /// <summary>
  /// Create the (shortest arc) quaternion that rotates a source direction to
  /// coincide with the target. This function is a function wrapper to the quat
  /// constructor.
  ///
  /// @see ext_quaternion_transform
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR qua<T, Q> rotateFromTo(vec<3, T, Q> const &sourceDirection, vec<3, T, Q> const &targetDirection) {
    return qua<T, Q>(normalize(sourceDirection), normalize(targetDirection));
  }

  /// <summary>
  /// Create a right-handed spherical billboard that rotates around a specified 'object' position.
  /// @see ext_quaternion_transform
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER qua<T, Q> quatbillboardRH(vec<3, T, Q> const &object, vec<3, T, Q> const &camPos, vec<3, T, Q> const &camUp, vec<3, T, Q> const &camFwd) {
    return toQuat(billboardRH<T, Q, 3, 3>(object, camPos, camUp, camFwd));
  }

  /// <summary>
  /// Create a left-handed spherical billboard that rotates around a specified 'object' position.
  /// @see ext_quaternion_transform
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER qua<T, Q> quatbillboardLH(vec<3, T, Q> const &object, vec<3, T, Q> const &camPos, vec<3, T, Q> const &camUp, vec<3, T, Q> const &camFwd) {
    return toQuat(billboardLH<T, Q, 3, 3>(object, camPos, camUp, camFwd));
  }

  /// <summary>
  /// Create a spherical billboard that rotates around a specified 'object' position;  uses default handedness.
  /// @see ext_quaternion_transform
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER qua<T, Q> quatbillboard(vec<3, T, Q> const &object, vec<3, T, Q> const &pos, vec<3, T, Q> const &up, vec<3, T, Q> const &forward) {
#if (GLM_CONFIG_CLIP_CONTROL & GLM_CLIP_CONTROL_LH_BIT)
    return quatbillboardLH<T, Q>(object, pos, up, forward);
#else
    return quatbillboardRH<T, Q>(object, pos, up, forward);
#endif
  }

  /// <summary>
  /// Given an axis, return the portion of the rotation that accounts for the twist about that axis.
  /// @see gtc_quaternion
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER qua<T, Q> twist(qua<T, Q> const &q, vec<3, T, Q> const &ref) {
    const vec<3, T, Q> xyz = dot(vec<3, T, Q>(q.x, q.y, q.z), ref) * ref;
    const qua<T, Q> twist(q.w, xyz);
    const T twist_len = length2(twist);
    return !detail::exactly_zero(twist_len) ? (twist / sqrt(twist_len)) : identity<qua<T, Q>>();
  }

  /// <summary>
  /// Decompose a quaternion into two concatenated rotations: swing (Y/Z axes) and twist (X axis).
  /// @see gtc_quaternion
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER void swingtwist(qua<T, Q> const &q, qua<T, Q> &outSwing, qua<T, Q> &outTwist) {
    const T s = sqrt(q.w * q.w + q.x * q.x);
    if (!detail::exactly_zero(s)) {
      outTwist = qua<T, Q>(q.w / s, q.x / s, T(0), T(0));
      outSwing = qua<T, Q>(s, T(0), (q.w * q.y - q.x * q.z) / s, (q.w * q.z + q.x * q.y) / s);
    }
    else {
      outTwist = identity<qua<T, Q>>();
      outSwing = q;
    }
  }

  /* API unification */
  /* Explicit support for all rotation matrices */

  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER qua<T, Q> quat_cast(qua<T, Q> const &q) {
    return q;
  }

  /// Converts a pure rotation 3 * 4 matrix to a quaternion.
  ///
  /// @tparam T Floating-point scalar types.
  /// @see gtc_quaternion
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER qua<T, Q> quat_cast(mat<3, 4, T, Q> const &m) {
    return quat_cast(mat<3, 3, T, Q>(m));
  }

  /// Converts a pure rotation 4 * 3 matrix to a quaternion.
  ///
  /// @tparam T Floating-point scalar types.
  /// @see gtc_quaternion
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER qua<T, Q> quat_cast(mat<4, 3, T, Q> const &m) {
    return quat_cast(mat<3, 3, T, Q>(m));
  }

  /* glm/gtx/quaternion.hpp */

  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER qua<T, Q> toQuat(qua<T, Q> const &q) {
    return q;
  }

  /// Converts a 3 * 4 matrix to a quaternion.
  /// @see gtx_quaternion
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER qua<T, Q> toQuat(mat<3, 4, T, Q> const &m) {
    return quat_cast(m);
  }

  /// Converts a 4 * 3 matrix to a quaternion.
  /// @see gtx_quaternion
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER qua<T, Q> toQuat(mat<4, 3, T, Q> const &m) {
    return quat_cast(m);
  }

  /* vector_query.hpp unification */

  /// Check whether a quaternion is normalized.
  /// @see gtx_vector_query
  /// @see gtx_matrix_query
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER bool isNormalized(qua<T, Q> const &q, T eps = epsilon<T>()) {
    GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'isNormalized' only accept floating-point inputs");
    return abs(length(q) - static_cast<T>(1)) <= static_cast<T>(2) * eps;
  }

  /// Check whether a quaternion is null.
  /// @see gtx_vector_query
  /// @see gtx_matrix_query
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER bool isNull(qua<T, Q> const &q, T eps = epsilon<T>()) {
    GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'isNull' only accept floating-point inputs");
    return length(q) <= eps;
  }

  /* fast_square_root.hpp unification */

  /// Faster than the common normalize function but less accurate.
  /// @see gtx_fast_square_root (dependence)
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER qua<T, Q> fastNormalize(qua<T, Q> const &x) {
    return x * fastInverseSqrt<T>(dot(x, x));
  }

  /// Normalize parameters and returns the dot product of x and y.
  /// @see gtx_normalize_dot
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER T normalizeDot(qua<T, Q> const &x, qua<T, Q> const &y) {
    return dot(x, y) * inversesqrt(dot(x, x) * dot(y, y));
  }

  /// Normalize parameters and returns the dot product of x and y.
  /// @see gtx_normalize_dot
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER T fastNormalizeDot(qua<T, Q> const &x, qua<T, Q> const &y) {
    return dot(x, y) * fastInverseSqrt(dot(x, x) * dot(y, y));
  }

  /* matrix_extensions.hpp unification */

  /// @private
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER qua<T, Q> inverseTransform(qua<T, Q> const &q) {
    return inverse(q);
  }

  /// @private
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<3, T, Q> extractScale(qua<T, Q> const &q) {
    ((void)(q));
    return vec<3, T, Q>(T(1));
  }

  /// @private
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER bool hasUniformScale(qua<T, Q> const &q, T eps = epsilon<T>()) {
    ((void)(q));
    ((void)(eps));
    return true;
  }

  /// @private
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR vec<3, T, Q> transformPos(qua<T, Q> const &q, vec<3, T, Q> const &v) {
    return operator*(q, v);
  }

  /// @private
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR vec<3, T, Q> transformDir(qua<T, Q> const &q, vec<3, T, Q> const &v) {
    return q * v;
  }

  /*
  ** {======================================================
  ** Fixes
  ** =======================================================
  */

  /// <summary>
  /// @private
  /// @GLMFix: genTypeTrait glm::qualifier support
  /// </summary>
  namespace detail {
    template<typename T, glm::qualifier Q>
    struct genTypeTrait<qua<T, Q>> {
      static const genTypeEnum GENTYPE = GENTYPE_QUAT;
    };
  }

  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER T _angle(qua<T, Q> const &q) {
    const T n = length(vec<3, T, Q>(q.x, q.y, q.z));
    return detail::approx_zero(n) ? T(0) : T(2) * atan2(n, abs(q.w));
  }

  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER T _angle(qua<T, Q> const &x, qua<T, Q> const &y) {
    return _angle(y * conjugate(x));
  }

  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER void _axisAngle(qua<T, Q> const &q, vec<3, T, Q> &out_axis, T &out_angle) {
    out_axis = axis(q);
    out_angle = angle(q);
  }

  /// <summary>
  /// Return the oriented angle between two quaternions based on a reference axis.
  /// @see gtx_vector_angle
  /// @private
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER T _orientedAngle(qua<T, Q> const &x, qua<T, Q> const &y, vec<3, T, Q> const &ref) {
    const qua<T, Q> rot = y * conjugate(x);
    return _angle(rot) * sign(dot(ref, axis(rot)));
  }

  /// <summary>
  /// @private
  /// @GLMFix: _slerp in vector_extensions
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER qua<T, Q> _slerp(qua<T, Q> const &x, qua<T, Q> const &y, T const &a) {
    return slerp(x, y, a);
  }

  template<typename T, typename S, qualifier Q>
  GLM_FUNC_QUALIFIER qua<T, Q> _slerp(qua<T, Q> const &x, qua<T, Q> const &y, T a, S k) {
    return slerp(x, y, a, k);
  }

  /* }====================================================== */

  /// @}
}
