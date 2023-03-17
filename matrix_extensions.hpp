/// @ref gtx_matrix_extensions
/// @file matrix_extensions.hpp
///
/// @defgroup gtx_matrix_extensions GLM_GTX_matrix_extensions
/// @ingroup gtx
///
/// GLM matrix extensions
///  1. API unifying functions, usually to support mat<3, 4> and mat<4, 3>
///  2. Functions emulated/ported from other popular vector-math libraries
///
/// @see core_func_common
/// @see ext_matrix_common

#pragma once
#if !defined(GLM_ENABLE_EXPERIMENTAL)
  #define GLM_ENABLE_EXPERIMENTAL
#endif

#include <limits>

#include <glm/glm.hpp>
#include <glm/gtx/euler_angles.hpp>
#include <glm/ext/matrix_clip_space.hpp>
#include <glm/ext/matrix_relational.hpp>
#include <glm/ext/matrix_transform.hpp>
#include <glm/ext/scalar_constants.hpp>
#if GLM_MESSAGES == GLM_ENABLE && !defined(GLM_EXT_INCLUDED)
  #pragma message("GLM: GLM_GTX_matrix_ext extension included")
#endif

namespace glm {

  /// @addtogroup gtx_matrix_extensions
  /// @{

  /// <summary>
  /// all(equal(x, y)) shorthand
  /// @private @see ext_matrix_relational
  /// </summary>
  template<length_t C, length_t R, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR bool all_equal(mat<C, R, T, Q> const &a, mat<C, R, T, Q> const &b) {
    return all(equal(a, b));
  }

  /// <summary>
  /// all(equal(x, y, eps)) shorthand
  /// @private @see ext_matrix_relational
  /// </summary>
  template<length_t C, length_t R, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR bool all_equal(mat<C, R, T, Q> const &a, mat<C, R, T, Q> const &b, T eps) {
    return all(equal(a, b, eps));
  }

  /// <summary>
  /// all(equal(x, y, eps)) shorthand
  /// @private @see ext_matrix_relational
  /// </summary>
  template<length_t C, length_t R, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR bool all_equal(mat<C, R, T, Q> const &a, mat<C, R, T, Q> const &b, vec<C, T, Q> const &eps) {
    return all(equal(a, b, eps));
  }

  /// <summary>
  /// all(equal(x, y, MaxULPs)) shorthand
  /// @private @see ext_matrix_relational
  /// </summary>
  template<length_t C, length_t R, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR bool all_equal(mat<C, R, T, Q> const &a, mat<C, R, T, Q> const &b, int MaxULPs) {
    return all(equal(a, b, MaxULPs));
  }

  /// <summary>
  /// all(equal(x, y, MaxULPs)) shorthand
  /// @private @see ext_matrix_relational
  /// </summary>
  template<length_t C, length_t R, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR bool all_equal(mat<C, R, T, Q> const &a, mat<C, R, T, Q> const &b, vec<C, int, Q> const &MaxULPs) {
    return all(equal(a, b, MaxULPs));
  }

  /// <summary>
  /// any(notEqual(x, y)) shorthand
  /// @private @see ext_matrix_relational
  /// </summary>
  template<length_t C, length_t R, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR bool any_notEqual(mat<C, R, T, Q> const &a, mat<C, R, T, Q> const &b) {
    return any(notEqual(a, b));
  }

  /// <summary>
  /// any(notEqual(x, y, eps)) shorthand
  /// @private @see ext_matrix_relational
  /// </summary>
  template<length_t C, length_t R, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR bool any_notEqual(mat<C, R, T, Q> const &a, mat<C, R, T, Q> const &b, T eps) {
    return any(notEqual(a, b, eps));
  }

  /// <summary>
  /// any(notEqual(x, y, eps)) shorthand
  /// @private @see ext_matrix_relational
  /// </summary>
  template<length_t C, length_t R, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR bool any_notEqual(mat<C, R, T, Q> const &a, mat<C, R, T, Q> const &b, vec<C, T, Q> const &eps) {
    return any(notEqual(a, b, eps));
  }

  /// <summary>
  /// any(notEqual(x, y, MaxULPs)) shorthand
  /// @private @see ext_matrix_relational
  /// </summary>
  template<length_t C, length_t R, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR bool any_notEqual(mat<C, R, T, Q> const &a, mat<C, R, T, Q> const &b, int MaxULPs) {
    return any(notEqual(a, b, MaxULPs));
  }

  /// <summary>
  /// any(notEqual(x, y, MaxULPs)) shorthand
  /// @private @see ext_matrix_relational
  /// </summary>
  template<length_t C, length_t R, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR bool any_notEqual(mat<C, R, T, Q> const &a, mat<C, R, T, Q> const &b, vec<C, int, Q> const &MaxULPs) {
    return any(notEqual(a, b, MaxULPs));
  }

  /// <summary>
  /// Return the diagonal-vector of the given matrix.
  /// @see gtx_matrix_operation
  /// </summary>
  template<length_t C, length_t R, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<(C < R ? C : R), T, Q> diagonal(mat<C, R, T, Q> const &m) {
    vec<(C < R ? C : R), T, Q> result(T(0));
    for (length_t i = 0; i < (C < R ? C : R); ++i)
      result[i] = m[i][i];
    return result;
  }

  /// <summary>
  /// Test if the matrix has an inverse (up to a given epsilon).
  /// @see core_func_matrix
  /// </summary>
  template<length_t C, length_t R, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER bool invertible(mat<C, R, T, Q> const &m, T eps = epsilon<T>()) {
    GLM_STATIC_ASSERT(C == R, "Symmetric Matrices");
    return !detail::approx_zero(determinant(m), eps);
  }

  /// <summary>
  /// Transforms the given point vector by this matrix M, i.e, computes M * (x, y, z, 1).
  /// This function does not divide by w, or output it, so it cannot have a projection.
  /// @see ext_matrix_transform
  /// </summary>
  template<length_t C, length_t R, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<3, T, Q> transformPos(mat<C, R, T, Q> const &m, vec<3, T, Q> const &v) {
    GLM_STATIC_ASSERT(C >= 4 && R >= 3, "invalid position transform");
    const typename mat<C, R, T, Q>::col_type &result = m * vec<4, T, Q>(v, T(1));
    return vec<3, T, Q>(result.x, result.y, result.z);
  }

  /// <summary>
  /// Functional operator*(matrix, vector) wrapper.
  /// @see ext_matrix_transform
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR vec<3, T, Q> transformPos(mat<3, 3, T, Q> const &m, vec<3, T, Q> const &v) {
    return operator*(m, v);
  }

  /// <summary>
  /// Functional operator*(matrix, vector) wrapper.
  /// @see ext_matrix_transforms
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR vec<3, T, Q> transformPos(mat<3, 4, T, Q> const &m, vec<3, T, Q> const &v) {
    return operator*(m, v);
  }

  /// <summary>
  /// Transforms a position by a mat4x4 with a perspective divide.
  /// @see ext_matrix_transform
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR vec<3, T, Q> transformPosPerspective(mat<4, 4, T, Q> const &m, vec<3, T, Q> const &v) {
    const vec<3, T, Q> res = transformPos(m, v);
    const T w = m[0].w * v.x + m[1].w * v.y + m[2].w * v.z + m[3].w;
    return res * (T(1) / w);
  }

  /// <summary>
  /// Transforms the given normal by this matrix m, i.e, computes M * (x, y, z, 0).
  /// This function does not divide by w or output it and cannot have a projection.
  /// @see ext_matrix_transform
  /// </summary>
  template<length_t C, length_t R, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER GLM_CONSTEXPR vec<3, T, Q> transformDir(mat<C, R, T, Q> const &m, vec<3, T, Q> const &v) {
    GLM_STATIC_ASSERT(C >= 3 && R >= 3, "invalid direction transform");
    const typename mat<C, R, T, Q>::col_type result = m * vec<4, T, Q>(v, T(0));
    return vec<3, T, Q>(result.x, result.y, result.z);
  }

  /// <summary>
  /// Create a translation, rotation, and scaling matrix.
  /// @see ext_matrix_transform
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER mat<4, 4, T, Q> trs(vec<3, T, Q> const &translation, qua<T, Q> const &rotation, vec<3, T, Q> const &scale) {
    const mat<3, 3, T, Q> r = toMat3(rotation);
    return mat<4, 4, T, Q>(
      vec<4, T, Q>(r[0] * scale.x, T(0)),
      vec<4, T, Q>(r[1] * scale.y, T(0)),
      vec<4, T, Q>(r[2] * scale.z, T(0)),
      vec<4, T, Q>(translation, T(1))
    );
  }

  /// <summary>
  /// Return the scaling components of the given mat3x3.
  /// @see gtx_matrix_query
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<3, T, Q> extractScale(mat<3, 3, T, Q> const &m) {
    return vec<3, T, Q>(length(m[0]), length(m[1]), length(m[2]));
  }

  /// <summary>
  /// Return the scaling components of the given mat4x3.
  /// @see gtx_matrix_query
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<3, T, Q> extractScale(mat<4, 3, T, Q> const &m) {
    return vec<3, T, Q>(length(m[0]), length(m[1]), length(m[2]));
  }

  /// <summary>
  /// Return the scaling components of the given mat3x4.
  /// @see gtx_matrix_query
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<3, T, Q> extractScale(mat<3, 4, T, Q> const &m) {
    return vec<3, T, Q>(
      length(vec<3, T, Q>(m[0].x, m[0].y, m[0].z)),
      length(vec<3, T, Q>(m[1].x, m[1].y, m[1].z)),
      length(vec<3, T, Q>(m[2].x, m[2].y, m[2].z))
    );
  }

  /// <summary>
  /// Return the scaling components of the given mat4x4.
  /// @see gtx_matrix_query
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<3, T, Q> extractScale(mat<4, 4, T, Q> const &m) {
    return vec<3, T, Q>(
      length(vec<3, T, Q>(m[0].x, m[0].y, m[0].z)),
      length(vec<3, T, Q>(m[1].x, m[1].y, m[1].z)),
      length(vec<3, T, Q>(m[2].x, m[2].y, m[2].z))
    );
  }

  /// <summary>
  /// Returns true if the matrix contains only uniform scaling (up to a given epsilon).
  /// @see gtx_matrix_query
  /// </summary>
  template<length_t C, length_t R, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER bool hasUniformScale(mat<C, R, T, Q> const &m, T eps = epsilon<T>()) {
    const vec<3, T, Q> scale = extractScale(m);
    return epsilonEqual(scale.x, scale.y, eps) && epsilonEqual(scale.x, scale.z, eps);
  }

  /// <summary>
  /// Returns true if the matrix contains a projection, i.e., the last row
  /// differs (up to an epsilon) of [0, 0, 0, 1].
  /// @see ext_matrix_projection
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER bool containsProjection(mat<4, 4, T, Q> const &m, T eps = epsilon<T>()) {
    const vec<4, T, Q> v = row(m, 3);
    return all(epsilonEqual(v, vec<4, T, Q>(T(0), T(0), T(0), T(1)), eps));
  }

  /// <summary>
  /// Create a matrix that mirrors the given plane normal.
  /// @see ext_matrix_projection
  /// </summary>
  template<length_t C, length_t R, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER mat<C, R, T, Q> planeMirror(vec<3, T, Q> const &n, T d = T(0)) {
    GLM_STATIC_ASSERT(C >= 3 && R >= 3, "invalid affine plane mirror");
    GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'planeMirror' only accept floating-point inputs");
    mat<C, R, T, Q> m(T(0));
    m[0].x = T(1) - T(2) * n.x * n.x;
    m[0].y = T(-2) * n.x * n.y;
    m[0].z = T(-2) * n.x * n.z;
    m[1].x = T(-2) * n.y * n.x;
    m[1].y = T(1) - T(2) * n.y * n.y;
    m[1].z = T(-2) * n.y * n.z;
    m[2].x = T(-2) * n.z * n.x;
    m[2].y = T(-2) * n.y * n.z;
    m[2].z = T(1) - T(2) * n.z * n.z;
    GLM_IF_CONSTEXPR(C > 3) {
      m[3].x = T(2) * d * n.x;
      m[3].y = T(2) * d * n.y;
      m[3].z = T(2) * d * n.z;
    }
    return m;
  }

  /// <summary>
  /// Create an affine transformation matrix that projects orthographically onto the provided plane.
  /// @see ext_matrix_projection
  /// </summary>
  template<length_t C, length_t R, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER mat<C, R, T, Q> orthoProjection(vec<3, T, Q> const &n, T d = T(0)) {
    GLM_STATIC_ASSERT(C >= 3 && R >= 3, "invalid affine plane projection");
    GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'orthoProjection' only accept floating-point inputs");
    mat<C, R, T, Q> m(T(0));
    m[0].x = T(1) - n.x * n.x;
    m[0].y = -n.x * n.y;
    m[0].z = -n.x * n.z;
    m[1].x = -n.y * n.x;
    m[1].y = T(1) - n.y * n.y;
    m[1].z = -n.y * n.z;
    m[2].x = -n.z * n.x;
    m[2].y = -n.y * n.z;
    m[2].z = T(1) - n.z * n.z;
    GLM_IF_CONSTEXPR(C > 3) {
      m[3].x = d * n.x;
      m[3].y = d * n.y;
      m[3].z = d * n.z;
    }
    return m;
  }

  /// <summary>
  /// Mouse coordinates must be scaled to: [-1, 1]
  /// @see ext_matrix_projection
  /// </summary>
  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER vec<3, T, Q> rayPicking(vec<3, T, Q> const &camForward, vec<3, T, Q> const &camUp, T fov, T aspectRatio, T zNear, T zFar, T mouseX, T mouseY) {
    const mat<4, 4, T, Q> proj = perspective(fov, aspectRatio, zNear, zFar);
    const mat<4, 4, T, Q> view = lookAt(vec<3, T, Q>(T(0)), camForward, camUp);
    const mat<4, 4, T, Q> invpv = inverse(proj * view);
    const vec<4, T, Q> screenPos = vec<4, T, Q>(mouseX, -mouseY, T(1), T(1));
    return normalize(vec<3, T, Q>(invpv * screenPos));  // Direction of the ray; originating at the camera position.
  }

  /// <summary>
  /// Returns the right-handed rotation matrix for a given forward and up-vector,
  /// assuming the vectors are normalized and not collinear.
  /// @see gtc_quaternion
  /// </summary>
  template<typename T, qualifier Q, length_t C = 4, length_t R = 4>
  GLM_FUNC_QUALIFIER mat<C, R, T, Q> lookRotationRH(vec<3, T, Q> const &fwd, vec<3, T, Q> const &up) {
    GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'lookRotationRH' only accept floating-point inputs");
    const vec<3, T, Q> f = -fwd;
    const vec<3, T, Q> s(normalize(cross(up, f)));
    const vec<3, T, Q> u(cross(f, s));
    mat<C, R, T, Q> Result(1);
    Result[0].x = s.x;
    Result[1].x = s.y;
    Result[2].x = s.z;
    Result[0].y = u.x;
    Result[1].y = u.y;
    Result[2].y = u.z;
    Result[0].z = f.x;
    Result[1].z = f.y;
    Result[2].z = f.z;
    return Result;
  }

  /// <summary>
  /// Create a left-handed rotation matrix for a given forward and up-vector,
  /// assuming the vectors are normalized and not collinear.
  /// @see gtc_quaternion
  /// </summary>
  template<typename T, qualifier Q, length_t C = 4, length_t R = 4>
  GLM_FUNC_QUALIFIER mat<C, R, T, Q> lookRotationLH(vec<3, T, Q> const &fwd, vec<3, T, Q> const &up) {
    GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'lookRotationLH' only accept floating-point inputs");
    const vec<3, T, Q> s(normalize(cross(up, fwd)));
    const vec<3, T, Q> u(cross(fwd, s));
    mat<C, R, T, Q> Result(1);
    Result[0].x = s.x;
    Result[1].x = s.y;
    Result[2].x = s.z;
    Result[0].y = u.x;
    Result[1].y = u.y;
    Result[2].y = u.z;
    Result[0].z = fwd.x;
    Result[1].z = fwd.y;
    Result[2].z = fwd.z;
    return Result;
  }

  /// <summary>
  /// Build a look at matrix based on the default handedness.
  /// @see gtc_quaternion
  /// </summary>
  template<typename T, qualifier Q, length_t C = 4, length_t R = 4>
  GLM_FUNC_QUALIFIER mat<C, R, T, Q> lookRotation(vec<3, T, Q> const &fwd, vec<3, T, Q> const &up) {
#if (GLM_CONFIG_CLIP_CONTROL & GLM_CLIP_CONTROL_LH_BIT)
    return lookRotationLH<T, Q, C, R>(fwd, up);
#else
    return lookRotationRH<T, Q, C, R>(fwd, up);
#endif
  }

  /// <summary>
  /// Creates a right-handed spherical billboard that rotates around a specified 'object' position.
  /// @see ext_matrix_transform
  /// </summary>
  template<typename T, qualifier Q, length_t C = 4, length_t R = 4>
  GLM_FUNC_QUALIFIER mat<C, R, T, Q> billboardRH(vec<3, T, Q> const &object, vec<3, T, Q> const &camPos, vec<3, T, Q> const &camUp, vec<3, T, Q> const &camFwd) {
    GLM_STATIC_ASSERT(C >= 3 && R >= 3 && C == R, "invalid billboard matrix");
    GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'billboardRH' only accept floating-point inputs");

    vec<3, T, Q> look = object - camPos;
    const T lengthSq = length2(look);
    if (detail::approx_zero(lengthSq))
      look = -camFwd;
    else
      look *= (T(1) / sqrt(lengthSq));

    const vec<3, T, Q> right = normalize(cross(camUp, look));
    const vec<3, T, Q> up = cross(look, right);
    mat<C, R, T, Q> Result(1);
    Result[0].x = right.x;
    Result[0].y = right.y;
    Result[0].z = right.z;
    Result[1].x = up.x;
    Result[1].y = up.y;
    Result[1].z = up.z;
    Result[2].x = look.x;
    Result[2].y = look.y;
    Result[2].z = look.z;
    GLM_IF_CONSTEXPR(R > 3) {
      Result[3][0] = object.x;
      Result[3][1] = object.y;
      Result[3][2] = object.z;
    }
    return Result;
  }

  /// <summary>
  /// Creates a left-handed spherical billboard that rotates around a specified 'object' position.
  /// @see ext_matrix_transform
  /// </summary>
  template<typename T, qualifier Q, length_t C = 4, length_t R = 4>
  GLM_FUNC_QUALIFIER mat<C, R, T, Q> billboardLH(vec<3, T, Q> const &object, vec<3, T, Q> const &camPos, vec<3, T, Q> const &camUp, vec<3, T, Q> const &camFwd) {
    GLM_STATIC_ASSERT(C >= 3 && R >= 3 && C == R, "invalid billboard matrix");
    GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'billboardLH' only accept floating-point inputs");

    vec<3, T, Q> look = camPos - object;
    const T lengthSq = length2(look);
    if (detail::approx_zero(lengthSq))
      look = -camFwd;
    else
      look *= (T(1) / sqrt(lengthSq));

    const vec<3, T, Q> right = normalize(cross(camUp, look));
    const vec<3, T, Q> up = cross(look, right);
    mat<C, R, T, Q> Result(1);
    Result[0].x = right.x;
    Result[0].y = right.y;
    Result[0].z = right.z;
    Result[1].x = up.x;
    Result[1].y = up.y;
    Result[1].z = up.z;
    Result[2].x = look.x;
    Result[2].y = look.y;
    Result[2].z = look.z;
    GLM_IF_CONSTEXPR(R > 3) {
      Result[3][0] = object.x;
      Result[3][1] = object.y;
      Result[3][2] = object.z;
    }
    return Result;
  }

  /// <summary>
  /// Create a spherical billboard that rotates around a specified 'object' position; uses default handedness.
  /// @see ext_matrix_transform
  /// </summary>
  template<typename T, qualifier Q, length_t C = 4, length_t R = 4>
  GLM_FUNC_QUALIFIER mat<C, R, T, Q> billboard(vec<3, T, Q> const &object, vec<3, T, Q> const &pos, vec<3, T, Q> const &up, vec<3, T, Q> const &forward) {
#if (GLM_CONFIG_CLIP_CONTROL & GLM_CLIP_CONTROL_LH_BIT)
    return billboardLH<T, Q, C, R>(object, pos, up, forward);
#else
    return billboardRH<T, Q, C, R>(object, pos, up, forward);
#endif
  }

  /* EulerAngle Extraction for all matrices with rotation parts */

  /// Extracts the (X * Y * Z) Euler angles from the rotation matrix M
  /// @see gtx_euler_angles
  template<length_t C, length_t R, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER void extractEulerAngleXYZ(mat<C, R, T, Q> const &M, T &t1, T &t2, T &t3) {
    GLM_STATIC_ASSERT(C >= 3 && R >= 3, "invalid extraction dimensions");
    GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'extractEulerAngle' only accept floating-point inputs");
    const T T1 = atan2(M[2][1], M[2][2]);
    const T C2 = sqrt(M[0][0] * M[0][0] + M[1][0] * M[1][0]);
    const T T2 = atan2(-M[2][0], C2);
    const T S1 = sin(T1);
    const T C1 = cos(T1);
    const T T3 = atan2(S1 * M[0][2] - C1 * M[0][1], C1 * M[1][1] - S1 * M[1][2]);
    t1 = -T1;
    t2 = -T2;
    t3 = -T3;
  }

  /// Extracts the (Y * X * Z) Euler angles from the rotation matrix M
  /// @see gtx_euler_angles
  template<length_t C, length_t R, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER void extractEulerAngleYXZ(mat<C, R, T, Q> const &M, T &t1, T &t2, T &t3) {
    GLM_STATIC_ASSERT(C >= 3 && R >= 3, "invalid extraction dimensions");
    GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'extractEulerAngle' only accept floating-point inputs");
    const T T1 = atan2(M[2][0], M[2][2]);
    const T C2 = sqrt(M[0][1] * M[0][1] + M[1][1] * M[1][1]);
    const T T2 = atan2(-M[2][1], C2);
    const T S1 = sin(T1);
    const T C1 = cos(T1);
    const T T3 = atan2(S1 * M[1][2] - C1 * M[1][0], C1 * M[0][0] - S1 * M[0][2]);
    t1 = T1;
    t2 = T2;
    t3 = T3;
  }

  /// Extracts the (X * Z * X) Euler angles from the rotation matrix M
  /// @see gtx_euler_angles
  template<length_t C, length_t R, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER void extractEulerAngleXZX(mat<C, R, T, Q> const &M, T &t1, T &t2, T &t3) {
    GLM_STATIC_ASSERT(C >= 3 && R >= 3, "invalid extraction dimensions");
    GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'extractEulerAngle' only accept floating-point inputs");
    const T T1 = atan2(M[0][2], M[0][1]);
    const T S2 = sqrt(M[1][0] * M[1][0] + M[2][0] * M[2][0]);
    const T T2 = atan2(S2, M[0][0]);
    const T S1 = sin(T1);
    const T C1 = cos(T1);
    const T T3 = atan2(C1 * M[1][2] - S1 * M[1][1], C1 * M[2][2] - S1 * M[2][1]);
    t1 = T1;
    t2 = T2;
    t3 = T3;
  }

  /// Extracts the (X * Y * X) Euler angles from the rotation matrix M
  /// @see gtx_euler_angles
  template<length_t C, length_t R, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER void extractEulerAngleXYX(mat<C, R, T, Q> const &M, T &t1, T &t2, T &t3) {
    GLM_STATIC_ASSERT(C >= 3 && R >= 3, "invalid extraction dimensions");
    GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'extractEulerAngle' only accept floating-point inputs");
    const T T1 = atan2(M[0][1], -M[0][2]);
    const T S2 = sqrt(M[1][0] * M[1][0] + M[2][0] * M[2][0]);
    const T T2 = atan2(S2, M[0][0]);
    const T S1 = sin(T1);
    const T C1 = cos(T1);
    const T T3 = atan2(-C1 * M[2][1] - S1 * M[2][2], C1 * M[1][1] + S1 * M[1][2]);
    t1 = T1;
    t2 = T2;
    t3 = T3;
  }

  /// Extracts the (Y * X * Y) Euler angles from the rotation matrix M
  /// @see gtx_euler_angles
  template<length_t C, length_t R, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER void extractEulerAngleYXY(mat<C, R, T, Q> const &M, T &t1, T &t2, T &t3) {
    GLM_STATIC_ASSERT(C >= 3 && R >= 3, "invalid extraction dimensions");
    GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'extractEulerAngle' only accept floating-point inputs");
    const T T1 = atan2(M[1][0], M[1][2]);
    const T S2 = sqrt(M[0][1] * M[0][1] + M[2][1] * M[2][1]);
    const T T2 = atan2(S2, M[1][1]);
    const T S1 = sin(T1);
    const T C1 = cos(T1);
    const T T3 = atan2(C1 * M[2][0] - S1 * M[2][2], C1 * M[0][0] - S1 * M[0][2]);
    t1 = T1;
    t2 = T2;
    t3 = T3;
  }

  /// Extracts the (Y * Z * Y) Euler angles from the rotation matrix M
  /// @see gtx_euler_angles
  template<length_t C, length_t R, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER void extractEulerAngleYZY(mat<C, R, T, Q> const &M, T &t1, T &t2, T &t3) {
    GLM_STATIC_ASSERT(C >= 3 && R >= 3, "invalid extraction dimensions");
    GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'extractEulerAngle' only accept floating-point inputs");
    const T T1 = atan2(M[1][2], -M[1][0]);
    const T S2 = sqrt(M[0][1] * M[0][1] + M[2][1] * M[2][1]);
    const T T2 = atan2(S2, M[1][1]);
    const T S1 = sin(T1);
    const T C1 = cos(T1);
    const T T3 = atan2(-S1 * M[0][0] - C1 * M[0][2], S1 * M[2][0] + C1 * M[2][2]);
    t1 = T1;
    t2 = T2;
    t3 = T3;
  }

  /// Extracts the (Z * Y * Z) Euler angles from the rotation matrix M
  /// @see gtx_euler_angles
  template<length_t C, length_t R, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER void extractEulerAngleZYZ(mat<C, R, T, Q> const &M, T &t1, T &t2, T &t3) {
    GLM_STATIC_ASSERT(C >= 3 && R >= 3, "invalid extraction dimensions");
    GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'extractEulerAngle' only accept floating-point inputs");
    const T T1 = atan2(M[2][1], M[2][0]);
    const T S2 = sqrt(M[0][2] * M[0][2] + M[1][2] * M[1][2]);
    const T T2 = atan2(S2, M[2][2]);
    const T S1 = sin(T1);
    const T C1 = cos(T1);
    const T T3 = atan2(C1 * M[0][1] - S1 * M[0][0], C1 * M[1][1] - S1 * M[1][0]);
    t1 = T1;
    t2 = T2;
    t3 = T3;
  }

  /// Extracts the (Z * X * Z) Euler angles from the rotation matrix M
  /// @see gtx_euler_angles
  template<length_t C, length_t R, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER void extractEulerAngleZXZ(mat<C, R, T, Q> const &M, T &t1, T &t2, T &t3) {
    GLM_STATIC_ASSERT(C >= 3 && R >= 3, "invalid extraction dimensions");
    GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'extractEulerAngle' only accept floating-point inputs");
    const T T1 = atan2(M[2][0], -M[2][1]);
    const T S2 = sqrt(M[0][2] * M[0][2] + M[1][2] * M[1][2]);
    const T T2 = atan2(S2, M[2][2]);
    const T S1 = sin(T1);
    const T C1 = cos(T1);
    const T T3 = atan2(-C1 * M[1][0] - S1 * M[1][1], C1 * M[0][0] + S1 * M[0][1]);
    t1 = T1;
    t2 = T2;
    t3 = T3;
  }

  /// Extracts the (X * Z * Y) Euler angles from the rotation matrix M
  /// @see gtx_euler_angles
  template<length_t C, length_t R, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER void extractEulerAngleXZY(mat<C, R, T, Q> const &M, T &t1, T &t2, T &t3) {
    GLM_STATIC_ASSERT(C >= 3 && R >= 3, "invalid extraction dimensions");
    GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'extractEulerAngle' only accept floating-point inputs");
    const T T1 = atan2(M[1][2], M[1][1]);
    const T C2 = sqrt(M[0][0] * M[0][0] + M[2][0] * M[2][0]);
    const T T2 = atan2(-M[1][0], C2);
    const T S1 = sin(T1);
    const T C1 = cos(T1);
    const T T3 = atan2(S1 * M[0][1] - C1 * M[0][2], C1 * M[2][2] - S1 * M[2][1]);
    t1 = T1;
    t2 = T2;
    t3 = T3;
  }

  /// Extracts the (Y * Z * X) Euler angles from the rotation matrix M
  /// @see gtx_euler_angles
  template<length_t C, length_t R, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER void extractEulerAngleYZX(mat<C, R, T, Q> const &M, T &t1, T &t2, T &t3) {
    GLM_STATIC_ASSERT(C >= 3 && R >= 3, "invalid extraction dimensions");
    GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'extractEulerAngle' only accept floating-point inputs");
    const T T1 = atan2(-M[0][2], M[0][0]);
    const T C2 = sqrt(M[1][1] * M[1][1] + M[2][1] * M[2][1]);
    const T T2 = atan2(M[0][1], C2);
    const T S1 = sin(T1);
    const T C1 = cos(T1);
    const T T3 = atan2(S1 * M[1][0] + C1 * M[1][2], S1 * M[2][0] + C1 * M[2][2]);
    t1 = T1;
    t2 = T2;
    t3 = T3;
  }

  /// Extracts the (Z * Y * X) Euler angles from the rotation matrix M
  /// @see gtx_euler_angles
  template<length_t C, length_t R, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER void extractEulerAngleZYX(mat<C, R, T, Q> const &M, T &t1, T &t2, T &t3) {
    GLM_STATIC_ASSERT(C >= 3 && R >= 3, "invalid extraction dimensions");
    GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'extractEulerAngle' only accept floating-point inputs");
    const T T1 = atan2(M[0][1], M[0][0]);
    const T C2 = sqrt(M[1][2] * M[1][2] + M[2][2] * M[2][2]);
    const T T2 = atan2(-M[0][2], C2);
    const T S1 = sin(T1);
    const T C1 = cos(T1);
    const T T3 = atan2(S1 * M[2][0] - C1 * M[2][1], C1 * M[1][1] - S1 * M[1][0]);
    t1 = T1;
    t2 = T2;
    t3 = T3;
  }

  /// Extracts the (Z * X * Y) Euler angles from the rotation matrix M
  /// @see gtx_euler_angles
  template<length_t C, length_t R, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER void extractEulerAngleZXY(mat<C, R, T, Q> const &M, T &t1, T &t2, T &t3) {
    GLM_STATIC_ASSERT(C >= 3 && R >= 3, "invalid extraction dimensions");
    GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'extractEulerAngle' only accept floating-point inputs");
    const T T1 = atan2(-M[1][0], M[1][1]);
    const T C2 = sqrt(M[0][2] * M[0][2] + M[2][2] * M[2][2]);
    const T T2 = atan2(M[1][2], C2);
    const T S1 = sin(T1);
    const T C1 = cos(T1);
    const T T3 = atan2(C1 * M[2][0] + S1 * M[2][1], C1 * M[0][0] + S1 * M[0][1]);
    t1 = T1;
    t2 = T2;
    t3 = T3;
  }

  /// <summary>
  /// Build a rotation matrix created from a normalized axis and an angle.
  ///
  /// @param m Input matrix multiplied by this rotation matrix.
  /// @param angle Rotation angle expressed in radians.
  /// @param v Rotation axis, must be normalized.
  /// @tparam T Value type used to build the matrix.
  /// @see gtx_rotate_normalized_axis
  /// </summary>
  template<length_t C, length_t R, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER mat<C, R, T, Q> rotateNormalizedAxis(mat<C, R, T, Q> const &m, T const &angle, vec<3, T, Q> const &v) {
    GLM_STATIC_ASSERT(C >= 3 && R >= 3, "invalid rotation matrix");
    GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'rotateNormalizedAxis' only accept floating-point inputs");

    const T a = angle;
    const T c = cos(a);
    const T s = sin(a);
    const vec<3, T, Q> axis(v);
    const vec<3, T, Q> temp((static_cast<T>(1) - c) * axis);

    mat<3, 3, T, Q> Rotate;
    Rotate[0].x = c + temp[0] * axis[0];
    Rotate[0].y = 0 + temp[0] * axis[1] + s * axis[2];
    Rotate[0].z = 0 + temp[0] * axis[2] - s * axis[1];
    Rotate[1].x = 0 + temp[1] * axis[0] - s * axis[2];
    Rotate[1].y = c + temp[1] * axis[1];
    Rotate[1].z = 0 + temp[1] * axis[2] + s * axis[0];
    Rotate[2].x = 0 + temp[2] * axis[0] + s * axis[1];
    Rotate[2].y = 0 + temp[2] * axis[1] - s * axis[0];
    Rotate[2].z = c + temp[2] * axis[2];

    mat<C, R, T, Q> Result;
    Result[0] = m[0] * Rotate[0].x + m[1] * Rotate[0].y + m[2] * Rotate[0].z;
    Result[1] = m[0] * Rotate[1].x + m[1] * Rotate[1].y + m[2] * Rotate[1].z;
    Result[2] = m[0] * Rotate[2].x + m[1] * Rotate[2].y + m[2] * Rotate[2].z;
    GLM_IF_CONSTEXPR(C > 3)
      Result[3] = m[3];
    return Result;
  }

  /// <summary>
  /// glm::inverse that assumes the last row is (0, 0, 0, 1)
  /// @see core_func_matrix
  /// </summary>
  template<length_t C, length_t R, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER mat<C, R, T, Q> inverseTransform(mat<C, R, T, Q> const &m) {
    GLM_STATIC_ASSERT(C >= 3 && R >= 3, "invalid extraction dimensions");
    const T OneOverDeterminant = static_cast<T>(1) / (
      + m[0][0] * (m[1][1] * m[2][2] - m[2][1] * m[1][2])
      - m[1][0] * (m[0][1] * m[2][2] - m[2][1] * m[0][2])
      + m[2][0] * (m[0][1] * m[1][2] - m[1][1] * m[0][2])
    );

    mat<C, R, T, Q> Inverse(T(0));
    Inverse[0][0] = +(m[1][1] * m[2][2] - m[2][1] * m[1][2]) * OneOverDeterminant;
    Inverse[1][0] = -(m[1][0] * m[2][2] - m[2][0] * m[1][2]) * OneOverDeterminant;
    Inverse[2][0] = +(m[1][0] * m[2][1] - m[2][0] * m[1][1]) * OneOverDeterminant;
    Inverse[0][1] = -(m[0][1] * m[2][2] - m[2][1] * m[0][2]) * OneOverDeterminant;
    Inverse[1][1] = +(m[0][0] * m[2][2] - m[2][0] * m[0][2]) * OneOverDeterminant;
    Inverse[2][1] = -(m[0][0] * m[2][1] - m[2][0] * m[0][1]) * OneOverDeterminant;
    Inverse[0][2] = +(m[0][1] * m[1][2] - m[1][1] * m[0][2]) * OneOverDeterminant;
    Inverse[1][2] = -(m[0][0] * m[1][2] - m[1][0] * m[0][2]) * OneOverDeterminant;
    Inverse[2][2] = +(m[0][0] * m[1][1] - m[1][0] * m[0][1]) * OneOverDeterminant;
    GLM_IF_CONSTEXPR(R > 3 && C > 3)
      Inverse[3][3] = T(1);
    return Inverse;
  }

  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER mat<3, 3, T, Q> inverseTransform(mat<3, 3, T, Q> const &m) {
    return inverse(m);
  }

  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER mat<2, 2, T, Q> inverseTransform(mat<2, 2, T, Q> const &m) {
    return inverse(m);
  }

  /*
  ** {======================================================
  ** Fixes
  ** =======================================================
  */

  /// <summary>
  /// @private
  /// @GLMFix: genTypeTrait qualifier support
  /// </summary>
  namespace detail {
    template<length_t C, length_t R, typename T, glm::qualifier Q>
    struct genTypeTrait<mat<C, R, T, Q>> {
      static const genTypeEnum GENTYPE = GENTYPE_MAT;
    };
  }

  /// <summary>
  /// @private
  /// @GLMFix: glm::scaleBias that ensures the matrix is initialized.
  /// </summary>
  template<typename T, qualifier Q = qualifier::defaultp>
  GLM_FUNC_QUALIFIER mat<4, 4, T, Q> _scaleBias(T scale, T bias) {
    mat<4, 4, T, Q> result(T(0));
    result[3] = vec<4, T, Q>(vec<3, T, Q>(bias), static_cast<T>(1));
    result[0][0] = scale;
    result[1][1] = scale;
    result[2][2] = scale;
    return result;
  }

  template<typename T, qualifier Q>
  GLM_FUNC_QUALIFIER mat<4, 4, T, Q> _scaleBias(mat<4, 4, T, Q> const &m, T scale, T bias) {
    return m * _scaleBias<T, Q>(scale, bias);
  }

  /// <summary>
  /// @private
  /// @GLMFix: Generalized implementation
  /// </summary>
  template<length_t C, length_t R, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER bool _isNull(mat<C, R, T, Q> const &m, T eps = epsilon<T>()) {
    bool result = true;
    for (length_t i = 0; i < C; ++i)
      result &= isNull(m[i], eps);  // @RemoveSqrt
    return result;
  }

  /// <summary>
  /// @private
  /// @GLMFix: Generalized implementation
  /// </summary>
  template<length_t C, length_t R, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER bool _isNormalized(mat<C, R, T, Q> const &m, T eps = epsilon<T>()) {
    bool result = true;
    for (length_t i = 0; i < C; ++i)
      result &= isNormalized(m[i], eps);  // @RemoveSqrt

    for (length_t i = 0; i < R; ++i) {
      typename mat<C, R, T, Q>::row_type v(T(0));
      for (length_t j = 0; j < C; ++j)
        v[j] = m[j][i];
      result &= isNormalized(v, eps);  // @RemoveSqrt
    }
    return result;
  }

  /// <summary>
  /// @private
  /// @GLMFix: Generalized implementation.
  /// </summary>
  template<length_t C, length_t R, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER bool _isOrthogonal(mat<C, R, T, Q> const &m, T const &epsilon) {
    bool result = true;
    for (length_t i = 0; i < C - 1; ++i) {
      for (length_t j = i + 1; j < C; ++j) {
        result &= areOrthogonal(m[i], m[j], epsilon);
      }
    }

    if (result) {
      const mat<C, R, T, Q> tmp = transpose(m);
      for (length_t i = 0; i < C - 1; ++i) {
        for (length_t j = i + 1; j < C; ++j) {
          result &= areOrthogonal(tmp[i], tmp[j], epsilon);
        }
      }
    }
    return result;
  }

  /// <summary>
  /// @private
  /// @GLMFix: Generalized implementation.
  /// </summary>
  template<length_t C, length_t R, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER bool _isIdentity(mat<C, R, T, Q> const &m, T eps = epsilon<T>()) {
    bool result = true;
    for (length_t i = 0; i < C; ++i) {
      for (length_t j = 0; j < R; ++j) {
        result &= ((i == j) ? (abs(m[i][j] - 1) <= eps) : (abs(m[i][j]) <= eps));
      }
    }
    return result;
  }

  /// <summary>
  /// @private
  /// @GLMFix: Generic axisAngle support for all rotation matrices; code from glm/gtx/matrix_interpolation.inl.
  /// </summary>
  template<length_t C, length_t R, typename T, qualifier Q>
  GLM_FUNC_QUALIFIER void _axisAngle(mat<C, R, T, Q> const &m, vec<3, T, Q> &axis, T &angle) {
    GLM_STATIC_ASSERT(C >= 3 && R >= 3, "invalid rotation matrix");
    const T epsilon = std::numeric_limits<T>::epsilon() * static_cast<T>(1e2);

    const bool nearSymmetrical = abs(m[1][0] - m[0][1]) < epsilon && abs(m[2][0] - m[0][2]) < epsilon && abs(m[2][1] - m[1][2]) < epsilon;
    if (nearSymmetrical) {
      const bool nearIdentity = abs(m[1][0] + m[0][1]) < epsilon && abs(m[2][0] + m[0][2]) < epsilon && abs(m[2][1] + m[1][2]) < epsilon && abs(m[0][0] + m[1][1] + m[2][2] - T(3.0)) < epsilon;
      if (nearIdentity) {
        angle = static_cast<T>(0.0);
        axis = vec<3, T, Q>(static_cast<T>(1.0), static_cast<T>(0.0), static_cast<T>(0.0));
        return;
      }

      angle = pi<T>();
      const T xx = (m[0][0] + static_cast<T>(1.0)) * static_cast<T>(0.5);
      const T yy = (m[1][1] + static_cast<T>(1.0)) * static_cast<T>(0.5);
      const T zz = (m[2][2] + static_cast<T>(1.0)) * static_cast<T>(0.5);
      const T xy = (m[1][0] + m[0][1]) * static_cast<T>(0.25);
      const T xz = (m[2][0] + m[0][2]) * static_cast<T>(0.25);
      const T yz = (m[2][1] + m[1][2]) * static_cast<T>(0.25);
      if ((xx > yy) && (xx > zz)) {
        if (xx < epsilon) {
          axis.x = static_cast<T>(0.0);
          axis.y = static_cast<T>(0.7071);
          axis.z = static_cast<T>(0.7071);
        }
        else {
          axis.x = sqrt(xx);
          axis.y = xy / axis.x;
          axis.z = xz / axis.x;
        }
      }
      else if (yy > zz) {
        if (yy < epsilon) {
          axis.x = static_cast<T>(0.7071);
          axis.y = static_cast<T>(0.0);
          axis.z = static_cast<T>(0.7071);
        }
        else {
          axis.y = sqrt(yy);
          axis.x = xy / axis.y;
          axis.z = yz / axis.y;
        }
      }
      else {
        if (zz < epsilon) {
          axis.x = static_cast<T>(0.7071);
          axis.y = static_cast<T>(0.7071);
          axis.z = static_cast<T>(0.0);
        }
        else {
          axis.z = sqrt(zz);
          axis.x = xz / axis.z;
          axis.y = yz / axis.z;
        }
      }
      return;
    }

    const T angleCos = (m[0][0] + m[1][1] + m[2][2] - static_cast<T>(1)) * static_cast<T>(0.5);
    if (angleCos >= static_cast<T>(1.0)) {
      angle = static_cast<T>(0.0);
    }
    else if (angleCos <= static_cast<T>(-1.0)) {
      angle = pi<T>();
    }
    else {
      angle = acos(angleCos);
    }
    axis = normalize(vec<3, T, Q>(m[1][2] - m[2][1], m[2][0] - m[0][2], m[0][1] - m[1][0]));
  }

  /* }====================================================== */

  /// @}
}
