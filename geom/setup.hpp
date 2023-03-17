/// @ref geom
/// @file geom/setup.hpp
///
/// @defgroup geom Geometry extensions
/// @brief Additional geometric structures not specified by GLSL specification

#pragma once

#include <limits>
#include <type_traits>

#include <glm/glm.hpp>
#include <glm/common.hpp>
#include <glm/gtx/compatibility.hpp>
#include <glm/gtx/vector_query.hpp>
#include <glm/ext/scalar_constants.hpp>
#include <glm/ext/scalar_relational.hpp>
#include <glm/ext/vector_relational.hpp>
#if GLM_CONFIG_ALIGNED_GENTYPES == GLM_ENABLE
  #include <glm/simd/platform.h>
#endif

#include "../scalar_extensions.hpp"

/* @COMPAT: GLM_CONFIG_DEFAULTED_DEFAULT_CTOR introduced in 0.9.9.9 */
#if !defined(GLM_CONFIG_DEFAULTED_DEFAULT_CTOR)
  #define GLM_CONFIG_DEFAULTED_DEFAULT_CTOR GLM_CONFIG_DEFAULTED_FUNCTIONS
#endif

/* @COMPAT: GLM_COMPILER_HIP introduced in 0.9.9.9 */
#if !defined(GLM_COMPILER_HIP)
  #define GLM_COMPILER_HIP 0
#endif

/* gtx/string_cast.hpp is available */
#if 0 && !((GLM_COMPILER & GLM_COMPILER_CUDA) || (GLM_COMPILER & GLM_COMPILER_HIP))
  #define GLM_GEOM_TOSTRING 1
  #include <glm/gtx/string_cast.hpp>
#else
  #define GLM_GEOM_TOSTRING 0
#endif

/* Runtime preconditions on geometric structures */
#if !defined(NDEBUG)
  #define GLM_GEOM_ASSUME(x, onError) \
    do {                              \
      if (!(x))                       \
        return (onError);             \
    } while (0)
#else
  #define GLM_GEOM_ASSUME(x, onError)
#endif

#define GLM_GEOM_BEGIN_NAMESPACE namespace glm {
#define GLM_GEOM_END_NAMESPACE }

/* Redefinition of GLM_FUNC_QUALIFIER for 'geom' */
#define GLM_GEOM_DECL static GLM_FUNC_DECL
#define GLM_GEOM_QUALIFIER static GLM_FUNC_QUALIFIER
#define GLM_GEOM_QUALIFIER_OUTLINE static GLM_NEVER_INLINE

/// @addtogroup geom
/// @{
GLM_GEOM_BEGIN_NAMESPACE
template<length_t L, typename T, qualifier Q = defaultp> struct AABB;
template<length_t L, typename T, qualifier Q = defaultp> struct Line;
template<length_t L, typename T, qualifier Q = defaultp> struct Ray;
template<length_t L, typename T, qualifier Q = defaultp> struct Segment;
template<length_t L, typename T, qualifier Q = defaultp> struct Triangle;
template<length_t L, typename T, qualifier Q = defaultp> struct Sphere;
template<length_t L, typename T, qualifier Q = defaultp> struct Disc;
template<length_t L, typename T, qualifier Q = defaultp> struct Capsule;
template<length_t L, typename T, qualifier Q = defaultp> struct Plane;
template<length_t L, typename T, qualifier Q = defaultp> struct Polygon;
template<            typename T, qualifier Q = defaultp> struct OBB;
GLM_GEOM_END_NAMESPACE
/// @}

/******************************************************************************
 * See Copyright Notice in g-truc/glm/blob/master/copying.txt
 ******************************************************************************/
