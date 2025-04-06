#pragma once

#include <iostream>
#include <string>

#include "util/types.h"

namespace limap {

// Unique identifier for cameras.
typedef uint32_t camera_t;

// Unique identifier for images.
typedef uint32_t image_t;

// Index of a 2D point in an image.
typedef uint32_t point2d_t;

// Unique identifier for 3D points.
typedef uint64_t point3d_t;

// Index of a 2D line segment in an image.
// Note: In this project, "line" is short for "line segment".
//       Please not to be confused with "infinite line".
typedef uint32_t line2d_t;

// Unique identifier for 3D lines (3D line tracks).
typedef uint64_t line3d_t;

// Unique identifier for 2D line segments.
typedef std::pair<image_t, line2d_t> image_line2d_t;

const line3d_t kInvalidLine3dId = std::numeric_limits<line3d_t>::max();

////////////////////////////////////////////////////////////////////////////////
// debug tools
////////////////////////////////////////////////////////////////////////////////

// helper function for printing
template <typename T>
void PrintVarDetailsDebug(const char* file_name, const int line_no,
                          const char* var_name, const T& var_value,
                          const std::string& msg = std::string()) {
  std::cout << "[" << file_name << ":" << line_no << "] " << var_name << " = "
            << var_value << " " << msg << std::endl;
}

// macro definition for printing
#define PRINT(msg) \
  std::cout << "[" << __FILE__ << ":" << __LINE__ << "] " << msg << std::endl;
#define PRINT_VAR(var) PrintVarDetailsDebug(__FILE__, __LINE__, #var, var);
#define PRINT_VAR_MSG(var, msg) \
  PrintVarDetailsDebug(__FILE__, __LINE__, #var, var, msg);

#ifndef VERBOSE
#define VERBOSE 0
#endif

}  // namespace limap
