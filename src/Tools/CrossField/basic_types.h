//
// Created by su on 2022/3/13.
//

#ifndef CROSSFIELD__BASIC_TYPES_H_
#define CROSSFIELD__BASIC_TYPES_H_

#include <cstdint>
#include <vector>
#include <array>
#include <cmath>
#include <cfloat>
#include <unordered_map>
#include <unordered_set>

#include "basic_utils.h"

#if defined(HAVE_64BIT_SIZE_T)
#define FNV_PRIME 1099511628211UL
#define FNV_OFFSET_BASIS 14695981039346656037UL
#else
#define PRIME 16777619UL
#define OFFSET_BASIS 2166136261UL
#endif

#ifndef M_PI
#define M_PI 3.14159265358979
#endif // M_PI


using std::unordered_map;
using std::unordered_set;

constexpr double EPS = 1.e-14; /* to detect 0 */

using Vec2 = std::array<double, 2>;
using Vec3 = std::array<double, 3>;
using ID = uint32_t;
using sid = int64_t;
using ID2 = std::array<ID, 2>;
using ID3 = std::array<ID, 3>;
using ID4 = std::array<ID, 4>;
using SID2 = std::array<sid, 2>;
using SID3 = std::array<sid, 3>;
using SID4 = std::array<sid, 4>;

const ID NO_ID = (ID)-1;
constexpr sid NO_SID = (sid)NO_ID;
constexpr uint8_t NO_U8 = (uint8_t)-1;
constexpr size_t NO_SIZE_T = (size_t)-1;
constexpr uint32_t NO_U32 = (uint32_t)-1;

/* Vec3 math */
inline double Dot(const Vec3 &a, const Vec3 &b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
inline double Length2(const Vec3 &a) { return Dot(a, a); }
inline double Length(const Vec3 &a) { return sqrt(Length2(a)); }
inline Vec3 operator-(const Vec3 &a, const Vec3 &b) {
  return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}
inline Vec3 operator+(const Vec3 &a, const Vec3 &b) {
  return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}
inline Vec3 operator*(const double &a, const Vec3 &b) {
  return {a * b[0], a * b[1], a * b[2]};
}
inline Vec3 operator*(const Vec3 &a, const double &b) {
  return {a[0] * b, a[1] * b, a[2] * b};
}
inline Vec3 Cross(const Vec3 &a, const Vec3 &b) {
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
          a[0] * b[1] - a[1] * b[0]};
}

inline ID2 Sorted(ID v1, ID v2) {
  if (v1 < v2) {
    return {v1, v2};
  }
  else { return {v2, v1}; }
}

inline double MaxAbs(double x, double y, double z) {
  return std::max(std::abs(x), std::max(std::abs(y), std::abs(z)));
}

inline double MaxAbs(const Vec3 &a) { return MaxAbs(a[0], a[1], a[2]); }

inline void NormalizeFast(Vec3 &a) {
  a = 1. / std::sqrt(Length2(a)) * a;
} /* no check, not safe, not accurate */
inline void Normalize(Vec3 &a) {
  double amp = MaxAbs(a);
  if (amp == 0.) { return; }
  a = amp * a;
  NormalizeFast(a);
}

inline double Clamp(double x, double lower, double upper) {
  return std::min(upper,
                  std::max(x, lower));
}
inline double AngleNVectors(Vec3 a, Vec3 b) {
  return acos(Clamp(Dot(a, b),
                    -1.,
                    1.));
}
inline double AngleVectors(Vec3 a, Vec3 b) {
  Normalize(a);
  Normalize(b);
  return AngleNVectors(a, b);
}
inline double AngleVectorsAlreadyNormalized(const Vec3 &a,
                                            const Vec3 &b) {
  return acos(Clamp(Dot(a, b),
                    -1.,
                    1.));
}

inline double Angle3Vertices(const Vec3 p1, const Vec3 p2, const Vec3 p3) {
  Vec3 a = p1 - p2;
  Vec3 b = p3 - p2;
  Vec3 c = Cross(a, b);
  double sinA = Length(c);
  double cosA = Dot(a, b);
  return std::atan2(sinA, cosA);
}

inline double Angle2Pi(const Vec3 &u, const Vec3 &v, const Vec3 &n) {
  const double dp = Dot(u, v);
  const double tp = Dot(n, Cross(u, v));
  double angle = atan2(tp, dp);
  if (angle < 0) angle += 2. * M_PI;
  return angle;
}

/* Vec2 math */
inline double Dot(const Vec2 &a, const Vec2 &b) {
  return a[0] * b[0] + a[1] * b[1];
}
inline double Length2(const Vec2 &a) { return Dot(a, a); }
inline double Length(const Vec2 &a) { return sqrt(Length2(a)); }

/* other */


inline Vec3 TriangleNormal(const Vec3 &p0, const Vec3 &p1, const Vec3 &p2) {
  Vec3 N = Cross(p2 - p0, p1 - p0);
  if (MaxAbs(N) == 0.) return {{0., 0., 0.}};
  Normalize(N);
  return N;
}

inline double TriangleArea(const Vec3 &p0, const Vec3 &p1, const Vec3 &p2) {
  Vec3 N = Cross(p2 - p0, p1 - p0);
  return Length(N) / 2.;
}

inline Vec3 TriangleCenter(const Vec3 &p0, const Vec3 &p1, const Vec3 &p2) {
  return 1. / 3. * (p0 + p1 + p2);
}

inline double CurveLength(const std::vector<Vec3> &pts) {
  if (pts.size() <= 1) return 0.;
  double L = 0.;
  for(int i=0;i< (int)pts.size() - 1;++i) L += Length(pts[i + 1] - pts[i]);
  return L;
}

inline double BboxDiag(const std::vector<Vec3> &points) {
  Vec3 mi = {DBL_MAX, DBL_MAX, DBL_MAX};
  Vec3 ma = {-DBL_MAX, -DBL_MAX, -DBL_MAX};
  for(auto i:points) {
    for(int d=0;d<3;++d) {
      mi[d] = std::min(i[d], mi[d]);
      ma[d] = std::max(i[d], ma[d]);
    }
  }
  return Length(ma - mi);
}

template<class T>
struct ArrayHash {
  size_t operator()(const T &key) const {
    int len = key.size();
    int s = sizeof(key[0]);
    size_t uc = 0xff;
    size_t hash = OFFSET_BASIS;
    for (int n = 0; n < len; ++n) {
      //const unsigned char* p = static_cast<const unsigned char*>(key[n]);
      auto p = static_cast<size_t>(key[n]);
      for (int k = s; k--;) {
        hash = (hash ^ (p & uc)) * PRIME;
        p >>= 8;
      }
    }
    return hash;
  }
};





#endif //CROSSFIELD__BASIC_TYPES_H_
