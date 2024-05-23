//
// Created by su on 2022/3/13.
//

#ifndef CROSSFIELD__BASIC_UTILS_H_
#define CROSSFIELD__BASIC_UTILS_H_

#include <cstddef>
#include <ctime>
#include <ostream>
#include <string>
#include <cstring>
#include <sstream>
#include <fstream>
#include <iostream>
#include <istream>
#include <vector>
#include <array>
#include <algorithm>

#define DBG(...) fprintf(stdout, "(DBG) %s:%i: ", __FILE__,__LINE__); show(std::cout, #__VA_ARGS__, __VA_ARGS__); fflush(stdout)

inline double Cpu() { return (double)(clock()) / (double)(CLOCKS_PER_SEC); }

inline void SwapBytes(char *array, int size, int n) {
  char *x = new char[size];
  for (int i = 0; i < n; i++) {
    char *a = &array[i * size];
    memcpy(x, a, size);
    for (int c = 0; c < size; c++) a[size - 1 - c] = x[c];
  }
  delete[] x;
}

inline std::vector<std::string> SplitFileName(const std::string &fileName) {
  // JFR DO NOT CHANGE TO std::vector<std::string> s(3), it segfaults while
  // destructor si called
  std::vector<std::string> s;
  s.resize(3);
  if (!fileName.empty()) {
    // returns [path, baseName, extension]
    int pos_dot = (int)fileName.find_last_of('.');
    int is_lash = (int)fileName.find_last_of("/\\");
    if (pos_dot == (int)std::string::npos) pos_dot = -1;
    if (is_lash == (int)std::string::npos) is_lash = -1;
    if (pos_dot > 0) s[2] = fileName.substr(pos_dot);
    if (is_lash > 0) s[0] = fileName.substr(0, is_lash + 1);
    s[1] =
        fileName
            .substr(s[0].size(), fileName.size() - s[0].size() - s[2].size());
  }
  return s;
}

using std::size_t;
template<class T1, class T2>
std::ostream &operator<<(std::ostream &os, const std::pair<T1, T2> &val) {
  return os << "(" << val.first << "," << val.second << ")";
}

template<class T, size_t N>
std::ostream &operator<<(std::ostream &os, const std::array<T, N> &values) {
  os << "(";
  for (size_t i = 0; i < values.size(); ++i) {
    os << values[i];
    if (i != values.size() - 1) {
      os << ", ";
    }
  }
  os << ")";
  return os;
}

template<class T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> &values) {
  os << "[";
  for (size_t i = 0; i < values.size(); ++i) {
    const T &x = values[i];
    os << x;
    if (i != values.size() - 1) {
      os << ", ";
    }
  }
  os << "]";
  return os;
}

static inline void SFormat(std::ostream &out, const char *s) {
  while (*s) {
    if (*s == '{' && *++s != '}')
      throw std::runtime_error("invalid format string: missing arguments");
    out << *s++;
  }
}

template<typename T, typename... Args>
static void SFormat(std::ostream &out,
                    const char *s,
                    const T &value,
                    const Args &... args) {
  while (*s) {
    if (*s == '{' && *++s == '}') {
      out << value;
      return SFormat(out, ++s, args...);
    }
    out << *s++;
  }
  printf("! SFormat problem, input: %s\n", s);
  throw std::runtime_error("extra arguments provided to printf");
}

/************************************/
/* Formatting and Logging functions */
template<typename... Args>
void error(const char *file,
           const unsigned long line,
           const char *format,
           const Args &... args) {
  fprintf(stdout, " %s:%i: ", file, line);
  std::ostringstream stream;
  SFormat(stream, format, args...);
  printf("- error :%s\n", stream.str().c_str());
}

template<typename... Args>
void warn(const char *file,
          const unsigned long line,
          const char *format,
          const Args &... args) {
  fprintf(stdout, " %s:%i: ", file, line);
  std::ostringstream stream;
  SFormat(stream, format, args...);
  printf("- warn :%s\n", stream.str().c_str());
}

template<typename... Args>
void info(const char *file,
          const unsigned long line,
          const char *format,
          const Args &... args) {
  fprintf(stdout, " %s:%i: ", file, line);
  std::ostringstream stream;
  SFormat(stream, format, args...);
  printf("- info :%s\n", stream.str().c_str());
}

template<typename... Args>
void debug(const char *file,
           const unsigned long line,
           const char *format,
           const Args &... args) {
  fprintf(stdout, " %s:%i: ", file, line);
  std::ostringstream stream;
  SFormat(stream, format, args...);
  printf("- debug :%s\n", stream.str().c_str());
}

template<class T>
void sort_unique(std::vector<T> &vec) {
  std::sort(vec.begin(), vec.end());
  vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
}

template<class T1, class T2>
T2 sort_unique_with_perm(
    const std::vector<T1> &in,
    std::vector<T1> &uniques,
    std::vector<T2> &old2new) {

  std::vector<T2> ids(in.size());
  for (T2 k = 0; k != in.size(); ++k) ids[k] = k;

  std::sort(ids.begin(), ids.end(),
            [&in](const T2 &a, const T2 &b) { return (in[a] < in[b]); }
  );

  uniques.resize(in.size());
  old2new.resize(in.size());
  for (T2 k = 0; k != in.size(); ++k) uniques[k] = in[k];

  std::sort(uniques.begin(), uniques.end());
  uniques.erase(std::unique(uniques.begin(), uniques.end()),
                uniques.end());
  T2 ic = 0; // index current
  T2 ir = 0; // index representation
  T2 cur_rep = 0; // index of current representation
  while (ic < in.size()) {
    ic = ir;
    while (ic < in.size() && in[ids[ic]] == in[ids[ir]]) {
      old2new[ids[ic]] = cur_rep;
      ++ic;
    }
    ir = ic;
    ++cur_rep;
  }
  return (T2)uniques.size();
}

template<class T>
void compress(const std::vector<std::vector<T> > &vov,
              std::vector<size_t> &first, std::vector<T> &values) {
  first.resize(vov.size() + 1);
  size_t count = 0;
  for (size_t i = 0; i < vov.size(); ++i) {
    first[i] = count;
    count += vov[i].size();
  }
  first[vov.size()] = count;
  values.resize(count);
  for (size_t i = 0; i < vov.size(); ++i) {
    for (size_t j = 0; j < vov[i].size(); ++j) {
      values[first[i] + j] = vov[i][j];
    }
  }
  first.shrink_to_fit();
  values.shrink_to_fit();
}

template<class T1, class T2>
inline std::vector<T2> dynamic_cast_vector(const std::vector<T1> &pointers) {
  std::vector<T2> output(pointers.size(), NULL);
  for (size_t i = 0; i < pointers.size(); ++i) {
    output[i] = dynamic_cast<T2>(pointers[i]);
  }
  return output;
}

template<class T>
std::vector<T> intersection(const std::vector<T> &v1,
                            const std::vector<T> &v2) {
  std::vector<T> s1 = v1;
  std::vector<T> s2 = v2;
  sort_unique(s1);
  sort_unique(s2);
  std::vector<T> s3;
  set_intersection(s1.begin(),
                   s1.end(),
                   s2.begin(),
                   s2.end(),
                   std::back_inserter(s3));
  return s3;
}

template<class T>
std::vector<T> merge(const std::vector<T> &v1, const std::vector<T> &v2) {
  std::vector<T> s1 = v1;
  std::vector<T> s2 = v2;
  sort_unique(s1);
  sort_unique(s2);
  std::vector<T> s3;
  set_union(s1.begin(), s1.end(), s2.begin(), s2.end(), std::back_inserter(s3));
  sort_unique(s3);
  return s3;
}

template<class T>
std::vector<T> difference(const std::vector<T> &v1, const std::vector<T> &v2) {
  std::vector<T> s1 = v1;
  std::vector<T> s2 = v2;
  sort_unique(s1);
  sort_unique(s2);
  std::vector<T> s3;
  set_difference(s1.begin(),
                 s1.end(),
                 s2.begin(),
                 s2.end(),
                 std::inserter(s3, s3.begin()));
  return s3;
}

template<class T>
void append(std::vector<T> &v1, const std::vector<T> &v2) {
  v1.insert(v1.end(), v2.begin(), v2.end());
}

template<typename T>
bool inVector(const T &value, const std::vector<T> &vec) {
  if (vec.empty()) return false;
  auto it = std::find(vec.begin(), vec.end(), value);
  return it != vec.end();
}

template<typename T, size_t N>
bool inArray(const T &value, const std::array<T, N> &vec) {
  if (vec.empty()) return false;
  auto it = std::find(vec.begin(), vec.end(), value);
  return it != vec.end();
}

/* For debug prints */
template<typename H1>
std::ostream &show(std::ostream &out, const char *label, H1 &&value) {
  return out << label << "=" << std::forward<H1>(value) << '\n';
}


template<typename H1, typename ...T>
std::ostream &show(std::ostream &out,
                   const char *label,
                   H1 &&value,
                   T &&... rest) {
  const char *p_comma = strchr(label, ',');
  return show(out.write(label, p_comma - label) << "="
                                                << std::forward<H1>(value)
                                                << ',',
              p_comma + 1,
              std::forward<T>(rest)...);
}

template<typename T>
void DebugOut(std::string file,T&&value){
  std::ofstream outfile;
  outfile.open(file+".txt");
  show(outfile,file.c_str(),value);
  outfile.close();
}

inline std::string env_var(const std::string &key) {
  char *var;
  var = getenv(key.c_str());
  std::string str_var;
  if (var != nullptr) {
    str_var = std::string(var);
  }
  return str_var;
}

inline void read_from_env(const std::string &key, double &value) {
  std::string str = env_var(key);
  if (str.empty()) return;
  value = std::stod(str);
}

inline void read_from_env(const std::string &key, int &value) {
  std::string str = env_var(key);
  if (str.empty()) return;
  value = std::stoi(str);
}

inline void read_from_env(const std::string &key, bool &value) {
  std::string str = env_var(key);
  if (!str.empty()) {
    if (str == "0" || str == "false" || str == "False") {
      value = false;
    } else if (str == "1" || str == "true" || str == "True") {
      value = true;
    }
  }
}

inline void read_from_env(const std::string &key, std::string &value) {
  value = env_var(key);
}




#endif //CROSSFIELD__BASIC_UTILS_H_
