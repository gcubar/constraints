
#ifndef BASE__H
#define BASE__H

#define QUOTE(arg) <arg>

// GNU GCC/G++ between 4.6 and 4.8.2
#if (defined(__GNUC__) || defined(__GNUG__)) && ((__GNUC__ == 4) && (__GNUC_MINOR__ >= 6) && (__GNUC_MINOR__ < 8) && (__GNUC_PATCHLEVEL__ < 2))
#  define TR1IFY(arg) tr1/arg
namespace cstd = std::tr1;
#else
#  define TR1IFY(arg) arg
namespace cstd = std;
// GCC > 4.8.1: must check the compiler option: -std=c++11
#endif

#define TR1INCLUDE(arg) QUOTE(TR1IFY(arg))


#include TR1INCLUDE(memory)
#include <vector>
#include <iostream>


#ifndef M_PI
#define M_PI        3.14159265358979323846
#endif

#ifndef M_PI_2
#define M_PI_2      1.57079632679489661923
#endif

struct point
{
    point()
      : x_(0.0), y_(0.0) {}

    point(double x, double y)
      : x_(x), y_(y) {}

    point(const point& p)
      : x_(p.x()), y_(p.y()) {}

    double x() const { return x_; }
    double y() const { return y_; }

  private:
    double x_, y_;
};

struct vec2
{
    vec2()
      : x_(0.0), y_(0.0) {}

    vec2(double x, double y)
      : x_(x), y_(y) {}

    vec2(const point &a, const point &b)
      : x_(b.x() - a.x()), y_(b.y() - a.y()) {}

    double x() const { return x_; }
    double y() const { return y_; }

  private:
    double x_, y_;
};

template<class Char, class Traits>
std::basic_ostream<Char, Traits>&
operator << (std::basic_ostream<Char, Traits>& os, const point& p)
{
    std::cout << std::fixed;
    std::cout.precision(2);
    os << "(" << p.x() << ", " << p.y() << ")";
    return os;
}

#endif
