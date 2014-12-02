
#ifndef CONSTRAINTS__H
#define CONSTRAINTS__H

#include "base.h"

// gecode
#include <gecode/minimodel.hh>
#include <gecode/float.hh>
#include <gecode/search.hh>
#include <gecode/gist.hh>

//#define CONSTRAINTS_DEBUG

#define EPSILON 0.001f

enum relational_op
  {
    op_equal,
    op_lesser,
    op_greater,
    op_lesser_equal,
    op_greater_equal
  };

const char* str_relop[5] =
  {
    "==",
    "<",
    ">",
    "<=",
    ">="
  };

template<relational_op op>
struct relational_operator
  {
    Gecode::BoolExpr exec(Gecode::LinFloatExpr expr1, Gecode::LinFloatExpr expr2)
      {
        assert(false); // not supported
      }
  };

template<>
struct relational_operator<op_equal>
  {
    Gecode::BoolExpr exec(Gecode::LinFloatExpr expr1, Gecode::LinFloatExpr expr2)
      {
        return expr1 == expr2;
      }
  };

template<>
struct relational_operator<op_lesser>
  {
    Gecode::BoolExpr exec(Gecode::LinFloatExpr expr1, Gecode::LinFloatExpr expr2)
      {
        return expr1 < expr2;
      }
  };

template<>
struct relational_operator<op_greater>
  {
    Gecode::BoolExpr exec(Gecode::LinFloatExpr expr1, Gecode::LinFloatExpr expr2)
      {
        return expr1 > expr2;
      }
  };

template<>
struct relational_operator<op_lesser_equal>
  {
    Gecode::BoolExpr exec(Gecode::LinFloatExpr expr1, Gecode::LinFloatExpr expr2)
      {
        return expr1 <= expr2;
      }
  };

template<>
struct relational_operator<op_greater_equal>
  {
    Gecode::BoolExpr exec(Gecode::LinFloatExpr expr1, Gecode::LinFloatExpr expr2)
      {
        return expr1 >= expr2;
      }
  };


// forwards
class constraint;
class container;

// references
typedef cstd::shared_ptr<constraint>              Constraint;
typedef cstd::shared_ptr<container>               Container;


class constraint
  {
    public:
      virtual void create(Gecode::Space &s, Container &c) = 0;
  };

class container
  {
    public:
      container()
        : has_point_(false) 
        {
        }

      container(point p)
        : point_(p), has_point_(true)
        {
        }

      void process_result()
        {
          values_.clear();
          for (int i = 0; i < vars_.size(); ++i)
            {
              double v = vars_[i].val().med();
              values_.push_back(v);
            }

          if (values_.size() > 1)
            { 
              point_ = point(values_[0], values_[1]);
            }

          // clean variable arguments array
          // exists a better method?
          vars_ = Gecode::FloatVarArgs();
        }

      point get_point()
        {
          return point_;
        }

      bool has_point() const
        {
          return has_point_;
        }

      void print_values() const
        {
          std::cout << std::fixed;
          std::cout.precision(2);

          std::cout << "{";
          for (size_t i = 0; i < values_.size(); ++i)
            {
              std::cout << values_[i] << " ";
            }
          std::cout << "}" << std::endl;
        }

      void print_vars() const
        {
          std::cout << vars_ << std::endl;
        }

      void print() const
        {
          if (has_point_)
            printf("(%0.3f;%0.3f)\n", point_.x(), point_.y());
          else
            { 
              std::cout << "( ";
              for (int i = 0; i < vars_.size(); ++i)
                {
                  printf("%x ", &vars_[i]);
                }
              printf(")\n");
            }
        }

      Gecode::FloatVarArgs& vars()
        {
          return vars_;
        }

    private:
      point                 point_;
      bool                  has_point_;
      std::vector<double>   values_;
      Gecode::FloatVarArgs  vars_;
  };

// find a random point that is inside the segment determined
// by the two points p & q
class point_in_segment : public constraint
  {
    public:
      point_in_segment(Container &p, Container &q)
        : p_(p), q_(q)
        {
        }

    protected:
      void create(Gecode::Space &s, Container &r)
        {
          // prerequisites
          assert(r);

          Gecode::FloatVarArgs &vars = r->vars();
          vars = Gecode::FloatVarArgs(s, 3, Gecode::Float::Limits::min, Gecode::Float::Limits::max);
          Gecode::FloatVar x(vars[0]), y(vars[1]), d(vars[2]);

#ifdef CONSTRAINTS_DEBUG
          std::cout << "point_in_segment{" << std::endl;
          std::cout << "\t"; p_->print();
          std::cout << "\t"; q_->print();
          std::cout << "\t"; r->print();
          std::cout << "}" << std::endl;
#endif

          // exclude d from equal zero
          Gecode::rel(s, d != 0);

          if (p_->has_point() && q_->has_point())
            { 
              point p1 = p_->get_point();
              Gecode::LinFloatExpr x1(p1.x());
              Gecode::LinFloatExpr y1(p1.y());
              point p2 = q_->get_point();
              Gecode::LinFloatExpr x2(p2.x());
              Gecode::LinFloatExpr y2(p2.y());

              // declare ranges for component variables
              Gecode::dom(s, x, std::min(p1.x(), p2.x()), std::max(p1.x(), p2.x()));
              Gecode::dom(s, y, std::min(p1.y(), p2.y()), std::max(p1.y(), p2.y()));

              // create parametric equations from point (p) and vector (u)
              //Gecode::rel(s, x == p1.x() + d * (p2.x() - p1.x()));
              //Gecode::rel(s, y == p1.y() + d * (p2.y() - p1.y()));
              Gecode::rel(s, x == x1 + d * (x2 - x1));
              Gecode::rel(s, y == y1 + d * (y2 - y1));
            }
          else
          if (p_->has_point() && !q_->has_point())
            {
              Gecode::FloatVarArgs &v2 = q_->vars();
              Gecode::FloatVar x2(v2[0]), y2(v2[1]);
              point p1 = p_->get_point();
              Gecode::LinFloatExpr x1(p1.x());
              Gecode::LinFloatExpr y1(p1.y());

              // restrict ranges for component variables
              //Gecode::rel(s, x >= Gecode::min(x1, x2) && x <= Gecode::max(x1, x2));
              //Gecode::rel(s, y >= Gecode::min(y1, y2) && y <= Gecode::max(y1, y2));
              Gecode::dom(s, x, std::min(p1.x(), x2.min()), std::max(p1.x(), x2.max()));
              Gecode::dom(s, y, std::min(p1.y(), y2.min()), std::max(p1.y(), y2.max()));

              // create parametric equations from point (p) and vector (u)
              Gecode::rel(s, x == x1 + d * (x2 - x1));
              Gecode::rel(s, y == y1 + d * (y2 - y1));
            }
          else
          if (!p_->has_point() && q_->has_point())
            {
              Gecode::FloatVarArgs &v1 = p_->vars();
              Gecode::FloatVar x1(v1[0]), y1(v1[1]);
              point p2 = q_->get_point();
              Gecode::LinFloatExpr x2(p2.x());
              Gecode::LinFloatExpr y2(p2.y());

              // declare domine for component variables
              //Gecode::rel(s, x >= Gecode::min(x1, x2) && x <= Gecode::max(x1, x2));
              //Gecode::rel(s, y >= Gecode::min(y1, y2) && y <= Gecode::max(y1, y2));
              Gecode::dom(s, x, std::min(x1.min(), p2.x()), std::max(x1.max(), p2.x()));
              Gecode::dom(s, y, std::min(y1.min(), p2.y()), std::max(y1.max(), p2.y()));

              // create parametric equations from point (p) and vector (u)
              Gecode::rel(s, x == x1 + d * (x2 - x1));
              Gecode::rel(s, y == y1 + d * (y2 - y1));
            }
          else
          if (!p_->has_point() && !q_->has_point())
            {
              Gecode::FloatVarArgs &v1 = p_->vars();
              Gecode::FloatVarArgs &v2 = q_->vars();
              Gecode::FloatVar x1(v1[0]), y1(v1[1]);
              Gecode::FloatVar x2(v2[0]), y2(v2[1]);

              // declare domine for component variables
              //Gecode::rel(s, x >= Gecode::min(x1, x2) && x <= Gecode::max(x1, x2));
              //Gecode::rel(s, y >= Gecode::min(y1, y2) && y <= Gecode::max(y1, y2));
              Gecode::dom(s, x, std::min(x1.min(), x2.min()), std::max(x1.max(), x2.max()));
              Gecode::dom(s, y, std::min(y1.min(), y2.min()), std::max(y1.max(), y2.max()));

              // create parametric equations from point (p) and vector (u)
              Gecode::rel(s, x == x1 + d * (x2 - x1));
              Gecode::rel(s, y == y1 + d * (y2 - y1));
            }
        }

    protected:
      Container   p_, q_;
  };

// find a random point that is inside the segment determined
// by the two points p & q, and the distance from p to the new point
// is greater, lesser or equal than "dist"
template <relational_op op>
class distance : public constraint
  {
    public:
      distance(Container &p, Container &r1, double dist)
        : p_(p), r1_(r1), dist_(dist)
        {
        }

    protected:
      virtual void create(Gecode::Space &s, Container &r)
        {
          // prerequisites
          assert(r);
          assert(r1_->vars().size() > 1);

          // reuse values assigned from r1 (point_in_segment)
          Gecode::FloatVarArgs &vars = r1_->vars();
          Gecode::FloatVar x(vars[0]), y(vars[1]);

#ifdef CONSTRAINTS_DEBUG
          std::cout << "distance{" << std::endl;
          std::cout << "\t"; p_->print();
          std::cout << "\t"; r1_->print();
          std::cout << "\tdist " << str_relop[op] << " " << dist_ << std::endl;
          std::cout << "}" << std::endl;
#endif

          relational_operator<op> executor;
          Gecode::LinFloatExpr left;

          if (p_->has_point())
            {
              point p1 = p_->get_point();
              left = Gecode::sqrt(Gecode::pow(p1.x() - x, 2) + Gecode::pow(p1.y() - y, 2));
            }
          else
            {
              Gecode::FloatVarArgs &v1 = p_->vars();
              Gecode::FloatVar x1(v1[0]), y1(v1[1]);
              left = Gecode::sqrt(Gecode::pow(x1 - x, 2) + Gecode::pow(y1 - y, 2));
            }

          Gecode::LinFloatExpr right(dist_);
          Gecode::BoolExpr bexpr = executor.exec(left, right);

          // make restriction about distance between points
          Gecode::rel(s, bexpr);
        }

    protected:
      Container p_, r1_;
      double    dist_;
};

// find two random points r1 & r2 that is inside segments that 
// are determined by point_in_segment constraint, and its parallel 
// with the segment determined by the two points p & q
class parallel_segments : public constraint
  {
    public:
      parallel_segments(Container &p, Container &q, Container &r1, Container &r2)
        : p_(p), q_(q), r1_(r1), r2_(r2)
        {
        }

    protected:
      virtual void create(Gecode::Space &s, Container &r)
        {
          // prerequisites
          assert(r);
          assert(r1_->vars().size() > 1);
          assert(r2_->vars().size() > 1);

          // create only one variable: sita
          Gecode::FloatVarArgs &vars = r->vars();
          vars = Gecode::FloatVarArgs(s, 1, Gecode::Float::Limits::min, Gecode::Float::Limits::max);
          Gecode::FloatVar sita(vars[0]);

          // reuse values assigned from point_in_segment objects
          Gecode::FloatVarArgs &r1 = r1_->vars();
          Gecode::FloatVarArgs &r2 = r2_->vars();
          Gecode::FloatVar x1(r1[0]), y1(r1[1]);
          Gecode::FloatVar x2(r2[0]), y2(r2[1]);

          // two segments are parallels when theirs vectors have equal directions
          // that means that u = sita * v, where u is a vector of segment (p, q), &
          // v is a vector (r1, r2)
          // http://www.tec-digital.itcr.ac.cr/revistamatematica/cursos-linea/Algebra-Lineal/algebra-vectorial-geova-walter/node5.html

          // create parametric equations from point (p) and vector (u)
          if (p_->has_point() && q_->has_point())
            { 
              point p = p_->get_point();
              point q = q_->get_point();

              Gecode::rel(s, (q.x() - p.x()) == sita * (x2 - x1));
              Gecode::rel(s, (q.y() - p.y()) == sita * (y2 - y1));
            }
          else
          if (p_->has_point() && !q_->has_point())
            {
              point p = p_->get_point();
              Gecode::FloatVarArgs &q = q_->vars();
              Gecode::FloatVar q_x(q[0]), q_y(q[1]);

              Gecode::rel(s, (q_x - p.x()) == sita * (x2 - x1));
              Gecode::rel(s, (q_y - p.y()) == sita * (y2 - y1));
            }
          else
          if (!p_->has_point() && q_->has_point())
            {
              Gecode::FloatVarArgs &p = p_->vars();
              Gecode::FloatVar p_x(p[0]), p_y(p[1]);
              point q = q_->get_point();

              Gecode::rel(s, (q.x() - p_x) == sita * (x2 - x1));
              Gecode::rel(s, (q.y() - p_y) == sita * (y2 - y1));
            }
          else
          if (!p_->has_point() && !q_->has_point())
            {
              Gecode::FloatVarArgs &p = p_->vars();
              Gecode::FloatVarArgs &q = q_->vars();
              Gecode::FloatVar p_x(p[0]), p_y(p[1]);
              Gecode::FloatVar q_x(q[0]), q_y(q[1]);

              Gecode::rel(s, (q_x - p_x) == sita * (x2 - x1));
              Gecode::rel(s, (q_y - p_y) == sita * (y2 - y1));
            }
        }

    protected:
      Container p_, q_;
      Container r1_, r2_;
  };

// find two random points r1 & r2 that is inside segments that
// are determined by point_in_segment constraint, and its perpendicular 
// with the segment determined by the two points p & q
class perpendicular_segment : public constraint
  {
    public:
      perpendicular_segment(Container &p, Container &q, Container &r1, Container &r2)
        : p_(p), q_(q), r1_(r1), r2_(r2)
        {
        }

    protected:
      virtual void create(Gecode::Space &s, Container &r)
        {
          // prerequisites
          assert(r);
          assert(r1_->vars().size() > 1);
          assert(r2_->vars().size() > 1);

          // reuse values assigned from point_in_segment objects
          Gecode::FloatVarArgs &r1 = r1_->vars();
          Gecode::FloatVarArgs &r2 = r2_->vars();
          Gecode::FloatVar x1(r1[0]), y1(r1[1]);
          Gecode::FloatVar x2(r2[0]), y2(r2[1]);

          // two segments are perpendicular when the product of theirs vectors are 
          // equal to zero, that means that u * v = 0, where u is a vector of 
          // segment (p, q), & v is a vector (res1, res2)
          // https://ar.answers.yahoo.com/question/index?qid=20090407081543AACvaIx

          // create equation with scalar product of vectors u & v
          if (p_->has_point() && q_->has_point())
            { 
              point p = p_->get_point();
              point q = q_->get_point();

              Gecode::rel(s, (q.x() - p.x()) * (x2 - x1) + (q.y() - p.y()) * (y2 - y1) == 0);
            }
          else
          if (p_->has_point() && !q_->has_point())
            {
              point p = p_->get_point();
              Gecode::FloatVarArgs &q = q_->vars();
              Gecode::FloatVar q_x(q[0]), q_y(q[1]);

              Gecode::rel(s, (q_x - p.x()) * (x2 - x1) + (q_y - p.y()) * (y2 - y1) == 0);
            }
          else
          if (!p_->has_point() && q_->has_point())
            {
              Gecode::FloatVarArgs &p = p_->vars();
              Gecode::FloatVar p_x(p[0]), p_y(p[1]);
              point q = q_->get_point();

              Gecode::rel(s, (q.x() - p_x) * (x2 - x1) + (q.y() - p_y) * (y2 - y1) == 0);
            }
          else
          if (!p_->has_point() && !q_->has_point())
            {
              Gecode::FloatVarArgs &p = p_->vars();
              Gecode::FloatVarArgs &q = q_->vars();
              Gecode::FloatVar p_x(p[0]), p_y(p[1]);
              Gecode::FloatVar q_x(q[0]), q_y(q[1]);

              Gecode::rel(s, (q_x - p_x) * (x2 - x1) + (q_y - p_y) * (y2 - y1) == 0);
            }
        }

    protected:
      Container p_, q_;
      Container r1_, r2_;
  };

// find a random point r1 that is inside segment that is determined
// by point_in_segment and make a new segment with a point q;
// the points p & q determine a segment such that with a new segment
// form an angle with amplitude more than "angle"
template <relational_op op>
class fix_angle : public constraint
  {
    public:
      fix_angle(const point &p, const point &q, Container &r1, double ang)
        : p_(p), q_(q), u_(q, p), res1_(r1), angle_(ang)
        {
        }

    protected:
      virtual void create(Gecode::Space &s, Container &r)
        {
          // prerequisites
          assert(r);
          assert(res1_->vars_.size() > 1);

          // reuse values assigned from point_in_segment objects
          Gecode::FloatVar x(res1_->vars_[0]), y(res1_->vars_[1]);

          relational_operator<op> executor;
          Gecode::LinFloatExpr left(cos(angle_));
          Gecode::LinFloatExpr right =
                          (u_.x() * (x - q_.x()) + u_.y() * (y - q_.y())) 
                          /
                          sqrt(
                                (pow(u_.x(), 2) + pow(u_.y(), 2)) *
                                (pow(x - q_.x(), 2) + pow(y - q_.y(), 2))
                              );

          Gecode::BoolExpr bexpr = executor.exec(left, right);

          // apply equation to determine the angle between two vectors
          Gecode::rel(s, bexpr);
        }

    protected:
      point   p_, q_;
      vec2    u_;
      Container  res1_;
      double  angle_;
  };

// find a random point r1 that is inside segment that is determined
// by point_in_segment and make a new segment with a point q;
// the points p & q determine a segment such that with a new segment
// form an angle with amplitude more than "angle"
template <relational_op op>
class angle : public constraint
  {
    public:
      angle(Container &a, Container &b, Container &c, double ang)
        : a_(a), b_(b), c_(c), angle_(ang)
        {
        }

    protected:
      virtual void create(Gecode::Space &s, Container &r)
        {
          // prerequisites
          assert(r);

          // td: how to check if point "p" & point_in_segment "r1" should be points of the same segment

          // reuse values assigned from point_in_segment objects
          Gecode::FloatVar a_x, a_y;
          Gecode::FloatVar b_x, b_y;
          Gecode::FloatVar c_x, c_y;

          //uffsss... eight combinations!!
          if (a_->has_point())
            {
              point a = a_->get_point();
              a_x = Gecode::FloatVar(s, Gecode::Float::Limits::min, Gecode::Float::Limits::max);
              a_y = Gecode::FloatVar(s, Gecode::Float::Limits::min, Gecode::Float::Limits::max);
              Gecode::rel(s, a_x == a.x() && a_y == a.y());
            }
          else
            { 
              Gecode::FloatVarArgs &a = a_->vars();
              a_x = Gecode::FloatVar(a[0]);
              a_y = Gecode::FloatVar(a[1]);
            }

          if (b_->has_point())
            {
              point b = b_->get_point();
              b_x = Gecode::FloatVar(s, Gecode::Float::Limits::min, Gecode::Float::Limits::max);
              b_y = Gecode::FloatVar(s, Gecode::Float::Limits::min, Gecode::Float::Limits::max);
              Gecode::rel(s, b_x == b.x() && b_y == b.y());
            }
          else
            { 
              Gecode::FloatVarArgs &b = b_->vars();
              b_x = Gecode::FloatVar(b[0]);
              b_y = Gecode::FloatVar(b[1]);
            }

          if (c_->has_point())
            {
              point c = c_->get_point();
              c_x = Gecode::FloatVar(s, Gecode::Float::Limits::min, Gecode::Float::Limits::max);
              c_y = Gecode::FloatVar(s, Gecode::Float::Limits::min, Gecode::Float::Limits::max);
              Gecode::rel(s, c_x == c.x() && c_y == c.y());
            }
          else
            {
              Gecode::FloatVarArgs &c = c_->vars();
              c_x = Gecode::FloatVar(c[0]);
              c_y = Gecode::FloatVar(c[1]);
            }

          relational_operator<op> executor;
          Gecode::LinFloatExpr left(cos(angle_));
          Gecode::LinFloatExpr right =
                          ((c_x - b_x) * (a_x - b_x) + (c_y - b_y) * (a_y - b_y)) 
                          /
                          sqrt(
                                (pow(c_x - b_x, 2) + pow(c_y - b_y, 2)) *
                                (pow(a_x - b_x, 2) + pow(a_y - b_y, 2))
                              );

          Gecode::BoolExpr bexpr = executor.exec(left, right);

          // apply equation to determine the angle between two vectors
          Gecode::rel(s, bexpr);
        }

    protected:
      Container a_, b_, c_;
      double    angle_;
  };

class constraints_manager
  {
    public:
      constraints_manager()
        {
        }

      // determine the array of points with the solution of all constraints
      bool solve()
        {
          ConstraintsModel cm(new constraints_model());

          for (size_t i = 0; i < result_list.size(); ++i)
            {
              cm->add(constraint_list[i], result_list[i]);
            }

          // before search the first, post branching
          cm->post_branching();

#ifdef CONSTRAINTS_DEBUG
          cm->print(std::cout);
          (void)cm->status();
          cm->print(std::cout);

          Gecode::Gist::Print<constraints_model> p("Print solution");
          Gecode::Gist::Options o;
          o.inspect.click(&p);
          Gecode::Gist::dfs(cm.get(), o);
          //Gecode::Gist::bab(cm.get(), o);
#endif

          //td: custom this time; export interface
          Gecode::Search::TimeStop stop(2000);
          
          Gecode::Search::Options opts;
          opts.stop = &stop;

          Gecode::DFS<constraints_model> search(cm.get(), opts);
          //Gecode::BAB<constraints_model> search(cm.get(), opts);

          // search and return first solution
          if (constraints_model* s = search.next())
            {
              //s->print();
              s->update_results(cm);

              delete s;
              return true;
            }

          return false;
        }

      // add new a constraint
      void add(Constraint &cstrt, Container &res)
        {
          constraint_list.push_back(cstrt);
          result_list.push_back(res);
        }

    private:
      class constraints_model;

      typedef cstd::shared_ptr<constraints_model> ConstraintsModel;

      class constraints_model : public Gecode::Space
        {
          public:
            constraints_model()
              {
              }

            constraints_model(bool share, constraints_model& cm)
              : Space(share, cm)
              {
                fvars.update(*this, share, cm.fvars);
              }

            virtual Space* copy(bool share)
              {
                return new constraints_model(share, *this);
              }

            // add new constraints to space model
            void add(Constraint &cstrt, Container &cont)
              {
                // call virtual method with particular constraints
                cstrt->create(*this, cont);

                // push all local argument variables into
                // manager array of argument variable;
                // this is important doing, but not any variable will be branched
                fargs << cont->vars();
                vres.push_back(cont);
              }

            // post branching with variable array
            void post_branching()
              {
                // create variable array from all added variable arguments
                fvars = Gecode::FloatVarArray(*this, fargs);

                // create random selector by local time
                Gecode::Rnd r;
                r.time();

                // post the branch picking a random solution
                Gecode::branch(*this, fvars, Gecode::FLOAT_VAR_SIZE_MIN(), Gecode::FLOAT_VAL_SPLIT_RND(r));
                //Gecode::branch(*this, fvars, Gecode::FLOAT_VAR_SIZE_MAX(), Gecode::FLOAT_VAL_SPLIT_RND(r));
              }

            void update_results(ConstraintsModel cmodel)
              {
                size_t idx = 0;
                int counter = 0;

                for (int i = 0; i < fvars.size(); ++i)
                  {
                    Container r = cmodel->vres[idx];
                    Gecode::FloatVarArgs &vars = r->vars();
                    vars[counter++] = fvars[i];

                    if (counter == vars.size())
                      {
                        counter = 0;
                        r->process_result();

                        while (++idx < cmodel->vres.size() && cmodel->vres[idx]->vars().size() <= counter);
                      }
                  }
              }

            // print variable array
            // only for debug purposes 
            void print(std::ostream& os) const
              {
                os << fvars << std::endl;
              }

          private:
            Gecode::FloatVarArray   fvars;
            Gecode::FloatVarArgs    fargs;
            std::vector<Container>  vres;
        };

    private:
      std::vector<Constraint> constraint_list;
      std::vector<Container>  result_list;
  };

#endif
