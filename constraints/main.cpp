
#include <iostream>

#include "constraints.h"
#include "gtest/gtest.h"

struct tools
{
    // http://stackoverflow.com/questions/328107/how-can-you-determine-a-point-is-between-two-other-points-on-a-line-segment
    static bool   is_point_inside_segment(point p1, point p2, point o, double epsilon = 0.001f)
    {
        // check if the points p1, p2 & o are aligned 
        double crossproduct = (o.y() - p1.y()) * (p2.x() - p1.x()) - (o.x() - p1.x()) * (p2.y() - p1.y());
        if (abs(crossproduct) > epsilon)
            return false;

        // check if dot product of (p2 - p1) & (o - p1) is positive
        double dotproduct = (o.x() - p1.x()) * (p2.x() - p1.x()) + (o.y() - p1.y()) * (p2.y() - p1.y());
        if (dotproduct < 0)
            return false;

        // check if dot product is less than square distance between p1 & p2
        double squared_length_p2_p1 = (p2.x() - p1.x()) * (p2.x() - p1.x()) + (p2.y() - p1.y()) * (p2.y() - p1.y());
        if (dotproduct > squared_length_p2_p1)
            return false;

        return true;
    }

    static double distance(point p1, point p2)
    {
        double dx = p1.x() - p2.x();
        double dy = p1.y() - p2.y();
        return sqrt( pow(dx, 2) + pow(dy, 2) );
    }

    static bool   parallels(point p1, point p2, point p3, point p4, double epsilon = 0.001f)
    {
        double dx1 = p1.x() - p2.x();
        double dy1 = p1.y() - p2.y();
        double dx2 = p3.x() - p4.x();
        double dy2 = p3.y() - p4.y();

        if (dx2 == 0)
        {
            return dx2 == dx1;
        }

        if (dy2 == 0)
        {
            return dy2 == dy1;
        }

        return abs(dx1 / dx2 - dy1 / dy2) < epsilon;
    }

    static bool   perpendiculars(point p1, point p2, point p3, point p4, double epsilon = 0.001f)
    {
        double dx1 = p1.x() - p2.x();
        double dy1 = p1.y() - p2.y();
        double dx2 = p3.x() - p4.x();
        double dy2 = p3.y() - p4.y();

        return dx1 * dx2 + dy1 * dy2 < epsilon;
    }
};

const int iterations = 200;
Container r_empty(new container());

point _M(-2, 1);
point _N(-2, 4);
point _O(8, 4);
point _P(8, 1);
point _Q(3, 1);

point _A(1, 1);
point _B(3, 1);
point _C(-1, 3);
point _D(5, 4);

Container M(new container(_M));
Container N(new container(_N));
Container O(new container(_O));
Container P(new container(_P));
Container Q(new container(_Q));

Container A(new container(_A));
Container B(new container(_B));
Container C(new container(_C));
Container D(new container(_D));

/*
TEST(Constraints, PointInSegment)
{
    Container r1(new container());
    Container r2(new container());

    Constraint c1(new point_in_segment(N, O));
    Constraint c2(new point_in_segment(r1, Q));

    constraints_manager cm;
    cm.add(c1, r1);
    cm.add(c2, r2);

    bool result;
    for (int i = 0; i < iterations; ++i)
    {
        // any new call to solve(), take a random solution
        if (cm.solve())
        {
            //std::cout << r1->get_point() << std::endl;
            //std::cout << r2->get_point() << std::endl;

            result = tools::is_point_inside_segment(_N, _O, r1->get_point());
            EXPECT_TRUE(result);

            result = tools::is_point_inside_segment(r1->get_point(), _Q, r2->get_point());
            EXPECT_TRUE(result);
        }
        else
        {
            std::cout << "[No solution found!]" << std::endl;
        }
    }
}

TEST(Constraints, Distance)
{
    Container r1(new container());
    Container r2(new container());

    Constraint c1(new point_in_segment(N, O));
    Constraint c2(new point_in_segment(O, P));
    Constraint c3(new distance<op_lesser>(N, r1, 6));
    Constraint c4(new distance<op_lesser>(O, r1, 6));
    Constraint c5(new distance<op_greater_equal>(O, r2, 1.5f));
    Constraint c6(new distance<op_greater_equal>(P, r2, 1.5f));

    constraints_manager cm;
    cm.add(c1, r1);
    cm.add(c2, r2);
    cm.add(c3, r_empty);
    cm.add(c4, r_empty);
    cm.add(c5, r_empty);
    cm.add(c6, r_empty);

    bool result;
    float d1, d2;
    for (int i = 0; i < iterations; ++i)
    {
        // any new call to solve(), take a random solution
        if (cm.solve())
        {
            //std::cout << r1->get_point() << std::endl;
            //std::cout << r2->get_point() << std::endl;

            result = tools::is_point_inside_segment(_N, _O, r1->get_point());
            EXPECT_TRUE(result);

            result = tools::is_point_inside_segment(_O, _P, r2->get_point());
            EXPECT_TRUE(result);

            d1 = tools::distance(_N, r1->get_point());
            d2 = tools::distance(_O, r1->get_point());
            result = d1 < 6 && d2 < 6;
            EXPECT_TRUE(result);

            d1 = tools::distance(_O, r2->get_point());
            d2 = tools::distance(_P, r2->get_point());
            result = d1 >= 1.5f && d2 >= 1.5f;
            EXPECT_TRUE(result);
        }
        else
        {
            std::cout << "[No solution found!]" << std::endl;
        }
    }
}

TEST(Constraints, Parallels)
{
    Container r1(new container());
    Container r2(new container());

    Constraint c1(new point_in_segment(N, O));
    Constraint c2(new point_in_segment(M, P));
    Constraint c3(new parallel_segments(M, N, r1, r2));

    constraints_manager cm;
    cm.add(c1, r1);
    cm.add(c2, r2);
    cm.add(c3, r_empty);

    bool result;
    for (int i = 0; i < iterations; ++i)
    {
        // any new call to solve(), take a random solution
        if (cm.solve())
        {
            //std::cout << r1->get_point() << std::endl;
            //std::cout << r2->get_point() << std::endl;

            result = tools::is_point_inside_segment(_N, _O, r1->get_point());
            EXPECT_TRUE(result);

            result = tools::is_point_inside_segment(_M, _P, r2->get_point());
            EXPECT_TRUE(result);

            result = tools::parallels(_M, _N, r1->get_point(), r2->get_point());
            EXPECT_TRUE(result);
        }
        else
        {
            std::cout << "[No solution found!]" << std::endl;
        }
    }
}

TEST(Constraints, Perpendiculars)
{
    Container r1(new container());
    Container r2(new container());

    Constraint c1(new point_in_segment(N, O));
    Constraint c2(new point_in_segment(M, P));
    Constraint c3(new perpendicular_segment(N, O, r1, r2));

    constraints_manager cm;
    cm.add(c1, r1);
    cm.add(c2, r2);
    cm.add(c3, r_empty);

    bool result;
    for (int i = 0; i < iterations; ++i)
    {
        // any new call to solve(), take a random solution
        if (cm.solve())
        {
            //std::cout << r1->get_point() << std::endl;
            //std::cout << r2->get_point() << std::endl;

            result = tools::is_point_inside_segment(_N, _O, r1->get_point());
            EXPECT_TRUE(result);

            result = tools::is_point_inside_segment(_M, _P, r2->get_point());
            EXPECT_TRUE(result);

            result = tools::perpendiculars(_N, _O, r1->get_point(), r2->get_point());
            EXPECT_TRUE(result);
        }
        else
        {
            std::cout << "[No solution found!]" << std::endl;
        }
    }
}

TEST(Constraints, MultipleConstraints1)
{
    Container r1(new container());
    Container r2(new container());
    Container r3(new container());
    Container r4(new container());
    Container r5(new container());

    Constraint c1(new point_in_segment(B, D));
    Constraint c2(new distance<op_lesser>(D, r1, 1));
    Constraint c3(new point_in_segment(A, C));
    Constraint c4(new point_in_segment(B, D));
    Constraint c5(new parallel_segments(C, r1, r2, r3));
    Constraint c6(new distance<op_equal>(A, r2, 2));
    Constraint c7(new point_in_segment(C, D));
    Constraint c8(new perpendicular_segment(C, A, r2, r5));

    constraints_manager cm;
    cm.add(c1, r1);
    cm.add(c2, r_empty);
    cm.add(c3, r2);
    cm.add(c4, r3);
    cm.add(c5, r4);
    cm.add(c6, r_empty);
    cm.add(c7, r5);
    cm.add(c8, r_empty);

    bool result;
    float dist;
    for (int i = 0; i < iterations; ++i)
    {
        if (cm.solve())
        {
            //std::cout << r1->get_point() << std::endl;
            //std::cout << r2->get_point() << std::endl;
            //std::cout << r3->get_point() << std::endl;
            //r4->print_values();
            //std::cout << r5->get_point() << std::endl;
            //r_empty->print_values();

            result = tools::is_point_inside_segment(_B, _D, r1->get_point());
            EXPECT_TRUE(result);

            result = tools::is_point_inside_segment(_A, _C, r2->get_point());
            EXPECT_TRUE(result);

            result = tools::is_point_inside_segment(_B, _D, r3->get_point());
            EXPECT_TRUE(result);

            result = tools::is_point_inside_segment(_C, _D, r5->get_point());
            EXPECT_TRUE(result);

            dist = tools::distance(_D, r1->get_point());
            EXPECT_TRUE(dist < 1);

            dist = tools::distance(_A, r2->get_point());
            EXPECT_TRUE(dist - 2 < 0.001f);

            result = tools::parallels(_C, r1->get_point(), r2->get_point(), r3->get_point());
            EXPECT_TRUE(result);

            result = tools::perpendiculars(_C, _A, r2->get_point(), r5->get_point());
            EXPECT_TRUE(result);
        }
        else
        {
            std::cout << "[No solution found!]" << std::endl;
        }
    }
}
*/

TEST(Constraints, MultipleConstraints2)
{
    /*point _T(0, 0);
    point _U(4.244, 0);
    point _V(4.244, 4.113);
    point _W(0, 4.113);*/
    point _T(0, 0);
    point _U(4, 0);
    point _V(4, 4);
    point _W(0, 4);

    Container T(new container(_T));
    Container U(new container(_U));
    Container V(new container(_V));
    Container W(new container(_W));

    Container r1(new container());
    Container r2(new container());
    Container r3(new container());
    Container r4(new container());

    Constraint c1(new point_in_segment(T, U));
    Constraint c2(new point_in_segment(W, V));
    Constraint c3(new point_in_segment(U, V));
    Constraint c4(new point_in_segment(r1, r2));

    //Constraint c5(new distance<op_greater>(U, r1, 2));
    //Constraint c6(new distance<op_lesser>(V, r2, 2));
    Constraint c5(new distance<op_greater_equal>(U, r1, 2));
    Constraint c6(new distance<op_lesser_equal>(V, r2, 2));

    Constraint c7(new perpendicular_segment(r1, r2, r3, r4));

    constraints_manager cm;
    cm.add(c1, r1);
    cm.add(c2, r2);
    //cm.add(c3, r3);

    cm.add(c4, r4);

    cm.add(c5, r_empty);
    cm.add(c6, r_empty);

    //cm.add(c7, r_empty);

    bool result;
    float dist;
    for (int i = 0; i < iterations; ++i)
    {
        if (cm.solve())
        {
            //std::cout << r1->get_point() << std::endl;
            //std::cout << r2->get_point() << std::endl;
            //std::cout << r3->get_point() << std::endl;
            //std::cout << r4->get_point() << std::endl;
            //r_empty->print_values();

            result = tools::is_point_inside_segment(_T, _U, r1->get_point());
            EXPECT_TRUE(result);

            result = tools::is_point_inside_segment(_W, _V, r2->get_point());
            EXPECT_TRUE(result);

            //result = tools::is_point_inside_segment(_U, _V, r3->get_point());
            //EXPECT_TRUE(result);

            result = tools::is_point_inside_segment(r1->get_point(), r2->get_point(), r4->get_point());
            EXPECT_TRUE(result);

            dist = tools::distance(_U, r1->get_point());
            EXPECT_TRUE(dist > 2);

            dist = tools::distance(_V, r2->get_point());
            EXPECT_TRUE(dist < 2);

            //result = tools::perpendiculars(r1->get_point(), r2->get_point(), r3->get_point(), r4->get_point());
            //EXPECT_TRUE(result);
        }
        else
        {
            std::cout << "[No solution found!]" << std::endl;
        }
    }
}

/*
// main function
int main(int argc, char* argv[])
{
    try
    {
        //test0();
        //test1();
        //test2();
        test3();


        // case 3
        std::cout << "-[CASE 3]-------------------------" << std::endl;
        Container r_cd1(new container());
        Container r_cd2(new container());
        Container r_ac(new container());

        Constraint ctt_pis1(new point_in_segment(C, D));
        Constraint ctt_pis2(new point_in_segment(C, D));
        Constraint ctt_pis3(new point_in_segment(A, C));
        Constraint ctt_mina(new angle<op_lesser>(A, B, r_cd1, M_PI / 6));
        //Constraint ctt_perp(new perpendicular_segment(A, C, r_cd2, r_ac));

        constraints_manager cm3;
        cm3.add(ctt_pis1, r_cd1);
        cm3.add(ctt_pis2, r_cd2);
        cm3.add(ctt_pis3, r_ac);
        cm3.add(ctt_mina, r_empty);
        //cm3.add(ctt_perp, r_empty);

        if (cm3.solve())
        {
            r_cd1->print_values();
            r_cd2->print_values();
            r_ac->print_values();
            r_empty->print_values();
            std::cout << "--------------------------" << std::endl;
        }
        else
        {
            std::cout << "ERROR: Solution not found!" << std::endl;
        }
        // ...

        // case 4
        std::cout << "-[CASE 4]-------------------------" << std::endl;
        point M(1, 1);
        point N(1, 4);
        point O(5, 4);
        point P(5, 1);
        point Q(3, 1);

        Container r_no(new container());
        Container r_mp(new container());
        Container r_rq(new container());

        Constraint ctt_pis_NO(new point_in_segment(N, O));
        Constraint ctt_pis_MP(new point_in_segment(M, P));
        Constraint ctt_mina_MQR(new free_angle<op_greater>(M, r_mp, r_no, M_PI / 2));
        Constraint ctt_pis_RQ(new point_in_segment(r_no, r_mp));

        constraints_manager cm4;
        cm4.add(ctt_pis_NO, r_no);
        cm4.add(ctt_pis_MP, r_mp);
        cm4.add(ctt_mina_MQR, r_empty);
        cm4.add(ctt_pis_RQ, r_rq);

        if (cm4.solve())
        {
            r_no->print_values();
            r_mp->print_values();
            r_empty->print_values();
            r_rq->print_values();
            std::cout << "--------------------------" << std::endl;
        }
        else
        {
            std::cout << "ERROR: Solution not found!" << std::endl;
        }
        // ...
    }
    catch (Gecode::Exception e)
    {
        std::cerr << "Gecode exception: " << e.what() << std::endl;
    }

    return 0;
}
*/