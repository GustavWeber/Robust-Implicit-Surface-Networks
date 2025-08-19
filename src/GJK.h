
#ifndef ROBUST_IMPLICIT_SURFACE_NETWORKS_GJK
#define ROBUST_IMPLICIT_SURFACE_NETWORKS_GJK

#include <CGAL/Kd_tree_rectangle.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Kernel/global_functions_3.h>
#include <CGAL/Origin.h>
#include <CGAL/enum.h>
#include <vector>

template<typename Kernel>
typename Kernel::Point_3 support_point(typename Kernel::Tetrahedron_3 tet, typename Kernel::Vector_3 direction){
    typename Kernel::Point_3 point;
    typename Kernel::FT max_scalar_product;
    for(int i = 0; i<4; i++){
        typename Kernel::FT scalar_product = CGAL::scalar_product(tet[i]-CGAL::Origin(), direction);
        if(scalar_product > max_scalar_product){
            point = tet[i];
            max_scalar_product = scalar_product;
        }
    }
    return point;
}

template<typename Kernel>
typename Kernel::Point_3 support_point(CGAL::Kd_tree_rectangle<typename Kernel::FT> rect, typename Kernel::Vector_3 direction){
    typename Kernel::Point_3 point;
    typename Kernel::FT max_scalar_product;
    std::array<typename Kernel::FT, 2> x_vals = {rect.min_coord(0), rect.max_coord(0)};
    std::array<typename Kernel::FT, 2> y_vals = {rect.min_coord(1), rect.max_coord(1)};
    std::array<typename Kernel::FT, 2> z_vals = {rect.min_coord(2), rect.max_coord(2)};
    for(auto x: x_vals){
        for(auto y: y_vals){
            for(auto z: z_vals){
                typename Kernel::FT scalar_product = CGAL::scalar_product(direction, typename Kernel::Vector_3(x, y, z));
                if(scalar_product > max_scalar_product){
                    point = {x, y, z};
                    max_scalar_product = scalar_product;
                }
            }
        }
    }
}

template<typename Kernel>
typename Kernel::Point_3 get_minkowski_diff_point(CGAL::Kd_tree_rectangle<typename Kernel::FT> rect, typename Kernel::Tetrahedron_3 tet, typename Kernel::Vector_3 direction){
    return CGAL::Origin() + support_point(rect, direction) - support_point(tet, -direction);
}

template<typename Kernel>
bool simplex_contains_origin(const std::vector<typename Kernel::Point_3> simplex){
typedef typename Kernel::Vector_3 Vector_3;
typedef typename Kernel::Plane_3 Plane_3;
    switch (simplex.size()) {
        case 1:
        {
            return simplex[0] == CGAL::Origin();
        }
        case 2:
        {
            return Kernel::Segment_3(simplex[0], simplex[1]).has_on(CGAL::Origin());
        }
        case 3:
        {
            return Kernel::Triangle_3(simplex[0], simplex[1], simplex[2]).has_on(CGAL::Origin());
        }
        case 4:
        {
            CGAL::Bounded_side side = Kernel::Tetrahedron_3(simplex[0], simplex[1], simplex[2], simplex[3]).bounded_side(CGAL::Origin());
            return side != CGAL::Bounded_side::ON_UNBOUNDED_SIDE;
        }
        default:
            assert(false);
            return false; //this should never happen
    }
}

template<typename Kernel>
void get_new_simplex_and_direction(std::vector<typename Kernel::Point_3>& simplex, typename Kernel::Vector_3& direction){
typedef typename Kernel::Vector_3 Vector_3;
typedef typename Kernel::Point_3 Point_3;

    if(simplex.size() == 4){
        //remove all vertices, that are on the opposite site of the opposing face than the origin 
        std::vector<Point_3> new_simplex;
        for(int i = 0; i<4; i++)
        {
            typename Kernel::Plane_3 plane(simplex[(i+1)%4], simplex[(i+2)%4], simplex[(i+3)%4]);
            if(plane.oriented_side(CGAL::Origin()) == plane.oriented_side(simplex[i])){
                new_simplex.emplace_back(simplex[i]);
            }
        }
        simplex = new_simplex;

        //since the origin does not lie inside the tet, we must have removed at least one vertex
    }

    if(simplex.size() == 3){
        //for each vertex: build the plane spanned by the vector between the two other verts and the normal of the triangle. if the vertex is in the same half-space to that plane as the origin, keep it 
        std::vector<Point_3> new_simplex;
        Vector_3 plane_normal = CGAL::cross_product(simplex[1]-simplex[0], simplex[2]-simplex[0]);
        for(int i = 0; i<3; i++){
            typename Kernel::Plane_3 plane(simplex[(i+1)%3], CGAL::cross_product(plane_normal, simplex[(i+2)%3]-simplex[(i+1)%3]));
            if(plane.oriented_side(CGAL::Origin()) == plane.oriented_side(simplex[i])){
                new_simplex.emplace_back(simplex[i]);
            }
        }
        simplex = new_simplex;

        //if we didn't shrink the simplex, the projected point must lie inside the triangle.
        if(simplex.size() == 3){
            direction = CGAL::Origin() - typename Kernel::Plane_3(simplex[0], simplex[1], simplex[2]).projection(CGAL::Origin());
            return;
        }
    }

    if(simplex.size() == 2){
        //Project origin onto line, if it lands beyond either endpoint, remove the other
        Point_3 a = simplex[0];
        Point_3 b = simplex[1];
        Vector_3 a0 = CGAL::Origin() - a;
        Vector_3 ab = b - a;

        //Get scaling factor of projecting a0 onto ab 
        typename Kernel::FT scalar = CGAL::scalar_product(a0, ab) / ab.squared_length();
        if(scalar < 0){
            simplex = {a};
        } 
        else if(scalar > 1){
            simplex = {b};
        }
        else {
            //the projected point must lie on the line
            direction = CGAL::Origin() - (a+ab*scalar);
            return;
        }
    }
    
    // Only one vertex can be left.
    direction = CGAL::Origin() - simplex[0];
}

template<typename Kernel>
typename Kernel::FT GetDistance(typename Kernel::Tetrahedron_3 tet, CGAL::Kd_tree_rectangle<typename Kernel::FT> rect){
typedef typename Kernel::Point_3 Point;
typedef typename Kernel::Vector_3 Vector;
    std::vector<Point> simplex = {};
    Vector direction(1, 0, 0); // arbitrary starting direction
    simplex.emplace_back(get_minkowski_diff_point<Kernel>(rect, tet, direction));
    while(true){

        // if simplex contains origin return 0
        if(simplex_contains_origin<Kernel>(simplex)) return 0;
        // reduce simplex and get new direction
        get_new_simplex_and_direction<Kernel>(simplex, direction); 
        // use direction to get new support point 
        Point new_point = get_minkowski_diff_point<Kernel>(rect, tet, direction);
        // if that new support point is not further along the search direction than the old closest point, return distance between simplex and origin 
        if(CGAL::scalar_product(direction, new_point-CGAL::Origin()) <= 0){
            return direction.squared_length();
        }
        // add that new support point to simplex
        simplex.emplace_back(new_point);
    }
}

#endif  // ROBUST_IMPLICIT_SURFACE_NETWORKS_GJK
