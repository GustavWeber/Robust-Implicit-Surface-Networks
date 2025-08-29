#ifndef ROBUST_IMPLICIT_SURFACE_NETWORK_APPROX_TET
#define ROBUST_IMPLICIT_SURFACE_NETWORK_APPROX_TET

#include <CGAL/Kd_tree_rectangle.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Kernel/global_functions_3.h>
#include <CGAL/squared_distance_3.h>

namespace ApproxTetDist{

template<typename Kernel>
void GetBoundingSphere(const typename Kernel::Tetrahedron_3& tet, typename Kernel::Point_3& center, typename Kernel::FT& radius){
    typedef typename Kernel::Point_3 Point;
    typedef typename Kernel::FT FT;
    typedef typename Kernel::Plane_3 Plane;

    Point circumcenter = CGAL::circumcenter(tet);
    std::array<Point, 4> on_sphere;
    int active = 0;
    for(int i = 0; i<4; i++){
        Plane plane(tet[i+1], tet[i+2], tet[i+3]);
        if(plane.oriented_side(circumcenter) == plane.oriented_side(tet[i])){
            // If the circumcenter is on the inside of the opposing plane, the vertex is on the bounding sphere
            on_sphere[active++] = tet[i];
        }
    }
    switch (active) {
        case 2:
            center = CGAL::circumcenter(on_sphere[0], on_sphere[1]);
            break;
        case 3:
            center = CGAL::circumcenter(on_sphere[0], on_sphere[1], on_sphere[2]);
            break;
        case 4:
            center = CGAL::circumcenter(on_sphere[0], on_sphere[1], on_sphere[2], on_sphere[3]);
            break;
        default:
            assert(false);
            //making sure to not have uninitialized output data but this should never happen
            radius = -1000;
            center = CGAL::Origin();
            return;
    }
    radius = CGAL::sqrt(CGAL::squared_distance(center, on_sphere[0]));

}

template<typename Kernel>
typename Kernel::FT GetDistance_Box(const typename Kernel::Tetrahedron_3& tet, const CGAL::Kd_tree_rectangle<typename Kernel::FT, CGAL::Dimension_tag<3>>& rect){
typedef typename Kernel::FT FT;
    const CGAL::Bbox_3& bbox = tet.bbox();
    FT distance = 0;
    FT dim_dist;
    for(int dim = 0; dim<3; dim++){
        if(bbox.min(dim) > rect.max_coord(dim)){
            dim_dist = bbox.min(dim) - rect.max_coord(dim);
        } else if(bbox.max(dim) < rect.min_coord(dim)){
            dim_dist = rect.min_coord(dim) - bbox.max(dim);
        } else {
            dim_dist = 0;
        }
        distance += dim_dist*dim_dist;
    }
    return distance;
}

template<typename Kernel>
typename Kernel::FT GetDistance_Sphere(const typename Kernel::Tetrahedron_3& tet, const CGAL::Kd_tree_rectangle<typename Kernel::FT, CGAL::Dimension_tag<3>>& rect){
    typedef typename Kernel::FT FT;
    typedef typename Kernel::Point_3 Point;
    
    Point center;
    FT radius;

    GetBoundingSphere<Kernel>(tet, center, radius);
    
    // Get Distance between center and box 
    FT squared_distance = 0;
    for(int dim = 0; dim <3; dim++){
        FT coord_dist = center[dim] - std::clamp(center[dim], rect.min_coord(dim), rect.max_coord(dim));
        squared_distance+=coord_dist*coord_dist;
    }
    FT distance = CGAL::sqrt(squared_distance)-radius; // Sadly we cannot just use the squared radius
    return distance * distance;
}
};
#endif //ROBUST_IMPLICIT_SURFACE_NETWORK_APPROX_TET
