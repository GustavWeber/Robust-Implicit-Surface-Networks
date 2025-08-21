#ifndef ROBUST_IMPLICIT_SURFACE_NETWORK_APPROX_TET
#define ROBUST_IMPLICIT_SURFACE_NETWORK_APPROX_TET

#include <CGAL/Kd_tree_rectangle.h>
#include <CGAL/Bbox_3.h>

namespace ApproxTetDist{

template<typename Kernel>
typename Kernel::FT GetDistance(const typename Kernel::Tetrahedron_3& tet, const CGAL::Kd_tree_rectangle<typename Kernel::FT, CGAL::Dimension_tag<3>>& rect){
typedef typename Kernel::FT FT;
    const CGAL::Bbox_3& bbox = tet.bbox();
    FT distance;
    FT dim_dist;
    for(int dim = 0; dim<3; dim++){
        if(bbox.min(dim) > rect.max_coord(dim)){
            dim_dist = bbox.min(dim) - rect.max_coord(dim);
        } else if(bbox.max(dim) < rect.min_coord(dim)){
            dim_dist = bbox.max(dim) - rect.max_coord(dim);
        } else {
            dim_dist = 0;
        }
        distance += dim_dist*dim_dist;
    }
    return distance;
}
};
#endif //ROBUST_IMPLICIT_SURFACE_NETWORK_APPROX_TET
