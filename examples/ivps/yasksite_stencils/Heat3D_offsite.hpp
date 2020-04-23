#ifndef _HEAT_OFFSITE_H
#define _HEAT_OFFSITE_H
//Implement various kinds of stencils required in the SeASiTe project.
//If new stencils are required update this file
#include "yask_compiler_api.hpp"
using namespace std;
using namespace yask;

namespace {
    class HeatBase_offsite : public yc_solution_with_radius_base {
        protected:
            // Indices & dimensions.
            yc_index_node_ptr t = new_step_index("t");           // step in time dim.
            yc_index_node_ptr x = new_domain_index("x");         // spatial dim.
            yc_index_node_ptr y = new_domain_index("y");         // spatial dim.
            yc_index_node_ptr z = new_domain_index("z");         // spatial dim.
            // Vars.
            yc_var_proxy data = yc_var_proxy("data", get_soln(), { t, x, y, z }); // time-varying 3D var.
            yc_var_proxy h_inv = yc_var_proxy("h_inv", get_soln()); // time-varying 3D var.

            // Return a coefficient.  Note: This returns completely fabricated
            // values only for illustrative purposes; they have no mathematical
            // significance.
            virtual double coeff(int di, int dj, int dk) const {
                return 1;
            }

            virtual void points(yc_number_node_ptr& v, yc_var_proxy& grid, const yc_number_node_ptr& t, const yc_number_node_ptr& x, const yc_number_node_ptr& y, const yc_number_node_ptr& z) =0;
            virtual void addPoints(yc_number_node_ptr& v, yc_var_proxy& grid, const yc_number_node_ptr& t, const yc_number_node_ptr& x, const yc_number_node_ptr& y, const yc_number_node_ptr& z) =0;
            virtual void boundaryHandling(yc_var_proxy& grid, const yc_number_node_ptr& t, const yc_number_node_ptr& x, const yc_number_node_ptr& y, const yc_number_node_ptr& z, const yc_bool_node_ptr& at_first_x, const yc_bool_node_ptr& at_last_x, const yc_bool_node_ptr& at_first_y, const yc_bool_node_ptr& at_last_y, const yc_bool_node_ptr& at_first_z, const yc_bool_node_ptr& at_last_z) =0;

        public:
            HeatBase_offsite(const string& name, int radius=1) : yc_solution_with_radius_base(name, radius) {
                //can support max. upto 3D space dimension
            }

            // Define equation at t+1 based on values at t.
            virtual void define() {
                // A condition that is true when index 'x','y'or 'z' is at boundary.
                // Normally bounding box calculation doesn't go to halo, but I hack
                // it so that the same stencil can be used for all vector foldings
                // without the problem of unimplemented unaligned writes
                auto at_first_x = (x < first_domain_index(x));
                auto at_first_y = (y < first_domain_index(y));
                auto at_first_z = (z < first_domain_index(z));

                auto at_last_x = (x > last_domain_index(x));
                auto at_last_y = (y > last_domain_index(y));
                auto at_last_z = (z > last_domain_index(z));

                boundaryHandling(data, t, x, y, z, at_first_x, at_last_x, at_first_y, at_last_y, at_first_z, at_last_z);
                yc_number_node_ptr v;
                addPoints(v, data, t, x, y, z);
                // define the value at t+1 to be equivalent to v.
                data(t+1, x, y, z) EQUALS v;
            }
    };

    class Heat3dStencil_offsite : public HeatBase_offsite {
        protected:
            virtual double coeff(int di, int dj, int dk) const
            {
                if((dj==0&&dk==0)&&(di==0))
                {
                    return -6;
                }
                else
                {
                    return 1;
                }
            }

            // Add additional points to v.
            virtual void points(yc_number_node_ptr& v, yc_var_proxy& data, const yc_number_node_ptr& t, const yc_number_node_ptr& x, const yc_number_node_ptr& y, const yc_number_node_ptr& z)
            {
                // points in x-dimension
                for(int r=get_radius(); r>=1; --r)
                {
                    v += coeff(r, 0, 0) * data(t, x+r, y, z) + coeff(-r, 0, 0) * data(t, x-r, y, z);
                }

                // points in y-dimension
                for(int r=get_radius(); r>=1; --r)
                {
                    v += coeff(0, r, 0) * data(t, x, y+r, z) + coeff( 0,-r, 0) * data(t, x, y-r, z);
                }

                // points in z-dimension
                for(int r=get_radius(); r >= 1; --r)
                {
                    v += coeff(0, 0, r) * data(t, x, y, z+r) + coeff(0, 0, -r) * data(t, x, y, z-r);
                }

                v += coeff(0, 0, 0) * data(t, x, y , z);
            }

            virtual void addPoints(yc_number_node_ptr& v, yc_var_proxy& data, const yc_number_node_ptr& t, const yc_number_node_ptr& x, const yc_number_node_ptr& y, const yc_number_node_ptr& z)
            {
                points(v, data, t, x, y, z);
                v *= h_inv;
            }

            // Neumann boundary conditions
            virtual void boundaryHandling(yc_var_proxy& data, const yc_number_node_ptr& t, const yc_number_node_ptr& x, const yc_number_node_ptr& y, const yc_number_node_ptr& z, const yc_bool_node_ptr& at_first_x, const yc_bool_node_ptr& at_last_x, const yc_bool_node_ptr& at_first_y, const yc_bool_node_ptr& at_last_y, const yc_bool_node_ptr& at_first_z, const yc_bool_node_ptr& at_last_z)
            {
                // boundaries in z-dimension
                //I do a hack here because unary '+' is not overloaded, adding
                //another term data(t,x,y,z) is also a hack to show dependency to
                //stencil compiler
               /* data(t, x, y, z) EQUALS (data(t-1, x, y, z) * 1.0)    IF_DOMAIN at_first_z;
                //I do a hack here because unary '+' is not overloaded
                data(t, x, y, z) EQUALS (data(t-1, x, y, z) * 1.0)    IF_DOMAIN at_last_z;

                // boundaries in y-dimension
                data(t, x, y, z) EQUALS (data(t-1, x, y, z) * 1.0)    IF_DOMAIN at_first_y;
                data(t, x, y, z) EQUALS (data(t-1, x, y, z) * 1.0)    IF_DOMAIN at_last_y;

                // boundaries in x-dimension
                data(t, x, y, z) EQUALS (data(t-1, x, y, z) * 1.0)    IF_DOMAIN at_first_x;
                data(t, x, y, z) EQUALS (data(t-1, x, y, z) * 1.0)    IF_DOMAIN at_last_x;
                */
            }

        public:
            Heat3dStencil_offsite(int radius=1) :
                HeatBase_offsite("Heat3D_offsite", radius) { }
            Heat3dStencil_offsite(const string& name, int radius=4) :
                HeatBase_offsite(name, radius) { }
    };

    YS_REGISTER_STENCIL_RADIUS("Heat3D_offsite", Heat3dStencil_offsite);
}
#endif
