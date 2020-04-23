#ifndef _WAVE_OFFSITE_H
#define _WAVE_OFFSITE_H
//Implement various kinds of stencils required in the SeASiTe project.
//If new stencils are required update this file
#include "yask_compiler_api.hpp"
using namespace std;
using namespace yask;

namespace {
    class WaveBase_offsite : public yc_solution_with_radius_base {

        protected:

            // Indices & dimensions.
            yc_index_node_ptr t = new_step_index("t");           // step in time dim.
            yc_index_node_ptr x = new_domain_index("x");         // spatial dim.
            yc_index_node_ptr y = new_domain_index("y");         // spatial dim.
            yc_index_node_ptr z = new_domain_index("z");         // spatial dim.

            // Vars.
            yc_var_proxy data =
                yc_var_proxy("data", get_soln(), { t, x, y, z }); // time-varying 3D var.

            yc_var_proxy vel  =
                yc_var_proxy("vel", get_soln(), {t, x, y, z }); // 3D var.


            yc_var_proxy prev_data  =
                yc_var_proxy("prev_data", get_soln(), {t, x, y, z }); // 3D var.


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

            WaveBase_offsite(const string& name, int radius=4) : yc_solution_with_radius_base(name, radius) {
                //can support max. upto 3D space dimension
            }

            // Define equation at t+1 based on values at t.
            virtual void define() {

                // A condition that is true when index 'x','y'or 'z' is at boundary.
                // Normally bounding box calculation doesn't go to halo, but I hack
                // it so that the same stencil can be used for all vector foldings
                // without the problem of unimplemented unaligned writes
                const yc_bool_node_ptr& at_first_x = (x < first_domain_index(x));
                const yc_bool_node_ptr& at_first_y = (y < first_domain_index(y));
                const yc_bool_node_ptr& at_first_z = (z < first_domain_index(z));

                const yc_bool_node_ptr& at_last_x = (x > last_domain_index(x));
                const yc_bool_node_ptr& at_last_y = (y > last_domain_index(y));
                const yc_bool_node_ptr& at_last_z = (z > last_domain_index(z));

                boundaryHandling(data, t,x,y,z, at_first_x, at_last_x, at_first_y, at_last_y, at_first_z, at_last_z);
                yc_number_node_ptr v;
                points(v, data, t, x, y, z);
                // define the value at t+1 to be equivalent to v.
                double c_0 = 0.1;
                yc_number_node_ptr tmp;
                tmp = c_0*data(t,x,y,z) + v;
                data(t+1, x, y, z) EQUALS 2*data(t,x,y,z) - prev_data(t,x,y,z) + vel(t,x,y,z)*tmp;
            }
    };

    class Wave1dStencil_offsite : public WaveBase_offsite {
        protected:
            virtual double coeff(int dr) const
            {
                return 1/(double)dr;
            }

            // Add additional points to v.
            virtual void points(yc_number_node_ptr& v, yc_var_proxy& data, const yc_number_node_ptr& t, const yc_number_node_ptr& x, const yc_number_node_ptr& y, const yc_number_node_ptr& z)
            {
                for(int r=get_radius(); r >= 1; --r)
                {
                    v += coeff(r) * (data(t, x, y, z+r  ) +  data(t, x, y, z-r  ));
                }

            }

            virtual void addPoints(yc_number_node_ptr& v, yc_var_proxy& data, const yc_number_node_ptr& t, const yc_number_node_ptr& x, const yc_number_node_ptr& y, const yc_number_node_ptr& z)
            {
                points(v, data, t, x, y, z);
                double c_0 = 0.1;
                yc_number_node_ptr tmp;
                tmp = c_0*data(t,x,y,z) + v;
                v += 2*data(t,x,y,z) - prev_data(t,x,y,z) + vel(t,x,y,z)*tmp;

            }
            //Neumann boundary of first order
            virtual void boundaryHandling(yc_var_proxy& data, const yc_number_node_ptr& t, const yc_number_node_ptr& x, const yc_number_node_ptr& y, const yc_number_node_ptr& z, const yc_bool_node_ptr& at_first_x, const yc_bool_node_ptr& at_last_x, const yc_bool_node_ptr& at_first_y, const yc_bool_node_ptr& at_last_y, const yc_bool_node_ptr& at_first_z, const yc_bool_node_ptr& at_last_z)
            {
                //I do a hack here because unary '+' is not overloaded, adding
                //another term data(t,x,y,z) is also a hack to show dependency to
                //stencil compiler
                data(t, x, y, z) EQUALS (data(t-1, x, y, z+1)*1.0)
                    IF_DOMAIN at_first_z;

                //I do a hack here because unary '+' is not overloaded
                data(t, x, y, z) EQUALS (data(t-1, x, y, z-1)*1.0)
                    IF_DOMAIN at_last_z;
            }


        public:
            Wave1dStencil_offsite(int radius=4) :
                WaveBase_offsite("Wave1D_offsite", radius) { }
            Wave1dStencil_offsite(const string& name, int radius=4) :
                WaveBase_offsite(name, radius) { }

    };
    YS_REGISTER_STENCIL_RADIUS("Wave1D_offsite",Wave1dStencil_offsite);

    class Wave2dStencil_offsite : public Wave1dStencil_offsite {
        protected:
            virtual double coeff(int dr) const
            {
                return 1/(double)dr;
            }

            // Add additional points to v.
            virtual void points(yc_number_node_ptr& v, yc_var_proxy& data, const yc_number_node_ptr& t, const yc_number_node_ptr& x, const yc_number_node_ptr& y, const yc_number_node_ptr& z)
            {
                for(int r=get_radius(); r>=1; --r)
                {
                    v += coeff( r ) * (data(t, x, y+r, z  ) + data(t, x, y-r, z  ));
                }
                Wave1dStencil_offsite::points(v, data, t, x, y, z);
            }

            virtual void addPoints(yc_number_node_ptr& v, yc_var_proxy& data, const yc_number_node_ptr& t, const yc_number_node_ptr& x, const yc_number_node_ptr& y, const yc_number_node_ptr& z)
            {
                points(v, data, t, x, y, z);
                double c_0 = 0.1;
                yc_number_node_ptr tmp;
                tmp = c_0*data(t,x,y,z) + v;
                v += 2*data(t,x,y,z) - prev_data(t,x,y,z) + vel(t,x,y,z)*tmp;
            }

            //Neumann boundary of first order
            virtual void boundaryHandling(yc_var_proxy& data, const yc_number_node_ptr& t, const yc_number_node_ptr& x, const yc_number_node_ptr& y, const yc_number_node_ptr& z, const yc_bool_node_ptr& at_first_x, const yc_bool_node_ptr& at_last_x, const yc_bool_node_ptr& at_first_y, const yc_bool_node_ptr& at_last_y, const yc_bool_node_ptr& at_first_z, const yc_bool_node_ptr& at_last_z)
            {
                    Wave1dStencil_offsite::boundaryHandling(data, t,x,y,z,at_first_x,at_last_x,at_first_y,at_last_y,at_first_z,at_last_z);

                      data(t, x, y, z) EQUALS (data(t-1,x,y+1,z)*1.0)
                      IF_DOMAIN at_first_y;

                      data(t, x, y, z) EQUALS (data(t-1,x,y-1,z)*1.0)
                      IF_DOMAIN at_last_y;
            }


        public:
            Wave2dStencil_offsite(int radius=4) :
                Wave1dStencil_offsite("Wave2D_offsite", radius) { }
            Wave2dStencil_offsite(const string& name, int radius=4) :
                Wave1dStencil_offsite(name, radius) { }

    };
    YS_REGISTER_STENCIL_RADIUS("Wave2D_offsite",Wave2dStencil_offsite);


    class Wave3dStencil_offsite : public Wave2dStencil_offsite {
        protected:
            virtual double coeff(int dr) const
            {
                return 1/(double)dr;
            }

            // Add additional points to v.
            virtual void points(yc_number_node_ptr& v, yc_var_proxy& data, const yc_number_node_ptr& t, const yc_number_node_ptr& x, const yc_number_node_ptr& y, const yc_number_node_ptr& z)
            {
                for(int r=get_radius(); r>=1; --r)
                {
                    v += coeff( r ) * (data(t, x+r, y, z ) + data(t, x-r, y, z ));
                }
                Wave2dStencil_offsite::points(v, data, t, x, y, z);
            }

            virtual void addPoints(yc_number_node_ptr& v, yc_var_proxy& data, const yc_number_node_ptr& t, const yc_number_node_ptr& x, const yc_number_node_ptr& y, const yc_number_node_ptr& z)
            {
                points(v, data, t, x, y, z);
                double c_0 = 0.1;
                yc_number_node_ptr tmp;
                tmp = c_0*data(t,x,y,z) + v;
                v += 2*data(t,x,y,z) - prev_data(t,x,y,z) + vel(t,x,y,z)*tmp;
            }

            //Neumann boundary of first order
            virtual void boundaryHandling(yc_var_proxy& data, const yc_number_node_ptr& t, const yc_number_node_ptr& x, const yc_number_node_ptr& y, const yc_number_node_ptr& z, const yc_bool_node_ptr& at_first_x, const yc_bool_node_ptr& at_last_x, const yc_bool_node_ptr& at_first_y, const yc_bool_node_ptr& at_last_y, const yc_bool_node_ptr& at_first_z, const yc_bool_node_ptr& at_last_z)
            {
/*                   Wave2dStencil_offsite::boundaryHandling(data, t,x,y,z,at_first_x,at_last_x,at_first_y,at_last_y,at_first_z,at_last_z);
                   data(t, x, y, z) EQUALS (data(t-1,x+1,y,z)*1.0)
                   IF_DOMAIN at_first_x;
                   data(t, x, y, z) EQUALS (data(t-1,x-1,y,z)*1.0)
                   IF_DOMAIN at_last_x;
                   */
            }


        public:
            Wave3dStencil_offsite(int radius=4) :
                Wave2dStencil_offsite("Wave3D_offsite", radius) { }
            Wave3dStencil_offsite(const string& name, int radius=4) :
                Wave2dStencil_offsite(name, radius) { }

    };

    YS_REGISTER_STENCIL_RADIUS("Wave3D_offsite",Wave3dStencil_offsite);
}

#endif
