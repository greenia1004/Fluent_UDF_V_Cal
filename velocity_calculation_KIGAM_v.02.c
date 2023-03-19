#include "udf.h"
#include "metric.h"

DEFINE_ON_DEMAND(on_demand_calc)
{
    #define NMAX 1000          // Maximum array size
    #if !RP_HOST
    Domain *d;
    Thread *t;
    
    /* Integrate Velocity Magnitude */

    real sum_vol_vel;          // volume velocity sum
    real sum_surf_vel;         // surface velocity sum
    real sum_vol;              // volume sum
    real sum_surf;             // surface sum
    real avg_vol_velA[NMAX];   // array of the average of cell volume weighted velocity
    real avg_surf_velA[NMAX];  // array of the average of surface area weighted velocity
    
    real x[ND_ND];             // cartesian coordinate x[0] -> x, x[1] -> y, x[2] -> z
    real A[3];                 // Area vector
    real radius_i;             // radius of inner cylinder: 0.075 [m]
    real radius_o;             // radius of outer cylinder: 0.0755 [m]
    real length;               // lenght of cylinders: 0.01 [m]
    real gap;                  // gap distance between inner and outer cylinders: gap = radius_o - radius_i
    real delta_gap;            // gap interval distance delta_gap = gap / ndiv
    int ndiv_gap;              // number to divide the gap
    //int ndiv_len;              // number to divide the length
    
    real radiusA[NMAX];         // radius array of iso-surfaces
    real delta_r = 2.0e-6/2.0; // delta r of iso-surface: half of mesh defeature size
    
    real rho;                  // rho at cylindrical coordinates (rho, phi, z)
    real phi;                  // phi at cylindrical coordinates (rho, phi, z)
    real z;                    // z at cartesian and cylindrical coorinates
    
    /* Initializing the variables */
    radius_i = 0.075;
    radius_o = 0.0755;
    length = 0.01;

    ndiv_gap = 100;
    //ndiv_len = 100;
    gap = radius_o - radius_i;
    delta_gap = gap / ndiv_gap;
    radius[0] = radius_i;                  // for inner cylinder
    radius[ndiv_gap-1] = radius_o;         // for outer cylinder
    for(int i=1; i < (ndiv_gap-2); i++)    // from 1 to 98
    {
        radius[i] = radius[0] + delta_gap * (real) i;
    }
    
    //real radial_v;
    //real tangetial_v;
    cell_t c;
    d = Get_Domain(1);

    thread_loop_c(t, d)
    {
        if (FLUID_THREAD_P(t))
        {
            for(i=0; i< ndiv_gap; i++)  // from 0 to 99
            {
                radius = radius[i];
                
                sum_vol_vel = 0.0;
                sum_surf_vel = 0.0;
                sum_vol = 0.0;
                sum_surf = 0.0;
                
                begin_c_loop_int(c,t)
                {
                    C_CENTROID(x,c,t);
                    if (sqrt(pow(x[0],2) + pow(x[1],2)) >= radius && sqrt(pow(x[0],2) + pow(x[1],2)) < radius + delta_r )
                    {
                        //radial_v = sqrt(C_U(c,t)*C_U(c,t)+C_V(c,t)*C_V(c,t));
                        //tangential_v = atan2(C_W(c,t)/C_U(c,t));

                        sum_vol_vel += sqrt(C_U(c,t)*C_U(c,t)+C_V(c,t)*C_V(c,t))*C_VOLUME(c,t);
                        sum_vol += C_VOLUME(c,t);
                    }
                    else
                    {                        
                        sum_vol_vel += 0;
                        sum_vol += 0;
                    }
                }
                end_c_loop_int(c,t)
                avg_vol_velA[i] = sum_vol_vel/sum_vol;
                
                begin_f_loop_int(f,t)
                {
                    F_CENTROID(x,f,t);
                    F_AREA(A,f,t);
                    if (sqrt(pow(x[0],2) + pow(x[1],2)) >= radius && sqrt(pow(x[0],2) + pow(x[1],2)) < radius + delta_r )
                    {
                        //radial_v = sqrt(C_U(c,t)*C_U(c,t)+C_V(c,t)*C_V(c,t));
                        //tangential_v = atan2(C_W(c,t)/C_U(c,t));

                        sum_surf_vel += sqrt(C_U(c,t)*C_U(c,t)+C_V(c,t)*C_V(c,t))*NV_MAG(A);
                        sum_surf += NV_MAG(A);
                    }
                    else
                    {                        
                        sum_surf_vel += 0;
                        sum_surf += 0;
                    }
                }
                end_f_loop_int(f,t)
                avg_surf_velA[i] = sum_surf_vel/sum_surf;
            }
        }

        #if RP_NODE
        avg_vol_vel = PRF_GRSUM1(avg_vol_vel);
        avg_surf_vel = PRF_GRSUM1(avg_surf_vel);
        #endif
    }

    #if PARALLEL
    if(I_AM_NODE_ZERO_P)
    #endif
    printf("== nnnn ==  Volume Averaged Velocity ==\n");
    for(i=0; i< ndiv_gap; i++)  // from 0 to 99
    {
        printf("%4d\t%g\n", i, avg_vol_velA[i]);
    }
    fflush(stdout);

    #endif
}
