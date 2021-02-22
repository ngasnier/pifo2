#include <stdexcept>
#include <math.h>
#include "Regridding.h"

namespace pifo {
    void Regridding::bilinearRegrid(double* x_in, long in_width, 
                        double* y_in, long in_height,
                        double* data_in, long size_in, 
                        bool cyclic, 
                        double* x_out, double* y_out, 
                        double* data_out, long size_out)
    {
        long* tab_i_in1 = new long[size_out];
        long* tab_i_in2 = new long[size_out];
        double* tab_x_adj1 = new double[size_out];
        double* tab_x_adj2 = new double[size_out];
        long* tab_j_in1 = new long[size_out];
        long* tab_j_in2 = new long[size_out];
        double* tab_y_adj1 = new double[size_out];
        double* tab_y_adj2 = new double[size_out];
        
        double x, y;

        long i_in1, j_in1;
        long i_in2, j_in2;

        double x_in1, y_in1;
        double x_in2, y_in2;

        double v1, v2, v3, v4;
        double vv1, vv2;
        double alpha_x, alpha_y;

        Regridding::optimizeGridIndices(x_in, in_width, x_out, size_out, cyclic, tab_i_in1, tab_i_in2, tab_x_adj1, tab_x_adj2);
        Regridding::optimizeGridIndices(y_in, in_height, y_out, size_out, false, tab_j_in1, tab_j_in2, tab_y_adj1, tab_y_adj2);
        
        for (long i=0;i<size_out;i++)
        {
            x = x_out[i];
            i_in1 = tab_i_in1[i];
            i_in2 = tab_i_in2[i];
            x_in1 = tab_x_adj1[i];
            x_in2 = tab_x_adj2[i];
            alpha_x = (x-x_in1)/(x_in2-x_in1);

            y = y_out[i];
            j_in1 = tab_j_in1[i];
            j_in2 = tab_j_in2[i];
            y_in1 = tab_y_adj1[i];
            y_in2 = tab_y_adj2[i];
            alpha_y = (y-y_in1)/(y_in2-y_in1);

            v1 = data_in[i_in1+in_width*j_in1];
            v2 = data_in[i_in1+in_width*j_in2];
            v3 = data_in[i_in2+in_width*j_in1];
            v4 = data_in[i_in2+in_width*j_in2];

            vv1 = alpha_y*v2 + (1-alpha_y)*v1;
            vv2 = alpha_y*v4 + (1-alpha_y)*v3;

            data_out[i] = alpha_x*vv2 + (1-alpha_x)*vv1 ;
        }

        delete tab_i_in1;
        delete tab_i_in2;
        delete tab_x_adj1;
        delete tab_x_adj2;
        delete tab_j_in1;
        delete tab_j_in2;
        delete tab_y_adj1;
        delete tab_y_adj2;
    }    
    
    
    // [ 0 1 2 ]      x
    // -------->
    // [ 2 1 0 ]
    void Regridding::optimizeGridIndices(double* x_in, long size_in, double* x_out, long size_out, bool cyclic, long* tab_i_in1, long* tab_i_in2, double* tab_x_adj1, double* tab_x_adj2)
    {
        if (size_in>1)
        {                
            long di_in = x_in[0]>x_in[size_in-1] ? -1 : 1;
            long i_in1 = 0;
            long i_in2 = 0;
            //double dx_start = (di_in>=0 ? x_in[1] - x_in[0] : x_in[size_in-2]-x_in[size_in-1]);
            double dx_end = (di_in>=0 ? x_in[size_in-1] - x_in[size_in-2] : x_in[0]-x_in[1]);
            double x_min = (di_in>=0 ? x_in[0] : x_in[size_in-1]);
            double x_max = (di_in>=0 ? x_in[size_in-1] : x_in[0]);
            //double x_start_cycle = (di_in>=0 ? x_in[0]-dx_start : x_in[size_in-1]-dx_start);
            double x_end_cycle = (di_in>=0 ? x_in[size_in-1]+dx_end : x_in[0]+dx_end);
            double cycle_length = x_end_cycle - x_min;
            double renorm = 0;
            
            double x = 0;

            for (long k=0;k<size_out;k++)
            {
                i_in1 = (di_in>=0 ? 0 : size_in-1);
                i_in2 = (di_in>=0 ? 1 : size_in-2);

                x = x_out[k];
                renorm = 0;
                if (!cyclic && (x<x_min || x>x_max)) 
                {
                    if (x<x_min)
                    {
                        if (di_in>=0)
                        {
                            i_in1 = 0;
                            i_in2 = 1;
                            tab_i_in1[k] = i_in1;
                            tab_i_in2[k] = i_in1;
                        }
                        else
                        {
                            i_in1 = size_in-1;
                            i_in2 = i_in1-1;
                            tab_i_in1[k] = i_in1;
                            tab_i_in2[k] = i_in1;
                        }
                    }
                    else
                    {
                        if (di_in>=0)
                        {
                            i_in1 = size_in-2;
                            i_in2 = i_in1+1;
                            tab_i_in1[k] = i_in2;
                            tab_i_in2[k] = i_in2;
                        }
                        else
                        {
                            i_in1 = 1;
                            i_in2 = i_in1-1;
                            tab_i_in1[k] = i_in2;
                            tab_i_in2[k] = i_in2;
                        }
                    }
                    tab_x_adj1[k] = x_in[i_in1];
                    tab_x_adj2[k] = x_in[i_in2];
                }
                else
                {
                    if ((x<x_min) || (x>x_end_cycle))
                    {
                        renorm = floor((x-x_min)/cycle_length);
                    }
                    x = x - renorm*cycle_length;

                    if (x>x_max && x<=x_end_cycle)
                    {
                        if (di_in>=0)
                        {
                            i_in1 = size_in-1;
                            i_in2 = 0;
                        }
                        else
                        {
                            i_in1 = 0;
                            i_in2 = size_in-1;
                        }
                        tab_x_adj1[k] = x_in[i_in1]+renorm*cycle_length;
                        tab_x_adj2[k] = x_end_cycle+renorm*cycle_length;
                    }
                    else
                    {
                        while ((x<x_in[i_in1] || x>x_in[i_in2]) 
                                && (i_in1>=0 && i_in1<size_in && i_in2>=0 && i_in2<size_in)) 
                        {
                            i_in1+=di_in;
                            i_in2+=di_in;
                        }
                        if (i_in1<0) i_in1 = size_in-1;
                        if (i_in2<0) i_in2 = size_in-1;
                        if (i_in1>size_in) i_in1 = 0;
                        if (i_in2>size_in) i_in2 = 0;
                        tab_x_adj1[k] = x_in[i_in1]+renorm*cycle_length;
                        tab_x_adj2[k] = x_in[i_in2]+renorm*cycle_length;
                    }
                    tab_i_in1[k] = i_in1;
                    tab_i_in2[k] = i_in2;
                }
            }
        }
        else
        {
            throw std::runtime_error("not enough coordinates to regrid.");
        }
    }
}