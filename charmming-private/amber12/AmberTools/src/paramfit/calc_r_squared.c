/*****************************************************
 * AMBER Bond Angle and Dihedral Parameter Optimiser *
 *                                                   *
 *           Written by: Robin Betz  (2011)          *
 *                       Ross Walker (2004)          *
 *                   UC San Diego                    *
 *           San Diego Supercomputer Center          *
 *            La Jolla, California, 92092            *
 *                       USA                         *
 *****************************************************/

/*calc_r_squared.c*/

/*The routines here are used to calculate the R squared value
  for a data fit

  R^2 here is calculated from the Pearson correlation

                         N*Sum(x[i]*y[i]) - (Sum(x[i]) * Sum(y[i])
  R =  -------------------------------------------------------------------------
       Sqrt[N*Sum(x[i]^2) - (Sum(x[i]))^2] * Sqrt[N*Sum(y[i]^2) - (Sum(y[i]))^2]

  R^2 is defined as R*R

*/

#include <stdio.h>
#include <math.h>

#include "constants.h"
#include "function_def.h"

double calc_r_squared(global_options_struct *global_options, parm_struct *parm_data, coords_struct *coords_data)
{
  /*Here, x are our fit and y are our QM data*/
  double sum_xy;
  double sum_x;
  double sum_y;
  double x_value, y_value;
  double sum_of_x_squareds;
  double sum_of_y_squareds;
  double r;
  double r_squared;
  int point;

  x_value=0.0;
  y_value=0.0;
  sum_xy=0.0;
  sum_x=0.0;
  sum_y=0.0;
  sum_of_x_squareds=0.0;
  sum_of_y_squareds=0.0;
  r=0.0;
  r_squared=0.0;
    
  /*The first stage is to calculate the sum x (sum of our fitted points), sum y (sum of qm data) and sum x * y*/
  /*Loop over all data points = NSTRUCTURES*/
  /*QM energies are stored internally in kcal/mol*/
  for (point=0; point<global_options->NSTRUCTURES; ++point)
  {
    y_value=coords_data[point].energy;
    if (global_options->FUNC_TO_FIT==SUM_SQUARES_AMBER_STANDARD)
    {
      x_value=eval_amber_std_for_single_struct(global_options, parm_data, &coords_data[point]);
    }
    else /*UNKNOWN FUNCTION - for the moment just set the value to zero*/
    {
      x_value=0.0;
    }
    /*Note we have to add our offset K to the x_value*/
    x_value+=global_options->K;
    sum_xy+=(x_value*y_value);
    sum_x+=x_value;
    sum_y+=y_value;
    sum_of_x_squareds+=(x_value*x_value);
    sum_of_y_squareds+=(y_value*y_value);
  }
  
  if (global_options->VERBOSITY>=HIGH)
  {
     printf("Sum of QM data points x MM data points(-K)  = %12.8E kCal^2/mol^2\n",sum_xy);
     printf("Sum of QM data points   = %12.8E kCal/mol, MM data points(-K)   = %12.8E kCal/mol\n",sum_y,sum_x);
     printf("Sum of QM data points^2 = %12.8E kCal/mol, MM data points(-K)^2 = %12.8E kCal/mol\n",sum_of_y_squareds,sum_of_x_squareds);
  }

  /*We now have all the data to calculate R*/

  r = (global_options->NSTRUCTURES * sum_xy) - (sum_x * sum_y);

  r = r / (     sqrt( (global_options->NSTRUCTURES * sum_of_x_squareds) - (sum_x * sum_x) )
              * sqrt( (global_options->NSTRUCTURES * sum_of_y_squareds) - (sum_y * sum_y) )
           );
           
  r_squared = r * r;

  if (r < 0 )
  {
     printf("!  WARNING negative correlation detected - R = %12.8E\n",r);
  }
            
  return r_squared;

}
