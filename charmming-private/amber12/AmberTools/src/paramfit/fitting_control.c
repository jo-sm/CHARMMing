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

/* fitting_control.c

   contains routines to control the fitting process calling
   the relevant energy minimisation routines etc.
*/

#include <stdio.h>

#include "function_def.h"

int do_fit(global_options_struct *global_options, parm_struct *parm_data, coords_struct *coords_data)
{	
   double initial_function_value;
   double initial_r_squared_value;
   double initial_energy_struc_1;
   double final_r_squared_value;
   double final_function_value;
   double final_energy_struc_1;
   int retval;
   bounds_struct bounds_data;
   
   /* Get initial information about the structures and ensure that if any dihedral phases are being *
    * parameterised that the phases are well represented in the sample structures.                  */
   if (global_options->CHECK_BOUNDS != NO)
   {
     retval = calculate_structure_diversity(global_options, &bounds_data, parm_data, coords_data);
     if (retval != SUCCESS)
       return retval;
     retval = check_dihedrals(global_options, parm_data, &bounds_data);
     if (retval != SUCCESS)
       return retval;
   }
   
   /*The initial parameters and details on the fitting etc will have been listed by the options_summary
     routine. So the first thing we do here is do an initial evaluation of our function*/
   if (global_options->FUNC_TO_FIT==SUM_SQUARES_AMBER_STANDARD)
   {
      if (global_options->VERBOSITY >= MEDIUM || global_options->RANDOM_SEED != 0) // don't print for test cases
      {
        printf("      --------------------------------- INITIAL PARAMETERS ------------------------------\n");
        printf("   Initial ");
        print_parameter_summary(global_options,parm_data);
        printf("      -----------------------------------------------------------------------------------\n");

      }
      
      /* Calculate how many bonds, angles, and dihedrals will be fit. */
      global_options->BOND_PARAMS = calculate_no_fit_params(parm_data, BONDS);
      global_options->ANGLE_PARAMS = calculate_no_fit_params(parm_data, ANGLES);
      global_options->DIHEDRAL_PARAMS = calculate_no_fit_params(parm_data, DIHEDRALS);
      
      initial_function_value = eval_sum_squares_amber_std(global_options, parm_data, coords_data);

      /*Calculate initial R^2 value*/
      initial_r_squared_value = calc_r_squared(global_options, parm_data, coords_data);
      initial_energy_struc_1=eval_amber_std_for_single_struct(global_options, parm_data, &coords_data[0]);
      if (global_options->VERBOSITY >= MEDIUM)
      {
        printf("   Sum of squares for initial parameters = %15.10f kcal^2/mol^2\n",initial_function_value);
        printf("   R^2 value for initial parameters      = %10.6f\n",initial_r_squared_value);
        printf("   Calculated energy with initial parameters for structure 1 = %10.6f KCal/mol\n\n",initial_energy_struc_1);
      }
   }
   else
   {
     printf("   ERROR IN DO_FIT - FUNCTION %d IS NOT YET IMPLEMENTED\n",global_options->FUNC_TO_FIT);
     return NOT_IMPLEMENTED;
   }
   
   // If K is being fit along with other parameters, give an error
   if (global_options->K_FIT==YES && global_options->PARAMETERS_TO_FIT!=K_ONLY)
   {
     printf("   ERROR: K is being fit along with other parameters. This will produce inaccurate\n");
     printf("         results and is disallowed. Specify a value for K or set PARAMETERS_TO_FIT=K_ONLY\n");
   }
 
   if (global_options->FITTING_FUNCTION==SIMPLEX)
      retval=minimise_function_simplex(global_options, parm_data, coords_data);
   
   else if (global_options->FITTING_FUNCTION==GENETIC)
   {
     if (global_options->VERBOSITY>=HIGH)
      print_parameter_summary(global_options,parm_data);
     
      retval = minimise_function_genetic(global_options, parm_data, coords_data);
   }
   else if (global_options->FITTING_FUNCTION==BOTH)
   {
     retval = minimise_function_genetic(global_options, parm_data, coords_data);
     retval = minimise_function_simplex(global_options, parm_data, coords_data);
   }
   else if (global_options->FITTING_FUNCTION==NONE)
   {
     printf("! Using no fitting function- returning initial parameters.\n");
   }
   else
   {
     printf("   ERROR IN do_fit() - UNKNOWN FITTING FUNCTION: %d\n",global_options->FITTING_FUNCTION);
     return UNKNOWN_OPT;
   }
   if (retval!=SUCCESS )
   {
      /*We can trap one error here - that is if we quit due to exceeding maximum iterations*/
      /*Other errors are fatal*/
      if (retval==EXCEEDEDMAXITERATIONS || retval==MINSTATIC)
        printf("!  Warning - final parameters do NOT represent a converged fit.\n");
      else
        return(retval);
   }

/*Now we shall attempt to calculate and R^2 value for our fit*/
      final_function_value = eval_sum_squares_amber_std(global_options, parm_data, coords_data);

    /* Check that the results are reasonable based on the given data in the input structures.       *
    * Dihedrals are already checked at the beginning, but do bonds and angles here to make sure    *
    * the result is in a well-represented area in the input data.                                  */
    if (global_options->CHECK_BOUNDS)
    {
      retval = check_angles(global_options, parm_data, &bounds_data);
      retval = check_bonds(global_options, parm_data, &bounds_data);
      clean_up_bounds(&bounds_data);
      if (retval != SUCCESS) return FAILURE;
    }
    
    final_r_squared_value=calc_r_squared(global_options, parm_data, coords_data);
    final_energy_struc_1=eval_amber_std_for_single_struct(global_options, parm_data, &coords_data[0]);

    if (global_options->VERBOSITY>=MEDIUM)
    {
      printf("   ---------------------------------- FINAL PARAMETERS -------------------------------\n");
    printf("   Fitted ");
    print_parameter_summary(global_options,parm_data);
    printf("   -----------------------------------------------------------------------------------\n");
    printf("   Function value with fitted parameters  =  %12.4f, R^2 = %12.4f\n",final_function_value,final_r_squared_value);
    printf("   Calculated energy with fitted parameters for structure 1 = %11.4f KCal/mol\n",final_energy_struc_1);
    printf("\n");
    }
    
    // if desired, save the frcmod file
    if (global_options->WRITE_FRCMOD)
      write_frcmod(global_options, parm_data);
    
    // if desired, write the list of final energies for comparison
    if (global_options->WRITE_ENERGY)
      write_energy(global_options, parm_data, coords_data, -1);
    
    // if desired, save a new toplogy and parameter file
    if (global_options->WRITE_PRMTOP)
      write_prmtop(global_options, parm_data);

   return SUCCESS;
}


