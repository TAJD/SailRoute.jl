 package  -------------------------------------------------------
 
   Performs several verification calculations given a file of grid spacings 
   and some observed quantity corresponding to each grid spacing.
 
   Computes:
   - order of convergence
   - Richardson extrapolation to zero grid spacing
   - grid convergence indices (GCI)
 
 --------------------------------------------------------------------------
 
Adapted from:
    
    NPARC Alliance CFD Verification and Validation Web Site
    Examining Spatial (Grid) Convergence
    verify.f90
    URL: http://www.grc.nasa.gov/WWW/wind/valid/tutorial/spatconv.html
    
    Nov '11: Updated to reflect Celik et al 2008.
    Sept '18: Translated to Julia

    Credit to Matthew Topper.
"""

""" Calculate the order of convergence values generated with three
grids of reducing resolution (ie grid_1 is finest). The values of the grids
are needed along with the ratios between them.

An iterative method with under-relaxation is used to calculate the order
of convergence as the refinement ratio is not necessarily constant.

This has been modified to the method of Celik (2008).
"""
function order_of_convergence(value_1, value_2, value_3, ratio_21, ratio_32, omega=0.5, tol=1.E-6)
    # Set a maximum residual and number of iterations
    max_res = 1.E6
    # calculate the epsilons.
    epsilon32 = float(value_3 - value_2)
    epsilon21 = float(value_2 - value_1)
    # Calculate the fraction
    epfrac = epsilon32 / epsilon21
    # Get the signed unit, s
    s = epfrac / abs(epfrac)
    # Initial guess at order of convergence, p
    p1 = (1. / log(ratio_21)) * abs(log(abs(epfrac))) # start_p
    # Initialise the residual and number of iterations
    residual = 1.0
    iterations = 0
    while abs(residual) > tol
        # Break if it's all gone bad
        if float(iterations) > max_res && residual > max_res
            println("Residual out of range or too many iterations")
        end
        # Get the last value
        p0 = p1
        # Calculate q
        q = log((ratio_21^p0 - s) / (ratio_32^p0 - s))
        # Calculate the p iteration
        pnew = (1. / log(ratio_21)) * abs(log(abs(epfrac)) + q)
        # Calculate the relaxation step.
        p1 = (1. - omega) * p0 + omega * pnew
        residual = p1 - p0 
        iterations += 1
    end
    # if abs(p1) > 1.0
    #     return 1.0
    # else
    return abs(p1)
    # end
end


""" Estimate the zero grid spacing value using richardsons extrapolation and
two grids of reducing resolution (ie grid_1 is finest). The refinement ratio
is needed. The order of convergence, p, is also required.
"""
function richardson_extrapolate(value_1, value_2, ratio_21, p)
    f_exact = ( ratio_21^p * value_1 - value_2 ) / ( ratio_21^p - 1.0 )
    return f_exact
end


""" This routine returns the relative error and extrapolated 
relative error. The values of the grids and needed along 
with the extrapolated value.
"""
function error_estimates(value_1, value_2, f_exact)
    # Get the approximate relative error
    e21a =  abs( (value_1 - value_2) / value_1 )
    # Get the extrapolated relative error
    e21ext = abs( ( f_exact - value_1 ) / f_exact )
    return e21a, e21ext
end


""" Calculate the fine and coarse grid convergence index for two grids of 
reducing resolution (ie grid_1 is finest). The refinement ration between the 
grids is required along with the approximate relative error (e21_approx) and 
the order of convergence, p.
"""
function gci(ratio_21, e21_approx, p)
    # Using a fixed safety factor as per Celik (2008)
    safety_factor = 1.25
    # Calculate the gci
    gci_fine = safety_factor * e21_approx / (ratio_21^p - 1.0)
    gci_coarse = ratio_21^p * gci_fine
    return gci_fine, gci_coarse
end


""" Calculate the ratio in succesive Eps as defined at the bottom of page
129 in Roache. If the ration is close to one then the asymptotic range has
been reached.
"""
function asymptotic_ratio(gci_fine_21, gci_fine_32, ratio_21, p)
    ratio = ratio_21^p * ( gci_fine_21 / gci_fine_32)
    return ratio
end


function apply_routine(print=false)
    value_1 = 0.970500
    value_2 = 0.968540
    value_3 = 0.961780
    h1 = 1.0
    h2 = 2.0
    h3 = 4.0
    ratio_21 = h2/h1
    ratio_32 = h3/h2
    ooc = order_of_convergence(value_1, value_2, value_3, ratio_21, ratio_32)
    f_exact_21 = richardson_extrapolate(value_1, value_2, ratio_21, ooc)
    f_exact_32 = richardson_extrapolate(value_2, value_3, ratio_32, ooc)
    e21a, e21ext = error_estimates(value_1, value_2, f_exact_21)
    e32a, e32ext = error_estimates(value_2, value_3, f_exact_32)
    gci_fine_21, gci_coarse_21 = gci(ratio_21, e21a, ooc)
    gci_fine_32, gci_coarse_32 = gci(ratio_32, e32a, ooc)
    ratio = asymptotic_ratio(gci_fine_21, gci_fine_32, ratio_21, ooc)

    if print==true
        println("Order of convergence: ", ooc)
        println("Richardson extrapolation: ", f_exact_21)
        println("Approximate error:  ", e21a)
        println("Extrapolated error: ", e21ext)
        println("Asymptotic ratio", ratio)
    end
    return ooc, f_exact_21, e21a, e21ext, gci_fine_21, gci_coarse_21, ratio
end


"""Return the fine GCI value"""
function GCI_calc(value_1, value_2, value_3, h1, h2, h3)
    ratio_21 = h2/h1
    ratio_32 = h3/h2
    ooc = order_of_convergence(value_1, value_2, value_3, ratio_21, ratio_32)
    f_exact_21 = richardson_extrapolate(value_1, value_2, ratio_21, ooc)
    f_exact_32 = richardson_extrapolate(value_2, value_3, ratio_32, ooc)
    e21a, e21ext = error_estimates(value_1, value_2, f_exact_21)
    e32a, e32ext = error_estimates(value_2, value_3, f_exact_32)
    gci_fine_21, gci_coarse_21 = gci(ratio_21, e21a, ooc)
    gci_fine_32, gci_coarse_32 = gci(ratio_32, e32a, ooc)
    ratio = asymptotic_ratio(gci_fine_21, gci_fine_32, ratio_21, ooc)
    return abs(gci_fine_21)
end


"""Return the fine GCI value"""
function extrap_value(value_1, value_2, value_3, h1, h2, h3)
    ratio_21 = h2/h1
    ratio_32 = h3/h2
    ooc = order_of_convergence(value_1, value_2, value_3, ratio_21, ratio_32)
    f_exact_21 = richardson_extrapolate(value_1, value_2, ratio_21, ooc)
    return abs(f_exact_21)
end


"""Return the order of convergence."""
function ooc_value(value_1, value_2, value_3, h1, h2, h3)
    ratio_21 = h2/h1
    ratio_32 = h3/h2
    ooc = order_of_convergence(value_1, value_2, value_3, ratio_21, ratio_32)
    return ooc
end
