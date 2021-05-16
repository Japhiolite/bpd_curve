# -*- coding: utf-8 -*-
"""
Created on Thu May 13 10:27:29 2021

@author: t.hoerbrand
"""


import numpy as np
import pandas as pd
import iapws
from plotnine import *
import pint as pint
import phreeqpy.iphreeqc.phreeqc_com as phreeqc_mod # install .com module for 64 bit from IPHREEQC 
from collections import namedtuple
from warnings import warn

# setup pint
u = pint.UnitRegistry()


# ----------------------------------------------------  
# functions
# ----------------------------------------------------

def make_chemical_model(chem_ions, chem_gas_mol, t, p, chem_units = "mg/L", pH = 7):
    '''  
    Function to create a PHREEQC model for a brine to calculate the boiling point and density.

    Parameters
    ----------
    chem_ions_mol : pd.DataFrame
        Contains the ionic composition of the fluid with the chemical units defined in chem_units.
    chem_gas_mol : pd.DataFrame
        Contains the gas amount of the individual gases in mol/kgw.
    t : float
        Temperature [°C]
    p : float
        Pressure [bara]
    chem_units : string, optional
        Units in which the chemical composition is supplied. 
        Options are "mg/L" (equal to ppm), "mol/L" (molarity) and "mol/kgw" (molality). 
        The default is "mg/L" which is equal to ppm.
    pH : float, optional
        DESCRIPTION. The default is 7.
        
    Returns
    -------
    chemical model string
    
    '''
    p = p/1.01325 # convert to atm for PHREEQC
    
    chemical_model = """
        TITLE 
        SOLUTION 1
        """ + \
        "\n units " + chem_units + \
        "\n temp " + str(t) + \
        "\n pH " + str(pH) + \
        "\n pressure " + str(p) + \
        "\n Na " + str(chem_ions["Na"][0]) + \
        "\n K " + str(chem_ions["Na"][0]) + \
        "\n Ca " + str(chem_ions["Na"][0]) + \
        "\n S(6) " + str(chem_ions["Na"][0]) + \
        "\n HCO3 " + str(chem_ions["Na"][0]) + \
        "\n Cl " + str(chem_ions["Cl"][0]) + \
        "\n EQUILIBRIUM_PHASES"   + \
        "\n Mtg(g) 10 " + str(chem_gas_mol["CH4"][0]) + \
        "\n Ntg(g) 10 " + str(chem_gas_mol["N2"][0]) + \
        "\n Hdg(g) 10 " + str(chem_gas_mol["H2"][0]) + \
        "\n CO2(g) 10 " + str(chem_gas_mol["CO2"][0]) + """
        SELECTED_OUTPUT
        -reset false
        USER_PUNCH
        -headings CO2 N2 CH4 H2O H2 rho t
        10 PUNCH 10^SI("CO2(g)") * 1.01325 /PR_PHI("CO2(g)")
        20 PUNCH 10^SI("Ntg(g)") * 1.01325 /PR_PHI("Ntg(g)")
        30 PUNCH 10^SI("Mtg(g)") * 1.01325 /PR_PHI("Mtg(g)")
        40 PUNCH 10^SI("H2O(g)") * 1.01325 /PR_PHI("H2O(g)")
        50 PUNCH 10^SI("Hdg(g)") * 1.01325 /PR_PHI("Hdg(g)")
        60 PUNCH RHO*1000
        70 PUNCH TC
        END
        """
        
    return chemical_model
      
  

def run_chemcial_model(chemical_model, database_path):  
    """
    

    Parameters
    ----------
    chemical_model : string
        PHREEQC input string, generated with "make_chemical_model"
    database_path : string
        Path to thermodynamic database used with phreeqc, e.g. r"C:/phreeqc/database/phreeqc.dat"

    Returns
    -------
    output : tuple
        results of PHREEQC simulation (selected output)
        Contains header in first row [0], intialization in second row [1] and results in third row [2]

    """
    
    # Initialize IPhreeqc
    phreeqc = phreeqc_mod.IPhreeqc()
    
    # Load thermodynamic database provided by phreeqc
    phreeqc.load_database(database_path)
    
    # Run chemical model
    phreeqc.run_string(chemical_model)  

    # Get selected output    
    #components = phreeqc.get_component_list()
    output = phreeqc.get_selected_output_array()
    
    return output


def get_chemical_properties(output):
    """
    Format simulation results from "run_chemical_model".

    Parameters
    ----------
    output : tuple
        Results of function "run_chemical_model".

    Returns
    -------
    partial_pressure : float
        Sum of partial pressures [atm].
    rho_phreeqc : float
        Denisty in [kg/m**3].
    t_phreeqc : float
        Temperature in [°C]

    """
    
    # Indices / output format provided with the USER_PUNCH block in function "make_chemical_model"
    partial_pressure = sum(output[2][0:5]) # why the hell is it not 0:4 here? Whatever it works...
    partial_pressure = partial_pressure
    rho_phreeqc = output[2][5]
    t_phreeqc = output[2][6]
    
    return partial_pressure, rho_phreeqc, t_phreeqc

def iterate_chemical_model(chem_ions, chem_gas_mol, t, p_convergence, database_path, tolerance = 0.01, chem_units = "mg/L", pH = 7):
    """  
    Function to iterate over chemical model

    Parameters
    ----------
    chem_ions_mol : pd.DataFrame
        Contains the ionic composition of the fluid with the chemical units defined in chem_units.
    chem_gas_mol : pd.DataFrame
        Contains the gas amount of the individual gases in mol/kgw.
    t : float
        Temperature [°C]
    p_convergence : float
        Pressure [bara]. It is converted to atm for PHREEQC in make_chemical_model
    chem_units : string, optional
        Units in which the chemical composition is supplied. 
        Options are "mg/L" (equal to ppm), "mol/L" (molarity) and "mol/kgw" (molality). 
        The default is "mg/L" which is equal to ppm.
    pH : float, optional
        DESCRIPTION. The default is 7.
        
    Returns
    -------
    chemical model string
    
    """
    
    upper_range = p_convergence + tolerance
    lower_range = p_convergence - tolerance
    stop_crit = 1 # stop criterion in case no convergence is reached
    pp = 0
    
    # brute force solver
    # adjusts temperature until pp and p_convergence are equal
    # stop criterion is, that partial pressure is within tolerance of convergence pressure
    while not(pp < upper_range and pp > lower_range):
        # run chemical model
        chemical_model = make_chemical_model(chem_ions, chem_gas_mol, t, p_convergence, chem_units = chem_units, pH = 7) 
        output = run_chemcial_model(chemical_model, database_path)
        pp, rho, t_phreeqc = get_chemical_properties(output)
        
        # adjust temperature for convergence
        if pp > upper_range:
            factor = p_convergence/pp + (1-p_convergence/pp)/1.03 # create relative factor for temperature adjustment
            t = t*factor
        elif pp < lower_range:
            factor = (p_convergence/pp + (1-p_convergence/pp)/1.03) # create relative factor for temperature adjustment
            t = t*factor
        #print(pp, p_convergence, t, factor, pp > upper_range, pp < lower_range)
        # stop criterion in case no convergence is reached
        stop_crit = stop_crit +1
        if stop_crit > 1000:
            warn("No convergence achieved.")
            return
        
    return pp, rho, t_phreeqc
        
# units are still a mess, need to be pinted
def bpdc_gas(depth, p0, chem_ions, chem_gas_mol, database_path, tolerance = 0.01, chem_units = "mol/kgw", units='common'):
    '''
    Function to calculate the hydrostatic boiling point depth curve, similar to:
    Haas Jr., J.L., 1971. The effect of salinity on the maximum thermal gradient 
    of a hydrothermal system at hydrostatic pressure. Econ. Geol. 66, 940–946.

    Parameters
    ----------
    depth : array-like
        An array of depths at which to evaluate the function. 
        Depth intervals should not be bigger than 1 m!
    p : float, optional
        surface pressure of well. The default is 1.01325 bar.
    method : string
        Method to calculate brine density. The default is "iapws".
    units : string or tuple
        String units for depth, pressure, temperature, and density
        respectively. The following are allowed:
        - 'common' to use metres, bar, deg Celsius, and kg/m**3.
        - 'SI' to use metres, Pascal, Kelvin, and kg/m**3.
        - 'Imperial' to use feet, psia, deg Fahrenheit, and g/cm**3.
        - tuple like ('m', 'atm', 'degC', 'g/cm**3') for custom units.

    Returns
    -------
    namedtuple
    '''
    u = UnitRegistry()

    g = 9.81 * u.m / u.s**2

    # Assign units to everything.
    if units == 'SI':
        u_d, u_p, u_t, u_rho = u.m, u.pa, u.K, u.kg/u.m**3
    if units == 'common':
        u_d, u_p, u_t, u_rho = u.m, u.bar, u('degC'), u.kg/u.m**3
    elif units == 'Imperial':
        u_d, u_p, u_t, u_rho = u.ft, u.psi, u('degF'), u.g/u.cm**3
    else:
        u_d, u_p, u_t, u_rho = list(map(u, units))
    
    # Override units with pint Quantity, if possible.
    # And compute the diff array.
    if isinstance(depth, pint.Quantity):
        u_d = depth.units
    else:
        depth = np.asanyarray(depth) * u_d
    depth_diff = np.diff(depth)
    if max(depth_diff.magnitude) > 1:
        warn("Step size of depth interval > 1. This leads to an error propagation in calculated temperatures. Consider reducing the interval size.")
                
    # Override units with pint Quantity, if possible.
    # And deal with not getting a p0 pressure.
    if isinstance(p0, pint.Quantity):
        u_p = p0.units
    elif p0 is None:
        p0 = 101325 * u_p
    else:
        p0 = p0 * u_p  
    p0 = p0.to('bar')  

    # Compute using PHREEQC option.
    pressure = np.atleast_1d(p0)
    # select starting temperature for iterate_chemical_models (reduces number of iterations)
    t_model = 100
    #print(chem_ions, chem_gas_mol, t_model, p0.to('atm').m, database_path, tolerance, chem_units)
    # iterate chemical model until partial pressure of model is equal to pressure
    pp, rho, tsat = iterate_chemical_model(chem_ions = chem_ions, 
                                           chem_gas_mol = chem_gas_mol, 
                                           t = t_model, # temperature is calculated in iterate_chemical_model; this is a starting value
                                           p_convergence = p0.to('atm').m, 
                                           database_path = database_path, 
                                           tolerance = tolerance, 
                                           chem_units = chem_units) 
    #tsat = tsat * u.K
    #print(tsat)
    tsat = np.atleast_1d(tsat)
    rho = rho * u.kg / u.m**3
    density = np.atleast_1d(rho)
    
    for i in depth_diff:
        # Calculate new pressure for this step.
        new_p = pressure[-1] + rho * g * i
        pressure = np.append(pressure, new_p)  # Has units.
        
        # select starting temperature for iterate_chemical_models (reduces number of iterations)
        if 'tsat' in locals():
            t_model = tsat[-1]
        else:
            t_model = 100
            
        # iterate chemical model until partial pressure of model is equal to pressure
        pp, rho, t = iterate_chemical_model(chem_ions = chem_ions, 
                                               chem_gas_mol = chem_gas_mol, 
                                               t = t_model, # temperature is calculated in iterate_chemical_model; this is a starting value
                                               p_convergence = pressure[-1].m, 
                                               database_path = database_path, 
                                               tolerance = tolerance, 
                                               chem_units = chem_units) 
        
        # Calculate new temperature for this step.
        #t = t * u.K
        tsat = np.append(tsat, t)
        print(tsat)
        # Calculate new density for next step.
        rho = rho * u.kg / u.m**3
        density = np.append(density, rho)
    
    # Finalize units.
    pressure = pressure.to(u_p)
    #tsat = tsat.to(u_t)
    density = density.to(u_rho)
         
        
    # Return a namedtuple to the user to retain units.
    BPD = namedtuple('BPD', ['depth', 'pressure', 'tsat', 'rho'])
    return BPD(depth, pressure, tsat, density)



      
  
        



# ----------------------------------------------------  
# correct for gases
# ---------------------------------------------------- 

chem_ions_mol = pd.DataFrame([[0.5, 0.5 , 0.5, 0.5, 0.5, 0.5]], columns = ['Na', 'K', 'Ca', 'Cl', 'SO4', 'HCO3'])
chem_gas_mol = pd.DataFrame([[0.01, 0.01, 0.001, 0.001]], columns = ['CH4', 'CO2', 'N2', 'H2'])






gas_curve = bpdc_gas(depth = df["depth_m"][0:40],
                p0=df["pressure_bara"].min(),
                chem_ions = chem_ions_mol,
                chem_gas_mol = chem_gas_mol,
                database_path = r"C:/phreeqc/database/pitzer.dat",
                tolerance = 0.1,
                chem_units = "mol/kgw")

gas_df = pd.DataFrame(gas_curve._asdict())

 




# this is a mess from here, I'll clean up on another day

pp, rho, t = iterate_chemical_model(chem_ions_mol, chem_gas_mol, 
                                    t = 150, 
                                    p_convergence = 19.086, 
                                    database_path = r"C:/phreeqc/database/pitzer.dat", 
                                    tolerance = 0.01, 
                                    chem_units = "mol/kgw") 

chemical_model = make_chemical_model(chem_ions_mol, chem_gas_mol, 208.23, 19.086, chem_units = "mol/kgw", pH = 7) 
output = run_chemcial_model(chemical_model, r"C:/phreeqc/database/pitzer.dat")
pp, rho, tsat = get_chemical_properties(output)

iapws.iapws97._TSat_P(1.9086)-273.15
      
pd.DataFrame([output[2]], columns = output[0])   

# only salinity correction
# sal_curve = bpdc_salinity(depth = df["depth_m"],
#                 p0=df["pressure_bara"].min(),
#                 chem_ions = chem_ions_mol,
#                 chem_gas_mol = chem_gas_mol,
#                 database_path = r"C:/phreeqc/database/pitzer.dat",
#                 tolerance = 0.1,
#                 chem_units = "mol/kgw")

# sal_df = pd.DataFrame(sal_curve._asdict())


# ----------------------------------------------------  
# compare results
# ---------------------------------------------------- 

# ggplot(df, aes(y = "depth_m")) +\
#     geom_line(aes(x="tsat_degC", colour = "'well pressure + iapws'")) + \
#     geom_line(aes(x = "temp_degC", colour = "'well data only'")) + \
#     scale_y_reverse() + \
#     labs(x = "T [°C]", y = "depth [m]", colour = "legend")

# ggplot(df, aes(y = "depth_m")) +\
#     geom_line(aes(x="tsat_degC", colour = "'well pressure + iapws'")) + \
#     geom_line(aes(x = "temp_degC", colour = "'well data only'")) + \
#     geom_line(aes(x = "tsat", y = "depth", colour = "'hydrostatic + iapws'"), data = hyd_df)  + \
#     scale_y_reverse() + \
#     labs(x = "T [°C]", y = "depth [m]", colour = "legend")
    
# ggplot(df, aes(y = "depth_m")) +\
#     geom_line(aes(x="tsat_degC", colour = "'well pressure + iapws'")) + \
#     geom_line(aes(x = "temp_degC", colour = "'well data only'")) + \
#     geom_line(aes(x = "tsat", y = "depth", colour = "'hydrostatic + iapws'"), data = hyd_df) + \
#     geom_line(aes(x = "tsat", y = "depth", colour = "'salinity corrected'"), data = sal_df) + \
#     scale_y_reverse() + \
#     labs(x = "T [°C]", y = "depth [m]", colour = "legend")
    
ggplot(df, aes(y = "depth_m")) +\
    geom_line(aes(x="tsat_degC", colour = "'well pressure + iapws'")) + \
    geom_line(aes(x = "temp_degC", colour = "'well data only'")) + \
    geom_line(aes(x = "tsat", y = "depth", colour = "'hydrostatic + iapws'"), data = hyd_df) + \
    geom_line(aes(x = "tsat", y = "depth", colour = "'gas and salinity corrected'"), data = gas_df) + \
    scale_y_reverse() + \
    labs(x = "T [°C]", y = "depth [m]", colour = "legend")
    
# ggplot(df, aes(y = "pressure_bara")) +\
#     geom_line(aes(x="tsat_degC", colour = "'well pressure + iapws'")) + \
#     geom_line(aes(x = "temp_degC", colour = "'well data only'")) + \
#     geom_line(aes(x = "tsat", y = "pressure", colour = "'hydrostatic + iapws'"), data = hyd_df) + \
#     geom_line(aes(x = "tsat", y = "pressure", colour = "'gas and salinity corrected'"), data = gas_df) + \
#     scale_y_reverse() + \
#     labs(x = "T [°C]", y = "pressure [bara]", colour = "legend")
    
# ggplot(df, aes(y = "depth_m")) + \
#     geom_line(aes(x="pressure_bara", colour = "'well pressure + iapws'")) + \
#     geom_line(aes(x = "pressure_bara", colour = "'well data only'")) + \
#     geom_line(aes(x = "pressure", y = "depth", colour = "'hydrostatic + iapws'"), data = hyd_df) + \
#     scale_y_reverse() + \
#     labs(x = "T [°C]", y = "depth [m]", colour = "legend")

# ggplot(df, aes(y = "depth_m")) +\
#     geom_line(aes(y = "rho", x = "tsat", colour = "'hydrostatic + iapws'"), data = hyd_df) + \
#     geom_line(aes(y = "rho", x = "tsat", colour = "'gas and salinity corrected'"), data = gas_df) + \
#     labs(y = "rho", x = "T [°C]", colour = "legend")


# general observations / assumptions etc
# * if the pressure is supplied, the hydrostatic pressure does not have to be calculated
# * Arnorsson: Rising hot waters begin to boil when steam pressure plus total gas pressure become equal to the hydrostatic pressure.



        
    
