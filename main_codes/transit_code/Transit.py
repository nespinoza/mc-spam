import sys
import ajplanet
import numpy as np
import occultquad


def rot(a,b,theta):
    "rotates a,b by theta"
    theta_rad = (np.pi / 180.0) * theta
    return a*np.cos(theta_rad) - b*np.sin(theta_rad), a*np.sin(theta_rad) + b*np.cos(theta_rad)

def uncorTransit(t, P, e, omega, t0,ld_1,ld_2,delta,TT,tau):
    """
    Calculates a transit using less-correlated parameters, from which the other pars are
    determined
    ld_1 =  u_1*cos(40) - u_2*sin(40)
    ld_2 =  u_1*sin(40) + u_2*cos(40)
    
    delta, T, tau defined as in Carter et al 2008
    """
    u1, u2 = rot( ld_1, ld_2, -45 )
    
    p = np.sqrt( delta )
    r = p # to follow notation of Carter et al
    sq = 1-r*TT/tau
    #print r, TT, tau
    r_a = (2 * np.pi / P) * ( (1+e*np.sin(omega)) / np.sqrt(1-e**2) ) * np.sqrt( (TT*tau)/(4*r) )
    if (sq >= 0):
        arg = r_a * ( (1+e*np.sin(omega)) / (1-e**2) ) * np.sqrt(1-r*TT/tau)
        if (arg < 1):
            i = np.arccos( arg )
        else:
            i = np.pi/2
    else:
        i = 0

    return Transit(t, P, i, r_a, p, t0, u1, u2, e, omega)

def uncorTransit2(t, P, e, omega, t0,u1,u2,delta,TT,tau):
    """
    Calculates a transit using less-correlated parameters, from which the other pars are
    determined
    ld_1 =  u_1*cos(40) - u_2*sin(40)
    ld_2 =  u_1*sin(40) + u_2*cos(40)
    
    delta, T, tau defined as in Carter et al 2008
    """    
    p = np.sqrt( delta )
    r = p # to follow notation of Carter et al
    sq = 1-r*TT/tau
    #print r, TT, tau
    r_a = (2 * np.pi / P) * ( (1+e*np.sin(omega)) / np.sqrt(1-e**2) ) * np.sqrt( (TT*tau)/(4*r) )
    if (sq >= 0):
        arg = r_a * ( (1+e*np.sin(omega)) / (1-e**2) ) * np.sqrt(1-r*TT/tau)
        if (arg < 1):
            i = np.arccos( arg )
        else:
            i = np.pi/2
    else:
        i = 0

    return Transit(t, P, i, r_a, p, t0, u1, u2, e, omega)


def Transit_stdpar(t, P, e, omega, ld_1,ld_2,delta):
    """
    Calculates a transit using  standard parameters.
    Meant to be used with fixed inclination and t0 for out program

    ld_1 =  u_1*cos(40) - u_2*sin(40)
    ld_2 =  u_1*sin(40) + u_2*cos(40)
    
    """
    u1, u2 = rot( ld_1, ld_2,-45)
    
    p = np.sqrt( delta )

    return Transit(t, P, i_ONTRANSIT, r_a_ONTRANSIT,\
                       p, t0_ONTRANSIT, u1, u2, e, omega)

def Transit(t, P, i, r_a, p, t0, u1, u2, e, omega):
    """ 
    Calculates transit light curve, normalized
    non-obvious pars:
    r_a=R_*/a
    p = R_p/R_*
    """

    if (e > 0):
        f = np.pi/2 - omega
        E = 2.0 * np.arctan2(np.sqrt( (1.0-e)/(1.0+e)) * np.sin(f/2.0), np.cos(f/2.0) )
        n = 2 * np.pi / P
        tperi = t0 - (E - e*np.sin(E))/n
    else:
        tperi = t0 - P/4.0

    Omega = np.pi

    (X,Y,Z) = ajplanet.pl_Orbit_array(t,tperi,P,1,e,omega,i,Omega)
    rpsky = (1/r_a) * np.sqrt(X**2 + Y**2)
    flux,flux0 = ajplanet.occultquad(rpsky,u1,u2,p)

    return flux, flux0

def uncorTransit_nl(t, P, e, omega, t0,c1,c2,c3,c4,delta,TT,tau):
    """
    Calculates a transit using less-correlated parameters, from which the other pars are
    determined
    
    delta, T, tau defined as in Carter et al 2008
    """
    p = np.sqrt( delta )
    r = p # to follow notation of Carter et al
    sq = 1-r*TT/tau
    r_a = (2 * np.pi / P) * ( (1+e*np.sin(omega)) / np.sqrt(1-e**2) ) * np.sqrt( (TT*tau)/(4*r) )
    if (sq >= 0):
        arg = r_a * ( (1+e*np.sin(omega)) / (1-e**2) ) * np.sqrt(1-r*TT/tau)
        if (arg < 1):
            i = np.arccos( arg )
        else:
            i = np.pi/2
    else:
        i = 0

    return Transit_nl(t, P, i, r_a, p, t0, c1,c2,c3,c4, e, omega)


def Transit_nl(t, P, i, r_a, p, t0, c1,c2,c3,c4, e, omega):
    """ 
    Calculates transit light curve, normalized
    non-obvious pars:
    r_a=R_*/a
    p = R_p/R_*
    """

    if (e > 0):
        f = np.pi/2.0 - omega
        E = 2.0 * np.arctan2(np.sqrt( (1.0-e)/(1.0+e)) * np.sin(f/2.0), np.cos(f/2.0) )
        n = 2 * np.pi / P
        tperi = t0 - (E - e*np.sin(E))/n
    else:
        tperi = t0 - P/4.0

    Omega = np.pi
    
    (X,Y,Z) = ajplanet.pl_Orbit_array(t,tperi,P,1, e,omega,i,Omega)
    rpsky = (1/r_a) * np.sqrt(X**2 + Y**2)
    
    flux,incog = ajplanet.occultnl(p,c1,c2,c3,c4,rpsky)

    return flux, incog
