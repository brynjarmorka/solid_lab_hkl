import numpy as np

"""
This is the script for calculating possible hkl values based on the
lenght of two translational vectors and the angle between them. The
script is adjusted for FCC lattices with the condition that all hkl
values are even or odd.

The first lines of code are for selecting the input values. The format:
# [assumed orientation, i.e. which diffraction pattern] 
# the measured values, and some metadata. [R1, R2, theta, assumed 
# orientation, [zone axis from assumed orientation and simulation]]
# >>>> resulting R1- and R2-vectors  >>> zone axis direction -> family

"""

# [100]
R = np.array([4.94, 6.98, np.pi*45.5/180, 100, [0,0,-4]], dtype=object)     
# >>>> [2,0,0], [2,2,0]       >>> P = [0,0,-4] -> {100}

# [110]
# R = np.array([4.23, 4.91, np.pi*54.2/180, 110, [2,2,0]], dtype=object)      
# >>>> [0,0,2], [-1,1,1]      >>>> P = [2,2,0]

# [110], alternative measurement
# R = np.array([6.92, 4.23, np.pi*35.5/180, 110, [2,2,0]], dtype=object)        
# >>>> [-2,2,0], [-1,1,1]     >>>> P = [-2,-2,0] -> {110}


# Boolean values to skip certain if statements. Skip when TRUE
skip_zone_axis_filter = False

# Constants
R_ratio = (R[0]/R[1])**2 # ratio between R1 og R2. Left side of equation (3).
# Adjusting the max and min hkl give exponentially more possibilities
# which is possible through the equations, but with the zone axis 
# filter the result stays the same for higher and lower values than +-3
max_hkl = 4 # max value for hkl numbers
min_hkl = -3 # min value for hkl numbers
# The threshold values for ratio and the angle
ratio_threshold = 0.05 # terskel for ratio, likning 10 i TEM-labinfo
angle_threshold = 0.1 # terskel for cos-ratio, likning 13 i TEM-labinfo


#
# Helper functions start
#

def hkl_sqare(list): 
    """ Calculates the sum of the square of each hkl value. """
    return (list[0]**2 + list[1]**2 + list[2]**2)

def cos_check(l1,l2):
    """ Calculates the cosinus check, i.e. equation (4) in the report. """
    return ( l1[0]*l2[0] + l1[1]*l2[1] + l1[2]*l2[2] )/(np.sqrt(hkl_sqare(l1))*np.sqrt(hkl_sqare(l2)))

def is_odd(num):
    """ Checking if the last bit of a integer is 1, i.e. if odd. """
    return num & 0x1

def is_even(num):
    """ Checking if the last bit of a integer is 0, i.e. if even. """
    return is_odd(num+1)

#
# Helper functions end
#


# lage listen med alle mulige hkl-kombinasjoner fra 0 til max_hkl.
# sjekker at hkl ikke er [0,0,0].
# gjør FCC-sjekken, altså kun odde- eller partall i én hkl
def make_hkl():
    """ 
    This function makes a list of all the possible hkl combinations for
    FCC lattices, with lowest value min_hkl and highest max_hkl-1. This
    function does the FCC check.
    Returns list(hkl_combinations) and list(hkl_squared_sums)
    
    
    """
    hkl = [] # first return value
    hkl_square_sum = [] # second return value
    for h in range (min_hkl, max_hkl):
            for k in range(min_hkl, max_hkl):
                    for l in range(min_hkl, max_hkl):
                        midlertidig = [h, k, l] # temporary storage of the hkl
                        # Will not append or check if [h,k,l] = [0,0,0]
                        if midlertidig != [0,0,0]:
                            # This is the FCC check, i.e. all even or all odd.
                            if (is_odd(h) and is_odd(k) and is_odd(l)) or (is_even(h) and is_even(k) and is_even(l)):
                                hkl.append(midlertidig)
                                hkl_square_sum.append(hkl_sqare(midlertidig))
    return hkl, hkl_square_sum


# sjekker først om ratio R_1 / R_2 er innenfor ratio_threshold
# sjekker så om vinkelen er innenfor angle_threshold
# skriver ut mulige hkl-kombinasjoner

def possible_hkl_check(hkl, hkl_square_sum, ratio_threshold, angle_threshold):
    """
    This is the function that does the ratio check, angle check and the
    zone axis check. The function returns the whole list of hkl values
    that match with the ratio and the angle, and it returns a much 
    shorter list with the hkl values that match with the zone axis.
    The function now also print the zone axis matches.

    Returns list(many_hkl) and list(few_hkl)
    
    """
    possbile_hkl = []
    hkl_zone_axis_match = []
    for i in range (len(hkl)):
        for j in range (len(hkl)):
            ratio = hkl_square_sum[i]/hkl_square_sum[j] # Calculating the ij-th ratio
            # Comparing the ij-th ratio with the ratio in the image
            if abs(R_ratio-ratio) < ratio_threshold:
                # Comparing the cos equation value with the cos of the angle in the image
                if abs(np.cos(R[2])-(cos_check((hkl[i]), (hkl[j])))) < angle_threshold:
                    # Appending hkl combinations that pass the ratio and angle check
                    possbile_hkl.append([hkl[i], hkl[j]])
                    # Filtering out the hkl pairs which have the same zone axis as the simulation
                    if (str(np.cross(hkl[i],hkl[j])) == str(np.array(R[4])) or skip_zone_axis_filter):
                        print(f"[{int(R[3])}]: R1 = {hkl[i]} & R2 = {hkl[j]}, crossproduct: {np.cross(hkl[i],hkl[j])}")
                        # print(f"& & {hkl[i]} & & {hkl[j]} ") # printing out the filteret ones for Latex
                        hkl_zone_axis_match.append([hkl[i], hkl[j]])                    
    return possbile_hkl, hkl_zone_axis_match

# fist makes the possible hkl list and sum list
hkl_list, hkl_square_sum = make_hkl()

# then checks which are possibilities
possibilities, hkl_zone_axis_match = possible_hkl_check(hkl_list, hkl_square_sum, ratio_threshold, angle_threshold)

# Gives the number of possible matches with angle and ratio check
print(f"Number of all possbile matches without the zone axis filter: {len(possibilities)}")