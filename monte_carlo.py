from math import log as ln
from math import sqrt as sqrt
from math import cos as cos
from math import sin as sin
from math import acos as arccos
from math import pi as pi
import random
import time
INF = 1000000000000000000000

def linspace(start,stop,n):
	res = [start]
	h = (stop-start)/n
	for i in range(0,n):
		res.append(start+((i+1)*h))
	
	return res

#Geometry Functions
#-----------------------------------------------------------------------------------------

#Define a function to check if a point is in a region, but not on the boundary
#Only works for rectangular geometries
def in_region(position,region):
	if(((position[0] > region[0][0]) and (position[0] < region[1][0])) and ((position[1] > region[0][1]) and (position[1] < region[2][1]))):
		return True
	else:
		return False

#Define a function to check if a point is in a region, including the boundary
#Only works for rectangular geometries
def on_region(position,region):
	if(((position[0] >= region[0][0]) and (position[0] <= region[1][0])) and ((position[1] >= region[0][1]) and (position[1] <= region[2][1]))):
		return True
	else:
		return False

#Define a function that finds the intersect of a line and a vertical line segment defined by two corner points
#Note: corner1 must be the lower corner
#Note: pos1 must be the "last" position
def vertical(pos1,pos2,corner1,corner2):
	x1 = pos1[0]
	x2 = pos2[0]
	y1 = pos1[1]
	y2 = pos2[1]
	x = corner1[0]
	
	m = (y2-y1)/(x2-x1)
	
	y = m*(x - x1) + y1
	
	x = round(x,8)
	y = round(y,8)
	
	if((y >= corner1[1]) and (y < corner2[1])):
		return [x,y]
	else:
		return None

#Define a function that finds the intersect of a line and a horizontal line segment defined by two corner points
#Note: corner1 must be the left-most corner
#Note: pos1 must be the "last" position
def horizontal(pos1,pos2,corner1,corner2):
	x1 = pos1[0]
	x2 = pos2[0]
	y1 = pos1[1]
	y2 = pos2[1]
	y = corner1[1]
	
	m = (y2-y1)/(x2-x1)
	
	x = ((y-y1)/m) + x1
	
	x = round(x,8)
	y = round(y,8)
	
	if((x >= corner1[0]) and (x < corner2[0])):
		return [x,y]
	else:
		return None

#-----------------------------------------------------------------------------------------



#Materials
#-----------------------------------------------------------------------------------------

#Import material cross sections
materials = ['pb206','pb207','pb208','h1','o16','c12','b10','b11','cd106','cd110','cd111','cd112','cd113','cd114','cd116','fe56','si28','ca40']
LIBRARY = {'el':{},'g':{}}

#Make logarithmic bin for cross section energies (for computing speed)
binning = []
for i in [1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0]:
	for j in [-9.0,-8.0,-7.0,-6.0,-5.0,-4.0,-3.0,-2.0,-1.0,0.0]:
		binning.append(i * (10**j))
binning.append(10.0)
binning.sort()

for i in range(0,len(materials)):
	for j in ['el','g']:
		filename = '/Users/puptheblue/Dropbox/Berkeley_Courses/Spring_2015/NE_155/Project/'+materials[i]+'_'+j+'.txt'
		file = open(filename,'r')
		lines = file.readlines()
		lines = lines[11:len(lines)-3]
		
		sub_lib = {}
		num = {}
		xss = {}
		for en in binning:
			num[en] = 0
			xss[en] = 0.0
		for line in lines:
			[energy,xs] = line.split()
			energy = float(energy)
			if(energy > 15):
				continue
			xs = float(xs)
			dif = float('inf')
			best_en = None
			for en in binning:
				if(abs(1 - en/energy) < dif):
					dif = abs(1 - en/energy)
					best_en = en
			num[best_en] += 1
			xss[best_en] += xs
		for en in num.keys():
			if(num[en] != 0):
				sub_lib[en] = (xss[en]/num[en])*(10**-24)
		
		LIBRARY[j][materials[i]] = sub_lib

#Define material densities
#             lead,  water, B4C,cadmium,portland cement,iron
densities = [11.34,18.0153,2.51,8.56   ,3.15           ,7.874]
#molar_masses = [205.974465278,206.975896887,207.976652071,1.00782503207,15.99491461956,12,10.012936992,11.009305406,105.90645941,109.90300207,110.904178107,111.902757809,112.904401662,113.90335854,115.904755809,27.97692653246,39.962590983]
molar_masses = [207.2,1,55.255,112.411,(0.253*60.0843 + 0.747*56.077),55.85]
abundances = [0.241,0.221,0.524,1,1,1,0.199,0.801,0.0125,0.1249,0.1280,0.2413,0.1222,0.2873,0.0749,1,1,1]
Na = 6.022 * (10**23)

def get_Sig(compound,x,energy):
	pb206_energy = None
	dif = float('inf')
	for en in LIBRARY[x]['pb206'].keys():
		if(abs(1 - en/energy) < dif):
			dif = abs(1 - en/energy)
			pb206_energy = en
	
	pb207_energy = None
	dif = float('inf')
	for en in LIBRARY[x]['pb207'].keys():
		if(abs(1 - en/energy) < dif):
			dif = abs(1 - en/energy)
			pb207_energy = en
			
	pb208_energy = None
	dif = float('inf')
	for en in LIBRARY[x]['pb208'].keys():
		if(abs(1 - en/energy) < dif):
			dif = abs(1 - en/energy)
			pb208_energy = en
	
	h1_energy = None
	dif = float('inf')
	for en in LIBRARY[x]['h1'].keys():
		if(abs(1 - en/energy) < dif):
			dif = abs(1 - en/energy)
			h1_energy = en
			
	o16_energy = None
	dif = float('inf')
	for en in LIBRARY[x]['o16'].keys():
		if(abs(1 - en/energy) < dif):
			dif = abs(1 - en/energy)
			o16_energy = en
	
	c12_energy = None
	dif = float('inf')
	for en in LIBRARY[x]['c12'].keys():
		if(abs(1 - en/energy) < dif):
			dif = abs(1 - en/energy)
			c12_energy = en
	
	b10_energy = None
	dif = float('inf')
	for en in LIBRARY[x]['b10'].keys():
		if(abs(1 - en/energy) < dif):
			dif = abs(1 - en/energy)
			b10_energy = en
			
	b11_energy = None
	dif = float('inf')
	for en in LIBRARY[x]['b11'].keys():
		if(abs(1 - en/energy) < dif):
			dif = abs(1 - en/energy)
			b11_energy = en
			
	cd106_energy = None
	dif = float('inf')
	for en in LIBRARY[x]['cd106'].keys():
		if(abs(1 - en/energy) < dif):
			dif = abs(1 - en/energy)
			cd106_energy = en
			
	cd110_energy = None
	dif = float('inf')
	for en in LIBRARY[x]['cd110'].keys():
		if(abs(1 - en/energy) < dif):
			dif = abs(1 - en/energy)
			cd110_energy = en
			
	cd111_energy = None
	dif = float('inf')
	for en in LIBRARY[x]['cd111'].keys():
		if(abs(1 - en/energy) < dif):
			dif = abs(1 - en/energy)
			cd111_energy = en
			
	cd112_energy = None
	dif = float('inf')
	for en in LIBRARY[x]['cd112'].keys():
		if(abs(1 - en/energy) < dif):
			dif = abs(1 - en/energy)
			cd112_energy = en
			
	cd113_energy = None
	dif = float('inf')
	for en in LIBRARY[x]['cd113'].keys():
		if(abs(1 - en/energy) < dif):
			dif = abs(1 - en/energy)
			cd113_energy = en
			
	cd114_energy = None
	dif = float('inf')
	for en in LIBRARY[x]['cd114'].keys():
		if(abs(1 - en/energy) < dif):
			dif = abs(1 - en/energy)
			cd114_energy = en
			
	cd116_energy = None
	dif = float('inf')
	for en in LIBRARY[x]['cd116'].keys():
		if(abs(1 - en/energy) < dif):
			dif = abs(1 - en/energy)
			cd116_energy = en
			
	fe56_energy = None
	dif = float('inf')
	for en in LIBRARY[x]['fe56'].keys():
		if(abs(1 - en/energy) < dif):
			dif = abs(1 - en/energy)
			fe56_energy = en
	
	si28_energy = None
	dif = float('inf')
	for en in LIBRARY[x]['si28'].keys():
		if(abs(1 - en/energy) < dif):
			dif = abs(1 - en/energy)
			si28_energy = en
			
	ca40_energy = None
	dif = float('inf')
	for en in LIBRARY[x]['ca40'].keys():
		if(abs(1 - en/energy) < dif):
			dif = abs(1 - en/energy)
			ca40_energy = en
	
	Sigs = {}
	Sigs['lead'] = 0.241*((11.34*Na)/207.2)*(LIBRARY[x]['pb206'][pb206_energy]) + 0.221*((11.34*Na)/207.2)*(LIBRARY[x]['pb207'][pb207_energy]) + 0.524*((11.34*Na)/207.2)*(LIBRARY[x]['pb208'][pb208_energy])
	Sigs['water'] = 2*((1*Na)/18.0153)*LIBRARY[x]['h1'][h1_energy] + ((1*Na)/18.0153)*LIBRARY[x]['o16'][o16_energy]
	Sigs['B4C'] = 4*0.199*((2.51*Na)/55.255)*LIBRARY[x]['b10'][b10_energy] + 4*0.801*((2.51*Na)/55.255)*LIBRARY[x]['b11'][b11_energy] + ((2.51*Na)/55.255)*LIBRARY[x]['c12'][c12_energy]
	Sigs['cadmium'] = 0.0125*((8.56*Na)/112.411)*LIBRARY[x]['cd106'][cd106_energy] + 0.1249*((8.56*Na)/112.411)*LIBRARY[x]['cd110'][cd110_energy] + 0.1280*((8.56*Na)/112.411)*LIBRARY[x]['cd111'][cd111_energy] + 0.2413*((8.56*Na)/112.411)*LIBRARY[x]['cd112'][cd112_energy] + 0.1222*((8.56*Na)/112.411)*LIBRARY[x]['cd113'][cd113_energy] + 0.2873*((8.56*Na)/112.411)*LIBRARY[x]['cd114'][cd114_energy] + 0.0749*((8.56*Na)/112.411)*LIBRARY[x]['cd116'][cd116_energy]
	Sigs['portland cement'] = 0.253*((3.15*Na)/(0.253*60.0843 + 0.747*56.077))*(LIBRARY[x]['si28'][si28_energy] + 2*LIBRARY[x]['o16'][o16_energy]) + 0.747*((3.15*Na)/(0.253*60.0843 + 0.747*56.077))*(LIBRARY[x]['ca40'][ca40_energy] + LIBRARY[x]['o16'][o16_energy])
	Sigs['iron'] = ((7.874*Na)/55.85)*LIBRARY[x]['fe56'][fe56_energy]
	Sigs['iron poly'] = ((2.681*Na)/97.90)*LIBRARY[x]['fe56'][fe56_energy] + 3*((2.681*Na)/97.90)*LIBRARY[x]['c12'][c12_energy] + 6*((2.681*Na)/97.90)*LIBRARY[x]['h1'][h1_energy]
	Sigs['void'] = 10**-8
	Sigs['grave'] = 10000000000000000000000000000000000000000
	
	return Sigs[compound]

#-----------------------------------------------------------------------------------------



#Inverted CDFs
#-----------------------------------------------------------------------------------------

#CDF for neutron elastic scattering
def P_el(x,A):
	if(A == 1.0):
		A += 0.000001
	constants = ((32*A)/( (2*(A**3 + A)) + (A**4 - (2*(A**2)) +1)*(ln(A**2 -1) - (2*ln(A+1))) ))
	upper = ((x * sqrt(A**2 + x**2 -1) * (A**2 + 2*(x**2) -1)) - (((A**2 -1)**2) * ln(sqrt(A**2 + x**2 -1) + x)))/(16*A)
	x = 0
	lower = ((x * sqrt(A**2 + x**2 -1) * (A**2 + 2*(x**2) -1)) - (((A**2 -1)**2) * ln(sqrt(A**2 + x**2 -1) + x)))/(16*A)
	#upper = (x * sqrt((A**2) + (x**2) -1) * ((A**2) + (2 * (x**2)) -1) -((((A**2) -1)**2) * ln((sqrt((A**2)+(x**2)-1))+x)))/(2 * ((A**3)+A)+(((A**2)-1))**2 * ln(A-1)-(((A**2)-1)**2) * ln(A+1))
	#x = -1
	#lower = (x * sqrt((A**2) + (x**2) -1) * ((A**2) + (2 * (x**2)) -1) -((((A**2) -1)**2) * ln((sqrt((A**2)+(x**2)-1))+x)))/(2 * ((A**3)+A)+(((A**2)-1))**2 * ln(A-1)-(((A**2)-1)**2) * ln(A+1))
	
	#return upper-lower
	return constants*(upper-lower)

#Inverse CDF for neutron elastic scattering
def P_el_inv(r,n,A):
	x = linspace(-1.0,1.0,n)
	y = []
	for i in range(0,n+1):
		y.append(P_el(x[i],A))
	
	#Invert
	n = len(x)-1
	mu = 0
	for i in range(0,n):
		if((r >= y[i]) and (r < y[i+1])):
			mu = x[i]
	
	theta = None
	r = random.random()
	if((r >= 0) and (r < 0.5)):
		theta = -1*arccos(mu)
	else:
		theta = arccos(mu)
	
	time2 = time.time()
	return theta
	
#Inverse CDF for path length
def P_path_inv(r,Sig):
	l = -1*(ln(1-r)/Sig)
	
	return l

#Inverse CDF for interaction type
def P_int_inv(r,Sigs):
	x = 0
	probs = [0]
	for i in range(0,len(Sigs)):
		Sig = Sigs[i]
		p = Sig/(sum(Sigs))
		probs.append(p + sum(probs))
	
	for i in range(0,len(probs)-1):
		if((r >= probs[i]) and (r < probs[i+1])):
			x = i
	
	return x

#Inverse CDF for atom scattered off of
def P_A_inv(r,x,compound,energy):
	pb206_energy = None
	dif = float('inf')
	for en in LIBRARY[x]['pb206'].keys():
		if(abs(1 - en/energy) < dif):
			dif = abs(1 - en/energy)
			pb206_energy = en
	
	pb207_energy = None
	dif = float('inf')
	for en in LIBRARY[x]['pb207'].keys():
		if(abs(1 - en/energy) < dif):
			dif = abs(1 - en/energy)
			pb207_energy = en
			
	pb208_energy = None
	dif = float('inf')
	for en in LIBRARY[x]['pb208'].keys():
		if(abs(1 - en/energy) < dif):
			dif = abs(1 - en/energy)
			pb208_energy = en
	
	h1_energy = None
	dif = float('inf')
	for en in LIBRARY[x]['h1'].keys():
		if(abs(1 - en/energy) < dif):
			dif = abs(1 - en/energy)
			h1_energy = en
			
	o16_energy = None
	dif = float('inf')
	for en in LIBRARY[x]['o16'].keys():
		if(abs(1 - en/energy) < dif):
			dif = abs(1 - en/energy)
			o16_energy = en
	
	c12_energy = None
	dif = float('inf')
	for en in LIBRARY[x]['c12'].keys():
		if(abs(1 - en/energy) < dif):
			dif = abs(1 - en/energy)
			c12_energy = en
	
	b10_energy = None
	dif = float('inf')
	for en in LIBRARY[x]['b10'].keys():
		if(abs(1 - en/energy) < dif):
			dif = abs(1 - en/energy)
			b10_energy = en
			
	b11_energy = None
	dif = float('inf')
	for en in LIBRARY[x]['b11'].keys():
		if(abs(1 - en/energy) < dif):
			dif = abs(1 - en/energy)
			b11_energy = en
			
	cd106_energy = None
	dif = float('inf')
	for en in LIBRARY[x]['cd106'].keys():
		if(abs(1 - en/energy) < dif):
			dif = abs(1 - en/energy)
			cd106_energy = en
			
	cd110_energy = None
	dif = float('inf')
	for en in LIBRARY[x]['cd110'].keys():
		if(abs(1 - en/energy) < dif):
			dif = abs(1 - en/energy)
			cd110_energy = en
			
	cd111_energy = None
	dif = float('inf')
	for en in LIBRARY[x]['cd111'].keys():
		if(abs(1 - en/energy) < dif):
			dif = abs(1 - en/energy)
			cd111_energy = en
			
	cd112_energy = None
	dif = float('inf')
	for en in LIBRARY[x]['cd112'].keys():
		if(abs(1 - en/energy) < dif):
			dif = abs(1 - en/energy)
			cd112_energy = en
			
	cd113_energy = None
	dif = float('inf')
	for en in LIBRARY[x]['cd113'].keys():
		if(abs(1 - en/energy) < dif):
			dif = abs(1 - en/energy)
			cd113_energy = en
			
	cd114_energy = None
	dif = float('inf')
	for en in LIBRARY[x]['cd114'].keys():
		if(abs(1 - en/energy) < dif):
			dif = abs(1 - en/energy)
			cd114_energy = en
			
	cd116_energy = None
	dif = float('inf')
	for en in LIBRARY[x]['cd116'].keys():
		if(abs(1 - en/energy) < dif):
			dif = abs(1 - en/energy)
			cd116_energy = en
			
	fe56_energy = None
	dif = float('inf')
	for en in LIBRARY[x]['fe56'].keys():
		if(abs(1 - en/energy) < dif):
			dif = abs(1 - en/energy)
			fe56_energy = en
	
	si28_energy = None
	dif = float('inf')
	for en in LIBRARY[x]['si28'].keys():
		if(abs(1 - en/energy) < dif):
			dif = abs(1 - en/energy)
			si28_energy = en
			
	ca40_energy = None
	dif = float('inf')
	for en in LIBRARY[x]['ca40'].keys():
		if(abs(1 - en/energy) < dif):
			dif = abs(1 - en/energy)
			ca40_energy = en
	
	Sigs = {}
	
	Sigs['lead'] = {}
	tot = 0.241*((11.34*Na)/207.2)*(LIBRARY[x]['pb206'][pb206_energy]) + 0.221*((11.34*Na)/207.2)*(LIBRARY[x]['pb207'][pb207_energy]) + 0.524*((11.34*Na)/207.2)*(LIBRARY[x]['pb208'][pb208_energy])
	Sigs['lead'][206] = (0.241*((11.34*Na)/207.2)*(LIBRARY[x]['pb206'][pb206_energy]))/tot
	Sigs['lead'][207] = (0.221*((11.34*Na)/207.2)*(LIBRARY[x]['pb207'][pb207_energy]))/tot
	Sigs['lead'][208] = 0.524*((11.34*Na)/207.2)*(LIBRARY[x]['pb208'][pb208_energy])/tot
	
	Sigs['water'] = {}
	tot = 2*((1*Na)/18.0153)*LIBRARY[x]['h1'][h1_energy] + ((1*Na)/18.0153)*LIBRARY[x]['o16'][o16_energy]
	Sigs['water'][1] = (2*((1*Na)/18.0153)*LIBRARY[x]['h1'][h1_energy])/tot
	Sigs['water'][16] = (((1*Na)/18.0153)*LIBRARY[x]['o16'][o16_energy])/tot
	
	Sigs['B4C'] = {}
	tot = 4*0.199*((2.51*Na)/55.255)*LIBRARY[x]['b10'][b10_energy] + 4*0.801*((2.51*Na)/55.255)*LIBRARY[x]['b11'][b11_energy] + ((2.51*Na)/55.255)*LIBRARY[x]['c12'][c12_energy]
	Sigs['B4C'][10] = (4*0.199*((2.51*Na)/55.255)*LIBRARY[x]['b10'][b10_energy])/tot
	Sigs['B4C'][11] = (4*0.801*((2.51*Na)/55.255)*LIBRARY[x]['b11'][b11_energy])/tot
	Sigs['B4C'][12] = (((2.51*Na)/55.255)*LIBRARY[x]['c12'][c12_energy])/tot
	
	Sigs['cadmium'] = {}
	tot = 0.0125*((8.56*Na)/112.411)*LIBRARY[x]['cd106'][cd106_energy] + 0.1249*((8.56*Na)/112.411)*LIBRARY[x]['cd110'][cd110_energy] + 0.1280*((8.56*Na)/112.411)*LIBRARY[x]['cd111'][cd111_energy] + 0.2413*((8.56*Na)/112.411)*LIBRARY[x]['cd112'][cd112_energy] + 0.1222*((8.56*Na)/112.411)*LIBRARY[x]['cd113'][cd113_energy] + 0.2873*((8.56*Na)/112.411)*LIBRARY[x]['cd114'][cd114_energy] + 0.0749*((8.56*Na)/112.411)*LIBRARY[x]['cd116'][cd116_energy]
	Sigs['cadmium'][106] = (0.0125*((8.56*Na)/112.411)*LIBRARY[x]['cd106'][cd106_energy])/tot
	Sigs['cadmium'][110] = (0.1249*((8.56*Na)/112.411)*LIBRARY[x]['cd110'][cd110_energy])/tot
	Sigs['cadmium'][111] = (0.1280*((8.56*Na)/112.411)*LIBRARY[x]['cd111'][cd111_energy])/tot
	Sigs['cadmium'][112] = (0.2413*((8.56*Na)/112.411)*LIBRARY[x]['cd112'][cd112_energy])/tot
	Sigs['cadmium'][113] = (0.1222*((8.56*Na)/112.411)*LIBRARY[x]['cd113'][cd113_energy])/tot
	Sigs['cadmium'][114] = (0.2873*((8.56*Na)/112.411)*LIBRARY[x]['cd114'][cd114_energy])/tot
	Sigs['cadmium'][116] = (0.0749*((8.56*Na)/112.411)*LIBRARY[x]['cd116'][cd116_energy])/tot
	
	Sigs['portland cement'] = {}
	tot = 0.253*((3.15*Na)/(0.253*60.0843 + 0.747*56.077))*(LIBRARY[x]['si28'][si28_energy] + 2*LIBRARY[x]['o16'][o16_energy]) + 0.747*((3.15*Na)/(0.253*60.0843 + 0.747*56.077))*(LIBRARY[x]['ca40'][ca40_energy] + LIBRARY[x]['o16'][o16_energy])
	Sigs['portland cement'][16] = (((3.15*Na)/(0.253*60.0843 + 0.747*56.077))*(0.253*2*LIBRARY[x]['o16'][o16_energy] + 0.747*LIBRARY[x]['o16'][o16_energy]))/tot
	Sigs['portland cement'][28] = ((0.253*((3.15*Na)/(0.253*60.0843 + 0.747*56.077)))*LIBRARY[x]['si28'][si28_energy])/tot
	Sigs['portland cement'][40] = ((0.747*((3.15*Na)/(0.253*60.0843 + 0.747*56.077)))*LIBRARY[x]['ca40'][ca40_energy])/tot
	
	Sigs['iron'] = {}
	Sigs['iron'][56] = 1
	
	Sigs['iron poly'] = {}
	tot = ((2.681*Na)/97.90)*LIBRARY[x]['fe56'][fe56_energy] + 3*((2.681*Na)/97.90)*LIBRARY[x]['c12'][c12_energy] + 6*((2.681*Na)/97.90)*LIBRARY[x]['h1'][h1_energy]
	Sigs['iron poly'][56] = (((2.681*Na)/97.90)*LIBRARY[x]['fe56'][fe56_energy])/tot
	Sigs['iron poly'][12] = (3*((2.681*Na)/97.90)*LIBRARY[x]['c12'][c12_energy])/tot
	Sigs['iron poly'][1] = (6*((2.681*Na)/97.90)*LIBRARY[x]['h1'][h1_energy])/tot
	
	set = Sigs[compound]
	Sigs = []
	for i in set.keys():
		Sigs.append(set[i])
	
	x = 0
	probs = [0]
	for i in range(0,len(Sigs)):
		Sig = Sigs[i]
		p = Sig/(sum(Sigs))
		probs.append(p + sum(probs))
	
	for i in range(0,len(probs)):
		probs[i] = probs[i]/probs[len(probs)-1]
	
	for i in range(0,len(probs)-1):
		if((r >= probs[i]) and (r < probs[i+1])):
			x = i
	
	A = set.keys()[x]
	time2 = time.time()
	
	return A

#-----------------------------------------------------------------------------------------

#Define a function to calculate the scattering energy
def scattering_energy(init_energy,A,mu):
	energy = (init_energy/((1+A)**2))*((mu+sqrt((A**2)+(mu**2)-1))**2)
	
	return energy

#MAIN
#-----------------------------------------------------------------------------------------

#Define the tally
TALLY = []

#Define the source wall, bottom of wall is origin
wall_bottom = 0
wall_top = 10.0

#Define the geometry
#	GEOMETRY[x] = [[[lower left],[lower right],[upper left],[upper right]],compound]
GEOMETRY = {}

#Materials
#GEOMETRY[0] = [[[-36.0,0.0],[0.0,0.0],[-36.0,10.0],[0.0,10.0]],'void']
GEOMETRY[1] = [[[0.0,0.0],[1.0,0.0],[0.0,10.0],[1.0,10.0]],'B4C']
GEOMETRY[2] = [[[1.0,0.0],[10.0,0.0],[1.0,10.0],[10.0,10.0]],'lead']
GEOMETRY[3] = [[[-36.0,-20.0],[130.0,-20.0],[-36.0,0.0],[130.0,0.0]],'portland cement']
GEOMETRY[4] = [[[-36.0,10.0],[130.0,10.0],[-36.0,50.0],[130.0,50.0]],'portland cement']
GEOMETRY[5] = [[[-36.0,0.0],[0.0,0.0],[-36.0,10.0],[0.0,10.0]],'void']
GEOMETRY[6] = [[[10.0,0.0],[32.5,0.0],[10.0,3.0],[32.5,3.0]],'iron']
GEOMETRY[7] = [[[10.0,3.0],[32.5,3.0],[10.0,7.0],[32.5,7.0]],'void']
GEOMETRY[8] = [[[10.0,7.0],[32.5,7.0],[10.0,10.0],[32.5,10.0]],'iron']
GEOMETRY[9] = [[[32.5,0.0],[32.6,0.0],[32.5,10.0],[32.6,10.0]],'cadmium']
GEOMETRY[10] = [[[32.6,0.0],[65.0,0.0],[32.6,3.5],[65.0,3.5]],'iron']
GEOMETRY[11] = [[[32.6,3.5],[65.0,3.5],[32.6,6.5],[65.0,6.5]],'void']
GEOMETRY[12] = [[[32.6,6.5],[65.0,6.5],[32.6,10.0],[65.0,10.0]],'iron']
GEOMETRY[13] = [[[65.0,0.0],[65.1,0.0],[65.0,10.0],[65.1,10.0]],'cadmium']
GEOMETRY[14] = [[[65.1,0.0],[97.5,0.0],[65.1,4.0],[97.5,4.0]],'iron']
GEOMETRY[15] = [[[65.1,4.0],[97.5,4.0],[65.1,6.0],[97.5,6.0]],'void']
GEOMETRY[16] = [[[65.1,6.0],[97.5,6.0],[65.1,10.0],[97.5,10.0]],'iron']
GEOMETRY[17] = [[[97.5,0.0],[97.6,0.0],[97.5,10.0],[97.6,10.0]],'cadmium']
GEOMETRY[18] = [[[97.6,0.0],[130.0,0.0],[97.6,4.5],[130.0,4.5]],'iron']
GEOMETRY[19] = [[[97.6,4.5],[130.0,4.5],[97.6,5.5],[130.0,5.5]],'void']
GEOMETRY[20] = [[[97.6,5.5],[130.0,5.5],[97.6,10.0],[130.0,10.0]],'iron']
GEOMETRY[21] = [[[130.0,0.0],[138.0,0.0],[130.0,10.0],[138.0,10.0]],'void']
GEOMETRY[22] = [[[130.0,-20.0],[138.0,-20.0],[130.0,0.0],[138.0,0.0]],'iron']
GEOMETRY[23] = [[[130.0,10.0],[138.0,10.0],[130.0,50.0],[138.0,50.0]],'iron']
GEOMETRY[24] = [[[138.0,-20.0],[218.0,-20.0],[138.0,0.0],[218.0,0.0]],'iron poly']
GEOMETRY[25] = [[[138.0,10.0],[218.0,10.0],[138.0,50.0],[218.0,50.0]],'iron poly']
GEOMETRY[26] = [[[138.0,0.0],[148.0,0.0],[138.0,3.5],[148.0,3.5]],'iron']
GEOMETRY[27] = [[[138.0,3.5],[148.0,3.5],[138.0,6.5],[148.0,6.5]],'void']
GEOMETRY[28] = [[[138.0,6.5],[148.0,6.5],[138.0,10.0],[148.0,10.0]],'iron']
GEOMETRY[29] = [[[148.0,0.0],[162.0,0.0],[148.0,10.0],[162.0,10.0]],'void']
GEOMETRY[30] = [[[162.0,0.0],[218.0,0.0],[162.0,3.5],[218.0,3.5]],'iron poly']
GEOMETRY[31] = [[[162.0,3.5],[218.0,3.5],[162.0,6.5],[218.0,6.5]],'void']
GEOMETRY[32] = [[[162.0,6.5],[218.0,6.5],[162.0,10.0],[218.0,10.0]],'iron poly']
GEOMETRY[33] = [[[218.0,-20.0],[230.0,-20.0],[218.0,3.5],[230.0,3.5]],'lead']
GEOMETRY[34] = [[[218.0,3.5],[230.0,3.5],[218.0,6.5],[230.0,6.5]],'void']
GEOMETRY[35] = [[[218.0,6.5],[230.0,6.5],[218.0,50.0],[230.0,50.0]],'lead']

#Graveyard
GEOMETRY[36] = [[[-1*INF,-20.0],[-36.0,-20.0],[-1*INF,50.0],[-36.0,50.0]],'grave']
GEOMETRY[37] = [[[-1*INF,-1*INF],[INF,-1*INF],[-1*INF,-20.0],[INF,-20.0]],'grave']
GEOMETRY[38] = [[[-1*INF,50.0],[INF,50.0],[-1*INF,INF],[-1*INF,INF]],'grave']
GEOMETRY[39] = [[[230.0,6.5],[INF,6.5],[230.0,INF],[INF,INF]],'grave']
GEOMETRY[40] = [[[230.0,-1*INF],[INF,-1*INF],[230.0,3.5],[INF,3.5]],'grave']

#Tally
GEOMETRY[41] = [[[230.0,3.5],[INF,3.5],[230,6.5],[INF,6.5]],'tally'] #tally (right)

#Define energy groups and number of particles to run in each group
energy_groups = [0.25,0.5,0.75,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0] #MeV
n_particles = [0.0,(3.25 * 10**4),(4.0 * 10**4),(3.5 * 10**4),(3 * 10**4),(2.25 * 10**4),(1.75 * 10**4),(1.5 * 10**4),(1.0 * 10**4),(8.0 * 10**3),(5.5 * 10**3),(4.0 * 10**3),(2.75 * 10**3),(2.0 * 10**3),(1.5 * 10**3),(1.0 * 10**3),(6.0 * 10**2),(4.5 * 10**2)]
scale = 1
for i in range(0,len(n_particles)):
	n_particles[i] = int(scale*n_particles[i])
print('# of particles: ' + str(sum(n_particles)))

#Run each particle
print('PARTICLES RUNNING')
timeinit = time.time()
for i in range(0,len(energy_groups)):
	for j in range(0,n_particles[i]):
		print(str(i)+'--'+str(j))
		#Start the particle somewhere on the source wall
		energy = energy_groups[i]
		last_energy = energy
		last_position = [0,0]
		next_last_position = [0,0]
		next_position = [0,0]
		position = [0.00001,random.random()*(wall_top-wall_bottom)]
		last_angle = 0.000000001
		angle = 0.000000001
		last_region = 1
		region = 1
		last_compound = 'water'
		compound = 'water'
		
		#Move the particle continually until tally, capture, or graveyard
		while(True):
			#Move the particle
			#print(position)
			l = P_path_inv(random.random(),get_Sig(compound,'el',energy)+get_Sig(compound,'g',energy))
			next_last_position[0] = position[0]
			next_last_position[1] = position[1]
			next_position[0] = cos(angle)*l + position[0]
			next_position[1] = sin(angle)*l + position[1]
			
			#Determine what material and region the particle is in
			next_region = None
			next_compound = None
			keys = GEOMETRY.keys()
			for k in range(0,len(keys)):
				check_region = GEOMETRY[keys[k]][0]
				if(in_region(next_position,check_region) == True):
					next_region = keys[k]
					next_compound = GEOMETRY[keys[k]][1]
			
			#Check whether the particle has left the current region, but has not entered the graveyard
			if(region == next_region):
				#If so accept the new position sample interaction
				last_position[0] = position[0]
				last_position[1] = position[1]
				position[0] = next_position[0]
				position[1] = next_position[1]
				last_region = region
				region = next_region
				last_compound = compound
				compound = next_compound
				Sigs = [get_Sig(compound,'el',energy),get_Sig(compound,'g',energy)]
				x = P_int_inv(random.random(),Sigs)
				if(x == 0):
					#The particle elastically scattered
					#Determine scattering angle
					A = P_A_inv(random.random(),'el',compound,energy)
					ang = P_el_inv(random.random(),90,A)
					mu = cos(ang)
					last_energy = energy
					energy = scattering_energy(last_energy,A,mu)
					angle += ang
					if(angle > pi):
						angle -= 2*pi
					if(angle < -1*pi):
						angle += 2*pi
				elif(x == 1):
					#The particle is captured stop tracking
					print('ABSORBED')
					break
			else:
				#If not move the particle to the boundary and continue
				left = vertical(next_last_position,next_position,GEOMETRY[region][0][0],GEOMETRY[region][0][2])
				right = vertical(next_last_position,next_position,GEOMETRY[region][0][1],GEOMETRY[region][0][3])
				upper = horizontal(next_last_position,next_position,GEOMETRY[region][0][2],GEOMETRY[region][0][3])
				lower = horizontal(next_last_position,next_position,GEOMETRY[region][0][0],GEOMETRY[region][0][1])
				
				if((left != None) and (((angle <= pi) and (angle > pi/2)) or ((angle > -1*pi) and (angle < -1*pi/2)))):
					last_position[0] = position[0]
					last_position[1] = position[1]
					position[0] = left[0]
					position[1] = left[1]
				if((right != None) and ((angle < pi/2) and (angle > -1*pi/2))):
					last_position[0] = position[0]
					last_position[1] = position[1]
					position[0] = right[0]
					position[1] = right[1]
				if((lower != None) and ((angle < 0) and (angle > -1*pi))):
					last_position[0] = position[0]
					last_position[1] = position[1]
					position[0] = lower[0]
					position[1] = lower[1]
				if((upper != None) and ((angle > 0) and (angle < pi))):
					last_position[0] = position[0]
					last_position[1] = position[1]
					position[0] = upper[0]
					position[1] = upper[1]
				
				#Check which region the new position is in
				keys = GEOMETRY.keys()
				#print(region)
				for k in range(0,len(keys)):
					if(keys[k] != region):
						check_region = GEOMETRY[keys[k]][0]
						if(on_region(position,check_region) == True):
							#print(keys[k])
							last_region = region
							region = keys[k]
							last_compound = compound
							compound = GEOMETRY[keys[k]][1]
							break
				
				#If the next region is the graveyard, stop tracking
				if(compound == 'grave'):
					print('GRAVEYARD')
					break
				#If in the tally, record, then stop tracking
				if(compound == 'tally'):
					print('TALLIED')
					TALLY.append(energy)
					break
timeend = time.time()

print('Time Elapsed: ' + str(timeend-timeinit) + ' s')

file = open('/Users/puptheblue/Dropbox/Berkeley_Courses/Spring_2015/NE_155/Project/TALLY.txt','w')
for i in range(0,len(TALLY)):
	file.write(str(TALLY[i])+'\n')
file.close()

#-----------------------------------------------------------------------------------------



#TEST ZONE
#-----------------------------------------------------------------------------------------

