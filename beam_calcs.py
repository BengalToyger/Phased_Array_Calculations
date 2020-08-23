#!/usr/bin/python3

import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

# Assuming a linear array of elements, calculates the phase 
# shift difference in degrees between the elements 
# (d_phi) needed to 
# steer the beam theta_s degrees from a direction normal
# to the array. 
def phase_shift_calc(d_div_lam, theta_s):
	d_phi = 360*np.sin(theta_s*np.pi/180)*d_div_lam
	return d_phi

# Calculates the Euclidean distance between a source and destination
def calc_distance(source, destination):
	r = np.sqrt((destination[0]-source[0])**2\
	+(destination[1]-source[1])**2\
	+(destination[2]-source[2])**2)
	return r

# Calculates the phase and magnitude of a wave at a distance r from an isotropic 
# source assuming
# amplitude A and an initial phase . 
def calc_wave(source, destination, phase, k, A):
	r = calc_distance(source, destination)
	wave = A*np.exp((k*r*2*np.pi+phase*np.pi/180)*1j)#/r
	return wave

# Convert spherical coorinates to rectangular coordinates
# azi is the angle in the XY plane relative to the X axis
# elev is the angle above or below the XY
def sphere_to_rect_xy(azi, elev, R):
	x = R*np.cos(elev)*np.cos(azi)
	y = R*np.cos(elev)*np.sin(azi)
	z = R*np.sin(elev)
	dest = np.array([x, y, z])
	return dest

# Convert spherical coorinates to rectangular coordinates
# azi is the angle in the XY plane relative to the X axis
# elev is the angle relative to the z axis
def sphere_to_rect_std(azi, elev, R):
	x = R*np.sin(elev)*np.cos(azi)
	y = R*np.sin(elev)*np.sin(azi)
	z = R*np.cos(elev)
	dest = np.array([x, y, z])
	return dest

# Calculates the array factor at azi, elev, and a distance R from the origin
# using the given sources, which should have members A (wave amplitude),
# phase (source starting phase), and coords (np.array of x y z coords)
def calc_af_at_dest(azi, elev, R, sources):
	AF_c = np.zeros((len(azi),len(elev)),dtype=complex)
	AF = np.zeros((len(azi),len(elev)))
	dest = sphere_to_rect_xy(azi,elev,R)
	for w in sources:
		AF_c=AF_c+calc_wave(w['coords'],dest,w['phase'],k,\
		w['A'])
	AF=np.abs(AF_c)

	norm_f = np.amax(AF)
	AF = AF/norm_f
	
	return AF

# Takes in any real phase value and places it between 0 and 360 (exclusive)
def wrap_phase_deg(phi):
	phi = phi % 360 # Put in range -360 to 360
	if phi < 0:
		phi = 360 + phi # Wrap 0 to 360
	phi = phi % 360 # Small negative numbers result in 360 degree phase. 
			# Phase shifters can only take 0 to 353. This makes any
			# 360 a 0
	return phi

# Takes in a spherical coordinate system and patch antenna parameters and 
# calculates the the antenna pattern. elev is the elevation, azi is the azimuth,
# W is width, L is length, k is wavenumber in free space. This function appears
# to be correct. It is found on 
# http://www.antenna-theory.com/antennas/patches/antenna.php
def patch_ant_pat(azi, elev, W, L, k):
	e_elev = e_elev_calc(azi, elev, W, L, k)
	e_azi = e_azi_calc(azi, elev, W, L, k)
	patch_pattern = np.sqrt(np.power(e_elev,2)+np.power(e_azi,2))
	norm_p = np.amax(patch_pattern)
	patch_pattern = patch_pattern/norm_p
	return patch_pattern

def e_elev_calc(azi, elev, W, L, k):
	num = num_calc(azi, elev, W, L, k)
	denom = denom_calc(azi, elev, W, L, k)
	fact = fact_calc(azi, elev, W, L, k)*np.cos(azi)
	return num/denom*fact

def e_azi_calc(azi, elev, W, L, k):
	num = -1*num_calc(azi, elev, W, L, k)
	denom = denom_calc(azi, elev, W, L, k)
	fact = fact_calc(azi, elev, W, L, k)*np.cos(elev)*np.sin(azi)
	return num/denom*fact

def num_calc(azi, elev, W, L, k):
	return np.sin(k*W*np.sin(elev)*np.sin(azi)/2)

def denom_calc(azi, elev, W, L, k):
	return k*W*np.sin(elev)*np.sin(azi)/2

def fact_calc(azi, elev, W, L, k):
	return np.cos(k*L*np.sin(elev)*np.cos(azi)/2)

# If running the module as a script, run the test code
if __name__ == '__main__':
	# Set up consts
	R = 1000 # Far-field radius
	freq = 10*10**9
	v_p = 3*10**8
	wave_l = v_p/freq
	k = 1/wave_l
	d_e = wave_l/2 # Distance between elements

	main_beam_azi = 45
	main_beam_elev = 90

	main_beam_azi = input("Please enter the input main beam azimuth in degrees\n")
	main_beam_elev = input("Please enter the input main beam elevation in degrees\n")

	main_beam_azi = float(main_beam_azi)
	main_beam_elev = float(main_beam_elev)

	# Translate azi and elev to unit vector in rectangular coords
	main_beam_vec = sphere_to_rect_xy(main_beam_azi*np.pi/180, main_beam_elev*np.pi/180, 1)

	print(main_beam_vec)

	# Determine corresponding rotations in zx and xy planes
	theta_zx = np.arctan(np.array(main_beam_vec[0])/\
		(main_beam_vec[2]+10**(-10)))*180/np.pi
	theta_yz = np.arctan(np.array(main_beam_vec[1])/\
		(main_beam_vec[2]+10**(-10)))*180/np.pi

	print(theta_zx)
	print(theta_yz)

	# Calculate the phase shift required to produce the rotation in each plane
	d_phi_zx = phase_shift_calc(d_e/wave_l,theta_zx)
	print(d_phi_zx)

	d_phi_yz = phase_shift_calc(d_e/wave_l,theta_yz)
	print(d_phi_yz)
	
	# Sum the phase shifts in correct planes
	phi_11 = d_phi_zx+d_phi_yz
	phi_01 = d_phi_yz
	phi_10 = d_phi_zx

	# Set up antenna sources
	sources = []
	sources.append(dict(\
		coords=np.array([d_e/2,d_e/2,0]),\
		phase=phi_11,\
		A=1\
	))


	sources.append(dict(\
		coords=np.array([-d_e/2,d_e/2,0]),\
		phase=phi_01,\
		A=1\
	))


	sources.append(dict(\
		coords=np.array([d_e/2,-d_e/2,0]),\
		phase=phi_10,\
		A=1\
	))


	sources.append(dict(\
		coords=np.array([-d_e/2,-d_e/2,0]),\
		phase=0,\
		A=1\
	))


	print(sources)
	print(wave_l)
	print(k)

	# Set up coordinate system
	num_azi = 200
	num_elev = 200

	azi = np.linspace(-np.pi, np.pi, num_azi)
	elev = np.linspace(-np.pi/2, np.pi/2, num_elev)
	azi_v, elev_v = np.meshgrid(azi, elev)
	
	# Calculate array factor
	AF = calc_af_at_dest(azi_v, elev_v, R, sources)
	

	fig = plt.figure()
	ax = plt.axes(projection='3d')
	ax.plot_wireframe(azi_v/np.pi*180, elev_v/np.pi*180, AF, label='Array Factor')
	ax.set_xlabel('Azimuth (degrees)')
	ax.set_ylabel('Elevation (degrees)')
	ax.set_zlabel('Array Factor')

	factor = sphere_to_rect_xy(azi_v, elev_v, AF)

	fig = plt.figure()
	ax = plt.axes(projection='3d')
	ax.plot_wireframe(factor[0], factor[1], factor[2], label='Array Pattern')
	ax.set_xlabel('X')
	ax.set_ylabel('Y')
	ax.set_zlabel('Z')

	W = 0.008
	L = 0.00932
	
	azi = np.linspace(-np.pi, np.pi, num_azi)
	elev = np.linspace(0.001, np.pi, num_elev)
	azi_v, elev_v = np.meshgrid(azi, elev)
	
	# Calculate antenna pattern
	AP = patch_ant_pat(azi_v, elev_v, W, L, 2*np.pi*k)
	
	patch_pattern = sphere_to_rect_std(azi_v, elev_v, AP)

	fig = plt.figure()
	ax = plt.axes(projection='3d')
	ax.plot_wireframe(patch_pattern[0], patch_pattern[1], patch_pattern[2], label='Array Pattern')
	ax.set_xlabel('X')
	ax.set_ylabel('Y')
	ax.set_zlabel('Z')
	
	plt.show()

