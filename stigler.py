#!/usr/bin/python3

# Copyright (C) 2015-2020 Bhalchandra Digambar Thatte

# This file is a part of the project STIGLER. The program is released under the
# GPL V3 license.

# This is a program to select food ingrediants while satisfying nutritional
# requirements (RDA - recommended daily allowance) while minimising either cost,
# or minimising weight of ingrediants given cost constraint, or maximising
# diversity of ingredients given constraints on cost and weight of ingrediants.

# convex optimization library
from cvxopt import matrix, solvers
# solvers.options['show_progress'] = False
from cvxopt import printing
printing.options['dformat'] = '%.1f'
printing.options['width'] = -1

import sys
import numpy as np
np.set_printoptions(threshold=sys.maxsize)
import csv

import argparse
parser = argparse.ArgumentParser(description='A program to find minimum cost, minimum weight, maximum variety food composition with cost and weight constraints')
parser.add_argument("-c", "--cost_hi", type=float, default=100.0, help = "upper bound on cost (used in weight and variety optimisations)")
parser.add_argument("-w", "--weight_hi", type=float, default=100.0, help = "upper bound on weight (used in variety optimisations)")
parser.add_argument("-s", "--scale", type=float, default=1.0, help = "scaling factor for RDA values")
parser.add_argument("-p", "--profile", type=str,
                    help = "profile directory", required=True)

args = parser.parse_args()
profile = args.profile
scale = args.scale
cost_hi = args.cost_hi
wt_hi = args.weight_hi

# input files --------------------------------------------------------
# costs of foods
file_costs = "./"+profile+"/cost.csv"

# Nutritient data of foods: fxn - foods (rows) x nutrients (columns);
# nxf - nutrients (rows) x foods (columns);
file_nut_fxn = "./data/fxn.csv"
file_nut_nxf = "./data/nxf.csv"

# Weight constraints on ingredients and food groups

file_rda = "./"+profile+"/rda.csv"
file_food_weights = "./"+profile+"/wts.csv"
file_group_weights = "./"+profile+"/group-wts.csv"

# load files ---------------------------------------------------------
# numpy array: nutrients
nutrients = np.loadtxt(file_rda, dtype='str', delimiter=",", skiprows=2, usecols=(0,1))
# nn = nutrients.size
# print(nutrients[0:nn])

# numpy array: foods
foods = np.loadtxt(file_nut_fxn, dtype='str', delimiter=",", skiprows=2, usecols=(0))
nf=foods.size
# print(foods[0:nf])

# Cost vector
C = matrix(np.loadtxt(file_costs, delimiter=",", skiprows=1, usecols=(3)))

# Nutrition constraints (nh = upper bounds on nutrients, Nx <= nh; nl =
# lower bounds on nutrients, nl <= Nx, i.e., -Nx <= -nl; Read
# rda-(dairy/vegan).csv
N_np = np.loadtxt(file_nut_nxf, delimiter=",", skiprows=1, usecols=range(2,nf+2)) # as numpy array
N = matrix(N_np)
nh = matrix(np.loadtxt(file_rda, delimiter=",", skiprows=2, usecols=(3)))
nl = matrix(np.loadtxt(file_rda, delimiter=",", skiprows=2, usecols=(2)))

# Weight constraints on individual foods (wl = minimum weight, wl<= IX,
# i.e., -IX <= -wl; wh = maximum weight, IX <= wh)
I = matrix(np.eye(nf))
wl = matrix(np.loadtxt(file_food_weights, delimiter=",", skiprows=2, usecols=(1)))
wh = matrix(np.loadtxt(file_food_weights, delimiter=",", skiprows=2, usecols=(2)))

# Weight constraints for food groups (group_defs_matrix*X <= wgh = maximum weights by groups; minimum
# weights by groups = wgl <= group_defs_matrix*X, i.e., -group_defs_matrix*X <= -wg)

g = [4,4,7,8,9,4,14,17]
# incidence matrix with food groups on rows and ingredients on columns
W = np.zeros(shape=(len(g),nf),dtype=np.float64)
d = 0
for k in range(len(g)):
    for i in range(g[k]):
        W[k][d+i]=1
    d += g[k]
group_def_matrix = matrix(W)
# group_defs_matrix = matrix(np.loadtxt(file_group_defs, delimiter=",", skiprows=1, usecols=range(1,nf+1)))
wgh = matrix(np.loadtxt(file_group_weights, delimiter=",", skiprows=2, usecols=(2)))
wgl = matrix(np.loadtxt(file_group_weights, delimiter=",", skiprows=2, usecols=(1)))

# --------------------------------------------------------------------
# Find minimum cost food composition
# --------------------------------------------------------------------
# solve the linear program
# minimize: C'X
# subject to GX <= h, X >= 0

G = matrix([N,-N,-I,I])
h = scale*matrix([nh,-nl,-wl,wh])

# G = matrix([N,-N,-I,I,-group_def_matrix, group_def_matrix])
# h = matrix([nh,-nl,-wl,wh,-gwl,gwh])
X=solvers.lp(C,G,h)
# X=solvers.lp(C,G,h,solver='glpk')

# --------------------------------------------------------------------
# Find minimum weight food composition
# --------------------------------------------------------------------
# solve the linear probgam
# minimize: 1'*Y (i.e., total weight = y_1+y_2+ ...+y_{nf}, where 1 is a column vector of 1s, and 1' is its transpose)
# subject to: GY <= h, Y >= 0

D = matrix(np.ones((nf,1)))

G = matrix([N,-N,-I,I,-group_def_matrix,group_def_matrix,C.trans()])
h = scale*matrix([nh,-nl,-wl,wh,-wgl,wgh,cost_hi])
Y=solvers.lp(D,G,h)

# --------------------------------------------------------------------
# Find maximum variety food composition
# --------------------------------------------------------------------
# solve the quadratic program
# minimize: Z'PZ + QZ
# subject to GZ <= h, Z >= 0

P_np = np.zeros(shape=(nf,nf),dtype=np.float64)

# TODO - avoid hardcoding. Number of ingredients in each food group
# (dairy - 4, flours - 4, grains - 7, legumes - 8, nuts and seeds - 9,
# oils - 4, fruits - 14, vegetables - 17)
g = [4,4,7,8,9,4,14,17]
d = 0
for k in range(len(g)):
    for i in range(g[k]):
        for j in range(g[k]):
            if i == j:
                P_np[d+i][d+i] = N_np[0][d+i]*N_np[0][d+i]*(g[k]-1)/g[k]
            else:
                P_np[d+i][d+j] = -N_np[0][d+i]*N_np[0][d+j]/g[k]
    d += g[k]

# for i in range(nf):
#     for j in range(nf):
#         if i == j:
#             P_np[i][i] = N_np[0][i]*N_np[0][i]*(nf-1)/nf
#         else:
#             P_np[i][j] = -N_np[0][i]*N_np[0][j]/nf

P = matrix(P_np)

q_np = np.zeros(shape=(nf,1),dtype=np.float64)
q = matrix(q_np)

G = matrix([N,-N,-I,I,-group_def_matrix,group_def_matrix,C.trans(),D.trans()])
h = scale*matrix([nh,-nl,-wl,wh,-wgl,wgh,cost_hi,wt_hi])

# Z = solvers.qp(P,q,G,h,A=None,b=None)
Z = solvers.qp(P,q,G,h)
# --------------------------------------------------------------------

if not X['status'] == "optimal":
    print('Cost minimisation unsuccessful. Weight constraints may be restrictive')

if not Y['status'] == "optimal":
    print('Weight minimisation unsuccessful. Cost or weight constraints may be restrictive')

if not Z['status'] == "optimal":
    print('Variety maximisation unsuccessful. Cost or weight constraints may be restrictive')

# if not (X['status'] == "optimal" and Y['status'] == "optimal" and Z['status'] == "optimal"):
#     quit()

# reports
print("--------------------")
print("Ingredients in grams")
print("--------------------")
# print(X['x']*100.0)
print('%-20s'%"Food", '%-12s'%"Min. Cost", '%-12s'%"Min. Weight", '%-12s'%"Max. Variety" )
for name, v1, v2, v3 in zip(foods,X['x']*100.0,Y['x']*100.0,Z['x']*100.0): print ('%-20s'%name, '%-12.0f' %v1, '%-12.0f' %v2, '%-12.0f' %v3)
cx = C.T*X['x']
cy = C.T*Y['x']
cz = C.T*Z['x']
a = cx[0]
b = cy[0]
c = cz[0]
print('--------------------------------------------------------')
print('%-20s'%"Cost in Reais", '%-12.1f'%a, '%-12.1f'%b, '%-12.1f'%c)

D = matrix(np.ones((nf,1)))
dx = D.T*X['x']
dy = D.T*Y['x']
dz = D.T*Z['x']
a = dx[0]
b = dy[0]
c = dz[0]
print('%-20s'%"Total weight (x100g)", '%-12.1f'%a, '%-12.1f'%b, '%-12.1f'%c)

print("-----------------------")
print("Nutritional composition")
print("-----------------------")

MX = N*X['x']
MY = N*Y['x']
MZ = N*Z['x']
# print(N)
print('%-24s'%"Nutrient", '%-12s'%"RDA", '%-14s'%"Max Permitted", '%-14s'%"Min. Cost", '%-14s'%"Min. Weight", '%-14s'%"Max. Variety")
for nut, lo, hi, v1, v2, v3 in zip(nutrients, nl, nh, MX, MY, MZ): print ('%-24s'%nut, '%-12.1f' %lo, '%-14.1f' %hi, '%-14.1f' %v1, '%-14.1f' %v2, '%-12.1f' %v3 )

