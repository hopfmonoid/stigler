# STIGLER

This is a program to select food ingrediants while satisfying nutritional
requirements (RDA - recommended daily allowance) while minimising either cost,
or minimising weight of ingrediants given cost constraint, or maximising
diversity of ingredients given constraints on cost and weight of ingrediants.

## History

The program was inspired by George Stigler, who first solved (manually) in 1939
the problem of selecting food ingrediants to meet nutritional
requirements. Stigler solved the problem (Stigler's 1939 diet) for 5 foods
(wheat flour, evaporated milk, cabbage, spinach, dried navy beans) of 9
nutrients (calories, protein, calcium, iron, vitamin A, thiamine (vitamin B1),
riboflavin (vitamin B2), niacin, ascorbic acid (vitamin C)), and obtained anual
cost of $39.93. Later the problem was formulated and solved by George Dantzig
exactly using Simplex algorithm for linear programming, and obtained a solution
with anual cost $39.69 (using 1939 data). (See wikipedia entry on Stigler's
diet)


## User profile

  The program supports creating a user profile. It consists of the following
  files, which must be kept in a profile directory (myprofile is a sample
  profile):

- cost.csv - food ingrediants and their costs

- rda.csv - user's recommended daily allowance

- wt-groups.csv - here user can set constraints on the weights of ingrediants in
  each food group (the groups are vegetables, oils, fruits, beans, etc.)

- wt.csv - here user can set constraints on the weights of individual food
  ingrediants.

## Data directory:

- fxn.csv - food ingrediants (rows) vs nutritional composition (columns). The
  file nxf.csv contains the same data with rows and columns interchanged. This
  data is compiled from the data available on the USDA (United States Department
  of Agriculture) website. Currently 67 food ingrediants and 22 nutrients are
  supported. (The scripts used for compiling this data will be available in the
  future.)

- rda-vegan.csv, rda-dairy.csv - sample recommended daily allowance for vegan
  and vegetarian (with dairy) diet of adult male. This data is compiled from the
  NIH (National Institute of Health) website. You may want to copy one of the
  files in your profile, or ceate rda.csv of your own.

## Usage:

python3 stigler.py [-h] [-c COST_HI] [-w WEIGHT_HI] [-s SCALE] -p <profile_dir>
where

- COST_HI = upper limit on the cost

- WEIGHT_HI = upper limit on weight of
ingrediants

- SCALE = scaling factor for RDA values, with default 1.0 (for example, you may
want to plan the diet for a week, in which case use scale 7).

## Warning

The program was created for fun, and I do not recommend you to plan your diet
using this program, without colsulting your doctor. The author takes no
responsibility for the consequences of using the program.
