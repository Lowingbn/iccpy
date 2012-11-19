"""
A file to store data labels used in Gadget snapshots
Each dictionary maps a block name (e.g. "RHO " or "Z) to a number of 
entries per particle (e.g. 1 or 12) and identifies which particles they apply to.

The data type is not given as it is assumed to be the floating point type
given in the gadget header (i.e. 32 bit or 64 bit).

POS, VEL, MASS and ID are not given as they are assumed to be in all gadget files.
 
"""


gas_only = (1,0,0,0,0,0)
gas_stars = (1,0,0,0,1,0)
allparts = (1,1,1,1,1,1)
stars_only = (0,0,0,0,1,0)

cecilia_labels = {
"POT ":(1, allparts),
"U   ":(1, gas_only),
"RHO ":(1, gas_only),
"NE  ":(1, gas_only),
"NH  ":(1, gas_only),
"HSML":(1, gas_only),
"SFR ":(1, gas_only),
"TSTP":(1, allparts),
"ACCE":(3, allparts),
"BFLD":(3, allparts),
"BFSM":(3, allparts),
"ROTB":(3, allparts),
"SRTB":(3, allparts),
"Zs  ":(1, gas_stars),
"Z   ":(12, gas_stars),
"INFO":(1, (0,0,0,0,0,0)),
"AGE ":(1, stars_only),
"iM  ":(1, stars_only),
"SLg ":(1, stars_only),
"HSMS":(1, stars_only),
"SUB ":(1, stars_only),
"HSBP":(1, (0,0,1,1,0,0)),
"TIPS":(9, allparts),
"DIPS":(36, allparts),
"CACO":(1, allparts),
"FLDE":(1, allparts),
"STDE":(1, allparts),
"PSDE":(1, allparts),
"ANRA":(3, allparts),
"LACA":(20, allparts),
"SHIN":(3, allparts),
"INDE":(1, allparts),
"EGYP":(1, gas_stars),
"EGYC":(1, gas_stars)}
