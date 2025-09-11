

nx = 50
ny = nx
nz = ny
dx = 1.
total_cytokines = 13
total_celltypes = 14
boundaryat = 5
relaxationmcs = 1000

true_size =  1 #### was 5

fipy_duration = int(1) #in hours
s_mcs = 60.
h_mcs = 1/60.
true_mass = 1


#conversions 
lineconv = true_size/nx #cm per grid
areaconv = true_size**2/nx**2
volumeconv = (true_size**3*1)/(nx**3*1)
massconv = true_mass


# parameters for PAMPS
D_PAMPS = 3*10**-7*s_mcs/areaconv
mu_PAMPS = 0.6 * h_mcs
k_nec_PAMPS = 225 * 10**-5 * volumeconv * h_mcs


# parameters for DAMPS
D_DAMPS = 3*10**-7*s_mcs/areaconv
mu_DAMPS = 0.6 * h_mcs
k_nec_DAMPS = 300 * 10**-4.6 * volumeconv * h_mcs #k_nec_DAMPS = 300 * 10**-4.6 * volumeconv * h_mcs
k_neutro_nec_DAMPS = 225 * 10**-4.1 * volumeconv * h_mcs

# parameters for CCL2
D_ccl2 = 3*10**-7*s_mcs/areaconv
mu_ccl2 = 0.6 * h_mcs
k_endoccl2 = 225 * 10**-5 * volumeconv * h_mcs
k_kera_ccl2 = 225 * 10**-5 * volumeconv * h_mcs
k_fibro_ccl2 = 225 * 10**-5 * volumeconv * h_mcs

# parameters for il1a
D_il1a = 3*10**-7*s_mcs/areaconv
mu_il1a = 0.6 * h_mcs
k_m1_il1a = 225 * 10**-5 * volumeconv * h_mcs

# parameters for il1b
D_il1b = 3*10**-7*s_mcs/areaconv
mu_il1b = 0.6 * h_mcs
k_m1_il1b = 225 * 10**-5 * volumeconv * h_mcs
k_neut_il1b = 225 * 10**-5 * volumeconv * h_mcs

# parameters for il6
D_il6 = 8.49*10**-8*s_mcs/areaconv
mu_il6 = 0.5*h_mcs
k_neut_il6 = 250 * 10**-5 * volumeconv * h_mcs
k_mast_il6 = 250 * 10**-5 * volumeconv * h_mcs

# parameters for il8
D_il8 =  2.09*10**-6*s_mcs/areaconv
mu_il8 = 0.2*h_mcs
k_mastil8 = 234 * 10**-5 * volumeconv * h_mcs
k_neutroil8 = 1.46 * 10**-5 * volumeconv * h_mcs


# parameters for il10
D_il10 = 4.45*10**-8*s_mcs/areaconv
mu_il10 = 0.5*h_mcs
k_m2_il10 = 96 * 10**-5 * volumeconv * h_mcs


# parameters for il1RA
D_il1RA = 1.48*10**-8*s_mcs/areaconv
mu_il1RA = 0.73 *h_mcs
k_m2_il1RA = 24.25 * volumeconv * h_mcs

# parameters for tnf-a
D_tnfa = 4.07*10**-9*s_mcs/areaconv
mu_tnfa = 0.5*0.225*h_mcs
k_neut_tnfa = 250 * 10**-5 * volumeconv * h_mcs
k_mast_tnfa = 270 * 10**-5 * volumeconv * h_mcs

# parameters for tgf-b
D_tgfb = 1.3*10**-7*s_mcs/areaconv
mu_tgfb = 0.5 * 1/25 * h_mcs
k_m2_tgfb = 50 * 10**-5 * volumeconv * h_mcs
k_plat_tgfb = 60 * 10**-5 * volumeconv * h_mcs


# parameters for PDGF
D_PDGF = 2.2*10**-7*s_mcs/areaconv
mu_PDGF = 0.5 * 1/25 * h_mcs
k_m2_PDGF = 30 * 10**-5 * volumeconv * h_mcs
k_plat_PDGF = 30 * 10**-5 * volumeconv * h_mcs

# parameters for FGF
D_FGF = 2.1*10**-7*s_mcs/areaconv
mu_FGF = 0.5 * 1/25 * h_mcs
k_m2_FGF = 22 * 10**-5 * volumeconv * h_mcs
k_plat_FGF = 10 * 10**-5 * volumeconv * h_mcs







# cell activation probability/saturation coeff
lmrtnf = 0.5
lm1il10 = 1.0
lftgf = 0.5
lmrpamp = 0.5
lmril1ra = 0.5
cccl2 = 5*10**-9


# 280k total each hour is 1k mcs # for thesting out this system do the values /10 or /100
#this means that if a lifespan here is 240 this means that in 240 / 100 = 2.4 loops before this is initiated
#this is a bit unintuitive because the max amount of loops in the system is 280 which is 280 000 mcs



endothelial_lifespan = 1000000 # live months to years
kera_lifespan = 1000000 # live around 1.5 months
neutro_lifespan_6 = 600
neutro_lifespan_8 = 800
fibro_start_lifespan = 23000 # live weeks-months
fibro_end_lifespan = 28000 # live weeks-months


fibro_final_lifespan = 1000000

platelet_lifespan = 9600 # live 7-10 days # calculate basically 9600 = 9600/100 = 96 cycles where each cycle is 1000 long so 96000 mcs
mast_lifespan = 1000000 # dont die in the system
mono_lifespan = 1000000 # dont die in the system
m1_lifespan = 1000000 # dont die in the system
m2_lifespan = 1000000 # dont die in the system
necrotic_lifespan = 1000000 # they need to be removed by macrophages
necrotic_neut_lifspan = 1000000 #  they need to be removed by macrophages
apoptotic_neut_lifespan = 1000000 # they need to be removed by macrophages



base_necrotic_chance_neut = 0.01
apop_chance_neut = 0.01
base_cell_necro_chance = 0.01
clearance_probability_macro = 0.03
base_chance_kera_death = 0.01
base_kera_necro_chance = 0.01
cell_deletion_chance = 1

base_fibro_death_chance = 0.01
base_fibro_deletion_from_grid = 0.05

sigmoida = 1
sigmoidb = 4






cytokine_decay = 0.001
