# 2d SRD test: big + small particles

units           lj
#atom_style      atomic
atom_style      hybrid sphere dipole
atom_modify     first big
dimension       3
#kspace_style    ewald 1e-4

# create big particles with sigma 1.0




#lattice         sc 1.0
region          box block 0 32 0 32 0 32
#region          box1 block 16 32 16 32 16 32
region          box2 block 5 11 1 30 2 30
#create_box      2 box
#lattice         sc 0.016
#create_atoms    1 region box1
#oelete_atoms    porosity box 0.40 23529
#read_restart ../restart.equil
reset_timestep 0

pair_style lj/cut 4.8
#shell   gunzip data.srd.big2.gz
read_data     data.srd.$(v_replica) 
#read_data     ../data.srd.$(v_replica) 
#read_data     data.srd.big2
##change_box   all x scale 2 y scale 2 z scale 2  set
#read_data     data.srd add append shift 0 0 32
#read_data     data.srd add append shift 0 32 0
#read_data     data.srd add append shift 0 32 32
#read_data     data.srd add append shift 32 0 0
#read_data     data.srd add append shift 32 0 32
#read_data     data.srd add append shift 32 32 0
#read_data     data.srd add append shift 32 32 32
##create_atoms    2 random 163840 3928 box


group           big type 1
group           small type 2

# need both mass settings due to hybrid atom style

mass		1 1.0e-2
mass		2 4e-4
set             type 1 mass 1.0e-2
set             type 2 mass 4e-4

set             type 1 diameter $(v_dia)
set             type 2 diameter 0.0
variable      lambda  equal "sqrt(v_dia^3*v_lbd)"
#set             type 1 dipole/random $(v_ran1) $(v_lambda)
#delete_atoms   group small

print "hello"

velocity        small create 1.0 $(v_ran1) loop geom dist gaussian

velocity        big create 1.0 $(v_ran1) loop geom dist gaussian
# SRD run
#
#
reset_timestep 0
neighbor        0.0 bin
neigh_modify    include big  delay 0 every 1 check no 
#comm_modify     mode multi group big vel yes
comm_modify  group big vel yes


#variables
variable bigmass equal mass(big)
variable smallmass equal mass(small)
variable cut equal "v_dia*2^(1./6.)"

# no pairwise interactions with small particles

#pair_style      lj4824/cut/dipole/long 2.5
#pair_style lj/cut/dipole/cut 4.8 16
pair_style      lj/cut 5
pair_coeff      1 1 1. $(v_dia)  $(v_cut)

pair_coeff      2 2 0.0 1.0 0.0
pair_coeff      1 2 0 2.90226534287825 2.5

# must set 1-2 cutoff to non-zero value
delete_atoms overlap 2.5 small big
pair_modify     shift yes

delete_atoms overlap 0.5 small big
# use fix SRD to push small particles out from inside big ones

timestep        0.00005
#compute 15 big group/group small
compute rotTemp big temp/sphere dof all
compute 10 big property/atom tqx tqy tqz
#compute 11 big property/atom fx fy fz
compute 12 big reduce ave c_10[1] c_10[2] c_10[3]
compute 13 big reduce ave fx fy fz
variable sqToq equal "(c_12[1]*c_12[1]+c_12[2]*c_12[2]+c_12[3]*c_12[3])"
variable sqForce equal "(c_13[1]*c_13[1]+c_13[2]*c_13[2]+c_13[3]*c_13[3])"


fix     4 big ave/time 1 100 100 v_sqToq
fix     6 big ave/time 1 100 100 v_sqForce
variable D0r equal "dt^2*f_4/(6*dt)"
variable D0t equal "dt^2*f_6/(6*dt)"

fix             1 big nve/sphere update dipole
#fix             10 big  langevin 20 1  100.0 48279
fix             2 small srd/mcs 240 big 1.0 1 49894  radius 0.8  collision noslip inside warn overlap no exact yes shift yes $(v_ran1) tstat yes scale 0.1
#fix             2_3 small temp/rescale $(v_N) 1 1 0.001 1.01
fix             3 all momentum 100  linear 1 1 1 
##fix             2 big langevin 1 1 0.16666666  201409 zero no omega yes gjf yes

#fix             2 small srd 20 big 1.0 0.2 49894 radius 0.88 search 0.2 collision noslip
#fix             2 small srd 20 big 1.0 0.2 49894 radius 0.88 search 0.2 collision nonoslip inside warn
#fix             2 small srd 10 big 1.0 1. 49894  radius 0.8      search 0.2 collision nonoslip inside warn overlap yes exact no shift yes $(v_ran1) tstat yes rescale yes
#fix             F_output small print 100 "${bigmass} ${smallmass}"
#compute         srdMSD small msd
#compute         bigMSD big msd

#compute mu big property/atom mux muy muz
#compute magnet big reduce ave c_mu[1] c_mu[2] c_mu[3]
#variable magnet equal  "sqrt((c_magnet[1]*c_magnet[1]+c_magnet[2]*c_magnet[2]+c_magnet[3]*c_magnet[3]))/13.95"
#fix     magnet big ave/time 1 100 100  v_magnet
#compute 1 big group/group small pair yes kspace no
#compute 1 big group/group small pair yes kspace no
#diagnostics

compute   tbig big temp
compute   tsmall small temp
variable  pebig equal pe*atoms/count(big)
variable  ebig equal etotal*atoms/count(big)
#thermo_style  custom step temp  etotal v_pebig v_ebig press  c_bigMSD[4] time  v_D0r v_D0t f_magnet
#run 30000
#unfix 10
#fix             10 big  langevin 1 1  100.0 48279
#run 30000
#unfix 10
#run 10000
#thermo_style  custom time c_rotTemp  etotal v_pebig v_ebig press  c_bigMSD[4] v_D0r v_D0t f_magnet
thermo_style  custom time c_rotTemp  etotal v_pebig v_ebig press c_tsmall c_tbig f_2[9]

thermo_modify temp tbig
thermo    600


#thermo          100
#thermo_style custom step temp c_rotTemp  v_D0r v_D0t  ke pe etotal

#dump           1 all custom 100 fluid.lammpstrj id type xu yu zu
#dump           1 big custom 100 big.lammpstrj id type xu yu zu
print "dump D_all"
dump            D_all big custom 600 big.lammpstrj.$(v_replica) id type xu yu zu vx vy vz
#dump D_small small custom 10000 small.lammpstrj.$(v_replica) id type xu yu zu vx vy vz
print "pair coeff"
pair_coeff      1 2 0 2.90226534287825 2.5
#variable make loop 50
run             60000   every 120000 &
	"write_data data.srd.$(v_replica)"
print "done"
write_data data.srd.$(v_replica)

