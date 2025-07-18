#!/bin/bash
##########################

Out=results_Fe_BCC.txt
#Test temperaturess
TestTemp=(10.0 20.0 50.0 100.0 150.0 200.0 300.0 400.0 500.0 600.0 700.0 800.0 900.0 1000)
clear

echo "# Results of lattice parameters predictions at different temperatures" | tee -a $Out
echo "# In-Temp, lx, ly, lz, Reported-Temp" | tee -a $Out
for i in "${TestTemp[@]}"
do
  echo -n ${i}", ">> $Out
  cat > in_bcc.temp << EOF
clear
log log.py append
# ------------------------ INITIALIZATION ----------------------------
# ------------------------ PRE-AMBLE--------------------------
units		metal
dimension 	3
boundary	p	p	p
atom_style		atomic
atom_modify map array
neighbor		1.0 bin
neigh_modify	every 1 delay 2 check yes
variable	nl	equal	10
# ------------------------ DOMAIN DEFINITION--------------------------
lattice	bcc	2.866
region	region1	prism	0	\${nl}	0	\${nl}	0	\${nl}	 0 0 0
create_box	1	region1
create_atoms	1	region region1
# ------------------------ FORCE FIELDS ------------------------------
# Interatomic potential
pair_style		eam/alloy
pair_coeff		* * Fe_Ack.eam Fe
neighbor        		2 bin
neigh_modify    		every 10 delay 0 check no
mass            		1 56
# ------------------------- EQUILIBRATION-------------------------------
variable myTemp  equal temp
variable myPress equal press 
variable myVol   equal vol
variable natoms equal "count(all)" 
variable teng equal "pe"
variable ecoh equal "v_teng/v_natoms"
variable lengthx equal "lx/v_nl"
variable lengthy equal "ly/v_nl"
variable lengthz equal "lz/v_nl"

fix           avgmyTemp  all ave/time 1 5 10 v_myTemp ave running start 10000
fix           avgmyPress  all ave/time 1 5 10 v_myPress ave running start 10000
fix           avgmyVol  all ave/time 1 5 10 v_myVol ave running start 10000
fix           avglx  all ave/time 1 5 10 v_lengthx mode scalar ave running start 10000
fix           avgly  all ave/time 1 5 10 v_lengthy mode scalar ave running start 10000
fix           avglz  all ave/time 1 5 10 v_lengthz mode scalar ave running start 10000
variable      T equal "f_avgmyTemp"
variable      P equal "f_avgmyPress"
variable      V equal "f_avgmyVol"
variable      LX equal "f_avglx"
variable      LY equal "f_avgly"
variable      LZ equal "f_avglz"

timestep      0.001
thermo 100
thermo_style custom step pe lx ly lz temp press pxx pyy pzz vol


velocity all create $i 4928459 dist gaussian
fix 1 all npt temp $i $i 0.1 aniso 0.0 0.0 1 drag 1 
run 10000

min_style       cg
minimize        1.0e-25  1.0e-25  5000  10000

print "flag: lx = \${LX} ly = \${LY} lz = \${LZ} Energy = \${ecoh}"
print "avarage temperature = \${T}"
print "lx = \${LX}"
print "ly = \${LY}"
print "lz = \${LZ}"
print "Total energy (eV) = \${ecoh}"
print "Lattice constantx (Angstoms) = \${lengthx}"
print "Lattice constanty (Angstoms) = \${lengthy}"
print "Lattice constantz (Angstoms) = \${lengthz}"
print "Temperature (K) = \${myTemp}"
# write_data Fe_$i.data
EOF

  mpirun -np 8 lmp_mpi -in in_bcc.temp >> in_bcc.out

  otemp=$(grep "flag:" in_bcc.out | cut -d' ' -f 4)
	echo -n  $otemp", ">> $Out
  otemp=$(grep "flag:" in_bcc.out | cut -d' ' -f 7)
	echo -n  $otemp", ">> $Out
  otemp=$(grep "flag:" in_bcc.out | cut -d' ' -f 10)
	echo -n  $otemp", ">> $Out
  grep "flag:" in_bcc.out | cut -d' ' -f 13 >> $Out
  rm in_bcc.out
done
###################################
