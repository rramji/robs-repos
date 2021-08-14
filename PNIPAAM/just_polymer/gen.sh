#!/bin/bash

save=$2
echo $save

get_qeq_charges()
{
	${scripts_dir}/createLammpsInput.pl -b $bgf -f "DREIDING" -o pqeq -s test
	${scripts_dir}/../codes/bin/lmp_expanse -in in.test_singlepoint -screen none
	${scripts_dir}/convertLammpsTrj.pl -b $bgf -l test.lammps -t last -o bgf -s $bgf
	q=$(${scripts_dir}/bgfcharge.pl $bgf | awk '{print $NF}')
	if [ $(echo $q | awk '{if($1==0)print 0; else print 1}') -gt 0 ]; then
		${scripts_dir}/modifyAtomData.pl -s $bgf -w $save -a "index==1" -f "CHARGE:${q}"
	fi
	rm -fr in.test* test.* log.* data.test
}

if [ $# -lt 2 ]; then
	echo "usage: $0 polymer_length(#monomers)  save_name"
	exit 1
fi
scripts_dir=/home/tpascal/scripts/
module purge
module load cpu/0.15.4  gcc/10.2.0 intel-mkl openbabel gsl

echo $1 $2

echo "Step 2: Creating polymer structure"
#create polymer
${scripts_dir}/polymerize.pl -b nipaam.corrected.bgf -d x -n $1 -h "1 4" -t "2 7" -q 1 -s pnipaam.bgf -f DREIDING

bgf="pnipaam.bgf"
get_qeq_charges

#echo "Step 6: Solvating"
##solvate
#${scripts_dir}/addSolvent.pl -i SiH.pnipaam.dry.bgf -f "DREIDING spcew" -n "density:1.0" -w spc -s SiH.pnipaam.solv.bgf
##${scripts_dir}/addSolvent.pl -i SiH.pnipaam.dry.bgf -f "DREIDING spcew" -n "density:1.0" -w spc -s SiH.pnipaam.solv.bgf > /dev/null
#bot=$(${scripts_dir}/getBounds.pl -b SiH.pnipaam.solv.bgf -o "fftype eq 'Si3'" | grep '^Z ' | awk '{print $2}')
#top=$(${scripts_dir}/getBounds.pl -b SiH.pnipaam.solv.bgf -o "resname ne 'WAT'" | grep '^Z ' | awk '{print $3+5}')
#${scripts_dir}/getBGFAtoms.pl -b SiH.pnipaam.solv.bgf -s SiH.pnipaam.solv.bgf -o "resname ne 'WAT' or (zcoord>${bot} and zcoord<${top})" -m 1 > /dev/null
#${scripts_dir}/splitAtomsByMol.pl -b SiH.pnipaam.solv.bgf -s SiH.pnipaam.solv.bgf -f resnum 
#${scripts_dir}/centerBGF.pl -b SiH.pnipaam.solv.bgf -s SiH.pnipaam.solv.bgf -c com_center
echo "Done"
