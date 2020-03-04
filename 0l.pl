
##### should add the above code for debugging

#@remove =(1..20);

#foreach (@remove){
#$temp=$_*50;
system ("rmdir /s/q 00MC_00050_fcc");
#system ("rmdir /s/q SiNWoutput");

#}
#unlink ("00PropertyT.dat");
#system ('rmdir /s/q DPD_denopt');
#system ('rmdir /s/q testnvt');
#system ("lmp_mpi -in ZnO_box_MD.in"); # to generate a designed ZnONW with the size matching graphene box in Z dimension
#system ("lmp_mpi -in graphene_zx.in");
#system ("lmp_mpi -in NPT2.lmps");
#system ("lmp_mpi -sf omp -pk omp 8 -in DOPE_SiNW.in");
#system ("lmp_mpi -sf omp -pk omp 8 -in SiNW_DOPE.in");

#system ("lmp_mpi -sf omp -pk omp 16 -in ZnO_box_MD_readdump.in");

#system ('"C:\Program Files\MPICH2\bin\mpiexec.exe" lmp_mpi -sf omp -pk omp 16 -in 00Tension_300.in');

 
#system ('"C:\Program Files\MPICH2\bin\mpiexec.exe" -np 8 lmp_mpi -in DPD_denopt.in');
system ('lmp_mpi -in NPT2_1.lmps');
#system ('lmp_serial -in BHH-tfmc.in');
print "all done";

