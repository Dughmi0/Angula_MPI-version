#you should add angula_mpi.for to your ANGULA\ code directory in before compilation as it uses the same modules angula.for uses.
#COMPILATION:
#you can easily compile it using mpif90 using the command:
      mpif90 -ffixed-line-length-none -std=legacy -fallow-argument-mismatch angula_mpi_optionB.for -o angula_mpi.exe
#SPOILER ALERT:
#you will get a warnning but it wont matter as I compilation will be successful; as my fortran knowledge is very limited I welcome any editing suggestions to fix unharming issue :').


#to change the maximum number if atoms per molecule and number of molecules per molecule type you should change it inside all files in the traductor directory and angula directory that contains dimat and dimmol and also diEat
#dimat is the number of maximum atoms per molecule
#dimmol is the maximum number of molecules per molecule type
#diEat is the maximum number of atomic energies affiliated with dimat
