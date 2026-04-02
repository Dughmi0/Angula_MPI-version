you should add angula_mpi.for to your ANGULA\ code directory in before compilation as it uses the same modules angula.for uses.
COMPILATION:
you can easily compile it using mpif90 using the command:
      mpif90 -ffixed-line-length-none -std=legacy -fallow-argument-mismatch angula_mpi_optionB.for -o angula_mpi.exe
SPOILER ALERT:
you will get a warnning but it wont matter as I compilation will be successful; as my fortran knowledge is very limited I welcome any editing suggestions to fix unharming issue :').
