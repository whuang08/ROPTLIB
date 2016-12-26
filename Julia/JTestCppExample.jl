# To run any test file in ROPTLIB/test in Julia,
# first generate a shared library by "make JuliaROPTLIB TP=name_of_a_test_file"
# in the terminal, e.g., "make JuliaROPTLIB TP=TestStieBrockett" for the problem
# defined in "TestStieBrockett.h". The shared library "TestStieBrockett.so" will
# be generated.

# add the shared library in ROPTLIB
Libdl.dlopen(path_to_lib * "/TestStieBrockett.so", Libdl.RTLD_GLOBAL)

# include the head file.
cxxinclude("TestStieBrockett.h")

# call the function declared in the head file.
@cxx main()


