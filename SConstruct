env = Environment()

env['CXX'] = 'g++'
env.Append(CCFLAGS=['-pthread', '-O3', '-Wno-inline', '-std=c++0x'],
           CPPPATH=['/usr/include', '/usr/include/eigen3', './src' , './Liquids'],
           LIBS=['cantera', 'sundials_cvodes', 'sundials_ida', 'sundials_nvecserial', 'lapack', 'blas'],
           LIBPATH=['/usr/lib'],
           LINKFLAGS=['-pthread'],
           FRAMEWORKS=[])

program = env.Program('./sprayFlamelet2way', ['./sprayFlamelet2way.cpp', './src/Lagrangian.cpp', './src/SprayStFlow.cpp', './src/Sim1D.cpp', './src/OneDim.cpp','./src/boundaries1D.cpp', './Liquids/Ethanol.cpp'])
Default(program)
