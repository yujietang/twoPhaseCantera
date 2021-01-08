env = Environment()

env['CXX'] = 'g++'
env.Append(CCFLAGS=['-pthread', '-O3', '-Wno-inline', '-std=c++0x'],
           CPPPATH=['/opt/cantera/include', './src'],
           LIBS=['cantera', 'sundials_cvodes', 'sundials_ida', 'sundials_nvecserial', 'lapack', 'blas'],
           LIBPATH=['/opt/cantera/lib'],
           LINKFLAGS=['-pthread'],
           FRAMEWORKS=[])

program = env.Program('./sprayFlamelet', ['./sprayFlamelet.cpp', './src/Lagrangian.cpp', './src/SprayStFlow.cpp', './src/Sim1D.cpp', './src/OneDim.cpp','./src/boundaries1D.cpp'])
Default(program)
