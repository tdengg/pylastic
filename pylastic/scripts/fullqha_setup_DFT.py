from pylastic.thermo.fullqh import Setup
a0=3.14500546899
setup = Setup()
setup.read_POS('POSCAR')
setup.setup_DFT_Cij(0.1,11)
