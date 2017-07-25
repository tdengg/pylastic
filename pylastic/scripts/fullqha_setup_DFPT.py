from pylastic.thermo.fullqh import Setup
setup = Setup()
setup.read_POS('POSCAR')
setup.setup_DFT(0.1,11)
setup.generate_supercells()
