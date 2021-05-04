
import sys,os,unittest
from lammps import lammps

class PythonCommand(unittest.TestCase):

    def setUp(self):
        machine=None
        if 'LAMMPS_MACHINE_NAME' in os.environ:
            machine=os.environ['LAMMPS_MACHINE_NAME']
        self.lmp=lammps(name=machine,
                        cmdargs=['-nocite',
                                 '-log','none',
                                 '-echo','screen',
                                 '-var','zpos','1.5',
                                 '-var','x','2'])
        # create demo input strings and files
        # a few commands to set up a box with a single atom
        self.demo_input="""
region       box block 0 $x 0 2 0 2
create_box 1 box
create_atoms 1 single 1.0 1.0 ${zpos}
"""
        # another command to add an atom and use a continuation line
        self.cont_input="""
create_atoms 1 single &
            0.2 0.1 0.1
"""
        self.demo_file='in.test'
        with open(self.demo_file,'w') as f:
            f.write(self.demo_input)
        self.cont_file='in.cont'
        with open(self.cont_file,'w') as f:
            f.write(self.cont_input)

    # clean up temporary files
    def tearDown(self):
        if os.path.exists(self.demo_file):
            os.remove(self.demo_file)
        if os.path.exists(self.cont_file):
            os.remove(self.cont_file)

    ##############################
    def testFile(self):
        """Test reading commands from a file"""
        natoms = self.lmp.get_natoms()
        self.assertEqual(natoms,0)
        self.lmp.file(self.demo_file)
        natoms = self.lmp.get_natoms()
        self.assertEqual(natoms,1)
        self.lmp.file(self.cont_file)
        natoms = self.lmp.get_natoms()
        self.assertEqual(natoms,2)

    def testNoFile(self):
        """Test (not) reading commands from no file"""
        self.lmp.file(None)
        natoms = self.lmp.get_natoms()
        self.assertEqual(natoms,0)

    def testCommand(self):
        """Test executing individual commands"""
        natoms = self.lmp.get_natoms()
        self.assertEqual(natoms,0)
        cmds = self.demo_input.splitlines()
        for cmd in cmds:
            self.lmp.command(cmd)
        natoms = self.lmp.get_natoms()
        self.assertEqual(natoms,1)

    def testCommandsList(self):
        """Test executing commands from list of strings"""
        natoms = self.lmp.get_natoms()
        self.assertEqual(natoms,0)
        cmds = self.demo_input.splitlines()+self.cont_input.splitlines()
        self.lmp.commands_list(cmds)
        natoms = self.lmp.get_natoms()
        self.assertEqual(natoms,2)

    def testCommandsString(self):
        """Test executing block of commands from string"""
        natoms = self.lmp.get_natoms()
        self.assertEqual(natoms,0)
        self.lmp.commands_string(self.demo_input+self.cont_input)
        natoms = self.lmp.get_natoms()
        self.assertEqual(natoms,2)

    def testNeighborList(self):
        self.lmp.command("units lj")
        self.lmp.command("atom_style atomic")
        self.lmp.command("atom_modify map array")
        self.lmp.command("boundary f f f")
        self.lmp.command("region box block 0 2 0 2 0 2")
        self.lmp.command("create_box 1 box")

        x = [
          1.0, 1.0, 1.0,
          1.0, 1.0, 1.5
        ]

        types = [1, 1]

        self.assertEqual(self.lmp.create_atoms(2, id=None, type=types, x=x), 2)
        nlocal = self.lmp.extract_global("nlocal")
        self.assertEqual(nlocal, 2)

        self.lmp.command("mass 1 1.0")
        self.lmp.command("velocity all create 3.0 87287")
        self.lmp.command("pair_style lj/cut 2.5")
        self.lmp.command("pair_coeff 1 1 1.0 1.0 2.5")
        self.lmp.command("neighbor 0.1 bin")
        self.lmp.command("neigh_modify every 20 delay 0 check no")

        self.lmp.command("run 0")

        self.assertEqual(self.lmp.find_pair_neighlist("lj/cut"), 0)
        nlist = self.lmp.get_neighlist(0)
        self.assertEqual(len(nlist), 2)
        atom_i, numneigh_i, neighbors_i = nlist[0]
        atom_j, numneigh_j, _ = nlist[1]

        self.assertEqual(atom_i, 0)
        self.assertEqual(atom_j, 1)

        self.assertEqual(numneigh_i, 1)
        self.assertEqual(numneigh_j, 0)

        self.assertEqual(1, neighbors_i[0])


##############################
if __name__ == "__main__":
    unittest.main()
