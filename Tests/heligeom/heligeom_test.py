import unittest

from ptools import *

#deactivate auto-diff for now:
def surreal(r):
    return r


class TestHeligeom(unittest.TestCase):
    def setUp(self):
        unittest.TestCase.setUp(self)
        
    def tearDown(self):
        pass
    
    def testGenerateBDNA(self):
        dna = DNA ("bp.ato.pdb","AAAAAATCGATCTATC",ADNA()) #tout atom
        result= dna.PrintPDB()
        expected = open("generate_B_DNA.expected").read()
        self.assertEqual(result, expected)

        
    def testBasicManipulation(self):
        #translate the DNA in coarse grain representation
        dna = DNA ("bp.red.pdb","generate_B_DNA.expected")

        #add a base Pair
        bp = BasePair (dna[0].GetRigidBody())
        dna.Add(bp)

        #add itself
        new = DNA (dna)
        dna.Add(new,BDNA())

        #change the type of a base
        dna.ChangeType (0,"A","bp.red.pdb")

        #turn the center base 
        dna.ApplyLocal (Roll(30), dna.Size()/2)

        #trim the extremities
        dna = dna.SubDNA (2, dna.Size()-3)

        #change to a All Atom representation
        dna.ChangeRepresentation ("bp.ato.pdb")

        result = dna.PrintPDB() 
        expected = open("basicManipulation.expected").read()
        self.assertEqual(result, expected)
        

    def testCGfromPDBFile(self):
        dna = DNA ("bp.ato.pdb","generate_B_DNA.expected") #gros grain

        result = dna.PrintPDB()
        expected = open("CGfromPDBFile.expected").read()
        self.assertEqual(result, expected)
        
    
    def testEnergy(self):
        prot = AttractRigidbody ("1A74_prot.red")
        dna = AttractRigidbody ("1A74_opt.red")

        cutoff=1000
  
        forcefield = AttractForceField1("aminon.par",surreal(cutoff))
        ener = forcefield.nonbon8(prot,dna, AttractPairList (prot, dna,cutoff))
        
        self.assertAlmostEqual(ener, -51.6955215854)
        self.assertEqual(prot.Size(), 706)
        self.assertEqual(dna.Size(), 231)

    
    
    
    
if __name__ == "__main__":    
    unittest.main()