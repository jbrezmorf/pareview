import unittest
from eclipse_reader import *


class TestEclipseIO(unittest.TestCase):
    def setUp(self):
        self.eio=EclipseIO()

    def test_read_egrid_mesh(self):
        pdo=[]
        self.eio.read_egrid_mesh(pdo,"test.egrid")
        
    def test_skip_after_keyword(self):
        f=open("xx.tmp","wb")
        f.write("jedna567\n9dva345678\ntri3456789"+(4*1024-30-4)*'x'+"milioneBINGO")
        f.close()
        f=open("xx.tmp","rb")
        self.failUnless( -1==self.eio.skip_after_keyword(f,"yy") )
        f.seek(0)
        self.failUnless( 1==self.eio.skip_after_keyword(f,"jedna") )
        self.assertEqual( 5, f.tell())
        self.failUnless( "567\n9"==f.read(5))
        f.seek(0)
        self.failUnless( 1==self.eio.skip_after_keyword(f,"dva3") )
        self.assertEqual( 14, f.tell())
        self.failUnless( "45678"==f.read(5))
        f.seek(0)
        self.failUnless( 1==self.eio.skip_after_keyword(f,"tri") )
        self.assertEqual( 23, f.tell())
        self.failUnless( "3456789xxx"==f.read(10))
        f.seek(0)
        self.failUnless( 1==self.eio.skip_after_keyword(f,"milione") )
        self.assertEqual( 4*1024+3, f.tell())
        self.assertEqual( "BINGO", f.read(5))
        f.close()
        
def main():
    unittest.main()

if __name__ == '__main__':
    main()
