import unittest
from eclipse_reader import *

class Utils(unittest.TestCase):
    def test_skip_after_keyword(self):
        f=open("xx.tmp","wb")
        f.write("")
        f.close()
        f=open("xx.tmp","rb")
        self.failUnless( -1==skip_after_keyword(f,"xx") )
        f.close()
        
def main():
    unittest.main()

if __name__ == '__main__':
    main()
