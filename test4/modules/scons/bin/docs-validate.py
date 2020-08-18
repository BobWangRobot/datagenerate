#!/usr/bin/env python
#
# Searches through the whole source tree and validates all
# documentation files against our own XSD in docs/xsd.
#
from __future__ import print_function

import sys,os
import SConsDoc

if __name__ == "__main__":
    if len(sys.argv)>1:
        if SConsDoc.validate_all_xml((sys.argv[1],)):
            print("OK")
        else:
            print("Validation failed! Please correct the errors above and try again.")
    else:
        if SConsDoc.validate_all_xml(['src',
                                      os.path.join('doc','design'),
                                      os.path.join('doc','developer'),
                                      os.path.join('doc','man'),
                                      os.path.join('doc','python10'),
                                      os.path.join('doc','reference'),
                                      os.path.join('doc','user')
                                      ]):
            print("OK")
        else:
            print("Validation failed! Please correct the errors above and try again.")
            sys.exit(1)
