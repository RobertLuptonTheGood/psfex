# -*- python -*-
import os
from lsst.sconsUtils import scripts
env = scripts.BasicSConstruct("meas_extensions_psfex")

if os.environ.has_key("PSFEX_DIR"):
    env.Append(CPPPATH = [os.environ["PSFEX_DIR"]]) # needed for config.h.  N.b. doesn't work with _wrap.cc
#env.Append(CCFLAGS = ['-DHAVE_CONFIG_H=1'])
