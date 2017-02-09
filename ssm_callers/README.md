## SSM callers

The Matlab (wrapper) functions to invoke the simulation of SSM's.
SSM can be any model producing structural data set to fit to the QSM's data sets.
Moreover, SSM can be implemented in any environment as long as this can be reached
from within Matlab (as, for example, the case with LPFG/C++ simulator based SSM's).

Currently we employ only LPFG-based simulators. We have also developed a simple
meta-language interface to define the macro-like (#define's) parameters 
for the simulators. See the description of this meta-languate in the wrapper functions.
