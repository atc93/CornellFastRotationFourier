#=== COMPILATION ===#

Being defined in the 'makefile'

To compile code, run

	make 'name of the code'

or to compile all the codes, just run

	make


#=== CODE ===#


--> FastRotationToyModel

src/FastRotationToyModel.cpp

Toy model that generates a Fast Rotation Signal for a detector located
at the end of the rin (360 deg). The parameters are:

 Momentum spread (1 sigma) in %
     Beam length (1 sigma) in ns
 Number of simulated muons
               Fill length in us

By default only save the Fast Rotation Signal in a histogram in the
ROOT file. To also save the tree, use the 'fillTree' switch. Outputs
are in 'plot' and 'root' directories
