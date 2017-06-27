This set of MATLAB files runs the SAEM-SL simulation for the first example in the paper (section 5.1).

It produces 30 estimations of parameters from the same dataset, starting at 30
random parameters location.

The user might want to reduce the value of the "numattempts" from 30 to, say 5-10.


Perhaps more interesting plots will be obtained by enlarging the mean value of the starting parameter values,
from bigtheta_start = [0,     1.39,       1.39] to bigtheta_start = [0,     1.61,       1.61]; 

Notice parameters are defined on log-scale.

# files:
- nonlingauss_run                    the main run file
- nonlingauss_statemodel             defines the X coordinate of the model
- nonlingauss_errormodel             defines the coordinate Y of the model, i.e. Y = X + noise
- nonlingauss_summaries              plug here user defined summaries for X and Y
- saem_synlik (*)                    the SAEM-SL algorithm  (*)


(*) this is customised for this specific problem and is not the same one found in the "theophylline" folder.