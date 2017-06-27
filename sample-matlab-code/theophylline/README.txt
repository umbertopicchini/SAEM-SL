This set of MATLAB files runs the SAEM-SL simulation for the second example in the paper (section 5.2).

It produces 3 estimation procedures of parameters from three randomly generated datasets, starting at the same parameters location.

The user might want to enlarge the "numattempts" value (in the paper we use numattempts=100, i.e. 100 independent estimation procedures).


Perhaps more interesting plots will be obtained by enlarging the mean value of the starting parameter values,
from bigtheta_start = [0,     1.39,       1.39] to bigtheta_start = [0,     1.61,       1.61]; 

Notice parameters are defined on log-scale.

# files:
- theophylline_run                    the main run file
- theophylline_statemodel             defines the X coordinate of the model
- theophylline_errormodel             defines the coordinate Y of the model, i.e. Y = X + noise
- theophylline_summaries              plug here user defined summaries for X and Y
- saem_synlik  (*)                    the SAEM-SL algorithm  (*)


(*) this is customised for this specific problem and is not the same one found in the "nonlingauss" folder.