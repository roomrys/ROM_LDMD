


## lagrangian_dmd.m

### Overview
The basic function of *lagrangian_dmd.m* is to use a FOM with lagrangian data to create a ROM using Lagrangian DMD. The entire *lagrangian_dmd.m* code can be separated into 7 distinct procedures as described below.
> **0. Load Data**: the FOM data is loaded into a data structure *gauges_struct* using the input *sim_number*. The runtime of the FOM is then extracted, saved, and deleted from *gauges_struct*. Then the field names of the data structure, which correspond to the gauge numbers, are extracted and stored in *gfields*.
> 
> **1. Data Analysis**: The Data Anlaysis step utilizes the *mins_AMR.m* function to discover which gauges use an Adaptive Mesh Refinement Level of 2 (where AMR L1 is the coarsest refinement and AMR L2 is the finest refinement; some gauges use both levels of refinement).
> 
> **2. Prep Data**: The Prep Data step utilized the *prep_data.m* function to return the scaled or unscaled FOM observable matrix *g_3Da* along with the *scale* used depending on if the input *isScale* is True or False, respectively. The cut off index, which is calculated using input *train_pct*, is used to shorten the observable *g_3Da* matrix, with the truncated matrix stored as *g_3D*, to only include data up to a certain time step.
> 
> **3. LDMD Offline**: The LDMD Offline step utilized the *lDMD.m* function to compute the eigendecomposition of the linear operator used in DMD.
> 
> **4. LDMD Online**: The LDMD Online step utilizes the *predict_DMD.m* function along with the variables returned from the LDMD Offline step to predict future observables. If *isScaled* is True, then the predicted observables will be unscaled prior to Error Analysis.
> 
> **5. Error Analysis**: The Error Analysis step utilizes the function *error_analysis.m* to compute and create error plots saved to the current working directory + ".\Plots\Error[*err*]\" + file name. Furthermore, if *r* is a 1D array with more than one entry, then the Error Analysis step will create a max error per rank plot, a runtime vs. rank plot, and a speed-up factor of FOM/ROM vs. rank plot.
> 
> **6. Video Creation**: The Video Creation step utilizes the functions *create_comp_video.m* and *create_video.m* to create videos of the FOM and ROM waves plotted together and individually, respectively.
> 
> **7. Program Runtime Analysis**: during each step of the overview, a timer was set to record the runtime of that step and saved to a data structure *runtimes* with fields corresponding to each step. The Program Runtime Analysis creates a histogram of the runtimes of each step. The plot has base-name "program_runtime" and is saved to the current working directory + ".\Plots\Runtime\" + file name.

### Inputs
> **sim_number (integer)**: "simulation number", determines which Full-Order Model data structure to load. All FOM .mat files have a base-name of "matlab_matrix"; the integer *sim_number* is appended to the end of the base-name to determine the full name of the data structure to load ("matlab_matrix"[*sim_number*]).
> 
> **err (integer: 1,2,3)**: determines which type of error will be computed each type of error will take the 2-norm of the difference between the FOM and ROM observables w.r.t. a specific variable and then divide this norm by:
>> *err = 2* -- the 2-norm of the FOM observable w.r.t. the aforementioned variable. 
>> *err = 3* -- the maximum absolute value of the FOM observable w.r.t. the aformentioned variable.
>
> **showFigs (string: 'off', 'on')**: shows figures being saved while code is running if *showFigs* = 'on'.
> 
> **train_pct (double: (0,1])**: specifies the percent of FOM data to use for training, i.e., to create the linear operator used in LDMD
> 
> **r (1D array of integers, or string: 'Compute')**: if a 1D array of integers is given, then *lagrangian_dmd.m* will run LDMD at each rank truncation given in r. If *r* = 'Compute', then the LDMD algorithm will calculate the rank truncation to use based of the Singular Value decay of the FOM data.
> 
>**isScaled (bool)**: boolean variable that will scale FOM observables along each variable between [-1, 1] if isScaled is true, and then unscale the data after LDMD.

### Outputs
> **g_3Da (matrix: ([gauges] * [observables/gauge]) X [time steps])**: the FOM observables matrix.
> 
> **Y_pred (matrix: ([gauges] * [observables/gauge]) X [time steps])**: the ROM observables matrix.
> 
> **Phi (matrix: ([gauges] * [observables/gauge]) X [rank trunc])**: the truncated eigenvectors matrix of the linear operator approximated with DMD.
> 
> **D (matrix: [rank trunc] X [rank trunc])**: the truncated (tall rectangular) eigenvalue matrix of the linear operator approximated with DMD.
> 
> **b (vector: [rank trunc] X 1)**: vector used in LDMD Online step where *b* = psuedoinv(*Phi*) * *g_3Da*(:, 1).

### Functions (separate files)
> #### categorize_AMR2.m (similar to categorize_AMR.m which returns both AMR_1 and AMRb gauges)
> The function *AMR2 = **categorize_AMR2**(gauges_struct)* iterates through each gauge matrix stored in the data structure *gauges_struct* and determines whether the gauge uses AMR L1, L2, or a mixture of both. The times for which each gauge is active is plotted in a histogram with base-name "lifetime_gauges" which is then saved to the current working directory + "\Plots\Gagues Lifetime\" + file name. The initial position of each gauges is also plotted where different colors correspond to the gauge using a certain AMR level. This plot is saved
>>  **Inputs**: 
>>> **gauges_struct (data structure containing matrices)**: data structure containing gauge ID's or gauge numbers as fieldnames and the states of the corresponding gauge in a matrix where the rows represent a time-step and the columns are organized as shown in the **Table 1** below.

|     01    |    02    |           03          |     04     |     05     |       06       | 07  |
|:---------:|:--------:|:---------------------:|:----------:|:----------:|:--------------:|-----|
| AMR Level | Time (s) | Topography Height (m) | x-coor (m) | y-coor (m) | Surface Height | Aux |
**Table 1**: Columns of matrices in *gauges_struct* are organized as shown above.

>> **Outputs**: 
>>> **AMR2.gauge_numbers (data structure)**: the gauge numbers of gauges who only use Adaptive Mesh Refinement (AMR) level 2
>>> **AMR2.dt_final (double)**: the minimum number of time steps where all AMR level 2 gauges are active
>>
>>  **Functions**:
>>>
>>
> #### prep_data.m
> The function *[g_3D, scale] = **prep_data**(gauges_struct, AMR2, isScaled)* extracts the individual matrices for each AMR L2 gauge from the data structure *gauges_struct* using the list of AMR L2 gauges in *AMR2.gauge_numbers*. The individual matrices are concatenated into a 3D matrix which is then sliced to only contain the observables we want (x-coor, y-coor, water height). The 3D matrix *g_3D* is flattened into 2D where all observables for all gauges at a single time step are stacked into a column; thus, the rows denote gauge data and the columns denote time-steps. Then, each observable in the 2D matrix *g_3D* is scaled between [-1, 1] if *isScaled* is true.
>> **Inputs**:
>>> **gauges_struct (data structure containing matrices)**: data structure containing gauge ID's or gauge numbers as fieldnames and the states of the corresponding gauge in a matrix where the rows represent a time-step and the columns are organized as shown in the **Table 1** above.
>>> **AMR2 (data structure)**: data structure containing the data structure *AMR2.gauge_numbers* which has the gauge numbers of gauges that use only AMR L2 and the double *AMR2.dt_final* which is the final time step when all AMR2 gauges are active.
>>> **isScaled (bool)**: boolean variable that will scale FOM observables along each variable between [-1, 1] if isScaled is true, and then unscale the data after LDMD.
>> 
>> **Outputs**: 
>>> **g_3D (matrix: ([gauges] * [observables/gauge]) X [time steps])**: matrix containing observables for the FOM where the rows are organized as shown in **Table 2** and each column is a time step.

|  Row                 | g[gauge #]: [data]      |
|:--------------------:|-------------------------|
|   01                 | g1: x-coor              |
|   ...                |           ...           |
|   num gauges         | g[last]:x-coor          |
|   num gauges + 1     | g1: y-coor              |
|   ...                |           ...           |
|   (num gauges)X2     | g[last]:y-coor          |
|   (num gauges)X2 + 1 | g1: surface height      |
|   ...                |           ...           |
|   (num gauges)X3     | g[last]:surface height  |
**Table 2**: The rows in the *g_3D* matrix are organized as shown above by stacking first the x-coor of all gauges, then the y-coor, and then the surface height.

>>> **scale (double vector: 1 X 3)**: the scale used on each observables (x-coor, y-coor, surface height) to scale the FOM data between [-1, 1].
>>
>>  **Functions**:
>>> **[U_dmd, S_dmd, V_dmd] = sv_analysis(S_dmd)**:
>>> **[U, S, V, K_tilde, W, D, r] = check_rank_trunc(U_dmd, S_dmd, V_dmd, K_tilde, W, D, r, Y2)**:
>>
> #### lDMD.m
> The function *[Phi, D, b, t_svd] = **lDMD**(Y, r)* runs through the Dynamical Mode Decomposition (DMD) Offline steps, i.e., reconstructing the eigendecomposition of the linear operator used in DMD. 
>> **Inputs**:
>>> **Y (double matrix: ([gauges] * [observables/gauge]) X [time steps])**: observables from FOM.
>>> **r (positive integer)**: the rank trunction to use during LDMD.
>> 
>> **Outputs**: 
>>> **Phi (double matrix: ([gauges] * [observables/gauge]) X [rank trunc])**: the truncated eigenvector matrix of the linear operator approximated by LDMD.
>>> **D (double matrix: [rank trunc] X [rank trunc])**: the truncated (tall rectangular) eigenvalue matrix of the linear operator approximated with DMD.
>>> **b (vector: [rank trunc] X 1)**: vector used in LDMD Online step where *b* = psuedoinv(*Phi*) * *Y*(:, 1).
>>> **t_svd (double)**: the time to compute the Singular Value Decomposition (SVD) of Y(:, 1:(end-1)).
>>
>>  **Functions**:
>>>
>>
> #### predict_DMD.m
> The function *Y_pred = **predict_DMD**(g_3D, Phi, D, b, params)* takes care of the online step of DMD (makes the predictions). The predictions are made as follows: Y = Phi * ((D .^ i) * b). This method yields imaginary terms in Y which we ignore in our final answer of Y_pred.
>> **Inputs**:
>>> **g_3D (matrix: ([gauges] * [observables/gauge]) X [time steps])**: matrix containing observables for the FOM where the rows are organized as shown in **Table 2** and each column is a time step.
>>> **Phi (double matrix: ([gauges] * [observables/gauge]) X [rank trunc])**: the truncated eigenvector matrix of the linear operator approximated by LDMD.
>>> **D (double matrix: [rank trunc] X [rank trunc])**: the truncated (tall rectangular) eigenvalue matrix of the linear operator approximated with DMD.
>>> **b (vector: [rank trunc] X 1)**: vector used in LDMD Online step where *b* = psuedoinv(*Phi*) * *Y*(:, 1).
>>> **params (data structure)**: extra parameters that define the problem at hand. In this instance, used for the integer denoting number of predictions wanted *params.num_pred*.
>>
>> **Outputs**: 
>>> **Y_pred (double matrix: ([gauges] * [observables/gauge]) X [num predictions in params.num_pred])**: The real parts of the ROM observables predicted.
>>
> #### error_analysis.m
> The function *function err_struct = **error_analysis**(G, Y_pred, chosen, AMR, gauges_struct, params)* takes care of all the error computation of DMD predicted observables vs. FOM observables.
>> **Inputs**:
>>> **G (matrix: ([gauges] * [observables/gauge]) X [num predictions])**: matrix containing observables for the FOM where the rows are organized as shown in **Table 2** and each column is a time step.
>>> **Y_pred (double matrix: ([gauges] * [observables/gauge]) X [num predictions in params.num_pred])**: The real parts of the ROM observables predicted.
>>> **chosen (string matrix: 1 X [number variables to use in LDMD])**: lists the variables to use in LDMD as a matrix of strings i.e. ["x", "y", "h"].
>>> **AMR (data structure)**: contains information on the adaptive mesh refinement used to obtain the FOM data.
>>> **gauges_struct (data structure)**: separates each observable by variables i.e. gauges_struct.x contains all observables of the x-coor for each gauge.
>>> **params (data structure)**: extra parameters that define the problem at hand. In this instance, used for the integer denoting the cut-off index for the training data stored in *params.cutoff_idx*.
>> 
>> **Outputs**: 
>>> **err_struct (data structure)**: contains information on the error along each variable, i.e. err_struct.x contains all error information for each gauge along the x-coor.
>>
>>  **Functions**:
>>>
>>
> #### create_comp_video.m
> The function ***create_comp_video**(ttl, Y_pred, Y, params, AMR, gauges_struct)* creates video of the ROM and FOM overlaid.
>> **Inputs**:
>>>
>> 
>> **Outputs**: 
>>>
>>
>>  **Functions**:
>>>
>>
> #### create_video.m
> The function ***create_video**(ttl, Y_pred, params, AMR, gauges_struct)* creates video of just one set of observables (ROM or FOM).
>> **Inputs**:
>>>
>> 
>> **Outputs**: 
>>>
>>
>>  **Functions**:
>>>
>>