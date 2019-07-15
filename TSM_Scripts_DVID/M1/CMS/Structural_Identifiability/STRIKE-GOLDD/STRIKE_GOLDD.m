%==========================================================================
% StrIkE-GOLDD (Structural Identifiability taken as Extended-Generalized
% Observability with Lie Derivatives and Decomposition)
%--------------------------------------------------------------------------
%
% A Matlab toolbox for structural identifiability analysis of nonlinear models
%
%--------------------------------------------------------------------------
%
% User input can be specified in file 'options.m'
%
%--------------------------------------------------------------------------
% Version v2.0.3 
% Last modified: 26/04/2018
% Alejandro Fernandez Villaverde (afvillaverde@iim.csic.es)
%==========================================================================

function STRIKE_GOLDD(varargin)

fprintf('\n\n ------------------------- \n');
fprintf(' >>> StrIkE-GOLDD TOOLBOX \n');
fprintf(' ------------------------- \n');

%==========================================================================
% Read options and add folders to path:
tStart = tic;
global x f p maxstates u unidflag 
if nargin > 0
    [modelname,paths,opts,submodels,prev_ident_pars] = run(varargin{1});
else
    [modelname,paths,opts,submodels,prev_ident_pars] = options;
end
maxstates = opts.maxstates;
addpath(genpath(paths.meigo));
addpath(genpath(paths.models));
addpath(genpath(paths.results));
addpath(genpath(paths.functions));

%==========================================================================
% Load model:
load(modelname)
fprintf('\n Analyzing identifiability of %s ... \n', modelname);

tic
%==========================================================================
% Initialize variables:
identifiables  = [];       % identifiable parameters.
nonidentif     = [];       % unidentifiable parameters.
obs_states     = [];       % observable states.
unobs_states   = [];       % unobservable states.
lastrank       = NaN;
decomp_flag    = 0;
unidflag       = 0;
skip_elim      = 0;
    
% ==========================================================================
% Remove parameters that have already been classified as identifiable:
if exist('prev_ident_pars','var') == 1
    for np=1:numel(prev_ident_pars)
        [~, original_index] = ismember(prev_ident_pars(np),p); 
        p(original_index)=[];                
    end
end
    
% ==========================================================================
% Remove states that have already been classified as observable:
if exist('prev_obs_states','var') == 1
    for nx=1:numel(prev_obs_states)
        [~, original_index] = ismember(prev_obs_states(nx),x); 
        x(original_index)=[];          
        f(original_index)=[];  
    end
end

%==========================================================================
% Dimensions of the problem: 
m = numel(h);                     % number of outputs
n = numel(x);                     % number of states
q = numel(p);                     % number of unknown parameters
r = n+q;                          % number of states + parameters
xaug = [x;p];                     % states and parameters 
faug = [f; zeros(numel(p),1)];    % augmented dynamics

% fprintf('\n >>> The model contains:\n %d state(s):\n %s',numel(x),char(x));
% fprintf('\n %d output(s):\n %s',numel(h),char(h));
% fprintf('\n %d input(s):\n %s',numel(u),char(u));
% fprintf('\n %d parameter(s):\n %s',q,char(p));
fprintf('\n >>> The model contains:\n%d state(s):\n',numel(x));
disp(x)
fprintf('%d output(s):\n ',numel(h));
disp(h)
fprintf('%d input(s):\n ',numel(u));
disp(u)
fprintf('%d parameter(s):\n',q);
disp(p)

%==========================================================================
% Check which states are directly measured, if any: 
saidas     = cell(m,1);
ismeasured = zeros(1,n);
for i=1:m
    saidas{i} = char(h(i));
end
for i=1:n
    ismeasured(i) = sum(strcmpi(char(x(i)),saidas));
end
meas_x_indices = find(ismeasured);       % indices of the measured states
unmeas_x_indices = find((1-ismeasured)); % indices of the unmeasured states
meas_x = x(meas_x_indices);              % names of the measured states

%==========================================================================
if opts.forcedecomp == 0 

    %==========================================================================
    % Build Identifiability-Observability matrix, Oi(x)
    nd = ceil(r/m-1);  % minimum number of Lie derivatives
    increaseLie = 1;
    fprintf('>>> Building the observability-identifiability matrix requires at least %d Lie derivatives',nd);
    fprintf('\n     Calculating derivatives: ');
    tic
    onx = zeros(m*(1+nd),r);
    onx = sym(onx);
    onx(1:m,:) = jacobian(h,xaug); % 1st block
    totaltime = toc;
    ind = 2;
    lasttime = 0;
    %--- Input derivatives ------------------------------------------------
    max_num_der = max(opts.nnzDerIn);
    if numel(u)>0
        for num_in=1:numel(u) % create array of derivatives of the inputs
            input_der(num_in,:) = [u(num_in),sym(strcat(char(u(num_in)),sprintf('_d')),[1 max_num_der])]; 
        end                   
        % replace zero derivatives of the inputs:
        input_der_prov = [input_der,zeros(numel(u),1)];
        for ind_u=1:numel(u)
            input_der_prov(ind_u,(opts.nnzDerIn(ind_u)+2):(max_num_der+1))=0;
        end
        input_der = input_der_prov(:,1:end-1);
    else
        input_der = [];
    end
    syms zero_input_der_dummy_name
    %----------------------------------------------------------------------
    past_Lie = h;        
    extra_term = 0;       
    while ind < (nd+2) && lasttime < opts.maxLietime % 2nd and subsequent blocks
        tic
        Lieh = onx(((ind-2)*m+1):(ind-1)*m,:)*faug;
        if ind>2 
            if numel(u) > 0
                extra_term = 0; % reset for each new Lie derivative
                for i=1:ind-2 
                    if i < size(input_der,2) 
                        lo_u_der   = input_der(:,i);   
                        hi_u_der   = input_der(:,i+1);
                        lo_u_der   = subs(lo_u_der,0,zero_input_der_dummy_name);
                        extra_term = extra_term + jacobian(past_Lie,lo_u_der)*hi_u_der;
                    end                        
                end
            end
        end       
        ext_Lie = Lieh + extra_term;
        past_Lie = ext_Lie;         
        onx(((ind-1)*m+1):(ind*m),:) = jacobian(ext_Lie,xaug); 
        lasttime = toc;
        totaltime = totaltime + lasttime;
        ind = ind+1;
        fprintf('%d ',ind-2); 
    end
    if ind == (nd+2)
        obsidentmatrix = sprintf('obs_ident_matrix_%s_%d_Lie_deriv',modelname,nd);
        obsidentfile = strcat(pwd,filesep,'results',filesep,obsidentmatrix);
        save(obsidentfile);
        while increaseLie == 1
            fprintf('\n >>> Observability-Identifiability matrix built with %d Lie derivatives',nd);
            fprintf('\n     (calculated in %d seconds)',totaltime);
            
            %==========================================================================
            % Check identifiability by calculating rank:
            fprintf('\n >>> Calculating rank...');
            tic
            if opts.replaceICs == 1
                xind = find(known_ics);
                if size(ics) ~= size(x)
                    ics = transpose(ics);
                end
                onx = subs(onx,x(xind),ics(xind));
            end  
            if opts.numeric == 1
                allvariables = symvar(onx);
                numeros = vpa(0.1+rand(size(allvariables)));
                numonx = subs(onx,allvariables,numeros); 
            else
                numonx = onx;
            end                 
            rango = double(rank(numonx));
            fprintf('\n     Rank = %d (calculated in %d seconds)',rango,toc);         
            if rango == r 
                obs_states = x;
                identifiables = p;
                increaseLie = 0;
            else% With that number of Lie derivatives the array is not full rank.                       
                % If possible, calculate one more Lie derivative and retry:
                if nd < r  && lasttime < opts.maxLietime && rango ~= lastrank
                    tic
                    nd = nd + 1;
                    ind = 1+nd;   
                    if numel(u) > 0
                        extra_term = 0; % reset for each new Lie derivative
                        for i=1:ind-2 
                            if i < size(input_der,2) 
                                lo_u_der   = input_der(:,i);   
                                hi_u_der   = input_der(:,i+1);
                                lo_u_der   = subs(lo_u_der,0,zero_input_der_dummy_name);
                                extra_term = extra_term + jacobian(past_Lie,lo_u_der)*hi_u_der;
                            end 
                        end
                    end
                    newLie   = onx(((ind-2)*m+1):(ind-1)*m,:)*faug; 
                    ext_Lie  = newLie + extra_term;
                    past_Lie = ext_Lie;        
                    newOnx   = jacobian(ext_Lie,xaug); 
                    onx      = [onx; newOnx];   
                    clear newLie newOnx
                    lasttime  = toc;
                    totaltime = totaltime + lasttime;
                    obsidentmatrix = sprintf('obs_ident_matrix_%s_%d_Lie_deriv',modelname,nd);
                    obsidentfile = strcat(pwd,filesep,'results',filesep,obsidentmatrix);
                    save(obsidentfile);
                    lastrank = rango;
                % If that is not possible, there are several possible causes:
                else 
                    if nd >= r % The maximum number of Lie derivatives has been reached
                        unidflag = 1; 
                        fprintf('\n    The model is structurally unidentifiable as a whole');
                    else
                        if rango == lastrank
                            numonx = numonx(1:(end-m),:);
                            nd = nd - 1;
                            unidflag = 1; % note that the pars may still be identifiable (rank deficiency may be due to initial conditions)
                        else
                            if lasttime >= opts.maxLietime
                                fprintf('\n    => More Lie derivatives would be needed to see if the model is structurally unidentifiable as a whole.');
                                fprintf('\n    However, the maximum computation time allowed for calculating each of them has been reached.');
                                fprintf('\n    You can increase it by changing <<opts.maxLietime>> (currently opts.maxLietime = %d)',opts.maxLietime);
                                unidflag = 0; 
                                if opts.unidentif == 1
                                    skip_elim = 1;
                                    increaseLie = 0;
                                end
                            end
                        end
                    end  
                    if skip_elim == 0 
                        % Eliminate columns one by one to check identifiability of the associated parameters: 
                        [identifiables,nonidentif,obs_states,unobs_states] = ...
                             elim_and_recalc(unmeas_x_indices,rango,numonx,opts);
                         increaseLie = 0;
                    end
                end
            end
        end        
    else% If the maxLietime has been reached, but the minimum of Lie derivatives has not been calculated:
        if opts.unidentif == 0
            if opts.decomp == 0
                fprintf('\n >>> Calculating rank...');
                tic
                if opts.replaceICs == 1
                    xind = find(known_ics);
                    if size(ics) ~= size(x)
                        ics = transpose(ics);
                    end
                    onx = subs(onx,x(xind),ics(xind));
                end  
                if opts.numeric == 1
                    allvariables = symvar(onx);
                    numeros = vpa(0.1+rand(size(allvariables)));
                    numonx = subs(onx,allvariables,numeros); 
                else
                    numonx = onx;
                end
                rango = double(rank(numonx));
                fprintf('\n     Rank = %d (calculated in %d seconds)',rango,toc); 
                [identifiables,nonidentif,obs_states,unobs_states] = ...
                elim_and_recalc(unmeas_x_indices,rango,numonx,opts);
            else
                decomp_flag = 1; 
            end
        end
    end
end % end of 'if opts.forcedecomp == 0'

if decomp_flag == 1 || opts.forcedecomp == 1 
    fprintf('\n The model will be decomposed. \n');
    if opts.decomp_user == 1
        model_results_folder = sprintf('results/decomp_user_%s_%s_maxLietime_%d',modelname,date,opts.maxLietime);
        if exist('submodels','var') == 0 || numel(submodels) == 0
            warning('You have selected to enter the submodels manually. However, no submodels have been specified => STRIKE-GOLDD will decompose the model using optimization. If you want to avoid it, create a non-empty submodels vector in the options file.') 
            opts.decomp_user = 0;
        end
    end
    if opts.decomp_user == 0
        model_results_folder = sprintf('results/decomp_%s_%s_maxstates_%d_maxLietime_%d',modelname,date,maxstates,opts.maxLietime);
        submodels = [];
    end
    mkdir(model_results_folder);
    [identifiables,nonidentif,obs_states,unobs_states] = decomp(modelname,opts,model_results_folder,submodels);
end

%==========================================================================
% Build the vectors of identifiable / non identifiable parameters and
% initial conditions (may be repeated if decomposition was used):
p_id         = symvar(identifiables);        
p_un         = symvar(nonidentif);          
obs_states   = symvar(obs_states);       
unobs_states = symvar(unobs_states); 

%==========================================================================
% Report results:
fprintf('\n\n ------------------------ \n');
fprintf(' >>> RESULTS SUMMARY:\n');
fprintf(' ------------------------ \n');
load(modelname,'p','x')

if numel(p_id) == numel(p)
    fprintf('\n >>> The model is structurally identifiable:');
    fprintf('\n     All its parameters are structurally identifiable. \n');
else    
    if unidflag == 1 
        fprintf('\n >>> The model is structurally unidentifiable.');
        fprintf('\n >>> These parameters are identifiable:\n      %s ',char(p_id));
        fprintf('\n >>> These parameters are unidentifiable:\n      %s \n',char(p_un));
    else             
        fprintf('\n >>> These parameters are identifiable:\n      %s ',char(p_id));
    end  
    %==========================================================================
    % Search for identifiable parameter combinations (combos):
    if opts.findcombos == 1 && exist('onx','var') == 1
        % Save results first, just in case the user kills the process:
        resultsname = sprintf('id_results_%s',modelname);
        fullresultsname = strcat(pwd,filesep,'results',filesep,resultsname);
        save(fullresultsname)      
        [parpde,stringpde] = combos(p_un,onx,n);
        if numel(parpde) == 0
            fprintf('\n\n >>> No identifiable combinations of parameters could be found');
        else
            fprintf('\n\n >>> There are identifiable combinations of parameters.');
            fprintf('\n     They can be found by solving the following PDE(s):');
            for numpdes = 1:numel(stringpde)
                fprintf('\n     %s \n',char(stringpde));
            end
        end       
    end
end
    
if(numel(obs_states))==numel(x)
    fprintf('\n >>> The model is observable:');
    fprintf('\n     All its states are observable. \n');
else
    if numel(obs_states)>0,    fprintf('\n >>> These states are observable (and their initial conditions, if considered unknown, are identifiable):\n      %s ',char(obs_states)); end
    if numel(unobs_states)>0,  fprintf('\n >>> These states are unobservable (and their initial conditions, if considered unknown, are unidentifiable):\n      %s ',char(unobs_states)); end
    if numel(meas_x)>0,        fprintf('\n >>> These states are directly measured:\n      %s ',char(meas_x)); end
end

%==========================================================================
totaltime = toc(tStart);
fprintf('\n Total execution time: %d \n\n',totaltime);

%==========================================================================
% Save results:
if decomp_flag == 0
    resultsname = sprintf('id_results_%s_%s',modelname,date);
else
    resultsname = sprintf('id_results_%s_maxstates_%d_%s',modelname,maxstates,date);
    delete VNS_report.mat
end    
fullresultsname = strcat(pwd,filesep,'results',filesep,resultsname);
save(fullresultsname);

end    
