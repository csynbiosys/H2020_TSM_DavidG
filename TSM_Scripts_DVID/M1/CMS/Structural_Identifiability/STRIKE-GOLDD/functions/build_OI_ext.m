%--------------------------------------------------------------------------
% Function that builds the generalized observability-identifiability matrix
% of the model specified in the 'options' script. It takes as input the
% number of Lie derivatives, provided by the user. The resulting array is
% stored in a MAT file.
%--------------------------------------------------------------------------

function build_OI_ext(num_der)

%==========================================================================
% Read options, add folders to path:
global x f p maxstates u %%% new
[modelname,paths,opts,submodels,prev_ident_pars] = options;
maxstates = opts.maxstates;
addpath(genpath(paths.models));
addpath(genpath(paths.results));
addpath(genpath(paths.functions));

%==========================================================================
% Load model:
load(modelname)
fprintf('\n Analyzing identifiability of %s ... \n', modelname);
    
% ==========================================================================
% Remove parameters that have already been classified as identifiable:
if exist('prev_ident_pars','var') == 1
    for np=1:numel(prev_ident_pars)
        [~, original_index] = ismember(prev_ident_pars(np),p); 
        p(original_index)=[];                
    end
end

%==========================================================================
% Dimensions of the problem: 
m = numel(h);                     % number of outputs
n = numel(x);                     % number of states
q = numel(p);                     % number of unknown parameters
xaug = [x;p];                     % states and parameters 
faug = [f; zeros(numel(p),1)];    % augmented dynamics

fprintf('\n >>> The model contains:\n %d states:\n %s',numel(x),char(x));
fprintf('\n %d outputs:\n %s',numel(h),char(h));
fprintf('\n %d inputs:\n %s',numel(u),char(u));%%% new
fprintf('\n %d parameters:\n %s',q,char(p));

%==========================================================================
% Build Identifiability-Observability matrix, Oi(x)
fprintf('\n >>> Building Identifiability-Observability matrix: ');

% 1st block:
tic
onx = jacobian(h,xaug); 
totaltime = toc;
ind = 2;

%--- Input derivatives ------------------------------------------------
max_num_der = max(opts.nnzDerIn);     
for num_in=1:numel(u) % create array of derivatives of the inputs
	input_der(num_in,:) = [u(num_in),sym(strcat(char(u(num_in)),sprintf('_d')),[1 max_num_der])]; 
end                   
% replace zero derivatives of the inputs:
input_der_prov = [input_der,zeros(numel(u),1)];
for ind_u=1:numel(u)
	input_der_prov(ind_u,(opts.nnzDerIn(ind_u)+2):(max_num_der+1))=0;
end
input_der = input_der_prov(:,1:end-1);
syms zero_input_der_dummy_name
%----------------------------------------------------------------------
past_Lie = h;         %%% new
extra_term = 0;       %%% new
    
% Lie derivatives:
while ind < (num_der+2) 
    tic
    newLie = onx(((ind-2)*m+1):(ind-1)*m,:)*faug;
    %--- add the extra term of the extended Lie derivative:
    if ind>2 %%% new
        if numel(u) > 0
            extra_term = 0; % reset for each new Lie derivative
            for i=1:ind-2 %%% new
				if i < size(input_der,2) 
					lo_u_der   = input_der(:,i);   
					hi_u_der   = input_der(:,i+1);
					lo_u_der   = subs(lo_u_der,0,zero_input_der_dummy_name);
					extra_term = extra_term + jacobian(past_Lie,lo_u_der)*hi_u_der;
				end 
            end%%% new
        end
    end %%% new 
    ext_Lie = newLie + extra_term;%%% new
    past_Lie = ext_Lie; %%% new       
    %---
    newOnx = jacobian(ext_Lie,xaug);
    onx = [onx; newOnx];   
    lasttime = toc;
    totaltime = totaltime + lasttime;
    clear newLie newOnx
    ind = ind+1;
    obsidentmatrix = sprintf('obs_ident_matrix_%s_%d_Lie_deriv',modelname,ind-2);
    obsidentfile = strcat(pwd,filesep,'results',filesep,obsidentmatrix);
    save(obsidentfile);
    fprintf('\n >>> Observability-Identifiability matrix built with %d Lie derivatives',ind-2);
    fprintf('\n     (calculated in %d seconds)\n',totaltime);
end
end    