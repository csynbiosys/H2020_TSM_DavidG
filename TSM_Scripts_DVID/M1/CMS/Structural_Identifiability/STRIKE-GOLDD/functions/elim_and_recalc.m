%--------------------------------------------------------------------------
% Function that determines identifiability of individual parameters one 
% by one, by successive elimination of its column in the identifiability 
% matrix and recalculation of its rank
%--------------------------------------------------------------------------
function [new_ident_pars,new_nonid_pars,new_obs_states,new_unobs_states] = elim_and_recalc(unmeas_xred_indices,rangoinicial,numonx,opts,varargin)

global p x unidflag

if nargin == 4    
    pred = p;
    xred = x; 
    q = numel(pred);
    n = numel(xred);
    qreal = q;
    identifiables = [];
    obs_states = [];
else
    % In case of decomposition
    % varargin: qreal,pred,xred,identifiables,obs_states
    qreal         = varargin{1};
    pred          = varargin{2};
    xred          = varargin{3};
    identifiables = varargin{4};
    obs_states    = varargin{5};
    q = numel(pred);
    n = numel(xred);
end

r  = n+q;
new_ident_pars = [];
new_nonid_pars = [];
new_obs_states  = [];
new_unobs_states  = [];
    

%==========================================================================
% ELIMINATE A PARAMETER:
%==========================================================================
% At each iteration we try removing a different column from onx:
for ind=1:qreal % only the first 'qreal' elements of pred are parameters; the following 'q-qreal' are states considered as parameters--we are not interested in their identifiability 
    isidentifiable = ismember(pred(ind),identifiables);
    if isidentifiable
        fprintf('\n Parameter %s has already been classified as identifiable.',char(pred(ind)))
    else 
        indices = 1:r;
        indices(n+ind) = [];
        num_rank = rank(numonx(:,indices));
        if num_rank == rangoinicial
            if (opts.unidentif == 1) || (opts.forcedecomp == 0 && opts.decomp == 0 && unidflag == 1) %%%
            	fprintf('\n    => Parameter %s is structurally unidentifiable',char(pred(ind)));
            	new_nonid_pars = [new_nonid_pars; pred(ind)];
                unidflag = 1;
            else
            	fprintf('\n    => We cannot decide about identifiability of parameter %s at the moment',char(pred(ind)));
            end  
        else
            if opts.unidentif == 0
            	fprintf('\n    => Parameter %s is structurally identifiable',char(pred(ind)));
            	new_ident_pars = [new_ident_pars; pred(ind)];  
            else
            	fprintf('\n    => We cannot decide about unidentifiability of parameter %s at the moment',char(pred(ind)));
            end      
        end
    end
end

%==========================================================================
% ELIMINATE A STATE:
%==========================================================================
% At each iteration we try removing a different state from 'xred':
if opts.checkObser == 1
    for ind=1:numel(unmeas_xred_indices) % for each unmeasured state
        original_index = unmeas_xred_indices(ind); % in this script, 'original_index' refers to xred
        isobservable = ismember(xred(original_index),obs_states);
        if isobservable
            fprintf('\n State %s has already been classified as observable.',char(xred(original_index)))
        else              
            indices = 1:r;
            indices(original_index) = []; %indices(original_index==unmeas_xred_indices) = [];
            num_rank = rank(numonx(:,indices));
            if num_rank == rangoinicial
                if (opts.unidentif == 1) || (opts.forcedecomp == 0 && opts.decomp == 0 && unidflag == 1) %%%
                    fprintf('\n    => State %s is unobservable',char(xred(original_index)));
                    new_unobs_states = [new_unobs_states; xred(original_index)];
                else
                    fprintf('\n    => We cannot decide about state %s at the moment',char(xred(original_index)));
                end 
            else
                if opts.unidentif == 0
                    fprintf('\n    => State %s is observable',char(xred(original_index)));
                    new_obs_states = [new_obs_states; xred(original_index)];  
                else
                    fprintf('\n    => We cannot decide about state %s at the moment',char(xred(original_index)));
                end   
            end
        end
    end
end

end
