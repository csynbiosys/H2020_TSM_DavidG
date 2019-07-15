function [] = runOED_M1_DVID_V1( resultBase, numLoops, numExperiments )

% cd ('../../');
% AMIGO_Startup();
% 
% cd ('\Examples\TSM\model_calibration\M1\OED');

% Selected boundaries for the parameters

theta_min = [0.03*0.1, 10*0.1, 10*0.1, 0, 0, 0,...
    0.1*0.1, 1*0.1, 30*0.1, 0.01, 0, 0,...
    (5e-2)*0.1, (1e-1)*0.1]; 
theta_max = [0.03*10, 10*10, 10*10, 100, 5, 5,...
    0.1*10, 1*10, 30*10, 1, 5, 5,...
    (5e-2)*10, (1e-1)*10];

% Create a matrix of initial guesses for the parameters, having as many
% rows as the number of PE iterations (numExperiments) 
% Each vector is passed as input to the computing function


M = zeros(numExperiments,length(theta_min));
for i=1:numExperiments
   
    p1 = lognrnd(-3.50655789731998,1.15129254649702);
    p2 = lognrnd(2.30258509299405,1.15129254649702);
    p3 = lognrnd(3.40119738166216,1.15129254649702);
    p4 = lognrnd(2.30258509299405,1.15129254649702);
    p5 = normrnd(2.5,1.25);
    p6 = normrnd(2.5,1.25);
    p7 = lognrnd(-2.30258509299405,1.15129254649702);
    p8 = lognrnd(2.22044604925031e-16,1.15129254649702);
    p9 = lognrnd(3.40119738166216,1.15129254649702);
    p10 = lognrnd(-2.30258509299405,1.15129254649702);
    p11 = normrnd(2.5,1.25);
    p12 = normrnd(2.5,1.25);
    p13 = lognrnd(-3.2188758248682,1.15129254649702);
    p14 = lognrnd(-2.30258509299405,1.15129254649702);
    
    M(i,:)=[p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14];
    for j=1:length(theta_min)
        
       if M(i,j)<theta_min(j)
           M(i,j)=theta_min(j);           
       end
       
       if M(i,j)>theta_max(j)
           M(i,j)=theta_max(j);
       end
       
    end
end

ParFull = M; % in this case I am fitting all the values
save('MatrixParameters_InputComparison_M1.mat','ParFull');
save(strcat('MatrixParameters_InputComparison_M1_',resultBase,'.mat'),'ParFull');

load('MatrixParameters_InputComparison_M1.mat');   

parfor epcc_exps=1:numExperiments
        stepd = 180;
        epccNumLoops = numLoops;
        try
            global_theta_guess = ParFull(epcc_exps,:);
            epccOutputResultFileNameBase = [resultBase,'-','OptstepseSS-',num2str(numLoops),'_loops-',num2str(epcc_exps)];
            [out]=OED_M1_DVID(epccOutputResultFileNameBase,epccNumLoops,stepd,epcc_exps,global_theta_guess);

        catch err
            %open file
            errorFile = [resultBase,'-','OptstepseSS-',num2str(numLoops),'_loops-',num2str(epcc_exps),'.errorLog'];
            fid = fopen(errorFile,'a+');
            fprintf(fid, '%s', err.getReport('extended', 'hyperlinks','off'));
            % close file
            fclose(fid);
        end
end