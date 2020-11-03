% ---------------------------------------------------------------------------------------------------------------------------------
%  (1) Please install CVX (http://cvxr.com/cvx/doc/install.html) in Matlab before running the following code.
%  (2) For details of the code, please refer to “Optimal beam pattern design for hybrid beamforming in millimeter wave communications,” IEEE Trans. Veh. Tech., vol. 69, no. 7, pp. 7987-7991, Jul. 2020.
% ---------------------------------------------------------------------------------------------------------------------------------

clc;
clear; 

Nt = 16;
 
N_rf = 5; %% RF chain number is critical to Benchmark 1

SteeringVector=@(theta) 1/sqrt(Nt)*exp(j*pi*(0:1:Nt-1)'*theta);
 

AoD_range = [-0.5  0.5]; 

 
%% Proposed Method 

%%% Basis Search and Dimensionality Reduction 
epsilon =0.99999;
[Fmatrix BasisNum]= BasisSearch(AoD_range, Nt, epsilon);
 
%%% Passband Flattening
TransitionBand  = 1/Nt;
AngleStep = 0.002;
[x_vector SDR_PerIndicator] = PassbandFlattening(AoD_range, Fmatrix, Nt, AngleStep, TransitionBand);

%%% Generate Reference Beamforming Vector 
BeamformingVectorRef = Fmatrix*x_vector;
BeamformingVectorRef =  BeamformingVectorRef/norm(BeamformingVectorRef,2);


%%% Hybrid Decomposition,  M_matrix -- Analog beamforming matrix, w_vector -- Digital beamforming vector
[M_matrix w_vector] = HybridStructuralization(BeamformingVectorRef);
BeamformingVector = M_matrix*w_vector;
BeamformingVector = BeamformingVector/norm(BeamformingVector,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 

%% Benchmark 1: Main idea is from "Channel Estimation and  Hybrid Precoding for Millimeter Wave Cellular Systems"
%%% f_approx_rect is the beamforming vector that can be implemented with N_rf RF chains 
%%% f_ref_rect is the reference beamforming vector for f_approx_rect
[f_ref_rect f_approx_rect M_matrix  w_vector] = LSOMP(AoD_range, N_rf, Nt, AngleStep, TransitionBand);
 



%% Benchmark 2: “Hierarchical codebook design for beamforming training in millimeter-wave communication,”
load HierCodebook.mat  
% Hierarchical codebook structure is  codebook{1}-->[16x1 double]     codebook{2}-->[16x2 double]     codebook{3}-->[16x4 double]     codebook{4}-->[16x8 double]     codebook{5}-->[16x16 double]

% f_subarray = codebook{3}(:,3).*SteeringVector(-0.25); %% AoD_range = [-0.25  0.25]
f_subarray = codebook{2}(:,2).*SteeringVector(-0.5); %% AoD_range = [-0.5  0.5]

f_subarray = f_subarray*diag(1./sqrt(diag(f_subarray'*f_subarray)));
 


%% Beam Pattern Plot, Calculation of Ripple Factor & Energy Ratio
cnt = 0;
amp_ref_rect_pass = [];
amp_approx_rect_pass = [];
amp_subarray_pass = [];
amp_BeamformingVector_pass = [];

for iitheta = -1:0.001:1
    cnt = cnt + 1;
    amp_approx_rect(cnt) =  norm(SteeringVector(iitheta)'*f_approx_rect,2);
    amp_subarray(cnt) =  norm(SteeringVector(iitheta)'*f_subarray,2);
	amp_BeamformingVector(cnt) =  norm(SteeringVector(iitheta)'*BeamformingVector,2);
    
	if  iitheta> min(AoD_range)   && iitheta < max(AoD_range) 
       amp_approx_rect_pass = [amp_approx_rect_pass amp_approx_rect(cnt) ];
       amp_subarray_pass = [amp_subarray_pass amp_subarray(cnt)];
       amp_BeamformingVector_pass = [amp_BeamformingVector_pass amp_BeamformingVector(cnt)];
	end
   
end

%%% Ripple Factor 
RF1 = sqrt(mean((amp_approx_rect_pass-mean(amp_approx_rect_pass)).^2))/mean(amp_approx_rect_pass);
RF2 = sqrt(mean((amp_subarray_pass-mean(amp_subarray_pass)).^2))/mean(amp_subarray_pass);
RF3 = sqrt(mean((amp_BeamformingVector_pass-mean(amp_BeamformingVector_pass)).^2))/mean(amp_BeamformingVector_pass);

display(['Ripple Factor of Benchmark 1 is ',  num2str(RF1),])
display(['Ripple Factor of Benchmark 2 is ',  num2str(RF2)])
display(['Ripple Factor of Proposed Method is ', num2str(RF3)])
display(' ')

%%% Energy Ratio
ER1 = sum(amp_approx_rect_pass.^2)/sum(amp_approx_rect.^2);
ER2 = sum(amp_subarray_pass.^2)/sum(amp_subarray.^2);
ER3 = sum(amp_BeamformingVector_pass.^2)/sum(amp_BeamformingVector.^2);

display(['Energy Ratio of Benchmark 1 is ',  num2str(ER1),])
display(['Energy Ratio of Benchmark 2 is ',  num2str(ER2)])
display(['Energy Ratio of Proposed Method is ', num2str(ER3)])




%%% Beam Pattern
figure(1)
plot(-1:0.001:1, amp_approx_rect,'r')
hold on 
plot(-1:0.001:1, amp_subarray )
plot(-1:0.001:1, amp_BeamformingVector,'g')
hold off
title('Beam Pattern (where the angle is Cosine AoD)')
legend('Benchmark 1 with N_{rf}= 5', 'Benchmark 2', 'Proposed method')
 
figure(2)
subplot(2,2,1)
polarplot(acos(-1:0.001:1),amp_approx_rect,'r')
thetalim([0 180])
title('Benchmark 1 with N_{rf}= 5')
subplot(2,2,2)
polarplot(acos(-1:0.001:1),amp_subarray)
thetalim([0 180])
title('Benchmark 2')
subplot(2,2,3)
polarplot(acos(-1:0.001:1),amp_BeamformingVector,'g')
thetalim([0 180])
title('Proposed method')
hold off
suptitle('Beam Pattern (where the angle is AoD)')