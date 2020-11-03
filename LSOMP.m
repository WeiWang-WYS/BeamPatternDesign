function [f_ref f_approx M_matrix  w_vector] = LSOMP(AoD_range, N_rf, Nt, AngleStep, TransitionBand)

SteeringVector=@(theta) 1/sqrt(Nt)*exp(theta*j*pi*(0:1:Nt-1)');
MatrixPass = [];
G_vector =[];


% Least Square (LS) method to derive reference beamforming vector
for itheta = -1:AngleStep:1
    if itheta<min(AoD_range)+TransitionBand||itheta>max(AoD_range)-TransitionBand
        G_vector = [G_vector' 0]';
    else
        G_vector = [G_vector' 1]';
    end
    MatrixPass = [MatrixPass SteeringVector(itheta)];
end

f_ref =  inv(MatrixPass*MatrixPass')*MatrixPass*G_vector;
f_ref = f_ref/norm(f_ref,2);

%%%%%% OMP to approximate the reference beamforming vector f_ref 
f_temp = f_ref;
M_matrix=[];

for iRF = 1:1:N_rf
%%%% Optimization  of vector under constant modulus constraint, refer to "Spatial Modulation for Uplink Multi-User mmWave MIMO Systems With Hybrid Structure"
   g_current = exp(j*2*pi*rand(Nt,1));
   g_next = exp(j*2*pi*rand(Nt,1));
   Gmatrix = f_temp*f_temp';
   while norm(g_current-g_next,2)>0.00001
        g_current = g_next;
        g = Gmatrix*g_current;
        g_next = exp(j*phase(g));        
   end
   M_matrix(:,iRF) = g_next;
   f_temp = f_temp - M_matrix*inv(M_matrix'*M_matrix)*M_matrix'*f_temp;
end


w_vector = inv(M_matrix'*M_matrix)*M_matrix'*f_ref;

f_approx = M_matrix*w_vector;

f_approx = f_approx/norm(f_approx,2);
