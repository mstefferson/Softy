% Function returns parameters used for recording something a distinct time
% intervals and aligns a final run time and recording time
% to fit with a given time step

% Inputs:
% dt: time step size
% t_tot: total run time
% t_rec: record after this much time

% Outputs:
% Nt: Number of time points. Includes zero
% N_rec: Number of recorded time points, including zero
% N_count: Number of time steps to count before recording. 

% Guide: t == timestep
% time    = t*dt
% t_tot   = Nt * dt
% t_rec   = N_count * dt  

function [t_tot,Nt,t_rec,N_rec,N_count] = TimeStepRecMaker(dt,t_tot,t_rec)

% Fix the recording time to be divisible by the time step
if mod(t_rec,dt) ~= 0
  t_rec = floor(t_rec/dt)*dt;
end

% Fix the total run time to be divisible by time step 
if mod(t_tot,dt) ~= 0
  t_tot = floor(t_tot/dt)*dt;
end

% Calculate the outputs
Nt = round(t_tot/dt);           % Number of time steps
N_rec = round(t_tot/t_rec) + 1; % Number of recorded points. +1 includes 0
N_count = round(t_rec/dt);      % Number of time steps before recording

end
