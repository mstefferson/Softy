% HR2DrotDenEvolverFTBodyIdC.M
%
% Description:
% Code is the main body  for the code that
% uses discrete Fourier transforms to solve the diffusion equation of rods in
% 2 spatial directions. These Rods are allowed to diffuse rotationally.
%
% Includes hard rod interactions, assuming in the interaction that the rods
% are infinitely thin.
%
%
%
% Density matrix is set up as rho(x,y,phi)-----> RHO(kx,ky,km)
%
% Isotropic diffusion. The propagaotr is a cube.
%


function [DenRecObj]  = ...
    HR2DrotDenEvolverFTBodyIdC(wfid,lfid,rho,ParamObj,...
    TimeObj,GridObj,DiffMobObj,feq)

fprintf(lfid,'In body of code\n');
% Create a text file that tells user what percent of the program has
% finished

PrgmTrackNm = sprintf('PrgmTracker_%i.txt',ParamObj.trial);
tfid        = fopen(PrgmTrackNm,'w');

%Set N since it used so frequently
Nx  = ParamObj.Nx;
Ny  = ParamObj.Ny;
Nm  = ParamObj.Nm;
N3 = Nx*Ny*Nm;

% FT initial density and max density
TotalDensity = sum(sum(sum(rho)));
rho_FT = fftshift(fftn(rho));

global Density_rec
global DensityFT_rec
%Initialize matrices that change size the +1 is to include initial density
if ParamObj.SaveMe
    Density_rec       = zeros( Nx, Ny, Nm, TimeObj.N_record + 1 );      % Store density amplitudes
    DensityFT_rec      = zeros( Nx, Ny, Nm, TimeObj.N_record + 1 );      % Store k-space amplitdues
    
    %Initialize records
    DensityFT_rec(:,:,:,1)   = rho_FT;
    Density_rec(:,:,:,1)     = rho;
else
    Density_rec = 0;
    DensityFT_rec = 0;
end

% keyboard
j_record = 2;     %Record holder

%Set up Diffusion operator, discrete k-space propagator, and interaction
%Set up Diffusion operator in cube form
[Lop] = DiffOpBuilderIsoDiffCube(DiffMobObj,GridObj);
Prop = exp(Lop .* TimeObj.delta_t);   % Exponentiate the elements


%%%%%%%%%%%%%%%%%%%Mayer function stuff%%%%%%%%%%%%%%%%%%%%%%%%%%
Fm_FT = fftshift(fftn( MayerFncDiffBtwPntsCalc(...
    Nx, Ny, Nm, ParamObj.Lx, ParamObj.Ly, ParamObj.L_rod) ));


%Hard rod interactions
if ParamObj.Interactions
    GammaExCube_FT = dRhoIntCalcVcFtId(rho,rho_FT,Fm_FT,ParamObj,...
        GridObj,DiffMobObj);
else
    GammaExCube_FT = zeros(Nx,Ny,Nm);
    fprintf(wfid,'Interacts are off my Lord\n');
end

%Driven Term
if ParamObj.Drive
    GammaDrCube_FT  = ...
        dRhoDriveCalcFtId(rho,ParamObj.vD,...
        GridObj.phi3D,GridObj.kx3D,GridObj.ky3D);
else
    GammaDrCube_FT = zeros(Nx,Ny,Nm);
    fprintf(wfid,'Driving is off my Lord\n');
end

%Total
GammaCube_FT = GammaDrCube_FT + GammaExCube_FT ;

% Take the first step- Euler. Element by element mulitplication
  % Take a step
if( ParamObj.StepMeth == 0 ) 
  NlPf =  TimeObj.delta_t;
 %[rho_FTnext] = DenStepperAB1c( Prop, rho_FT, GammaCube_FT, TimeObj.delta_t );
 [rho_FTnext] = DenStepperAB1cPf( Prop, rho_FT, GammaCube_FT, NlPf );
elseif( ParamObj.StepMeth == 1 )
  NlPf = 3 * TimeObj.delta_t / 2;
  NlPrevPf = TimeObj.delta_t / 2;
%   [rho_FTnext] = DenStepperAB1c( Prop, rho_FT, GammaCube_FT,TimeObj.delta_t );
  [rho_FTnext] = DenStepperAB1cPf( Prop, rho_FT, GammaCube_FT,NlPf );
elseif( ParamObj.StepMeth == 2 ) 
  NlPf = TimeObj.delta_t .* Prop;
%   [rho_FTnext] = DenStepperHAB1c( Prop, rho_FT, GammaCube_FT,TimeObj.delta_t );
  [rho_FTnext] = DenStepperHAB1cPf( Prop, rho_FT, GammaCube_FT, NlPf);
elseif( ParamObj.StepMeth == 3 )
  NlPf = 3 * TimeObj.delta_t / 2 .* Prop;
  NlPrevPf = TimeObj.delta_t / 2 .* Prop .* Prop;
%   [rho_FTnext] = DenStepperHAB1c( Prop, rho_FT, GammaCube_FT,TimeObj.delta_t );
  [rho_FTnext] = DenStepperHAB1cPf( Prop, rho_FT, GammaCube_FT, NlPf );
elseif( ParamObj.StepMeth == 4 )
  NlPf = TimeObj.delta_t / 2 .* ( 1 + Prop);
%   [rho_FTnext] = DenStepperBHAB1c( Prop, rho_FT, GammaCube_FT,TimeObj.delta_t );
  [rho_FTnext] = DenStepperBHAB1cPf( Prop, rho_FT, GammaCube_FT, NlPf );
elseif( ParamObj.StepMeth == 5 )
  NlPf = TimeObj.delta_t / 2 * ( 2 + Prop );
  NlPrevPf = TimeObj.delta_t / 2;
%   [rho_FTnext] = DenStepperBHAB1c( Prop, rho_FT, GammaCube_FT,TimeObj.delta_t );
  [rho_FTnext] = DenStepperBHAB1cPf( Prop, rho_FT, GammaCube_FT, NlPf);
elseif( ParamObj.StepMeth == 6 ) 
  GamProp = ( Prop - 1 ) ./ Lop;
  [rho_FTnext] = DenStepperEEM1c( Prop, GamProp, rho_FT,GammaCube_FT);
%   keyboard
else
  fprintf('No stepping method selected');
end


tic
ShitIsFucked = 0;
SteadyState  = 0;

% keyboard
fprintf(lfid,'Starting master time loop\n');
for t = 1:TimeObj.N_time-1
    %Save the previous and take one step forward.
    % Save the old drho
    GammaCube_FTprev = GammaCube_FT;
    rho_FTprev  = rho_FT;
    
    %Need to update rho!!!
    rho_FT      = rho_FTnext;
    
    % Calculate rho if there is driving or interactions
    if ParamObj.Interactions || ParamObj.Drive
        rho    = real(ifftn(ifftshift(rho_FT)));
    end
    
    %Hard rod interactions
    if ParamObj.Interactions
        GammaExCube_FT = dRhoIntCalcVcFtId(rho,rho_FT,Fm_FT,ParamObj,...
            GridObj,DiffMobObj);
    end
    
    %Driven Term
    if ParamObj.Drive
        GammaDrCube_FT  = dRhoDriveCalcFtId(...
            rho,ParamObj.vD,GridObj.phi3D,GridObj.kx3D,GridObj.ky3D);
    end
    
    
    GammaCube_FT = GammaDrCube_FT + GammaExCube_FT ;
    % Take a step
    if( ParamObj.StepMeth == 0 ) 
%         [rho_FTnext] = DenStepperAB1c( Prop, rho_FT, GammaCube_FT, TimeObj.delta_t  );
        [rho_FTnext] = DenStepperAB1cPf( Prop, rho_FT, GammaCube_FT, NlPf );
    elseif( ParamObj.StepMeth == 1 )
%       [rho_FTnext] = DenStepperAB2c( ...
%         Prop, rho_FT,GammaCube_FT,GammaCube_FTprev,TimeObj.delta_t );
        [rho_FTnext] = DenStepperAB2cPf( ...
             Prop, rho_FT,GammaCube_FT,GammaCube_FTprev,NlPf, NlPrevPf );
    elseif( ParamObj.StepMeth == 2 ) 
%       [rho_FTnext] = DenStepperHAB1c( Prop, rho_FT, GammaCube_FT,TimeObj.delta_t );
      [rho_FTnext] = DenStepperHAB1cPf( Prop, rho_FT, GammaCube_FT, NlPf);
    elseif( ParamObj.StepMeth == 3 )
%       [rho_FTnext] = DenStepperHAB2c( ...
%         Prop, rho_FT, GammaCube_FT,GammaCube_FTprev, TimeObj.delta_t );
        [rho_FTnext] = DenStepperHAB2cPf( ...
            Prop, rho_FT, GammaCube_FT,GammaCube_FTprev, NlPf, NlPrevPf );
    elseif( ParamObj.StepMeth == 4 )
%       [rho_FTnext] = DenStepperBHAB1c( Prop, rho_FT, GammaCube_FT,TimeObj.delta_t );
      [rho_FTnext] = DenStepperBHAB1cPf( Prop, rho_FT, GammaCube_FT, NlPf );
    elseif( ParamObj.StepMeth == 5 )
%       [rho_FTnext] = DenStepperBHAB2c( ...
%         Prop, rho_FT, GammaCube_FT,GammaCube_FTprev,TimeObj.delta_t );
      [rho_FTnext] = DenStepperBHAB2cPf( ...
        Prop, rho_FT, GammaCube_FT,GammaCube_FTprev, NlPf, NlPrevPf );
    elseif( ParamObj.StepMeth == 6 ) 
      [rho_FTnext] = DenStepperEEM1c( Prop, GamProp, rho_FT,GammaCube_FT);
    end

   
    %Save everything (this includes the initial state)
    if (mod(t,TimeObj.N_count)== 0)
        if ParamObj.SaveMe
            [SteadyState,ShitIsFucked,MaxReldRho] = ...
                VarRecorderTrackerCube(wfid,tfid,TimeObj,t,...
                rho_FT,rho_FTprev,TotalDensity ,j_record);
        else
            [SteadyState,ShitIsFucked,MaxReldRho] = ...
                VarRecorderTrackerNoSaveCube(wfid,tfid,TimeObj,t,...
                rho_FT,rho_FTprev,TotalDensity);
        end %if save
        if ShitIsFucked == 1 || SteadyState == 1
%             keyboard
            break
        end
        j_record = j_record+1;
    end %end recording
    
end %end time loop
fprintf(lfid,'Finished master time loop\n');
%  keyboard

% Update last rho
if ParamObj.SaveMe
    if ShitIsFucked == 0 && SteadyState == 0
        t =  t + 1;
        rho_FT      = rho_FTnext;
        if (mod(t,TimeObj.N_count)==0)
            if ParamObj.SaveMe
                [SteadyState,ShitIsFucked,MaxReldRho] = ...
                    VarRecorderTrackerCube(wfid,tfid,TimeObj,t,...
                    rho_FT,rho_FTprev,TotalDensity ,j_record);
            else
                [SteadyState,ShitIsFucked,MaxReldRho] = ...
                    VarRecorderTrackerNoSaveCube(wfid,tfid,TimeObj,t,...
                    rho_FT,rho_FTprev,TotalDensity);
            end %if save   [SteadyState,ShitIsFucked,MaxReldRho] = ...
            
        end % End recording
    end
end %end if save

%If something broke, return zeros. Else, return the goods
if ShitIsFucked
    fprintf(wfid,'Density is either negative or not conserved.\n');
    fprintf(wfid,'I have done %i steps out of %i.\n',t, TimeObj.N_time);
    
elseif SteadyState
    fprintf(wfid,'Things are going steady if you know what I mean.\n');
    fprintf(wfid,'I have done %i steps out of %i.\n',t, TimeObj.N_time);
end

% Get rid of zeros in record matrices
Record_hold   = 1:j_record;
TimeRecVec    = (0:j_record-1) * TimeObj.t_record;
if ParamObj.SaveMe
    Density_rec   = Density_rec(:,:,:,Record_hold);
    DensityFT_rec = DensityFT_rec(:,:,:,Record_hold);
else
    Density_rec   = rho;
    DensityFT_rec = rho_FT;
end %end if save

trun = toc;

%Save the structure
DenRecObj = struct('DidIBreak', ShitIsFucked,'SteadyState', SteadyState,...
    'j_record', j_record,...
    'TimeRecVec',TimeRecVec,...
    'RunTime', trun, ...
    'bc',ParamObj.bc,...
    'Density_rec',Density_rec,'DensityFT_rec', DensityFT_rec,...
    'MaxReldRho',MaxReldRho,'feq',feq);


fclose(wfid); %Close warning statement file
fclose(tfid); %Close program tracker file
% keyboard
end %function
