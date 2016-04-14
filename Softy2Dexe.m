% Input creater for HR2DrotMainIdC
% Isotropic diffusion

CurrentDir = pwd;
addpath( genpath( CurrentDir) );

% Now can change number of grid points in the x, y, phi direction
Run  = 1; % Run main from here
Move = 0; % Move files to a nice location

%%%%%%%% Trial %%%%%%%%%%%%
trial    = 7;

%%%%%% Turn on/off interactions%%%%%%%%%
Interactions = 1;
SaveMe       = 1;

%%%%%%%%%%%%% Box and Rod Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
Nx      = 128;
Ny      = 128;

%%%%%%%%% Initial density parameters%%%%%%%%%%%%%%%%%%
% Dimensionless  scaled concentration bc > 1.501 or bc < 1.499 if
% perturbing about equilbrum
bc  = 4.0;     % Scaled concentration
R   = 1;        % Disk raidus
Rs  = 1.855;        % Soft shoulder distance
Lx  = 10*R;     % Box length
Ly  = 10*R;     % Box length

%%%%%%%%%%%%%%%Time recording %%%%%%%%%%%%%%%%%%%%%%%%%%
dt          = 1e-3; %time step
t_rec       = 0.5e-1; %time interval for recording dynamics
t_tot       = 10;   %total time
ss_epsilon  = 1e-8;                          %steady state condition

NumModesX   = 4;
NumModesY   = 4;

% Weight of the spatial sinusoidal perturbation. %
% Perturbations added to rho(i,j,k) = 1. Must be small
WeightPos   = 1e-2;
Random      = 0;       % Random perturbation coeffs

% Stepping method
StepMeth = 0;  % Not in yet

%%%%%%%%%%%%%%%%%%%%% Physical Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tmp    = 1;            % Temperature
eps    = 1.0;
a      = 0.1;

% mobility
Mob = 1;
Mob_pos  = Mob;
Mob_rot  = Mob;
D = Mob / Tmp;
% Concentration and rod stuff
b       = R ^2;            % Average excluded a guess for now.
c       = bc / b;                   % Concentration
Norm    = c * Lx * Ly;              % number of particles

ParamObj = struct( 'SaveMe',SaveMe,'Nx', Nx, 'Ny', Ny, 'Lx', Lx, 'Ly', Ly, 'eps', eps, ...
    'bc',bc,'a', a,'Tmp', Tmp, 'R', R, 'Rs', Rs) ;
GridObj = GridMakerPBCxk(Nx, Ny, 0, Lx, Ly);
[t_tot,Nt,t_rec,N_rec,N_count] = TimeStepRecMaker(dt,t_tot,t_rec);

TimeObj = struct( 'dt', dt, 't_rec', t_rec, 't_tot', t_tot, 'ss_epsilon',ss_epsilon,...
    'N_time', Nt, 'N_rec',N_rec,'N_count',N_count);

%% Initialize rho. Very rough
rho =  c .* ones(Nx,Ny);
kx0 = Nx/2+1;
ky0 = Ny/2+1;
for i = 1:NumModesX
    for j = 1:NumModesY
       rho = rho +  rand() * WeightPos * c * ( ...
            cos( GridObj.kx(kx0 + i) .* GridObj.x2D + GridObj.ky(ky0 + j) * GridObj.y2D ) + ...
            sin( GridObj.kx(kx0 + i) .* GridObj.x2D + GridObj.ky(ky0 + j) * GridObj.y2D ) );
    end
end
% keyboard
rho_FT = fftshift(fftn(rho));

global Density_rec
global DensityFT_rec
%Initialize matrices that change size the +1 is to include initial density
if ParamObj.SaveMe
    Density_rec       = zeros( Nx, Ny, TimeObj.N_rec + 1 );      % Store density amplitudes
    DensityFT_rec      = zeros( Nx, Ny, TimeObj.N_rec + 1 );      % Store k-space amplitdues
    
    %Initialize records
    DensityFT_rec(:,:,1)   = rho_FT;
    Density_rec(:,:,1)     = rho;
else
    Density_rec = 0;
    DensityFT_rec = 0;
end

% keyboard
j_record = 2;     %Record holder

%Set up Diffusion operator, discrete k-space propagator, and interaction
%Set up Diffusion operator in cube form
[Lop] = - D .* ( GridObj.kx2D .^ 2 + GridObj.ky2D .^ 2);
Prop = exp(Lop .* TimeObj.dt);   % Exponentiate the elements

%%%%%%%%%%%%%%%%%%% Interaction stuff%%%%%%%%%%%%%%%%%%%%%%%%%%
V    = SSpotential(Nx,Ny,Lx,Ly,R,Rs,eps, a);
V_FT = fftshift(fftn( V ) );

%Hard rod interactions
if Interactions
    GammaCube_FT  = dRhoIntCalc2D(rho,rho_FT,V_FT,ParamObj,GridObj,Mob);
else
    GammaCube_FT = zeros(Nx,Ny);
end

NlPf = TimeObj.dt / 2 .* ( 1 + Prop);
 
%NlPf = TimeObj.dt;
tic
ShitIsFucked = 0;
SteadyState  = 0;

% First step
[rho_FTnext] = DenStepperAB1cPf( Prop, rho_FT, GammaCube_FT, NlPf );
% [rho_FTnext] = DenStepperBHAB1cPf( Prop, rho_FT, GammaCube_FT, NlPf );

 NlPf = TimeObj.dt / 2 * ( 2 + Prop );
 NlPrevPf = TimeObj.dt / 2;
% keyboard
for t = 1:TimeObj.N_time-1
    %Save the previous and take one step forward.
    % Save the old drho
    GammaCube_FTprev = GammaCube_FT;
    rho_FTprev  = rho_FT;

    %Need to update rho!!!
    rho_FT      = rho_FTnext;
    
    %Hard rod interactions
    if Interactions
        rho    = real(ifftn(ifftshift(rho_FT)));
        GammaCube_FT  = dRhoIntCalc2D(rho,rho_FT,V_FT,ParamObj,GridObj,Mob);
    end
    
%[rho_FTnext] = DenStepperAB1cPf( Prop, rho_FT, GammaCube_FT, NlPf );
%[rho_FTnext] = DenStepperBHAB1cPf( Prop, rho_FT, GammaCube_FT, NlPf );
 [rho_FTnext] = DenStepperBHAB2cPf( Prop, rho_FT, GammaCube_FT,...
     GammaCube_FTprev, NlPf, NlPrevPf );
    %Save everything (this includes the initial state)
    if (mod(t,TimeObj.N_count)== 0)
        if SaveMe
             if ~Interactions
                 rho    = real(ifftn(ifftshift(rho_FT)));
         end
            Density_rec(:,:,j_record) = rho;
            DensityFT_rec(:,:,j_record) = rho_FT;
        end %if save
        if isnan(rho)
            fprintf('Shit is NAN-fucked\n');
            ShitIsFucked = 1;
        end
        if min(min(min(rho))) < 0
            fprintf('Shit is neg-fucked\n');
            ShitIsFucked = 1;
        end
        if ShitIsFucked == 1 || SteadyState == 1
            %keyboard
            break
        end
        j_record = j_record+1;
    end %end recording
    
end %end time loop


% Update last rho
if SaveMe
    t =  t + 1;
    rho_FT      = rho_FTnext;
    rho = real( ifftn( ifftshift( rho_FT) ) );
    if (mod(t,TimeObj.N_count)==0)
        Density_rec(:,:,j_record) = rho;
        DensityFT_rec(:,:,j_record) = rho_FT;
    end % End recording
end
%end if save

Record_hold   = 1:j_record;
TimeRecVec    = (0:j_record-1) * TimeObj.t_rec;
if SaveMe
    Density_rec   = Density_rec(:,:,Record_hold);
    DensityFT_rec = DensityFT_rec(:,:,Record_hold);
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
    'Density_rec',Density_rec,'DensityFT_rec', DensityFT_rec);

%%
if ShitIsFucked == 0
Fig = figure();
set(Fig, 'WindowStyle', 'normal');
%PosVec = [680 558 1200 800];
%Fig.Position = PosVec;

nFrames = length(TimeRecVec);
% keyboard
set(gcf,'renderer','zbuffer')

axh1 = subplot(1,2,1); % Save the handle of the subplot
colorbar('peer',axh1);
axpos1 = get(axh1,'position'); % Save the position as ax
set(axh1,'NextPlot','replaceChildren',...
    'CLim', [0 max(max(max(Density_rec)))],...
    'YDir','normal');
set(axh1,'position',axpos1); % Manually setting this holds the position with colorbar
axh1.XTick = 0:Lx/5:Lx; axh1.YTick = 0:Ly/5:Ly;
xlabel('x'); ylabel('y')
axis square



axh2 = subplot(1,2,2); % Save the handle of the subplot
colorbar('peer',axh2);
axpos2 = get(axh2,'position'); % Save the position as ax
set(axh2,'NextPlot','replaceChildren',...
    'CLim',[min(min(min(log( Density_rec ) ))) max(max(max(log( Density_rec ) )))],...
    'YDir','normal');
set(axh2,'position',axpos2); % Manually setting this holds the position with colorbar
axh2.XTick = 0:Lx/5:Lx; axh2.YTick = 0:Ly/5:Ly;
xlabel('x'); ylabel('y')
axis square


%Initialize the movie structure
MovStr = sprintf('SoftCrys%.d.avi',trial);
Mov = VideoWriter(MovStr);
Mov.FrameRate = 4;
open(Mov);
%F( length(TimeRecVec) ) = struct('cdata',[],'colormap',[]);

for i = 1:length(TimeRecVec)
    titstr = sprintf('c(x,y) t = %.1f', TimeRecVec(i) );
    title(axh1,titstr);
    pcolor( axh1,GridObj.x,GridObj.y,Density_rec(:,:,i) );
    shading(axh1,'interp')

    titstr = sprintf('log( c(x,y) ) t = %.1f', TimeRecVec(i) );
    title(axh2,titstr);
    pcolor( axh2,GridObj.x,GridObj.y, log(Density_rec(:,:,i))  );
    shading(axh2,'interp')
   

   % F = getframe(Fig,[0 0 PosVec(3) PosVec(4)]);
   F = getframe(Fig);
    writeVideo(Mov,F);
%     keyboard
end
    
close(Mov);

end %% ShitIsFucked plotter
