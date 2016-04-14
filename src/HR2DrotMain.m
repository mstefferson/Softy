% HR2DrotMain.m
%
% Program is the main() for running the diffusion of 2D hard rods with
% orientation. Program handles interactions using DDFT

function  [DenFinal, DenFTFinal, GridObj, ParamObj,TimeObj,...
    DidIBreak,SteadyState,MaxReldRho] = ...
    HR2DrotMain(InputFile)
% Add paths (this should already be added, but just to be careful)
% Save error messages in file
try
    EvolvedDen = 0;DenFinal = 0;DenFTFinal = 0;GridObj = 0;ParamObj = 0;
    TimeObj = 0;DidIBreak = 0;SteadyState = 0;MaxReldRho = 0;
    
    %     keyboard
    tMainID  = tic;
    %Grab the parameters
    % keyboard
    DataTemp    = importdata(InputFile);
    ParamVec    = DataTemp.data(1,:);
    TimeVec     = DataTemp.data(2,~isnan(DataTemp.data(2,:)));    %Pull out the time stuff
    FileNameMat = DataTemp.textdata(1);
    Path2Save   = DataTemp.textdata(2);
    IntDenType  = DataTemp.textdata(3);
    LoadName    = DataTemp.textdata(4);
    
    % Make some  objects
    %     ParamNmVec = {'trial' 'Interactions' 'Nx' 'Ny' 'Nm' 'Lx' 'Ly' 'L_rod' 'Diam' 'Eta_visc'...
    %         'Tmp' 'Norm' 'WeightPos' 'WeightAng' 'NumModesX' 'NumModesY' 'NumModesM' 'bc'};
    ParamNmVec = {'trial' 'Interactions' 'Drive' 'StepMeth' 'IntCond' 'Movies' 'SaveMe'...
        'Nx' 'Ny' 'Nm' 'Lx' 'Ly' 'L_rod' 'Diam' 'Eta_visc'...
        'kB' 'Tmp' 'Norm' 'WeightPos' 'WeightAng' 'NumModesX' 'NumModesY' ...
        'NumModesM' 'bc' 'Mob_pos','Mob_rot'  };
    TimeNmVec  = {'delta_t' 't_record' 't_tot' 'ss_epsilon'};
    
    ParamObj   = struct('NmVec',{ParamNmVec},'ValVec',...
        ParamVec,'trial',ParamVec(1),...
        'Interactions',ParamVec(2), 'Drive',ParamVec(3),...
        'StepMeth', ParamVec(4), 'IntCond', ParamVec(5),...
        'MakeOP',ParamVec(6),'MakeMovies',ParamVec(7),'SaveMe',ParamVec(8),...
        'Nx', ParamVec(9),'Ny', ParamVec(10),'Nm', ParamVec(11),...
        'Lx', ParamVec(12),'Ly', ParamVec(13),'L_rod', ParamVec(14), ...
        'Tmp',ParamVec(15), 'Norm',ParamVec(16), 'WeightPos',ParamVec(17), ...
        'WeightAng',ParamVec(18), 'Random',ParamVec(19),...
        'NumModesX',ParamVec(20), 'NumModesY',ParamVec(21), ...
        'NumModesM', ParamVec(22),'bc',ParamVec(23), ...
        'c', ParamVec(24), ...
        'Mob_par',ParamVec(25),'Mob_perp',ParamVec(26),...
        'Mob_rot',ParamVec(27), 'vD', ParamVec(28) );
    
    fprintf('Read input and made ParamObj\n')
    
    %     keyboard
    % Create a file that holds warning print statements
    WarningStmtString = sprintf('WarningStmts_%i.txt',ParamObj.trial);
    wfid  = fopen(WarningStmtString,'a+');    % a+ allows to append data
    
    LocString = sprintf('Location_%i.txt',ParamObj.trial);
    lfid      = fopen(LocString,'a+');    % a+ allows to append data
    fprintf(lfid,'Starting main, current code\n');
    
    %Time Recording
    N_time   = ceil(TimeVec(3)/TimeVec(1)); %number of time steps
    N_record = ceil(TimeVec(3)/TimeVec(2)); %number of time points to record. Does not include initial density
    N_count  = ceil(TimeVec(2)/TimeVec(1)); %spacing between times to record
    
    TimeObj = struct('TimeNmVec',{TimeNmVec},'TimeVecOrg',{TimeVec},...
        'delta_t',TimeVec(1), 't_record',TimeVec(2), 't_tot',TimeVec(3), 'ss_epsilon',TimeVec(4),...
        'N_time', N_time, 'N_record',N_record,'N_count',N_count);
    
    % Fix the time
    [TimeObj.t_tot,TimeObj.N_time,TimeObj.t_rec,TimeObj.N_rec,TimeObj.N_count]= ...
        TimeStepRecMaker(TimeObj.delta_t,TimeObj.t_tot,TimeObj.t_record);
    fprintf(lfid,'Made Time Obj\n');
    fprintf('Made Time Obj\n')
    %%%Make all the grid stuff%%%%%%%%%%%%%%
    [GridObj] = GridMakerPBCxk(...
        ParamObj.Nx,ParamObj.Ny,ParamObj.Nm,ParamObj.Lx,ParamObj.Ly);
    fprintf(lfid,'Made Grid\n');
    fprintf('Made Grid\n')
    
    %Make diffusion coeff (send smallest dx dy for stability
    [DiffMobObj] =  DiffMobCoupCoeffCalc( wfid,ParamObj.Tmp,...
        ParamObj.Mob_par,ParamObj.Mob_perp,ParamObj.Mob_rot,...
        TimeObj.delta_t, min(GridObj.dx,GridObj.dy),...
        GridObj.dphi,GridObj.kx2D, GridObj.ky2D,ParamObj.vD);
    
    fprintf(lfid,'Made diffusion object\n');
    fprintf('Made diffusion object\n');
    
    %Initialze density
    [rho] = MakeConc(GridObj,ParamObj);
    Nc    = 20;
    % Equilib distribution
    [Coeff_best,~] = CoeffCalcExpCos2D(Nc,GridObj.phi,ParamObj.bc); % Calculate coeff
    feq = DistBuilderExpCos2Dsing(Nc,GridObj.phi,Coeff_best);        % Build equil distribution
    fprintf(lfid,'Made initial density\n');
    fprintf('Made initial density\n');
    
    %     keyboard
    % Run the main code
    tBodyID      = tic;
    
    [DenRecObj]  = HR2DrotDenEvolverFTBody(wfid,lfid,rho,ParamObj, TimeObj,GridObj,DiffMobObj,feq);
    %     keyboard
    EvolvedDen = 1;
    BodyRunTime  = toc(tBodyID);
    fprintf(lfid,'Made density object\n');
    fprintf(lfid,'Body Run Time = %f\n\n', BodyRunTime);
    fprintf('Made density object\n');
    fprintf('Body Run Time = %f\n\n', BodyRunTime);
    
    % Store final density and transform
    DenFinal   = DenRecObj.Density_rec(:,:,:,end);
    DenFTFinal = DenRecObj.DensityFT_rec(:,:,:,end);
    DidIBreak  = DenRecObj.DidIBreak;
    SteadyState = DenRecObj.SteadyState;
    MaxReldRho  = DenRecObj.MaxReldRho;
    
    %         keyboard
    
    % Run movies if you want
    if ParamObj.MakeOP  == 1
        tOpID           = tic ;
        %                 keyboard
        if  DenRecObj.DidIBreak == 0
            [OrderParamObj] = CPNrecMaker(...
                ParamObj.Nx,ParamObj.Ny,DenRecObj.TimeRecVec,...
                GridObj,DenRecObj.Density_rec,feq);
        else %Don't incldue the blowed up denesity for movies. They don't like it.
            TimeRecVecTemp = DenRecObj.TimeRecVec(1:end-1);
            [OrderParamObj] = CPNrecMaker(ParamObj.Nx,ParamObj.Ny,...
                TimeRecVecTemp,GridObj,...
                DenRecObj.Density_rec(:,:,:,1:length(TimeRecVecTemp)),...
                feq);
        end
        OpRunTime       = toc(tOpID);
        fprintf(lfid,'Made interaction order paramater object\n');
        fprintf('Made interaction order paramater object\n');
        
        if ParamObj.MakeMovies == 1
            % Build OP records
            
            % Make matlab movies
            tMovID       = tic;
            %         keyboard
            HoldX = ParamObj.Nx /2 + 1;
            HoldY = ParamObj.Ny /2 + 1;
            
            if DenRecObj.DidIBreak == 0
          
                OPMovieMakerTgtherAvi(ParamObj.trial,GridObj.x,GridObj.y, GridObj.phi, ...
                OrderParamObj.C_rec, OrderParamObj.NOP_rec,OrderParamObj.POP_rec,...
                reshape( DenRecObj.Density_rec(HoldX, HoldY, : , :), [ParamObj.Nm length(DenRecObj.TimeRecVec)] ),...
                DenRecObj.TimeRecVec)
            
            else
                
                OPMovieMakerTgtherAvi(ParamObj.trial,GridObj.x,GridObj.y, GridObj.phi, ...
                OrderParamObj.C_rec, OrderParamObj.NOP_rec,OrderParamObj.POP_rec,...
                reshape( DenRecObj.Density_rec(HoldX, HoldY, : ,1 :end - 1), [ParamObj.Nm length(DenRecObj.TimeRecVec) - 1] ),...
                DenRecObj.TimeRecVec(1:end - 1 ) )
                
            end


            %        MovieObj = OPMovieMakerTgtherMat(GridObj,ParamObj,OrderParamObj,...
            %              DenRecObj.Density_rec,feq);
            %                 keyboard
            MovRunTime   = toc(tMovID);
            %         keyboard
            fprintf(lfid,'Made movies\n');
            % Record how long it took
            fprintf(lfid,'OrderParam Run time = %f\n', OpRunTime);
            fprintf(lfid,'Make Mov Run Time = %f\n',  MovRunTime);
            
            
        end % End if movies
        
        kx0 = ParamObj.Nx / 2 + 1;
        ky0 = ParamObj.Ny / 2 + 1;
        km0 = ParamObj.Nm / 2 + 1;
        Nrec = length( DenRecObj.TimeRecVec);

        FTind2plot = zeros( 8, 3 );
        FTmat2plot = zeros( 8, Nrec );
        
        FTind2plot(1,:) = [kx0     ky0     km0 + 1];
        FTind2plot(2,:) = [kx0 + 1 ky0     km0 + 1];
        FTind2plot(3,:) = [kx0     ky0 + 1 km0 + 1];
        FTind2plot(4,:) = [kx0 + 1 ky0 + 1 km0 + 1];
        FTind2plot(5,:) = [kx0     ky0     km0 + 2];
        FTind2plot(6,:) = [kx0 + 1 ky0     km0 + 2];
        FTind2plot(7,:) = [kx0     ky0 + 1 km0 + 2];
        FTind2plot(8,:) = [kx0 + 1 ky0 + 1 km0 + 2];
        
        for i = 1:8 
            FTmat2plot(i,:) =  reshape(... 
            DenRecObj.DensityFT_rec( FTind2plot(i,1), FTind2plot(i,2), FTind2plot(i,3),: ),...   
            [ 1, Nrec ]  );   
        end
        
%         keyboard
        ampPlotterFT(FTmat2plot, FTind2plot, DenRecObj.TimeRecVec, ParamObj.Nx, ParamObj.Ny,...
            ParamObj.Nm, DenRecObj.bc,ParamObj.vD,  ParamObj.SaveMe, ParamObj.trial)
   
    end % if OP
    
    if ParamObj.SaveMe
        MemObj = 0;
        % Save all parameters
        
        % Save everything. Save seperately for big files
        DenStr = sprintf('DenRec_%i',ParamObj.trial);
        TimeStr = sprintf('TimeObj_%i',ParamObj.trial);
        ParamStr = sprintf('ParamObj_%i',ParamObj.trial);
        GridStr = sprintf('GridObj_%i',ParamObj.trial);
        
        save(DenStr,'DenRecObj','-v7.3')
        save(TimeStr,'GridObj','-v7.3')
        save(ParamStr,'ParamObj','-v7.3')
        save(GridStr,'GridObj','-v7.3')
        
        if ParamObj.MakeOP
            OpStr = sprintf('OP_%i',ParamObj.trial);
            save(OpStr,'OrderParamObj','-v7.3')
        end
    end
    % Save how long everything took
    fprintf(lfid,'Everything saved\n');
    TotRunTime = toc(tMainID);
    fprintf(lfid,'Total Run time = %f\n', TotRunTime);
    
    fclose('all');
    % keyboard
    %     if ParamObj.SaveMe
    % Move everything
    %         MoveStrTxt = sprintf('*%i.txt', ParamObj.trial);
    %         MoveStrMat = sprintf('*%i.mat', ParamObj.trial);
    %         movefile(MoveStrTxt,Path2Save{1});
    %         movefile(MoveStrMat,Path2Save{1});
    %         cd /home/mws/Documents/Research/BG/DDFT/Outputs
    %     end
    
catch err %Catch errors
    
    
    ErrFileNmStr = sprintf('errFile%i.txt',ParamObj.trial);
    efid         = fopen(ErrFileNmStr,'a+');
    % write the error to file and to screen
    % first line: message
    %     fprintf(efid,'%s', err.getReport('extended', 'hyperlinks','off')) ;
    fprintf('%s', err.getReport('extended')) ;
    disp(err.message);
    fclose(efid);
    fclose('all');
    
    keyboard
    %    keyboard
    if ParamObj.SaveMe
        
        TimeStr = sprintf('TimeObj_%i',ParamObj.trial);
        ParamStr = sprintf('ParamObj_%i',ParamObj.trial);
        GridStr = sprintf('GridObj_%i',ParamObj.trial);
        
        save(TimeStr,'GridObj','-v7.3')
        save(ParamStr,'ParamObj','-v7.3')
        save(GridStr,'GridObj','-v7.3')
        if EvolvedDen
            DenStr = sprintf('DenRec_%i',ParamObj.trial);
            save(DenStr,'DenRecObj','-v7.3');
        end
    end
    
end %End try and catch

% clc
%close all
end % End HR2DrotVgrExeMain.m
