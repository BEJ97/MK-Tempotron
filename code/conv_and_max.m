function [ AllVec,Time_Chnl_Lbl ] = conv_and_max( CIN,numScales,numRot,numFilters,size1,size2,...
            Timings,paramCell,stateCell, timeslice, ...
            sqfilters,filtersize, doMAXforEachEvent, post_proc_th,  ...
            showOnTheFlyConv,showLayerFigures,rot,circularRField,conv_disp_interval,...
            CONV_refresh_window_only, useCentroidShift)

time_last_disp = -inf;
cnt_disp = 0;

N = size(CIN,1);    
S1 = cell(numScales,numRot);
C1 = cell(numScales,numRot);
for i = 1:numScales
    for j = 1:numRot
        S1{i,j} = stateCell{i,j}.J;     % all zeros at the very beginning.
        C1{i,j} = stateCell{i,j}.J_max; % all zeros at the very beginning.
    end
end
AllVec = zeros( size1*size2*numFilters, size(Timings,2) );%Timings 是超过阀值时的时间点

if showLayerFigures,
    S1C1Edge_fig = figure('position',[10 60 1000 500]);
else
    S1C1Edge_fig = [];
end

if showOnTheFlyConv
    ontheflyconv_fig = figure('position',[10 500 700 300]);
else
    ontheflyconv_fig = [];
end


for iEvt = 1:N
    
    event_in = CIN(iEvt,:);
    time = CIN(iEvt,1);    
    if iEvt>1, timelast = CIN(iEvt-1,1); else timelast = -inf; end

    % -------- on the fly convolution.--------
    for i = 1:numScales
        for j = 1:numRot
            [stateCell{i,j}] = conv_integNoFire(event_in, ...
                paramCell{i,j},stateCell{i,j},CONV_refresh_window_only);
        end
    end
    
    if doMAXforEachEvent,
        % --------MAX operation. each event affects a region --------
        [ stateCell ] = myMAXoperation_win( stateCell,filtersize, iEvt,event_in,sqfilters,post_proc_th);
    end
    
    if (showOnTheFlyConv)
        if (time-time_last_disp >= conv_disp_interval) || time == CIN(end,1)
            cnt_disp = cnt_disp+1;
            show_S1_C1(ontheflyconv_fig,stateCell);
            time_last_disp = time;
        end
    end
    
    
    % ====== at timings selected by motion symbol detector : =============
    e = find(time == Timings(1,:));
    if (~isempty(e)) && (time~=timelast )

        if CONV_refresh_window_only, 
            % conv only affected different small windows. now we need to 
            % refresh S1 whole map before whole-map MAX (for doMAXforEachEvent==0).
            % Actually, here it also refreshes C1 maps (for doMAXforEachEvent==1)
            for i = 1:numScales
                for j = 1:numRot
                    [stateCell{i,j}]= leakage_refresh_whole_map(paramCell{i,j},stateCell{i,j},time);
                end
            end
        end
        
        if ~doMAXforEachEvent,
            % --------MAX operation. whole map --------
            [ stateCell ] = myMAXoperation( stateCell,filtersize,post_proc_th);
        end
        
        for i = 1:numScales
            for j = 1:numRot
                S1{i,j} = abs( stateCell{i,j}.J );
                C1{i,j} = stateCell{i,j}.J_max;
                % thresholding on the refreshed C1 map: ******
                C1{i,j}=C1{i,j}.*(C1{i,j}>=post_proc_th);
            end
        end
 
        if showLayerFigures,
            % input frame recon.
            tmp = ( CIN(:,1)<=CIN(iEvt,1) ) & ( CIN(:,1)>=CIN(iEvt,1)-timeslice ) ;
            I1 = reconstaer_back(CIN(tmp,:), size1, size2);
            
            show_S1_C1_I_edge_figs(S1C1Edge_fig,rot,circularRField,I1,S1,C1,filtersize);
        end

        if (useCentroidShift),
            if size(Timings,1)<5
                warning('No Centroid information! centroidShift disabled.');
            else
                centroid = Timings(4:5,e);
                shiftHorizontal = 1;
                shiftVertical = 0;
                [ C1 ] = centroidShift(C1, centroid, shiftHorizontal, shiftVertical);
            end
        end
            
        [ Vec] = myRasterScan(C1);
        for inde = 1:length(e)
            AllVec(:,e(inde)) = Vec;
        end
    end
    
    
end
    
Time_Chnl_Lbl = Timings;

