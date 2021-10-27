function show_S1_C1(h_fig,stateCell)

[numScales, numRot] = size(stateCell);

s1 = cell(numScales,numRot);
c1 = cell(numScales,numRot);

for i = 1:numScales
    for j = 1:numRot
        s1{i,j} = abs(stateCell{i,j}.J);
        c1{i,j} = stateCell{i,j}.J_max;
    end
end


figure(h_fig),
cont = 0;
nr = numScales;
nc = numRot*2;

for idxScale = 1:numScales
    for idxRot = 1:numRot*2
        cont = cont+1;
        if idxRot<= numRot
            subplot(nr,nc,cont),
            imshow(s1{idxScale,idxRot},[]);
            if (idxScale==1)&&(idxRot==1), title('s1'); end
            drawnow;
        elseif idxRot<=numRot*2
            subplot(nr,nc,cont),
            imshow(c1{idxScale,idxRot-numRot}~=0,[]); 
            if (idxScale==1)&&(idxRot==numRot+1), title('c1 b4Th'); end
            drawnow;
        end
    end
end







