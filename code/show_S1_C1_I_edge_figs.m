function show_S1_C1_I_edge_figs(S1C1Edge_fig,rot,circularRField,I1,s1,c1,filtersize)


[numScales, numRot] = size(s1);

% extract edges, no merge
edge=[];
edgeidx=0;
for idxRot=1:numRot
    for idxScale=1:numScales  
        tmp=c1{idxScale,idxRot};
        [M,N] = size(tmp);
        fsz = filtersize((idxScale-1)*numRot+idxRot);
        for rr=1:M
            for cc=1:N
                if (tmp(rr,cc)~=0) % a surviving neuron
                    edgeidx=edgeidx+1;
                    switch rot(idxRot)
                        case 0
                            rot_type = 1;
                        case 45
                            rot_type = 2;
                        case 90
                            rot_type = 3;
                        case 135
                            rot_type = 4;
                    end
                    edge(edgeidx,:)=[rr,cc,fsz,rot_type];
                    %[center_row, center_col, length, orientation_type]
                end
            end
        end
    end
end

% lines reconstruction
recon=myrecon(edge,M,N,circularRField);

figure(S1C1Edge_fig),
cont = 0;

for idxScale = 1:numScales
    for idxRot = 1:numRot*2+1
        cont = cont+1;
        if idxRot<= numRot
            subplot(numScales,numRot*2+1,cont),
            imshow(s1{idxScale,idxRot},[]);
            if (idxScale==1)&&(idxRot==1), title('s1'); end
            drawnow;
        elseif idxRot<=numRot*2
            subplot(numScales,numRot*2+1,cont),
            imshow(c1{idxScale,idxRot-numRot}~=0,[]); 
            if (idxScale==1)&&(idxRot==numRot+1), title('c1'); end
            drawnow;
        else
            if idxScale==1
                subplot(numScales,numRot*2+1,cont),
                imshow(I1,[]); 
                title('input'); 
                drawnow;
            end
            if idxScale==2
                subplot(numScales,numRot*2+1,cont),
                imshow(recon,[]); 
                title('edge recon'); 
                drawnow;
            end
        end
    end
end





