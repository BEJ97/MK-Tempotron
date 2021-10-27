function recon=myrecon(edge,M,N,circularRField)
if(nargin<4),
    circularRField=0;
end

recon=zeros(M,N);
for edgeidx=1:size(edge,1)  % for each line
    rr=edge(edgeidx,1); % center point's row coordinate
    cc=edge(edgeidx,2); % center point's col coordinate
    szidx=(edge(edgeidx,3)-1)/2; % szidx
    diridx=edge(edgeidx,4); % diridx
    switch diridx
        case 1  % horizontal
            startpoint_col=max(1,cc-szidx);
            endpoint_col=min(N,cc+szidx);
            recon(rr, startpoint_col:endpoint_col)=1;
        case 2  % 45
			recon(rr,cc)=1;
            if(circularRField),limit=floor(szidx/sqrt(2));else limit=szidx; end
            % circularRField limit LUT: [if limit=floor(szidx/sqrt(2)) ]
            % szidx : 1  2  3  4  5  6 
            % limit : 0  1  2  2  3  4
            for i=1:limit
                if (rr+i>=1)&&(rr+i<=M)&& (cc-i>=1)&&(cc-i<=N)
                    recon(rr+i,cc-i)=1;
                end
                if (rr-i>=1)&&(rr-i<=M)&& (cc+i>=1)&&(cc+i<=N)
                    recon(rr-i,cc+i)=1;
                end
            end
        case 3  % vertical
            startpoint_row=max(1,rr-szidx);
            endpoint_row=min(M,rr+szidx);
            recon(startpoint_row:endpoint_row, cc)=1;
        case 4  % 135
            recon(rr,cc)=1;
            if(circularRField),limit=floor(szidx/sqrt(2));else limit=szidx; end
            for i=1:limit
                if (rr+i>=1)&&(rr+i<=M)&& (cc+i>=1)&&(cc+i<=N)
                    recon(rr+i,cc+i)=1;
                end
                if (rr-i>=1)&&(rr-i<=M)&& (cc-i>=1)&&(cc-i<=N)
                    recon(rr-i,cc-i)=1;
                end
            end
    end
end

