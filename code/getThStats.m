function [ THs ] = getThStats(CINsfile,V0,dt,use_single_exponential, lut1,lut2, nGroup,datafolder)

if ~exist('datafolder','var')
    datafolder = pwd;
end

timedLog('getting Th stats...');

load (CINsfile,'NumCINs','CINs','Labels');
N = NumCINs;
th = zeros(1,N);

for i = 1:N
    
    CIN = CINs{i};
    
    ts = CIN(:,1);
    [ V ]= edInteg (ts,V0,dt,use_single_exponential,lut1,lut2); %%%%%%%%%V即为MSD的膜电压   
    th(i) = max(V);   %每个样本的最大膜电位	
end
save([datafolder,'/','Vmsd'],'th');
                  %%%%%
%THs2 = zeros(1200,nGroup); 
%k = 1;                         %%%%%%%%
%for i = 1:nGroup                                    %%%%%%%%%%%%
%    ind = Labels == i-1;
%    THs2(ind,i) = th(ind);                          
%end
%THs3 = reshape(th, [7000,10]);                                                 %%%%%%%%%%%
%save([datafolder,'/','Vmax'],'THs3');


THs = zeros(1,nGroup);
%for i = 1:nGroup
%    ind = Labels == i-1;
%    THs(i) = min(th(ind))*0.3;	%改  原THs(i) = min(th(ind))*0.3       %min(th(ind))该类别某个最小peak的样本
	% if THs(i) < 150            %改  原没
		% THs(i) = 150;           %
	% end                          %
%end
% THs(1) = 220;
% THs(2) = 155;
% THs(3) = 190;
% THs(4) = 170;
% THs(5) = 110;
% THs(6) = 290;
% THs(7) = 280;
% THs(8) = 250;
% THs(9) = 220;
% THs(10) = 220;


save ([datafolder,'/','THs'], 'THs')

    