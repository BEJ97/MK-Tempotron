
% Tempotron: a neuron that learns spike timing-based decisions
% Rober Gutig 2006 Nature Neuroscience
NumImages=2;%改
% for i=1:NumImages
  % ImageMatrix=zeros(10,1);%10*1的矩阵
	% for c=1:10
		% ImageMatrix(c)=randperm(255,1);
	% end
  % TrainPtns(:,i)=ImageMatrix;
% end
% save('TrainPtns','TrainPtns');
load('TrainPtns2');
TrainPtns=TrainPtns*1e-3;  % scale to ms
nAfferents = size(TrainPtns,1);%改
nPtns = NumImages;
nOutputs = 1;    %%%%%%%%%%   1改
%nOutputs = 5;
 
loadData=1;% 是否载入已保存的模型
 
V_thr = 1; V_rest = 0;
T = 256e-3;         % pattern duration ms
dt = 1e-3;
tau_m = 20e-3; % tau_m = 15e-3;？？？
tau_s = tau_m/4;
% K(t?ti)=V0(exp[?(t?ti)/τ]–exp[?(t?ti)/τs])
aa = exp(-(0:dt:3*tau_m)/tau_m)-exp(-(0:dt:3*tau_m)/tau_s);
 
V0 = 1/max(exp(-(0:dt:3*tau_m)/tau_m)-exp(-(0:dt:3*tau_m)/tau_s))+0;
lmd = 2e-2;%1e-2/V0;   % optimal performance lmd=3e-3*T/(tau_m*nAfferents*V0)  1e-4/V0;学习率

maxEpoch = 50;

mu = 0.99 %0.99;  % momentum factor动量
% generate patterns (each pattern consists one spik-e per afferent)
 
if loadData ==0 %初始化网络
    weights = 1e-1*randn(nAfferents,nOutputs);  % 1e-3*randn(nAfferents,1);
    save('weights0','weights');
else
    load('weights0','weights');
end
%Class = logical(eye(nOutputs));     % desired class label for each pattern
Class = false(1,2); Class(1)=true;
%Class = de2bi(1:2,'left-msb'); Class=Class';
 
correctRate=zeros(1,maxEpoch);
dw_Past=zeros(nAfferents,nPtns,nOutputs);  % momentum for accelerating learning.上一个权重的更新，用于动量计算
a = 0.5;
b = 0.5;
for epoch=1:maxEpoch    
    
    Class_Tr = false(nOutputs,nPtns);  % actual outputs of training
    for pp=1:nPtns%nPtns 
 %       Class_Tr = false(nOutputs,1);  % actual outputs of training
                
        for neuron=1:nOutputs
            Vmax=0; tmax=0;
            fired=false;        
            Vm1=zeros(1,256); indx1= 1; % trace pattern 1
            for t=dt:dt:T
                Vm = 0; 
                if fired==false
                    Tsyn=find(TrainPtns(:,pp)<=t+0.1*dt);    % no cut window
                else
                    Tsyn=find(TrainPtns(:,pp)<=t_fire+0.1*dt); % shut down inputs
                end
				
				if ~isempty(Tsyn)                    
                    A1=TrainPtns(:,pp);
                    A2=A1(Tsyn);
					
					if Class(neuron,pp) == true
						if Class_Tr(neuron,pp)==false
							K =(V0-a)*(exp(-(t-A2)/tau_m)-exp(-(t-A2)/tau_s)); % the kernel value for each fired afferent
						else
							K =(V0)*(exp(-(t-A2)/tau_m)-exp(-(t-A2)/tau_s));
						end
					elseif Class(neuron,pp) == false
						if Class_Tr(neuron,pp)==true
							K =(V0+b)*(exp(-(t-A2)/tau_m)-exp(-(t-A2)/tau_s));
						else
							K =(V0)*(exp(-(t-A2)/tau_m)-exp(-(t-A2)/tau_s));
						end
					else
						K =(V0)*(exp(-(t-A2)/tau_m)-exp(-(t-A2)/tau_s));
					end
                    A1=weights(:,neuron);
                    firedWeights=A1(Tsyn);
                    Vm = Vm + firedWeights'*K ;
				end
				
				
 
                Vm = Vm + V_rest;
                if Vm>=V_thr && fired==false % fire
                    fired=true;
                    t_fire=t;
                    Class_Tr(neuron,pp)=true;
                end
				
				if Vm>Vmax
                    Vmax=Vm; tmax=t;
                end
                
				  
 
                %if pp==1
                    Vm1(indx1)=Vm;
                    indx1=indx1+1;
                %end
            end
 
            if pp==1
					Vm2 = Vm1;
					figure(1); plot(dt:dt:T,Vm2);
					title(strcat('Image ',char('A'+pp-1),'; neuron: ',num2str(neuron))); drawnow;
            end
			if pp==2
					Vm3 = Vm1;
					figure(2); plot(dt:dt:T,Vm3);
					title(strcat('Image ',char('A'+pp-1),'; neuron: ',num2str(neuron))); drawnow;
			end
            if Vmax<=0
                tmax=max(TrainPtns(:,pp));
            end
            
            if Class_Tr(neuron,pp)~=Class(neuron,pp) % error
                
                Tsyn=find(TrainPtns(:,pp)<=tmax+0.1*dt); 
                if ~isempty(Tsyn)                    
                    A1=TrainPtns(:,pp);
                    A2=A1(Tsyn);
                    K =V0*(exp(-(tmax-A2)/tau_m)-exp(-(tmax-A2)/tau_s)); % the kernel value for each fired afferent
                    A1=weights(:,neuron);
                    dwPst=dw_Past(:,pp,neuron);
                    if fired==false    % LTP
                        Dw=lmd*K;
                    else           % LTD
                        Dw=-1.1*lmd*K;
                    end
                    A1(Tsyn) = A1(Tsyn) + Dw + mu*dwPst(Tsyn);
                    weights(:,neuron)=A1;
                    dwPst(Tsyn)=Dw;
                    dw_Past(:,pp,neuron) = dwPst;
                end                
            end            
            
        end  % end of one neuron computation
        
   end % end of one image
   %CC=isequal(Class,Class_Tr);
   %correctRate(epoch)=sum(Class==Class_Tr)/length(Class);
   CC = bi2de(Class_Tr','left-msb');
end
save('TrainedWt','weights');

% figure(3); plot(1:maxEpoch,correctRate,'-b.');




