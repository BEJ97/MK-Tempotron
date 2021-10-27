function [ V ]= edInteg (ts,V0,dt,use_single_exponential,lut1,lut2)
% event-driven leaky integration cell

NumEvt = length(ts);  %�¼�����

K1 = 0;
K2 = 0;
V = zeros(NumEvt,1); 
t_last = -1;

for i = 1:NumEvt
    ti = ts(i);% ÿ���¼���ʱ��
    delta_t = ti-t_last;%��ǰʱ������һʱ��Ĳ�ֵ
    
    lut_addr = round(delta_t /dt)+1;
    if lut_addr<=length(lut1)
        Sc1 = lut1(lut_addr);%lut1�����m��Ӧ��ָ��˥�����㲿��
    else
        Sc1 = 0;
    end
    K1 = Sc1*K1;
    K1 = K1 + V0;
    
    if ~use_single_exponential,
        if lut_addr<=length(lut2)
            Sc2 = lut2(lut_addr);
        else
            Sc2 = 0;
        end
        K2 = Sc2*K2;
        K2 = K2 + V0;
        K = K1-K2;
    else
        K = K1;
    end
    
    V(i) = K;             % total potential. (i.e. the potential in Motion symbol detector)��i���¼�����ʱ��Ĥ��λ
    
    t_last = ti;
end

