function [ V ]= edInteg (ts,V0,dt,use_single_exponential,lut1,lut2)
% event-driven leaky integration cell

NumEvt = length(ts);  %事件数量

K1 = 0;
K2 = 0;
V = zeros(NumEvt,1); 
t_last = -1;

for i = 1:NumEvt
    ti = ts(i);% 每个事件的时间
    delta_t = ti-t_last;%当前时间与上一时间的差值
    
    lut_addr = round(delta_t /dt)+1;
    if lut_addr<=length(lut1)
        Sc1 = lut1(lut_addr);%lut1代表τm对应的指数衰减计算部分
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
    
    V(i) = K;             % total potential. (i.e. the potential in Motion symbol detector)第i个事件进入时的膜电位
    
    t_last = ti;
end

