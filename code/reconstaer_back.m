function [J]=reconstaer_back(CIN, size1, size2)
% AER--> image, discard sign of events from the very beginning.
% size1: #rows;  size2: #cols

J = zeros(size1,size2);
if size(CIN,1)>0
    x=CIN(:,4);
    y=CIN(:,5);
    signo=CIN(:,6);
	
    for i=1:length(x)
		J(x(i),y(i)) = J(x(i),y(i)) + abs(signo(i));
    end
end
