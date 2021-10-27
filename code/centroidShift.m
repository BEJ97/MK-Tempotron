function [ C1 ] = centroidShift(C1, centroid, shiftHorizontal, shiftVertical)

% by default, only shift horizontally.
if ~exist('shiftHorizontal','var')
    shiftHorizontal = 1;
end
if ~exist('shiftVertical','var')
    shiftVertical = 0;
end

xc = round(centroid(1));
yc = round(centroid(2));

[numScales, numRot] = size(C1);

for i = 1:numScales
    for j = 1:numRot
        A = C1{i,j};
        [M,N] = size(A);
        
        n1 = M/2-xc;
        if shiftVertical
            B = zeros(M,N);
            if n1>=0,   % shift downwards
                B(n1+1:M, :) = A(1:M-n1, :);
            else        % shift upwards
                n1 = abs(n1);
                B(1:M-n1, :) = A(n1+1:M, :);
            end
            A = B;
        end
        
            
        n2 = N/2-yc;
        if shiftHorizontal
            B = zeros(M,N);
            if n2>=0,   % shift rightwards
                B(:, n2+1:N) = A(:, 1:N-n2);
            else        % shift leftwards
                n2 = abs(n2);
                B(:, 1:N-n2) = A(:, n2+1:N);
            end
            A = B;
        end
        
        C1{i,j} = A;
    end
end

