function [A, rerr] = dchol2(A0)

% Preliminaries
tol = 1e-20; % Tolerance
A = triu(A0);
n = length(A0(:, 1));

% Compute LL' factorization
for i = 1 : n
    % Find diagonal elements
    s = 0;
    for j = 1 : i-1
        s = s + A(j,i) * A(j,i);
    end
    t = A(i,i) - s;

    if t <= tol
        % Diagonal element is algorithmically 0
        A(i,i) = 0;
        A(i, i+1:n) = 0;
    else
        % Find off-diagonal elements
        A(i,i) = sqrt(t);
        for j = i+1 : n
            s = 0;
            for k = 1 : i-1
                s = s + A(k,i) * A(k,j);
            end
            A(i,j) = A(i,j) / A(i,i) - s / A(i,i);
        end
    end
end

% Estimate error
rerr = norm(A'*A - A0);


