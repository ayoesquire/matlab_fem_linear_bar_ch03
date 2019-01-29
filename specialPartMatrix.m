function [x, y] = specialPartMatrix(rem, FUnknown, U, F, K)

%%
    % Function to perform special matrix partition when at least one
    % non-zero nodal displacement value(s) {U(i)} is specified as a
    % boundary condition
%%

k0 = zeros(rem,length(U));           % Initialize the new destination matrix
fP = zeros(rem,1);         % (Ditto)

rP = 1;

for i = 1:length(U)% Populate partition matrix
    if FUnknown(i) == 0
        k0(rP,:) = K(i,:); % Assign element corresponding to global matrix
        disp(FUnknown(i)); disp(i);

        % Repeat for fP:
        fP(rP,1) = F(i,1);

        rP = rP + 1;       % index for row is incremented
    end
end

% Extract sub-matrix by way of product (i.e. the matrix multiplication)
% of kP and U then,

% Super-impose the produc and solve the equation:

f0 = fP - (k0*U);

% Not forgetting to evaluate kP by calling normalPartMatrix() function:
kP = normalPartMatrix(rem, FUnknown, F, K);

% Display result of partition (Note: no need to display kP and fP as the
% normalPartMatrix() does this already):
disp('k0 = '); disp(k0);
disp('f0 = '); disp(f0);


% Return values:
x = kP;
y = f0;
