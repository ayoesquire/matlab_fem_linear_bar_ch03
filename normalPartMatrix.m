function [x, y] = normalPartMatrix(rem, FUnknown, F, K)

%%
    % Function to perform regular matrix partition when zero or no
    % nodal displacement value(s) {U(i)} is specified as a boundary
    % condition
%%

kP = zeros(rem);                        % Initialize and pre-allocate new
fP = zeros(rem,1);                      % partition matrices

rP = 1; cP = 1;
for i = 1:length(FUnknown)
    if FUnknown(i) == 0
        % disp(i);
        for j = 1:length(FUnknown)
            if FUnknown(j) == 0
                disp([num2str(i),',',num2str(j)]);
                kP(rP,cP) = K(i,j);   % Assign element corresponding to
                cP = cP + 1;          % global matrix
            end
        end
        % Repeat for fP:
        fP(rP,1) = F(i,1);
        cP = 1;           % Reset index for column while
        rP = rP + 1;       % index for row is incremented
    end
end

% Display result of partition:
disp('kP = '); disp(kP);
disp('fP = '); disp(fP);

% Return values:
x = kP;
y = fP;