% example_3_1_03_280119.m
% Consider the structure composed of two linear bars as showm
% in Fig. 3.2.Given E = 210GPa, A = 0.003m2, P = 10kN, and node
% 3 is displaced to the right by 0.002m, determine:
% 1. The global stiffness matrix for the structure
% 2. The displacement at node 2
% 3. The reactions at nodes 1 and 3
% 4. The stress in e ach bar
% NOTE: Modified Steps 4, 5 & 6 for example_3_1_01_011218.m using the
% algorithm for example_4_2_01_230119.m
%%
clear
clc


% Declare Modulus of Elasticity:
E = input('Enter Modulus of Elasticity of Bar, E = ');
% If user entry is empty, it assigns a default value:
if isempty(E)
    E = 210e6;
end

% Declare no. of elements, n:
n = input('Enter of Elements, n = ');
% If user entry is empty, it assigns a default value:
if isempty(n)
    n = 1;
end

% Declare vector for varying length of bar, L, as vector:
L = ones(1,n);
for i = 1:n
    L(i) = input(['Enter Length of element ', num2str(i),' = ']);
end

% Declare cross-sectional area, A, as vector:
A = input('Enter Cross-sectional area, A (meter-sqr) = ');

% Display parameters:
disp('********************');               % New-line
disp(['Modulus of Elasticity, E = ',num2str(E),' Pa']);
disp(['Number of Elements, n = ', num2str(n)]);
disp(['Length of element(s), L(i) = ',num2str(L),' m']);
disp(['Cross-sectional Area, A = ', num2str(A),' m^2']);

%----------STEP 2: Generating the Element Stiffness Matrix, k:-------------

k = zeros(2,2*n);   % Intialize elemental stiffness matrix that will hold
                    % the subsequently generated square stiffness matrices
                    % for each element. However, note that this process is
                    % optional.
for i = 1:2:2*n
    k(1:2,i:i+1) = linearBarElementStiffness(E,A,L((i+1)/2));
end

%disp(['Element Stiffness Matrices, k = ',num2str(k)]);
disp(k);

%-----------STEP 3: Assembling the Global Stiffness Matrix---------
% The size of the global stiffness matrix is (n+1) x (n+1). Thus, intialize
% a zero matrix and make calls to the appropriate assemble function
K = zeros(n+1,n+1);
r = 1; c = 2;
for i = 1:n
    K = linearBarAssemble(K,k(1:2,r:c),i,i+1);
    r = r + 2; c = c + 2;
end

disp('K'); disp(K);

%-----------STEP 4: Applying the Boundary Conditions:----------------

% Define the Global Displacement matrix, U:
U = zeros(n+1,1);
UUnknown = U;
% Define the Global Nodal Force matrix, F:
F = zeros(n+1,1);
FUnknown = F;

% Iter
countU = 0; countF = 0; countZeros = 0;

% Prompt for known variables:

for i = 1:length(U)
    val = input(['Enter value of Global Displ, U(',num2str(i),')'...
        ' if known, or simply press "Enter" if unknown: ']);
    if isempty(val)
        UUnknown(i) = i;
        countU = countU + 1;
    else
        U(i) = val;
        if U(i) ~= 0
            countZeros = countZeros + 1;    % Flags for non-zero values
        end                                 % present
    end
end
for i = 1:length(F)
    val = input(['Enter value of Global Force, F(',num2str(i),')'...
        ' if known, or simply press "Enter" if unknown: ']);
    if isempty(val)
        FUnknown(i) = i;
        countF = countF + 1;
    else
        F(i) = val;
    end
end

fprintf('\tU\t\tF\n');
disp([U,F]);
disp('U unknowns    F unknowns')
disp([UUnknown,FUnknown]);

%-------------STEP 5: Solving the Equations:-------------------
% Using partitioning
% Test to see which parameter should be used between UUnknown and FUnknown:
% Becareful how you replace either 'UUnknown' and 'FUnknown' in subsequent code
% as 'UUnknown' and 'FUnknown' may still be relevant rather than using 'unknown'

% Initialize and populate the partition matrix:

%rem = length(unknown) - count;  % Calculates the remainder
rem = length(FUnknown) - countF;  % Calculates the remainder
disp(['No. of Remainder = ', num2str(rem)]);
fprintf('\n');     % Display newline \n Search Help keyword is 'metacharacter'

if (countF > 0) && (countZeros > 0)
    disp('Special partitioning with super-imposition required:');
    [kP,f0] = specialPartMatrix(rem, FUnknown, U, F, K);
    % Compute the corresponding partition displacement vector, uP:
    uP = kP\f0;
elseif (count > 0)
    disp('Normal parititioning required:');
    [kP,fP] = normalPartMatrix(rem, FUnknown, F, K);
    % Compute the corresponding partition displacement vector, uP:
    uP = kP\fP;
else
    disp('Error in specified boundary conditions!');
end

% Finally, obtain the corresponding partition displacement vector, :
disp('uP = '); disp(uP);

% ------------ STEP 6: Post-processing ----------------:

% ------Re-construct the original displacement vector dimension and
% integrate with the newly computed results accordingly:
uR = zeros(n+1,1);
rowP = 1;           % If you don't do this and use 'rowR' as your iter, it
                    % return an error stating that the index dimensions do
                    % not match. It is for this same reason I employed this
                    % tactic in the preceding the for loop

for i = 1:length(uR)
    if UUnknown(i) == 0
        uR(i) = U(i);
    elseif UUnknown(i) ~= 0
        uR(i) = uP(rowP);
        rowP = rowP + 1;
    end
end

disp(uR);

% ------ Complete the re-construstion by replacing U with uR:
U = uR;
disp('Global Displacement, U = '); disp(U);

% ------ Now, repeat the same to re-construct the original Global Force
% vector, by simply computing the following:

F = K*U;
disp('Global Force, F = '); disp(F);

% ------ Setup the element nodal displacement vector by making calls to
% the function linearBarElementForces().m
% NOTE: This node positions follow the i, j, m, where m represents the
% node at the middle.

c = 1;           % Reset column range index
r = 1;           % Reset iter

u = zeros(2,n);   % Initilize matrix to store elemental node variables
                  % NOTE: Quadratic elements always have 3 nodes hence the
                  % dimension of the row of the elemental matrix u is set
                  % as 3. Please take note.

% Initialize the element force vectors (f)
f = zeros(2,n);    % Same as above

% Initialize the element stress vectors (sigma)
sigma = zeros(2,n); % Same as above

% Generate the element nodal vectors (u, f) and simultaneously compute the
% element stress vectors (sigma):
disp(['Element Nodal Displacement (u), Force (f) & Stress (sigma) '...
    'vectors are given as follows: ']);

for i = 1:n
    % Compute vectors:
    u(:,i) = [U(r); U(r+1)];
    f(:,i) = linearBarElementForces(k(1:2,c:c+1),u(:,i));
    sigma(:,i) = f(:,i)/A;
    % NOTE: Alternatively, the elemental stresses may be computed by making
    % calls to the appropriate function as indicated in the commented line
    % below (Uncomment to test):
    % sigma(:,i) = linearBarElementStresses(k(1:2,c:c+1),u(:,i),A);

    % Simultaneously display results:
    disp(['u',num2str(i),' = ']); disp(u(:,i));
    disp(['f',num2str(i),' = ']); disp(f(:,i));
    disp(['sigma',num2str(i),' = ']); disp(sigma(:,i));

    c = c + 2; r = r + 1;
end
