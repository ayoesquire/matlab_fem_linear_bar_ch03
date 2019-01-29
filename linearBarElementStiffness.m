function y = linearBarElementStiffness(E,A,L)
% linearBarElementStiffness     This function returns the element
%                               stiffness matrix for a linear bar with
%                               modulus of elasticity E, cross-sectional
%                               area A, and length L. The sizeof the
%                               element stiffness matrix is 2 x 2.
y = [E*A/L -E*A/L; -E*A/L E*A/L];