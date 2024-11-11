%{
||
|| @file 	Solve_EvenOdd_to_ForwardBack.m
|| @version	2.0
|| @author	Nathaniel Furman
|| @contact	furmann@uci.edu
|| @credit
|| | 
|| #
||
|| @description
|| | Sets up and solves the symbolic equations for converting between even
|| | and odd CST [S] parameters to forward and backward [T] parameters for
|| | 'Converstion_Seo_to_Tfb.m' script.
|| #
||
|| @license
|| | This library is free software; you can redistribute it and/or
|| | modify it under the terms of the GNU Lesser General Public
|| | License as published by the Free Software Foundation; version
|| | 2.1 of the License.
|| |
|| | This library is distributed in the hope that it will be useful,
|| | but WITHOUT ANY WARRANTY; without even the implied warranty of
|| | MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
|| | Lesser General Public License for more details.
|| |
|| | You should have received a copy of the GNU Lesser General Public
|| | License along with this library; if not, write to the Free Software
|| | Foundation Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110 USA
|| #
||
%}

clear;clc;

debug = true;

% Problem setup (coped from 'Conversion_Seo_to_Tfb.m' on 2022-03-21)
%{
b = S a, where e/o mean even and odd modes and 0/d are the locations of
each port
[ b_e(0)  = [S] [ a_e(0)
  b_o(0)          a_o(0)
  b_e(d)          a_e(d)
  b_o(d) ]        a_o(d) ]

-----

We need state vector as:
[E1(d); E2(d)] = [T_coupler] [E1(0); E2(0)], where 1/2 mean the respective
waveguide and +/- are the forward and backward directions of propagation
[ E1+(d)   = [T] [ E1+(0)
  E1-(d)           E1-(0)
  E2+(d)           E2+(0)
  E2-(d) ]         E2-(0) ]

-----

We also know the relationship between even/odd modes and fields in the
individual waveguides:
a_e(0) = E1+(0) + E2+(0)
a_o(0) = E1+(0) - E2+(0)

And similarly for the other three combination sets

%}

solvingLossless = true;

% Set up the symbolic variables used in the calculation
% For scattering symbolic variables:
%  - a/b represents the incident/reflected waves
%  - e/o represents the even and odd modes
%  - 0/D represents the port locations
% For waveguide electric fields:
%  - 1/2 represents the respective waveguide
%  - p/m represents the forward (plus) or backward (minus) direction
%  - 0/D represents the port locations
syms s11 s12 s13 s14 s21 s22 s23 s24 s31 s32 s33 s34 s41 s42 s43 s44 ...
    ae0 ao0 aeD aoD be0 bo0 beD boD e1p_0 e1m_0 e2p_0 e2m_0 e1p_1 e1m_1 e2p_1 e2m_1

% Given scattering matrix
s = [s11 s12 s13 s14; s21 s22 s23 s24; s31 s32 s33 s34; s41 s42 s43 s44];

% Scattering matrix equation
eqnS = [be0; bo0; beD; boD] == s * [ae0; ao0; aeD; aoD];

% Substitue relationship between even/odd modes and forward/backward modes
eqnT = subs(eqnS, [ae0, ao0, aeD, aoD, be0, bo0, beD, boD], ...
    [(e1p_0 + e1m_1), (e1p_0 - e1m_1), (e2p_0 + e2m_1), (e2p_0 - e2m_1), ...
     (e1m_0 + e1p_1), (e1m_0 - e1p_1), (e2m_0 + e2p_1), (e2m_0 - e2p_1)]);
 
% Solve for the electric fields at 'D'
mySols = solve(eqnT, [e1p_1 e1m_1 e2p_1 e2m_1]);

% Combine for general solution
outGeneral = [ ...
    simplify(mySols.e1p_1); ...
    simplify(mySols.e1m_1); ...
    simplify(mySols.e2p_1); ...
    simplify(mySols.e2m_1)];

% If using ideal case, no need for common divisor
if solvingLossless
    common_divisor = 1;
else
    common_divisor = (s11*s33 - s13*s31 - s11*s34 - s12*s33 + s13*s32 + s14*s31 + s12*s34 - s14*s32 + s11*s43 - s13*s41 + s21*s33 - s23*s31 - s11*s44 - s12*s43 + s13*s42 + s14*s41 - s21*s34 - s22*s33 + s23*s32 + s24*s31 + s12*s44 - s14*s42 + s22*s34 - s24*s32 + s21*s43 - s23*s41 - s21*s44 - s22*s43 + s23*s42 + s24*s41 + s22*s44 - s24*s42);
end

outGeneralsimplified = simplify(subs(outGeneral(1),[e1p_0 e1m_0 e2p_0 e2m_0]));

%
% Simplify and factor out common denominator
outGeneralIdeal = simplify(outGeneral.*common_divisor);

if solvingLossless
    outGeneralSimplified = simplify(subs(outGeneralIdeal, ...
        [s11 s12 s14 s22 s24 s31 s33  s42 s44], ...
        [0   0   0   0   0   0   0   0]));
else
    outGeneralSimplified = outGeneralIdeal;
end

% Collect the sixteen transfer matrix values per the input/output
% Only the terms in the first set of parenthesis correspond to T(i,j)
% outGeneral = [T] * [e1p0; e1m0; e2p0; e2m0]
generalPattern = "(" + wildcardPattern + ")*";

% T(1,1)
tmpVar = string(collect(outGeneralSimplified(1), e1p_0));
if debug, disp(tmpVar), end
tmpVar = extract(tmpVar, generalPattern);
if ~isempty(tmpVar)
    t11 = extractBetween(tmpVar(1), 2, strlength(tmpVar(1))-2);
else
    t11 = "0";
end
fprintf("t11 = %s;\n", t11)
% T(2,1)
tmpVar = string(collect(outGeneralSimplified(2), e1p_0));
if debug, disp(tmpVar), end
tmpVar = extract(tmpVar, generalPattern);
if ~isempty(tmpVar)
    t21 = extractBetween(tmpVar(1), 2, strlength(tmpVar(1))-2);
else
    t21 = "0";
end
fprintf("t21 = %s;\n", t21)
% T(3,1)
tmpVar = string(collect(outGeneralSimplified(3), e1p_0));
if debug, disp(tmpVar), end
tmpVar = extract(tmpVar, generalPattern);
if ~isempty(tmpVar)
    t31 = extractBetween(tmpVar(1), 2, strlength(tmpVar(1))-2);
else
    t31 = "0";
end
fprintf("t31 = %s;\n", t31)
% T(4,1)
tmpVar = string(collect(outGeneralSimplified(4), e1p_0));
if debug, disp(tmpVar), end
tmpVar = extract(tmpVar, generalPattern);
if ~isempty(tmpVar)
    t41 = extractBetween(tmpVar(1), 2, strlength(tmpVar(1))-2);
else
    t41 = "0";
end
fprintf("t41 = %s;\n", t41)

% T(1,2)
tmpVar = string(collect(outGeneralSimplified(1), e1m_0));
if debug, disp(tmpVar), end
tmpVar = extract(tmpVar, generalPattern);
if ~isempty(tmpVar)
    t12 = extractBetween(tmpVar(1), 2, strlength(tmpVar(1))-2);
else
    t12 = "0";
end
fprintf("t12 = %s;\n", t12)
% T(2,2)
tmpVar = string(collect(outGeneralSimplified(2), e1m_0));
if debug, disp(tmpVar), end
tmpVar = extract(tmpVar, generalPattern);
if ~isempty(tmpVar)
    t22 = extractBetween(tmpVar(1), 2, strlength(tmpVar(1))-2);
else
    t22 = "0";
end
fprintf("t22 = %s;\n", t22)
% T(3,2)
tmpVar = string(collect(outGeneralSimplified(3), e1m_0));
if debug, disp(tmpVar), end
tmpVar = extract(tmpVar, generalPattern);
if ~isempty(tmpVar)
    t32 = extractBetween(tmpVar(1), 2, strlength(tmpVar(1))-2);
else
    t32 = "0";
end
fprintf("t32 = %s;\n", t32)
% T(4,2)
tmpVar = string(collect(outGeneralSimplified(4), e1m_0));
if debug, disp(tmpVar), end
tmpVar = extract(tmpVar, generalPattern);
if ~isempty(tmpVar)
    t42 = extractBetween(tmpVar(1), 2, strlength(tmpVar(1))-2);
else
    t42 = "0";
end
fprintf("t42 = %s;\n", t42)

% And repeat the above pattern
tmpVar = string(collect(outGeneralSimplified(1), e2p_0));
if debug, disp(tmpVar), end
tmpVar = extract(tmpVar, generalPattern);
if ~isempty(tmpVar)
    t13 = extractBetween(tmpVar(1), 2, strlength(tmpVar(1))-2);
else
    t13 = "0";
end
fprintf("t13 = %s;\n", t13)
tmpVar = string(collect(outGeneralSimplified(2), e2p_0));
if debug, disp(tmpVar), end
tmpVar = extract(tmpVar, generalPattern);
if ~isempty(tmpVar)
    t23 = extractBetween(tmpVar(1), 2, strlength(tmpVar(1))-2);
else
    t23 = "0";
end
fprintf("t23 = %s;\n", t23)
tmpVar = string(collect(outGeneralSimplified(3), e2p_0));
if debug, disp(tmpVar), end
tmpVar = extract(tmpVar, generalPattern);
if ~isempty(tmpVar)
    t33 = extractBetween(tmpVar(1), 2, strlength(tmpVar(1))-2);
else
    t33 = "0";
end
fprintf("t33 = %s;\n", t33)
tmpVar = string(collect(outGeneralSimplified(4), e2p_0));
if debug, disp(tmpVar), end
tmpVar = extract(tmpVar, generalPattern);
if ~isempty(tmpVar)
    t43 = extractBetween(tmpVar(1), 2, strlength(tmpVar(1))-2);
else
    t43 = "0";
end
fprintf("t43 = %s;\n", t43)

tmpVar = string(collect(outGeneralSimplified(1), e2m_0));
if debug, disp(tmpVar), end
tmpVar = extract(tmpVar, generalPattern);
if ~isempty(tmpVar)
    t14 = extractBetween(tmpVar(1), 2, strlength(tmpVar(1))-2);
else
    t14 = "0";
end
fprintf("t14 = %s;\n", t14)
tmpVar = string(collect(outGeneralSimplified(2), e2m_0));
if debug, disp(tmpVar), end
tmpVar = extract(tmpVar, generalPattern);
if ~isempty(tmpVar)
    t24 = extractBetween(tmpVar(1), 2, strlength(tmpVar(1))-2);
else
    t24 = "0";
end
fprintf("t24 = %s;\n", t24)
tmpVar = string(collect(outGeneralSimplified(3), e2m_0));
if debug, disp(tmpVar), end
tmpVar = extract(tmpVar, generalPattern);
if ~isempty(tmpVar)
    t34 = extractBetween(tmpVar(1), 2, strlength(tmpVar(1))-2);
else
    t34 = "0";
end
fprintf("t34 = %s;\n", t34)
tmpVar = string(collect(outGeneralSimplified(4), e2m_0));
if debug, disp(tmpVar), end
tmpVar = extract(tmpVar, generalPattern);
if ~isempty(tmpVar)
    t44 = extractBetween(tmpVar(1), 2, strlength(tmpVar(1))-2);
else
    t44 = "0";
end
fprintf("t44 = %s;\n", t44)

%% Footer
%{
|| @changelog
|| | 1.0 2022-03-21 - Nathaniel Furman : Initial Release
|| | 2.0 2022-04-19 - Nathaniel Furman : Simplifed for ideal case
|| #
%}