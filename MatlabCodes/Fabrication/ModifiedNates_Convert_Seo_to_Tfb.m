%{
||
|| @file 	Convert_Seo_to_Tfb.m
|| @version	1.0
|| @author	Nathaniel Furman
|| @contact	furmann@uci.edu
|| @credit
|| | 
|| #
||
|| @description
|| | Converts 4-port (4x4xM) scattering matrix in terms of even and odd 
|| | modes to the transfer matrix in terms of forward and backward waves 
|| | in each waveguide. 'M' is the length of the scattering matrix
|| | corresponding to the number of frequency points included.
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

% Problem setup
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

function T_Parameters = ModifiedNates_Convert_Seo_to_Tfb(S_Parameters)

S_size = size(S_Parameters);
T_Parameters = zeros(S_size);

% Determine if single frequency or multiple (prevent throwing error)
numElements = size(S_size);
if numElements(2) < 3
    % Means transfer matrix is only for a single frequency
    freq_bin = 1;
else
    % Otherwise calculate for frequency range
    freq_bin = S_size(3);
end

% Convert for each frequency
for jj = 1:freq_bin
    % Convert to compacted notation
    s11 = S_Parameters(1,1,jj);
    s12 = S_Parameters(1,2,jj);
    s13 = S_Parameters(1,3,jj);
    s14 = S_Parameters(1,4,jj);

    s21 = S_Parameters(2,1,jj);
    s22 = S_Parameters(2,2,jj);
    s23 = S_Parameters(2,3,jj);
    s24 = S_Parameters(2,4,jj);

    s31 = S_Parameters(3,1,jj);
    s32 = S_Parameters(3,2,jj);
    s33 = S_Parameters(3,3,jj);
    s34 = S_Parameters(3,4,jj);

    s41 = S_Parameters(4,1,jj);
    s42 = S_Parameters(4,2,jj);
    s43 = S_Parameters(4,3,jj);
    s44 = S_Parameters(4,4,jj);
    
    % Solutions calculated using Symbolic maths  
    t11 = 0;
    t21 = -(s31 + s42)/(s31 - s42);
    t31 = -(2*s31*s42)/(s31 - s42);
    t41 = 0;
    t12 = (s13 + s24)/(s13 - s24);
    t22 = 0;
    t32 = 0;
    t42 = 2/(s13 - s24);
    t13 = -(2*s13*s24)/(s13 - s24);
    t23 = 0;
    t33 = 0;
    t43 = -(s13 + s24)/(s13 - s24);
    t14 = 0;
    t24 = 2/(s31 - s42);
    t34 = (s31 + s42)/(s31 - s42);
    t44 = 0;
    
    T_Parameters(1,1,jj) = t11;
    T_Parameters(1,2,jj) = t12;
    T_Parameters(1,3,jj) = t13;
    T_Parameters(1,4,jj) = t14;
    
    T_Parameters(2,1,jj) = t21;
    T_Parameters(2,2,jj) = t22;
    T_Parameters(2,3,jj) = t23;
    T_Parameters(2,4,jj) = t24;
    
    T_Parameters(3,1,jj) = t31;
    T_Parameters(3,2,jj) = t32;
    T_Parameters(3,3,jj) = t33;
    T_Parameters(3,4,jj) = t34;
    
    T_Parameters(4,1,jj) = t41;
    T_Parameters(4,2,jj) = t42;
    T_Parameters(4,3,jj) = t43;
    T_Parameters(4,4,jj) = t44;
    
%     T_Parameters(:,:,jj) = T_Parameters(:,:,jj)./(s11*s33 - s13*s31 - s11*s34 - s12*s33 + s13*s32 + s14*s31 + s12*s34 - s14*s32 + s11*s43 - s13*s41 + s21*s33 - s23*s31 - s11*s44 - s12*s43 + s13*s42 + s14*s41 - s21*s34 - s22*s33 + s23*s32 + s24*s31 + s12*s44 - s14*s42 + s22*s34 - s24*s32 + s21*s43 - s23*s41 - s21*s44 - s22*s43 + s23*s42 + s24*s41 + s22*s44 - s24*s42);

end

end

%% Footer
%{
|| @changelog
|| | 1.0 2022-03-21 - Nathaniel Furman : Initial Release
|| #
%}