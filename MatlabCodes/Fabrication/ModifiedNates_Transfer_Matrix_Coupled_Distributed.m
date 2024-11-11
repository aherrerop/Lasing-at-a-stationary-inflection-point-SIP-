%{
||
|| @file 	Transfer_Matrix_Coupled_Distributed.m
|| @version	1.0
|| @author	Nathaniel Furman
|| @contact	furmann@uci.edu
|| @credit
|| | 
|| #
||
|| @description
|| | Calculates the transfer matrix for a distributed directional coupler
|| | given the even/odd mode wavenumbers and the coupler length.
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

function [T] = ModifiedNates_Transfer_Matrix_Coupled_Distributed(Length, k_even, k_odd)

% Even and odd mode propagations
Even_Propagation = exp(1i*k_even*Length);
Odd_Propagation  = exp(1i*k_odd *Length);

% Build [S] matrix
S = zeros(4);

% b_e(0) = S13 * a_e(d)
S(1,3) = exp(-1i*k_even*Length);
% b_o(0) = S24 * a_o(d)
S(2,4) = exp(-1i*k_odd *Length);

% b_e(d) = S31 * a_e(0)
S(3,1) = exp(1i*k_even*Length);
% b_o(d) = S42 * a_o(0)
S(4,2) = exp(1i*k_odd *Length);

% We see [S] is symmetric and unitary here

% Convert to proper transfer matrix
T = ModifiedNates_Convert_Seo_to_Tfb(S);

end

%% Footer
%{
|| @changelog
|| | 1.0 2022-03-22 - Nathaniel Furman : Initial Release
|| #
%}