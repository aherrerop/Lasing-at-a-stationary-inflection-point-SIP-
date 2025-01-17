%{
||
|| @file 	returnEffIndexOdd.m
|| @version	1.0
|| @author	Nathaniel Furman
|| @contact	furmann@uci.edu
|| @credit
|| | 
|| #
||
|| @description
|| | Returns the even mode effective refractive index for a coupled
|| | 450x220 nm Si with SiO2 cladding OWG and a 150 nm gap based on given
|| | frequency parameters.
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

function n_odd = returnEffIndexOdd(freq_min, freq_max, freq_bin)

if freq_min < 190 || freq_max > 196 % In THz
    error("Frequency not in range, must select between 190-196 THz")
end

% Copied from "CST-Results/Coupled-WG-Wideband-Epsilon.txt"

freq_sim = [ ...
    190.00000000000; ...
    190.58950000000; ...
    191.27110000000; ...
    192.08110000000; ...
    192.79990000000; ...
    192.99970000000; ...
    193.00000000000; ...
    193.00120000000; ...
    193.03240000000; ...
    193.94710000000; ...
    194.50000000000; ...
    195.04000000000; ...
    195.33820000000; ...
    195.59620000000; ...
    195.82510000000; ...
    195.97000000000; ...
    195.99970000000; ...
    196.00000000000; ...
];

%{
% Port one odd mode
epsilon_12 = [ ...
    5.3015638484620; ...
	5.3297431206802; ...
	5.3621181233797; ...
	5.4003052865962; ...
	5.4339338873317; ...
	5.4432383818855; ...
	5.4432523385629; ...
	5.4433081648717; ...
	5.4447594127670; ...
	5.4871046275048; ...
	5.5125125999687; ...
	5.5371918954786; ...
	5.5507630694854; ...
	5.5624719629413; ...
	5.5728348205537; ...
	5.5793824857912; ...
	5.5807233769806; ...
	5.5807464030428; ...
];

% Port two odd mode
epsilon_22 = [ ...
    5.3017276145412; ...
	5.3299060227115; ...
	5.3622800346922; ...
	5.4004660320699; ...
	5.4340936086439; ...
	5.4433978202433; ...
	5.4434117764991; ...
	5.4434676011093; ...
	5.4449188048907; ...
	5.4872627343686; ...
	5.5126699375086; ...
	5.5373484871048; ...
	5.5509192515067; ...
	5.5626277918912; ...
	5.5729903372851; ...
	5.5795378053693; ...
	5.5808786561899; ...
	5.5809035462660; ...
];
%}

% Calculate average odd mode refractive index
%n_odd_sim = sqrt((epsilon_12+epsilon_22)/2);
n_odd_sim = [ ...
   2.302530288943362; ...
   2.308641282593693; ...
   2.315642260591206; ...
   2.323872986919262; ...
   2.331097112517580; ...
   2.333091961553252; ...
   2.333094952532151; ...
   2.333106916322203; ...
   2.333417902740281; ...
   2.342473837833989; ...
   2.347890812780409; ...
   2.353140495442569; ...
   2.356022317486838; ...
   2.358505856981544; ...
   2.360701713245323; ...
   2.362088090139792; ...
   2.362371904799337; ...
   2.362376975559659; ...
];

n_odd = interp1(freq_sim, n_odd_sim, ...
    linspace(freq_min, freq_max, freq_bin), 'spline');

end

%% Footer
%{
|| @changelog
|| | 1.0 2022-04-03 - Nathaniel Furman : Initial Release
|| #
%}