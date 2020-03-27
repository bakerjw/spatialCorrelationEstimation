function [fun] = vario_fun_form(options)

% Specify the variogram functional form to be fit
%
% Created by Yilin Chen 02/20/2020
% Modified by Jack Baker, 3/17/2020, to update input parameter
% specification
%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Variables
%   options.funcForm = 1 for sill*(1-exp(-3 * h / range)); 
%                    = 2 for sill*(1-exp(-h^0.55 / range))
%   options.fixedSill = 1 for assuming sill is a constant of 1
%                   = 0 for assuming sill is a variable
%
% Output Variables
%   fun = variogram function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if options.funcForm == 1 && options.fixedSill == 1
    fun = @(b,x) (1 - exp(-3 .* x ./ b(1)));

elseif options.funcForm == 2 && options.fixedSill == 1
    fun = @(b,x) (1 - exp(-x.^0.55 ./ b(1)));

elseif options.funcForm == 1 && options.fixedSill == 0
    fun = @(b,x)  b(1).*(1 - exp(-3 .* x ./ b(2)));

elseif options.funcForm == 2 && options.fixedSill == 0
    fun = @(b,x)  b(1).*(1 - exp(-x.^0.55 ./ b(2)));
    
else
    warning('Error, incorrect inputs to variogram fitting function')
end

end
