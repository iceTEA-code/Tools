%
% Ptotal = PR_Z(z_in,pp,sf,cp,nuclide)
%
% Calculates total production rates (from spallation and muon pathways) for
% given depths based on computed sample parameters, using CronusCalc 
% functions.
%
% z_in should be a depth or array of depths (in g/cm^2).
%
% pp is a struct of physical constant, generated using get_pars.m.
%
% sf is a struct of scaling factors, generated for each sample using 
% get_pars.m.
%
% cp is a struct of computed parameters, generated for each sample using 
% get_pars.m.
%
% nuclide should 10 or 26, corresponding to 10Be or 26Al. 
%
% Output is the total production rate at given depths for the given 
% nuclide, for all samples.
%
% Written by Richard Selwyn Jones, Durham University
% richard.s.jones@durham.ac.uk
% Part of the iceTEA tools suite, which is built on versions of 
% CRONUS-Earth and CRONUScalc code.
%
%
%%

function Ptotal = PR_Z(z_in,pp,sf,cp,nuclide)

  % Check inputs
  if (nargin ~= 5)
      error('PR_Z has wrong number of inputs!');
  end
  if nuclide ~= 10 && nuclide ~= 26
      error('nuclide should be "10" or "26"!');
  end

  z = z_in;

  if nuclide == 10 || nuclide == 26
  
      [PtotalBe,PtotalAl,~,~,~,~] = prodz1026(z,pp,sf,cp);
      
      if nuclide == 10
          Ptotal = PtotalBe;
      elseif nuclide == 26
          Ptotal = PtotalAl;
      end
  
  end
  
end
