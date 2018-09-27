%
% out = get_jags()
%
% Installs JAGS (Just Another Gibbs Sampler). To be used with
% transect_regress_spline.m, before the JAGS model is run.
%
% No inputs are required.
%
% Returns 1 if JAGS is installed, or 0 if it is not installed.
%
% Written by Richard Selwyn Jones, Durham University
% richard.s.jones@durham.ac.uk
% Part of the iceTEA tools suite.
%
%
%%

function out = get_jags()
  
  warning('JAGS (Just Another Gibbs Sampler) is required for statistical analysis.');
  choice1 = questdlg('Unable to find the JAGS program. Should I try to download it now from the internet?',...
      'Download JAGS now?', ...
      'Yeah, please download.','No thanks.','Yeah, please download.');
  switch choice1
      case 'Yeah, please download.'
          
          if ispc
              try
                  websave('JAGS-4.3.0.exe','https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/Windows/JAGS-4.3.0.exe');
                  disp('JAGS-4.3.0.exe has downloaded to your current directory. Run the JAGS executable, follow the prompts and install in the iceTEA Functions directory.')
                  out = 1;
              catch
                  warning('Download failed. This may be caused by firewall issues, or some other reason. Manually download the program here: https://sourceforge.net/projects/mcmc-jags/. Once downloaded, follow the prompts and install in the iceTEA Functions directory.')
                  out = 0;
              end
          elseif ismac
              try
                  websave('JAGS-4.3.0.dmg','https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/Mac%20OS%20X/JAGS-4.3.0.dmg');
                  disp('JAGS-4.3.0.dmg has downloaded to your current directory. Run the JAGS disk image, follow the prompts and install on the system.')
                  out = 1;
              catch
                  warning('Download failed. This may be caused by firewall issues, or some other reason. Manually download the program here: https://sourceforge.net/projects/mcmc-jags/. Once downloaded, follow the prompts and install on the system.')
                  out = 0;
              end
          else
              disp('Not PC or Mac?  You can try downloading other JAGS versions here: http://mcmc-jags.sourceforge.net/.')
              out = 0;
          end
          
      otherwise
          disp('Unable to continue without JAGS.  Try downloading manually from here: https://sourceforge.net/projects/mcmc-jags/. Once downloaded, follow the prompts and install in the iceTEA Functions directory.')
          out = 0;
  end

end
