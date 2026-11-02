% Runs the simulator in GNU Octave: constants.m, part1.m and part2.m, unchanged.
% From the repository root:  octave --no-gui tools/octave/run_all.m
% Like the MATLAB scripts, it writes constants.mat, part1.mat and part2.mat to the
% current folder.

addpath(fileparts(mfilename('fullpath')));   % stand-ins for xline and yline
save_default_options('-v7');                  % .mat files MATLAB and SciPy can read
more off;

disp('>> constants'); constants
disp('>> part1'); part1
disp('>> part2');
try
  part2
catch err
  % Octave cannot write the anonymous function `equations` to a MAT file, and
  % saving is the last line of part2.m: every result exists by then.
  if isempty(strfind(err.message, 'equations'))
    rethrow(err);
  end
  clear equations err
  save('part2.mat');
end
