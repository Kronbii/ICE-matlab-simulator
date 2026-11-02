function varargout = yline(y, spec, varargin)
  % Stand-in for MATLAB's yline, which GNU Octave lacks: a horizontal line at y
  % across the current x-limits. Only on the path when running in Octave.
  if nargin < 2
    spec = '-';
  end
  xl = xlim();
  held = ishold();
  hold on;
  h = plot(xl, [y y], spec, varargin{:});
  xlim(xl);
  if ~held
    hold off;
  end
  if nargout > 0
    varargout{1} = h;
  end
end
