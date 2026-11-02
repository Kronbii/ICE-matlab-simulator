function varargout = xline(x, spec, varargin)
  % Stand-in for MATLAB's xline, which GNU Octave lacks: a vertical line at x
  % across the current y-limits. Only on the path when running in Octave.
  if nargin < 2
    spec = '-';
  end
  yl = ylim();
  held = ishold();
  hold on;
  h = plot([x x], yl, spec, varargin{:});
  ylim(yl);
  if ~held
    hold off;
  end
  if nargout > 0
    varargout{1} = h;
  end
end
