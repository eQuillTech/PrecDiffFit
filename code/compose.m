% Octave compatibility shim for MATLAB's compose(): returns a char string.
function s=compose(fmt,varargin)
    s=sprintf(fmt,varargin{:});
end
