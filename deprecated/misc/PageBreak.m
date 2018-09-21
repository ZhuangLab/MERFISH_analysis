function breakString = PageBreak(varargin)
    breakString = '-------------------------------------------------------------------------';
    if nargin < 1 || ~strcmp(varargin{1}, 'nodisplay')
        disp(breakString);
    end        
end