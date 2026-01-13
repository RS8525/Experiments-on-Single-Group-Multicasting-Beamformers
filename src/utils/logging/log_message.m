function log_message(message, varargin)
% LOG_MESSAGE Print timestamped log message
%
% Inputs:
%   message  - String message to log
%   varargin - Optional arguments for sprintf formatting
%
% Example:
%   log_message('Starting experiment...');
%   log_message('Iteration %d: error = %.6f', iter, error);

timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');

if nargin > 1
    message = sprintf(message, varargin{:});
end

fprintf('[%s] %s\n', timestamp, message);

end
