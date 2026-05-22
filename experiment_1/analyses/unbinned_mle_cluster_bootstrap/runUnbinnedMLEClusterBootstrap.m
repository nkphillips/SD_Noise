function varargout = runUnbinnedMLEClusterBootstrap(varargin)
% runUnbinnedMLEClusterBootstrap  Deprecated compatibility wrapper.
% Use runSerialDependenceMLEBootstrap instead.

    warning('runUnbinnedMLEClusterBootstrap:deprecated', ...
        ['runUnbinnedMLEClusterBootstrap has moved to ' ...
         'runSerialDependenceMLEBootstrap. Forwarding this call.']);

    wrapper_dir = fileparts(mfilename('fullpath'));
    analysis_dir = fullfile(wrapper_dir, '..');
    addpath(fullfile(analysis_dir, 'functions', 'serial_dependence_mle'));
    addpath(fullfile(analysis_dir, 'plotting', 'serial_dependence_mle'));

    if nargout == 0
        runSerialDependenceMLEBootstrap(varargin{:});
    else
        [varargout{1:nargout}] = runSerialDependenceMLEBootstrap(varargin{:});
    end
end
