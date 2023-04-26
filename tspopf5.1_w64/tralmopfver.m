function rv = tralmopfver(varargin)
%TRALMOPFVER  Prints or returns TRALMOPF version info for current installation.
%   V = TRALMOPFVER returns the current TRALMOPF version number.
%   V = TRALMOPFVER('all') returns a struct with the fields Name, Version,
%   Release and Date (all strings). Calling TRALMOPFVER without assigning the
%   return value prints the version and release date of the current
%   installation of TRALMOPF.
%
%   See also MPVER.

%   TSPOPF for MATPOWER
%   $Id: tralmopfver.m 2645 2015-03-11 20:53:01Z ray $
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2007-2014 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/tspopf/ for more info.

v = struct( 'Name',     'TRALMOPF', ... 
            'Version',  '5.1', ...
            'Release',  '', ...
            'Date',     '12-Mar-2015' );
if nargout > 0
    if nargin > 0
        rv = v;
    else
        rv = v.Version;
    end
else
    fprintf('%-22s Version %-9s  %11s\n', v.Name, v.Version, v.Date);
end
