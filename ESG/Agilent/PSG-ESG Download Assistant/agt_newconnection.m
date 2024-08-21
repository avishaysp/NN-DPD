function SGConnection = agt_newconnection(varargin)
% PSG/ESG Download Assistant, Version 1.2
% Copyright (C) 2003 Agilent Technologies, Inc.
%
% function SGConnection = agt_newconnection(varargin)
% This function generates connection structures according to the different interface types.
% The valid interface types are: 'GPIB' or 'TCPIP'.
% 
% The function call for each interface is:
% 
% For GPIB connection:
%     agt_newconnection('GPIB',boardindex, primaryaddress, secondaryaddress)
%     boardindex       integer default 0 
%     primaryaddress   integer default 19
%     secondaryaddress integer default 0.  (optional parameter)
%     sample: connection = agt_newconnection('gpib',0,19);
% For TCPIP connection:
%     agt_newconnection('TCPIP',IPAddress)
%     IPAddress  string 
%     sample: connection = agt_newconnection('TCPIP','11.11.11.11');

if (nargin < 1) 
    SGConnection = default;
    return;
end

if (~ischar(varargin{1}))
   error('the first parameter has to be a string showing the connection interface type: ''GPIB'' or ''TCPIP''');
end

interface = varargin{1};

if (strcmpi(interface,'GPIB'))
    SGConnection = GPIBConnection(varargin{2:end});
    return;
end

if (strcmpi(interface,'TCPIP'))
    SGConnection =  TCPIPConnection(varargin{2:end});
    return;
end

msg = sprintf('%s%s%s','Invalid interface type: ''',interface,'''. The interface type should be one of ''GPIB'' or ''TCPIP''');
error(msg);


function SGConnection = default()
    SGConnection.type = 'GPIB';
    SGConnection.GPIB.board   = int32 ( 0 );
    SGConnection.GPIB.primary  = int32 ( 19 );
    SGConnection.GPIB.secondary = int32 ( 0 );
    
function SGConnection = GPIBConnection(varargin)
    
    SGConnection.type = 'GPIB';
    SGConnection.GPIB.board   = int32 ( 0 );
    SGConnection.GPIB.primary  = int32 ( 19 );
    SGConnection.GPIB.secondary = int32 ( 0 );

    if (nargin ==0) return; end
    if (nargin >0) 
        SGConnection.GPIB.board = int32 ( checkPositiveInteger(varargin{1},'GPIB board number') );
    end
    
    if (nargin >1) 
        SGConnection.GPIB.primary = int32 ( checkPositiveInteger(varargin{2},'GPIB primary address') );
    end
    
    if (nargin >2) 
        SGConnection.GPIB.secondary = int32 ( checkPositiveInteger(varargin{3},'GPIB secondary address') );
    end
            
function SGConnection = TCPIPConnection(varargin)
    SGConnection.type = 'TCPIP';
    SGConnection.TCPIP.board = int32 ( 0 );
    SGConnection.TCPIP.address = '0.0.0.0';
    if (nargin ==0) return; end
    
    if (nargin >0) 
        SGConnection.TCPIP.address = checkString(varargin{1},'TCPIP address');
    end

function val = checkPositiveInteger(input, errormsgHeader)

    val = input;
    if (~isnumeric(input))
        msg = strcat(errormsgHeader, ' should be an integer');
        error(msg);
    end
    val = abs(floor(val(1)));

function val = checkString(input,errormsgHeader)
    val = input;
    if (~ischar(input))
        msg = strcat(errormsgHeader, ' should be a string');
        error(msg);
    end