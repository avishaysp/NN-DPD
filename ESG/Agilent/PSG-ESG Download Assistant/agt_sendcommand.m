function [status, status_description] = agt_sendcommand(connection,SCPI)
% PSG/ESG Download Assistant, Version 1.2
% Copyright (C) 2003 Agilent Technologies, Inc.
%
% function [status, status_description] = agt_sendcommand(connection,SCPI)
% The function sends SCPI commands to the instruemnt.
%
% Output:        
%   status               integer     status of sending a SCPI string. 0:succeeded, -1:failed
%                                    NOTE:  status only returns a status on whether or not the SCPI
%                                    string was sent to the signal generator.  It does not report
%                                    whether or not the SCPI string was valid.
%   status_description   string      if status is < 0, status_description contains an error message.
% Input:
%   connection    a structure generated by the agt_newconnection function.
%   SCPI          a SCPI command. Note: use agt_query to send query command.
%
[status, status_description] = agt_sgIOmx( int32(2), connection, SCPI );