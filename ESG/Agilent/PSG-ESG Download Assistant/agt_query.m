function [status, status_description, query_result] = agt_query(connection,SCPI)
% PSG/ESG Download Assistant, Version 1.2
% Copyright (C) 2003 Agilent Technologies, Inc.
%
% function [status, status_description] = agt_query(connection,SCPI)
% The function sends SCPI commands to the instrument.
%
% Output:        
%   status                  integer     status of download. 0:succeeded, -1:failed
%   status_description      string      if status is < 0, status_description contains an error message.
%   query_result            string      if status is 0, query_result contains the result of the query.
% Input:
%   connection    a structure generated by the agt_newconnection function.
%   SCPI          a SCPI command. Note: use agt_sendcommand to send other SCPI commands.
%
[status, response] = agt_sgIOmx( int32(3), connection, SCPI );
if (status <0) 
   status_description = response;
   query_result = [];
else
   status_description = 'succeeded';
   query_result = response;
end
