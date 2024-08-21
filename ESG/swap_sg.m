function [a,b,c,d]=swap_sg(sg1,sg2,aa,bb,cc,dd)
% swap_sg.m
%
% PURPOSE: swap between two signal generators .
% PARAMETERS: sp1 - GPIB Instrument object already initialized and open.
%             sp2 - GPIB Instrument object already initialized and open.

if (nargin~=2 && nargin~=6)
    help swap_sg
    error('Wrong input arguments');
end

aaa=sg2.Tag;
sg2.Tag=sg1.Tag;
sg1.Tag=aaa;
a=cc;b=dd;
c=aa;d=bb;